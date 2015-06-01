// Xerus - A General Purpose Tensor Library
// Copyright (C) 2014-2015 Benjamin Huber and Sebastian Wolf. 
// 
// Xerus is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
// 
// Xerus is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Affero General Public License for more details.
// 
// You should have received a copy of the GNU Affero General Public License
// along with Xerus. If not, see <http://www.gnu.org/licenses/>.
//
// For further information on Xerus visit https://libXerus.org 
// or contact us at contact@libXerus.org.

#pragma once

#include "contractionHeuristic.h"
#include "fullTensor.h"
#include "index.h"
#include "misc/missingFunctions.h"
#include "misc/test.h"

namespace xerus {

    template<bool isOperator>
	/// The TTNetwork class is used to represent TTTensor and TToperators (depending on the template argument) and is a special kind of TensorNetwork.
    class TTNetwork : public TensorNetwork {
	protected:
		static constexpr size_t N=isOperator?2:1;
		
		enum cannonicalization_t {
			NONE, LEFT, RIGHT
		};
		cannonicalization_t cannonicalization; 
    public:
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Constructors - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
		/// Constructs an order zero TTNetwork.
        explicit TTNetwork() = default;
        
        /// Constructs an zero initialized TTNetwork with the given degree and ranks all equal to one. Naturally for TTOperators the degree must be even.
        explicit TTNetwork(const size_t _degree);
        
		/// Constructs a TTNetwork from the given FullTensor, using the higher order SVD algorithm. Opionally an accuracy can be given.
        explicit TTNetwork(const FullTensor& _full, const double _eps=1e-15); //TODO no magic numbers
        
		/// Copy constructor for TTNetworks.
        implicit TTNetwork(const TTNetwork & _cpy);
        
		/// Move constructor for TTNetworks.
        implicit TTNetwork(      TTNetwork&& _mov);
        
		/// Transforms a given TensorNetwork to a TTNetwork. NOTE this is not yet implemented different from casting to FullTensor and then using a HOSVD.
        explicit TTNetwork(const TensorNetwork &_cpy, double _eps=1e-14);
        
		/// Transforms a given TensorNetwork to a TTNetwork. NOTE this is not yet implemented different from casting to FullTensor and then using a HOSVD.
        explicit TTNetwork(TensorNetwork &&_mov, double _eps=1e-14);
        
		/// Random constructs a TTNetwork with the given dimensions and ranks. The entries of the componend tensors are sampled independendly using the provided random generator and distribution.
        template<class generator, class distribution>
        static TTNetwork construct_random(const std::vector<size_t>& _dimensions, const std::vector<size_t> &_ranks, generator& _rnd, distribution& _dist) {
            const size_t N = isOperator?2:1;
            REQUIRE(_ranks.size() == _dimensions.size()/N-1,"non-matching amount of ranks given to TTNetwork::construct_random");
            #ifndef DISABLE_RUNTIME_CHECKS_
                for (size_t d : _dimensions) {
                    REQUIRE(d > 0, "trying to construct random TTTensor with dimension 0");
                }
                for (size_t d : _ranks) {
                    REQUIRE(d > 0, "trying to construct random TTTensor with rank 0");
                }
            #endif
            
            TTNetwork result;
            if(_dimensions.size() == 0) {
    // 			std::shared_ptr<Tensor> newT(new FullTensor(0));
    // 			newT->data.get()[0] = dis(rnd);
    // 			result.nodes.emplace_back(std::move(newT));
                result.factor = _dist(_rnd); //TODO create testcase for degree-0 TTtensors
                result.dimensions.clear();
            } else {
                result.nodes.clear();
                result.dimensions = _dimensions;
                result.factor = 1.0;
                
                // externalLinks
                size_t numNodes = _dimensions.size()/N;
                result.externalLinks.emplace_back(0, 0, result.dimensions[0], false);
                for(size_t i = 1; i < numNodes; ++i) {
                    result.externalLinks.emplace_back(i, 1, result.dimensions[i], false);
                }
                if(N == 2) {
                    result.externalLinks.emplace_back(0, 1, result.dimensions[numNodes], false);
                    for(size_t i = 1; i < numNodes; ++i) {
                        result.externalLinks.emplace_back(i, 2, result.dimensions[numNodes+i], false);
                    }
                }
                
                size_t maxDim1 = misc::product(_dimensions);
                size_t maxDim2 = 1;
                std::vector<size_t> constructVector;
                std::vector<TensorNode::Link> neighbors;
                for (size_t i=0; i<_dimensions.size()/N; ++i) {
                    constructVector.clear();
                    neighbors.clear();
                    if (i!=0) {
                        size_t r = std::min(_ranks[i-1],std::min(maxDim1, maxDim2));
                        constructVector.emplace_back(r);
                        neighbors.emplace_back(i-1, N+(i>=2?1:0), r, false);
                    }
                    for (size_t j=0; j< N; ++j) { 
                        size_t r = _dimensions[i+j*numNodes];
                        constructVector.emplace_back(r); 
                        neighbors.emplace_back(-1, i+j*numNodes, r, true);
                    }
                    if (i!=_dimensions.size()/N-1) {
                        size_t r = std::min(_ranks[i],std::min(maxDim1/_dimensions[i], maxDim2*_dimensions[i]));
                        constructVector.emplace_back(r);
                        neighbors.emplace_back(i+1, 0, r, false);
                    }
                    
                    result.nodes.emplace_back(
                        std::unique_ptr<Tensor>(new FullTensor(FullTensor::construct_random(constructVector, _rnd, _dist))), 
                        std::move(neighbors)
                    );
                    
                    maxDim1 /= _dimensions[i];
                    maxDim2 *= _dimensions[i];
                }
            }
            REQUIRE(result.nodes.size() == _dimensions.size()/N,"ie");
            result.cannonicalize_right();
            REQUIRE(result.is_valid_tt(), "ie");
            return result;
        }
        
		/// Random constructs a TTNetwork with the given dimensions and ranks. The entries of the componend tensors are sampled independendly using the provided random generator and distribution.
        template<class generator, class distribution>
        static TTNetwork construct_random(const std::vector<size_t>& _dimensions, size_t _rank, generator& _rnd, distribution& _dist) {
            const size_t N = isOperator?2:1;
            return construct_random(_dimensions, std::vector<size_t>(_dimensions.size()/N-1, _rank), _rnd, _dist);
        }
        
        /// Construct a TTOperator with the given dimensions representing the identity. (Only applicable for TTOperators, i.e. not for TTtensors).
        static TTNetwork construct_identity(const std::vector<size_t>& _dimensions);
        
        TTNetwork& operator=(const TTNetwork& _other) = default;
        
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Internal helper functions - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    protected:
        void round_train(const std::vector<size_t>& _maxRanks, const double _eps);

		static void construct_train_from_full(TensorNetwork& _out, const FullTensor& _A, const double _eps);
        
        static void contract_stack(const IndexedTensorWritable<TensorNetwork> &_me);
            
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Miscellaneous - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    public:
		/// tests whether the network resembles that of a TTTensor and checks consistency with the underlying tensor objects
        /// @note will not check for orthogonality
        bool is_valid_tt() const;
		
        static TTNetwork dyadic_product(const TTNetwork &_lhs, const TTNetwork &_rhs);
        
        static TTNetwork dyadic_product(const std::vector<std::reference_wrapper<TTNetwork>> &_tensors);
        
		/// returns a reference to the component tensor at position @a _idx in [0..numComponents]
		/// @note the first and last component tensor have virtual indices of dimension 1 added so that all components are of the same degree
		const Tensor &get_component(size_t _idx) const;
		
		/// sets a single component tensor at position @a _idx to be equal to @a _T
		/// updates meta-data but might leave the TT tensor in an illegal state. it is the callers responsibility to update the component tensors consistently
		void set_component(size_t _idx, const Tensor &_T);
		void set_component(size_t _idx, std::unique_ptr<Tensor> &&_T);
		
        /// Splits the TTNetwork into two parts by removing the node at @a _position.
        std::pair<TensorNetwork, TensorNetwork> chop(const size_t _position) const;
        
        void round(value_t _eps);

        void round(size_t _maxRank);

        void round(const std::vector<size_t> &_maxRanks);
        
        void round(int _maxRank);

        std::vector<size_t> ranks() const;
        
        size_t rank(size_t _i) const;
        
        size_t datasize() const;
        
        /// moves core to the left
        void cannonicalize_left();
        
        /// moves core to the right
        void cannonicalize_right();
            
        /// swaps all external indices to create the transposed operator
        template<bool B = isOperator, typename std::enable_if<B, int>::type = 0>
        void transpose() {
            Index i,r,l,j;
            for (size_t n=0; n<nodes.size(); ++n) {
                size_t ld = n==0?0:1;
                size_t rd = n==nodes.size()-1?0:1;
                REQUIRE(nodes[n].neighbors.size() == ld + rd + 2, "ie " << ld << " " << rd << " " << nodes[n].neighbors.size() << " n= " << n);
                (*nodes[n].tensorObject)(l^ld,j,i,r^rd) = (*nodes[n].tensorObject)(l^ld,i,j,r^rd);
                TensorNode::Link &l1 = nodes[n].neighbors[ld];
                TensorNode::Link &l2 = nodes[n].neighbors[ld+1];
                REQUIRE(l1.external, "ie "<< l1.other << " " << ld << " " << n << " " << nodes[n].neighbors.size() << " " << nodes.size());
                REQUIRE(l2.external, "ie "<< l2.other);
                std::swap(externalLinks[l1.indexPosition].dimension, externalLinks[l2.indexPosition].dimension);
                std::swap(l1.dimension, l2.dimension);
            }
        }
        
        
        virtual TensorNetwork* get_copy() const override;
        
        virtual value_t frob_norm() const override;
        
        virtual bool is_in_expected_format() const override;
        
        /*- - - - - - - - - - - - - - - - - - - - - - - - - -  Basic arithmetics - - - - - - - - - - - - - - - - - - - - - - - - - - */
        TTNetwork& operator+=(const TTNetwork& _other);
        TTNetwork  operator+(const TTNetwork& _other) const;
        
        TTNetwork& operator-=(const TTNetwork& _other);
        TTNetwork  operator-(const TTNetwork& _other) const;
        
        TTNetwork& operator*=(const value_t _prod);
        TTNetwork  operator*(const value_t _prod) const;
        
        TTNetwork& operator/=(const value_t _div);
        TTNetwork  operator/(const value_t _div) const;
        
        
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Operator specializations - - - - - - - - - - - - - - - - - - - - - - - - - - */
        virtual bool specialized_contraction(IndexedTensorWritable<TensorNetwork> &_out, const IndexedTensorReadOnly<TensorNetwork> &_me, const IndexedTensorReadOnly<TensorNetwork> &_other) const override;
        virtual bool specialized_sum(IndexedTensorWritable<TensorNetwork> &_out, const IndexedTensorReadOnly<TensorNetwork> &_me, const IndexedTensorReadOnly<TensorNetwork> &_other) const override;
        virtual void specialized_evaluation(const IndexedTensorWritable<TensorNetwork> &_me, const IndexedTensorReadOnly<TensorNetwork> &_other) override;
        
    };

    template<bool isOperator>
    static _inline_ TTNetwork<isOperator> operator*(const value_t _lhs, const TTNetwork<isOperator>& _rhs) { return _rhs*_lhs; }

    /// Returns the frobenius norm of the given tensor
    template<bool isOperator>
    static _inline_ value_t frob_norm(const TTNetwork<isOperator>& _tensor) { return _tensor.frob_norm(); }


    typedef TTNetwork<false> TTTensor;
    typedef TTNetwork<true> TTOperator;

    namespace internal {

        template<bool isOperator>
        /// Internal class used to represent stacks consiting of (possibly multiply) applications of TTOperators to either a TTTensor or TTOperator.
        class TTStack : public TTNetwork<isOperator> {
        public:
            /*- - - - - - - - - - - - - - - - - - - - - - - - - - Operator specializations - - - - - - - - - - - - - - - - - - - - - - - - - - */
            virtual void specialized_evaluation(const IndexedTensorWritable<TensorNetwork> &_me _unused_ , const IndexedTensorReadOnly<TensorNetwork> &_other _unused_) override {
                LOG(fatal, "TTStack not supported as a storing type");
            }
            
            virtual TensorNetwork* get_copy() const override {
                LOG(fatal, "forbidden");
                return nullptr;
            }
            
            virtual value_t frob_norm() const override {
                Index i;
                TTNetwork<isOperator> tmp(this->degree());
                tmp(i&0) = (*this)(i&0);
                return tmp.frob_norm();
            }
        };
    }
}
