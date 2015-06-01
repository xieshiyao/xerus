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
	public:
		static constexpr const size_t N=isOperator?2:1;
		
		bool cannonicalized;
		size_t corePosition;		
		
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
            REQUIRE(_ranks.size() == _dimensions.size()/N-1,"Non-matching amount of ranks given to TTNetwork::construct_random");
            #ifndef DISABLE_RUNTIME_CHECKS_
                for (const size_t d : _dimensions) {
                    REQUIRE(d > 0, "Trying to construct random TTTensor with dimension 0 is illegal.");
                }
                for (const size_t d : _ranks) {
                    REQUIRE(d > 0, "Trying to construct random TTTensor with rank 0 is illegal.");
                }
            #endif
            
            TTNetwork result(_dimensions.size());
			const size_t numComponents = _dimensions.size()/N;
			size_t maxDim1 = 1;
			size_t maxDim2 = misc::product(_dimensions);
			
			for(size_t i = 0; i < numComponents; ++i) {
				size_t oldMaxDim = std::min(maxDim1, maxDim2);
				maxDim1 *= _dimensions[i];
				maxDim2 /= _dimensions[i];
				size_t maxDim = std::min(maxDim1, maxDim2);
				result.set_component(i, FullTensor::construct_random({i==0?1:std::min(oldMaxDim, _ranks[i-1]), _dimensions[i], i==numComponents-1?1:std::min(maxDim, _ranks[i])}, _rnd, _dist));
			}
            result.cannonicalize_left();
            REQUIRE(result.is_valid_tt(), "Internal Error.");
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
        
		/// @brief moves the core to @a _position
		/// all components left of @a _position will be left-orthogonal, those to the right will be right-orthogonal
		void move_core(size_t _position);
		
        /// moves core to the left
        void cannonicalize_left();
        
        /// moves core to the right
        void cannonicalize_right();
            
        /// swaps all external indices to create the transposed operator
        template<bool B = isOperator, typename std::enable_if<B, int>::type = 0>
        void transpose() {
            Index i,r,l,j;
            for (size_t n=0; n < degree(); ++n) {
				std::unique_ptr<Tensor> At(nodes[n].tensorObject->construct_new());
                (*At)(l,j,i,r) = get_component(n)(l,i,j,r);
				set_component(n, std::move(At));
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
                LOG(fatal, "Forbidden");
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
