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

namespace xerus {

template<bool isOperator>
class TTNetwork : public TensorNetwork {    
public:
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Constructors - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    explicit TTNetwork();
    
    explicit TTNetwork(const size_t _degree);
    
    explicit TTNetwork(const FullTensor& _full, const double _eps=1e-15) { //TODO no magic numbers
		construct_train_from_full(*this, _full, _eps);
	}
    
	explicit TTNetwork(const TTNetwork & _cpy);
    
    implicit TTNetwork(      TTNetwork&& _mov);
	
	explicit TTNetwork(const TensorNetwork &_cpy, double _eps=1e-14);
	
	explicit TTNetwork(TensorNetwork &&_mov, double _eps=1e-14);
	
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
			
            size_t maxDim1 = product(_dimensions);
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
					std::shared_ptr<Tensor>(new FullTensor(FullTensor::construct_random(constructVector, _rnd, _dist))), 
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
    
    template<class generator, class distribution>
	_inline_ static TTNetwork construct_random(const std::vector<size_t>& _dimensions, size_t _rank, generator& _rnd, distribution& _dist) {
		const size_t N = isOperator?2:1;
		return construct_random(_dimensions, std::vector<size_t>(_dimensions.size()/N-1, _rank), _rnd, _dist);
	}
    
    static TTNetwork construct_identity(const std::vector<size_t>& _dimensions);
    
    //TODO does this make sense?
    TTNetwork& operator=(const TTNetwork& _other) = default;

    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Internal helper functions - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    protected:
    static void construct_train_from_full(TensorNetwork& _out, const FullTensor& _A, const double _eps);
    
    static void round_train(TensorNetwork& _me, const std::vector<size_t>& _maxRanks, const double _eps);
    
    static void contract_stack(const IndexedTensorWritable<TensorNetwork> &_me);
        
    /// tests whether the network resembles that of a TTTensor and checks consistency with the udnerlying tensor objects
    /// @note will not check for orthogonality
    bool is_valid_tt() const;
    
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Miscellaneous - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
public:
    static TTNetwork dyadic_product(const TTNetwork &_lhs, const TTNetwork &_rhs);
	
	static TTNetwork dyadic_product(const std::vector<std::reference_wrapper<TTNetwork>> &_tensors);
	
	void round(value_t _eps) {
		cannonicalize_left();
		const size_t N = isOperator?2:1;
		round_train(*this, std::vector<size_t>(degree()/N-1, size_t(-1)), _eps);
	}

	void round(size_t _maxRank) {
		cannonicalize_left();
		const size_t N = isOperator?2:1;
		round_train(*this, std::vector<size_t>(degree()/N-1 ,_maxRank), 1e-15);
	}

	void round(const std::vector<size_t> &_maxRanks) {
		cannonicalize_left();
		round_train(*this, _maxRanks, 1e-15);
	}
    
    _inline_ void round(int _maxRank) {
        REQUIRE( _maxRank > 0, "MaxRank must be positive");
        round(size_t(_maxRank));
    }
	
	
	std::vector<size_t> ranks() const {
		std::vector<size_t> res;
		for (size_t n=0; n<nodes.size()-1; ++n) {
			res.push_back(nodes[n].neighbors.back().dimension);
		}
		return res;
	}
	
	size_t rank(size_t _i) const {
		REQUIRE(_i < nodes.size()-1, "requested illegal rank");
		return nodes[_i].neighbors.back().dimension;
	}
	
	size_t datasize() const {
		size_t result = 0;
		for (const TensorNode &n : nodes) {
			result += n.tensorObject->size;
		}
		return result;
	}

	
	/// moves core to the left (eg. to perform a round operation that then moves it to the right)
	void cannonicalize_left() {
		Index i,r,j;
		FullTensor core(2);
		for (size_t n=nodes.size()-1; n>0; --n) {
            REQUIRE(!nodes[n].erased, "ie n="<<n);
			Tensor &currTensor = *nodes[n].tensorObject;
			( core(j,r), currTensor(r,i&1) ) = RQ(currTensor(j,i&1));  //TODO we want a rank-detecting QR at this point?
			Tensor &nextTensor = *nodes[n-1].tensorObject;
			nextTensor(j&1,i) = nextTensor(j&1,r) * core(r,i);
			if (currTensor.dimensions[0] != nodes[n].neighbors.front().dimension) {
				nodes[n].neighbors.front().dimension = nodes[n-1].neighbors.back().dimension = currTensor.dimensions[0];
			}
		}
	}
	
	/// moves core to the right
	void cannonicalize_right() {
		Index i,r,j;
		FullTensor core(2);
		for (size_t n=0; n<nodes.size()-1; ++n) {
			Tensor &currTensor = *nodes[n].tensorObject;
			( currTensor(i&1,r), core(r,j) ) = QR(currTensor(i&1,j)); //TODO we want a rank-detecting QR at this point?
			Tensor &nextTensor = *nodes[n+1].tensorObject;
			nextTensor(i,j&1) = core(i,r) * nextTensor(r,j&1);
			if (nextTensor.dimensions[0] != nodes[n+1].neighbors.front().dimension) {
				nodes[n+1].neighbors.front().dimension = nodes[n].neighbors.back().dimension = nextTensor.dimensions[0];
			}
		}
	}
	    
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
	
	/*- - - - - - - - - - - - - - - - - - - - - - - - - -  Basic arithmetics - - - - - - - - - - - - - - - - - - - - - - - - - - */
	TTNetwork& operator+=(const TTNetwork& _other) {
		Index i;
		(*this)(i&0) = (*this)(i&0) + _other(i&0);
		return *this;
	}
	TTNetwork  operator+(const TTNetwork& _other) const {
		TTNetwork cpy(*this);
		cpy += _other;
		return cpy;
	}
	
	TTNetwork& operator-=(const TTNetwork& _other) {
		Index i;
		(*this)(i&0) = (*this)(i&0) - _other(i&0);
		return *this;
	}
	TTNetwork  operator-(const TTNetwork& _other) const {
		TTNetwork cpy(*this);
		cpy -= _other;
		return cpy;
	}
	
	TTNetwork& operator*=(const value_t _prod) {
		factor *= _prod;
		return *this;
	}
	TTNetwork  operator*(const value_t _prod) const {
		TTNetwork result(*this);
		result *= _prod;
		return result;
	}
	
	TTNetwork& operator/=(const value_t _div) {
		factor /= _div;
		return *this;
	}
	TTNetwork  operator/(const value_t _div) const {
		TTNetwork result(*this);
		result /= _div;
		return result;
	}
	
	
	/*- - - - - - - - - - - - - - - - - - - - - - - - - - Operator specializations - - - - - - - - - - - - - - - - - - - - - - - - - - */
	virtual bool specialized_contraction(IndexedTensorWritable<TensorNetwork> &_out, const IndexedTensorReadOnly<TensorNetwork> &_me, const IndexedTensorReadOnly<TensorNetwork> &_other) const override;
	virtual bool specialized_sum(IndexedTensorWritable<TensorNetwork> &_out, const IndexedTensorReadOnly<TensorNetwork> &_me, const IndexedTensorReadOnly<TensorNetwork> &_other) const override;
	virtual void specialized_evaluation(const IndexedTensorWritable<TensorNetwork> &_me, const IndexedTensorReadOnly<TensorNetwork> &_other) override;
	
	virtual TensorNetwork* get_copy() const override {
		return new TTNetwork(*this);
	}
	
	virtual value_t frob_norm() const override {
		REQUIRE(is_valid_tt(), "frob_norm of illegal TT");
		return nodes.back().tensorObject->frob_norm();
	}
	
	virtual bool is_in_expected_format() const override {
		return is_valid_tt();
	}
};

template<bool isOperator>
_inline_ TTNetwork<isOperator> operator*(const value_t _lhs, const TTNetwork<isOperator>& _rhs) { return _rhs*_lhs; }

/// Returns the frobenius norm of the given tensor
template<bool isOperator>
_inline_ value_t frob_norm(const TTNetwork<isOperator>& _tensor) { return _tensor.frob_norm(); }


typedef TTNetwork<false> TTTensor;
typedef TTNetwork<true> TTOperator;

namespace internal {

template<bool isOperator>
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


