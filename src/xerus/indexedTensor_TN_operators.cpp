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

#include <xerus/indexedTensor_TN_operators.h>
#include <xerus/indexedTensor_tensor_operators.h>
#include <xerus/index.h>
#include <xerus/fullTensor.h>
#include <xerus/tensorNetwork.h>
#include <xerus/misc/test.h>
#include <xerus/misc/missingFunctions.h>

namespace xerus {
    
    template<> 
    void IndexedTensorWritable<Tensor>::operator=(const IndexedTensorReadOnly<TensorNetwork> &_rhs) {
        std::vector<Index> rightIndices = _rhs.get_assigned_indices();
//         LOG(LetsTrySomething, _rhs.tensorObjectReadOnly->factor);
		TensorNetwork cpy(*_rhs.tensorObjectReadOnly);
		TensorNetwork::trace_out_double_indices(rightIndices, cpy(rightIndices));
        
        std::set<size_t> all;
        for (size_t i=0; i < cpy.nodes.size(); ++i) {
            if (!cpy.nodes[i].erased) all.insert(i);
        }
        
        if (all.empty()) {
			tensorObject->reset({}, DONT_SET_ZERO());
            tensorObject->operator[]({}) = cpy.factor;
		} else {
			size_t res = cpy.contract(all);
			
			std::vector<Index> externalOrder;
			for(size_t i = 0; i < cpy.nodes[res].neighbors.size(); ++i) { externalOrder.emplace_back(); }
			
			std::vector<Index> internalOrder;
			for(const TensorNode::Link& link: cpy.nodes[res].neighbors) {
				REQUIRE(link.external, "Internal Error " << link.other << " " << link.indexPosition);
				internalOrder.emplace_back(externalOrder[link.indexPosition]);
			}
			
			const std::vector<Index> myIndices = get_assigned_indices(get_eval_degree(rightIndices));
			std::vector<Index> outOrder;
			for (const Index &idx : myIndices) {
				REQUIRE(misc::contains(rightIndices, idx), "Every index on the LHS must appear somewhere on the RHS");
				size_t spanSum = 0;
				for (size_t j = 0; rightIndices[j] != idx; ++j) {
					spanSum += rightIndices[j].span;
				}
				
				for (size_t i=0; i < idx.span; ++i) {
					outOrder.push_back(externalOrder[spanSum+i]);
				}
			}
		
			(*tensorObject)(outOrder) = (*cpy.nodes[res].tensorObject)(internalOrder);
		}
    }

    
    template<> void IndexedTensorWritable<TensorNetwork>::operator=(const IndexedTensorReadOnly<Tensor>& _rhs) {
		tensorObject->specialized_evaluation(*this, static_cast<IndexedTensorMoveable<TensorNetwork>>(_rhs)); // TODO change this to not cast
    }
    
    
    /// shuffles the external links of _lhs according to the indices of the indexedTensors
    /// lhs contains a copy of rhs, thus we have to swap the rhs.indices to resemble those of the lhs
    void TensorNetwork::shuffle_indices(std::vector<Index> &_currentIndices, const IndexedTensorWritable<TensorNetwork> &_lhs) {
        const std::vector<Index> lhsIndices = _lhs.get_assigned_indices();
        // writeable copy of the left indices
        size_t passedDegree1=0;
        for (size_t i=0; i<_currentIndices.size(); passedDegree1+=_currentIndices[i].span, ++i) {
            if (_currentIndices[i]!= lhsIndices[i]) {
				// find correct index
				size_t j=i+1;
				size_t passedDegree2=passedDegree1+_currentIndices[i].span;
				for (; j<_currentIndices.size(); passedDegree2+=_currentIndices[j].span, ++j) {
					if (_currentIndices[j] == lhsIndices[i]) break;
				}
				if (j < _currentIndices.size()) {
					std::swap(_currentIndices[i],_currentIndices[j]);
					
					for (size_t n=0; n<_currentIndices[i].span; ++n) {
						_lhs.tensorObject->swap_external_links(passedDegree1+n,passedDegree2+n);
					}
				} else {
					// expand...
					LOG(fatal, "TN expansion marked as >won't fix<");
				}
			}
            REQUIRE(_currentIndices[i].span == lhsIndices[i].span, "Index span mismatch");
        }
    }

    template<>
    void IndexedTensorWritable<TensorNetwork>::operator=(const IndexedTensorReadOnly<TensorNetwork>& _rhs) {
		tensorObject->specialized_evaluation(*this, _rhs);
    }
    
    
    
    IndexedTensorMoveable<TensorNetwork> operator*(value_t _factor, const IndexedTensorReadOnly<TensorNetwork>  &  _rhs) {
        TensorNetwork *res = _rhs.tensorObjectReadOnly->get_copy();
        res->factor *= _factor;
        IndexedTensorMoveable<TensorNetwork> result(res, _rhs.indices);
        return result;
    }
    
    IndexedTensorMoveable<TensorNetwork> operator*(value_t _factor, IndexedTensorMoveable<TensorNetwork> &&  _rhs) {
        IndexedTensorMoveable<TensorNetwork> result(std::move(_rhs));
        result.tensorObject->factor *= _factor;
        return result;
    }
    
    



    IndexedTensorMoveable<TensorNetwork> operator*(const IndexedTensorReadOnly<TensorNetwork> &  _lhs, const IndexedTensorReadOnly<TensorNetwork> &  _rhs) {
        IndexedTensorMoveable<TensorNetwork> result;
        if(!_lhs.tensorObjectReadOnly->specialized_contraction(result, _lhs, _rhs) && !_rhs.tensorObjectReadOnly->specialized_contraction(result, _rhs, _lhs)) {
            result.tensorObject = new TensorNetwork(*_lhs.tensorObjectReadOnly);
            result.tensorObjectReadOnly = result.tensorObject;
            result.indices = _lhs.get_assigned_indices();
            result.deleteTensorObject = true;
            TensorNetwork::add_network_to_network(result, _rhs);
        }
        return result;
    }


    IndexedTensorMoveable<TensorNetwork> operator*(IndexedTensorMoveable<TensorNetwork> &&  _lhs, const IndexedTensorReadOnly<TensorNetwork>  &  _rhs) {
        IndexedTensorMoveable<TensorNetwork> result;
        if(!_lhs.tensorObjectReadOnly->specialized_contraction(result, std::move(_lhs), _rhs) && !_rhs.tensorObjectReadOnly->specialized_contraction(result, _rhs, std::move(_lhs))) {
            result.tensorObject = _lhs.tensorObject;
            result.tensorObjectReadOnly = _lhs.tensorObjectReadOnly;
            result.indices = _lhs.get_assigned_indices();
            result.deleteTensorObject = true;
            _lhs.deleteTensorObject = false;
            TensorNetwork::add_network_to_network(result, _rhs);
        }
        return result;
    }

    
    IndexedTensorMoveable<TensorNetwork> operator+(const IndexedTensorReadOnly<TensorNetwork>  &  _lhs, const IndexedTensorReadOnly<TensorNetwork>  &  _rhs) {
        IndexedTensorMoveable<TensorNetwork> result;
        if(!_lhs.tensorObjectReadOnly->specialized_sum(result, _lhs, _rhs) && !_rhs.tensorObjectReadOnly->specialized_sum(result, _rhs, _lhs)) {
            LOG(warning, "Using FullTensor fallback for TensorNetwork sum!");
            std::unique_ptr<FullTensor> tmpResult(new FullTensor(_lhs.degree())); //TODO sparse
            (*tmpResult)(_lhs.indices) = IndexedTensorMoveable<Tensor>(_lhs) + IndexedTensorMoveable<Tensor>(_rhs);
            result.assign(IndexedTensorMoveable<TensorNetwork>(new TensorNetwork(std::move(tmpResult)), _lhs.indices));
        }
        return result;
    }
    
    
    IndexedTensorMoveable<TensorNetwork> operator-(const IndexedTensorReadOnly<TensorNetwork>  & _lhs, const IndexedTensorReadOnly<TensorNetwork>  &  _rhs) {
        return _lhs+(-1*_rhs);
    }
    
    IndexedTensorMoveable<TensorNetwork> operator-(      IndexedTensorMoveable<TensorNetwork> && _lhs, const IndexedTensorReadOnly<TensorNetwork> &  _rhs) {
		_lhs.tensorObject->factor *= -1;
		return (-1)*operator+(std::move(_lhs), _rhs);
	}
    IndexedTensorMoveable<TensorNetwork> operator-(const IndexedTensorReadOnly<TensorNetwork> &  _lhs,       IndexedTensorMoveable<TensorNetwork> && _rhs) {
		_rhs.tensorObject->factor*=-1;
		return operator+(std::move(_rhs), _lhs);
	} 
    IndexedTensorMoveable<TensorNetwork> operator-(      IndexedTensorMoveable<TensorNetwork> && _lhs,       IndexedTensorMoveable<TensorNetwork> && _rhs) {
		_rhs.tensorObject->factor*=-1;
		return operator+(std::move(_rhs), std::move(_lhs));
	}
}