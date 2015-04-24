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

#include "xerus.h"

namespace xerus {

    template<>
    void IndexedTensorWritable<Tensor>::operator=(const IndexedTensorReadOnly<Tensor>&  _rhs)  {
        if(!_rhs.uses_tensor(tensorObject)) {
            // If LHS and RHS object don't coincide we can directly evaluate
            this->tensorObject->reset(_rhs.get_evaluated_dimensions(indices), DONT_SET_ZERO());
            evaluate(*this, _rhs);
        } else {
            // If the tensors in fact coincide we have to use a tmp object
            IndexedTensorMoveable<Tensor> tmpTensor(_rhs);
            this->tensorObject->reset(_rhs.get_evaluated_dimensions(indices), DONT_SET_ZERO());
            evaluate(*this, tmpTensor);
        }
    }

    namespace internal {
        template<bool PLUS>
        _inline_ IndexedTensorMoveable<Tensor> plus_minus(const IndexedTensorReadOnly<Tensor> & _lhs, const IndexedTensorReadOnly<Tensor> & _rhs) {
            const std::vector<Index> lhsIndices = _lhs.get_assigned_indices();
            const std::vector<Index> rhsIndices = _rhs.get_assigned_indices();
            
            // Get open rhs indices that are not contained in lhs
            std::vector<Index> lhsNotRhs;
            std::vector<Index> lhsOpen;
            size_t lhsOpenDim = 0;
            for(size_t i = 0; i < lhsIndices.size(); ++i) {
                if(_lhs.is_open(lhsIndices[i])) {
                    lhsOpen.emplace_back(lhsIndices[i]);
                    lhsOpenDim += lhsIndices[i].span;
                    if(!contains(rhsIndices, lhsIndices[i])) {
                        lhsNotRhs.emplace_back(lhsIndices[i]);
                    }
                }
            }
            
            // Get open lhs indices that are not contained in rhs
            std::vector<Index> rhsNotLhs;
            std::vector<Index> rhsOpen;
            size_t rhsOpenDim = 0;
            for(size_t i = 0; i < rhsIndices.size(); ++i) {
                if(_rhs.is_open(rhsIndices[i])) {
                    rhsOpen.emplace_back(rhsIndices[i]);
                    rhsOpenDim += rhsIndices[i].span;
                    if(!contains(lhsIndices, rhsIndices[i])) {
                        rhsNotLhs.emplace_back(rhsIndices[i]);
                    }
                }
            }
             
            if(rhsNotLhs.empty() && lhsNotRhs.empty()) {
                IndexedTensorMoveable<Tensor> tmp1(_lhs.tensorObjectReadOnly->construct_new(std::vector<size_t>(lhsOpenDim, 1)), lhsOpen); //TODO get rid of vector
                static_cast<IndexedTensorWritable<Tensor>&>(tmp1) = _lhs;
                IndexedTensorMoveable<Tensor> tmp2(_rhs.tensorObjectReadOnly->construct_new(std::vector<size_t>(lhsOpenDim, 1)), lhsOpen); //TODO get rid of vector
                static_cast<IndexedTensorWritable<Tensor>&>(tmp2) = _rhs;
                
                if(!tmp1.tensorObjectReadOnly->is_sparse() && !tmp2.tensorObjectReadOnly->is_sparse()) { // Both Full
                    if(PLUS) {
                        static_cast<FullTensor&>(*tmp1.tensorObject) += static_cast<FullTensor&>(*tmp2.tensorObject);
                    } else {
                        static_cast<FullTensor&>(*tmp1.tensorObject) -= static_cast<FullTensor&>(*tmp2.tensorObject);
                    }
                    return tmp1;
                } else if(tmp1.tensorObjectReadOnly->is_sparse() && tmp2.tensorObjectReadOnly->is_sparse()) { // Both Sparse
                    if(PLUS) {
                        static_cast<SparseTensor&>(*tmp1.tensorObject) += static_cast<SparseTensor&>(*tmp2.tensorObject);
                    } else {
                        static_cast<SparseTensor&>(*tmp1.tensorObject) -= static_cast<SparseTensor&>(*tmp2.tensorObject);
                    }
                    return tmp1;
                } else if(!tmp1.tensorObjectReadOnly->is_sparse() && tmp2.tensorObjectReadOnly->is_sparse()) { // Full + Sparse
                    if(PLUS) {
                        static_cast<FullTensor&>(*tmp1.tensorObject) += static_cast<SparseTensor&>(*tmp2.tensorObject);
                    } else {
                        static_cast<FullTensor&>(*tmp1.tensorObject) -= static_cast<SparseTensor&>(*tmp2.tensorObject);
                    }
                    return tmp1;
                } else { // Sparse + Full (switch lhs and rhs)
                    if(PLUS) {
                        static_cast<FullTensor&>(*tmp2.tensorObject) += static_cast<SparseTensor&>(*tmp1.tensorObject);
                    } else {
                        static_cast<FullTensor&>(*tmp2.tensorObject) -= static_cast<SparseTensor&>(*tmp1.tensorObject);
                        tmp2.tensorObject->factor *= -1;
                    }
                    return tmp2;
                }
            } else {
                LOG(fatal, "Sum/Diff with non common indices not yet implemented.");
                return IndexedTensorMoveable<Tensor>(_lhs.tensorObjectReadOnly->construct_new({}), lhsOpen);
            }
        }
    }
    
    IndexedTensorMoveable<Tensor> operator+(const IndexedTensorReadOnly<Tensor> & _lhs, const IndexedTensorReadOnly<Tensor> & _rhs) {
        return internal::plus_minus<true>(_lhs, _rhs);
    }
        
    IndexedTensorMoveable<Tensor> operator-(const IndexedTensorReadOnly<Tensor> & _lhs, const IndexedTensorReadOnly<Tensor> & _rhs) {
        return internal::plus_minus<false>(_lhs, _rhs);
    }
    
    IndexedTensorMoveable<Tensor> operator*(const value_t _lhs, const IndexedTensorReadOnly<Tensor>& _rhs) {
        IndexedTensorMoveable<Tensor> result(_rhs);
        result.tensorObject->factor *= _lhs;
        return result;
    }

    IndexedTensorMoveable<Tensor> operator*(const IndexedTensorReadOnly<Tensor> & _lhs, const value_t _rhs) {
        return _rhs * _lhs;
    }

    IndexedTensorMoveable<Tensor> operator/(const IndexedTensorReadOnly<Tensor> & _lhs, const value_t _rhs) {
        return (1/_rhs)*_lhs;
    }
}
