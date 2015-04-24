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

#include "indexedTensorWritable.h"

namespace xerus {
	template<class tensor_type, tensor_type_restrictions>

	/**
	 * usercreated indices for tensor equations, e.g. @f$ A_{i,j} = B_{i,k,l} C_{k,l,j} @f$
	 */
    class IndexedTensor final : public IndexedTensorWritable<tensor_type> {
    public:     
        /// There is no usefull copy constructor, because the handling of the tensorObject is unclear
        IndexedTensor(const IndexedTensor &_other ) = delete;
        
        /// Move constructor
        IndexedTensor(IndexedTensor &&_other ) : IndexedTensorWritable<tensor_type>(std::move(_other)) { }
        
        /// Constructs an IndexedTensor with the given tensor and indices and if ordered to do so takes owership of the tensorObject
        ALLOW_MOVE(std::vector<Index>, T)
        IndexedTensor(tensor_type* const _tensorObject, T&& _indices, const bool _takeOwnership) :
            IndexedTensorWritable<tensor_type>(_tensorObject, std::forward<T>(_indices), _takeOwnership) {}
            
            
        // Assignment operators -- Used for tensor assignment WITH indices (i.e. in general the LHS and RHS indexTensors do NOT have the same indices)
        // NOTE: The following would be deleted due to move constructor
        _inline_ void operator=(const IndexedTensor<tensor_type>& _rhs) {
            operator=(static_cast<const IndexedTensorReadOnly<tensor_type> &>(_rhs));
        }
        
        /// Assignment operators
        void operator=(const IndexedTensorReadOnly<Tensor>&         _rhs) {
            static_cast<IndexedTensorWritable<tensor_type>&>(*this) = _rhs;
        }
        
        void operator=(const IndexedTensorReadOnly<TensorNetwork>&  _rhs) {
            static_cast<IndexedTensorWritable<tensor_type>&>(*this) = _rhs;
        }
    };
}
