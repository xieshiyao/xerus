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
	template<class tensor_type>
	/**
	 * @brief Internal representation of an readable and writeable indexed Tensor or TensorNetwork.
	 * @details This class appears inplicitly by indexing any Tensor or TensorNetwork. It is not recommended to use
	 * it explicitly or to store variables of this type (unless you really know what you are doing).
	 */
    class IndexedTensor final : public IndexedTensorWritable<tensor_type> {
    public:     
        /// There is no usefull copy constructor, because the handling of the tensorObject is unclear
        IndexedTensor(const IndexedTensor &_other ) = delete;
        
        /// Move constructor
        IndexedTensor(IndexedTensor &&_other );
        
        /// Constructs an IndexedTensor with the given tensor and indices and if ordered to do so takes owership of the tensorObject
        IndexedTensor(tensor_type* const _tensorObject, const std::vector<Index>& _indices, const bool _takeOwnership);
        
        /// Constructs an IndexedTensor with the given tensor and indices and if ordered to do so takes owership of the tensorObject
        IndexedTensor(tensor_type* const _tensorObject, std::vector<Index>&& _indices, const bool _takeOwnership);
            
            
        // Assignment operators -- Used for tensor assignment WITH indices (i.e. in general the LHS and RHS indexTensors do NOT have the same indices)
        // NOTE: The following would be deleted due to move constructor
        void operator=(const IndexedTensor<tensor_type>& _rhs);
        
        /// Assignment operator
        void operator=(const IndexedTensorReadOnly<Tensor>&         _rhs);
        
        /// Assignment operator
        void operator=(const IndexedTensorReadOnly<TensorNetwork>&  _rhs);
    };
}
