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

#include "indexedTensor.h"

namespace xerus {
    /// IndexedTensor variant that signifies, that the tensorObject is moveable. meant for internal use only (strange things happen when storing these)
    template<class tensor_type>
    class IndexedTensorMoveable final : public IndexedTensorWritable<tensor_type> {
    public:
        /// Creates an empty indexed Tensor, should only be used internally
        IndexedTensorMoveable();
        
        /// There is no usefull copy constructor, because the handling of the tensorObject is unclear
        IndexedTensorMoveable(const IndexedTensorMoveable &_other ) = delete;
        
        /// Move constructor
        IndexedTensorMoveable(IndexedTensorMoveable &&_other );
        
        /// Constructs an IndexedTensorMoveable with the given tensor and indices and if ordered to do so takes ownership of the tensorObject
        IndexedTensorMoveable(tensor_type* const _tensorObject, const std::vector<Index>& _indices);
        
        /// Constructs an IndexedTensorMoveable with the given tensor and indices and if ordered to do so takes ownership of the tensorObject
        IndexedTensorMoveable(tensor_type* const _tensorObject, std::vector<Index>&& _indices);
        
        /// Allow conversions from indexed TensorNetworks to indexed Tensors
        explicit IndexedTensorMoveable(const IndexedTensorReadOnly<TensorNetwork> &  _other );
        explicit IndexedTensorMoveable(      IndexedTensorReadOnly<TensorNetwork> && _other );
        
        //TODO do we want implicit conversions from IndexedTensorReadOnly<Tensor> to IndexedTensorMoveable<Tensor> ?
        /// Allow implicit conversions from indexed Tensors to indexed TensorNetworks
        IndexedTensorMoveable(const IndexedTensorReadOnly<Tensor> &  _other);
        IndexedTensorMoveable(      IndexedTensorReadOnly<Tensor> && _other);
    };
}
