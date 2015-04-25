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
    IndexedTensorMoveable<TensorNetwork>::IndexedTensorMoveable(const IndexedTensorReadOnly<Tensor> &  _other) : 
        IndexedTensorWritable<TensorNetwork>(new TensorNetwork(*_other.tensorObjectReadOnly), _other.indices, true) { }
    
    template<>
    IndexedTensorMoveable<TensorNetwork>::IndexedTensorMoveable(      IndexedTensorReadOnly<Tensor> && _other) : 
        IndexedTensorWritable<TensorNetwork>(new TensorNetwork(*_other.tensorObjectReadOnly), std::move(_other.indices), true) { }
    
    // Allow casts from tensor Network to FullTensor
    IndexedTensorMoveable<Tensor>::IndexedTensorMoveable(const IndexedTensorReadOnly<TensorNetwork> &  _other ) : 
        IndexedTensorWritable<Tensor>(_other.tensorObjectReadOnly->fully_contracted_tensor().release(), _other.indices, true) { }

    IndexedTensorMoveable<Tensor>::IndexedTensorMoveable(      IndexedTensorReadOnly<TensorNetwork> && _other ) : 
        IndexedTensorWritable<Tensor>(_other.tensorObjectReadOnly->fully_contracted_tensor().release(), std::move(_other.indices), true) { }
}
