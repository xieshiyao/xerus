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

/**
 * @file
 * @brief Header file for the basic operators of indexed tensor networks.
 */

#pragma once

#include "indexedTensorMoveable.h"
#include "indexedTensor_operators.h"

namespace xerus {
    IndexedTensorMoveable<TensorNetwork> operator*(const IndexedTensorReadOnly<TensorNetwork> &  _lhs, const IndexedTensorReadOnly<TensorNetwork> &  _rhs); // creates a new network
    
    IndexedTensorMoveable<TensorNetwork> operator*(      IndexedTensorMoveable<TensorNetwork> && _lhs, const IndexedTensorReadOnly<TensorNetwork> &  _rhs); // copy right into left network
    
    static _inline_ IndexedTensorMoveable<TensorNetwork> operator*(const IndexedTensorReadOnly<TensorNetwork> &  _lhs,       IndexedTensorMoveable<TensorNetwork> && _rhs) {
		return operator*(std::move(_rhs), _lhs);
	}
	
    static _inline_ IndexedTensorMoveable<TensorNetwork> operator*(      IndexedTensorMoveable<TensorNetwork> && _lhs,       IndexedTensorMoveable<TensorNetwork> && _rhs) {
		return operator*(std::move(_lhs), _rhs);
	}

}
