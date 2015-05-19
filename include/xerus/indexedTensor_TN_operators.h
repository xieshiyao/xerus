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

#include "indexedTensorMoveable.h"

namespace xerus {
    // Defined in the IndexedTensor classes: 
    template<> void IndexedTensorWritable<Tensor>::operator=(const IndexedTensorReadOnly<TensorNetwork>&  _rhs);
    
    template<> void IndexedTensorWritable<TensorNetwork>::operator=(const IndexedTensorReadOnly<Tensor>&  _rhs);
    
    template<> void IndexedTensorWritable<TensorNetwork>::operator=(const IndexedTensorReadOnly<TensorNetwork>&  _rhs);
    
    
    
    IndexedTensorMoveable<TensorNetwork> operator*(value_t _factor, const IndexedTensorReadOnly<TensorNetwork>& _network);
	
	IndexedTensorMoveable<TensorNetwork> operator*(value_t _factor, IndexedTensorMoveable<TensorNetwork> &&_network);

    _inline_ IndexedTensorMoveable<TensorNetwork> operator*(const IndexedTensorReadOnly<TensorNetwork> & _network, value_t _factor) {
        return _factor*_network;
    }
    
    _inline_ IndexedTensorMoveable<TensorNetwork> operator*(IndexedTensorMoveable<TensorNetwork> &&_network, value_t _factor) {
        return _factor*std::move(_network);
    }

    
    
    IndexedTensorMoveable<TensorNetwork> operator*(const IndexedTensorReadOnly<TensorNetwork> &  _lhs, const IndexedTensorReadOnly<TensorNetwork> &  _rhs); // creates a new network
    
    IndexedTensorMoveable<TensorNetwork> operator*(      IndexedTensorMoveable<TensorNetwork> && _lhs, const IndexedTensorReadOnly<TensorNetwork> &  _rhs); // copy right into left network
    
    _inline_ IndexedTensorMoveable<TensorNetwork> operator*(const IndexedTensorReadOnly<TensorNetwork> &  _lhs,       IndexedTensorMoveable<TensorNetwork> && _rhs) {
		return operator*(std::move(_rhs), _lhs);
	}
	
    _inline_ IndexedTensorMoveable<TensorNetwork> operator*(      IndexedTensorMoveable<TensorNetwork> && _lhs,       IndexedTensorMoveable<TensorNetwork> && _rhs) {
		return operator*(std::move(_lhs), _rhs);
	}

	
	
    IndexedTensorMoveable<TensorNetwork> operator+(const IndexedTensorReadOnly<TensorNetwork> &  _lhs, const IndexedTensorReadOnly<TensorNetwork> &  _rhs); // creates a new network
    
    _inline_ IndexedTensorMoveable<TensorNetwork> operator+(      IndexedTensorMoveable<TensorNetwork> && _lhs, const IndexedTensorReadOnly<TensorNetwork> &  _rhs) { // TODO own function
		return operator+(_lhs, _rhs);
	}
	
    _inline_ IndexedTensorMoveable<TensorNetwork> operator+(const IndexedTensorReadOnly<TensorNetwork> &  _lhs,       IndexedTensorMoveable<TensorNetwork> && _rhs) {
		return operator+(std::move(_rhs), _lhs);
	}
	
    _inline_ IndexedTensorMoveable<TensorNetwork> operator+(      IndexedTensorMoveable<TensorNetwork> && _lhs,       IndexedTensorMoveable<TensorNetwork> && _rhs) {
		return operator+(std::move(_lhs), _rhs);
	}
	
	
	
	IndexedTensorMoveable<TensorNetwork> operator-(const IndexedTensorReadOnly<TensorNetwork> &  _lhs, const IndexedTensorReadOnly<TensorNetwork> &  _rhs); // creates a new network
    IndexedTensorMoveable<TensorNetwork> operator-(      IndexedTensorMoveable<TensorNetwork> && _lhs, const IndexedTensorReadOnly<TensorNetwork> &  _rhs);
    IndexedTensorMoveable<TensorNetwork> operator-(const IndexedTensorReadOnly<TensorNetwork> &  _lhs,       IndexedTensorMoveable<TensorNetwork> && _rhs);
    IndexedTensorMoveable<TensorNetwork> operator-(      IndexedTensorMoveable<TensorNetwork> && _lhs,       IndexedTensorMoveable<TensorNetwork> && _rhs); 
}
