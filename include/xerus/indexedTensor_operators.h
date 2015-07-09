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
 * @brief Header file for the basic operators of general indexedTensors.
 */

#pragma once

#include "indexedTensorMoveable.h"

namespace xerus {
	template<class tensor_type>
    static inline IndexedTensorMoveable<tensor_type> operator*(const value_t _factor, const IndexedTensorReadOnly<tensor_type>& _tensor) {
		IndexedTensorMoveable<tensor_type> result(_tensor);
		*result.tensorObject *= _factor;
		return result;
	}
	
	template<class tensor_type>
	static inline IndexedTensorMoveable<tensor_type> operator*(const IndexedTensorReadOnly<tensor_type>& _tensor, const value_t _factor) {
		return operator*(_factor, _tensor);
	}
	
	
    template<class tensor_type>
	static inline IndexedTensorMoveable<tensor_type> operator*(const value_t _factor, IndexedTensorMoveable<tensor_type>&& _tensor) {
		IndexedTensorMoveable<tensor_type> result(std::move(_tensor));
		*result.tensorObject *= _factor;
		return result;
	}
	
	template<class tensor_type>
    static inline IndexedTensorMoveable<tensor_type> operator*(IndexedTensorMoveable<tensor_type>&& _tensor, const value_t _factor) {
		return operator*(_factor, std::move(_tensor));
	}
	
    template<class tensor_type>
	static inline IndexedTensorMoveable<tensor_type> operator/(const IndexedTensorReadOnly<tensor_type> & _tensor, const value_t _divisor) {
		IndexedTensorMoveable<tensor_type> result(_tensor);
		*result.tensorObject /= _divisor;
		return result;
	}
	
    template<class tensor_type>
    static inline IndexedTensorMoveable<tensor_type> operator/(IndexedTensorMoveable<tensor_type>&& _tensor, const value_t _divisor) {
		IndexedTensorMoveable<tensor_type> result(std::move(_tensor));
		*result.tensorObject /= _divisor;
		return result;
	}
	
	void operator+=(IndexedTensorWritable<Tensor> &  _lhs, const IndexedTensorReadOnly<Tensor> &  _rhs);
	void operator-=(IndexedTensorWritable<Tensor> &  _lhs, const IndexedTensorReadOnly<Tensor> &  _rhs);
	
	IndexedTensorMoveable<Tensor> operator+(const IndexedTensorReadOnly<Tensor> &  _lhs, const IndexedTensorReadOnly<Tensor> &  _rhs);
	IndexedTensorMoveable<Tensor> operator+(      IndexedTensorMoveable<Tensor> && _lhs, const IndexedTensorReadOnly<Tensor> &  _rhs);
	IndexedTensorMoveable<Tensor> operator+(const IndexedTensorReadOnly<Tensor> &  _lhs,       IndexedTensorMoveable<Tensor> && _rhs);
	IndexedTensorMoveable<Tensor> operator+(      IndexedTensorMoveable<Tensor> && _lhs,       IndexedTensorMoveable<Tensor> && _rhs);
	
	IndexedTensorMoveable<Tensor> operator-(const IndexedTensorReadOnly<Tensor> &  _lhs, const IndexedTensorReadOnly<Tensor> &  _rhs);
	IndexedTensorMoveable<Tensor> operator-(      IndexedTensorMoveable<Tensor> && _lhs, const IndexedTensorReadOnly<Tensor> &  _rhs);
	IndexedTensorMoveable<Tensor> operator-(const IndexedTensorReadOnly<Tensor> &  _lhs,       IndexedTensorMoveable<Tensor> && _rhs);
	IndexedTensorMoveable<Tensor> operator-(      IndexedTensorMoveable<Tensor> && _lhs,       IndexedTensorMoveable<Tensor> && _rhs);
	
	
	IndexedTensorMoveable<Tensor> operator+(const IndexedTensorReadOnly<Tensor> &  _lhs, const IndexedTensorReadOnly<TensorNetwork> &  _rhs);
	IndexedTensorMoveable<Tensor> operator+(const IndexedTensorReadOnly<TensorNetwork> &  _lhs, const IndexedTensorReadOnly<Tensor> &  _rhs);
	IndexedTensorMoveable<Tensor> operator+(      IndexedTensorMoveable<Tensor> && _lhs, const IndexedTensorReadOnly<TensorNetwork> &  _rhs);
	IndexedTensorMoveable<Tensor> operator+(const IndexedTensorReadOnly<TensorNetwork> &  _lhs,       IndexedTensorMoveable<Tensor> && _rhs);
	
	IndexedTensorMoveable<Tensor> operator-(const IndexedTensorReadOnly<Tensor> &  _lhs, const IndexedTensorReadOnly<TensorNetwork> &  _rhs);
	IndexedTensorMoveable<Tensor> operator-(const IndexedTensorReadOnly<TensorNetwork> &  _lhs, const IndexedTensorReadOnly<Tensor> &  _rhs);
	IndexedTensorMoveable<Tensor> operator-(      IndexedTensorMoveable<Tensor> && _lhs, const IndexedTensorReadOnly<TensorNetwork> &  _rhs);
	IndexedTensorMoveable<Tensor> operator-(const IndexedTensorReadOnly<TensorNetwork> &  _lhs,       IndexedTensorMoveable<Tensor> && _rhs);
	
	
	IndexedTensorMoveable<TensorNetwork> operator+(const IndexedTensorReadOnly<TensorNetwork> &  _lhs, const IndexedTensorReadOnly<TensorNetwork> &  _rhs);
	IndexedTensorMoveable<TensorNetwork> operator-(const IndexedTensorReadOnly<TensorNetwork> &  _lhs, const IndexedTensorReadOnly<TensorNetwork> &  _rhs);
	
	void operator+=(IndexedTensorWritable<TensorNetwork> &  _lhs, const IndexedTensorReadOnly<TensorNetwork> &  _rhs);
	void operator-=(IndexedTensorWritable<TensorNetwork> &  _lhs, const IndexedTensorReadOnly<TensorNetwork> &  _rhs);
}
