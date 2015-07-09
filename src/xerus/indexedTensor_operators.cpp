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
 * @brief Implementation of the basic indexedTensor operators.
 */

#include <xerus/indexedTensor_TN_operators.h>
#include <xerus/indexedTensor_tensor_operators.h>
#include <xerus/index.h>
#include <xerus/fullTensor.h>
#include <xerus/sparseTensor.h>
#include <xerus/tensorNetwork.h>
#include <xerus/misc/check.h>
#include <xerus/misc/missingFunctions.h>

namespace xerus {
    
	void operator+=(IndexedTensorWritable<Tensor> & _lhs, const IndexedTensorReadOnly<Tensor> & _rhs) {
		std::unique_ptr<Tensor> reorderedRhs(_rhs.tensorObjectReadOnly->construct_new());
		(*reorderedRhs)(_lhs.indices) = _rhs;
		
		if(_lhs.tensorObjectReadOnly->is_sparse()) {
			REQUIRE(reorderedRhs->is_sparse(), "Cannot calculate Sparse += Dense");
			static_cast<SparseTensor&>(*_lhs.tensorObject) += static_cast<SparseTensor&>(*reorderedRhs);
		} else {
			static_cast<FullTensor&>(*_lhs.tensorObject) += *reorderedRhs;
		}		
	}
	
	void operator-=(IndexedTensorWritable<Tensor> & _lhs, const IndexedTensorReadOnly<Tensor> & _rhs) {
		std::unique_ptr<Tensor> reorderedRhs(_rhs.tensorObjectReadOnly->construct_new());
		(*reorderedRhs)(_lhs.indices) = _rhs;
		
		if(_lhs.tensorObjectReadOnly->is_sparse()) {
			REQUIRE(reorderedRhs->is_sparse(), "Cannot calculate Sparse += Dense");
			static_cast<SparseTensor&>(*_lhs.tensorObject) -= static_cast<SparseTensor&>(*reorderedRhs);
		} else {
			static_cast<FullTensor&>(*_lhs.tensorObject) -= *reorderedRhs;
		}
	}
    
    
    IndexedTensorMoveable<Tensor> operator+(const IndexedTensorReadOnly<Tensor> & _lhs, const IndexedTensorReadOnly<Tensor> & _rhs) {
		if(!_lhs.tensorObjectReadOnly->is_sparse() || _rhs.tensorObjectReadOnly->is_sparse()) {
			IndexedTensorMoveable<Tensor> result(_lhs);
			result.perform_traces();
			operator+=(result, _rhs);
			return result;
		} else {
			IndexedTensorMoveable<Tensor> result(_rhs);
			result.perform_traces();
			operator+=(result, _lhs);
			return result;
		}
	}
	
	IndexedTensorMoveable<Tensor> operator+(      IndexedTensorMoveable<Tensor> && _lhs, const IndexedTensorReadOnly<Tensor> &  _rhs) {
		if(!_lhs.tensorObjectReadOnly->is_sparse() || _rhs.tensorObjectReadOnly->is_sparse()) {
			IndexedTensorMoveable<Tensor> result(std::move(_lhs));
			result.perform_traces();
			operator+=(result, _rhs);
			return result;
		} else {
			IndexedTensorMoveable<Tensor> result(_rhs);
			result.perform_traces();
			operator+=(result, _lhs);
			return result;
		}
	}
	
	IndexedTensorMoveable<Tensor> operator+(const IndexedTensorReadOnly<Tensor> &  _lhs,       IndexedTensorMoveable<Tensor> && _rhs){
		if(!_rhs.tensorObjectReadOnly->is_sparse() || _lhs.tensorObjectReadOnly->is_sparse()) {
			IndexedTensorMoveable<Tensor> result(std::move(_rhs));
			result.perform_traces();
			operator+=(result, _lhs);
			return result;
		} else {
			IndexedTensorMoveable<Tensor> result(_lhs);
			result.perform_traces();
			operator+=(result, _rhs);
			return result;
		}
	}
	
	IndexedTensorMoveable<Tensor> operator+(      IndexedTensorMoveable<Tensor> && _lhs,       IndexedTensorMoveable<Tensor> && _rhs){
		if(!_lhs.tensorObjectReadOnly->is_sparse() || _rhs.tensorObjectReadOnly->is_sparse()) {
			IndexedTensorMoveable<Tensor> result(std::move(_lhs));
			result.perform_traces();
			operator+=(result, _rhs);
			return result;
		} else {
			IndexedTensorMoveable<Tensor> result(std::move(_rhs));
			result.perform_traces();
			operator+=(result, _lhs);
			return result;
		}
	}
	
	
	
	IndexedTensorMoveable<Tensor> operator-(const IndexedTensorReadOnly<Tensor> & _lhs, const IndexedTensorReadOnly<Tensor> & _rhs) {
		if(!_lhs.tensorObjectReadOnly->is_sparse() || _rhs.tensorObjectReadOnly->is_sparse()) {
			IndexedTensorMoveable<Tensor> result(_lhs);
			result.perform_traces();
			operator-=(result, _rhs);
			return result;
		} else {
			IndexedTensorMoveable<Tensor> result(_rhs);
			result.perform_traces();
			operator-=(result, _lhs);
			*result.tensorObject *= -1.0;
			return result;
		}
	}
	
	IndexedTensorMoveable<Tensor> operator-(      IndexedTensorMoveable<Tensor> && _lhs, const IndexedTensorReadOnly<Tensor> &  _rhs) {
		if(!_lhs.tensorObjectReadOnly->is_sparse() || _rhs.tensorObjectReadOnly->is_sparse()) {
			IndexedTensorMoveable<Tensor> result(std::move(_lhs));
			result.perform_traces();
			operator-=(result, _rhs);
			return result;
		} else {
			IndexedTensorMoveable<Tensor> result(_rhs);
			result.perform_traces();
			operator-=(result, _lhs);
			*result.tensorObject *= -1.0;
			return result;
		}
	}
	
	IndexedTensorMoveable<Tensor> operator-(const IndexedTensorReadOnly<Tensor> &  _lhs,       IndexedTensorMoveable<Tensor> && _rhs){
		if(!_rhs.tensorObjectReadOnly->is_sparse() || _lhs.tensorObjectReadOnly->is_sparse()) {
			IndexedTensorMoveable<Tensor> result(std::move(_rhs));
			result.perform_traces();
			operator-=(result, _lhs);
			*result.tensorObject *= -1.0;
			return result;
		} else {
			IndexedTensorMoveable<Tensor> result(_lhs);
			result.perform_traces();
			operator-=(result, _rhs);
			return result;
		}
	}
	
	IndexedTensorMoveable<Tensor> operator-(      IndexedTensorMoveable<Tensor> && _lhs,       IndexedTensorMoveable<Tensor> && _rhs){
		if(!_lhs.tensorObjectReadOnly->is_sparse() || _rhs.tensorObjectReadOnly->is_sparse()) {
			IndexedTensorMoveable<Tensor> result(std::move(_lhs));
			result.perform_traces();
			operator-=(result, _rhs);
			return result;
		} else {
			IndexedTensorMoveable<Tensor> result(std::move(_rhs));
			result.perform_traces();
			operator-=(result, _lhs);
			*result.tensorObject *= -1.0;
			return result;
		}
	}
	
	
	IndexedTensorMoveable<Tensor> operator+(const IndexedTensorReadOnly<Tensor> &  _lhs, const IndexedTensorReadOnly<TensorNetwork> &  _rhs) {
		return _lhs+IndexedTensorMoveable<Tensor>(_rhs);
	}
	
	IndexedTensorMoveable<Tensor> operator+(const IndexedTensorReadOnly<TensorNetwork> &  _lhs, const IndexedTensorReadOnly<Tensor> &  _rhs) {
		return IndexedTensorMoveable<Tensor>(_lhs)+_rhs;
	}
	
	IndexedTensorMoveable<Tensor> operator-(const IndexedTensorReadOnly<Tensor> &  _lhs, const IndexedTensorReadOnly<TensorNetwork> &  _rhs) {
		return _lhs-IndexedTensorMoveable<Tensor>(_rhs);
	}
	
	IndexedTensorMoveable<Tensor> operator-(const IndexedTensorReadOnly<TensorNetwork> &  _lhs, const IndexedTensorReadOnly<Tensor> &  _rhs) {
		return IndexedTensorMoveable<Tensor>(_lhs)-_rhs;
	}
	
	
    IndexedTensorMoveable<TensorNetwork> operator+(const IndexedTensorReadOnly<TensorNetwork>  &  _lhs, const IndexedTensorReadOnly<TensorNetwork>  &  _rhs) {
        IndexedTensorMoveable<TensorNetwork> result;
        if(!_lhs.tensorObjectReadOnly->specialized_sum(result, _lhs, _rhs) && !_rhs.tensorObjectReadOnly->specialized_sum(result, _rhs, _lhs)) {
            result.assign(IndexedTensorMoveable<TensorNetwork>(IndexedTensorMoveable<Tensor>(_lhs) + IndexedTensorMoveable<Tensor>(_rhs)));
        }
        return result;
    }
    
    IndexedTensorMoveable<TensorNetwork> operator-(const IndexedTensorReadOnly<TensorNetwork>  & _lhs, const IndexedTensorReadOnly<TensorNetwork>  &  _rhs) {
        return _lhs+(-1*_rhs);
    }
    
    
	void operator+=(IndexedTensorWritable<TensorNetwork> &  _lhs, const IndexedTensorReadOnly<TensorNetwork> &  _rhs) {
		_lhs = _lhs + _rhs;
	}
	
	void operator-=(IndexedTensorWritable<TensorNetwork> &  _lhs, const IndexedTensorReadOnly<TensorNetwork> &  _rhs) {
		_lhs = _lhs - _rhs;
	}
}