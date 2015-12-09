// // Xerus - A General Purpose Tensor Library
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
 * @brief Implementation of the FullTensor class.
 */

#include <xerus/tensor.h>
#include <xerus/fullTensor.h>
#include <xerus/sparseTensor.h>
#include <xerus/tensorNetwork.h>
#include <xerus/blasLapackWrapper.h>
#include <xerus/selectedFunctions.h>
#include <cstring>

namespace xerus {
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Constructors - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    FullTensor::FullTensor() : FullTensor(std::vector<size_t>({})) { }

    FullTensor::FullTensor(const FullTensor&  _other) : Tensor(_other) { }

    FullTensor::FullTensor(      FullTensor&& _other) : Tensor(std::move(_other)){ }
    
    FullTensor::FullTensor(const Tensor&  _other) : Tensor(_other) {
		if(_other.is_sparse()) {
			denseData.reset(new value_t[size], internal::array_deleter_vt);
			misc::array_set_zero(denseData.get(), size);
			for(const std::pair<size_t, value_t>& entry : *static_cast<const SparseTensor&>(_other).sparseData) {
				denseData.get()[entry.first] = entry.second;
			}
			sparseData.reset();
			representation = Representation::Dense;
		}
	}

    FullTensor::FullTensor(      Tensor&& _other) : Tensor(std::move(_other)) {
		if(_other.is_sparse()) {
			denseData.reset(new value_t[size], internal::array_deleter_vt);
			misc::array_set_zero(denseData.get(), size);
			for(const std::pair<size_t, value_t>& entry : *sparseData) {
				denseData.get()[entry.first] = entry.second;
			}
			sparseData.reset();
			representation = Representation::Dense;
		}
	}
    
    FullTensor::FullTensor(const TensorNetwork& _other) : FullTensor(*_other.fully_contracted_tensor()) { }

    FullTensor::FullTensor(const size_t _degree) : FullTensor(std::vector<size_t>(_degree, 1)) { }
    
    FullTensor::FullTensor(const std::vector<size_t>&  _dimensions, _unused_ DONT_SET_ZERO) : Tensor(_dimensions) { }
        
    FullTensor::FullTensor(      std::vector<size_t>&& _dimensions, _unused_ DONT_SET_ZERO) : Tensor(std::move(_dimensions)) { }
    
    FullTensor::FullTensor(const std::vector<size_t>&  _dimensions) : FullTensor(_dimensions, DONT_SET_ZERO()) {
        misc::array_set_zero(denseData.get(), size);
    }
    
    FullTensor::FullTensor(      std::vector<size_t>&& _dimensions) : FullTensor(std::move(_dimensions), DONT_SET_ZERO()) {
        misc::array_set_zero(denseData.get(), size);
    }
    
    FullTensor::FullTensor(const std::vector<size_t> & _dimensions, std::unique_ptr<value_t[]>&& _data) : Tensor(_dimensions) {
		denseData.reset(_data.release(), internal::array_deleter_vt);
	}
        
    FullTensor::FullTensor(      std::vector<size_t>&& _dimensions, std::unique_ptr<value_t[]>&& _data) : Tensor(std::move(_dimensions)) {
		denseData.reset(_data.release(), internal::array_deleter_vt);
	}
    
    

        
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Standard operators - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    FullTensor& FullTensor::operator=(const FullTensor& _other) {
		assign(_other); // Delegate to Tensor
		denseData = _other.denseData;
        return *this;
    }

    FullTensor& FullTensor::operator=(FullTensor&& _other) {
        std::swap(dimensions, _other.dimensions);
        std::swap(size, _other.size);
        factor = _other.factor;
		std::swap(denseData, _other.denseData);
        return *this;
    }
    
    FullTensor& FullTensor::operator=(const Tensor& _other) {
		if(_other.is_sparse()) {
			reset(_other.dimensions);
			operator+=(_other);
		} else {
			operator=(static_cast<const FullTensor&>(_other));
		}
        return *this;
    }

    FullTensor& FullTensor::operator=(Tensor&& _other) {
        if(_other.is_sparse()) {
			reset(_other.dimensions);
			operator+=(_other);
		} else {
			operator=(std::move(static_cast<FullTensor&&>(_other)));
		}
        return *this;
    }



}


