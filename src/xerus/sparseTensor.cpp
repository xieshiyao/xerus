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
 * @brief Implementation of the SparseTensor class.
 */

#include <xerus/sparseTensor.h>

namespace xerus {
	SparseTensor::SparseTensor() : Tensor(Representation::Sparse) { }

    SparseTensor::SparseTensor( const Tensor&  _other) : Tensor(_other) {
		use_sparse_representation();
	}
    
    SparseTensor::SparseTensor(       Tensor&& _other) : Tensor(std::move(_other)) {
		use_sparse_representation();
	}
    
    SparseTensor::SparseTensor(const std::vector<size_t> & _dimensions) : Tensor(_dimensions, Representation::Sparse) { }
        
	SparseTensor::SparseTensor(      std::vector<size_t>&& _dimensions) : Tensor(std::move(_dimensions), Representation::Sparse) { }
    
    SparseTensor::SparseTensor(std::initializer_list<size_t>&& _dimensions) : Tensor(std::move(_dimensions), Representation::Sparse) { }
    
    SparseTensor::SparseTensor(const FullTensor& _full, const double _eps) : Tensor(_full) {
		denseData.reset();
		representation = Representation::Sparse;
		sparseData.reset(new std::map<size_t, value_t>());
        for(size_t i = 0; i < _full.size; ++i) {
            if(std::abs(_full[i]) >= _eps) {
                sparseData->insert({i, _full[i]});
//                 sparseData->emplace_hint(sparseData->end(), i, fullDataPtr[i]); TODO use this instead
            }
        }
    }
    
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Standard operators - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
//     SparseTensor& SparseTensor::operator=(const SparseTensor&  _other) {
//         assign(_other);
// 		sparseData = _other.sparseData;
//         return *this;
//     }
//     
//     SparseTensor& SparseTensor::operator=(      SparseTensor&& _other) {
//         assign(std::move(_other));
// 		sparseData = _other.sparseData;
//         return *this;
//     }
    
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Access - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    
    value_t SparseTensor::at(const size_t _position) const {
        const std::map<size_t, value_t>::const_iterator entry = sparseData->find(_position);
        if(entry == sparseData->end()) {
            return 0;
        } else {
            return factor*entry->second;
        }
    }
    
    value_t SparseTensor::at(const std::vector<size_t>& _indices) const {
        CHECK(_indices.size() == dimensions.size(), fatal, "Wrong number of indices given");
        size_t finalIndex = 0;
        for(size_t i = 0; i < _indices.size(); ++i) {
            REQUIRE(_indices[i] < dimensions[i], "Index "<< i <<" out of bounds " << _indices[i] << " >=! " << dimensions[i]);
            finalIndex *= dimensions[i];
            finalIndex += _indices[i];
        }
        
        const std::map<size_t, value_t>::const_iterator entry = sparseData->find(finalIndex);
        if(entry == sparseData->end()) {
            return 0;
        } else {
            return factor*entry->second;
        }
    }
    
    
}
