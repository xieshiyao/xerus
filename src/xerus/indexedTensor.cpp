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
 * @brief Implementation of the IndexedTensor class.
 */

#include <xerus/indexedTensor.h>

#include <xerus/index.h>
#include <xerus/misc/missingFunctions.h>
#include <xerus/tensor.h>
#include <xerus/tensorNetwork.h>

namespace xerus {
    template<class tensor_type>
    IndexedTensor<tensor_type>::IndexedTensor(IndexedTensor &&_other ) : IndexedTensorWritable<tensor_type>(std::move(_other)) { }

    template<class tensor_type>
    IndexedTensor<tensor_type>::IndexedTensor(tensor_type* const _tensorObject, const std::vector<Index>& _indices, const bool _takeOwnership) :
        IndexedTensorWritable<tensor_type>(_tensorObject, _indices, _takeOwnership) {}
        
    template<class tensor_type>
    IndexedTensor<tensor_type>::IndexedTensor(tensor_type* const _tensorObject, std::vector<Index>&& _indices, const bool _takeOwnership) :
        IndexedTensorWritable<tensor_type>(_tensorObject, std::move(_indices), _takeOwnership) {}
        
        
    template<class tensor_type>
    void IndexedTensor<tensor_type>::operator=(IndexedTensor<tensor_type>&& _rhs) {
        operator=(static_cast<IndexedTensorReadOnly<tensor_type>&&>(_rhs));
    }
    
    template<class tensor_type>
    void IndexedTensor<tensor_type>::operator=(IndexedTensorReadOnly<Tensor>&& _rhs) {
        static_cast<IndexedTensorWritable<tensor_type>&&>(*this) = std::move(_rhs);
    }
    
    template<class tensor_type>
    void IndexedTensor<tensor_type>::operator=(IndexedTensorReadOnly<TensorNetwork>&& _rhs) {
		static_cast<IndexedTensorWritable<tensor_type>&&>(*this) = std::move(_rhs);
    }
    
    // IndexedTensorReadOnly may be instanciated as
    template class IndexedTensor<Tensor>;
    template class IndexedTensor<TensorNetwork>;
    
}
