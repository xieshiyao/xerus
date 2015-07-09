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
 * @brief Implementation of the IndexedTensorWritable class.
 */

#include <xerus/indexedTensorWritable.h>

#include <xerus/index.h>
#include <xerus/misc/missingFunctions.h>
#include <xerus/tensor.h>
#include <xerus/tensorNetwork.h>

namespace xerus {
    template<class tensor_type>
    IndexedTensorWritable<tensor_type>::IndexedTensorWritable() : IndexedTensorReadOnly<tensor_type>(), tensorObject(nullptr), deleteTensorObject(false) {}
    
    template<class tensor_type>
    IndexedTensorWritable<tensor_type>::IndexedTensorWritable(IndexedTensorWritable &&_other ) : IndexedTensorReadOnly<tensor_type>(std::move(_other)), tensorObject(_other.tensorObject), deleteTensorObject(_other.deleteTensorObject) {
        // Take ownership
        _other.deleteTensorObject = false;
    }
    
    template<class tensor_type>
    IndexedTensorWritable<tensor_type>::IndexedTensorWritable(tensor_type* const _tensorObject, const std::vector<Index>& _indices, const bool _takeOwnership) :
        IndexedTensorReadOnly<tensor_type>(_tensorObject, _indices), tensorObject(_tensorObject), deleteTensorObject(_takeOwnership) {}
        
        
    template<class tensor_type>
    IndexedTensorWritable<tensor_type>::IndexedTensorWritable(tensor_type* const _tensorObject, std::vector<Index>&& _indices, const bool _takeOwnership) :
        IndexedTensorReadOnly<tensor_type>(_tensorObject, std::move(_indices)), tensorObject(_tensorObject), deleteTensorObject(_takeOwnership) {}

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Destructor - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    template<class tensor_type>
    IndexedTensorWritable<tensor_type>::~IndexedTensorWritable() { 
        if(deleteTensorObject) {
            delete this->tensorObject;
        }
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Other - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    template<class tensor_type>
    bool IndexedTensorWritable<tensor_type>::is_owner() const {
        return deleteTensorObject;
    }
    
    template<class tensor_type>
    void IndexedTensorWritable<tensor_type>::delete_if_owner() {
        if(deleteTensorObject) {
            delete this->tensorObject;
        }
    }
    
    template<class tensor_type>
    void IndexedTensorWritable<tensor_type>::assign(      IndexedTensorWritable&& _other) {
        this->tensorObject = _other.tensorObject;
        this->tensorObjectReadOnly = _other.tensorObjectReadOnly;
        this->indices = std::move(_other.indices);
        this->deleteTensorObject = _other.deleteTensorObject;
        _other.deleteTensorObject = false;
    }
    
    template<class tensor_type>
    void IndexedTensorWritable<tensor_type>::reset(tensor_type* const _tensorObject, const std::vector<Index>& _indices, const bool _takeOwnership) {
        // Delete old tensorObject
        delete_if_owner();
        
        // Set new parameters
        this->tensorObject = _tensorObject;
        this->tensorObjectReadOnly = _tensorObject;
        this->indices = _indices;
        this->deleteTensorObject = _takeOwnership;
    }
    
    template<class tensor_type>
    void IndexedTensorWritable<tensor_type>::reset(tensor_type* const _tensorObject, std::vector<Index>&& _indices, const bool _takeOwnership) {
        // Delete old tensorObject
        delete_if_owner();
        
        // Set new parameters
        this->tensorObject = _tensorObject;
        this->tensorObjectReadOnly = _tensorObject;
        this->indices = std::move(_indices);
        this->deleteTensorObject = _takeOwnership;
    }
    
    template<class tensor_type>
    void IndexedTensorWritable<tensor_type>::operator=(const IndexedTensorWritable<tensor_type>& _rhs) {
        operator=(static_cast<const IndexedTensorReadOnly<tensor_type> &>(_rhs));
    }
    
    template<>
    void IndexedTensorWritable<Tensor>::perform_traces() {
		REQUIRE(deleteTensorObject, "IndexedTensorMoveable must own its tensor object");
		const std::vector<Index> assIndices = this->get_assigned_indices();
		std::vector<Index> openIndices;
		bool allOpen = true;
		for(const Index& idx : assIndices) {
			if(idx.open()) {
				openIndices.push_back(idx);
			} else {
				allOpen = false;
			}
		}
		if(!allOpen) { 
			(*this->tensorObject)(openIndices) = *this;
			this->indices = openIndices;
		}
	}
    
    // IndexedTensorReadOnly may be instanciated as
    template class IndexedTensorWritable<Tensor>;
    template class IndexedTensorWritable<TensorNetwork>;
    
}
