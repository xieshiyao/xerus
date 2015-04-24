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

#include "indexedTensorReadOnly.h"

namespace xerus {
    template<class tensor_type, tensor_type_restrictions>
    class IndexedTensorWritable : public IndexedTensorReadOnly<tensor_type> {
    public:
        // Non-const pointer to the tensor object. Must always coincide with tensorObject
        tensor_type* tensorObject;
        
        /// Flag indicating, that the indexedTensor is the owner of the tensor object
        bool deleteTensorObject;
        
        /// Creates an empty indexed Tensor, should only be used internally
    protected:
        IndexedTensorWritable() : IndexedTensorReadOnly<tensor_type>(), tensorObject(nullptr), deleteTensorObject(false) {};
        
    public:
        /// There is no usefull copy constructor, because the handling of the tensorObject is unclear
        IndexedTensorWritable(const IndexedTensorWritable &_other ) = delete;
        
        /// Move constructor
        IndexedTensorWritable(IndexedTensorWritable &&_other ) : IndexedTensorReadOnly<tensor_type>(std::move(_other)), tensorObject(_other.tensorObject), deleteTensorObject(_other.deleteTensorObject) {
            // Take ownership
            _other.deleteTensorObject = false;
        }
        
        /// Constructs an IndexedTensorWritable with the given tensor and takes owership of the tensorObject if requested
        ALLOW_MOVE(std::vector<Index>, IdxVec)
        IndexedTensorWritable(tensor_type* const _tensorObject, IdxVec&& _indices, const bool _takeOwnership) :
            IndexedTensorReadOnly<tensor_type>(_tensorObject, std::forward<IdxVec>(_indices)), tensorObject(_tensorObject), deleteTensorObject(_takeOwnership) {}

        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Destructor - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
        
        virtual ~IndexedTensorWritable() { 
            if(deleteTensorObject) {
                delete this->tensorObject;
            }
        }
        
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Other - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
        
        bool is_owner() const {
            return deleteTensorObject;
        }
        
        void delete_if_owner() {
            if(deleteTensorObject) {
                delete this->tensorObject;
            }
        }
        
        /// Replacement for the original = operator
        void copy_assign(const IndexedTensorWritable&  _other) {
            if(_other.deleteTensorObject) {
                tensorObject = _other.tensorObject.deep_copy();
            } else {
                tensorObject = _other.tensorObject;
            }
            this->tensorObjectReadOnly = tensorObject;
            this->indices = _other.indices;
            this->deleteTensorObject = _other.deleteTensorObject;
        }
        
        /// Replacement for the original = operator
        void assign(      IndexedTensorWritable&& _other) {
            this->tensorObject = _other.tensorObject;
            this->tensorObjectReadOnly = _other.tensorObjectReadOnly;
            this->indices = std::move(_other.indices);
            this->deleteTensorObject = _other.deleteTensorObject;
            _other.deleteTensorObject = false;
        }
        
        /// Resets the IndexedTensorWritable
        ALLOW_MOVE(std::vector<Index>, IdxVec)
        void reset(tensor_type* const _tensorObject, IdxVec&& _indices, const bool _takeOwnership = false) {
            // Delete old tensorObject
            delete_if_owner();
            
            // Set new parameters
            this->tensorObject = _tensorObject;
            this->tensorObjectReadOnly = _tensorObject;
            this->indices = std::forward<IdxVec>(_indices);
            this->deleteTensorObject = _takeOwnership;
        }
        
        // Assignment operators -- Used for tensor assignment WITH indices (i.e. in general the LHS and RHS indexTensors do NOT have the same indices)
        // Implementation is provided in the corresponding *_specialization.hpp files
        void operator=(const IndexedTensorReadOnly<Tensor>&         _rhs);
        void operator=(const IndexedTensorReadOnly<TensorNetwork>&  _rhs);
        
        // NOTE: The following would be deleted due to move constructor and is therefore implemented here
        _inline_ void operator=(const IndexedTensorWritable<tensor_type>& _rhs) {
            operator=(static_cast<const IndexedTensorReadOnly<tensor_type> &>(_rhs));
        }
    };
}
