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
    // Necessary forward declaritons
    class Tensor;
    class TensorNetwork;
    
    template<class tensor_type>
    class IndexedTensorWritable : public IndexedTensorReadOnly<tensor_type> {
    public:
        // Non-const pointer to the tensor object. Must always coincide with tensorObject
        tensor_type* tensorObject;
        
        /// Flag indicating, that the indexedTensor is the owner of the tensor object
        bool deleteTensorObject;
        
        /// Creates an empty indexed Tensor, should only be used internally
    protected:
        IndexedTensorWritable();
        
    public:
        /// There is no usefull copy constructor, because the handling of the tensorObject is unclear
        IndexedTensorWritable(const IndexedTensorWritable &_other ) = delete;
        
        /// Move constructor
        IndexedTensorWritable(IndexedTensorWritable &&_other );
        
        /// Constructs an IndexedTensorWritable with the given tensor and takes owership of the tensorObject if requested
        IndexedTensorWritable(tensor_type* const _tensorObject, const std::vector<Index>& _indices, const bool _takeOwnership);
        
        /// Constructs an IndexedTensorWritable with the given tensor and takes owership of the tensorObject if requested
        IndexedTensorWritable(tensor_type* const _tensorObject, std::vector<Index>&& _indices, const bool _takeOwnership);
        
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Destructor - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
        
        virtual ~IndexedTensorWritable();
        
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Other - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
        
        bool is_owner() const;
        
        void delete_if_owner();
        
        /// Replacement for the original = operator
        void copy_assign(const IndexedTensorWritable&  _other);
        
        /// Replacement for the original = operator
        void assign(      IndexedTensorWritable&& _other);
        
        /// Resets the IndexedTensorWritable
        void reset(tensor_type* const _tensorObject, const std::vector<Index>& _indices, const bool _takeOwnership = false);
        
        /// Resets the IndexedTensorWritable
        void reset(tensor_type* const _tensorObject, std::vector<Index>&& _indices, const bool _takeOwnership = false);
        
        // Assignment operators -- Used for tensor assignment WITH indices (i.e. in general the LHS and RHS indexTensors do NOT have the same indices)
        // Implementation is provided in the corresponding *_specialization.hpp files
        void operator=(const IndexedTensorReadOnly<Tensor>&         _rhs);
        void operator=(const IndexedTensorReadOnly<TensorNetwork>&  _rhs);
        
        // NOTE: The following would be deleted due to move constructor and is therefore implemented here
        void operator=(const IndexedTensorWritable<tensor_type>& _rhs);
    };
}
