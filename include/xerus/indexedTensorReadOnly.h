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

#include "basic.h"
#include <vector>

namespace xerus {
    // Necessary forward declaritons
    class Index;

    
	/**
	 * Class representing any Tensor or TensorNetwork object equipped with an Index order, that can at least be read (i.e. is not nessecarily writeable)
	 */
    template<class tensor_type>
    class IndexedTensorReadOnly {
    public:
        /// Pointer to the associated Tensor/TensorNetwork object
        const tensor_type* tensorObjectReadOnly;
         
        /// Vector of the associates indices 
        std::vector<Index> indices;
        
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Constructors - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    protected:
        /// Creates an empty IndexedTensorReadOnly, should only be used internally
        IndexedTensorReadOnly();
        
    public:
        /// There is no usefull copy constructor.
        IndexedTensorReadOnly(const IndexedTensorReadOnly & _other ) = delete;
        
        /// Move-constructor
        IndexedTensorReadOnly(IndexedTensorReadOnly<tensor_type> && _other );
        
        /// Constructs an IndexedTensorReadOnly using the given pointer and indices
        IndexedTensorReadOnly(const tensor_type* const _tensorObjectReadOnly, const std::vector<Index>& _indices);
        
        /// Constructs an IndexedTensorReadOnly using the given pointer and indices
        IndexedTensorReadOnly(const tensor_type* const _tensorObjectReadOnly, std::vector<Index>&& _indices);
        
    public:
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Destructor - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
        
        /// Destructor must be virtual
        virtual ~IndexedTensorReadOnly();
        
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Others - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
        bool uses_tensor(const tensor_type *otherTensor) const;
        
        size_t degree() const;
        
        bool is_open(const Index& idx) const;
        
        bool is_contained_and_open(const Index& idx) const;
        
        std::vector<size_t> get_evaluated_dimensions(const std::vector<Index>& _indexOrder) const;
        
        std::vector<Index> get_assigned_indices() const;
        
        std::vector<Index> get_assigned_indices(const size_t _futureDegree, const bool _assignDimensions = false) const;
        
        #ifndef DISABLE_RUNTIME_CHECKS_
            /// Checks whether the indices are usefull in combination with the current dimensions
            void check_indices(const bool _allowNonOpen = true) const;
        #endif
    };
    
    template<class tensor_type>
    value_t frob_norm(const IndexedTensorReadOnly<tensor_type>& _idxTensor);
    
    size_t get_eval_degree(const std::vector<Index>& _indices);
}

