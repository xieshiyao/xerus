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

#include "assignedIndices.h"

#define tensor_type_restrictions typename std::enable_if< \
               std::is_same<Tensor,        typename std::decay<tensor_type>::type>{}  \
            || std::is_same<TensorNetwork, typename std::decay<tensor_type>::type>{},  \
        int>::type = 0

namespace xerus {
    // Forward declaration of the two allowed types
    class Tensor;
    class TensorNetwork;
    
	/**
	 * Class representing any Tensor or TensorNetwork object equipped with an Index order, that can at least be read (i.e. is not nessecarily writeable)
	 */
    template<class tensor_type, tensor_type_restrictions>
    class IndexedTensorReadOnly {
    public:
        /// Pointer to the associated Tensor/TensorNetwork object
        const tensor_type* tensorObjectReadOnly;
         
        /// Vector of the associates indices 
        std::vector<Index> indices;
        
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Constructors - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    protected:
        /// Creates an empty IndexedTensorReadOnly, should only be used internally
        IndexedTensorReadOnly() : tensorObjectReadOnly(nullptr) {};
        
    public:
        /// There is no usefull copy constructor.
        IndexedTensorReadOnly(const IndexedTensorReadOnly & _other ) = delete;
        
        /// Move-constructor
        IndexedTensorReadOnly(IndexedTensorReadOnly<tensor_type> && _other ) :
            tensorObjectReadOnly(_other.tensorObjectReadOnly),
            indices(std::move(_other.indices))
            { }
        
        /// Constructs an IndexedTensorReadOnly using the given pointer and indices
        ALLOW_MOVE(std::vector<Index>, VECTOR)
        IndexedTensorReadOnly(const tensor_type* const _tensorObjectReadOnly, VECTOR&& _indices) : tensorObjectReadOnly(_tensorObjectReadOnly), indices(std::forward<VECTOR>(_indices)) { }
        
    public:
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Destructor - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
        
        /// Destructor must be virtual
        virtual ~IndexedTensorReadOnly() { }
        
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Others - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
        bool uses_tensor(const tensor_type *otherTensor) const {
            return otherTensor == tensorObjectReadOnly;
        }
        
        size_t degree() const {
            return tensorObjectReadOnly->degree();
        }
        
        bool is_open(const Index& idx) const {
            REQUIRE(idx.fixed() || contains(indices, idx), "Index " << idx << " not contained in indices: " << indices); // TODO would be nice to check that fixed indices are also contained...
            return !idx.fixed() && count(indices, idx) == 1;
        }
        
        bool is_contained_and_open(const Index& idx) const {
            return !idx.fixed() && count(indices, idx) == 1;
        }
        
        std::vector<size_t> get_evaluated_dimensions(const std::vector<Index>& _indexOrder) const {
            std::vector<size_t> evalDimensions;
            evalDimensions.reserve(_indexOrder.size());
            for(const Index& idx : _indexOrder) {
                REQUIRE(count(indices, idx) == 1, "All indices of evaluation target must appear exactly once.");
                
                // Find index
                size_t indexPos = 0, dimCount = 0;
                while(indices[indexPos] != idx) { dimCount += indices[indexPos].flags[Index::Flag::INVERSE_SPAN] ? degree()-indices[indexPos].span : indices[indexPos].span; ++indexPos; }
                
                // Calculate span
                size_t span = indices[indexPos].flags[Index::Flag::INVERSE_SPAN] ? degree()-indices[indexPos].span : indices[indexPos].span;
                
                REQUIRE(dimCount+span <= tensorObjectReadOnly->dimensions.size(), "Order determined by Indices is to large.");
                
                // Insert dimensions
                for(size_t i = 0; i < span; ++i) {
                    evalDimensions.emplace_back(tensorObjectReadOnly->dimensions[dimCount+i]);
                }
            }
            return evalDimensions;
        }
        
        
        std::vector<Index> get_assigned_indices() const {
            return get_assigned_indices(this->degree(), true);
        }
        
        std::vector<Index> get_assigned_indices(const size_t _futureDegree, const bool _assignDimensions = false) const {
            std::vector<Index> assignedIndices;
            assignedIndices.reserve(indices.size());
            
            size_t dimensionCount = 0;
            for(const Index& idx : indices) {
                // We don't look at span zero indices
                if((   !idx.flags[Index::Flag::INVERSE_SPAN] && idx.span == 0) 
                    || (idx.flags[Index::Flag::INVERSE_SPAN] && idx.span == _futureDegree) 
                    || (idx.flags[Index::Flag::FRACTIONAL_SPAN] && _futureDegree/idx.span == 0)) 
                {
                    REQUIRE(!idx.fixed(), "Fixed index must have span == 1");
                    continue;
                }
                
                // Fixed indices want special treatnment
                if(idx.fixed()) {
                    REQUIRE(!idx.flags[Index::Flag::INVERSE_SPAN], "Fixed index must not have inverse span.");
                    REQUIRE(!idx.flags[Index::Flag::FRACTIONAL_SPAN], "Fixed index must not have fractional span.");
                    REQUIRE(idx.span == 1, "Fixed index must have span one.");
                    if(_assignDimensions) {
                        assignedIndices.emplace_back(idx.valueId, 1, tensorObjectReadOnly->dimensions[dimensionCount++], Index::Flag::OPEN, Index::Flag::FIXED, false);
                    } else {
                        assignedIndices.emplace_back(idx.valueId, 1, 0, Index::Flag::OPEN, Index::Flag::FIXED, false);
                        dimensionCount++;
                    }
                } else {
                    // Set span
                    size_t span;
                    if(idx.flags[Index::Flag::INVERSE_SPAN]) {
                        REQUIRE(idx.span < _futureDegree, "Index used with variable span (e.g. i&3) would have negative span " << _futureDegree << " - " << idx.span << " = " << _futureDegree - idx.span << "!");
                        span = _futureDegree-idx.span;
                    } else if(idx.flags[Index::Flag::FRACTIONAL_SPAN]) {
                        CHECK(_futureDegree%idx.span == 0, warning, "Fractional span used in tensor with an degree that is not divisable by the given fraction.");
                        span = _futureDegree/idx.span;
                    } else {
                        span = idx.span;
                    }
                    
                    // Calculate multDimension
                    size_t multDimension = 1;
                    if(_assignDimensions) {
                        REQUIRE(dimensionCount+span <= tensorObjectReadOnly->dimensions.size(), "Order determined by Indices is to large.");
                        for(size_t i = 0; i < span; ++i) {
                            multDimension *= tensorObjectReadOnly->dimensions[dimensionCount++];
                        }
                    } else {
                        dimensionCount += span;
                    }
                    
                    // Determine whether index is open
                    bool open = true;
                    for(size_t i = 0; i < assignedIndices.size(); ++i) {
                        if(idx == assignedIndices[i]) {
                            REQUIRE(assignedIndices[i].open(), "An index must not appere more than twice!");
                            assignedIndices[i].open(false);
                            open = false;
                            break;
                        }
                    }
                    
                    assignedIndices.emplace_back(idx.valueId, span, multDimension, Index::Flag::OPEN, open);
                }
            }
            
            REQUIRE(dimensionCount >= _futureDegree, "Order determined by Indices is to small. Order according to the indices " << dimensionCount << ", according to the tensor " << _futureDegree);
            REQUIRE(dimensionCount <= _futureDegree, "Order determined by Indices is to large. Order according to the indices " << dimensionCount << ", according to the tensor " << _futureDegree);
            return assignedIndices;
        }
        
        AssignedIndices assign_indices() const {
            AssignedIndices assignedIndices;
            
            assignedIndices.allIndicesOpen = true;
            assignedIndices.numIndices = 0;
            assignedIndices.indices.reserve(indices.size());
            assignedIndices.indexDimensions.reserve(indices.size());
            
            size_t dimensionCount = 0;
            for(const Index& idx : indices) {
                // We don't look at span zero indices
                if((!idx.flags[Index::Flag::INVERSE_SPAN] && idx.span == 0) || (idx.flags[Index::Flag::INVERSE_SPAN] && idx.span == degree()) || (idx.flags[Index::Flag::FRACTIONAL_SPAN] && degree()/idx.span == 0)) {
                    REQUIRE(!idx.fixed(), "Fixed index must have span == 1");
                    continue;
                }
                
                // Fixed indices want special treatnment
                if(idx.fixed()) {
                    REQUIRE(!idx.flags[Index::Flag::INVERSE_SPAN], "Fixed index must not have inverse span.");
                    REQUIRE(!idx.flags[Index::Flag::FRACTIONAL_SPAN], "Fixed index must not have fractional span.");
                    REQUIRE(idx.span == 1, "Fixed index must have span one.");
                    
                    assignedIndices.numIndices++;
                    assignedIndices.indices.emplace_back(idx.valueId, 1, tensorObjectReadOnly->dimensions[dimensionCount], Index::Flag::FIXED, Index::Flag::OPEN, true, false);
                    assignedIndices.indexDimensions.emplace_back(tensorObjectReadOnly->dimensions[dimensionCount++]);
                    assignedIndices.allIndicesOpen = false;
                } else {
                    // Set span
                    size_t span;
                    if(idx.flags[Index::Flag::INVERSE_SPAN]) {
                        REQUIRE(idx.span < degree(), "Index used with variable span (e.g. i&3) would have negative span " << degree() << " - " << idx.span << " = " << degree() - idx.span << "!");
                        span = degree()-idx.span;
                    } else if(idx.flags[Index::Flag::FRACTIONAL_SPAN]) {
                        CHECK(degree()%idx.span == 0, warning, "Fractional span used in tensor with an degree that is not divisable by the given fraction.");
                        span = degree()/idx.span;
                    } else {
                        span = idx.span;
                    }
                    
                    // Calculate multDimension
                    REQUIRE(dimensionCount+span <= tensorObjectReadOnly->dimensions.size(), "Order determined by Indices is to large.");
                    size_t multDimension = 1;
                    for(size_t i=0; i < span; ++i) {
                        multDimension *= tensorObjectReadOnly->dimensions[dimensionCount++];
                    }
                    
                    
                    // Determine whether index is open
                    bool open = true;
                    for(size_t i = 0; i < assignedIndices.indices.size(); ++i) {
                        if(idx == assignedIndices.indices[i]) {
                            REQUIRE(assignedIndices.indices[i].open(), "An index must not appere more than twice!");
                            assignedIndices.indices[i].open(false);
                            open = false;
                            break;
                        }
                    }
                    
                    assignedIndices.numIndices++;
                    assignedIndices.indices.emplace_back(idx.valueId, span, multDimension, Index::Flag::OPEN, open);
                    assignedIndices.indexDimensions.emplace_back(multDimension);
                    assignedIndices.allIndicesOpen = assignedIndices.allIndicesOpen && open; 
                }
            }
            
            REQUIRE(assignedIndices.numIndices == assignedIndices.indices.size() && assignedIndices.numIndices == assignedIndices.indexDimensions.size(), "Internal Error"); 
            REQUIRE(dimensionCount >= degree(), "Order determined by Indices is to small. Order according to the indices " << dimensionCount << ", according to the tensor " << degree());
            REQUIRE(dimensionCount <= degree(), "Order determined by Indices is to large. Order according to the indices " << dimensionCount << ", according to the tensor " << degree());
            return assignedIndices;
        }
        
        #ifndef DISABLE_RUNTIME_CHECKS_
            /// Checks whether the indices are usefull in combination with the current dimensions
            void check_indices(const bool _allowNonOpen = true) const {
                size_t dimensionCount = 0;
                for(const Index& idx : indices) {
                    REQUIRE(_allowNonOpen || !idx.fixed(), "Fixed indices are not allowed here.");
                    REQUIRE(!idx.fixed() || idx.span == 1, "Fixed index must have span == 1");
                    REQUIRE(_allowNonOpen || count(indices, idx) == 1, "Traces are not allowed here.");
                    REQUIRE(idx.fixed() || count(indices, idx) <= 2, "An index must not appere more than twice!");
                    
                    if(idx.flags[Index::Flag::INVERSE_SPAN]) {
                        REQUIRE(idx.span <= degree(), "Index used with variable span (e.g. i&3) would have negative span " << (long)(degree() - idx.span) << "!");
                        REQUIRE(!idx.fixed(), "Fixed index must not have inverse span");
                        dimensionCount += degree()-idx.span;
                    } else if(idx.flags[Index::Flag::FRACTIONAL_SPAN]) {
                        CHECK(degree()%idx.span == 0, warning, "Fractional span used in tensor with an degree that is not divisable by the given fraction.");
                        dimensionCount += degree()/idx.span;
                    } else {
                        dimensionCount += idx.span;
                    }
                }
                REQUIRE(dimensionCount >= degree(), "Order determined by Indices is to small. Order according to the indices " << dimensionCount << ", according to the tensor " << degree());
                REQUIRE(dimensionCount <= degree(), "Order determined by Indices is to large. Order according to the indices " << dimensionCount << ", according to the tensor " << degree());
            }
        #endif
    };

    template<class T>
    _inline_ value_t frob_norm(const IndexedTensorReadOnly<T>& _A) {
        return _A.tensorObjectReadOnly->frob_norm(); 
    }
    
    _inline_ size_t get_eval_degree(const std::vector<Index>& _indices) {
        size_t degree = 0;
        for(const Index& idx : _indices) {
            REQUIRE(!idx.flags[Index::Flag::INVERSE_SPAN], "Internal Error");
            if(!idx.fixed() && count(_indices, idx) != 2) { degree += idx.span; }
        }
        return degree;
    }
}

