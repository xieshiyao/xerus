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
 * @brief Implementation of the IndexedTensorReadOnly class.
 */

#include <xerus/indexedTensorReadOnly.h>
#include <xerus/indexedTensorMoveable.h>

#include <xerus/index.h>
#include <xerus/misc/missingFunctions.h>
#include <xerus/misc/check.h>
#include <xerus/tensor.h>
#include <xerus/tensorNetwork.h>

namespace xerus {
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Constructors - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    template<class tensor_type>
    IndexedTensorReadOnly<tensor_type>::IndexedTensorReadOnly() : tensorObjectReadOnly(nullptr) { }

    template<class tensor_type>
    IndexedTensorReadOnly<tensor_type>::IndexedTensorReadOnly(IndexedTensorReadOnly<tensor_type> && _other ) :
        tensorObjectReadOnly(_other.tensorObjectReadOnly),
        indices(std::move(_other.indices))
        { }
    
    
    template<class tensor_type>
    IndexedTensorReadOnly<tensor_type>::IndexedTensorReadOnly(const tensor_type* const _tensorObjectReadOnly, const std::vector<Index>& _indices)
        : tensorObjectReadOnly(_tensorObjectReadOnly), indices(_indices) { }
        
    template<class tensor_type>
    IndexedTensorReadOnly<tensor_type>::IndexedTensorReadOnly(const tensor_type* const _tensorObjectReadOnly, std::vector<Index>&& _indices)
        : tensorObjectReadOnly(_tensorObjectReadOnly), indices(std::move(_indices)) { }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Destructor - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    template<class tensor_type>
    IndexedTensorReadOnly<tensor_type>::~IndexedTensorReadOnly() { }
    
    
	/*- - - - - - - - - - - - - - - - - - - - - - - - - - Aritmetic Operators - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
//     template<class tensor_type>
//     IndexedTensorMoveable<tensor_type> IndexedTensorReadOnly<tensor_type>::operator*(const value_t _factor) const {
// 		IndexedTensorMoveable<tensor_type> result(*this);
//         *result.tensorObject *= _factor;
//         return result;
// 	}
// 	
// 	template<class tensor_type>
//     IndexedTensorMoveable<tensor_type> IndexedTensorReadOnly<tensor_type>::operator/(const value_t _divisor) const {
// 		IndexedTensorMoveable<tensor_type> result(*this);
//         *result.tensorObject /= _divisor;
//         return result;
// 	}
	
	template<class tensor_type>
    IndexedTensorReadOnly<tensor_type>::operator value_t() const {
		REQUIRE(degree() == 0, "cannot cast tensors of degree > 0 to value_t. did you mean frob_norm() or similar?");
		return (*tensorObjectReadOnly)[0];
	}
	
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Others - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    template<class tensor_type>
    bool IndexedTensorReadOnly<tensor_type>::uses_tensor(const tensor_type *otherTensor) const {
        return otherTensor == tensorObjectReadOnly;
    }
    
    template<class tensor_type>
    size_t IndexedTensorReadOnly<tensor_type>::degree() const {
        return tensorObjectReadOnly->degree();
    }
    
    template<class tensor_type>
    bool IndexedTensorReadOnly<tensor_type>::is_contained_and_open(const Index& idx) const {
        return !idx.fixed() && misc::count(indices, idx) == 1;
    }
    
    template<class tensor_type>
    std::vector<size_t> IndexedTensorReadOnly<tensor_type>::get_evaluated_dimensions(const std::vector<Index>& _indexOrder) const {
        std::vector<size_t> evalDimensions;
        evalDimensions.reserve(_indexOrder.size());
		
        for(const Index& idx : _indexOrder) {
            REQUIRE(misc::count(indices, idx) == 1, "All indices of evaluation target must appear exactly once.");
            
            // Find index
            size_t indexPos = 0, dimCount = 0;
            while(indices[indexPos] != idx) {
				dimCount += indices[indexPos++].actual_span(degree());
			}
            
            // Calculate span
            const size_t span = indices[indexPos].actual_span(degree());
            
            REQUIRE(dimCount+span <= tensorObjectReadOnly->dimensions.size(), "Order determined by Indices is to large. Tensor has " << tensorObjectReadOnly->dimensions.size() << " indices at least " << dimCount+span);
            
            // Insert dimensions
            for(size_t i = 0; i < span; ++i) {
                evalDimensions.emplace_back(tensorObjectReadOnly->dimensions[dimCount+i]);
            }
        }
        return evalDimensions;
    }
    
    template<class tensor_type>
    std::vector<Index> IndexedTensorReadOnly<tensor_type>::get_assigned_indices() const {
        return get_assigned_indices(this->degree(), true);
    }
    
    template<class tensor_type>
    std::vector<Index> IndexedTensorReadOnly<tensor_type>::get_assigned_indices(const size_t _futureDegree, const bool _assignDimensions) const {
        std::vector<Index> assignedIndices;
        assignedIndices.reserve(indices.size());
        
        size_t dimensionCount = 0;
        for(const Index& idx : indices) {
            // We don't look at span zero indices
            if(idx.actual_span(_futureDegree) == 0) { continue; }
            
            // Fixed indices want special treatnment
            if(idx.fixed()) {
                if(_assignDimensions) {
                    assignedIndices.emplace_back(idx.valueId, 1, tensorObjectReadOnly->dimensions[dimensionCount++], Index::Flag::OPEN, Index::Flag::FIXED, false);
                } else {
                    assignedIndices.emplace_back(idx.valueId, 1, 0, Index::Flag::OPEN, Index::Flag::FIXED, false);
                    dimensionCount++;
                }
            } else {
                // Set span
                const size_t span = idx.actual_span(_futureDegree);
                
                // Calculate multDimension
                size_t multDimension = 1;
                if(_assignDimensions) {
                    REQUIRE(dimensionCount+span <= tensorObjectReadOnly->dimensions.size(), "Order determined by Indices is to large: " << dimensionCount+span << " > " << tensorObjectReadOnly->dimensions.size());
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
    
    #ifndef DISABLE_RUNTIME_CHECKS_
        /// Checks whether the indices are usefull in combination with the current dimensions
        template<class tensor_type>
        void IndexedTensorReadOnly<tensor_type>::check_indices(const bool _allowNonOpen) const {
            size_t dimensionCount = 0;
            for(const Index& idx : indices) {
                REQUIRE(_allowNonOpen || !idx.fixed(), "Fixed indices are not allowed here.");
                REQUIRE(_allowNonOpen || misc::count(indices, idx) == 1, "Traces are not allowed here.");
                REQUIRE(misc::count(indices, idx) <= 2, "An index must not appere more than twice!");
				dimensionCount += idx.actual_span(degree());
            }
            REQUIRE(dimensionCount >= degree(), "Order determined by Indices is to small. Order according to the indices " << dimensionCount << ", according to the tensor " << degree());
            REQUIRE(dimensionCount <= degree(), "Order determined by Indices is to large. Order according to the indices " << dimensionCount << ", according to the tensor " << degree());
        }
    #endif
    
    // IndexedTensorReadOnly may be instanciated as
    template class IndexedTensorReadOnly<Tensor>;
    template class IndexedTensorReadOnly<TensorNetwork>;
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - External functions - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    template<class tensor_type>
    value_t frob_norm(const IndexedTensorReadOnly<tensor_type>& _idxTensor) {
        return _idxTensor.tensorObjectReadOnly->frob_norm(); 
    }
    
    template value_t frob_norm<Tensor>(const IndexedTensorReadOnly<Tensor>& _idxTensor);
    template value_t frob_norm<TensorNetwork>(const IndexedTensorReadOnly<TensorNetwork>& _idxTensor);
    
    size_t get_eval_degree(const std::vector<Index>& _indices) {
        size_t degree = 0;
        for(const Index& idx : _indices) {
            REQUIRE(idx.flags[Index::Flag::ASSINGED], "Internal Error");
            if(!idx.fixed() && misc::count(_indices, idx) != 2) { degree += idx.span; }
        }
        return degree;
    }
    
}
