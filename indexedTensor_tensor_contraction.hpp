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

#include "xerus.h"

namespace xerus {

    #ifdef CHECK_
    /// Check if common and open indices defined by _lhs and _rhs coincide with the ones defined by _result
    void check_for_index_compatability(const AssignedIndices& _resultAssIndices, const AssignedIndices& _lhsAssIndices, const AssignedIndices& _rhsAssIndices) {
        LOG(ContractionDebug, "Checking input indices...");
        REQUIRE(_resultAssIndices.allIndicesOpen, "Result of contraction must not contain traces or fixed indices!");
        
        // Check that every index in result appears exactly once on lhs or rhs and that span and dimension conincide
        for(size_t i = 0; i < _resultAssIndices.numIndices; ++i) {
            size_t j;
            // Look for the index in lhs
            for(j = 0; j < _lhsAssIndices.numIndices; ++j) {
                if(_resultAssIndices.indices[i] == _lhsAssIndices.indices[j]) {
                    REQUIRE(_resultAssIndices.indices[i].span == _lhsAssIndices.indices[j].span, "Span of indices in result and lhs must conincide.");
                    REQUIRE(_resultAssIndices.indexDimensions[i] == _lhsAssIndices.indexDimensions[j], "Dimensions of indices in result and lhs must conincide.");
                    REQUIRE(_lhsAssIndices.indexOpen[j], "Index appearing in result of contraction must be open in lhs.");
                    REQUIRE(!contains(_rhsAssIndices.indices, _resultAssIndices.indices[i]), "Index appearing in result of contraction must not appear in lhs AND rhs.");
                    break;
                }
            }
            
            // If the index is not found in the lhs, search rhs
            if(j == _lhsAssIndices.numIndices) {
                for(j = 0; j < _rhsAssIndices.numIndices; ++j) {
                    if(_resultAssIndices.indices[i] == _rhsAssIndices.indices[j]) {
                        REQUIRE(_resultAssIndices.indices[i].span == _rhsAssIndices.indices[j].span, "Span of indices in result and rhs must conincide.");
                        REQUIRE(_resultAssIndices.indexDimensions[i] == _rhsAssIndices.indexDimensions[j], "Dimensions of indices in result and rhs must conincide.");
                        REQUIRE(_rhsAssIndices.indexOpen[j], "Index appearing in result of contraction must be open in rhs.");
                        break;
                    }
                }
                
                REQUIRE(j < _rhsAssIndices.numIndices, "Index appearing in the result of contraction must appear on either lhs or rhs side.");
            }
        }
        
        // Check that every index in lhs is either non-open, appears in rhs with right span and dimension, or is contained in result
        for(size_t i = 0; i < _lhsAssIndices.numIndices; ++i) {
            if(!_lhsAssIndices.indexOpen[i]) {
                REQUIRE(_lhsAssIndices.indices[i].is_fixed() || !contains(_rhsAssIndices.indices, _lhsAssIndices.indices[i]), "Index that part of a trace in lhs, must not appear in rhs");
                // It cannot be contained in result because of previous checks
                continue;
            }
            
            // Look for the index in rhs
            size_t j;
            for(j = 0; j < _rhsAssIndices.numIndices; ++j) {
                if(_lhsAssIndices.indices[i] == _rhsAssIndices.indices[j]) {
                    REQUIRE(_lhsAssIndices.indices[i].span == _rhsAssIndices.indices[j].span, "Span of indices in lhs and rhs of contraction must conincide.");
                    REQUIRE(_lhsAssIndices.indexDimensions[i] == _rhsAssIndices.indexDimensions[j], "Dimensions of indices in lhs and rhs of contraction  must conincide.");
                    REQUIRE(_rhsAssIndices.indexOpen[j], "Index appearing open in lhs of contraction must also be open in rhs.");
                    break;
                }
            }
            
            REQUIRE(j < _rhsAssIndices.numIndices || contains(_resultAssIndices.indices, _lhsAssIndices.indices[i]), "Index appearing in the lhs of contraction must appear on either rhs or result.");
        }
        
        // Check that every index in rhs is either non-open, or appears in lhs or result
        for(size_t i = 0; i < _rhsAssIndices.numIndices; ++i) {
            REQUIRE(!_rhsAssIndices.indexOpen[i] || contains(_lhsAssIndices.indices, _rhsAssIndices.indices[i]) || contains(_resultAssIndices.indices, _rhsAssIndices.indices[i]), "Every index appearing open in rhs of contraction must either appear in lhs or result.");
        }
        
        LOG(ContractionDebug, "Input indices look right.");
    }
    #endif
    
    /// Tests whether the indices of @a _other are seperated in @a _candidate and wether the order of the common indices coincides. Also gives @a _needsReshuffle = true, if @a _candidate contains non-open indices
    void test_seperation_and_order(bool& _needsReshuffle, bool& _isOrdered, const AssignedIndices& _candidate, const AssignedIndices& _other) {
        // If there are non-open indices the tensor has to be reshuffeled
        if(!_candidate.allIndicesOpen) {_needsReshuffle = true; _isOrdered = true; return; }
        
        // A set of zero indices is trivially seperated and ordered, the same applies if the result has no indices (and all indices are open).
        if(_candidate.indices.empty() || _other.indices.empty()) { _needsReshuffle = false; _isOrdered = true; return; }

        _needsReshuffle = false;
        _isOrdered = true;
        size_t resultStartIndex = 0;
        if(contains(_other.indices, _candidate.indices[0])) { // First index is common
            bool switched = false;
            while(_other.indices[resultStartIndex] != _candidate.indices[0]) { // Find the position of the first index in the other tensor
                ++resultStartIndex;
                REQUIRE(resultStartIndex < _other.indices.size(), "Internal Error");
            }
            for(size_t i=1; i < _candidate.numIndices; ++i) {
                if(!contains(_other.indices, _candidate.indices[i])) { 
                    switched = true; 
                } else if(switched) {
                    _needsReshuffle = true;
                    _isOrdered = true; // If the the tensor is reshuffeled it is also ordered afterwards
                    break; 
                } else if(resultStartIndex+i >= _other.numIndices || _other.indices[resultStartIndex+i] != _candidate.indices[i]) { 
                    _isOrdered = false;
                }
            }
        } else { // First index is un-common
            bool openIndexFound = false;
            for(size_t i=1; i < _candidate.numIndices; ++i){
                if(contains(_other.indices, _candidate.indices[i]) ) {
                    if(!openIndexFound) {
                        openIndexFound = true;
                        while(_other.indices[resultStartIndex] != _candidate.indices[i]) {  // Find this first common index in the result tensor
                            ++resultStartIndex;
                            REQUIRE(resultStartIndex < _other.indices.size(), "Internal Error");
                        }
                        resultStartIndex -= i; // We want to add i in the following
                    } else {
                        if(resultStartIndex+i >= _other.numIndices || _other.indices[resultStartIndex+i] != _candidate.indices[i]) { 
                            _isOrdered = false; 
                        }
                    }
                } else if(openIndexFound) {
                    _needsReshuffle = true;
                    _isOrdered = true; // If the the tensor is reordered it is also ordered afterwards
                    break;
                }
            }
        }
    }

    /// Check whether a perfect order is possible, i.e. the indices belonging to lhs and rhs are separated in the result
    bool perfect_order_possible(const AssignedIndices& _result, const AssignedIndices& _lhs, const AssignedIndices& _rhs) {
        LOG(ContractionDebug, "Checking whether a perfect order is possible...");
        if(!_lhs.indices.empty() && !_rhs.indices.empty() && !_result.indices.empty()) { // If any index set is empty perfect order is always possible
            bool switched = false;
            if(contains(_lhs.indices, _result.indices[0])){ // First index is contained in LHS
                for(size_t i=1; i< _result.numIndices; ++i) {
                    if(!contains(_lhs.indices, _result.indices[i])) {
                        switched = true;
                    } else if(switched) {
                        return false;
                    }
                }
            } else { // == first index is contained in RHS
                for(size_t i=1; i< _result.numIndices; ++i) {
                    if(!contains(_rhs.indices, _result.indices[i])) {
                        switched = true;
                    } else if(switched) {
                        return false;
                    }
                }
            }
        }
        return true;
    }

    /** Performs a contraction of all common indices of _lhs and _rhs and saves the resulting tensor in _result
    * precondition: compatible indices for _result = _lhs * _rhs; and correct dimensions in result tensorObject
    * postcondition: _result.tensorObject updated to correcly contain the product
    */
    void contract(const IndexedTensorWritable<Tensor>& _result, const IndexedTensorReadOnly<Tensor>& _lhs, const IndexedTensorReadOnly<Tensor>& _rhs) {
        // Get the assigned indices, we assume that result is already of right dimension
        const AssignedIndices lhsAssIndices = _lhs.assign_indices();
        const AssignedIndices rhsAssIndices = _rhs.assign_indices();
        const AssignedIndices resultAssIndices = _result.assign_indices();

        #ifdef CHECK_
            check_for_index_compatability(resultAssIndices, lhsAssIndices, rhsAssIndices);
        #endif
        
        if(_lhs.tensorObjectReadOnly->is_sparse() && _rhs.tensorObjectReadOnly->is_sparse() && _result.tensorObjectReadOnly->is_sparse()) {
            // We have to propagate the common factors
            _result.tensorObject->factor = _lhs.tensorObjectReadOnly->factor*_rhs.tensorObjectReadOnly->factor;
            
            std::vector<Index> lhsOpen, rhsOpen, common;
            
            for(const Index& idx : lhsAssIndices.indices) {
                if(contains(resultAssIndices.indices, idx)) {
                    lhsOpen.push_back(idx);
                } else if(contains(rhsAssIndices.indices, idx)) {
                    common.push_back(idx);
                }
            }
            for(const Index& idx : rhsAssIndices.indices) {
                if(contains(resultAssIndices.indices, idx)) {
                    rhsOpen.push_back(idx);
                }
            }
            
            CsUniquePtr lhsCS = to_cs_format(_lhs, lhsOpen, common);
            
            CsUniquePtr rhsCS = to_cs_format(_rhs, common, rhsOpen);
            
            CsUniquePtr resultCS = matrix_matrix_product(lhsCS, rhsCS);
            
            lhsOpen.insert(lhsOpen.end(), rhsOpen.begin(), rhsOpen.end());
            
            evaluate(_result, from_cs_format(resultCS, _result.get_evaluated_dimensions(lhsOpen))(lhsOpen));
            return;
        }
        
        // We don't want the result to have any factor that might interfere
        _result.tensorObject->factor = 1.0;
        
        // Check whether both sides are separated and whether their open indices are ordered
        LOG(ContractionDebug, "Checking Index uniqueness, seperation and order...");
        bool lhsNeedsReshuffle, lhsIsOrdered;
        bool rhsNeedsReshuffle, rhsIsOrdered;
        test_seperation_and_order(lhsNeedsReshuffle, lhsIsOrdered, lhsAssIndices, resultAssIndices);
        test_seperation_and_order(rhsNeedsReshuffle, rhsIsOrdered, rhsAssIndices, resultAssIndices);
        
        
        // Check compatibility
        bool bothSidesAreCompatible = true;
        if(!lhsNeedsReshuffle && !rhsNeedsReshuffle) { // If either needs a reshuffle, both sides will be compatible afterwards
            LOG(ContractionDebug, "Both sides are separated => Checking whether they are compatible");
            std::vector<Index>::const_iterator lhsItr = lhsAssIndices.indices.begin();
            std::vector<Index>::const_iterator rhsItr = rhsAssIndices.indices.begin();
            while (lhsItr != lhsAssIndices.indices.end() && contains(_result.indices, *lhsItr)) { lhsItr++; } // skip open indices (lhs)
            while (rhsItr != rhsAssIndices.indices.end() && contains(_result.indices, *rhsItr)) { rhsItr++; } // skip open indices (rhs)
            // Note that indices are separated at this point
            while (lhsItr != lhsAssIndices.indices.end() && rhsItr != rhsAssIndices.indices.end() && !contains(_result.indices, *lhsItr)) {
                if(*lhsItr != *rhsItr) {
                    LOG(ContractionDebug, "Common index order does not coincide => Not compatible.");
                    bothSidesAreCompatible = false;
                    break;
                }
                lhsItr++;
                rhsItr++;
            }
        }

        //Check whether it is cheaper to perform pre- or post-ordering
        bool reorderResult = true;
        if(perfect_order_possible(resultAssIndices, lhsAssIndices, rhsAssIndices)) {
            LOG(ContractionDebug, "Perfect order is possible. Calculating preorder costs");
            size_t preOrderCosts = 0;
            if(!lhsIsOrdered) { preOrderCosts += _lhs.tensorObjectReadOnly->size; }
            if(!rhsIsOrdered) { preOrderCosts += _rhs.tensorObjectReadOnly->size; }
            // if not compatible, the smaller one is no "cost" as it has to be reordered anyway
            if(!bothSidesAreCompatible) { preOrderCosts -= std::min(std::min(_lhs.tensorObjectReadOnly->size, _rhs.tensorObjectReadOnly->size), preOrderCosts); }
            
            
            if(preOrderCosts < _result.tensorObjectReadOnly->size) {
                LOG(ContractionDebug, "Cost for preorder are less than postOrder costs => Preorder factors to achieve perfect order.");
                if(!lhsIsOrdered) { lhsNeedsReshuffle = true; }
                if(!rhsIsOrdered) { rhsNeedsReshuffle = true; }
                reorderResult = false;
            }
        }

        if(!lhsNeedsReshuffle && !rhsNeedsReshuffle && !bothSidesAreCompatible) {
            LOG(ContractionDebug, "Common indices not compatible =>  reorder the smaller factor.");
            if(_lhs.tensorObjectReadOnly->size < _rhs.tensorObjectReadOnly->size) {
                lhsNeedsReshuffle = true;
            } else {
                rhsNeedsReshuffle = true;
            }
        }
        
        // Actually reorder LHS and RHS if necessary and create temporary result if reordering is necessary and calculate the dimensions of the matrification
        std::unique_ptr<IndexedTensor<Tensor>> lhsSaveSlot, rhsSaveSlot, workingResultSaveSlot;
        const IndexedTensorReadOnly<Tensor> *actualLhs, *actualRhs;
        const IndexedTensorWritable<Tensor> *workingResult;
        size_t leftDim = 1, rightDim=1, midDim=1;
        
        if(lhsNeedsReshuffle) {
            LOG(ContractionDebug, "Reordering LHS");
            std::vector<Index> newIndexOrder;
            
            // Add open indices in the order as they appear in the result
            for(size_t i = 0; i < resultAssIndices.numIndices; ++i) {
                if(contains(lhsAssIndices.indices, resultAssIndices.indices[i])) {
                    newIndexOrder.push_back(resultAssIndices.indices[i]);
                    leftDim *= resultAssIndices.indexDimensions[i];
                }
            }
            
            // Add common indices in the order as they appear in RHS
            for(size_t i = 0; i < rhsAssIndices.numIndices; ++i) {
                if(contains(lhsAssIndices.indices, rhsAssIndices.indices[i])) {
                    newIndexOrder.push_back(rhsAssIndices.indices[i]);
                    midDim *= rhsAssIndices.indexDimensions[i];
                }
            }
            
            lhsSaveSlot.reset(new IndexedTensor<Tensor>(new FullTensor(_lhs.get_evaluated_dimensions(newIndexOrder), DONT_SET_ZERO()), newIndexOrder, true));
            evaluate(*lhsSaveSlot, _lhs);
            actualLhs = lhsSaveSlot.get();
        } else {
            LOG(ContractionDebug, "LHS does not need reordering");
            for(size_t i = 0; i < lhsAssIndices.numIndices; ++i) {
                if(contains(resultAssIndices.indices, lhsAssIndices.indices[i])) {
                    leftDim *= lhsAssIndices.indexDimensions[i];
                } else {
                    midDim *= lhsAssIndices.indexDimensions[i];
                }
            }
            actualLhs = &_lhs;
        }
        
        if(rhsNeedsReshuffle) {
            LOG(ContractionDebug, "Reordering RHS");
            std::vector<Index> newIndexOrder;
            
            // Add common indices in the order of as they appear in the acutal LHS
            for(const Index& idx : actualLhs->indices) {
                if(contains(_rhs.indices, idx)) {
                    newIndexOrder.push_back(idx);
                }
            }
            
            // Add open indices in the order as they appear in the result
            for(size_t i = 0; i < resultAssIndices.numIndices; ++i) {
                if(contains(_rhs.indices, resultAssIndices.indices[i])) {
                    newIndexOrder.push_back(resultAssIndices.indices[i]);
                    rightDim *= resultAssIndices.indexDimensions[i];
                }
            }
            
            rhsSaveSlot.reset(new IndexedTensor<Tensor>(new FullTensor(_rhs.get_evaluated_dimensions(newIndexOrder), DONT_SET_ZERO()), newIndexOrder, true));
            evaluate(*rhsSaveSlot, _rhs);
            actualRhs = rhsSaveSlot.get();
        } else {
            LOG(ContractionDebug, "RHS does not need reordering");
            for(size_t i = 0; i < rhsAssIndices.numIndices; ++i) {
                if(contains(resultAssIndices.indices, rhsAssIndices.indices[i])) {
                    rightDim *= rhsAssIndices.indexDimensions[i];
                }
            }
            actualRhs = &_rhs;
        }
        
        if(reorderResult) {
            LOG(ContractionDebug, "Creating temporary result tensor");
            std::vector<Index> newIndexOrder;
            
            // Add LHS indices in the order as they appear
            for(const Index& idx : actualLhs->indices) {
                if(contains(_result.indices, idx)) {
                    newIndexOrder.push_back(idx);
                }
            }
            
            // Add RHS indices in the order of as they appear
            for(const Index& idx : actualRhs->indices) {
                if(contains(_result.indices, idx)) {
                    newIndexOrder.push_back(idx);
                }
            }
            
            workingResultSaveSlot.reset(new IndexedTensor<Tensor>(new FullTensor(_result.get_evaluated_dimensions(newIndexOrder), DONT_SET_ZERO()), newIndexOrder, true));
            workingResult = workingResultSaveSlot.get();
        } else {
            workingResult = &_result;
            static_cast<FullTensor*>(_result.tensorObject)->ensure_own_data_no_copy();
            
            // Check whether sides have to be switched to achieve perfect order
            if(!_result.indices.empty() && !contains(_lhs.indices, _result.indices[0])) {
                std::swap(actualLhs, actualRhs);
                std::swap(leftDim, rightDim);
            }
        }
        
        // Check if Matrices have to be transposed
        bool lhsTrans = !(actualLhs->indices.empty() || contains(resultAssIndices.indices, actualLhs->indices[0]));
        bool rhsTrans = !(actualRhs->indices.empty() || !contains(resultAssIndices.indices, actualRhs->indices[0]));
        
        LOG(ContractionDebug, "Performing Matrix multiplication of " << leftDim << "x" << midDim << " * " << midDim << "x" << rightDim << ".");
        
        const value_t commonFactor = _lhs.tensorObjectReadOnly->factor * _rhs.tensorObjectReadOnly->factor;
        blasWrapper::matrix_matrix_product(static_cast<FullTensor*>(workingResult->tensorObject)->data.get(), leftDim, rightDim, commonFactor, static_cast<const FullTensor*>(actualLhs->tensorObjectReadOnly)->data.get(), lhsTrans, midDim, static_cast<const FullTensor*>(actualRhs->tensorObjectReadOnly)->data.get(), rhsTrans);
        
        if(reorderResult) {
            LOG(ContractionDebug, "Reordering result");
            evaluate(_result, *workingResultSaveSlot);
        }
        
        return; // This was fun wasn't it?
    }

    IndexedTensorMoveable<Tensor> contract(const IndexedTensorReadOnly<Tensor>& _lhs, const IndexedTensorReadOnly<Tensor>& _rhs) {
        const std::vector<Index> lhsIndices = _lhs.get_assigned_indices();
        const std::vector<Index> rhsIndices = _rhs.get_assigned_indices();
        std::vector<Index> outIndices;
        std::vector<size_t> outDimensions;
                
        size_t dimensionCount = 0;
        for(const Index& idx : lhsIndices) {
            if(_lhs.is_open(idx) && !contains(rhsIndices, idx)) {
                outIndices.emplace_back(idx.valueId, idx.inverseSpan ? _lhs.degree()-idx.span : idx.span, false);
                for(size_t i=0; i < outIndices.back().span; ++i) {
                    outDimensions.push_back(_lhs.tensorObjectReadOnly->dimensions[dimensionCount++]);
                }
            } else {
                dimensionCount += idx.inverseSpan ? _lhs.degree()-idx.span : idx.span;
            }
        }
        
        dimensionCount = 0;
        for(const Index& idx : rhsIndices) {
            if(_rhs.is_open(idx) && !contains(lhsIndices, idx)) {
                outIndices.emplace_back(idx.valueId, idx.inverseSpan ? _rhs.degree()-idx.span : idx.span, false);
                for(size_t i=0; i < outIndices.back().span; ++i) {
                    outDimensions.push_back(_rhs.tensorObjectReadOnly->dimensions[dimensionCount++]);
                }
            } else {
                dimensionCount += idx.inverseSpan ? _rhs.degree()-idx.span : idx.span;
            }
        }
        
        
        // TODO More intelligence when to use a sparse result
        Tensor* resultTensor;
        if(_lhs.tensorObjectReadOnly->is_sparse() && _rhs.tensorObjectReadOnly->is_sparse()) {
            resultTensor = new SparseTensor(outDimensions);
        } else {
            resultTensor = new FullTensor(outDimensions, DONT_SET_ZERO());
        }
        IndexedTensorMoveable<Tensor> result(resultTensor, outIndices);
        contract(result, _lhs, _rhs);
        return result;
    }
}
