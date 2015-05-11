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

    #ifndef DISABLE_RUNTIME_CHECKS_
    /// Check if common and open indices defined by _lhs and _rhs coincide with the ones defined by _result
    void check_input_validity(const AssignedIndices& _resultAssIndices, const AssignedIndices& _lhsAssIndices, const AssignedIndices& _rhsAssIndices) {
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
                    REQUIRE(_lhsAssIndices.indices[j].open(), "Index appearing in result of contraction must be open in lhs.");
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
                        REQUIRE(_rhsAssIndices.indices[j].open(), "Index appearing in result of contraction must be open in rhs.");
                        break;
                    }
                }
                
                REQUIRE(j < _rhsAssIndices.numIndices, "Index appearing in the result of contraction must appear on either lhs or rhs side.");
            }
        }
        
        // Check that every index in lhs is either non-open, appears in rhs with right span and dimension, or is contained in result
        for(size_t i = 0; i < _lhsAssIndices.numIndices; ++i) {
            if(!_lhsAssIndices.indices[i].open()) {
                REQUIRE(_lhsAssIndices.indices[i].fixed() || !contains(_rhsAssIndices.indices, _lhsAssIndices.indices[i]), "Index that part of a trace in lhs, must not appear in rhs");
                // It cannot be contained in result because of previous checks
                continue;
            }
            
            // Look for the index in rhs
            size_t j;
            for(j = 0; j < _rhsAssIndices.numIndices; ++j) {
                if(_lhsAssIndices.indices[i] == _rhsAssIndices.indices[j]) {
                    REQUIRE(_lhsAssIndices.indices[i].span == _rhsAssIndices.indices[j].span, "Span of indices in lhs and rhs of contraction must conincide.");
                    REQUIRE(_lhsAssIndices.indexDimensions[i] == _rhsAssIndices.indexDimensions[j], "Dimensions of indices in lhs and rhs of contraction must conincide.");
                    REQUIRE(_rhsAssIndices.indices[j].open(), "Index appearing open in lhs of contraction must also be open in rhs.");
                    break;
                }
            }
            
            REQUIRE(j < _rhsAssIndices.numIndices || contains(_resultAssIndices.indices, _lhsAssIndices.indices[i]), "Index appearing in the lhs of contraction must appear on either rhs or result.");
        }
        
        // Check that every index in rhs is either non-open, or appears in lhs or result
        for(size_t i = 0; i < _rhsAssIndices.numIndices; ++i) {
            REQUIRE(!_rhsAssIndices.indices[i].open() || contains(_lhsAssIndices.indices, _rhsAssIndices.indices[i]) || contains(_resultAssIndices.indices, _rhsAssIndices.indices[i]), "Every index appearing open in rhs of contraction must either appear in lhs or result.");
        }
        
        LOG(ContractionDebug, "Input indices look alright.");
    }
    #endif
    
    /// Tests whether the indices of @a _other are seperated in @a _candidate and whether the order of the common indices coincides. Also gives @a _needsReshuffle = true, if @a _candidate contains non-open indices
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
            
            // Find the position of the first index in the other tensor
            while(_other.indices[resultStartIndex] != _candidate.indices[0]) { ++resultStartIndex; }
            
            // Iterate over all _candidate indices and check for correct order
            for(size_t i=1; i < _candidate.numIndices; ++i) {
                if(!contains(_other.indices, _candidate.indices[i])) {
                    switched = true; 
                } else if(switched) { // We are not seperated -> everything is lost.
                    _needsReshuffle = true; 
                    _isOrdered = true; // If the the tensor is reshuffeled it is also ordered afterwards
                    break; 
                } else if(resultStartIndex+i >= _other.numIndices || _other.indices[resultStartIndex+i] != _candidate.indices[i]) { // The order does not coincide
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
                } else if(openIndexFound) { // We are not seperated -> everything is lost.
                    _needsReshuffle = true;
                    _isOrdered = true; // If the the tensor is reordered it is also ordered afterwards
                    break;
                }
            }
        }
    }
    
    /// Check whether lhsIndices and rhsIndices are compatible in the sense that the common indices are in the same order
    bool check_compatibility(const AssignedIndices& _resultAssIndices, const AssignedIndices& _lhsAssIndices, const AssignedIndices& _rhsAssIndices) {
        bool bothSidesAreCompatible = true;
        LOG(ContractionDebug, "Both sides are separated => Checking whether they are compatible");
        std::vector<Index>::const_iterator lhsItr = _lhsAssIndices.indices.begin();
        std::vector<Index>::const_iterator rhsItr = _rhsAssIndices.indices.begin();
        while (lhsItr != _lhsAssIndices.indices.end() && contains(_resultAssIndices.indices, *lhsItr)) { lhsItr++; } // skip open indices (lhs)
        while (rhsItr != _rhsAssIndices.indices.end() && contains(_resultAssIndices.indices, *rhsItr)) { rhsItr++; } // skip open indices (rhs)
        // Note that indices are separated at this point
        while (lhsItr != _lhsAssIndices.indices.end() && rhsItr != _rhsAssIndices.indices.end() && !contains(_resultAssIndices.indices, *lhsItr)) {
            if(*lhsItr != *rhsItr) {
                LOG(ContractionDebug, "Common index order does not coincide => Not compatible.");
                bothSidesAreCompatible = false;
                break;
            }
            lhsItr++;
            rhsItr++;
        }
        return bothSidesAreCompatible;
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

    /** Performs the contraction of all common indices of _lhs and _rhs and saves the resulting tensor in _result
    * precondition: compatible indices for _result = _lhs * _rhs; and correct dimensions in result tensorObject
    * postcondition: _result.tensorObject updated to correcly contain the product
    */
    void contract(const IndexedTensorWritable<Tensor>& _result, const IndexedTensorReadOnly<Tensor>& _lhs, const IndexedTensorReadOnly<Tensor>& _rhs) {
        // Get the assigned indices (we assume that result is already of the right dimensions)
        const AssignedIndices lhsAssIndices = _lhs.assign_indices(); // TODO get rid of AssignedIndices
        const AssignedIndices rhsAssIndices = _rhs.assign_indices();
        const AssignedIndices resultAssIndices = _result.assign_indices();
        
        const std::vector<Index> lhsIndices = _lhs.get_assigned_indices();
        const std::vector<Index> rhsIndices = _rhs.get_assigned_indices();
        const std::vector<Index> resultIndices = _result.get_assigned_indices();

        #ifndef DISABLE_RUNTIME_CHECKS_
            check_input_validity(resultAssIndices, lhsAssIndices, rhsAssIndices);
        #endif
        
        // Check whether both sides are separated and whether their open indices are ordered
        LOG(ContractionDebug, "Checking Index uniqueness, seperation and order...");
        bool lhsNeedsReshuffle, lhsIsOrdered;
        bool rhsNeedsReshuffle, rhsIsOrdered;
        test_seperation_and_order(lhsNeedsReshuffle, lhsIsOrdered, lhsAssIndices, resultAssIndices);
        test_seperation_and_order(rhsNeedsReshuffle, rhsIsOrdered, rhsAssIndices, resultAssIndices);
        
        // Check compatibility, if either needs a reshuffle, both sides will be compatible afterwards
        const bool bothSidesAreCompatible = lhsNeedsReshuffle || rhsNeedsReshuffle || check_compatibility(resultAssIndices, lhsAssIndices, rhsAssIndices);
        
        // If both sides ar enot compatible we have to reorder the smaller one
        if(!bothSidesAreCompatible) {
            LOG(ContractionDebug, "Common indices not compatible =>  reorder the smaller factor.");
            if(_lhs.tensorObjectReadOnly->reorder_costs() < _rhs.tensorObjectReadOnly->reorder_costs()) {
                lhsNeedsReshuffle = true;
                lhsIsOrdered = true;
            } else {
                rhsNeedsReshuffle = true;
                rhsIsOrdered = true;
            }
        }
        
        // Check whether it is cheaper to perform pre- or post-ordering
        bool reorderResult = true;
        if(perfect_order_possible(resultAssIndices, lhsAssIndices, rhsAssIndices)) {
            LOG(ContractionDebug, "Perfect order is possible. Calculating preorder costs");
            size_t preOrderCosts = 0;
            if(!lhsIsOrdered) { preOrderCosts += _lhs.tensorObjectReadOnly->reorder_costs(); }
            if(!rhsIsOrdered) { preOrderCosts += _rhs.tensorObjectReadOnly->reorder_costs(); }
            
            if(preOrderCosts < _result.tensorObjectReadOnly->size) { // TODO approximate the cost for reordering result (non-zero entries not yet known)
                LOG(ContractionDebug, "Cost for preorder are less than postOrder costs => Preorder factors to achieve perfect order.");
                if(!lhsIsOrdered) { lhsNeedsReshuffle = true; }
                if(!rhsIsOrdered) { rhsNeedsReshuffle = true; }
                reorderResult = false;
            }
        }

        // Construct the vectors for reordering
        std::vector<Index> lhsOpenIndices, rhsOpenIndices, commonIndices;
        size_t leftDim = 1, midDim = 1, rightDim = 1;
        
        // Set lhs and common indices
        if(lhsNeedsReshuffle) {
            // Add open indices in the order as they appear in the result
            for(const Index& idx : resultIndices) {
                if(contains(lhsIndices, idx)) {
                    lhsOpenIndices.push_back(idx);
                    leftDim *= idx.dimension;
                }
            }
            
            // Add common indices in the order as they appear in RHS
            for(const Index& idx : rhsIndices) {
                if(contains(lhsIndices, idx)) {
                    commonIndices.push_back(idx);
                    midDim *= idx.dimension;
                }
            }
        } else {
            // Add indices in the order they also appear in LHS
            for(const Index& idx : lhsIndices) {
                if(contains(resultIndices, idx)) {
                    lhsOpenIndices.push_back(idx);
                    leftDim *= idx.dimension;
                } else if(contains(rhsIndices, idx)) {
                    commonIndices.push_back(idx);
                    midDim *= idx.dimension;
                }
            }
        }
        
        // Set lhs and common indices
        if(rhsNeedsReshuffle) {
            // Add open indices in the order as they appear in the result
            for(const Index& idx : resultIndices) {
                if(contains(rhsIndices, idx)) {
                    rhsOpenIndices.push_back(idx);
                    rightDim *= idx.dimension;
                }
            }
        } else {
            // Add indices in the order they also appear in RHS
            for(const Index& idx : rhsIndices) {
                if(contains(resultIndices, idx)) {
                    rhsOpenIndices.push_back(idx);
                    rightDim *= idx.dimension;
                }
            }
        }
        
        std::vector<Index> workingResultIndices(lhsOpenIndices);
        workingResultIndices.insert(workingResultIndices.end(), rhsOpenIndices.begin(), rhsOpenIndices.end());
        
        // - - - - - - - - - - - - - - - - - - - - - - - SparseTensor * SparseTensor ==> FullTensor - - - - - - - - - - - - - - - - - - - - - - -
        if(_lhs.tensorObjectReadOnly->is_sparse() && _rhs.tensorObjectReadOnly->is_sparse() && _result.tensorObjectReadOnly->is_sparse()) {
            // We have to propagate the common factors
            _result.tensorObject->factor = _lhs.tensorObjectReadOnly->factor*_rhs.tensorObjectReadOnly->factor;
            
            CsUniquePtr lhsCS = to_cs_format(_lhs, lhsOpenIndices, commonIndices);
            
            CsUniquePtr rhsCS = to_cs_format(_rhs, commonIndices, rhsOpenIndices);
            
            CsUniquePtr resultCS = matrix_matrix_product(lhsCS, rhsCS);
            
            evaluate(_result, from_cs_format(resultCS, _result.get_evaluated_dimensions(workingResultIndices))(workingResultIndices)); // TODO this evaluation should be part of from_cs_format
            
        // - - - - - - - - - - - - - - - - - - - - - - - Other - - - - - - - - - - - - - - - - - - - - - - - 
        } else { 
            // We don't want the result to have any factor since we take the factors into account during the calculation.
            _result.tensorObject->factor = 1.0;
            
            // Actually reorder LHS and RHS if necessary and create temporary result if reordering is necessary and calculate the dimensions of the matrification
            std::unique_ptr<IndexedTensor<Tensor>> lhsSaveSlot, rhsSaveSlot, workingResultSaveSlot;
            const IndexedTensorReadOnly<Tensor> *actualLhs, *actualRhs;
            const IndexedTensorWritable<Tensor> *workingResult;
            
            if(lhsNeedsReshuffle) {
                LOG(ContractionDebug, "Reordering LHS");
                
                // Add the common indices to the open ones
                lhsOpenIndices.insert(lhsOpenIndices.end(), commonIndices.begin(), commonIndices.end());
                
                lhsSaveSlot.reset(new IndexedTensor<Tensor>(_lhs.tensorObjectReadOnly->construct_new(_lhs.get_evaluated_dimensions(lhsOpenIndices), DONT_SET_ZERO()), std::move(lhsOpenIndices), true));
                evaluate(*lhsSaveSlot, _lhs);
                actualLhs = lhsSaveSlot.get();
            } else {
                actualLhs = &_lhs;
            }
            
            if(rhsNeedsReshuffle) {
                LOG(ContractionDebug, "Reordering RHS");
                
                // Add the open indices to the common ones
                commonIndices.insert(commonIndices.end(), rhsOpenIndices.begin(), rhsOpenIndices.end());
               
                
                rhsSaveSlot.reset(new IndexedTensor<Tensor>(_rhs.tensorObjectReadOnly->construct_new(_rhs.get_evaluated_dimensions(commonIndices), DONT_SET_ZERO()), std::move(commonIndices), true));
                evaluate(*rhsSaveSlot, _rhs);
                actualRhs = rhsSaveSlot.get();
            } else {
                actualRhs = &_rhs;
            }
            
            if(reorderResult) {
                LOG(ContractionDebug, "Creating temporary result tensor");
                
                workingResultSaveSlot.reset(new IndexedTensor<Tensor>(_result.tensorObjectReadOnly->construct_new(_result.get_evaluated_dimensions(workingResultIndices), DONT_SET_ZERO()), std::move(workingResultIndices), true));
                workingResult = workingResultSaveSlot.get();
            } else {
                workingResult = &_result;
                _result.tensorObject->ensure_own_data_no_copy();
                
                // Check whether both sides have to be swaped to achieve perfect order
                if(!_result.indices.empty() && !contains(lhsAssIndices.indices, _result.indices[0])) {
                    std::swap(actualLhs, actualRhs);
                    std::swap(leftDim, rightDim);
                }
            }
            
            // Check if Matrices have to be transposed
            const bool lhsTrans = !(actualLhs->indices.empty() || contains(resultAssIndices.indices, actualLhs->indices[0]));
            const bool rhsTrans = !(actualRhs->indices.empty() || !contains(resultAssIndices.indices, actualRhs->indices[0]));
            
            LOG(ContractionDebug, "Performing Matrix multiplication of " << leftDim << "x" << midDim << " * " << midDim << "x" << rightDim << ".");
            
            const value_t commonFactor = _lhs.tensorObjectReadOnly->factor * _rhs.tensorObjectReadOnly->factor;
            
            const bool lhsSparse = _lhs.tensorObjectReadOnly->is_sparse();
            const bool rhsSparse = _rhs.tensorObjectReadOnly->is_sparse();
            const bool resultSparse = _result.tensorObjectReadOnly->is_sparse();
            
            LOG(bla, "WHAT: " << lhsSparse << rhsSparse << resultSparse);
            // Select actual case
            if(!lhsSparse && !rhsSparse && !resultSparse) { // Full * Full => Full
                blasWrapper::matrix_matrix_product(static_cast<FullTensor*>(workingResult->tensorObject)->data.get(), leftDim, rightDim, commonFactor, static_cast<const FullTensor*>(actualLhs->tensorObjectReadOnly)->data.get(), lhsTrans, midDim, static_cast<const FullTensor*>(actualRhs->tensorObjectReadOnly)->data.get(), rhsTrans);
            } else if(lhsSparse && !rhsSparse && !resultSparse) { // Sparse * Full => Full
                LOG(bla, "WHAT");
                matrix_matrix_product(static_cast<FullTensor*>(workingResult->tensorObject)->data.get(), leftDim, rightDim, commonFactor, *static_cast<const SparseTensor*>(actualLhs->tensorObjectReadOnly)->entries.get(), lhsTrans, midDim, static_cast<const FullTensor*>(actualRhs->tensorObjectReadOnly)->data.get(), rhsTrans);
            } else if(!lhsSparse && rhsSparse && !resultSparse) { // Full * Sparse => Full
                LOG(bla, "WHAT");
                matrix_matrix_product(static_cast<FullTensor*>(workingResult->tensorObject)->data.get(), leftDim, rightDim, commonFactor, static_cast<const FullTensor*>(actualLhs->tensorObjectReadOnly)->data.get(), lhsTrans, midDim, *static_cast<const SparseTensor*>(actualRhs->tensorObjectReadOnly)->entries.get(), rhsTrans);
            } else if(lhsSparse && !rhsSparse && resultSparse) { // Sparse * Full => Sparse
//                 matrix_matrix_product(static_cast<SparseTensor*>(workingResult->tensorObject)->entries.get(), leftDim, rightDim, commonFactor, *static_cast<const SparseTensor*>(actualLhs->tensorObjectReadOnly)->entries.get(), lhsTrans, midDim, static_cast<const FullTensor*>(actualRhs->tensorObjectReadOnly)->data.get(), rhsTrans);
            } else if(!lhsSparse && rhsSparse && resultSparse) { // Full * Sparse => Sparse
//                 matrix_matrix_product(static_cast<SparseTensor*>(workingResult->tensorObject)->entries.get(), leftDim, rightDim, commonFactor, static_cast<const FullTensor*>(actualLhs->tensorObjectReadOnly)->data.get(), lhsTrans, midDim, *static_cast<const SparseTensor*>(actualRhs->tensorObjectReadOnly)->entries.get(), rhsTrans);
            } // TODO
            
            if(reorderResult) {
                LOG(ContractionDebug, "Reordering result");
                evaluate(_result, *workingResultSaveSlot);
            }
        }
    }

    
    IndexedTensorMoveable<Tensor> contract(const IndexedTensorReadOnly<Tensor>& _lhs, const IndexedTensorReadOnly<Tensor>& _rhs) {
        const std::vector<Index> lhsIndices = _lhs.get_assigned_indices();
        const std::vector<Index> rhsIndices = _rhs.get_assigned_indices();
        std::vector<Index> outIndices;
        std::vector<size_t> outDimensions;
                
        size_t lhsOpenDim = 1, rhsOpenDim = 1;
        size_t dimensionCount = 0;
        for(const Index& idx : lhsIndices) {
            if(_lhs.is_open(idx) && !contains(rhsIndices, idx)) {
                outIndices.emplace_back(idx.valueId, idx.flags[Index::Flag::INVERSE_SPAN] ? _lhs.degree()-idx.span : idx.span);
                for(size_t i=0; i < outIndices.back().span; ++i) {
                    lhsOpenDim *= _lhs.tensorObjectReadOnly->dimensions[dimensionCount];
                    outDimensions.push_back(_lhs.tensorObjectReadOnly->dimensions[dimensionCount]);
                    dimensionCount++;
                }
            } else {
                dimensionCount += idx.flags[Index::Flag::INVERSE_SPAN] ? _lhs.degree()-idx.span : idx.span;
            }
        }
        
        dimensionCount = 0;
        for(const Index& idx : rhsIndices) {
            if(_rhs.is_open(idx) && !contains(lhsIndices, idx)) {
                outIndices.emplace_back(idx.valueId, idx.flags[Index::Flag::INVERSE_SPAN] ? _rhs.degree()-idx.span : idx.span);
                for(size_t i=0; i < outIndices.back().span; ++i) {
                    rhsOpenDim *= _rhs.tensorObjectReadOnly->dimensions[dimensionCount];
                    outDimensions.push_back(_rhs.tensorObjectReadOnly->dimensions[dimensionCount]);
                    dimensionCount++;
                }
            } else {
                dimensionCount += idx.flags[Index::Flag::INVERSE_SPAN] ? _rhs.degree()-idx.span : idx.span;
            }
        }
        
        Tensor* resultTensor;
        if(   (_lhs.tensorObjectReadOnly->is_sparse() && _rhs.tensorObjectReadOnly->is_sparse())
           || (_lhs.tensorObjectReadOnly->is_sparse() && _lhs.tensorObjectReadOnly->count_non_zero_entries() < lhsOpenDim)
           || (_rhs.tensorObjectReadOnly->is_sparse() && _rhs.tensorObjectReadOnly->count_non_zero_entries() < rhsOpenDim)
        ) {
            resultTensor = new SparseTensor(outDimensions);
        } else {
            resultTensor = new FullTensor(outDimensions, DONT_SET_ZERO());
        }
        
        IndexedTensorMoveable<Tensor> result(resultTensor, outIndices);
        contract(result, _lhs, _rhs);
        return result;
    }
}
