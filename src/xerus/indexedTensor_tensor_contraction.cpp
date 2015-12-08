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
 * @brief Implementation of the (indexed) FullTensor contraction.
 */

#include <xerus/indexedTensor_tensor_operators.h>
#include <xerus/basic.h>
#include <xerus/index.h>
#include <xerus/misc/check.h>
#include <xerus/fullTensor.h>
#include <xerus/sparseTensor.h>
#include <xerus/cs_wrapper.h>
#include <xerus/sparseTimesFullContraction.h>
#include <xerus/blasLapackWrapper.h>
#include <memory>

namespace xerus {
    
    #ifndef DISABLE_RUNTIME_CHECKS_
    /// Check if common and open indices defined by _lhs and _rhs coincide with the ones defined by _result
    void check_input_validity(const std::vector<Index>& _resultIndices, const std::vector<Index>& _lhsIndices, const std::vector<Index>& _rhsIndices) {
        LOG(ContractionDebug, "Checking input indices...");
        REQUIRE(Index::all_open(_resultIndices), "Result of contraction must not contain traces or fixed indices!");
        
        // Check that every index in result appears exactly once on lhs or rhs and that span and dimension conincide
        for(size_t i = 0; i < _resultIndices.size(); ++i) {
            size_t j;
            // Look for the index in lhs
            for(j = 0; j < _lhsIndices.size(); ++j) {
                if(_resultIndices[i] == _lhsIndices[j]) {
                    REQUIRE(_resultIndices[i].span == _lhsIndices[j].span, "Span of indices in result and lhs must conincide.");
                    REQUIRE(_resultIndices[i].dimension() == _lhsIndices[j].dimension(), "Dimensions of indices in result and lhs must conincide.");
                    REQUIRE(_lhsIndices[j].open(), "Index appearing in result of contraction must be open in lhs.");
                    REQUIRE(!misc::contains(_rhsIndices, _resultIndices[i]), "Index appearing in result of contraction must not appear in lhs AND rhs.");
                    break;
                }
            }
            
            // If the index is not found in the lhs, search rhs
            if(j == _lhsIndices.size()) {
                for(j = 0; j < _rhsIndices.size(); ++j) {
                    if(_resultIndices[i] == _rhsIndices[j]) {
                        REQUIRE(_resultIndices[i].span == _rhsIndices[j].span, "Span of indices in result and rhs must conincide.");
                        REQUIRE(_resultIndices[i].dimension() == _rhsIndices[j].dimension(), "Dimensions of indices in result and rhs must conincide.");
                        REQUIRE(_rhsIndices[j].open(), "Index appearing in result of contraction must be open in rhs.");
                        break;
                    }
                }
                
                REQUIRE(j < _rhsIndices.size(), "Index appearing in the result of contraction must appear on either lhs or rhs side.");
            }
        }
        
        // Check that every index in lhs is either non-open, appears in rhs with right span and dimension, or is contained in result
        for(size_t i = 0; i < _lhsIndices.size(); ++i) {
            if(!_lhsIndices[i].open()) {
                REQUIRE(_lhsIndices[i].fixed() || !misc::contains(_rhsIndices, _lhsIndices[i]), "Index that part of a trace in lhs, must not appear in rhs");
                // It cannot be contained in result because of previous checks
                continue;
            }
            
            // Look for the index in rhs
            size_t j;
            for(j = 0; j < _rhsIndices.size(); ++j) {
                if(_lhsIndices[i] == _rhsIndices[j]) {
                    REQUIRE(_lhsIndices[i].span == _rhsIndices[j].span, "Span of indices in lhs and rhs of contraction must conincide.");
                    REQUIRE(_lhsIndices[i].dimension() == _rhsIndices[j].dimension(), "Dimensions of indices in lhs and rhs of contraction must conincide.");
                    REQUIRE(_rhsIndices[j].open(), "Index appearing open in lhs of contraction must also be open in rhs.");
                    break;
                }
            }
            
            REQUIRE(j < _rhsIndices.size() || misc::contains(_resultIndices, _lhsIndices[i]), "Index appearing in the lhs of contraction must appear on either rhs or result.");
        }
        
        // Check that every index in rhs is either non-open, or appears in lhs or result
        for(size_t i = 0; i < _rhsIndices.size(); ++i) {
            REQUIRE(!_rhsIndices[i].open() || misc::contains(_lhsIndices, _rhsIndices[i]) || misc::contains(_resultIndices, _rhsIndices[i]), "Every index appearing open in rhs of contraction must either appear in lhs or result.");
        }
        
        LOG(ContractionDebug, "Input indices look alright.");
    }
    #endif
    
    /// Tests whether the indices of @a _other are seperated in @a _candidate and whether the order of the common indices coincides. Also gives @a _needsReshuffle = true, if @a _candidate contains non-open indices
    void test_seperation_and_order(bool& _needsReshuffle, bool& _isOrdered, const std::vector<Index>& _candidate, const std::vector<Index>& _other) {
        // If there are non-open indices the tensor has to be reshuffeled
        if(!Index::all_open(_candidate)) {_needsReshuffle = true; _isOrdered = true; return; }
        
        // A set of zero indices is trivially seperated and ordered, the same applies if the result has no indices (and all indices are open).
        if(_candidate.empty() || _other.empty()) { _needsReshuffle = false; _isOrdered = true; return; }

        _needsReshuffle = false;
        _isOrdered = true;
        size_t resultStartIndex = 0;
        if(misc::contains(_other, _candidate[0])) { // First index is common
            bool switched = false;
            
            // Find the position of the first index in the other tensor
            while(_other[resultStartIndex] != _candidate[0]) { ++resultStartIndex; }
            
            // Iterate over all _candidate indices and check for correct order
            for(size_t i=1; i < _candidate.size(); ++i) {
                if(!misc::contains(_other, _candidate[i])) {
                    switched = true; 
                } else if(switched) { // We are not seperated -> everything is lost.
                    _needsReshuffle = true; 
                    _isOrdered = true; // If the the tensor is reshuffeled it is also ordered afterwards
                    break; 
                } else if(resultStartIndex+i >= _other.size() || _other[resultStartIndex+i] != _candidate[i]) { // The order does not coincide
                    _isOrdered = false;
                }
            }
        } else { // First index is un-common
            bool openIndexFound = false;
            for(size_t i=1; i < _candidate.size(); ++i){
                if(misc::contains(_other, _candidate[i]) ) {
                    if(!openIndexFound) {
                        openIndexFound = true;
                        while(_other[resultStartIndex] != _candidate[i]) {  // Find this first common index in the result tensor
                            ++resultStartIndex;
                            REQUIRE(resultStartIndex < _other.size(), "Internal Error");
                        }
                        resultStartIndex -= i; // We want to add i in the following
                    } else {
                        if(resultStartIndex+i >= _other.size() || _other[resultStartIndex+i] != _candidate[i]) { 
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
    bool check_compatibility(const std::vector<Index>& _resultAssIndices, const std::vector<Index>& _lhsAssIndices, const std::vector<Index>& _rhsAssIndices) {
        LOG(ContractionDebug, "Checking whether they are compatible");
        bool bothSidesAreCompatible = true;
        std::vector<Index>::const_iterator lhsItr = _lhsAssIndices.begin();
        std::vector<Index>::const_iterator rhsItr = _rhsAssIndices.begin();
        while (lhsItr != _lhsAssIndices.end() && misc::contains(_resultAssIndices, *lhsItr)) { lhsItr++; } // skip open indices (lhs)
        while (rhsItr != _rhsAssIndices.end() && misc::contains(_resultAssIndices, *rhsItr)) { rhsItr++; } // skip open indices (rhs)
        // Note that indices are separated at this point
        while (lhsItr != _lhsAssIndices.end() && rhsItr != _rhsAssIndices.end() && !misc::contains(_resultAssIndices, *lhsItr)) {
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
    bool perfect_order_possible(const std::vector<Index>& _result, const std::vector<Index>& _lhs, const std::vector<Index>& _rhs) {
        LOG(ContractionDebug, "Checking whether a perfect order is possible...");
        if(!_lhs.empty() && !_rhs.empty() && !_result.empty()) { // If any index set is empty perfect order is always possible
            bool switched = false;
            if(misc::contains(_lhs, _result[0])){ // First index is contained in LHS
                for(size_t i=1; i< _result.size(); ++i) {
                    if(!misc::contains(_lhs, _result[i])) {
                        switched = true;
                    } else if(switched) {
                        return false;
                    }
                }
            } else { // == first index is contained in RHS
                for(size_t i=1; i< _result.size(); ++i) {
                    if(!misc::contains(_rhs, _result[i])) {
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
    void contract(const IndexedTensorWritable<Tensor>& _result, const IndexedTensorReadOnly<Tensor>& _lhs, const std::vector<Index>& _lhsIndices, const IndexedTensorReadOnly<Tensor>& _rhs, const std::vector<Index>& _rhsIndices) {
        // Get the assigned indices for result (we assume that result is already of the right dimensions)
        const std::vector<Index> resultIndices = _result.get_assigned_indices();

        #ifndef DISABLE_RUNTIME_CHECKS_
            check_input_validity(resultIndices, _lhsIndices, _rhsIndices);
        #endif
        
        // Check whether both sides are separated and whether their open indices are ordered
        LOG(ContractionDebug, "Checking Index uniqueness, seperation and order...");
        bool lhsNeedsReshuffle, lhsIsOrdered;
        bool rhsNeedsReshuffle, rhsIsOrdered;
        test_seperation_and_order(lhsNeedsReshuffle, lhsIsOrdered, _lhsIndices, resultIndices);
        test_seperation_and_order(rhsNeedsReshuffle, rhsIsOrdered, _rhsIndices, resultIndices);
        
        // Check compatibility, if either needs a reshuffle, both sides will be compatible afterwards
        const bool bothSidesAreCompatible = lhsNeedsReshuffle || rhsNeedsReshuffle || check_compatibility(resultIndices, _lhsIndices, _rhsIndices);
        
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
        if(perfect_order_possible(resultIndices, _lhsIndices, _rhsIndices)) {
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
		lhsOpenIndices.reserve(_lhsIndices.size());
		rhsOpenIndices.reserve(_rhsIndices.size());
		commonIndices.reserve(std::min(_lhsIndices.size(), _rhsIndices.size()));
        size_t leftDim = 1, midDim = 1, rightDim = 1;
        
        // Set lhs and common indices
        if(lhsNeedsReshuffle) {
            // Add open indices in the order as they appear in the result
            for(const Index& idx : resultIndices) {
                if(misc::contains(_lhsIndices, idx)) {
                    lhsOpenIndices.push_back(idx);
                    leftDim *= idx.dimension();
                }
            }
            
            // Add common indices in the order as they appear in RHS
            for(const Index& idx : _rhsIndices) {
                if(misc::contains(_lhsIndices, idx)) {
                    commonIndices.push_back(idx);
                    midDim *= idx.dimension();
                }
            }
        } else {
            // Add indices in the order they also appear in LHS
            for(const Index& idx : _lhsIndices) {
                if(misc::contains(resultIndices, idx)) {
                    lhsOpenIndices.push_back(idx);
                    leftDim *= idx.dimension();
                } else if(misc::contains(_rhsIndices, idx)) {
                    commonIndices.push_back(idx);
                    midDim *= idx.dimension();
                }
            }
        }
        
        // Set lhs and common indices
        if(rhsNeedsReshuffle) {
            // Add open indices in the order as they appear in the result
            for(const Index& idx : resultIndices) {
                if(misc::contains(_rhsIndices, idx)) {
                    rhsOpenIndices.push_back(idx);
                    rightDim *= idx.dimension();
                }
            }
        } else {
            // Add indices in the order they also appear in RHS
            for(const Index& idx : _rhsIndices) {
                if(misc::contains(resultIndices, idx)) {
                    rhsOpenIndices.push_back(idx);
                    rightDim *= idx.dimension();
                }
            }
        }
        
        std::vector<Index> workingResultIndices(lhsOpenIndices);
        workingResultIndices.insert(workingResultIndices.end(), rhsOpenIndices.begin(), rhsOpenIndices.end());
        
        // - - - - - - - - - - - - - - - - - - - - - - - SparseTensor * SparseTensor ==> SparseTensor - - - - - - - - - - - - - - - - - - - - - - -
        if(_lhs.tensorObjectReadOnly->is_sparse() && _rhs.tensorObjectReadOnly->is_sparse() && _result.tensorObjectReadOnly->is_sparse()) {
            // We have to propagate the common factors
            _result.tensorObject->factor = _lhs.tensorObjectReadOnly->factor*_rhs.tensorObjectReadOnly->factor;
            
            CsUniquePtr lhsCS = to_cs_format(_lhs, lhsOpenIndices, commonIndices);
            CsUniquePtr rhsCS = to_cs_format(_rhs, commonIndices, rhsOpenIndices);
            CsUniquePtr resultCS = matrix_matrix_product(lhsCS, rhsCS);
            
            evaluate(_result, from_cs_format(resultCS, _result.get_evaluated_dimensions(workingResultIndices))(workingResultIndices));
            
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
                if(!_result.indices.empty() && !misc::contains(_lhsIndices, resultIndices[0])) {
                    std::swap(actualLhs, actualRhs);
                    std::swap(leftDim, rightDim);
                }
            }
            
            // Check if Matrices have to be transposed
            const bool lhsTrans = !(actualLhs->indices.empty() || misc::contains(resultIndices, actualLhs->indices[0]));
            const bool rhsTrans = !(actualRhs->indices.empty() || !misc::contains(resultIndices, actualRhs->indices[0]));
            
            LOG(ContractionDebug, "Performing Matrix multiplication of " << leftDim << "x" << midDim << " * " << midDim << "x" << rightDim << ".");
            
            const value_t commonFactor = _lhs.tensorObjectReadOnly->factor * _rhs.tensorObjectReadOnly->factor;
            
            const bool lhsSparse = actualLhs->tensorObjectReadOnly->is_sparse();
            const bool rhsSparse = actualRhs->tensorObjectReadOnly->is_sparse();
            const bool resultSparse = workingResult->tensorObjectReadOnly->is_sparse();
            
            // Select actual case
            if(!lhsSparse && !rhsSparse && !resultSparse) { // Full * Full => Full
                blasWrapper::matrix_matrix_product(static_cast<FullTensor*>(workingResult->tensorObject)->unsanitized_data_pointer(), leftDim, rightDim, commonFactor, static_cast<const FullTensor*>(actualLhs->tensorObjectReadOnly)->unsanitized_data_pointer(), lhsTrans, midDim, static_cast<const FullTensor*>(actualRhs->tensorObjectReadOnly)->unsanitized_data_pointer(), rhsTrans);
            } else if(lhsSparse && !rhsSparse && !resultSparse) { // Sparse * Full => Full
                matrix_matrix_product(static_cast<FullTensor*>(workingResult->tensorObject)->unsanitized_data_pointer(), leftDim, rightDim, commonFactor, *static_cast<const SparseTensor*>(actualLhs->tensorObjectReadOnly)->sparseData.get(), lhsTrans, midDim, static_cast<const FullTensor*>(actualRhs->tensorObjectReadOnly)->unsanitized_data_pointer(), rhsTrans);
            } else if(!lhsSparse && rhsSparse && !resultSparse) { // Full * Sparse => Full
                matrix_matrix_product(static_cast<FullTensor*>(workingResult->tensorObject)->unsanitized_data_pointer(), leftDim, rightDim, commonFactor, static_cast<const FullTensor*>(actualLhs->tensorObjectReadOnly)->unsanitized_data_pointer(), lhsTrans, midDim, *static_cast<const SparseTensor*>(actualRhs->tensorObjectReadOnly)->sparseData.get(), rhsTrans);
            } else if(lhsSparse && !rhsSparse && resultSparse) { // Sparse * Full => Sparse
                matrix_matrix_product(*static_cast<SparseTensor*>(workingResult->tensorObject)->sparseData.get(), leftDim, rightDim, commonFactor, *static_cast<const SparseTensor*>(actualLhs->tensorObjectReadOnly)->sparseData.get(), lhsTrans, midDim, static_cast<const FullTensor*>(actualRhs->tensorObjectReadOnly)->unsanitized_data_pointer(), rhsTrans);
            } else if(!lhsSparse && rhsSparse && resultSparse) { // Full * Sparse => Sparse
                matrix_matrix_product(*static_cast<SparseTensor*>(workingResult->tensorObject)->sparseData.get(), leftDim, rightDim, commonFactor, static_cast<const FullTensor*>(actualLhs->tensorObjectReadOnly)->unsanitized_data_pointer(), lhsTrans, midDim, *static_cast<const SparseTensor*>(actualRhs->tensorObjectReadOnly)->sparseData.get(), rhsTrans);
            } else {
                LOG(fatal, "Invalid combiantion of sparse/dense tensors in contraction");
            }
            
            if(reorderResult) {
                LOG(ContractionDebug, "Reordering result");
                evaluate(_result, *workingResultSaveSlot);
            }
        }
    }
    
    void contract(const IndexedTensorWritable<Tensor>& _result, const IndexedTensorReadOnly<Tensor>& _lhs,  const IndexedTensorReadOnly<Tensor>& _rhs) {
        contract(_result, _lhs, _lhs.get_assigned_indices(), _rhs, _rhs.get_assigned_indices());
    }

    
    IndexedTensorMoveable<Tensor> contract(const IndexedTensorReadOnly<Tensor>& _lhs, const IndexedTensorReadOnly<Tensor>& _rhs) {
        const std::vector<Index> lhsIndices = _lhs.get_assigned_indices();
        const std::vector<Index> rhsIndices = _rhs.get_assigned_indices();
        std::vector<Index> outIndices;
		outIndices.reserve(_lhs.indices.size() + _rhs.indices.size());
        std::vector<size_t> outDimensions;
		outDimensions.reserve(_lhs.degree()+_rhs.degree());
                
        size_t lhsOpenDim = 1, rhsOpenDim = 1;
        size_t dimensionCount = 0;
        for(const Index& idx : lhsIndices) {
            if(idx.open() && !misc::contains(rhsIndices, idx)) {
                outIndices.emplace_back(idx);
                for(size_t i=0; i < idx.span; ++i) {
                    lhsOpenDim *= _lhs.tensorObjectReadOnly->dimensions[dimensionCount];
                    outDimensions.push_back(_lhs.tensorObjectReadOnly->dimensions[dimensionCount]);
                    dimensionCount++;
                }
            } else {
                dimensionCount += idx.span;
            }
        }
        
        dimensionCount = 0;
        for(const Index& idx : rhsIndices) {
            if(idx.open() && !misc::contains(lhsIndices, idx)) {
                outIndices.emplace_back(idx);
                for(size_t i=0; i < idx.span; ++i) {
                    rhsOpenDim *= _rhs.tensorObjectReadOnly->dimensions[dimensionCount];
                    outDimensions.push_back(_rhs.tensorObjectReadOnly->dimensions[dimensionCount]);
                    dimensionCount++;
                }
            } else {
                dimensionCount += idx.span;
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
        contract(result, _lhs, lhsIndices, _rhs, rhsIndices);
        return result;
    }
}
