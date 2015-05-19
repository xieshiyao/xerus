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

#include <xerus/indexedTensor_tensor_operators.h>
#include <xerus/basic.h>
#include <xerus/index.h>
#include <xerus/misc/test.h>
#include <xerus/fullTensor.h>
#include <xerus/sparseTensor.h>
#include <memory>
#include <xerus/misc/selectedFunctions.h>
#include <xerus/misc/missingFunctions.h>

namespace xerus {

    /** "increases" the imaginary indices that lead to the current pointer position of oldPosition to match the index i
    * precondition: oldPosition corresponding to index (i-1)
    * postcondition: oldPosition corresponding to index i
    */
    void increase_indices(const size_t _i, const value_t*& _oldPosition, const size_t _numIndices, const size_t* const _steps, const size_t* const _multDimensions) {
        size_t index = _numIndices-1;
        _oldPosition += _steps[index];
        size_t multStep = _multDimensions[index];
        while(_i%multStep == 0) {
            _oldPosition -= _multDimensions[index]*_steps[index]; // "reset" current index to 0
            --index; // Advance to next index
            _oldPosition += _steps[index]; // increase next index
            multStep *= _multDimensions[index]; // next stepSize
        }
    }

    void sum_traces(   value_t* const _newPosition,
                                const value_t* _oldPosition,
                                const size_t* const _doubleSteps,
                                const size_t* const _doubleMultDimensions,
                                const size_t _numDoubleIndexPairs,
                                const size_t _numSummations) {
        *_newPosition = *_oldPosition;
        for(size_t k = 1; k < _numSummations; ++k) {
            increase_indices(k, _oldPosition, _numDoubleIndexPairs, _doubleSteps, _doubleMultDimensions);
            *_newPosition += *_oldPosition;
        }
    }
    void sum_traces(   value_t* const _newPosition,
                                const value_t* _oldPosition,
                                const size_t* const _doubleSteps,
                                const size_t* const _doubleMultDimensions,
                                const size_t _numDoubleIndexPairs,
                                const size_t _numSummations,
                                const size_t _orderedIndicesMultDim ) {
        array_copy(_newPosition, _oldPosition, _orderedIndicesMultDim);
        for(size_t k = 1; k < _numSummations; ++k) {
            increase_indices(k, _oldPosition, _numDoubleIndexPairs, _doubleSteps, _doubleMultDimensions);
            array_add(_newPosition, 1.0, _oldPosition, _orderedIndicesMultDim);
        }
    }

    size_t get_position(const std::pair<size_t, value_t>& _entry,
                        const size_t* const _baseIndexDim,
                        const size_t* const _divisors,
                        const size_t* const _attributes,
                        const size_t numIndices ) {
        size_t position = 0;
        for(size_t i = 0; i < numIndices; ++i) {
            position += ((_entry.first/_divisors[i])%_baseIndexDim[i])*_attributes[i];
        }
        return position;
    }

    bool check_position(size_t & _position,
                        const std::pair<size_t, value_t>& _entry,
                        const size_t* const _baseIndexDim,
                        const size_t* const _divisors,
                        const size_t* const _attributes,
                        const bool* const _fixedFlags,
                        const bool* const _traceFlags,
                        const size_t numIndices ) {
        _position = 0;
        // Process each index
        for(size_t i = 0; i < numIndices; ++i) {
            const size_t indexPosition = (_entry.first/_divisors[i])%_baseIndexDim[i];
            
            // The position is only changed if the index is not fixed...
            if(_fixedFlags[i]) {
                if(indexPosition != _attributes[i]) {
                    return false; // The indexPosition differs from the fixed position => Do not add this value
                }
            // ... and not part of a trace.
            } else if(_traceFlags[i]) {
                if(indexPosition != (_entry.first/_divisors[_attributes[i]])%_baseIndexDim[_attributes[i]]) {
                    return false; // This entry is not on the diagonal of the trace => Do not add this value
                }
            } else {
                _position += indexPosition*_attributes[i];
            }
        }
        return true;
    }

    std::unique_ptr<const size_t[]> get_dimension_array(const std::vector<Index> _indices) {
        if(_indices.size() == 0) {
            return std::unique_ptr<const size_t[]>();
        } else {
            size_t* stepSizes = new size_t[_indices.size()];
            for(size_t i = 0; i < _indices.size(); ++i) {
                stepSizes[i] = _indices[i].dimension();
            }
            return std::unique_ptr<const size_t[]>(stepSizes);
        }
    }
    
    std::unique_ptr<const size_t[]> get_step_sizes(const std::vector<Index> _indices) {
        if(_indices.size() == 0) {
            return std::unique_ptr<const size_t[]>();
        } else {
            size_t* stepSizes = new size_t[_indices.size()];
            stepSizes[_indices.size()-1] = 1;
            for(size_t i = _indices.size()-1; i > 0; --i) {
                stepSizes[i-1] = stepSizes[i]*_indices[i].dimension();
            }
            return std::unique_ptr<const size_t[]>(stepSizes);
        }
    }


    void evaluate(const IndexedTensorWritable<Tensor>& _out, const IndexedTensorReadOnly<Tensor>& _base) {
        // Get the assigned indices
        const std::vector<Index> baseIndices = _base.get_assigned_indices();
        const std::vector<Index> outIndices = _out.get_assigned_indices();
        
        // Extract base index dimensions
        const std::unique_ptr<const size_t[]> baseIndexDimensions = get_dimension_array(baseIndices);
        
        #ifndef DISABLE_RUNTIME_CHECKS_ // Performe complete check whether the input is valid
            REQUIRE(_out.tensorObjectReadOnly != _base.tensorObjectReadOnly, "Target of evaluation must not conincide with base!");
            REQUIRE(!_out.tensorObjectReadOnly->is_sparse() || _base.tensorObjectReadOnly->is_sparse(), "Evaluation of FullTensor to SparseTensor not implemented and probably not useful.");
            
            // Check base indices
            for(size_t i = 0; i < baseIndices.size(); ++i) {
                // If the index is fixed we don't expect to find it anywhere
                if(baseIndices[i].fixed()) {
                    REQUIRE(baseIndices[i].span == 1, "Fixed indices must have span one. Indices are: " << baseIndices << ", total should be " << baseIndices.size() << ". The problem is: " << baseIndices[i] << " -- " << baseIndices[i].fixed());
                    continue;
                }
                
                // Try to find index in _out
                size_t j = 0;
                while(j < outIndices.size() && baseIndices[i] != outIndices[j]) { ++j; }
                if(j < outIndices.size()) {
                    REQUIRE(baseIndices[i].dimension() == outIndices[j].dimension(), "The indexDimensions in the target and base of evaluation must coincide. Here " << baseIndices[i].dimension() << "!=" << outIndices[j].dimension() << ". For index " << baseIndices[i] << " == " << outIndices[j]);
                    REQUIRE(baseIndices[i].span == outIndices[j].span, "The indexSpans in the target and base of evaluation must coincide.");
                    REQUIRE(baseIndices[i].open(), "Indices appearing in the target of evaluation must not be part of a trace nor be fixed. Base: " << baseIndices << " Out: " << outIndices);
                    continue;
                }
                
                // Try to find index a second time in base
                j = 0;
                while(j < baseIndices.size() && (i == j || baseIndices[i] != baseIndices[j])) { ++j; }
                REQUIRE(j < baseIndices.size(), "All indices of evalutation base must either be fixed, appear in the target or be part of a trace. Base: " << baseIndices << " Out: " << outIndices);
                REQUIRE(count(baseIndices, baseIndices[i]) == 2, "Indices must appear at most two times. Base: " << baseIndices << " Out: " << outIndices);
                REQUIRE(baseIndices[i].dimension() == baseIndices[j].dimension(), "The indexDimensions of two traced indices must conince.");
                REQUIRE(baseIndices[i].span == 1 && baseIndices[j].span == 1, "The indexSpans of traced indices must be one (It is ambigious what a trace of span 2 indices is meant to be).");
            }
            
            // Check out indices
            for(size_t i = 0; i < outIndices.size(); ++i) {
                REQUIRE(outIndices[i].open(),  "Traces and fixed indices are not allowed in the target of evaluation. Base: " << baseIndices << " Out: " << outIndices);
                REQUIRE(count(baseIndices, outIndices[i]) == 1, "Every index of the target must appear exactly once in the base of evaluation. Base: " << baseIndices << " Out: " << outIndices);
            }
        #endif
        
        // If there is no index reshuffling, we can simplify a lot
        if(baseIndices == outIndices) {
            if(!_out.tensorObjectReadOnly->is_sparse() && !_base.tensorObjectReadOnly->is_sparse()) { // Full => Full
                _out.tensorObject->factor = _base.tensorObjectReadOnly->factor;
                static_cast<FullTensor*>(_out.tensorObject)->data = static_cast<const FullTensor*>(_base.tensorObjectReadOnly)->data;
                
            } else if(_out.tensorObjectReadOnly->is_sparse() && _base.tensorObjectReadOnly->is_sparse()) { // Sparse => Sparse
                _out.tensorObject->factor = _base.tensorObjectReadOnly->factor;
                static_cast<SparseTensor*>(_out.tensorObject)->entries = static_cast<const SparseTensor*>(_base.tensorObjectReadOnly)->entries;
            
            } else if(_base.tensorObjectReadOnly->is_sparse()) { // Sparse => Full
                FullTensor& outTensor = *static_cast<FullTensor*>(_out.tensorObject);
                value_t* const outData = outTensor.data.get();
                outTensor.ensure_own_data_no_copy();
                array_set_zero(outData, outTensor.size);
                for(const std::pair<size_t, value_t>& entry : *static_cast<const SparseTensor*>(_base.tensorObjectReadOnly)->entries) {
                    outData[entry.first] = entry.second;
                }
            }
            return; // We are finished here
        }
        
        // We need the step sizes of the base indices
        const std::unique_ptr<const size_t[]> baseIndexStepSizes = get_step_sizes(baseIndices);
        
        
        // In every case we have to ensure that _out has its own data, since we gonna rewrite it.
        _out.tensorObject->ensure_own_data_no_copy();
        
        //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - Full => Full   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if(!_out.tensorObjectReadOnly->is_sparse() && !_base.tensorObjectReadOnly->is_sparse()) {
            // Extract out index dimensions
            const std::unique_ptr<const size_t[]> outIndexDimensions = get_dimension_array(outIndices);
        
            // Propagate the constant factor, since we won't apply it for FullTensors
            _out.tensorObject->factor = _base.tensorObjectReadOnly->factor;
            
            // Count how many indices in the back are already ordered (We know that base has at least as many indices as out)
            size_t numOrderedIndices;
            for(numOrderedIndices = 0; numOrderedIndices < outIndices.size() && baseIndices[baseIndices.size()-1-numOrderedIndices] == outIndices[outIndices.size()-1-numOrderedIndices]; ++numOrderedIndices) { }
            
            const size_t orderedIndexDim = baseIndexStepSizes[baseIndices.size()-numOrderedIndices-1];
            
            VLA(size_t[outIndices.size()-numOrderedIndices], stepSizes); // How much we have to move in base when the i-th index of _out increases
            
            size_t fixedIndexOffset = 0; // Fixed offset in base due to fixed indices
            
            std::vector<size_t> traceStepSizes; // How much we have to move in base if the index of the i-th trace is increased
            std::vector<size_t> traceDimensions; // The dimensions of the traces
            size_t totalTraceDim = 1; // Total number of summantions, i.e. Product of all trace dimensions
            
            // Calculate stepSizes for our tensor. We will march in steps of orderedIndicesMultDim in _out.data
            for(size_t i = 0; i < baseIndices.size()-numOrderedIndices; ++i) {
                // First try to find the index in _out (if it is not contained it must be fixed or part of a trace)
                size_t outPos = 0;
                while(outPos < outIndices.size() && baseIndices[i] != outIndices[outPos]) { ++outPos; }
                
                if(outPos < outIndices.size()) { // If we found it we are basically finished
                    stepSizes[outPos] = baseIndexStepSizes[i];
                } else if(baseIndices[i].fixed()) { // One reason for an index to not be in _out is to be fixed
                    fixedIndexOffset += size_t(baseIndices[i].valueId)*baseIndexStepSizes[i];
                } else { // If the Index is not fixed, then it has to be part of a trace
                    for(size_t j = i+1; j < baseIndices.size()-numOrderedIndices; ++j) {
                        if(baseIndices[i] == baseIndices[j]) {
                            traceStepSizes.emplace_back(baseIndexStepSizes[i]+baseIndexStepSizes[j]);
                            traceDimensions.emplace_back(baseIndices[i].dimension());
                            totalTraceDim *= baseIndices[i].dimension();
                            break;
                        }
                    }
                }
            }
            
            // Get pointers to the data
            const value_t* oldPosition = static_cast<const FullTensor*>(_base.tensorObjectReadOnly)->data.get()+fixedIndexOffset;
            value_t* const newPosition = static_cast<FullTensor*>(_out.tensorObject)->data.get();
            
            if(orderedIndexDim == 1) { // We have to copy/add single entries
                if(totalTraceDim == 1) { // We don't need to sum any traces
                    newPosition[0] = *oldPosition;
                    for(size_t i = 1; i < _out.tensorObject->size; ++i) {
                        increase_indices(i, oldPosition, outIndices.size()-numOrderedIndices, stepSizes, outIndexDimensions.get());
                        newPosition[i] = *oldPosition;
                    }
                } else { // We have to add traces
                    sum_traces(newPosition, oldPosition, traceStepSizes.data(), traceDimensions.data(), traceDimensions.size(), totalTraceDim);
                    for(size_t i = 1; i < _out.tensorObject->size; ++i) {
                        increase_indices(i, oldPosition, outIndices.size()-numOrderedIndices, stepSizes, outIndexDimensions.get());
                        sum_traces(newPosition+i, oldPosition, traceStepSizes.data(), traceDimensions.data(), traceDimensions.size(), totalTraceDim);
                    }
                }
            } else { // We can copy/add larger blocks
                if(totalTraceDim == 1) { // We don't need to sum any traces
                    array_copy(newPosition, oldPosition, orderedIndexDim);
                    for(size_t i = 1; i < _out.tensorObject->size/orderedIndexDim; ++i) {
                        increase_indices(i, oldPosition, outIndices.size()-numOrderedIndices, stepSizes, outIndexDimensions.get());
                        array_copy(newPosition + i*orderedIndexDim, oldPosition, orderedIndexDim);
                    }
                } else { // We have to add traces
                    sum_traces(newPosition, oldPosition, traceStepSizes.data(), traceDimensions.data(), traceDimensions.size(), totalTraceDim, orderedIndexDim);
                    for(size_t i = 1; i < _out.tensorObject->size/orderedIndexDim; ++i) {
                        increase_indices(i, oldPosition, outIndices.size()-numOrderedIndices, stepSizes, outIndexDimensions.get());
                        sum_traces(newPosition + i*orderedIndexDim, oldPosition, traceStepSizes.data(), traceDimensions.data(), traceDimensions.size(), totalTraceDim, orderedIndexDim);
                    }
                }
            }
        }
        //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - Sparse => Both  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        else if(_base.tensorObjectReadOnly->is_sparse()) {
            VLA(bool[baseIndices.size()]  , fixedFlags); // Flag for each index indicating whether the index is fixed
            VLA(bool[baseIndices.size()]  , traceFlags); // Flag for each index indicating whether the index is part of a trace
            VLA(size_t[baseIndices.size()], attributes);  // Either the factor in _out, the value of an fixed index or the position of the other part of a trace
            bool peacefullIndices = true;
            
            // Process base indices
            const std::unique_ptr<const size_t[]> outIndexStepSizes = get_step_sizes(outIndices);
            for(size_t i = 0; i < baseIndices.size(); ++i) {
                // First try to find the index in out
                size_t outPos = 0;
                while(outPos < outIndices.size() && outIndices[outPos] != baseIndices[i]) { ++outPos; }
                if(outPos < outIndices.size()) {
                    fixedFlags[i] = false;
                    traceFlags[i] = false;
                    attributes[i] = outIndexStepSizes[outPos];
                    continue;
                }
                
                // Check whether the index is fixed
                if(baseIndices[i].fixed()) {
                    fixedFlags[i] = true;
                    traceFlags[i] = false;
                    attributes[i] = (size_t) baseIndices[i].valueId;
                    peacefullIndices = false;
                    continue;
                }
                
                // If the index is not in out and not fixed then it has to be part of a trace
                fixedFlags[i] = false;
                traceFlags[i] = true;
                peacefullIndices = false;
                for(attributes[i] = 0; baseIndices[i] != baseIndices[attributes[i]] || attributes[i] == i; ++attributes[i]) { }
            }
            
            // Get the entries and the factor of our base tensor
            const std::map<size_t, value_t>& baseEntries = *static_cast<const SparseTensor*>(_base.tensorObjectReadOnly)->entries;
            const value_t factor = _base.tensorObjectReadOnly->factor;
            
            
            // Check whether _out is sparse
            if(_out.tensorObjectReadOnly->is_sparse()) {
                // Ensure that _out is empty
                std::map<size_t, value_t>& outEntries = *static_cast<SparseTensor*>(_out.tensorObject)->entries;
                outEntries.clear();
                
                if(peacefullIndices) {
                    for(const std::pair<size_t, value_t>& entry : baseEntries) {
                        outEntries.insert(std::pair<size_t, value_t>(get_position(entry, baseIndexDimensions.get(), baseIndexStepSizes.get(), attributes, baseIndices.size()), factor*entry.second));
                    }
                } else {
                    for(const std::pair<size_t, value_t>& entry : baseEntries) {
                        size_t newPosition;
                        if(check_position(newPosition, entry, baseIndexDimensions.get(), baseIndexStepSizes.get(), attributes, fixedFlags, traceFlags, baseIndices.size())) {
                            outEntries[newPosition] += factor*entry.second;
                        }
                    }
                }
            } else {
                // Ensure that _out is empty
                value_t* const dataPointer = static_cast<FullTensor*>(_out.tensorObject)->data.get();
                array_set_zero(dataPointer, _out.tensorObject->size);
                
                if(peacefullIndices) {
                    for(const std::pair<size_t, value_t>& entry : baseEntries) {
                        dataPointer[get_position(entry, baseIndexDimensions.get(), baseIndexStepSizes.get(), attributes, baseIndices.size())] = factor*entry.second;
                    }
                } else {
                    for(const std::pair<size_t, value_t>& entry : baseEntries) {
                        size_t newPosition;
                        if(check_position(newPosition, entry, baseIndexDimensions.get(), baseIndexStepSizes.get(), attributes, fixedFlags, traceFlags, baseIndices.size())) {
                            dataPointer[newPosition] += factor*entry.second;
                        }
                    }
                }
            }
        }
    }
}
