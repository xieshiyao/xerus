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

#include "../../include/xerus.h"

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

    std::unique_ptr<const size_t[]> get_step_sizes(const AssignedIndices _assIndices) {
        if(_assIndices.numIndices == 0) {
            return std::unique_ptr<const size_t[]>();
        } else {
            size_t* stepSizes = new size_t[_assIndices.numIndices];
            stepSizes[_assIndices.numIndices-1] = 1;
            for(size_t i = _assIndices.numIndices-1; i > 0; --i) {
                stepSizes[i-1] = stepSizes[i]*_assIndices.indexDimensions[i];
            }
            return std::unique_ptr<const size_t[]>(stepSizes);
        }
    }


    void evaluate(const IndexedTensorWritable<Tensor>& _out, const IndexedTensorReadOnly<Tensor>& _base) {
        // Get the assigned indices
        const AssignedIndices baseAssIdx = _base.assign_indices();
        const AssignedIndices outAssIdx = _out.assign_indices();
        
        #ifndef DISABLE_RUNTIME_CHECKS_ // Performe complete check whether the input is valid
            REQUIRE(_out.tensorObjectReadOnly != _base.tensorObjectReadOnly, "Target of evaluation must not conincide with base!");
            REQUIRE(!_out.tensorObjectReadOnly->is_sparse() || _base.tensorObjectReadOnly->is_sparse(), "Evaluation of SparseTensor to FullTensor not implemented and probably not useful.");
            
            // Check base indices
            for(size_t i = 0; i < baseAssIdx.numIndices; ++i) {
                // If the index is fixed we don't expect to find it anywhere
                if(baseAssIdx.indices[i].fixed()) {
                    REQUIRE(baseAssIdx.indices[i].span == 1, "Fixed indices must have span one. Indices are: " << baseAssIdx.indices << ", total should be " << baseAssIdx.numIndices << ". The problem is: " << baseAssIdx.indices[i] << " -- " << baseAssIdx.indices[i].fixed());
                    continue;
                }
                
                // Try to find index in _out
                size_t j = 0;
                while(j < outAssIdx.numIndices && baseAssIdx.indices[i] != outAssIdx.indices[j]) { ++j; }
                if(j < outAssIdx.numIndices) {
                    REQUIRE(baseAssIdx.indexDimensions[i] == outAssIdx.indexDimensions[j], "The indexDimensions in the target and base of evaluation must coincide. Here " << baseAssIdx.indexDimensions[i] << "!=" << outAssIdx.indexDimensions[j] << ". For index " << baseAssIdx.indices[i] << " == " << outAssIdx.indices[j]);
                    REQUIRE(baseAssIdx.indices[i].span == outAssIdx.indices[j].span, "The indexSpans in the target and base of evaluation must coincide.");
                    REQUIRE(baseAssIdx.indices[i].open(), "Indices appearing in the target of evaluation must not be part of a trace nor be fixed. Base: " << baseAssIdx.indices << " Out: " << outAssIdx.indices);
                    continue;
                }
                
                // Try to find index a second time in base
                j = 0;
                while(j < baseAssIdx.numIndices && (i == j || baseAssIdx.indices[i] != baseAssIdx.indices[j])) { ++j; }
                REQUIRE(j < baseAssIdx.numIndices, "All indices of evalutation base must either be fixed, appear in the target or be part of a trace. Base: " << baseAssIdx.indices << " Out: " << outAssIdx.indices);
                REQUIRE(count(baseAssIdx.indices, baseAssIdx.indices[i]) == 2, "Indices must appear at most two times. Base: " << baseAssIdx.indices << " Out: " << outAssIdx.indices);
                REQUIRE(baseAssIdx.indexDimensions[i] == baseAssIdx.indexDimensions[j], "The indexDimensions of two traced indices must conince.");
                REQUIRE(baseAssIdx.indices[i].span == 1 && baseAssIdx.indices[j].span == 1, "The indexSpans of traced indices must be one (It is ambigious what a trace of span 2 indices is meant to be).");
            }
            
            // Check out indices
            for(size_t i = 0; i < outAssIdx.numIndices; ++i) {
                REQUIRE(outAssIdx.indices[i].open(),  "Traces and fixed indices are not allowed in the target of evaluation. Base: " << baseAssIdx.indices << " Out: " << outAssIdx.indices);
                REQUIRE(count(baseAssIdx.indices, outAssIdx.indices[i]) == 1, "Every index of the target must appear exactly once in the base of evaluation. Base: " << baseAssIdx.indices << " Out: " << outAssIdx.indices);
            }
        #endif
        
        // If there is no index reshuffling, we can simplify a lot
        if(baseAssIdx.indices == outAssIdx.indices) {
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
        const std::unique_ptr<const size_t[]> baseIndexStepSizes = get_step_sizes(baseAssIdx);
        
        
        // In every case we have to ensure that _out has its own data, since we gonna rewrite it.
        _out.tensorObject->ensure_own_data_no_copy();
        
        //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - Full => Full   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if(!_out.tensorObjectReadOnly->is_sparse() && !_base.tensorObjectReadOnly->is_sparse()) {
            // Propagate the constant factor, since we won't apply it for FullTensors
            _out.tensorObject->factor = _base.tensorObjectReadOnly->factor;
            
            // Count how many indices in the back are already ordered (We know that base has at least as many indices as out)
            size_t numOrderedIndices;
            for(numOrderedIndices = 0; numOrderedIndices < outAssIdx.numIndices && baseAssIdx.indices[baseAssIdx.numIndices-1-numOrderedIndices] == outAssIdx.indices[outAssIdx.numIndices-1-numOrderedIndices]; ++numOrderedIndices) { }
            
            const size_t orderedIndexDim = baseIndexStepSizes[baseAssIdx.numIndices-numOrderedIndices-1];
            
            VLA(size_t[outAssIdx.numIndices-numOrderedIndices], stepSizes); // How much we have to move in base when the i-th index of _out increases
            
            size_t fixedIndexOffset = 0; // Fixed offset in base due to fixed indices
            
            std::vector<size_t> traceStepSizes; // How much we have to move in base if the index of the i-th trace is increased
            std::vector<size_t> traceDimensions; // The dimensions of the traces
            size_t totalTraceDim = 1; // Total number of summantions, i.e. Product of all trace dimensions
            
            // Calculate stepSizes for our tensor. We will march in steps of orderedIndicesMultDim in _out.data
            for(size_t i = 0; i < baseAssIdx.numIndices-numOrderedIndices; ++i) {
                // First try to find the index in _out (if it is not contained it must be fixed or part of a trace)
                size_t outPos = 0;
                while(outPos < outAssIdx.numIndices && baseAssIdx.indices[i] != outAssIdx.indices[outPos]) { ++outPos; }
                
                if(outPos < outAssIdx.numIndices) { // If we found it we are basically finished
                    stepSizes[outPos] = baseIndexStepSizes[i];
                } else if(baseAssIdx.indices[i].fixed()) { // One reason for an index to not be in _out is to be fixed
                    fixedIndexOffset += size_t(baseAssIdx.indices[i].valueId)*baseIndexStepSizes[i];
                } else { // If the Index is not fixed, then it has to be part of a trace
                    for(size_t j = i+1; j < baseAssIdx.numIndices-numOrderedIndices; ++j) {
                        if(baseAssIdx.indices[i] == baseAssIdx.indices[j]) {
                            traceStepSizes.emplace_back(baseIndexStepSizes[i]+baseIndexStepSizes[j]);
                            traceDimensions.emplace_back(baseAssIdx.indexDimensions[i]);
                            totalTraceDim *= baseAssIdx.indexDimensions[i];
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
                        increase_indices(i, oldPosition, outAssIdx.numIndices-numOrderedIndices, stepSizes, outAssIdx.indexDimensions.data());
                        newPosition[i] = *oldPosition;
                    }
                } else { // We have to add traces
                    sum_traces(newPosition, oldPosition, traceStepSizes.data(), traceDimensions.data(), traceDimensions.size(), totalTraceDim);
                    for(size_t i = 1; i < _out.tensorObject->size; ++i) {
                        increase_indices(i, oldPosition, outAssIdx.numIndices-numOrderedIndices, stepSizes, outAssIdx.indexDimensions.data());
                        sum_traces(newPosition+i, oldPosition, traceStepSizes.data(), traceDimensions.data(), traceDimensions.size(), totalTraceDim);
                    }
                }
            } else { // We can copy/add larger blocks
                if(totalTraceDim == 1) { // We don't need to sum any traces
                    array_copy(newPosition, oldPosition, orderedIndexDim);
                    for(size_t i = 1; i < _out.tensorObject->size/orderedIndexDim; ++i) {
                        increase_indices(i, oldPosition, outAssIdx.numIndices-numOrderedIndices, stepSizes, outAssIdx.indexDimensions.data());
                        array_copy(newPosition + i*orderedIndexDim, oldPosition, orderedIndexDim);
                    }
                } else { // We have to add traces
                    sum_traces(newPosition, oldPosition, traceStepSizes.data(), traceDimensions.data(), traceDimensions.size(), totalTraceDim, orderedIndexDim);
                    for(size_t i = 1; i < _out.tensorObject->size/orderedIndexDim; ++i) {
                        increase_indices(i, oldPosition, outAssIdx.numIndices-numOrderedIndices, stepSizes, outAssIdx.indexDimensions.data());
                        sum_traces(newPosition + i*orderedIndexDim, oldPosition, traceStepSizes.data(), traceDimensions.data(), traceDimensions.size(), totalTraceDim, orderedIndexDim);
                    }
                }
            }
        }
        //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - Sparse => Both  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        else if(_base.tensorObjectReadOnly->is_sparse()) {
            VLA(bool[baseAssIdx.numIndices]  , fixedFlags); // Flag for each index indicating whether the index is fixed
            VLA(bool[baseAssIdx.numIndices]  , traceFlags); // Flag for each index indicating whether the index is part of a trace
            VLA(size_t[baseAssIdx.numIndices], attributes);  // Either the factor in _out, the value of an fixed index or the position of the other part of a trace
            bool peacefullIndices = true;
            
            // Process base indices
            const std::unique_ptr<const size_t[]> outIndexStepSizes = get_step_sizes(outAssIdx);
            for(size_t i = 0; i < baseAssIdx.numIndices; ++i) {
                // First try to find the index in out
                size_t outPos = 0;
                while(outPos < outAssIdx.numIndices && outAssIdx.indices[outPos] != baseAssIdx.indices[i]) { ++outPos; }
                if(outPos < outAssIdx.numIndices) {
                    fixedFlags[i] = false;
                    traceFlags[i] = false;
                    attributes[i] = outIndexStepSizes[outPos];
                    continue;
                }
                
                // Check whether the index is fixed
                if(baseAssIdx.indices[i].fixed()) {
                    fixedFlags[i] = true;
                    traceFlags[i] = false;
                    attributes[i] = (size_t) baseAssIdx.indices[i].valueId;
                    peacefullIndices = false;
                    continue;
                }
                
                // If the index is not in out and not fixed then it has to be part of a trace
                fixedFlags[i] = false;
                traceFlags[i] = true;
                peacefullIndices = false;
                for(attributes[i] = 0; baseAssIdx.indices[i] != baseAssIdx.indices[attributes[i]] || attributes[i] == i; ++attributes[i]) { }
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
                        outEntries.insert(std::pair<size_t, value_t>(get_position(entry, baseAssIdx.indexDimensions.data(), baseIndexStepSizes.get(), attributes, baseAssIdx.numIndices), factor*entry.second));
                    }
                } else {
                    for(const std::pair<size_t, value_t>& entry : baseEntries) {
                        size_t newPosition;
                        if(check_position(newPosition, entry, baseAssIdx.indexDimensions.data(), baseIndexStepSizes.get(), attributes, fixedFlags, traceFlags, baseAssIdx.numIndices)) {
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
                        dataPointer[get_position(entry, baseAssIdx.indexDimensions.data(), baseIndexStepSizes.get(), attributes, baseAssIdx.numIndices)] = factor*entry.second;
                    }
                } else {
                    for(const std::pair<size_t, value_t>& entry : baseEntries) {
                        size_t newPosition;
                        if(check_position(newPosition, entry, baseAssIdx.indexDimensions.data(), baseIndexStepSizes.get(), attributes, fixedFlags, traceFlags, baseAssIdx.numIndices)) {
                            dataPointer[newPosition] += factor*entry.second;
                        }
                    }
                }
            }
        }
    }
}
