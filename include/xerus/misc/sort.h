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



START_MISC_NAMESPACE

// NOT TESTED -- NOTE -- NOT TESTED

static _inline_ void merge(const size_t* const _leftPositions,
            const value_t* const _leftValues,
            const size_t _leftSize, 
            const size_t* const _rightPositions,
            const value_t* const _rightValues,
            const size_t _rightSize,
            size_t* _outPositions,
            value_t* const _outValues,
            size_t& _outSize
            ) {
    size_t leftIndex = 0, rightIndex = 0, outIndex = 0;
    
    while(leftIndex < _leftSize && rightIndex < _rightSize) {
        if(_leftPositions[leftIndex] == _rightPositions[rightIndex]) {
            _outValues[outIndex] = _leftValues[leftIndex] + _rightValues[rightIndex];
            _outPositions[outIndex] = _leftPositions[leftIndex] + _rightPositions[rightIndex];
            ++leftIndex; ++rightIndex;
        } else if(_leftPositions[leftIndex] < _rightPositions[rightIndex]) {
            _outValues[outIndex] = _leftValues[leftIndex];
            _outPositions[outIndex] = _leftPositions[leftIndex];
            ++leftIndex;
        } else {
            _outValues[outIndex] = _rightValues[leftIndex];
            _outPositions[outIndex] = _rightPositions[leftIndex];
            ++rightIndex;
        }
        ++outIndex;
    }
    _outSize = outIndex;
}

void merge_sort(const size_t* const _positions,
                const value_t* const _values,
                size_t* const _outPositions,
                value_t* const _outValues,
                size_t _size ) {
    // End of the recursion
    if(_size == 1) {
        _outPositions[0] = _positions[0];
        _outValues[0] = _values[0];
    }
    
    // Recursive call for the left half
    size_t leftSize = _size/2;
    std::unique_ptr<size_t[]> leftPositions(new size_t[leftSize]);
    std::unique_ptr<value_t[]> leftValues(new value_t[leftSize]);
    merge_sort(_positions, _values, leftPositions.get(), leftValues.get(), leftSize);
    
    // Recursive call for the right half
    size_t rightSize = _size/2+_size%2;
    std::unique_ptr<size_t[]> rightPositions(new size_t[leftSize]);
    std::unique_ptr<value_t[]> rightValues(new value_t[leftSize]);
    merge_sort(_positions+leftSize, _values+leftSize, rightPositions.get(), rightValues.get(), rightSize);
    
    // Merge both halfes
    merge(leftPositions.get(), leftValues.get(), leftSize, rightPositions.get(), rightValues.get(), rightSize, _outPositions, _outValues, _size);
}

END_MISC_NAMESPACE
    