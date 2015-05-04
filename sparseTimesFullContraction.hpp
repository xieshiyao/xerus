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
    void matrix_matrix_product( double* const _C,
                                const size_t _leftDim,
                                const size_t _rightDim,
                                const double _alpha,
                                const double* const _A,
                                const bool _transposeA,
                                const size_t _midDim,
                                const std::map<size_t, double>& _B,
                                const bool _transposeB) {
        // Prepare output array
        array_set_zero(_C, _leftDim*_rightDim);
        
        for(const std::pair<size_t, double>& entry : _B) {
            const size_t i = entry.first/_rightDim;
            const size_t j = entry.first%_rightDim;
            
            for(size_t k = 0; k < _leftDim; ++k) {
                _C[k*_leftDim+j] += entry.second*_A[k*_midDim+i];
            }
        }
        
        array_scale(_C, _alpha, _leftDim*_rightDim);
    }
    
    void matrix_matrix_product( double* const _C,
                                const size_t _leftDim,
                                const size_t _rightDim,
                                const double _alpha,
                                const std::map<size_t, double>& _A,
                                const bool _transposeA,
                                const size_t _midDim,
                                const double* const _B,
                                const bool _transposeB) {
        // Prepare output array
        array_set_zero(_C, _leftDim*_rightDim);
        
        for(const std::pair<size_t, double>& entry : _A) {
            const size_t i = entry.first/_rightDim;
            const size_t j = entry.first%_rightDim;
            
            array_add(_C+i*_rightDim, entry.second, _B+j*_rightDim, _rightDim);
        }
        
        array_scale(_C, _alpha, _leftDim*_rightDim);
    }
}
