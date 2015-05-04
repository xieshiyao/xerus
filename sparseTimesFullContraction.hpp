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
    
    //TODO this is most likely not efficent
    std::unique_ptr<double[]> transpose(const double* const _A, const size_t _leftDim, const size_t _rightDim) {
        double* const AT = new double[_leftDim*_rightDim];
        for(size_t i = 0; i < _leftDim; ++i) {
            for(size_t j = 0; j < _rightDim; ++j) {
                AT[j*_leftDim+i] = _A[i*_rightDim+j];
            }
        }
        return std::unique_ptr<double[]>(AT);
    }
    
    //TODO this is most likely not efficent
    void transpose_inplace(double* const _A, const size_t _leftDim, const size_t _rightDim) {
        for(size_t i = 0; i < _leftDim; ++i) {
            for(size_t j = i+1; j < _rightDim; ++j) {
                const double x = _A[i*_rightDim+j];
                _A[i*_rightDim+j] = _A[j*_leftDim+i];
                _A[j*_leftDim+i] = x;
            }
        }
    }
    
    // - - - - - - - - - - - - - - - - - - - - - - - - - Mix to Full - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    
    void matrix_matrix_product( double* const _C,
                                const size_t _leftDim,
                                const size_t _rightDim,
                                const double _alpha,
                                const std::map<size_t, double>& _A,
                                const bool _transposeA,
                                const size_t _midDim,
                                const double* const _B) {
        // Prepare output array
        array_set_zero(_C, _leftDim*_rightDim);
        
        // Transposition of A only changes how i and j are calculated
        if(!_transposeA) {
            for(const std::pair<size_t, double>& entry : _A) {
                const size_t i = entry.first/_rightDim;
                const size_t j = entry.first%_rightDim;
                array_add(_C+i*_rightDim, _alpha*entry.second, _B+j*_rightDim, _rightDim);
            }
        } else {
            for(const std::pair<size_t, double>& entry : _A) {
                const size_t i = entry.first%_rightDim;
                const size_t j = entry.first/_rightDim;
                array_add(_C+i*_rightDim, _alpha*entry.second, _B+j*_rightDim, _rightDim);
            }
        }
    }
    
    void matrix_matrix_product( double* const _C,
                                const size_t _leftDim,
                                const size_t _rightDim,
                                const double _alpha,
                                const double* const _A,
                                const bool _transposeA,
                                const size_t _midDim,
                                const std::map<size_t, double>& _B,
                                const bool _transposeB) {
        // It is significantly faster to calculate (B^T * A*T)^T
        if(!_transposeA) {
            const std::unique_ptr<double[]> AT = transpose(_A, _leftDim, _midDim);
            matrix_matrix_product(_C, _leftDim, _rightDim, _alpha, _B, !_transposeB, _midDim, AT.get());
        } else {
            matrix_matrix_product(_C, _leftDim, _rightDim, _alpha, _B, !_transposeB, _midDim, _A);
        }
        transpose_inplace(_C, _rightDim, _leftDim);
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
        if(_transposeB) {
            const std::unique_ptr<double[]> BT = transpose(_B, _midDim, _rightDim);
            matrix_matrix_product(_C, _leftDim, _rightDim, _alpha, _A, _transposeA, _midDim, BT.get());
        } else {
            matrix_matrix_product(_C, _leftDim, _rightDim, _alpha, _A, _transposeA, _midDim, _B);
        }
    }
    
    
    // - - - - - - - - - - - - - - - - - - - - - - - - - Mix to Sparse - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    
    
    void matrix_matrix_product( std::map<size_t, double>& _C,
                                const size_t _leftDim,
                                const size_t _rightDim,
                                const double _alpha,
                                const std::map<size_t, double>& _A,
                                const bool _transposeA,
                                const size_t _midDim,
                                const double* const _B) {
        std::map<size_t, double*> rows;
        
        // Transposition of A only changes how i and j are calculated
        if(!_transposeA) {
            for(const std::pair<size_t, double>& entry : _A) {
                const size_t i = entry.first/_rightDim;
                const size_t j = entry.first%_rightDim;
                auto row = rows.find(i);
                if(row != rows.end()) {
                    array_add(row->second, _alpha*entry.second, _B+j*_rightDim, _rightDim);
                } else {
                    double* const tmpPtr = new double[_rightDim];
                    rows.insert(std::pair<size_t, double*>(i, tmpPtr));
                    array_scaled_copy(tmpPtr, _alpha*entry.second, _B+j*_rightDim, _rightDim);
                }
            }
        } else {
            for(const std::pair<size_t, double>& entry : _A) {
                const size_t i = entry.first%_rightDim;
                const size_t j = entry.first/_rightDim;
                auto row = rows.find(i);
                if(row != rows.end()) {
                    array_add(row->second, _alpha*entry.second, _B+j*_rightDim, _rightDim);
                } else {
                    double* const tmpPtr = new double[_rightDim];
                    rows.insert(std::pair<size_t, double*>(i, tmpPtr));
                    array_scaled_copy(tmpPtr, _alpha*entry.second, _B+j*_rightDim, _rightDim);
                }
            }
        }
        
        for(const std::pair<size_t, double*>& row : rows) {
            for(size_t i = 0; i < _rightDim; ++i) {
                #pragma GCC diagnostic push
                #pragma GCC diagnostic ignored "-Wfloat-equal"
                if(row.second[i] != 0) {
                    _C.insert(std::pair<size_t, double>(row.first*_rightDim+i, row.second[i]));
                }
                #pragma GCC diagnostic pop
            }
            delete[] row.second;
        }
    }
    
}
