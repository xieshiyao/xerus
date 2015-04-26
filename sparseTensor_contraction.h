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

#include "sparseTensor.h"

#include <suitesparse/cs.h>

// CSparse sparse matrix representation
// typedef struct cs_sparse {
//     int nzmax; // maximum number of entries
//     int m; // number of rows
//     int n; // number of columns
//     int *p; // column pointers (size n+1) or col indices (size nzmax)
//     int *i; // row indices, size nzmax
//     double *x; // numerical values, size nzmax
//     int nz; // # of entries in triplet matrix, -1 for compressed-col
// };


// #include "xerus.h"

namespace xerus {

    // Converts an Indexed SparseTensor and an given matrification to the CSparse sparse matrix format
    cs_di_sparse to_cs_format(const IndexedTensorReadOnly<Tensor>& _tensor, const std::vector<Index>& _lhsIndices, const std::vector<Index>& _rhsIndices);
    
    SparseTensor from_cs_format(const cs_di_sparse& _cs_format, const std::vector<size_t>& _dimensions);
    
    void print_cs(const cs_di_sparse& _cs_format);

    void matrix_matrix_product( std::set<value_t>& _C,
                                const size_t _leftDim,
                                const size_t _rightDim,
                                const double _alpha,
                                const std::set<value_t>& _A,
                                const bool _transposeA,
                                const size_t _middleDim,
                                const std::set<value_t>& _B,
                                const bool _transposeB);
    
}
