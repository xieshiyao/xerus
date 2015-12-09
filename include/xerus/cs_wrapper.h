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
 * @brief Header file for suitesparse wrapper functions.
 */

#pragma once

#ifdef __has_include
    #if __has_include(<cs.h>)
        #include <cs.h>
    #elif __has_include(<suitesparse/cs.h>)
        #include <suitesparse/cs.h>
	#else
		#pragma error no SuiteSparse found
    #endif
#else
    #include <cs.h>
#endif

#include "tensor.h"


namespace xerus {
    /// Unique_ptr wrapper that should always be used to encapsulate the CS sparse matrix format.
    typedef std::unique_ptr<cs_di, cs_di*(*)(cs_di*)> CsUniquePtr;
    
    /** @brief: Function that estimates whether a contraction result is sparse.
     * @param _lhsDim the first dimension of the left matrix and the result, i.e. N in an NxM * MxK contraction.
     * @param _midDim the second dimension of the left matrix and the first dimension of the right matrix, i.e. M in an NxM * MxK contraction.
     * @param _rhsDim the second dimension of the right matrix and the result, i.e. K in an NxM * MxK contraction.
     * @param _lhsEntries number of non-zero entries of the left matrix.
     * @param _rhsEntries number of non-zero entries of the right matrix.
     */
    bool sparse_result(const size_t _lhsDim, const size_t _midDim, const size_t _rhsDim, const size_t _lhsEntries, const size_t _rhsEntries);
    
    /// Allocates a CS sparse matrix with given dimensions and number of entries
    CsUniquePtr create_cs(const size_t _m, const size_t _n, const size_t _N); 

    // Converts an Indexed Tensor and a given matrification to the CSparse sparse matrix format
    CsUniquePtr to_cs_format(const IndexedTensorReadOnly<Tensor>& _tensor, const std::vector<Index>& _lhsIndices, const std::vector<Index>& _rhsIndices);
    
    /// Calculates the Matrix Matrix product between to CS sparse matrices
    CsUniquePtr matrix_matrix_product(const CsUniquePtr& _lhs, const CsUniquePtr& _rhs);
    
    /// Retransforms a CS sparse matrix to sparse Tensor format
    Tensor from_cs_format(const CsUniquePtr& _cs_format, const std::vector<size_t>& _dimensions);
    
    /// Prints a matrix in cs format
    void print_cs(const CsUniquePtr& _cs_format);
}
