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
	namespace internal {
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
		
		///@brief Allocates a CS sparse matrix with given dimensions and number of entries.
		CsUniquePtr create_cs(const size_t _m, const size_t _n, const size_t _N); 

		///@brief: Converts the given @a _tensor to the CS format using the given matrification.
		CsUniquePtr to_cs_format(const Tensor& _tensor, const size_t _splitPos, const bool _transpose);
		
		///@brief: Retransforms a CS sparse matrix to sparse Tensor format
		void from_cs_format(Tensor& _output, const CsUniquePtr& _cs_format);
		
		///@brief: Calculates the Matrix Matrix product between to sparse matrices
		void matrix_matrix_product( std::map<size_t, double>& _C,
									const size_t _leftDim,
							const size_t _rightDim,
							const double _alpha,
							const std::map<size_t, double>& _A,
							const bool _transposeA,
							const size_t _midDim,
							const std::map<size_t, double>& _B,
							const bool _transposeB );
		
		/// Prints a matrix in cs format
		void print_cs(const CsUniquePtr& _cs_format);
	}
}
