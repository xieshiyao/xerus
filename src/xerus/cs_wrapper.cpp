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
* @brief Implementation of the suitesparse wrapper functions.
*/

#include <xerus/cs_wrapper.h>
#include <xerus/index.h>
#include <xerus/tensor.h>
#include <xerus/misc/performanceAnalysis.h>

#include <xerus/misc/check.h>

namespace xerus {
	namespace internal {
		
		CsUniquePtr create_cs(const size_t _m, const size_t _n, const size_t _N) {
			REQUIRE(_m < std::numeric_limits<int>::max() && _n < std::numeric_limits<int>::max() && _N < std::numeric_limits<int>::max(), "Sparse Tensor is to large for SuiteSparse (" << _m << " x " << _n << ", " << _N << ")");
			return CsUniquePtr(cs_spalloc((int) _m, (int) _n, (int) _N, 1, 0), &cs_spfree);
		}
		
		// Converts an Indexed SparseTensor and an given matrification to the CSparse sparse matrix format
		CsUniquePtr to_cs_format(const std::map<size_t, double>& _input, const size_t _m, const size_t _n, const bool _transpose) {
			CsUniquePtr cs_format = create_cs(_m, _n, _input.size());
			size_t entryPos = 0;
			cs_format->i[0] = 0;
			
			if(_transpose) {
				// Easy setting: We want A^T (m x n) in compressed coloum storage. We have A (n x m) and transform
				// it to compressed row storage. This results in A^T in compressed coloum storage as demanded.
				
				int currRow = -1;
				
				for(const std::pair<size_t, value_t>& entry : _input) {
					cs_format->x[entryPos] = entry.second;
					cs_format->i[entryPos] = (int) (entry.first%_m);
					while(currRow < (int) (entry.first/_m)) {
						cs_format->p[++currRow] = int(entryPos);
					}
					entryPos++;
				}
				
				REQUIRE(size_t(currRow) < _n && entryPos == _input.size(), "Internal Error " << currRow << ", " << _n << " | " << entryPos << ", " <<  _input.size());
				
				while(currRow < (int) _n) {
					cs_format->p[++currRow] = int(entryPos);
				}
				
			} else {
				
				// We need to get the ordering of the transpose first
				std::map<std::pair<int, int>, value_t> transposedData;
				for(const std::pair<size_t, value_t>& entry : _input) {
					transposedData.emplace(std::pair<int, int>(int(entry.first%_n), int(entry.first/_n)), entry.second);
				}
				
				int currCol = -1;
				
				for(const std::pair<std::pair<int, int>, value_t>& entry : transposedData) {
					cs_format->x[entryPos] = entry.second;
					cs_format->i[entryPos] = entry.first.second;
					while(currCol < entry.first.first) {
						cs_format->p[++currCol] = int(entryPos);
					}
					entryPos++;
				}
				
				REQUIRE(currCol < (int) _n && entryPos == _input.size(), "Internal Error " << currCol << ", " <<  _n << " | " << entryPos << ", " << _input.size());
				
				while(currCol < (int) _n) {
					cs_format->p[++currCol] = int(entryPos);
				}
			}
			
			return cs_format;
		}
		
		void from_cs_format(std::map<size_t, double>& _output, const CsUniquePtr& _cs_format, const double _alpha) {
			REQUIRE(_cs_format, "NullPtr cannot be converted to Tensor.");
			
			for(int i = 0; i < _cs_format->n; ++i) {
				for(int j = _cs_format->p[i]; j < _cs_format->p[i+1]; ++j) {
					IF_CHECK( auto ret = ) _output.emplace(_cs_format->i[j]*_cs_format->n+i, _alpha*_cs_format->x[j]);
					REQUIRE(ret.second, "Internal Error");
				}
			}
		}
		
		void matrix_matrix_product( std::map<size_t, double>& _C,
								const size_t _leftDim,
								const size_t _rightDim,
								const double _alpha,
								const std::map<size_t, double>& _A,
								const bool _transposeA,
								const size_t _midDim,
								const std::map<size_t, double>& _B,
								const bool _transposeB ) {
			const CsUniquePtr lhsCs = to_cs_format(_A, _leftDim, _midDim, _transposeA);
			const CsUniquePtr rhsCs = to_cs_format(_B, _midDim, _rightDim, _transposeB);
			const CsUniquePtr resultCs(cs_multiply(lhsCs.get(), rhsCs.get()), &cs_spfree);
			from_cs_format(_C, resultCs, _alpha);
		}
		
		void print_cs(const CsUniquePtr& _cs_format) {
			std::cout << "Sparse Matrix parameters: n = " << _cs_format->n << ", m = " << _cs_format->m << ", Max Entries = " << _cs_format->nzmax << std::endl;
			std::cout << "Values =  {";
			std::cout << _cs_format->x[0];
			for(int i = 1; i < _cs_format->nzmax; ++i) {
				std::cout << ", " << _cs_format->x[i];
			}
			std::cout << "}" << std::endl;
			
			std::cout << "Row idx = {";
			std::cout << _cs_format->i[0];
			for(int i = 1; i < _cs_format->nzmax; ++i) {
				std::cout << ", " << _cs_format->i[i];
			}
			std::cout << "}" << std::endl;
			
			std::cout << "Col Pos = {";
			std::cout << _cs_format->p[0];
			for(int i = 1; i < _cs_format->n + 1; ++i) {
				std::cout << ", " << _cs_format->p[i];
			}
			std::cout << "}" << std::endl;
		}
	}
}
