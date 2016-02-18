// Xerus - A General Purpose Tensor Library
// Copyright (C) 2014-2016 Benjamin Huber and Sebastian Wolf. 
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

#include <memory>
#include <map>
#include <mutex>

#include <suitesparse/cholmod.h>


namespace xerus {
	class Tensor;
	
	namespace internal {
		///@brief wrapper object for the cholmod_common struct to automatically call the constructor and destructor
		struct CholmodCommon final {
			struct RestrictedAccess final {
				cholmod_common* c;
				std::mutex &lock;
				RestrictedAccess(cholmod_common* _c, std::mutex &_lock);
				operator cholmod_common*() const;
				~RestrictedAccess();
			};
			
			std::unique_ptr<cholmod_common> c;
			std::mutex lock;
			CholmodCommon();
			~CholmodCommon();
			RestrictedAccess operator&();
			std::function<void(cholmod_sparse*)> get_deleter();
		};
		
		///@brief wrapper class for the cholmod sparse matrix objects
		class CholmodSparse final {
		public:
			using ptr_type = std::unique_ptr<cholmod_sparse, std::function<void(cholmod_sparse*)>>;
			ptr_type matrix;
			
			///@brief stores the given cholmod_sparse object in this wrapper. NOTE: the object must have been created by the same thread that uses this constructor!
			CholmodSparse(cholmod_sparse *_matrix);
			
			///@brief Creates a cholmod sparse matrix with given dimensions and number of entries.
			CholmodSparse(const size_t _m, const size_t _n, const size_t _N);
			
			///@brief Converts the given @a _tensor to the cholmod sparse format using the given matrification.
			CholmodSparse(const std::map<size_t, double>& _input, const size_t _m, const size_t _n, const bool _transpose);
			
			///@brief Transforms a cholmod sparse matrix to sparse Tensor format.
			std::map<size_t, double> to_map(double _alpha=1.0) const;
			
			///@brief Calculates the Matrix Matrix product with another sparse matrix.
			CholmodSparse operator*(const CholmodSparse &_rhs) const;
			
			///@brief Calculates the Matrix Matrix product between two sparse matrices.
			static void matrix_matrix_product(std::map<size_t, double>& _C,
										const size_t _leftDim,
										const size_t _rightDim,
										const double _alpha,
										const std::map<size_t, double>& _A,
										const bool _transposeA,
										const size_t _midDim,
										const std::map<size_t, double>& _B,
										const bool _transposeB);
		};
	}
}
