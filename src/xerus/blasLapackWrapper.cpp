// Xerus - A General Purpose Tensor Library
// Copyright (C) 2014-2017 Benjamin Huber and Sebastian Wolf. 
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
* @brief Implementation of the blas and lapack wrapper functions.
*/


#include <complex.h>
// fix for non standard-conform complex implementation
#undef I

// workaround for broken Lapack
#define lapack_complex_float    float _Complex
#define lapack_complex_double   double _Complex
extern "C"
{
	#include <cblas.h> 
}

#ifdef __has_include
	#if __has_include(<lapacke.h>)
		#include <lapacke.h>
	#elif __has_include(<lapacke/lapacke.h>)
		#include <lapacke/lapacke.h>
	#else
		#pragma error no lapacke found
	#endif
#else
	#include <lapacke.h>
#endif


#include <memory>
#include <xerus/misc/standard.h>
#include <xerus/misc/performanceAnalysis.h>
#include <xerus/misc/check.h>

#include <xerus/misc/stringUtilities.h>
#include <xerus/basic.h>

#include <xerus/blasLapackWrapper.h>
#include <xerus/misc/basicArraySupport.h>
#include <xerus/misc/math.h>
#include <xerus/misc/internal.h>



namespace xerus {
	namespace blasWrapper {
		// the following routines use a work array as temporary storage
		// to avoid the overhead of repeated reallocation for many small calls, every thread pre-allocates a small scratch-space
// 		const size_t DEFAULT_WORKSPACE_SIZE = 1024*1024;
// 		thread_local value_t defaultWorkspace[DEFAULT_WORKSPACE_SIZE]; // NOTE recheck compatibility with eigen (dolfin) when reinserting this!
		
		
		//----------------------------------------------- LEVEL I BLAS ----------------------------------------------------------
		
		double one_norm(const double* const _x, const size_t _n) {
			REQUIRE(_n <= static_cast<size_t>(std::numeric_limits<int>::max()), "Dimension to large for BLAS/Lapack");
			
			XERUS_PA_START;
			
			const double result = cblas_dasum(static_cast<int>(_n), _x, 1);
			
			XERUS_PA_END("Dense BLAS", "One Norm", misc::to_string(_n));
			
			return result;
		}
		
		double two_norm(const double* const _x, const size_t _n) {
			REQUIRE(_n <= static_cast<size_t>(std::numeric_limits<int>::max()), "Dimension to large for BLAS/Lapack");
			
			XERUS_PA_START;
			
			const double result = cblas_dnrm2(static_cast<int>(_n), _x, 1);
			
			XERUS_PA_END("Dense BLAS", "Two Norm", misc::to_string(_n));
			
			return result;
		}
		
		double dot_product(const double* const _x, const size_t _n, const double* const _y) {
			REQUIRE(_n <= static_cast<size_t>(std::numeric_limits<int>::max()), "Dimension to large for BLAS/Lapack");
			
			XERUS_PA_START;
			
			const double result = cblas_ddot(static_cast<int>(_n), _x, 1, _y, 1);
			
			XERUS_PA_END("Dense BLAS", "Dot Product", misc::to_string(_n)+"*"+misc::to_string(_n));
			
			return result;
		}
		
		
		//----------------------------------------------- LEVEL II BLAS ---------------------------------------------------------
		
		void matrix_vector_product(double* const _x, const size_t _m, const double _alpha, const double* const _A, const size_t _n, const bool _transposed, const double* const _y) {
			// Delegate call if appropriate
			if(_m == 1) { *_x = _alpha*dot_product(_A, _n, _y); return;}
			
			REQUIRE(_m <= static_cast<size_t>(std::numeric_limits<int>::max()), "Dimension to large for BLAS/Lapack");
			REQUIRE(_n <= static_cast<size_t>(std::numeric_limits<int>::max()), "Dimension to large for BLAS/Lapack");
			
			XERUS_PA_START;
			if(!_transposed) {
				cblas_dgemv(CblasRowMajor, CblasNoTrans, static_cast<int>(_m), static_cast<int>(_n), _alpha, _A, static_cast<int>(_n), _y, 1, 0.0, _x, 1);
			} else {
				cblas_dgemv(CblasRowMajor, CblasTrans, static_cast<int>(_n), static_cast<int>(_m), _alpha, _A, static_cast<int>(_m) , _y, 1, 0.0, _x, 1);
			}
			
			XERUS_PA_END("Dense BLAS", "Matrix Vector Product", misc::to_string(_m)+"x"+misc::to_string(_n)+" * "+misc::to_string(_n));
		}
		
		void dyadic_vector_product(double* _A, const size_t _m, const size_t _n, const double _alpha, const double*const  _x, const double* const _y) {
			REQUIRE(_m <= static_cast<size_t>(std::numeric_limits<int>::max()), "Dimension to large for BLAS/Lapack");
			REQUIRE(_n <= static_cast<size_t>(std::numeric_limits<int>::max()), "Dimension to large for BLAS/Lapack");
			
			XERUS_PA_START;
			
			//Blas wants to add the product to A, but we don't.
			misc::set_zero(_A, _m*_n);
			
			cblas_dger(CblasRowMajor, static_cast<int>(_m), static_cast<int>(_n), _alpha, _x, 1, _y, 1, _A, static_cast<int>(_n));
			
			XERUS_PA_END("Dense BLAS", "Dyadic Vector Product", misc::to_string(_m)+" o "+misc::to_string(_n));
		}
		
		
		//----------------------------------------------- LEVEL III BLAS --------------------------------------------------------
		/// Performs the Matrix-Matrix product c = a * b
		void matrix_matrix_product( double* const _C,
									const size_t _leftDim,
									const size_t _rightDim,
									const double _alpha,
									const double* const _A,
									const size_t _lda,
									const bool _transposeA,
									const size_t _middleDim,
									const double* const _B,
									const size_t _ldb,
									const bool _transposeB) {
			//Delegate call if appropriate
			if(_leftDim == 1) {
				matrix_vector_product(_C, _rightDim, _alpha, _B, _middleDim, !_transposeB, _A);
			} else if(_rightDim == 1) {
				matrix_vector_product(_C, _leftDim, _alpha, _A, _middleDim, _transposeA, _B);
			} else if(_middleDim == 1) { 
				dyadic_vector_product(_C, _leftDim, _rightDim, _alpha, _A, _B);
			} else {
			
				REQUIRE(_leftDim <= static_cast<size_t>(std::numeric_limits<int>::max()), "Dimension to large for BLAS/Lapack");
				REQUIRE(_middleDim <= static_cast<size_t>(std::numeric_limits<int>::max()), "Dimension to large for BLAS/Lapack");
				REQUIRE(_rightDim <= static_cast<size_t>(std::numeric_limits<int>::max()), "Dimension to large for BLAS/Lapack");
				REQUIRE(_lda <= static_cast<size_t>(std::numeric_limits<int>::max()), "Dimension to large for BLAS/Lapack");
				REQUIRE(_ldb <= static_cast<size_t>(std::numeric_limits<int>::max()), "Dimension to large for BLAS/Lapack");
				
				XERUS_PA_START;
				
				cblas_dgemm( CblasRowMajor,                             // Array storage format
						_transposeA ? CblasTrans : CblasNoTrans,        // LHS transposed?
						_transposeB ? CblasTrans : CblasNoTrans,        // RHS transposed?
						static_cast<int>(_leftDim),                     // Left dimension
						static_cast<int>(_rightDim),                    // Right dimension
						static_cast<int>(_middleDim),                   // Middle dimension
						_alpha,                                         // Factor to the Product
						_A,                                             // Pointer to LHS
						static_cast<int>(_lda),                         // LDA
						_B,                                             // Pointer to RHS
						static_cast<int>(_ldb),                         // LDB
						0.0,                                            // Factor of C (Zero if only the product is required)
						_C,                                             // Pointer to result
						static_cast<int>(_rightDim)                     // LDC
				);
				
				XERUS_PA_END("Dense BLAS", "Matrix-Matrix-Multiplication", misc::to_string(_leftDim)+"x"+misc::to_string(_middleDim)+" * "+misc::to_string(_middleDim)+"x"+misc::to_string(_rightDim));
			}
		}
		
		
		
		//----------------------------------------------- LAPACK ----------------------------------------------------------------
		
		void svd( double* const _U, double* const _S, double* const _Vt, const double* const _A, const size_t _m, const size_t _n) {
			//Create copy of A
			const std::unique_ptr<double[]> tmpA(new double[_m*_n]);
			misc::copy(tmpA.get(), _A, _m*_n);
			
			svd_destructive(_U, _S, _Vt, tmpA.get(), _m, _n);
		}
		
		
		void svd_destructive( double* const _U, double* const _S, double* const _Vt, double* const _A, const size_t _m, const size_t _n) {
			REQUIRE(_m <= static_cast<size_t>(std::numeric_limits<int>::max()), "Dimension to large for BLAS/Lapack");
			REQUIRE(_n <= static_cast<size_t>(std::numeric_limits<int>::max()), "Dimension to large for BLAS/Lapack");
			
			XERUS_PA_START;
			
            int lapackAnswer = LAPACKE_dgesdd(LAPACK_ROW_MAJOR, 'S', static_cast<int>(_m), static_cast<int>(_n), _A, static_cast<int>(_n), _S, _U, static_cast<int>(std::min(_m, _n)), _Vt, static_cast<int>(_n));
			CHECK(lapackAnswer == 0, warning, "Lapack failed to compute SVD. Answer is: " << lapackAnswer);
			CHECK(lapackAnswer == 0, warning, "Call was: LAPACKE_dgesdd(LAPACK_ROW_MAJOR, 'S', " << static_cast<int>(_m) << ", " << static_cast<int>(_n) << ", " << _A << ", " << static_cast<int>(_n) <<", " 
			<< _S <<", " << _U << ", " << static_cast<int>(std::min(_m, _n)) << ", " << _Vt << ", " << static_cast<int>(_n) << ");");
			if(lapackAnswer != 0) {
				LOG(warning, "SVD failed ");
			}
			
			XERUS_PA_END("Dense LAPACK", "Singular Value Decomposition", misc::to_string(_m)+"x"+misc::to_string(_n));
		}
		
		
		std::tuple<std::unique_ptr<double[]>, std::unique_ptr<double[]>, size_t> qc(const double* const _A, const size_t _m, const size_t _n) {
			const std::unique_ptr<double[]> tmpA(new double[_m*_n]);
			misc::copy(tmpA.get(), _A, _m*_n);
			
			return qc_destructive(tmpA.get(), _m, _n);
		}
		
		
		std::tuple<std::unique_ptr<double[]>, std::unique_ptr<double[]>, size_t> qc_destructive(double* const _A, const size_t _m, const size_t _n) {
			REQUIRE(_m <= static_cast<size_t>(std::numeric_limits<int>::max()), "Dimension to large for BLAS/Lapack");
			REQUIRE(_n <= static_cast<size_t>(std::numeric_limits<int>::max()), "Dimension to large for BLAS/Lapack");
			
			REQUIRE(_n > 0, "Dimension n must be larger than zero");
			REQUIRE(_m > 0, "Dimension m must be larger than zero");
			
			XERUS_PA_START;
			
			// Maximal rank is used by Lapacke
			const size_t maxRank = std::min(_m, _n);
			
			// Tmp Array for Lapacke
			const std::unique_ptr<double[]> tau(new double[maxRank]);
			
			const std::unique_ptr<int[]> permutation(new int[_n]());
			misc::set_zero(permutation.get(), _n); // Lapack requires the entries to be zero.
			
			// Calculate QR factorisations with column pivoting
			IF_CHECK(int lapackAnswer = ) LAPACKE_dgeqp3(LAPACK_ROW_MAJOR, static_cast<int>(_m), static_cast<int>(_n), _A, static_cast<int>(_n), permutation.get(), tau.get());
			REQUIRE(lapackAnswer == 0, "Unable to perform QC factorisaton (dgeqp3). Lapacke says: " << lapackAnswer );
			
			
			// Determine the actual rank
			size_t rank;
			for (rank = 1; rank <= maxRank; ++rank) {
				if (rank == maxRank || std::abs(_A[rank+rank*_n]) < 16*std::numeric_limits<double>::epsilon()*_A[0]) {
					break;
				}
			}
			
			
			// Create the matrix C
			std::unique_ptr<double[]> C(new double[rank*_n]);
			misc::set_zero(C.get(), rank*_n); 
			
			// Copy the upper triangular Matrix C (rank x _n) into position
			for (size_t col = 0; col < _n; ++col) {
				const size_t targetCol = static_cast<size_t>(permutation[col]-1); // For Lapack numbers start at 1 (instead of 0).
				for(size_t row = 0; row < rank && row < col+1; ++row) {
					C[row*_n + targetCol] = _A[row*_n + col];
				}
			}
			
			
			// Create orthogonal matrix Q
			IF_CHECK(lapackAnswer = ) LAPACKE_dorgqr(LAPACK_ROW_MAJOR, static_cast<int>(_m), static_cast<int>(maxRank), static_cast<int>(maxRank), _A, static_cast<int>(_n), tau.get());
			CHECK(lapackAnswer == 0, error, "Unable to reconstruct Q from the QC factorisation. Lapacke says: " << lapackAnswer);
			
			// Copy the newly created Q into position
			std::unique_ptr<double[]> Q(new double[_m*rank]);
			if(rank == _n) {
				misc::copy(Q.get(), _A, _m*rank);
			} else {
				for(size_t row = 0; row < _m; ++row) {
					misc::copy(Q.get()+row*rank, _A+row*_n, rank);
				}
			}
			
			XERUS_PA_END("Dense LAPACK", "QRP Factorisation", misc::to_string(_m)+"x"+misc::to_string(rank)+" * "+misc::to_string(rank)+"x"+misc::to_string(_n));
			
			return std::make_tuple(std::move(Q), std::move(C), rank);
		}
		
		
		std::tuple<std::unique_ptr<double[]>, std::unique_ptr<double[]>, size_t> cq(const double* const _A, const size_t _m, const size_t _n) {
			const std::unique_ptr<double[]> tmpA(new double[_m*_n]);
			misc::copy(tmpA.get(), _A, _m*_n);
			
			return cq_destructive(tmpA.get(), _m, _n);
		}
		
		
		// We use that in col-major we get At = Qt * Ct => A = C * Q, i.e. doing the calculation in col-major and switching Q and C give the desired result.
		std::tuple<std::unique_ptr<double[]>, std::unique_ptr<double[]>, size_t> cq_destructive(double* const _A, const size_t _m, const size_t _n) {
			REQUIRE(_n <= static_cast<size_t>(std::numeric_limits<int>::max()), "Dimension to large for BLAS/Lapack");
			REQUIRE(_m <= static_cast<size_t>(std::numeric_limits<int>::max()), "Dimension to large for BLAS/Lapack");
			
			REQUIRE(_m > 0, "Dimension m must be larger than zero");
			REQUIRE(_n > 0, "Dimension n must be larger than zero");
			
			XERUS_PA_START;
			
			// Maximal rank is used by Lapacke
			const size_t maxRank = std::min(_n, _m);
			
			// Tmp Array for Lapacke
			const std::unique_ptr<double[]> tau(new double[maxRank]);
			
			const std::unique_ptr<int[]> permutation(new int[_m]());
			misc::set_zero(permutation.get(), _m); // Lapack requires the entries to be zero.
			
			// Calculate QR factorisations with column pivoting
			IF_CHECK(int lapackAnswer = ) LAPACKE_dgeqp3(LAPACK_COL_MAJOR, static_cast<int>(_n), static_cast<int>(_m), _A, static_cast<int>(_n), permutation.get(), tau.get());
			REQUIRE(lapackAnswer == 0, "Unable to perform QC factorisaton (dgeqp3). Lapacke says: " << lapackAnswer );
			
			
			// Determine the actual rank
			size_t rank;
			for (rank = 1; rank <= maxRank; ++rank) {
				if (rank == maxRank || std::abs(_A[rank+rank*_n]) < 16*std::numeric_limits<double>::epsilon()*_A[0]) {
					break;
				}
			}
			
			
			// Create the matrix C
			std::unique_ptr<double[]> C(new double[rank*_m]);
			misc::set_zero(C.get(), rank*_m); 
			
			// Copy the upper triangular Matrix C (rank x _m) into position
			for (size_t col = 0; col < _m; ++col) {
				const size_t targetCol = static_cast<size_t>(permutation[col]-1); // For Lapack numbers start at 1 (instead of 0).
				misc::copy(C.get()+targetCol*rank, _A+col*_n, std::min(rank, col+1));
			}
			
			
			// Create orthogonal matrix Q
			IF_CHECK(lapackAnswer = ) LAPACKE_dorgqr(LAPACK_COL_MAJOR, static_cast<int>(_n), static_cast<int>(maxRank), static_cast<int>(maxRank), _A, static_cast<int>(_n), tau.get());
			CHECK(lapackAnswer == 0, error, "Unable to reconstruct Q from the QC factorisation. Lapacke says: " << lapackAnswer);
			
			// Copy the newly created Q into position
			std::unique_ptr<double[]> Q(new double[_n*rank]);
			misc::copy(Q.get(), _A, _n*rank);
			
			XERUS_PA_END("Dense LAPACK", "QRP Factorisation", misc::to_string(_n)+"x"+misc::to_string(rank)+" * "+misc::to_string(rank)+"x"+misc::to_string(_m));
			
			return std::make_tuple(std::move(C), std::move(Q), rank);
		}
		
		
		void qr(double* const _Q, double* const _R, const double* const _A, const size_t _m, const size_t _n) {
			// Create tmp copy of A since Lapack wants to destroy it
			const std::unique_ptr<double[]> tmpA(new double[_m*_n]);
			misc::copy(tmpA.get(), _A, _m*_n);
			
			qr_destructive(_Q, _R, tmpA.get(), _m, _n);
		}
		
		
		void inplace_qr(double* const _AtoQ, double* const _R, const size_t _m, const size_t _n) {
			qr_destructive(_AtoQ, _R, _AtoQ, _m, _n);
		}
		
		
		void qr_destructive( double* const _Q, double* const _R, double* const _A, const size_t _m, const size_t _n) {
			REQUIRE(_m <= static_cast<size_t>(std::numeric_limits<int>::max()), "Dimension to large for BLAS/Lapack");
			REQUIRE(_n <= static_cast<size_t>(std::numeric_limits<int>::max()), "Dimension to large for BLAS/Lapack");
			
			REQUIRE(_n > 0, "Dimension n must be larger than zero");
			REQUIRE(_m > 0, "Dimension m must be larger than zero");
			
			REQUIRE(_Q && _R && _A, "QR decomposition must not be called with null pointers: Q:" << _Q << " R: " << _R << " A: " << _A);
			REQUIRE(_A != _R, "_A and _R must be different, otherwise qr call will fail.");
			
			XERUS_PA_START;
			
			// Maximal rank is used by Lapacke
			const size_t rank = std::min(_m, _n); 
			
			// Tmp Array for Lapacke
			const std::unique_ptr<double[]> tau(new double[rank]);
			
			// Calculate QR factorisations
//             LOG(Lapacke, "Call to dorgqr with parameters: " << LAPACK_ROW_MAJOR << ", " << static_cast<int>(_m)  << ", " << static_cast<int>(_n)  << ", " << _A << ", " << static_cast<int>(_n)  << ", " << tau.get());
			IF_CHECK( int lapackAnswer = ) LAPACKE_dgeqrf(LAPACK_ROW_MAJOR, static_cast<int>(_m), static_cast<int>(_n), _A, static_cast<int>(_n), tau.get());
			CHECK(lapackAnswer == 0, error, "Unable to perform QR factorisaton. Lapacke says: " << lapackAnswer );
			
			// Copy the upper triangular Matrix R (rank x _n) into position
			for(size_t row =0; row < rank; ++row) {
				misc::set_zero(_R+row*_n, row); // Set starting zeros
				misc::copy(_R+row*_n+row, _A+row*_n+row, _n-row); // Copy upper triangular part from lapack result.
			}
			
			// Create orthogonal matrix Q (in tmpA)
	//         LOG(Lapacke, "Call to dorgqr with parameters: " << LAPACK_ROW_MAJOR << ", " << static_cast<int>(_m)  << ", " << static_cast<int>(rank)  << ", " << static_cast<int>(rank) << ", " << _A << ", " << static_cast<int>(_n)  << ", " << tau.get());
			IF_CHECK( lapackAnswer = ) LAPACKE_dorgqr(LAPACK_ROW_MAJOR, static_cast<int>(_m), static_cast<int>(rank), static_cast<int>(rank), _A, static_cast<int>(_n), tau.get());
			CHECK(lapackAnswer == 0, error, "Unable to reconstruct Q from the QR factorisation. Lapacke says: " << lapackAnswer);
			
			// Copy Q (_m x rank) into position
			if(_A != _Q) {
				if(_n == rank) {
					misc::copy(_Q, _A, _m*_n);
				} else {
					for(size_t row = 0; row < _m; ++row) {
						misc::copy(_Q+row*rank, _A+row*_n, rank);
					}
				}
			} else if(_n != rank) { // Note extra treatmeant to avoid memcpy overlap
				for(size_t row = 1; row < _m; ++row) {
					misc::copy_inplace(_Q+row*rank, _A+row*_n, rank);
				}
			}
			
			XERUS_PA_END("Dense LAPACK", "QR Factorisation", misc::to_string(_m)+"x"+misc::to_string(_n));
		}
		
		
		void rq( double* const _R, double* const _Q, const double* const _A, const size_t _m, const size_t _n) {
			// Create tmp copy of A since Lapack wants to destroy it
			const std::unique_ptr<double[]> tmpA(new double[_m*_n]);
			misc::copy(tmpA.get(), _A, _m*_n);
			
			rq_destructive(_R, _Q, tmpA.get(), _m, _n);
		}
		
		
		void inplace_rq( double* const _R, double* const _AtoQ, const size_t _m, const size_t _n) {
			rq_destructive(_R, _AtoQ, _AtoQ, _m, _n);
		}
		
		
		void rq_destructive( double* const _R, double* const _Q, double* const _A, const size_t _m, const size_t _n) {
			REQUIRE(_m <= static_cast<size_t>(std::numeric_limits<int>::max()), "Dimension to large for BLAS/Lapack");
			REQUIRE(_n <= static_cast<size_t>(std::numeric_limits<int>::max()), "Dimension to large for BLAS/Lapack");
			
			REQUIRE(_n > 0, "Dimension n must be larger than zero");
			REQUIRE(_m > 0, "Dimension m must be larger than zero");
			
			REQUIRE(_Q && _R && _A, "QR decomposition must not be called with null pointers: R " << _R << " Q: " << _Q << " A: " << _A);
			REQUIRE(_A != _R, "_A and _R must be different, otherwise qr call will fail.");
			
			XERUS_PA_START;
			
			// Maximal rank is used by Lapacke
			const size_t rank = std::min(_m, _n); 
			
			// Tmp Array for Lapacke
			const std::unique_ptr<double[]> tau(new double[rank]);
			
			IF_CHECK( int lapackAnswer = ) LAPACKE_dgerqf(LAPACK_ROW_MAJOR, static_cast<int>(_m), static_cast<int>(_n), _A, static_cast<int>(_n), tau.get());
			CHECK(lapackAnswer == 0, error, "Unable to perform QR factorisaton. Lapacke says: " << lapackAnswer << ". Call was: LAPACKE_dgerqf(LAPACK_ROW_MAJOR, "<<static_cast<int>(_m)<<", "<<static_cast<int>(_n)<<", "<<_A<<", "<<static_cast<int>(_n)<<", "<<tau.get()<<");" );
			
			
			// Copy the upper triangular Matrix R (_m x rank) into position.
			size_t row = 0;
			for( ; row < _m - rank; ++row) {
				misc::copy(_R+row*rank, _A+row*_n+_n-rank, rank);
			}
			for(size_t skip = 0; row < _m; ++row, ++skip) {
				misc::set_zero(_R+row*rank, skip); // Set starting zeros
				misc::copy(_R+row*rank+skip, _A+row*_n+_n-rank+skip, rank-skip); // Copy upper triangular part from lapack result.
			}
			
			// Create orthogonal matrix Q (in _A). Lapacke expects to get the last rank rows of A...
			IF_CHECK( lapackAnswer = ) LAPACKE_dorgrq(LAPACK_ROW_MAJOR, static_cast<int>(rank), static_cast<int>(_n), static_cast<int>(rank), _A+(_m-rank)*_n, static_cast<int>(_n), tau.get()); 
			CHECK(lapackAnswer == 0, error, "Unable to reconstruct Q from the RQ factorisation. Lapacke says: " << lapackAnswer << ". Call was: LAPACKE_dorgrq(LAPACK_ROW_MAJOR, "<<static_cast<int>(rank)<<", "<<static_cast<int>(_n)<<", "<<static_cast<int>(rank)<<", "<<_A+(_m-rank)*_n<<", "<<static_cast<int>(_n)<<", "<<tau.get()<<");");

			
			//Copy Q (rank x _n) into position
			if(_A != _Q) {
				misc::copy(_Q, _A+(_m-rank)*_n, rank*_n);
			}
			
			XERUS_PA_END("Dense LAPACK", "RQ Factorisation", misc::to_string(_m)+"x"+misc::to_string(_n));
		}
		
		
		static bool is_symmetric(const double* const _A, const size_t _n) {
			double max = 0;
			for (size_t i=0; i<_n*_n; ++i) {
				max = std::max(max, _A[i]);
			}
			
			for (size_t i=0; i<_n; ++i) {
				for (size_t j=i+1; j<_n; ++j) {
					if (std::abs(_A[i*_n + j] - _A[i + j*_n]) >= 4 * max * std::numeric_limits<double>::epsilon()) {
// 						LOG(aslkdjj, std::abs(_A[i*_n + j] - _A[i + j*_n]) << " / " << _A[i*_n + j] << " " << max);
						return false;
					}
				}
			}
			return true;
		}
		
		/// @brief checks whether the diagonal of @a _A is all positive or all negative. returns false otherwise
		static bool pos_neg_definite_diagonal(const double* const _A, const size_t _n) {
			bool positive = (_A[0] > 0);
			const size_t steps=_n+1;
			if (positive) {
				for (size_t i=1; i<_n; ++i) {
					if (_A[i*steps] < std::numeric_limits<double>::epsilon()) {
						return false;
					}
				}
				return true;
			} else {
				for (size_t i=1; i<_n; ++i) {
					if (_A[i*steps] > -std::numeric_limits<double>::epsilon()) {
						return false;
					}
				}
				return true;
			}
			
		}
		
		/// Solves Ax = b for x
		/// order of checks and solvers inspired by matlabs mldivide https://de.mathworks.com/help/matlab/ref/mldivide.html
		void solve(double* const _x, const double* const _A, const size_t _m, const size_t _n, const double* const _b, const size_t _nrhs) {
			REQUIRE(_m <= static_cast<size_t>(std::numeric_limits<int>::max()), "Dimension to large for BLAS/Lapack");
			REQUIRE(_n <= static_cast<size_t>(std::numeric_limits<int>::max()), "Dimension to large for BLAS/Lapack");
			REQUIRE(_nrhs <= static_cast<size_t>(std::numeric_limits<int>::max()), "Dimension to large for BLAS/Lapack");
			
			const std::unique_ptr<double[]> tmpA(new double[_m*_n]);
			misc::copy(tmpA.get(), _A, _m*_n);
			
			LOG(debug, "solving with...");
			
			// not rectangular -> fallback to least-squares (QR or SVD decomposition to solve)
			if (_m != _n) {
				LOG(debug, "SVD");
				const std::unique_ptr<double[]> tmpB(new double[_m*_nrhs]);
				misc::copy(tmpB.get(), _b, _m*_nrhs);
				solve_least_squares_destructive(_x, tmpA.get(), _m, _n, tmpB.get(), _nrhs);
				return;
			}
			
			// not symmetric -> LU solver
			if (!is_symmetric(_A, _n)) {
				LOG(debug, "LU");
				XERUS_PA_START;
				
				std::unique_ptr<int[]> pivot(new int[_n]);
				
				misc::copy(_x, _b, _n);
				
				IF_CHECK( int lapackAnswer = ) LAPACKE_dgesv(
					LAPACK_ROW_MAJOR,
					static_cast<int>(_n),		// Dimensions of A (nxn)
					static_cast<int>(_nrhs), 	// Number of rhs b
					tmpA.get(),					// input: A, output: L and U
					static_cast<int>(_n),		// LDA
					pivot.get(),				// output: permutation P
					_x,							// input: rhs b, output: solution x
					static_cast<int>(_nrhs)		// ldb
				);
				CHECK(lapackAnswer == 0, error, "Unable to solve Ax = b (PLU solver). Lapacke says: " << lapackAnswer);
				
				XERUS_PA_END("Dense LAPACK", "Solve (PLU)", misc::to_string(_n)+"x"+misc::to_string(_n)+"x"+misc::to_string(_nrhs));
				return;
			}
			
			// positive or negative diagonal -> try cholesky
			if (pos_neg_definite_diagonal(_A, _n)) {
				LOG(debug, "trying cholesky");
				int lapackAnswer = 0;
				{
					XERUS_PA_START;
					
					lapackAnswer = LAPACKE_dpotrf2(
						LAPACK_ROW_MAJOR,
						'U', 					// Upper triangle of A is read
						static_cast<int>(_n),	// dimensions of A
						tmpA.get(),				// input: A, output: cholesky factorisation
						static_cast<int>(_n)	// LDA
					);
					
					XERUS_PA_END("Dense LAPACK", "Cholesky decomposition", misc::to_string(_n)+"x"+misc::to_string(_n));
				}
				
				if (lapackAnswer == 0) {
					LOG(debug, "cholesky");
					XERUS_PA_START;
					
					misc::copy(_x, _b, _n);
					
					lapackAnswer = LAPACKE_dpotrs(
						LAPACK_ROW_MAJOR,
						'U',					// upper triangle of cholesky decomp is stored in tmpA
						static_cast<int>(_n),	// dimensions of A
						static_cast<int>(_nrhs),// number of rhs
						tmpA.get(),				// input: cholesky decomp
						static_cast<int>(_n), 	// lda
						_x,						// input: rhs b, output: solution x
						static_cast<int>(_nrhs)	// ldb
					);
					CHECK(lapackAnswer == 0, error, "Unable to solve Ax = b (cholesky solver). Lapacke says: " << lapackAnswer);
					
					XERUS_PA_END("Dense LAPACK", "Solve (Cholesky)", misc::to_string(_n)+"x"+misc::to_string(_n)+"x"+misc::to_string(_nrhs));
					
					return;
				} else {
					// restore tmpA
					misc::copy(tmpA.get(), _A, _m*_n);
				}
			}
			
			LOG(debug, "LDL");
			// non-definite diagonal or choleksy failed -> fallback to LDL^T decomposition
			XERUS_PA_START;
			
			misc::copy(_x, _b, _n);
			std::unique_ptr<int[]> pivot(new int[_n]);
			
			LAPACKE_dsysv(
				LAPACK_ROW_MAJOR,
				'U',					// upper triangular part of _A is accessed
				static_cast<int>(_n),	// dimension of A
				static_cast<int>(_nrhs),// number of rhs
				tmpA.get(), 			// input: A, output: part of the LDL decomposition
				static_cast<int>(_n),	// lda
				pivot.get(), 			// output: details of blockstructure of D
				_x,						// input: rhs b, output: solution x
				static_cast<int>(_nrhs)	// ldb
			);
			
			XERUS_PA_END("Dense LAPACK", "Solve (LDL)", misc::to_string(_n)+"x"+misc::to_string(_n)+"x"+misc::to_string(_nrhs));
		}
		
	
		void solve_least_squares( double* const _x, const double* const _A, const size_t _m, const size_t _n, const double* const _b, const size_t _p){
			const std::unique_ptr<double[]> tmpA(new double[_m*_n]);
			misc::copy(tmpA.get(), _A, _m*_n);
			
			const std::unique_ptr<double[]> tmpB(new double[_m*_p]);
			misc::copy(tmpB.get(), _b, _m*_p);
			
			solve_least_squares_destructive(_x, tmpA.get(), _m, _n, tmpB.get(), _p);
		}
		
		
		void solve_least_squares_destructive( double* const _x, double* const _A, const size_t _m, const size_t _n, double* const _b, const size_t _p){
			REQUIRE(_m <= static_cast<size_t>(std::numeric_limits<int>::max()), "Dimension to large for BLAS/Lapack");
			REQUIRE(_n <= static_cast<size_t>(std::numeric_limits<int>::max()), "Dimension to large for BLAS/Lapack");
			REQUIRE(_p <= static_cast<size_t>(std::numeric_limits<int>::max()), "Dimension to large for BLAS/Lapack");
			
			XERUS_PA_START;
			
			std::unique_ptr<int[]> pivot(new int[_n]);
			misc::set_zero(pivot.get(), _n);
			
			std::unique_ptr<double[]> signulars(new double[std::min(_n, _m)]);
			
			int rank;
			
			double* bOrX;
			if(_m >= _n) {
				bOrX = _b;
			} else {
				bOrX = _x;
				misc::copy(bOrX, _b, _m*_p);
				misc::set_zero(bOrX+_m*_p, (_n-_m)*_p); // Lapacke is unhappy if the array contains NANs...
			}
			
// 			IF_CHECK( int lapackAnswer = ) LAPACKE_dgelsy(
// 				LAPACK_ROW_MAJOR, 
// 				static_cast<int>(_m),   // Left dimension of A
// 				static_cast<int>(_n),   // Right dimension of A
// 				static_cast<int>(_p),	// Number of rhss
// 				_A,         			// Matrix A
// 				static_cast<int>(_n),   // LDA
// 				bOrX,       			// On input b, on output x
// 				static_cast<int>(_p),          			// LDB
// 				pivot.get(),			// Pivot, entries must be zero to allow pivoting
// 				xerus::EPSILON,      	// Used to determine the accuracy of the Lapacke call. Basically all singular values smaller than RCOND*s[0] are ignored. (s[0] is the largest signular value)
// 				&rank);     			// Outputs the rank of A
			
			IF_CHECK( int lapackAnswer = ) LAPACKE_dgelsd(
				LAPACK_ROW_MAJOR, 
				static_cast<int>(_m),   // Left dimension of A
				static_cast<int>(_n),   // Right dimension of A
				static_cast<int>(_p),	// Number of rhss
				_A,         			// Matrix A
				static_cast<int>(_n),   // LDA
				bOrX,       			// On input b, on output x
				static_cast<int>(_p),	// LDB
				signulars.get(),		// Pivot, entries must be zero to allow pivoting
				xerus::EPSILON,      	// Used to determine the accuracy of the Lapacke call. Basically all singular values smaller than RCOND*s[0] are ignored. (s[0] is the largest signular value)
				&rank);     			// Outputs the rank of A
			
			CHECK(lapackAnswer == 0, error, "Unable to solves min ||Ax - b||_2 for x. Lapacke says: " << lapackAnswer << " sizes are " << _m << " x " << _n << " * " << _p);
			
			if(_m >= _n) { // I.e. bOrX is _b
				misc::copy(_x, bOrX, _n*_p);
			}
			
			XERUS_PA_END("Dense LAPACK", "Solve Least Squares", misc::to_string(_m)+"x"+misc::to_string(_n)+" * "+misc::to_string(_p));
		}
		
	} // namespace blasWrapper

} // namespace xerus

