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

#include "../../../include/xerus/misc/blasLapackWrapper.h"
#include "../../../include/xerus/misc/selectedFunctions.h"

#include <complex.h>
// Fix for broken complex implementation
#undef I

// Workaround for brocken Lapack
#define lapack_complex_float    float _Complex
#define lapack_complex_double   double _Complex
extern "C"
{
    #include <cblas.h> 
}
#include <lapacke/lapacke.h>

#include <memory>
#include <xerus/misc/standard.h>
#include <xerus/misc/test.h>
 
#ifdef BLAS_ANALYSIS
    #include <map>
    #include <sstream>
    #include "timeMeasure.h"
    #include "stringUtilities.h"
#endif


START_MISC_NAMESPACE

#ifdef BLAS_ANALYSIS
    std::map<std::string, CallCounter> callCounter;
    
    std::string print_blas_analysis() {
        std::stringstream ss;
        ss << "| " << std::endl;
        ss << "| " << std::endl;
        for(const std::pair<std::string, CallCounter>& call : callCounter) {
            ss << "| ======================== " << call.first << " =======================" << std::endl;
            ss << "| Total time " << call.second.totalTime/1000 << " ms in " << call.second.totalCalls << " calls. In detail (showing only calls that contribute more than 1 ms in total)" << std::endl;
            for(const std::pair<std::string, std::pair<size_t, size_t>>& subCall : call.second.calls) {
                if(subCall.second.second/1000 >= 1) {
                    ss << "| | " << std::setfill (' ') << std::setw(5) << subCall.second.first 
                    << " calls to  " << std::setfill (' ') << std::setw(25) << subCall.first 
                    << " taking " << std::setfill (' ') << std::setw(8) << subCall.second.second/1000 
                    << " ms. In average that is " << std::setfill (' ') << std::setw(6) << subCall.second.second/(1000*subCall.second.first) << " ms per call" << std::endl;
                }
            }
            ss << "| " << std::endl;
        }

        ss << "| " << std::endl;
        ss << "| " << std::endl;
        ss << "| Total Blas / Lapack time: " << total_blas_time() << " ms " << std::endl;
        ss << "| " << std::endl;
        ss << "| ";
        return ss.str();
    }
    
    size_t total_blas_time() {
        size_t totalTime = 0;
        for(const std::pair<std::string, CallCounter>& call : callCounter) {
            totalTime += call.second.totalTime;
        }
        return totalTime/1000;
    }
    
    #define START_TIME size_t startTime = uTime()
    
    #define ADD_CALL(name, parameter) add_blas_call(name, parameter, startTime)
    
    void add_blas_call(const std::string& _callName, const std::string& _callParameter, const size_t _startTime) {
        size_t passedTime = uTime()-_startTime;
        CallCounter& call = callCounter[_callName];
        call.totalCalls++;
        call.totalTime += passedTime;
        call.calls[_callParameter].first++;
        call.calls[_callParameter].second += passedTime;
    }
#else 
    std::string print_blas_analysis() {
            return "Blas Analysis is deaktivated";
    }
    
    size_t total_blas_time() {
        return 1;
    }
    
    #define START_TIME void()
    
    #define ADD_CALL(name, parameter) void()
#endif

namespace blasWrapper {
    
    //----------------------------------------------- LEVEL I BLAS ----------------------------------------------------------
    
    double two_norm(const double* const _x, const size_t _n) {
        REQUIRE(_n <= (size_t) std::numeric_limits<int>::max(), "Dimension to large for BLAS/Lapack");
        
        START_TIME;
        
        double result = cblas_dnrm2((int) _n, _x, 1);
        
        ADD_CALL("Two Norm", to_string(_n));
        
        return result;
    }
    
    double dot_product(const double* const _x, const size_t _n, const double* const _y) {
        REQUIRE(_n <= (size_t) std::numeric_limits<int>::max(), "Dimension to large for BLAS/Lapack");
        
        START_TIME;
        
        double result = cblas_ddot((int) _n, _x, 1, _y, 1);
        
        ADD_CALL("Dot Product", to_string(_n)+"*"+to_string(_n));
        
        return result;
    }
    
    
    //----------------------------------------------- LEVEL II BLAS ---------------------------------------------------------
    
    void matrix_vector_product(double* const _x, const size_t _m, const double _alpha, const double* const _A, const size_t _n, const bool _transposed, const double* const _y) {
        // Delegate call if appropriate
        if(_m == 1) { *_x = _alpha*dot_product(_A, _n, _y); return;}
        
        REQUIRE(_m <= (size_t) std::numeric_limits<int>::max(), "Dimension to large for BLAS/Lapack");
        REQUIRE(_n <= (size_t) std::numeric_limits<int>::max(), "Dimension to large for BLAS/Lapack");
        
        START_TIME;
        if(!_transposed) {
            cblas_dgemv(CblasRowMajor, CblasNoTrans, (int)_m, (int)_n, _alpha, _A, (int)_n, _y, 1, 0.0, _x, 1);
        } else {
            cblas_dgemv(CblasRowMajor, CblasTrans, (int)_n, (int)_m, _alpha, _A, (int)_m , _y, 1, 0.0, _x, 1);
        }
        
        ADD_CALL("Matrix Vector Product", to_string(_m)+"x"+to_string(_n)+" * "+to_string(_n));
    }
    
    void dyadic_vector_product(double* _A, const size_t _m, const size_t _n, const double _alpha, const double*const  _x, const double* const _y) {
        REQUIRE(_m <= (size_t) std::numeric_limits<int>::max(), "Dimension to large for BLAS/Lapack");
        REQUIRE(_n <= (size_t) std::numeric_limits<int>::max(), "Dimension to large for BLAS/Lapack");
        
        START_TIME;
        
        //Blas wants to add the product to A, but we don't
        array_set_zero(_A, _m*_n);
        
        cblas_dger(CblasRowMajor, (int)_m, (int)_n, _alpha, _x, 1, _y, 1, _A, (int)_n);
        
        ADD_CALL("Dyadic Vector Product", to_string(_m)+" o "+to_string(_n));
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
        
            REQUIRE(_leftDim <= (size_t) std::numeric_limits<int>::max(), "Dimension to large for BLAS/Lapack");
            REQUIRE(_middleDim <= (size_t) std::numeric_limits<int>::max(), "Dimension to large for BLAS/Lapack");
            REQUIRE(_rightDim <= (size_t) std::numeric_limits<int>::max(), "Dimension to large for BLAS/Lapack");
            REQUIRE(_lda <= (size_t) std::numeric_limits<int>::max(), "Dimension to large for BLAS/Lapack");
            REQUIRE(_ldb <= (size_t) std::numeric_limits<int>::max(), "Dimension to large for BLAS/Lapack");
            
            START_TIME;
            
            cblas_dgemm(    CblasRowMajor,                                  // Array storage format
                            _transposeA ? CblasTrans : CblasNoTrans,        // LHS transposed?
                            _transposeB ? CblasTrans : CblasNoTrans,        // RHS transposed?
                            (int) _leftDim,                                 // Left dimension
                            (int) _rightDim,                                // Right dimension
                            (int) _middleDim,                               // Middle dimension
                            _alpha,                                         // Factor to the Product
                            _A,                                             // Pointer to LHS
                            (int) _lda,                                     // LDA
                            _B,                                             // Pointer to RHS
                            (int) _ldb,                                     // LDB
                            0.0,                                            // Factor of C (Zero if only the product is required)
                            _C,                                             // Pointer to result
                            (int) _rightDim                                 // LDC
                    );
            
            ADD_CALL("Matrix-Matrix-Multiplication", to_string(_leftDim)+"x"+to_string(_middleDim)+" * "+to_string(_middleDim)+"x"+to_string(_rightDim));
        }
    }
    
    
    
    //----------------------------------------------- LAPACK ----------------------------------------------------------------

    void svd( double* const _U, double* const _S, double* const _Vt, const double* const _A, const size_t _m, const size_t _n) {
        //Create copy of A
        const std::unique_ptr<double[]> tmpA(new double[_m*_n]);
        array_copy(tmpA.get(), _A, _m*_n);
        
        svd_destructive(_U, _S, _Vt, tmpA.get(), _m, _n);
    }
    
    void svd_destructive( double* const _U, double* const _S, double* const _Vt, double* const _A, const size_t _m, const size_t _n) {
        REQUIRE(_m <= (size_t) std::numeric_limits<int>::max(), "Dimension to large for BLAS/Lapack");
        REQUIRE(_n <= (size_t) std::numeric_limits<int>::max(), "Dimension to large for BLAS/Lapack");
        
        START_TIME;
        
        int lapackAnswer = LAPACKE_dgesdd(LAPACK_ROW_MAJOR, 'S', (int) _m, (int) _n, _A, (int) _n, _S, _U, (int) std::min(_m, _n), _Vt, (int) _n);
        CHECK(lapackAnswer == 0, error, "Lapack failed to compute SVD. Answer is: " << lapackAnswer);
        
        ADD_CALL("Singular Value Decomposition", to_string(_m)+"x"+to_string(_n));
    }
    
    
    
    void qr( double* const _Q, double* const _R, const double* const _A, const size_t _m, const size_t _n) {
        // Create tmp copy of A since Lapack wants to destroy it
        const std::unique_ptr<double[]> tmpA(new double[_m*_n]);
        array_copy(tmpA.get(), _A, _m*_n);
        
        qr_destructive(_Q, _R, tmpA.get(), _m, _n);
    }
    
    void inplace_qr(double* const _AtoQ, double* const _R, const size_t _m, const size_t _n) {
        qr_destructive(_AtoQ, _R, _AtoQ, _m, _n);
    }
    
    void qr_destructive( double* const _Q, double* const _R, double* const _A, const size_t _m, const size_t _n) {
        REQUIRE(_m <= (size_t) std::numeric_limits<int>::max(), "Dimension to large for BLAS/Lapack");
        REQUIRE(_n <= (size_t) std::numeric_limits<int>::max(), "Dimension to large for BLAS/Lapack");
        
        REQUIRE(_n > 0, "Dimension n must be larger than zero");
        REQUIRE(_m > 0, "Dimension m must be larger than zero");
        
        REQUIRE(_Q && _R && _A, "QR decomposition must not be called with null pointers: Q:" << _Q << " R: " << _R << " A: " << _A);
        REQUIRE(_A != _R, "_A and _R must be different, otherwise qr call will fail.");
        
        START_TIME;
        
        // Maximal rank is used by Lapacke
        const size_t rank = std::min(_m, _n); 
        
        // Tmp Array for Lapacke
        const std::unique_ptr<double[]> tau(new double[rank]);
        
        // Calculate QR factorisations
//         LOG(Lapacke, "Call to dorgqr with parameters: " << LAPACK_ROW_MAJOR << ", " << (int) _m  << ", " << (int) _n  << ", " << _A << ", " << (int) _n  << ", " << tau.get());
        int lapackAnswer = LAPACKE_dgeqrf(LAPACK_ROW_MAJOR, (int) _m, (int) _n, _A, (int) _n, tau.get());
        CHECK(lapackAnswer == 0, error, "Unable to perform QR factorisaton. Lapacke says: " << lapackAnswer );
        
        // Copy the upper triangular Matrix R (rank x _n) into position
        for(size_t row =0; row < rank; ++row) {
            array_set_zero(_R+row*_n, row); // Set starting zeros
            array_copy(_R+row*_n+row, _A+row*_n+row, _n-row); // Copy upper triangular part from lapack result.
        }
        
        // Create orthogonal matrix Q (in tmpA)
//         LOG(Lapacke, "Call to dorgqr with parameters: " << LAPACK_ROW_MAJOR << ", " << (int) _m  << ", " << (int) rank  << ", " << (int) rank << ", " << _A << ", " << (int) _n  << ", " << tau.get());
        lapackAnswer = LAPACKE_dorgqr(LAPACK_ROW_MAJOR, (int) _m, (int) rank, (int) rank, _A, (int) _n, tau.get());
        CHECK(lapackAnswer == 0, error, "Unable to reconstruct Q from the QR factorisation. Lapacke says: " << lapackAnswer);
        
        //Copy Q (_m x rank) into position, if Q is not to be constructed in place of A
        if(_A != _Q) {
            if(_m == _n) {
                array_copy(_Q, _A, _m*_n);
            } else {
                for(size_t row =0; row < _m; ++row) {
                    array_copy(_Q+row*rank, _A+row*_n, rank);
                }
            }
        }
        
        ADD_CALL("QR Factorisation", to_string(_m)+"x"+to_string(_n));
    }
    
    

    void rq( double* const _R, double* const _Q, const double* const _A, const size_t _m, const size_t _n) {
        // Create tmp copy of A since Lapack wants to destroy it
        const std::unique_ptr<double[]> tmpA(new double[_m*_n]);
        array_copy(tmpA.get(), _A, _m*_n);
        
        rq_destructive(_R, _Q, tmpA.get(), _m, _n);
    }
    
    void inplace_rq( double* const _R, double* const _AtoQ, const size_t _m, const size_t _n) {
        rq_destructive(_R, _AtoQ, _AtoQ, _m, _n);
    }
    
    
    void rq_destructive( double* const _R, double* const _Q, double* const _A, const size_t _m, const size_t _n) {
        REQUIRE(_m <= (size_t) std::numeric_limits<int>::max(), "Dimension to large for BLAS/Lapack");
        REQUIRE(_n <= (size_t) std::numeric_limits<int>::max(), "Dimension to large for BLAS/Lapack");
        
        REQUIRE(_n > 0, "Dimension n must be larger than zero");
        REQUIRE(_m > 0, "Dimension m must be larger than zero");
        
        REQUIRE(_Q && _R && _A, "QR decomposition must not be called with null pointers: R " << _R << " Q: " << _Q << " A: " << _A);
        REQUIRE(_A != _R, "_A and _R must be different, otherwise qr call will fail.");
        
        START_TIME;
        
        // Maximal rank is used by Lapacke
        const size_t rank = std::min(_m, _n); 
        
        // Tmp Array for Lapacke
        const std::unique_ptr<double[]> tau(new double[rank]);
        
        int lapackAnswer = LAPACKE_dgerqf(LAPACK_ROW_MAJOR, (int) _m, (int) _n, _A, (int) _n, tau.get());
        CHECK(lapackAnswer == 0, error, "Unable to perform QR factorisaton. Lapacke says: " << lapackAnswer );
        
        
        // Copy the upper triangular Matrix R (_m x rank) into position.
        size_t row = 0;
        for( ; row < _m - rank; ++row) {
            array_copy(_R+row*rank, _A+row*_n+_n-rank, rank);
        }
        for(size_t skip = 0; row < _m; ++row, ++skip) {
            array_set_zero(_R+row*rank, skip); // Set starting zeros
            array_copy(_R+row*rank+skip, _A+row*_n+_n-rank+skip, rank-skip); // Copy upper triangular part from lapack result.
        }
        
        // Create orthogonal matrix Q (in _A). Lapacke expects to get the last rank rows of A...
        lapackAnswer = LAPACKE_dorgrq(LAPACK_ROW_MAJOR, (int) rank, (int) _n, (int) rank, _A+(_m-rank)*_n, (int) _n, tau.get()); 
        CHECK(lapackAnswer == 0, error, "Unable to reconstruct Q from the RQ factorisation. Lapacke says: " << lapackAnswer);

        
        //Copy Q (rank x _n) into position
        if(_A != _Q) {
            array_copy(_Q, _A+(_m-rank)*_n, rank*_n);
        }
        
        ADD_CALL("RQ Factorisation", to_string(_m)+"x"+to_string(_n));
    }
    
    
/*  TODO we need test cases for these  
    /// Solves Ax = b for x
    void solve( double* const _x, const double* const _A, const size_t _n, const double* const _b) {
        const std::unique_ptr<double[]> tmpA(new double[_n*_n]);
        array_copy(tmpA.get(), _A, _n*_n);
        array_copy(_x, _b, _n);
        
        solve_destructive(_x, tmpA.get(), _n);
    }
    
    /// Solves Ax = b for x
    void solve_destructive( double* const _bToX, double* const _A, const size_t _n) {
        REQUIRE(_n <= (size_t) std::numeric_limits<int>::max(), "Dimension to large for BLAS/Lapack");
        
        START_TIME;
        
        std::unique_ptr<int[]> pivot(new int[_n]);
        
        int lapackAnswer = LAPACKE_dgesv(
            LAPACK_ROW_MAJOR,
            (int) _n,       // Dimensions of A (nxn)
            1,              // Number of b's, here always one
            _A,             // The input matrix A, will be destroyed
            (int) _n,       // LDA
            pivot.get(),    // Output of the pivot ordering
            _bToX,          // Input the vector b, output the vector x
            1 );
        CHECK(lapackAnswer == 0, error, "Unable to solves Ax = b. Lapacke says: " << lapackAnswer);
        
        ADD_CALL("Solve ", to_string(_n)+"x"+to_string(_n));
    }*/
    
    /// Solves min ||Ax - b||_2 for x
    void solve_least_squares( double* const _x, const double* const _A, const size_t _m, const size_t _n, const double* const _b){
        const std::unique_ptr<double[]> tmpA(new double[_m*_n]);
        array_copy(tmpA.get(), _A, _m*_n);
        
        const std::unique_ptr<double[]> tmpB(new double[_n]);
        array_copy(tmpB.get(), _b, _n);
        
        solve_least_squares_destructive(_x, tmpA.get(), _m, _n, tmpB.get());
    }
    
    /// Solves min ||Ax - b||_2 for x
    void solve_least_squares_destructive( double* const _x, double* const _A, const size_t _m, const size_t _n, double* const _b){
        REQUIRE(_m <= (size_t) std::numeric_limits<int>::max(), "Dimension to large for BLAS/Lapack");
        REQUIRE(_n <= (size_t) std::numeric_limits<int>::max(), "Dimension to large for BLAS/Lapack");
        
        START_TIME;
        
        std::unique_ptr<int[]> pivot(new int[_n]);
        array_set_zero(pivot.get(), _n);
        int rank;
        
        double* bOrX;
        if(_m >= _n) {
            bOrX = _x;
            array_copy(bOrX, _b, _n);
        } else {
            bOrX = _b;
        }
        
        int lapackAnswer = LAPACKE_dgelsy(
            LAPACK_ROW_MAJOR, 
            (int) _m,   // Left dimension of A
            (int) _n,   // Right dimension of A
            1,          // Number of b's, here always one
            _A,         // Matrix A
            (int) _n,   // LDA
            bOrX,       // On input b, on output x
            1,          // LDB, here always one
            pivot.get(),// Pivot, entries must be zero to allow pivoting
            1e-12,      // Used to determine the accuracy of the Lapacke call. Basically all singular values smaller than RCOND*s[0] are ignored. (s[0] is the largest signular value)
            &rank);     // Outputs the rank of A
        CHECK(lapackAnswer == 0, error, "Unable to solves min ||Ax - b||_2 for x. Lapacke says: " << lapackAnswer);
        
        if(_m < _n) {
            array_copy(_x, bOrX, _m);
        }
        
        ADD_CALL("Solve Least Squares", to_string(_m)+"x"+to_string(_n));
    }
    
}

END_MISC_NAMESPACE
