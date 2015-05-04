#pragma once

#include <cstring>
#include "blasLapackWrapper.h"


START_MISC_NAMESPACE

template <typename T>
_inline_ void array_set_zero(T* const __restrict _x, const size_t _n) {
    memset(_x, 0, _n*sizeof(T));
}

template <typename T>
_inline_ void array_copy(T* const __restrict _x, const T* const _y, const size_t _n) {
    memcpy(_x, _y, _n*sizeof(T));
}

 _inline_ void array_set_zero(double* const __restrict _x, const size_t _n) {
	for(size_t i=0; i<_n; ++i) { _x[i] = 0; }
}

_inline_ void array_copy(double* const __restrict _x, const double* const _y, const size_t _n) {
	memcpy(_x, _y, _n*sizeof(double));
}

_inline_ void array_scaled_copy(double* const __restrict _x, const double _alpha, const double* const _y, const size_t _n) {
	for(size_t i=0; i<_n; ++i) { _x[i] = _alpha*_y[i]; }
}

_inline_ void array_scale(double* const __restrict _x, const double _alpha, const size_t _n) {
	REQUIRE(_n <= (size_t) std::numeric_limits<int>::max(), "Dimension to large for blas/lapack");
	cblas_dscal((int) _n, _alpha, _x, 1);
}

_inline_ void array_add(double* const __restrict _x, const double _alpha, const double* const _y, const size_t _n) {
	REQUIRE(_n <= (size_t) std::numeric_limits<int>::max(), "Dimension to large for blas/lapack");
	cblas_daxpy((int) _n, _alpha, _y, 1, _x, 1);
}

_inline_ void array_scale_add(const double _alpha, double* const __restrict _x, const double _beta, const double* const _y, const size_t _n) {
	for(size_t i = 0; i < _n; i++ ) { _x[i] = _alpha*_x[i]+_beta*_y[i]; }
}

END_MISC_NAMESPACE

