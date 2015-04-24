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
