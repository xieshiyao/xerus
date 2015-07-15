#pragma once

#include "misc/standard.h"

namespace xerus {
    namespace misc {

        template <typename T>
        void array_set_zero(T* const __restrict _x, const size_t _n);

        template <typename T>
        void array_copy(T* const __restrict _x, const T* const _y, const size_t _n);

        void array_set_zero(double* const __restrict _x, const size_t _n);

        void array_copy(double* const __restrict _x, const double* const _y, const size_t _n);

        void array_scaled_copy(double* const __restrict _x, const double _alpha, const double* const _y, const size_t _n);

        void array_scale(double* const __restrict _x, const double _alpha, const size_t _n);

        void array_add(double* const __restrict _x, const double _alpha, const double* const _y, const size_t _n);

        void array_scale_add(const double _alpha, double* const __restrict _x, const double _beta, const double* const _y, const size_t _n);

    }
}

