#pragma once

#ifndef LAPACK_GLOBAL
#if defined(LAPACK_GLOBAL_PATTERN_LC) || defined(ADD_)
#define LAPACK_GLOBAL(lcname,UCNAME)  lcname##_
#elif defined(LAPACK_GLOBAL_PATTERN_UC) || defined(UPPER)
#define LAPACK_GLOBAL(lcname,UCNAME)  UCNAME
#elif defined(LAPACK_GLOBAL_PATTERN_MC) || defined(NOCHANGE)
#define LAPACK_GLOBAL(lcname,UCNAME)  lcname
#else
#define LAPACK_GLOBAL(lcname,UCNAME)  lcname##_
#endif
#endif

extern "C" {
/*
typedef int lapack_int;

void LAPACK_GLOBAL(dgeqp3, DGEQP3)( lapack_int* m, lapack_int* n, double* a, lapack_int* lda,
									lapack_int* jpvt, double* tau, double* work,
									lapack_int* lwork, lapack_int *info );
	*/
	
}
