/*****************************************************************************
  Copyright (c) 2010, Intel Corp.
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of Intel Corporation nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
  THE POSSIBILITY OF SUCH DAMAGE.
******************************************************************************
* Contents: Native C interface to LAPACK
* Author: Intel Corporation
* Generated November, 2011
* Modified June 2015 to meet the needs of the xerus library.
*****************************************************************************/

/**
 * @file
 * @brief Header file for all lapack functions defined in the lapack fortran library (version 3.5.0) to allow direct calls to these fortran functions.
 */

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

typedef int lapack_int;
typedef int lapack_logical;
typedef std::complex<double> lapack_complex_double;
typedef std::complex<float> lapack_complex_float;

extern "C" {

#define LAPACK_sgetrf LAPACK_GLOBAL(sgetrf,SGETRF)
#define LAPACK_dgetrf LAPACK_GLOBAL(dgetrf,DGETRF)
#define LAPACK_cgetrf LAPACK_GLOBAL(cgetrf,CGETRF)
#define LAPACK_zgetrf LAPACK_GLOBAL(zgetrf,ZGETRF)
#define LAPACK_sgbtrf LAPACK_GLOBAL(sgbtrf,SGBTRF)
#define LAPACK_dgbtrf LAPACK_GLOBAL(dgbtrf,DGBTRF)
#define LAPACK_cgbtrf LAPACK_GLOBAL(cgbtrf,CGBTRF)
#define LAPACK_zgbtrf LAPACK_GLOBAL(zgbtrf,ZGBTRF)
#define LAPACK_sgttrf LAPACK_GLOBAL(sgttrf,SGTTRF)
#define LAPACK_dgttrf LAPACK_GLOBAL(dgttrf,DGTTRF)
#define LAPACK_cgttrf LAPACK_GLOBAL(cgttrf,CGTTRF)
#define LAPACK_zgttrf LAPACK_GLOBAL(zgttrf,ZGTTRF)
#define LAPACK_spotrf LAPACK_GLOBAL(spotrf,SPOTRF)
#define LAPACK_dpotrf LAPACK_GLOBAL(dpotrf,DPOTRF)
#define LAPACK_cpotrf LAPACK_GLOBAL(cpotrf,CPOTRF)
#define LAPACK_zpotrf LAPACK_GLOBAL(zpotrf,ZPOTRF)
#define LAPACK_dpstrf LAPACK_GLOBAL(dpstrf,DPSTRF)
#define LAPACK_spstrf LAPACK_GLOBAL(spstrf,SPSTRF)
#define LAPACK_zpstrf LAPACK_GLOBAL(zpstrf,ZPSTRF)
#define LAPACK_cpstrf LAPACK_GLOBAL(cpstrf,CPSTRF)
#define LAPACK_dpftrf LAPACK_GLOBAL(dpftrf,DPFTRF)
#define LAPACK_spftrf LAPACK_GLOBAL(spftrf,SPFTRF)
#define LAPACK_zpftrf LAPACK_GLOBAL(zpftrf,ZPFTRF)
#define LAPACK_cpftrf LAPACK_GLOBAL(cpftrf,CPFTRF)
#define LAPACK_spptrf LAPACK_GLOBAL(spptrf,SPPTRF)
#define LAPACK_dpptrf LAPACK_GLOBAL(dpptrf,DPPTRF)
#define LAPACK_cpptrf LAPACK_GLOBAL(cpptrf,CPPTRF)
#define LAPACK_zpptrf LAPACK_GLOBAL(zpptrf,ZPPTRF)
#define LAPACK_spbtrf LAPACK_GLOBAL(spbtrf,SPBTRF)
#define LAPACK_dpbtrf LAPACK_GLOBAL(dpbtrf,DPBTRF)
#define LAPACK_cpbtrf LAPACK_GLOBAL(cpbtrf,CPBTRF)
#define LAPACK_zpbtrf LAPACK_GLOBAL(zpbtrf,ZPBTRF)
#define LAPACK_spttrf LAPACK_GLOBAL(spttrf,SPTTRF)
#define LAPACK_dpttrf LAPACK_GLOBAL(dpttrf,DPTTRF)
#define LAPACK_cpttrf LAPACK_GLOBAL(cpttrf,CPTTRF)
#define LAPACK_zpttrf LAPACK_GLOBAL(zpttrf,ZPTTRF)
#define LAPACK_ssytrf LAPACK_GLOBAL(ssytrf,SSYTRF)
#define LAPACK_dsytrf LAPACK_GLOBAL(dsytrf,DSYTRF)
#define LAPACK_csytrf LAPACK_GLOBAL(csytrf,CSYTRF)
#define LAPACK_zsytrf LAPACK_GLOBAL(zsytrf,ZSYTRF)
#define LAPACK_chetrf LAPACK_GLOBAL(chetrf,CHETRF)
#define LAPACK_zhetrf LAPACK_GLOBAL(zhetrf,ZHETRF)
#define LAPACK_ssptrf LAPACK_GLOBAL(ssptrf,SSPTRF)
#define LAPACK_dsptrf LAPACK_GLOBAL(dsptrf,DSPTRF)
#define LAPACK_csptrf LAPACK_GLOBAL(csptrf,CSPTRF)
#define LAPACK_zsptrf LAPACK_GLOBAL(zsptrf,ZSPTRF)
#define LAPACK_chptrf LAPACK_GLOBAL(chptrf,CHPTRF)
#define LAPACK_zhptrf LAPACK_GLOBAL(zhptrf,ZHPTRF)
#define LAPACK_sgetrs LAPACK_GLOBAL(sgetrs,SGETRS)
#define LAPACK_dgetrs LAPACK_GLOBAL(dgetrs,DGETRS)
#define LAPACK_cgetrs LAPACK_GLOBAL(cgetrs,CGETRS)
#define LAPACK_zgetrs LAPACK_GLOBAL(zgetrs,ZGETRS)
#define LAPACK_sgbtrs LAPACK_GLOBAL(sgbtrs,SGBTRS)
#define LAPACK_dgbtrs LAPACK_GLOBAL(dgbtrs,DGBTRS)
#define LAPACK_cgbtrs LAPACK_GLOBAL(cgbtrs,CGBTRS)
#define LAPACK_zgbtrs LAPACK_GLOBAL(zgbtrs,ZGBTRS)
#define LAPACK_sgttrs LAPACK_GLOBAL(sgttrs,SGTTRS)
#define LAPACK_dgttrs LAPACK_GLOBAL(dgttrs,DGTTRS)
#define LAPACK_cgttrs LAPACK_GLOBAL(cgttrs,CGTTRS)
#define LAPACK_zgttrs LAPACK_GLOBAL(zgttrs,ZGTTRS)
#define LAPACK_spotrs LAPACK_GLOBAL(spotrs,SPOTRS)
#define LAPACK_dpotrs LAPACK_GLOBAL(dpotrs,DPOTRS)
#define LAPACK_cpotrs LAPACK_GLOBAL(cpotrs,CPOTRS)
#define LAPACK_zpotrs LAPACK_GLOBAL(zpotrs,ZPOTRS)
#define LAPACK_dpftrs LAPACK_GLOBAL(dpftrs,DPFTRS)
#define LAPACK_spftrs LAPACK_GLOBAL(spftrs,SPFTRS)
#define LAPACK_zpftrs LAPACK_GLOBAL(zpftrs,ZPFTRS)
#define LAPACK_cpftrs LAPACK_GLOBAL(cpftrs,CPFTRS)
#define LAPACK_spptrs LAPACK_GLOBAL(spptrs,SPPTRS)
#define LAPACK_dpptrs LAPACK_GLOBAL(dpptrs,DPPTRS)
#define LAPACK_cpptrs LAPACK_GLOBAL(cpptrs,CPPTRS)
#define LAPACK_zpptrs LAPACK_GLOBAL(zpptrs,ZPPTRS)
#define LAPACK_spbtrs LAPACK_GLOBAL(spbtrs,SPBTRS)
#define LAPACK_dpbtrs LAPACK_GLOBAL(dpbtrs,DPBTRS)
#define LAPACK_cpbtrs LAPACK_GLOBAL(cpbtrs,CPBTRS)
#define LAPACK_zpbtrs LAPACK_GLOBAL(zpbtrs,ZPBTRS)
#define LAPACK_spttrs LAPACK_GLOBAL(spttrs,SPTTRS)
#define LAPACK_dpttrs LAPACK_GLOBAL(dpttrs,DPTTRS)
#define LAPACK_cpttrs LAPACK_GLOBAL(cpttrs,CPTTRS)
#define LAPACK_zpttrs LAPACK_GLOBAL(zpttrs,ZPTTRS)
#define LAPACK_ssytrs LAPACK_GLOBAL(ssytrs,SSYTRS)
#define LAPACK_dsytrs LAPACK_GLOBAL(dsytrs,DSYTRS)
#define LAPACK_csytrs LAPACK_GLOBAL(csytrs,CSYTRS)
#define LAPACK_zsytrs LAPACK_GLOBAL(zsytrs,ZSYTRS)
#define LAPACK_chetrs LAPACK_GLOBAL(chetrs,CHETRS)
#define LAPACK_zhetrs LAPACK_GLOBAL(zhetrs,ZHETRS)
#define LAPACK_ssptrs LAPACK_GLOBAL(ssptrs,SSPTRS)
#define LAPACK_dsptrs LAPACK_GLOBAL(dsptrs,DSPTRS)
#define LAPACK_csptrs LAPACK_GLOBAL(csptrs,CSPTRS)
#define LAPACK_zsptrs LAPACK_GLOBAL(zsptrs,ZSPTRS)
#define LAPACK_chptrs LAPACK_GLOBAL(chptrs,CHPTRS)
#define LAPACK_zhptrs LAPACK_GLOBAL(zhptrs,ZHPTRS)
#define LAPACK_strtrs LAPACK_GLOBAL(strtrs,STRTRS)
#define LAPACK_dtrtrs LAPACK_GLOBAL(dtrtrs,DTRTRS)
#define LAPACK_ctrtrs LAPACK_GLOBAL(ctrtrs,CTRTRS)
#define LAPACK_ztrtrs LAPACK_GLOBAL(ztrtrs,ZTRTRS)
#define LAPACK_stptrs LAPACK_GLOBAL(stptrs,STPTRS)
#define LAPACK_dtptrs LAPACK_GLOBAL(dtptrs,DTPTRS)
#define LAPACK_ctptrs LAPACK_GLOBAL(ctptrs,CTPTRS)
#define LAPACK_ztptrs LAPACK_GLOBAL(ztptrs,ZTPTRS)
#define LAPACK_stbtrs LAPACK_GLOBAL(stbtrs,STBTRS)
#define LAPACK_dtbtrs LAPACK_GLOBAL(dtbtrs,DTBTRS)
#define LAPACK_ctbtrs LAPACK_GLOBAL(ctbtrs,CTBTRS)
#define LAPACK_ztbtrs LAPACK_GLOBAL(ztbtrs,ZTBTRS)
#define LAPACK_sgecon LAPACK_GLOBAL(sgecon,SGECON)
#define LAPACK_dgecon LAPACK_GLOBAL(dgecon,DGECON)
#define LAPACK_cgecon LAPACK_GLOBAL(cgecon,CGECON)
#define LAPACK_zgecon LAPACK_GLOBAL(zgecon,ZGECON)
#define LAPACK_sgbcon LAPACK_GLOBAL(sgbcon,SGBCON)
#define LAPACK_dgbcon LAPACK_GLOBAL(dgbcon,DGBCON)
#define LAPACK_cgbcon LAPACK_GLOBAL(cgbcon,CGBCON)
#define LAPACK_zgbcon LAPACK_GLOBAL(zgbcon,ZGBCON)
#define LAPACK_sgtcon LAPACK_GLOBAL(sgtcon,SGTCON)
#define LAPACK_dgtcon LAPACK_GLOBAL(dgtcon,DGTCON)
#define LAPACK_cgtcon LAPACK_GLOBAL(cgtcon,CGTCON)
#define LAPACK_zgtcon LAPACK_GLOBAL(zgtcon,ZGTCON)
#define LAPACK_spocon LAPACK_GLOBAL(spocon,SPOCON)
#define LAPACK_dpocon LAPACK_GLOBAL(dpocon,DPOCON)
#define LAPACK_cpocon LAPACK_GLOBAL(cpocon,CPOCON)
#define LAPACK_zpocon LAPACK_GLOBAL(zpocon,ZPOCON)
#define LAPACK_sppcon LAPACK_GLOBAL(sppcon,SPPCON)
#define LAPACK_dppcon LAPACK_GLOBAL(dppcon,DPPCON)
#define LAPACK_cppcon LAPACK_GLOBAL(cppcon,CPPCON)
#define LAPACK_zppcon LAPACK_GLOBAL(zppcon,ZPPCON)
#define LAPACK_spbcon LAPACK_GLOBAL(spbcon,SPBCON)
#define LAPACK_dpbcon LAPACK_GLOBAL(dpbcon,DPBCON)
#define LAPACK_cpbcon LAPACK_GLOBAL(cpbcon,CPBCON)
#define LAPACK_zpbcon LAPACK_GLOBAL(zpbcon,ZPBCON)
#define LAPACK_sptcon LAPACK_GLOBAL(sptcon,SPTCON)
#define LAPACK_dptcon LAPACK_GLOBAL(dptcon,DPTCON)
#define LAPACK_cptcon LAPACK_GLOBAL(cptcon,CPTCON)
#define LAPACK_zptcon LAPACK_GLOBAL(zptcon,ZPTCON)
#define LAPACK_ssycon LAPACK_GLOBAL(ssycon,SSYCON)
#define LAPACK_dsycon LAPACK_GLOBAL(dsycon,DSYCON)
#define LAPACK_csycon LAPACK_GLOBAL(csycon,CSYCON)
#define LAPACK_zsycon LAPACK_GLOBAL(zsycon,ZSYCON)
#define LAPACK_checon LAPACK_GLOBAL(checon,CHECON)
#define LAPACK_zhecon LAPACK_GLOBAL(zhecon,ZHECON)
#define LAPACK_sspcon LAPACK_GLOBAL(sspcon,SSPCON)
#define LAPACK_dspcon LAPACK_GLOBAL(dspcon,DSPCON)
#define LAPACK_cspcon LAPACK_GLOBAL(cspcon,CSPCON)
#define LAPACK_zspcon LAPACK_GLOBAL(zspcon,ZSPCON)
#define LAPACK_chpcon LAPACK_GLOBAL(chpcon,CHPCON)
#define LAPACK_zhpcon LAPACK_GLOBAL(zhpcon,ZHPCON)
#define LAPACK_strcon LAPACK_GLOBAL(strcon,STRCON)
#define LAPACK_dtrcon LAPACK_GLOBAL(dtrcon,DTRCON)
#define LAPACK_ctrcon LAPACK_GLOBAL(ctrcon,CTRCON)
#define LAPACK_ztrcon LAPACK_GLOBAL(ztrcon,ZTRCON)
#define LAPACK_stpcon LAPACK_GLOBAL(stpcon,STPCON)
#define LAPACK_dtpcon LAPACK_GLOBAL(dtpcon,DTPCON)
#define LAPACK_ctpcon LAPACK_GLOBAL(ctpcon,CTPCON)
#define LAPACK_ztpcon LAPACK_GLOBAL(ztpcon,ZTPCON)
#define LAPACK_stbcon LAPACK_GLOBAL(stbcon,STBCON)
#define LAPACK_dtbcon LAPACK_GLOBAL(dtbcon,DTBCON)
#define LAPACK_ctbcon LAPACK_GLOBAL(ctbcon,CTBCON)
#define LAPACK_ztbcon LAPACK_GLOBAL(ztbcon,ZTBCON)
#define LAPACK_sgerfs LAPACK_GLOBAL(sgerfs,SGERFS)
#define LAPACK_dgerfs LAPACK_GLOBAL(dgerfs,DGERFS)
#define LAPACK_cgerfs LAPACK_GLOBAL(cgerfs,CGERFS)
#define LAPACK_zgerfs LAPACK_GLOBAL(zgerfs,ZGERFS)
#define LAPACK_dgerfsx LAPACK_GLOBAL(dgerfsx,DGERFSX)
#define LAPACK_sgerfsx LAPACK_GLOBAL(sgerfsx,SGERFSX)
#define LAPACK_zgerfsx LAPACK_GLOBAL(zgerfsx,ZGERFSX)
#define LAPACK_cgerfsx LAPACK_GLOBAL(cgerfsx,CGERFSX)
#define LAPACK_sgbrfs LAPACK_GLOBAL(sgbrfs,SGBRFS)
#define LAPACK_dgbrfs LAPACK_GLOBAL(dgbrfs,DGBRFS)
#define LAPACK_cgbrfs LAPACK_GLOBAL(cgbrfs,CGBRFS)
#define LAPACK_zgbrfs LAPACK_GLOBAL(zgbrfs,ZGBRFS)
#define LAPACK_dgbrfsx LAPACK_GLOBAL(dgbrfsx,DGBRFSX)
#define LAPACK_sgbrfsx LAPACK_GLOBAL(sgbrfsx,SGBRFSX)
#define LAPACK_zgbrfsx LAPACK_GLOBAL(zgbrfsx,ZGBRFSX)
#define LAPACK_cgbrfsx LAPACK_GLOBAL(cgbrfsx,CGBRFSX)
#define LAPACK_sgtrfs LAPACK_GLOBAL(sgtrfs,SGTRFS)
#define LAPACK_dgtrfs LAPACK_GLOBAL(dgtrfs,DGTRFS)
#define LAPACK_cgtrfs LAPACK_GLOBAL(cgtrfs,CGTRFS)
#define LAPACK_zgtrfs LAPACK_GLOBAL(zgtrfs,ZGTRFS)
#define LAPACK_sporfs LAPACK_GLOBAL(sporfs,SPORFS)
#define LAPACK_dporfs LAPACK_GLOBAL(dporfs,DPORFS)
#define LAPACK_cporfs LAPACK_GLOBAL(cporfs,CPORFS)
#define LAPACK_zporfs LAPACK_GLOBAL(zporfs,ZPORFS)
#define LAPACK_dporfsx LAPACK_GLOBAL(dporfsx,DPORFSX)
#define LAPACK_sporfsx LAPACK_GLOBAL(sporfsx,SPORFSX)
#define LAPACK_zporfsx LAPACK_GLOBAL(zporfsx,ZPORFSX)
#define LAPACK_cporfsx LAPACK_GLOBAL(cporfsx,CPORFSX)
#define LAPACK_spprfs LAPACK_GLOBAL(spprfs,SPPRFS)
#define LAPACK_dpprfs LAPACK_GLOBAL(dpprfs,DPPRFS)
#define LAPACK_cpprfs LAPACK_GLOBAL(cpprfs,CPPRFS)
#define LAPACK_zpprfs LAPACK_GLOBAL(zpprfs,ZPPRFS)
#define LAPACK_spbrfs LAPACK_GLOBAL(spbrfs,SPBRFS)
#define LAPACK_dpbrfs LAPACK_GLOBAL(dpbrfs,DPBRFS)
#define LAPACK_cpbrfs LAPACK_GLOBAL(cpbrfs,CPBRFS)
#define LAPACK_zpbrfs LAPACK_GLOBAL(zpbrfs,ZPBRFS)
#define LAPACK_sptrfs LAPACK_GLOBAL(sptrfs,SPTRFS)
#define LAPACK_dptrfs LAPACK_GLOBAL(dptrfs,DPTRFS)
#define LAPACK_cptrfs LAPACK_GLOBAL(cptrfs,CPTRFS)
#define LAPACK_zptrfs LAPACK_GLOBAL(zptrfs,ZPTRFS)
#define LAPACK_ssyrfs LAPACK_GLOBAL(ssyrfs,SSYRFS)
#define LAPACK_dsyrfs LAPACK_GLOBAL(dsyrfs,DSYRFS)
#define LAPACK_csyrfs LAPACK_GLOBAL(csyrfs,CSYRFS)
#define LAPACK_zsyrfs LAPACK_GLOBAL(zsyrfs,ZSYRFS)
#define LAPACK_dsyrfsx LAPACK_GLOBAL(dsyrfsx,DSYRFSX)
#define LAPACK_ssyrfsx LAPACK_GLOBAL(ssyrfsx,SSYRFSX)
#define LAPACK_zsyrfsx LAPACK_GLOBAL(zsyrfsx,ZSYRFSX)
#define LAPACK_csyrfsx LAPACK_GLOBAL(csyrfsx,CSYRFSX)
#define LAPACK_cherfs LAPACK_GLOBAL(cherfs,CHERFS)
#define LAPACK_zherfs LAPACK_GLOBAL(zherfs,ZHERFS)
#define LAPACK_zherfsx LAPACK_GLOBAL(zherfsx,ZHERFSX)
#define LAPACK_cherfsx LAPACK_GLOBAL(cherfsx,CHERFSX)
#define LAPACK_ssprfs LAPACK_GLOBAL(ssprfs,SSPRFS)
#define LAPACK_dsprfs LAPACK_GLOBAL(dsprfs,DSPRFS)
#define LAPACK_csprfs LAPACK_GLOBAL(csprfs,CSPRFS)
#define LAPACK_zsprfs LAPACK_GLOBAL(zsprfs,ZSPRFS)
#define LAPACK_chprfs LAPACK_GLOBAL(chprfs,CHPRFS)
#define LAPACK_zhprfs LAPACK_GLOBAL(zhprfs,ZHPRFS)
#define LAPACK_strrfs LAPACK_GLOBAL(strrfs,STRRFS)
#define LAPACK_dtrrfs LAPACK_GLOBAL(dtrrfs,DTRRFS)
#define LAPACK_ctrrfs LAPACK_GLOBAL(ctrrfs,CTRRFS)
#define LAPACK_ztrrfs LAPACK_GLOBAL(ztrrfs,ZTRRFS)
#define LAPACK_stprfs LAPACK_GLOBAL(stprfs,STPRFS)
#define LAPACK_dtprfs LAPACK_GLOBAL(dtprfs,DTPRFS)
#define LAPACK_ctprfs LAPACK_GLOBAL(ctprfs,CTPRFS)
#define LAPACK_ztprfs LAPACK_GLOBAL(ztprfs,ZTPRFS)
#define LAPACK_stbrfs LAPACK_GLOBAL(stbrfs,STBRFS)
#define LAPACK_dtbrfs LAPACK_GLOBAL(dtbrfs,DTBRFS)
#define LAPACK_ctbrfs LAPACK_GLOBAL(ctbrfs,CTBRFS)
#define LAPACK_ztbrfs LAPACK_GLOBAL(ztbrfs,ZTBRFS)
#define LAPACK_sgetri LAPACK_GLOBAL(sgetri,SGETRI)
#define LAPACK_dgetri LAPACK_GLOBAL(dgetri,DGETRI)
#define LAPACK_cgetri LAPACK_GLOBAL(cgetri,CGETRI)
#define LAPACK_zgetri LAPACK_GLOBAL(zgetri,ZGETRI)
#define LAPACK_spotri LAPACK_GLOBAL(spotri,SPOTRI)
#define LAPACK_dpotri LAPACK_GLOBAL(dpotri,DPOTRI)
#define LAPACK_cpotri LAPACK_GLOBAL(cpotri,CPOTRI)
#define LAPACK_zpotri LAPACK_GLOBAL(zpotri,ZPOTRI)
#define LAPACK_dpftri LAPACK_GLOBAL(dpftri,DPFTRI)
#define LAPACK_spftri LAPACK_GLOBAL(spftri,SPFTRI)
#define LAPACK_zpftri LAPACK_GLOBAL(zpftri,ZPFTRI)
#define LAPACK_cpftri LAPACK_GLOBAL(cpftri,CPFTRI)
#define LAPACK_spptri LAPACK_GLOBAL(spptri,SPPTRI)
#define LAPACK_dpptri LAPACK_GLOBAL(dpptri,DPPTRI)
#define LAPACK_cpptri LAPACK_GLOBAL(cpptri,CPPTRI)
#define LAPACK_zpptri LAPACK_GLOBAL(zpptri,ZPPTRI)
#define LAPACK_ssytri LAPACK_GLOBAL(ssytri,SSYTRI)
#define LAPACK_dsytri LAPACK_GLOBAL(dsytri,DSYTRI)
#define LAPACK_csytri LAPACK_GLOBAL(csytri,CSYTRI)
#define LAPACK_zsytri LAPACK_GLOBAL(zsytri,ZSYTRI)
#define LAPACK_chetri LAPACK_GLOBAL(chetri,CHETRI)
#define LAPACK_zhetri LAPACK_GLOBAL(zhetri,ZHETRI)
#define LAPACK_ssptri LAPACK_GLOBAL(ssptri,SSPTRI)
#define LAPACK_dsptri LAPACK_GLOBAL(dsptri,DSPTRI)
#define LAPACK_csptri LAPACK_GLOBAL(csptri,CSPTRI)
#define LAPACK_zsptri LAPACK_GLOBAL(zsptri,ZSPTRI)
#define LAPACK_chptri LAPACK_GLOBAL(chptri,CHPTRI)
#define LAPACK_zhptri LAPACK_GLOBAL(zhptri,ZHPTRI)
#define LAPACK_strtri LAPACK_GLOBAL(strtri,STRTRI)
#define LAPACK_dtrtri LAPACK_GLOBAL(dtrtri,DTRTRI)
#define LAPACK_ctrtri LAPACK_GLOBAL(ctrtri,CTRTRI)
#define LAPACK_ztrtri LAPACK_GLOBAL(ztrtri,ZTRTRI)
#define LAPACK_dtftri LAPACK_GLOBAL(dtftri,DTFTRI)
#define LAPACK_stftri LAPACK_GLOBAL(stftri,STFTRI)
#define LAPACK_ztftri LAPACK_GLOBAL(ztftri,ZTFTRI)
#define LAPACK_ctftri LAPACK_GLOBAL(ctftri,CTFTRI)
#define LAPACK_stptri LAPACK_GLOBAL(stptri,STPTRI)
#define LAPACK_dtptri LAPACK_GLOBAL(dtptri,DTPTRI)
#define LAPACK_ctptri LAPACK_GLOBAL(ctptri,CTPTRI)
#define LAPACK_ztptri LAPACK_GLOBAL(ztptri,ZTPTRI)
#define LAPACK_sgeequ LAPACK_GLOBAL(sgeequ,SGEEQU)
#define LAPACK_dgeequ LAPACK_GLOBAL(dgeequ,DGEEQU)
#define LAPACK_cgeequ LAPACK_GLOBAL(cgeequ,CGEEQU)
#define LAPACK_zgeequ LAPACK_GLOBAL(zgeequ,ZGEEQU)
#define LAPACK_dgeequb LAPACK_GLOBAL(dgeequb,DGEEQUB)
#define LAPACK_sgeequb LAPACK_GLOBAL(sgeequb,SGEEQUB)
#define LAPACK_zgeequb LAPACK_GLOBAL(zgeequb,ZGEEQUB)
#define LAPACK_cgeequb LAPACK_GLOBAL(cgeequb,CGEEQUB)
#define LAPACK_sgbequ LAPACK_GLOBAL(sgbequ,SGBEQU)
#define LAPACK_dgbequ LAPACK_GLOBAL(dgbequ,DGBEQU)
#define LAPACK_cgbequ LAPACK_GLOBAL(cgbequ,CGBEQU)
#define LAPACK_zgbequ LAPACK_GLOBAL(zgbequ,ZGBEQU)
#define LAPACK_dgbequb LAPACK_GLOBAL(dgbequb,DGBEQUB)
#define LAPACK_sgbequb LAPACK_GLOBAL(sgbequb,SGBEQUB)
#define LAPACK_zgbequb LAPACK_GLOBAL(zgbequb,ZGBEQUB)
#define LAPACK_cgbequb LAPACK_GLOBAL(cgbequb,CGBEQUB)
#define LAPACK_spoequ LAPACK_GLOBAL(spoequ,SPOEQU)
#define LAPACK_dpoequ LAPACK_GLOBAL(dpoequ,DPOEQU)
#define LAPACK_cpoequ LAPACK_GLOBAL(cpoequ,CPOEQU)
#define LAPACK_zpoequ LAPACK_GLOBAL(zpoequ,ZPOEQU)
#define LAPACK_dpoequb LAPACK_GLOBAL(dpoequb,DPOEQUB)
#define LAPACK_spoequb LAPACK_GLOBAL(spoequb,SPOEQUB)
#define LAPACK_zpoequb LAPACK_GLOBAL(zpoequb,ZPOEQUB)
#define LAPACK_cpoequb LAPACK_GLOBAL(cpoequb,CPOEQUB)
#define LAPACK_sppequ LAPACK_GLOBAL(sppequ,SPPEQU)
#define LAPACK_dppequ LAPACK_GLOBAL(dppequ,DPPEQU)
#define LAPACK_cppequ LAPACK_GLOBAL(cppequ,CPPEQU)
#define LAPACK_zppequ LAPACK_GLOBAL(zppequ,ZPPEQU)
#define LAPACK_spbequ LAPACK_GLOBAL(spbequ,SPBEQU)
#define LAPACK_dpbequ LAPACK_GLOBAL(dpbequ,DPBEQU)
#define LAPACK_cpbequ LAPACK_GLOBAL(cpbequ,CPBEQU)
#define LAPACK_zpbequ LAPACK_GLOBAL(zpbequ,ZPBEQU)
#define LAPACK_dsyequb LAPACK_GLOBAL(dsyequb,DSYEQUB)
#define LAPACK_ssyequb LAPACK_GLOBAL(ssyequb,SSYEQUB)
#define LAPACK_zsyequb LAPACK_GLOBAL(zsyequb,ZSYEQUB)
#define LAPACK_csyequb LAPACK_GLOBAL(csyequb,CSYEQUB)
#define LAPACK_zheequb LAPACK_GLOBAL(zheequb,ZHEEQUB)
#define LAPACK_cheequb LAPACK_GLOBAL(cheequb,CHEEQUB)
#define LAPACK_sgesv LAPACK_GLOBAL(sgesv,SGESV)
#define LAPACK_dgesv LAPACK_GLOBAL(dgesv,DGESV)
#define LAPACK_cgesv LAPACK_GLOBAL(cgesv,CGESV)
#define LAPACK_zgesv LAPACK_GLOBAL(zgesv,ZGESV)
#define LAPACK_dsgesv LAPACK_GLOBAL(dsgesv,DSGESV)
#define LAPACK_zcgesv LAPACK_GLOBAL(zcgesv,ZCGESV)
#define LAPACK_sgesvx LAPACK_GLOBAL(sgesvx,SGESVX)
#define LAPACK_dgesvx LAPACK_GLOBAL(dgesvx,DGESVX)
#define LAPACK_cgesvx LAPACK_GLOBAL(cgesvx,CGESVX)
#define LAPACK_zgesvx LAPACK_GLOBAL(zgesvx,ZGESVX)
#define LAPACK_dgesvxx LAPACK_GLOBAL(dgesvxx,DGESVXX)
#define LAPACK_sgesvxx LAPACK_GLOBAL(sgesvxx,SGESVXX)
#define LAPACK_zgesvxx LAPACK_GLOBAL(zgesvxx,ZGESVXX)
#define LAPACK_cgesvxx LAPACK_GLOBAL(cgesvxx,CGESVXX)
#define LAPACK_sgbsv LAPACK_GLOBAL(sgbsv,SGBSV)
#define LAPACK_dgbsv LAPACK_GLOBAL(dgbsv,DGBSV)
#define LAPACK_cgbsv LAPACK_GLOBAL(cgbsv,CGBSV)
#define LAPACK_zgbsv LAPACK_GLOBAL(zgbsv,ZGBSV)
#define LAPACK_sgbsvx LAPACK_GLOBAL(sgbsvx,SGBSVX)
#define LAPACK_dgbsvx LAPACK_GLOBAL(dgbsvx,DGBSVX)
#define LAPACK_cgbsvx LAPACK_GLOBAL(cgbsvx,CGBSVX)
#define LAPACK_zgbsvx LAPACK_GLOBAL(zgbsvx,ZGBSVX)
#define LAPACK_dgbsvxx LAPACK_GLOBAL(dgbsvxx,DGBSVXX)
#define LAPACK_sgbsvxx LAPACK_GLOBAL(sgbsvxx,SGBSVXX)
#define LAPACK_zgbsvxx LAPACK_GLOBAL(zgbsvxx,ZGBSVXX)
#define LAPACK_cgbsvxx LAPACK_GLOBAL(cgbsvxx,CGBSVXX)
#define LAPACK_sgtsv LAPACK_GLOBAL(sgtsv,SGTSV)
#define LAPACK_dgtsv LAPACK_GLOBAL(dgtsv,DGTSV)
#define LAPACK_cgtsv LAPACK_GLOBAL(cgtsv,CGTSV)
#define LAPACK_zgtsv LAPACK_GLOBAL(zgtsv,ZGTSV)
#define LAPACK_sgtsvx LAPACK_GLOBAL(sgtsvx,SGTSVX)
#define LAPACK_dgtsvx LAPACK_GLOBAL(dgtsvx,DGTSVX)
#define LAPACK_cgtsvx LAPACK_GLOBAL(cgtsvx,CGTSVX)
#define LAPACK_zgtsvx LAPACK_GLOBAL(zgtsvx,ZGTSVX)
#define LAPACK_sposv LAPACK_GLOBAL(sposv,SPOSV)
#define LAPACK_dposv LAPACK_GLOBAL(dposv,DPOSV)
#define LAPACK_cposv LAPACK_GLOBAL(cposv,CPOSV)
#define LAPACK_zposv LAPACK_GLOBAL(zposv,ZPOSV)
#define LAPACK_dsposv LAPACK_GLOBAL(dsposv,DSPOSV)
#define LAPACK_zcposv LAPACK_GLOBAL(zcposv,ZCPOSV)
#define LAPACK_sposvx LAPACK_GLOBAL(sposvx,SPOSVX)
#define LAPACK_dposvx LAPACK_GLOBAL(dposvx,DPOSVX)
#define LAPACK_cposvx LAPACK_GLOBAL(cposvx,CPOSVX)
#define LAPACK_zposvx LAPACK_GLOBAL(zposvx,ZPOSVX)
#define LAPACK_dposvxx LAPACK_GLOBAL(dposvxx,DPOSVXX)
#define LAPACK_sposvxx LAPACK_GLOBAL(sposvxx,SPOSVXX)
#define LAPACK_zposvxx LAPACK_GLOBAL(zposvxx,ZPOSVXX)
#define LAPACK_cposvxx LAPACK_GLOBAL(cposvxx,CPOSVXX)
#define LAPACK_sppsv LAPACK_GLOBAL(sppsv,SPPSV)
#define LAPACK_dppsv LAPACK_GLOBAL(dppsv,DPPSV)
#define LAPACK_cppsv LAPACK_GLOBAL(cppsv,CPPSV)
#define LAPACK_zppsv LAPACK_GLOBAL(zppsv,ZPPSV)
#define LAPACK_sppsvx LAPACK_GLOBAL(sppsvx,SPPSVX)
#define LAPACK_dppsvx LAPACK_GLOBAL(dppsvx,DPPSVX)
#define LAPACK_cppsvx LAPACK_GLOBAL(cppsvx,CPPSVX)
#define LAPACK_zppsvx LAPACK_GLOBAL(zppsvx,ZPPSVX)
#define LAPACK_spbsv LAPACK_GLOBAL(spbsv,SPBSV)
#define LAPACK_dpbsv LAPACK_GLOBAL(dpbsv,DPBSV)
#define LAPACK_cpbsv LAPACK_GLOBAL(cpbsv,CPBSV)
#define LAPACK_zpbsv LAPACK_GLOBAL(zpbsv,ZPBSV)
#define LAPACK_spbsvx LAPACK_GLOBAL(spbsvx,SPBSVX)
#define LAPACK_dpbsvx LAPACK_GLOBAL(dpbsvx,DPBSVX)
#define LAPACK_cpbsvx LAPACK_GLOBAL(cpbsvx,CPBSVX)
#define LAPACK_zpbsvx LAPACK_GLOBAL(zpbsvx,ZPBSVX)
#define LAPACK_sptsv LAPACK_GLOBAL(sptsv,SPTSV)
#define LAPACK_dptsv LAPACK_GLOBAL(dptsv,DPTSV)
#define LAPACK_cptsv LAPACK_GLOBAL(cptsv,CPTSV)
#define LAPACK_zptsv LAPACK_GLOBAL(zptsv,ZPTSV)
#define LAPACK_sptsvx LAPACK_GLOBAL(sptsvx,SPTSVX)
#define LAPACK_dptsvx LAPACK_GLOBAL(dptsvx,DPTSVX)
#define LAPACK_cptsvx LAPACK_GLOBAL(cptsvx,CPTSVX)
#define LAPACK_zptsvx LAPACK_GLOBAL(zptsvx,ZPTSVX)
#define LAPACK_ssysv LAPACK_GLOBAL(ssysv,SSYSV)
#define LAPACK_dsysv LAPACK_GLOBAL(dsysv,DSYSV)
#define LAPACK_csysv LAPACK_GLOBAL(csysv,CSYSV)
#define LAPACK_zsysv LAPACK_GLOBAL(zsysv,ZSYSV)
#define LAPACK_ssysvx LAPACK_GLOBAL(ssysvx,SSYSVX)
#define LAPACK_dsysvx LAPACK_GLOBAL(dsysvx,DSYSVX)
#define LAPACK_csysvx LAPACK_GLOBAL(csysvx,CSYSVX)
#define LAPACK_zsysvx LAPACK_GLOBAL(zsysvx,ZSYSVX)
#define LAPACK_dsysvxx LAPACK_GLOBAL(dsysvxx,DSYSVXX)
#define LAPACK_ssysvxx LAPACK_GLOBAL(ssysvxx,SSYSVXX)
#define LAPACK_zsysvxx LAPACK_GLOBAL(zsysvxx,ZSYSVXX)
#define LAPACK_csysvxx LAPACK_GLOBAL(csysvxx,CSYSVXX)
#define LAPACK_chesv LAPACK_GLOBAL(chesv,CHESV)
#define LAPACK_zhesv LAPACK_GLOBAL(zhesv,ZHESV)
#define LAPACK_chesvx LAPACK_GLOBAL(chesvx,CHESVX)
#define LAPACK_zhesvx LAPACK_GLOBAL(zhesvx,ZHESVX)
#define LAPACK_zhesvxx LAPACK_GLOBAL(zhesvxx,ZHESVXX)
#define LAPACK_chesvxx LAPACK_GLOBAL(chesvxx,CHESVXX)
#define LAPACK_sspsv LAPACK_GLOBAL(sspsv,SSPSV)
#define LAPACK_dspsv LAPACK_GLOBAL(dspsv,DSPSV)
#define LAPACK_cspsv LAPACK_GLOBAL(cspsv,CSPSV)
#define LAPACK_zspsv LAPACK_GLOBAL(zspsv,ZSPSV)
#define LAPACK_sspsvx LAPACK_GLOBAL(sspsvx,SSPSVX)
#define LAPACK_dspsvx LAPACK_GLOBAL(dspsvx,DSPSVX)
#define LAPACK_cspsvx LAPACK_GLOBAL(cspsvx,CSPSVX)
#define LAPACK_zspsvx LAPACK_GLOBAL(zspsvx,ZSPSVX)
#define LAPACK_chpsv LAPACK_GLOBAL(chpsv,CHPSV)
#define LAPACK_zhpsv LAPACK_GLOBAL(zhpsv,ZHPSV)
#define LAPACK_chpsvx LAPACK_GLOBAL(chpsvx,CHPSVX)
#define LAPACK_zhpsvx LAPACK_GLOBAL(zhpsvx,ZHPSVX)
#define LAPACK_sgeqrf LAPACK_GLOBAL(sgeqrf,SGEQRF)
#define LAPACK_dgeqrf LAPACK_GLOBAL(dgeqrf,DGEQRF)
#define LAPACK_cgeqrf LAPACK_GLOBAL(cgeqrf,CGEQRF)
#define LAPACK_zgeqrf LAPACK_GLOBAL(zgeqrf,ZGEQRF)
#define LAPACK_sgeqpf LAPACK_GLOBAL(sgeqpf,SGEQPF)
#define LAPACK_dgeqpf LAPACK_GLOBAL(dgeqpf,DGEQPF)
#define LAPACK_cgeqpf LAPACK_GLOBAL(cgeqpf,CGEQPF)
#define LAPACK_zgeqpf LAPACK_GLOBAL(zgeqpf,ZGEQPF)
#define LAPACK_sgeqp3 LAPACK_GLOBAL(sgeqp3,SGEQP3)
#define LAPACK_dgeqp3 LAPACK_GLOBAL(dgeqp3,DGEQP3)
#define LAPACK_cgeqp3 LAPACK_GLOBAL(cgeqp3,CGEQP3)
#define LAPACK_zgeqp3 LAPACK_GLOBAL(zgeqp3,ZGEQP3)
#define LAPACK_sorgqr LAPACK_GLOBAL(sorgqr,SORGQR)
#define LAPACK_dorgqr LAPACK_GLOBAL(dorgqr,DORGQR)
#define LAPACK_sormqr LAPACK_GLOBAL(sormqr,SORMQR)
#define LAPACK_dormqr LAPACK_GLOBAL(dormqr,DORMQR)
#define LAPACK_cungqr LAPACK_GLOBAL(cungqr,CUNGQR)
#define LAPACK_zungqr LAPACK_GLOBAL(zungqr,ZUNGQR)
#define LAPACK_cunmqr LAPACK_GLOBAL(cunmqr,CUNMQR)
#define LAPACK_zunmqr LAPACK_GLOBAL(zunmqr,ZUNMQR)
#define LAPACK_sgelqf LAPACK_GLOBAL(sgelqf,SGELQF)
#define LAPACK_dgelqf LAPACK_GLOBAL(dgelqf,DGELQF)
#define LAPACK_cgelqf LAPACK_GLOBAL(cgelqf,CGELQF)
#define LAPACK_zgelqf LAPACK_GLOBAL(zgelqf,ZGELQF)
#define LAPACK_sorglq LAPACK_GLOBAL(sorglq,SORGLQ)
#define LAPACK_dorglq LAPACK_GLOBAL(dorglq,DORGLQ)
#define LAPACK_sormlq LAPACK_GLOBAL(sormlq,SORMLQ)
#define LAPACK_dormlq LAPACK_GLOBAL(dormlq,DORMLQ)
#define LAPACK_cunglq LAPACK_GLOBAL(cunglq,CUNGLQ)
#define LAPACK_zunglq LAPACK_GLOBAL(zunglq,ZUNGLQ)
#define LAPACK_cunmlq LAPACK_GLOBAL(cunmlq,CUNMLQ)
#define LAPACK_zunmlq LAPACK_GLOBAL(zunmlq,ZUNMLQ)
#define LAPACK_sgeqlf LAPACK_GLOBAL(sgeqlf,SGEQLF)
#define LAPACK_dgeqlf LAPACK_GLOBAL(dgeqlf,DGEQLF)
#define LAPACK_cgeqlf LAPACK_GLOBAL(cgeqlf,CGEQLF)
#define LAPACK_zgeqlf LAPACK_GLOBAL(zgeqlf,ZGEQLF)
#define LAPACK_sorgql LAPACK_GLOBAL(sorgql,SORGQL)
#define LAPACK_dorgql LAPACK_GLOBAL(dorgql,DORGQL)
#define LAPACK_cungql LAPACK_GLOBAL(cungql,CUNGQL)
#define LAPACK_zungql LAPACK_GLOBAL(zungql,ZUNGQL)
#define LAPACK_sormql LAPACK_GLOBAL(sormql,SORMQL)
#define LAPACK_dormql LAPACK_GLOBAL(dormql,DORMQL)
#define LAPACK_cunmql LAPACK_GLOBAL(cunmql,CUNMQL)
#define LAPACK_zunmql LAPACK_GLOBAL(zunmql,ZUNMQL)
#define LAPACK_sgerqf LAPACK_GLOBAL(sgerqf,SGERQF)
#define LAPACK_dgerqf LAPACK_GLOBAL(dgerqf,DGERQF)
#define LAPACK_cgerqf LAPACK_GLOBAL(cgerqf,CGERQF)
#define LAPACK_zgerqf LAPACK_GLOBAL(zgerqf,ZGERQF)
#define LAPACK_sorgrq LAPACK_GLOBAL(sorgrq,SORGRQ)
#define LAPACK_dorgrq LAPACK_GLOBAL(dorgrq,DORGRQ)
#define LAPACK_cungrq LAPACK_GLOBAL(cungrq,CUNGRQ)
#define LAPACK_zungrq LAPACK_GLOBAL(zungrq,ZUNGRQ)
#define LAPACK_sormrq LAPACK_GLOBAL(sormrq,SORMRQ)
#define LAPACK_dormrq LAPACK_GLOBAL(dormrq,DORMRQ)
#define LAPACK_cunmrq LAPACK_GLOBAL(cunmrq,CUNMRQ)
#define LAPACK_zunmrq LAPACK_GLOBAL(zunmrq,ZUNMRQ)
#define LAPACK_stzrzf LAPACK_GLOBAL(stzrzf,STZRZF)
#define LAPACK_dtzrzf LAPACK_GLOBAL(dtzrzf,DTZRZF)
#define LAPACK_ctzrzf LAPACK_GLOBAL(ctzrzf,CTZRZF)
#define LAPACK_ztzrzf LAPACK_GLOBAL(ztzrzf,ZTZRZF)
#define LAPACK_sormrz LAPACK_GLOBAL(sormrz,SORMRZ)
#define LAPACK_dormrz LAPACK_GLOBAL(dormrz,DORMRZ)
#define LAPACK_cunmrz LAPACK_GLOBAL(cunmrz,CUNMRZ)
#define LAPACK_zunmrz LAPACK_GLOBAL(zunmrz,ZUNMRZ)
#define LAPACK_sggqrf LAPACK_GLOBAL(sggqrf,SGGQRF)
#define LAPACK_dggqrf LAPACK_GLOBAL(dggqrf,DGGQRF)
#define LAPACK_cggqrf LAPACK_GLOBAL(cggqrf,CGGQRF)
#define LAPACK_zggqrf LAPACK_GLOBAL(zggqrf,ZGGQRF)
#define LAPACK_sggrqf LAPACK_GLOBAL(sggrqf,SGGRQF)
#define LAPACK_dggrqf LAPACK_GLOBAL(dggrqf,DGGRQF)
#define LAPACK_cggrqf LAPACK_GLOBAL(cggrqf,CGGRQF)
#define LAPACK_zggrqf LAPACK_GLOBAL(zggrqf,ZGGRQF)
#define LAPACK_sgebrd LAPACK_GLOBAL(sgebrd,SGEBRD)
#define LAPACK_dgebrd LAPACK_GLOBAL(dgebrd,DGEBRD)
#define LAPACK_cgebrd LAPACK_GLOBAL(cgebrd,CGEBRD)
#define LAPACK_zgebrd LAPACK_GLOBAL(zgebrd,ZGEBRD)
#define LAPACK_sgbbrd LAPACK_GLOBAL(sgbbrd,SGBBRD)
#define LAPACK_dgbbrd LAPACK_GLOBAL(dgbbrd,DGBBRD)
#define LAPACK_cgbbrd LAPACK_GLOBAL(cgbbrd,CGBBRD)
#define LAPACK_zgbbrd LAPACK_GLOBAL(zgbbrd,ZGBBRD)
#define LAPACK_sorgbr LAPACK_GLOBAL(sorgbr,SORGBR)
#define LAPACK_dorgbr LAPACK_GLOBAL(dorgbr,DORGBR)
#define LAPACK_sormbr LAPACK_GLOBAL(sormbr,SORMBR)
#define LAPACK_dormbr LAPACK_GLOBAL(dormbr,DORMBR)
#define LAPACK_cungbr LAPACK_GLOBAL(cungbr,CUNGBR)
#define LAPACK_zungbr LAPACK_GLOBAL(zungbr,ZUNGBR)
#define LAPACK_cunmbr LAPACK_GLOBAL(cunmbr,CUNMBR)
#define LAPACK_zunmbr LAPACK_GLOBAL(zunmbr,ZUNMBR)
#define LAPACK_sbdsqr LAPACK_GLOBAL(sbdsqr,SBDSQR)
#define LAPACK_dbdsqr LAPACK_GLOBAL(dbdsqr,DBDSQR)
#define LAPACK_cbdsqr LAPACK_GLOBAL(cbdsqr,CBDSQR)
#define LAPACK_zbdsqr LAPACK_GLOBAL(zbdsqr,ZBDSQR)
#define LAPACK_sbdsdc LAPACK_GLOBAL(sbdsdc,SBDSDC)
#define LAPACK_dbdsdc LAPACK_GLOBAL(dbdsdc,DBDSDC)
#define LAPACK_ssytrd LAPACK_GLOBAL(ssytrd,SSYTRD)
#define LAPACK_dsytrd LAPACK_GLOBAL(dsytrd,DSYTRD)
#define LAPACK_sorgtr LAPACK_GLOBAL(sorgtr,SORGTR)
#define LAPACK_dorgtr LAPACK_GLOBAL(dorgtr,DORGTR)
#define LAPACK_sormtr LAPACK_GLOBAL(sormtr,SORMTR)
#define LAPACK_dormtr LAPACK_GLOBAL(dormtr,DORMTR)
#define LAPACK_chetrd LAPACK_GLOBAL(chetrd,CHETRD)
#define LAPACK_zhetrd LAPACK_GLOBAL(zhetrd,ZHETRD)
#define LAPACK_cungtr LAPACK_GLOBAL(cungtr,CUNGTR)
#define LAPACK_zungtr LAPACK_GLOBAL(zungtr,ZUNGTR)
#define LAPACK_cunmtr LAPACK_GLOBAL(cunmtr,CUNMTR)
#define LAPACK_zunmtr LAPACK_GLOBAL(zunmtr,ZUNMTR)
#define LAPACK_ssptrd LAPACK_GLOBAL(ssptrd,SSPTRD)
#define LAPACK_dsptrd LAPACK_GLOBAL(dsptrd,DSPTRD)
#define LAPACK_sopgtr LAPACK_GLOBAL(sopgtr,SOPGTR)
#define LAPACK_dopgtr LAPACK_GLOBAL(dopgtr,DOPGTR)
#define LAPACK_sopmtr LAPACK_GLOBAL(sopmtr,SOPMTR)
#define LAPACK_dopmtr LAPACK_GLOBAL(dopmtr,DOPMTR)
#define LAPACK_chptrd LAPACK_GLOBAL(chptrd,CHPTRD)
#define LAPACK_zhptrd LAPACK_GLOBAL(zhptrd,ZHPTRD)
#define LAPACK_cupgtr LAPACK_GLOBAL(cupgtr,CUPGTR)
#define LAPACK_zupgtr LAPACK_GLOBAL(zupgtr,ZUPGTR)
#define LAPACK_cupmtr LAPACK_GLOBAL(cupmtr,CUPMTR)
#define LAPACK_zupmtr LAPACK_GLOBAL(zupmtr,ZUPMTR)
#define LAPACK_ssbtrd LAPACK_GLOBAL(ssbtrd,SSBTRD)
#define LAPACK_dsbtrd LAPACK_GLOBAL(dsbtrd,DSBTRD)
#define LAPACK_chbtrd LAPACK_GLOBAL(chbtrd,CHBTRD)
#define LAPACK_zhbtrd LAPACK_GLOBAL(zhbtrd,ZHBTRD)
#define LAPACK_ssterf LAPACK_GLOBAL(ssterf,SSTERF)
#define LAPACK_dsterf LAPACK_GLOBAL(dsterf,DSTERF)
#define LAPACK_ssteqr LAPACK_GLOBAL(ssteqr,SSTEQR)
#define LAPACK_dsteqr LAPACK_GLOBAL(dsteqr,DSTEQR)
#define LAPACK_csteqr LAPACK_GLOBAL(csteqr,CSTEQR)
#define LAPACK_zsteqr LAPACK_GLOBAL(zsteqr,ZSTEQR)
#define LAPACK_sstemr LAPACK_GLOBAL(sstemr,SSTEMR)
#define LAPACK_dstemr LAPACK_GLOBAL(dstemr,DSTEMR)
#define LAPACK_cstemr LAPACK_GLOBAL(cstemr,CSTEMR)
#define LAPACK_zstemr LAPACK_GLOBAL(zstemr,ZSTEMR)
#define LAPACK_sstedc LAPACK_GLOBAL(sstedc,SSTEDC)
#define LAPACK_dstedc LAPACK_GLOBAL(dstedc,DSTEDC)
#define LAPACK_cstedc LAPACK_GLOBAL(cstedc,CSTEDC)
#define LAPACK_zstedc LAPACK_GLOBAL(zstedc,ZSTEDC)
#define LAPACK_sstegr LAPACK_GLOBAL(sstegr,SSTEGR)
#define LAPACK_dstegr LAPACK_GLOBAL(dstegr,DSTEGR)
#define LAPACK_cstegr LAPACK_GLOBAL(cstegr,CSTEGR)
#define LAPACK_zstegr LAPACK_GLOBAL(zstegr,ZSTEGR)
#define LAPACK_spteqr LAPACK_GLOBAL(spteqr,SPTEQR)
#define LAPACK_dpteqr LAPACK_GLOBAL(dpteqr,DPTEQR)
#define LAPACK_cpteqr LAPACK_GLOBAL(cpteqr,CPTEQR)
#define LAPACK_zpteqr LAPACK_GLOBAL(zpteqr,ZPTEQR)
#define LAPACK_sstebz LAPACK_GLOBAL(sstebz,SSTEBZ)
#define LAPACK_dstebz LAPACK_GLOBAL(dstebz,DSTEBZ)
#define LAPACK_sstein LAPACK_GLOBAL(sstein,SSTEIN)
#define LAPACK_dstein LAPACK_GLOBAL(dstein,DSTEIN)
#define LAPACK_cstein LAPACK_GLOBAL(cstein,CSTEIN)
#define LAPACK_zstein LAPACK_GLOBAL(zstein,ZSTEIN)
#define LAPACK_sdisna LAPACK_GLOBAL(sdisna,SDISNA)
#define LAPACK_ddisna LAPACK_GLOBAL(ddisna,DDISNA)
#define LAPACK_ssygst LAPACK_GLOBAL(ssygst,SSYGST)
#define LAPACK_dsygst LAPACK_GLOBAL(dsygst,DSYGST)
#define LAPACK_chegst LAPACK_GLOBAL(chegst,CHEGST)
#define LAPACK_zhegst LAPACK_GLOBAL(zhegst,ZHEGST)
#define LAPACK_sspgst LAPACK_GLOBAL(sspgst,SSPGST)
#define LAPACK_dspgst LAPACK_GLOBAL(dspgst,DSPGST)
#define LAPACK_chpgst LAPACK_GLOBAL(chpgst,CHPGST)
#define LAPACK_zhpgst LAPACK_GLOBAL(zhpgst,ZHPGST)
#define LAPACK_ssbgst LAPACK_GLOBAL(ssbgst,SSBGST)
#define LAPACK_dsbgst LAPACK_GLOBAL(dsbgst,DSBGST)
#define LAPACK_chbgst LAPACK_GLOBAL(chbgst,CHBGST)
#define LAPACK_zhbgst LAPACK_GLOBAL(zhbgst,ZHBGST)
#define LAPACK_spbstf LAPACK_GLOBAL(spbstf,SPBSTF)
#define LAPACK_dpbstf LAPACK_GLOBAL(dpbstf,DPBSTF)
#define LAPACK_cpbstf LAPACK_GLOBAL(cpbstf,CPBSTF)
#define LAPACK_zpbstf LAPACK_GLOBAL(zpbstf,ZPBSTF)
#define LAPACK_sgehrd LAPACK_GLOBAL(sgehrd,SGEHRD)
#define LAPACK_dgehrd LAPACK_GLOBAL(dgehrd,DGEHRD)
#define LAPACK_cgehrd LAPACK_GLOBAL(cgehrd,CGEHRD)
#define LAPACK_zgehrd LAPACK_GLOBAL(zgehrd,ZGEHRD)
#define LAPACK_sorghr LAPACK_GLOBAL(sorghr,SORGHR)
#define LAPACK_dorghr LAPACK_GLOBAL(dorghr,DORGHR)
#define LAPACK_sormhr LAPACK_GLOBAL(sormhr,SORMHR)
#define LAPACK_dormhr LAPACK_GLOBAL(dormhr,DORMHR)
#define LAPACK_cunghr LAPACK_GLOBAL(cunghr,CUNGHR)
#define LAPACK_zunghr LAPACK_GLOBAL(zunghr,ZUNGHR)
#define LAPACK_cunmhr LAPACK_GLOBAL(cunmhr,CUNMHR)
#define LAPACK_zunmhr LAPACK_GLOBAL(zunmhr,ZUNMHR)
#define LAPACK_sgebal LAPACK_GLOBAL(sgebal,SGEBAL)
#define LAPACK_dgebal LAPACK_GLOBAL(dgebal,DGEBAL)
#define LAPACK_cgebal LAPACK_GLOBAL(cgebal,CGEBAL)
#define LAPACK_zgebal LAPACK_GLOBAL(zgebal,ZGEBAL)
#define LAPACK_sgebak LAPACK_GLOBAL(sgebak,SGEBAK)
#define LAPACK_dgebak LAPACK_GLOBAL(dgebak,DGEBAK)
#define LAPACK_cgebak LAPACK_GLOBAL(cgebak,CGEBAK)
#define LAPACK_zgebak LAPACK_GLOBAL(zgebak,ZGEBAK)
#define LAPACK_shseqr LAPACK_GLOBAL(shseqr,SHSEQR)
#define LAPACK_dhseqr LAPACK_GLOBAL(dhseqr,DHSEQR)
#define LAPACK_chseqr LAPACK_GLOBAL(chseqr,CHSEQR)
#define LAPACK_zhseqr LAPACK_GLOBAL(zhseqr,ZHSEQR)
#define LAPACK_shsein LAPACK_GLOBAL(shsein,SHSEIN)
#define LAPACK_dhsein LAPACK_GLOBAL(dhsein,DHSEIN)
#define LAPACK_chsein LAPACK_GLOBAL(chsein,CHSEIN)
#define LAPACK_zhsein LAPACK_GLOBAL(zhsein,ZHSEIN)
#define LAPACK_strevc LAPACK_GLOBAL(strevc,STREVC)
#define LAPACK_dtrevc LAPACK_GLOBAL(dtrevc,DTREVC)
#define LAPACK_ctrevc LAPACK_GLOBAL(ctrevc,CTREVC)
#define LAPACK_ztrevc LAPACK_GLOBAL(ztrevc,ZTREVC)
#define LAPACK_strsna LAPACK_GLOBAL(strsna,STRSNA)
#define LAPACK_dtrsna LAPACK_GLOBAL(dtrsna,DTRSNA)
#define LAPACK_ctrsna LAPACK_GLOBAL(ctrsna,CTRSNA)
#define LAPACK_ztrsna LAPACK_GLOBAL(ztrsna,ZTRSNA)
#define LAPACK_strexc LAPACK_GLOBAL(strexc,STREXC)
#define LAPACK_dtrexc LAPACK_GLOBAL(dtrexc,DTREXC)
#define LAPACK_ctrexc LAPACK_GLOBAL(ctrexc,CTREXC)
#define LAPACK_ztrexc LAPACK_GLOBAL(ztrexc,ZTREXC)
#define LAPACK_strsen LAPACK_GLOBAL(strsen,STRSEN)
#define LAPACK_dtrsen LAPACK_GLOBAL(dtrsen,DTRSEN)
#define LAPACK_ctrsen LAPACK_GLOBAL(ctrsen,CTRSEN)
#define LAPACK_ztrsen LAPACK_GLOBAL(ztrsen,ZTRSEN)
#define LAPACK_strsyl LAPACK_GLOBAL(strsyl,STRSYL)
#define LAPACK_dtrsyl LAPACK_GLOBAL(dtrsyl,DTRSYL)
#define LAPACK_ctrsyl LAPACK_GLOBAL(ctrsyl,CTRSYL)
#define LAPACK_ztrsyl LAPACK_GLOBAL(ztrsyl,ZTRSYL)
#define LAPACK_sgghrd LAPACK_GLOBAL(sgghrd,SGGHRD)
#define LAPACK_dgghrd LAPACK_GLOBAL(dgghrd,DGGHRD)
#define LAPACK_cgghrd LAPACK_GLOBAL(cgghrd,CGGHRD)
#define LAPACK_zgghrd LAPACK_GLOBAL(zgghrd,ZGGHRD)
#define LAPACK_sggbal LAPACK_GLOBAL(sggbal,SGGBAL)
#define LAPACK_dggbal LAPACK_GLOBAL(dggbal,DGGBAL)
#define LAPACK_cggbal LAPACK_GLOBAL(cggbal,CGGBAL)
#define LAPACK_zggbal LAPACK_GLOBAL(zggbal,ZGGBAL)
#define LAPACK_sggbak LAPACK_GLOBAL(sggbak,SGGBAK)
#define LAPACK_dggbak LAPACK_GLOBAL(dggbak,DGGBAK)
#define LAPACK_cggbak LAPACK_GLOBAL(cggbak,CGGBAK)
#define LAPACK_zggbak LAPACK_GLOBAL(zggbak,ZGGBAK)
#define LAPACK_shgeqz LAPACK_GLOBAL(shgeqz,SHGEQZ)
#define LAPACK_dhgeqz LAPACK_GLOBAL(dhgeqz,DHGEQZ)
#define LAPACK_chgeqz LAPACK_GLOBAL(chgeqz,CHGEQZ)
#define LAPACK_zhgeqz LAPACK_GLOBAL(zhgeqz,ZHGEQZ)
#define LAPACK_stgevc LAPACK_GLOBAL(stgevc,STGEVC)
#define LAPACK_dtgevc LAPACK_GLOBAL(dtgevc,DTGEVC)
#define LAPACK_ctgevc LAPACK_GLOBAL(ctgevc,CTGEVC)
#define LAPACK_ztgevc LAPACK_GLOBAL(ztgevc,ZTGEVC)
#define LAPACK_stgexc LAPACK_GLOBAL(stgexc,STGEXC)
#define LAPACK_dtgexc LAPACK_GLOBAL(dtgexc,DTGEXC)
#define LAPACK_ctgexc LAPACK_GLOBAL(ctgexc,CTGEXC)
#define LAPACK_ztgexc LAPACK_GLOBAL(ztgexc,ZTGEXC)
#define LAPACK_stgsen LAPACK_GLOBAL(stgsen,STGSEN)
#define LAPACK_dtgsen LAPACK_GLOBAL(dtgsen,DTGSEN)
#define LAPACK_ctgsen LAPACK_GLOBAL(ctgsen,CTGSEN)
#define LAPACK_ztgsen LAPACK_GLOBAL(ztgsen,ZTGSEN)
#define LAPACK_stgsyl LAPACK_GLOBAL(stgsyl,STGSYL)
#define LAPACK_dtgsyl LAPACK_GLOBAL(dtgsyl,DTGSYL)
#define LAPACK_ctgsyl LAPACK_GLOBAL(ctgsyl,CTGSYL)
#define LAPACK_ztgsyl LAPACK_GLOBAL(ztgsyl,ZTGSYL)
#define LAPACK_stgsna LAPACK_GLOBAL(stgsna,STGSNA)
#define LAPACK_dtgsna LAPACK_GLOBAL(dtgsna,DTGSNA)
#define LAPACK_ctgsna LAPACK_GLOBAL(ctgsna,CTGSNA)
#define LAPACK_ztgsna LAPACK_GLOBAL(ztgsna,ZTGSNA)
#define LAPACK_sggsvp LAPACK_GLOBAL(sggsvp,SGGSVP)
#define LAPACK_dggsvp LAPACK_GLOBAL(dggsvp,DGGSVP)
#define LAPACK_cggsvp LAPACK_GLOBAL(cggsvp,CGGSVP)
#define LAPACK_zggsvp LAPACK_GLOBAL(zggsvp,ZGGSVP)
#define LAPACK_stgsja LAPACK_GLOBAL(stgsja,STGSJA)
#define LAPACK_dtgsja LAPACK_GLOBAL(dtgsja,DTGSJA)
#define LAPACK_ctgsja LAPACK_GLOBAL(ctgsja,CTGSJA)
#define LAPACK_ztgsja LAPACK_GLOBAL(ztgsja,ZTGSJA)
#define LAPACK_sgels LAPACK_GLOBAL(sgels,SGELS)
#define LAPACK_dgels LAPACK_GLOBAL(dgels,DGELS)
#define LAPACK_cgels LAPACK_GLOBAL(cgels,CGELS)
#define LAPACK_zgels LAPACK_GLOBAL(zgels,ZGELS)
#define LAPACK_sgelsy LAPACK_GLOBAL(sgelsy,SGELSY)
#define LAPACK_dgelsy LAPACK_GLOBAL(dgelsy,DGELSY)
#define LAPACK_cgelsy LAPACK_GLOBAL(cgelsy,CGELSY)
#define LAPACK_zgelsy LAPACK_GLOBAL(zgelsy,ZGELSY)
#define LAPACK_sgelss LAPACK_GLOBAL(sgelss,SGELSS)
#define LAPACK_dgelss LAPACK_GLOBAL(dgelss,DGELSS)
#define LAPACK_cgelss LAPACK_GLOBAL(cgelss,CGELSS)
#define LAPACK_zgelss LAPACK_GLOBAL(zgelss,ZGELSS)
#define LAPACK_sgelsd LAPACK_GLOBAL(sgelsd,SGELSD)
#define LAPACK_dgelsd LAPACK_GLOBAL(dgelsd,DGELSD)
#define LAPACK_cgelsd LAPACK_GLOBAL(cgelsd,CGELSD)
#define LAPACK_zgelsd LAPACK_GLOBAL(zgelsd,ZGELSD)
#define LAPACK_sgglse LAPACK_GLOBAL(sgglse,SGGLSE)
#define LAPACK_dgglse LAPACK_GLOBAL(dgglse,DGGLSE)
#define LAPACK_cgglse LAPACK_GLOBAL(cgglse,CGGLSE)
#define LAPACK_zgglse LAPACK_GLOBAL(zgglse,ZGGLSE)
#define LAPACK_sggglm LAPACK_GLOBAL(sggglm,SGGGLM)
#define LAPACK_dggglm LAPACK_GLOBAL(dggglm,DGGGLM)
#define LAPACK_cggglm LAPACK_GLOBAL(cggglm,CGGGLM)
#define LAPACK_zggglm LAPACK_GLOBAL(zggglm,ZGGGLM)
#define LAPACK_ssyev LAPACK_GLOBAL(ssyev,SSYEV)
#define LAPACK_dsyev LAPACK_GLOBAL(dsyev,DSYEV)
#define LAPACK_cheev LAPACK_GLOBAL(cheev,CHEEV)
#define LAPACK_zheev LAPACK_GLOBAL(zheev,ZHEEV)
#define LAPACK_ssyevd LAPACK_GLOBAL(ssyevd,SSYEVD)
#define LAPACK_dsyevd LAPACK_GLOBAL(dsyevd,DSYEVD)
#define LAPACK_cheevd LAPACK_GLOBAL(cheevd,CHEEVD)
#define LAPACK_zheevd LAPACK_GLOBAL(zheevd,ZHEEVD)
#define LAPACK_ssyevx LAPACK_GLOBAL(ssyevx,SSYEVX)
#define LAPACK_dsyevx LAPACK_GLOBAL(dsyevx,DSYEVX)
#define LAPACK_cheevx LAPACK_GLOBAL(cheevx,CHEEVX)
#define LAPACK_zheevx LAPACK_GLOBAL(zheevx,ZHEEVX)
#define LAPACK_ssyevr LAPACK_GLOBAL(ssyevr,SSYEVR)
#define LAPACK_dsyevr LAPACK_GLOBAL(dsyevr,DSYEVR)
#define LAPACK_cheevr LAPACK_GLOBAL(cheevr,CHEEVR)
#define LAPACK_zheevr LAPACK_GLOBAL(zheevr,ZHEEVR)
#define LAPACK_sspev LAPACK_GLOBAL(sspev,SSPEV)
#define LAPACK_dspev LAPACK_GLOBAL(dspev,DSPEV)
#define LAPACK_chpev LAPACK_GLOBAL(chpev,CHPEV)
#define LAPACK_zhpev LAPACK_GLOBAL(zhpev,ZHPEV)
#define LAPACK_sspevd LAPACK_GLOBAL(sspevd,SSPEVD)
#define LAPACK_dspevd LAPACK_GLOBAL(dspevd,DSPEVD)
#define LAPACK_chpevd LAPACK_GLOBAL(chpevd,CHPEVD)
#define LAPACK_zhpevd LAPACK_GLOBAL(zhpevd,ZHPEVD)
#define LAPACK_sspevx LAPACK_GLOBAL(sspevx,SSPEVX)
#define LAPACK_dspevx LAPACK_GLOBAL(dspevx,DSPEVX)
#define LAPACK_chpevx LAPACK_GLOBAL(chpevx,CHPEVX)
#define LAPACK_zhpevx LAPACK_GLOBAL(zhpevx,ZHPEVX)
#define LAPACK_ssbev LAPACK_GLOBAL(ssbev,SSBEV)
#define LAPACK_dsbev LAPACK_GLOBAL(dsbev,DSBEV)
#define LAPACK_chbev LAPACK_GLOBAL(chbev,CHBEV)
#define LAPACK_zhbev LAPACK_GLOBAL(zhbev,ZHBEV)
#define LAPACK_ssbevd LAPACK_GLOBAL(ssbevd,SSBEVD)
#define LAPACK_dsbevd LAPACK_GLOBAL(dsbevd,DSBEVD)
#define LAPACK_chbevd LAPACK_GLOBAL(chbevd,CHBEVD)
#define LAPACK_zhbevd LAPACK_GLOBAL(zhbevd,ZHBEVD)
#define LAPACK_ssbevx LAPACK_GLOBAL(ssbevx,SSBEVX)
#define LAPACK_dsbevx LAPACK_GLOBAL(dsbevx,DSBEVX)
#define LAPACK_chbevx LAPACK_GLOBAL(chbevx,CHBEVX)
#define LAPACK_zhbevx LAPACK_GLOBAL(zhbevx,ZHBEVX)
#define LAPACK_sstev LAPACK_GLOBAL(sstev,SSTEV)
#define LAPACK_dstev LAPACK_GLOBAL(dstev,DSTEV)
#define LAPACK_sstevd LAPACK_GLOBAL(sstevd,SSTEVD)
#define LAPACK_dstevd LAPACK_GLOBAL(dstevd,DSTEVD)
#define LAPACK_sstevx LAPACK_GLOBAL(sstevx,SSTEVX)
#define LAPACK_dstevx LAPACK_GLOBAL(dstevx,DSTEVX)
#define LAPACK_sstevr LAPACK_GLOBAL(sstevr,SSTEVR)
#define LAPACK_dstevr LAPACK_GLOBAL(dstevr,DSTEVR)
#define LAPACK_sgees LAPACK_GLOBAL(sgees,SGEES)
#define LAPACK_dgees LAPACK_GLOBAL(dgees,DGEES)
#define LAPACK_cgees LAPACK_GLOBAL(cgees,CGEES)
#define LAPACK_zgees LAPACK_GLOBAL(zgees,ZGEES)
#define LAPACK_sgeesx LAPACK_GLOBAL(sgeesx,SGEESX)
#define LAPACK_dgeesx LAPACK_GLOBAL(dgeesx,DGEESX)
#define LAPACK_cgeesx LAPACK_GLOBAL(cgeesx,CGEESX)
#define LAPACK_zgeesx LAPACK_GLOBAL(zgeesx,ZGEESX)
#define LAPACK_sgeev LAPACK_GLOBAL(sgeev,SGEEV)
#define LAPACK_dgeev LAPACK_GLOBAL(dgeev,DGEEV)
#define LAPACK_cgeev LAPACK_GLOBAL(cgeev,CGEEV)
#define LAPACK_zgeev LAPACK_GLOBAL(zgeev,ZGEEV)
#define LAPACK_sgeevx LAPACK_GLOBAL(sgeevx,SGEEVX)
#define LAPACK_dgeevx LAPACK_GLOBAL(dgeevx,DGEEVX)
#define LAPACK_cgeevx LAPACK_GLOBAL(cgeevx,CGEEVX)
#define LAPACK_zgeevx LAPACK_GLOBAL(zgeevx,ZGEEVX)
#define LAPACK_sgesvd LAPACK_GLOBAL(sgesvd,SGESVD)
#define LAPACK_dgesvd LAPACK_GLOBAL(dgesvd,DGESVD)
#define LAPACK_cgesvd LAPACK_GLOBAL(cgesvd,CGESVD)
#define LAPACK_zgesvd LAPACK_GLOBAL(zgesvd,ZGESVD)
#define LAPACK_sgesdd LAPACK_GLOBAL(sgesdd,SGESDD)
#define LAPACK_dgesdd LAPACK_GLOBAL(dgesdd,DGESDD)
#define LAPACK_cgesdd LAPACK_GLOBAL(cgesdd,CGESDD)
#define LAPACK_zgesdd LAPACK_GLOBAL(zgesdd,ZGESDD)
#define LAPACK_dgejsv LAPACK_GLOBAL(dgejsv,DGEJSV)
#define LAPACK_sgejsv LAPACK_GLOBAL(sgejsv,SGEJSV)
#define LAPACK_dgesvj LAPACK_GLOBAL(dgesvj,DGESVJ)
#define LAPACK_sgesvj LAPACK_GLOBAL(sgesvj,SGESVJ)
#define LAPACK_sggsvd LAPACK_GLOBAL(sggsvd,SGGSVD)
#define LAPACK_dggsvd LAPACK_GLOBAL(dggsvd,DGGSVD)
#define LAPACK_cggsvd LAPACK_GLOBAL(cggsvd,CGGSVD)
#define LAPACK_zggsvd LAPACK_GLOBAL(zggsvd,ZGGSVD)
#define LAPACK_ssygv LAPACK_GLOBAL(ssygv,SSYGV)
#define LAPACK_dsygv LAPACK_GLOBAL(dsygv,DSYGV)
#define LAPACK_chegv LAPACK_GLOBAL(chegv,CHEGV)
#define LAPACK_zhegv LAPACK_GLOBAL(zhegv,ZHEGV)
#define LAPACK_ssygvd LAPACK_GLOBAL(ssygvd,SSYGVD)
#define LAPACK_dsygvd LAPACK_GLOBAL(dsygvd,DSYGVD)
#define LAPACK_chegvd LAPACK_GLOBAL(chegvd,CHEGVD)
#define LAPACK_zhegvd LAPACK_GLOBAL(zhegvd,ZHEGVD)
#define LAPACK_ssygvx LAPACK_GLOBAL(ssygvx,SSYGVX)
#define LAPACK_dsygvx LAPACK_GLOBAL(dsygvx,DSYGVX)
#define LAPACK_chegvx LAPACK_GLOBAL(chegvx,CHEGVX)
#define LAPACK_zhegvx LAPACK_GLOBAL(zhegvx,ZHEGVX)
#define LAPACK_sspgv LAPACK_GLOBAL(sspgv,SSPGV)
#define LAPACK_dspgv LAPACK_GLOBAL(dspgv,DSPGV)
#define LAPACK_chpgv LAPACK_GLOBAL(chpgv,CHPGV)
#define LAPACK_zhpgv LAPACK_GLOBAL(zhpgv,ZHPGV)
#define LAPACK_sspgvd LAPACK_GLOBAL(sspgvd,SSPGVD)
#define LAPACK_dspgvd LAPACK_GLOBAL(dspgvd,DSPGVD)
#define LAPACK_chpgvd LAPACK_GLOBAL(chpgvd,CHPGVD)
#define LAPACK_zhpgvd LAPACK_GLOBAL(zhpgvd,ZHPGVD)
#define LAPACK_sspgvx LAPACK_GLOBAL(sspgvx,SSPGVX)
#define LAPACK_dspgvx LAPACK_GLOBAL(dspgvx,DSPGVX)
#define LAPACK_chpgvx LAPACK_GLOBAL(chpgvx,CHPGVX)
#define LAPACK_zhpgvx LAPACK_GLOBAL(zhpgvx,ZHPGVX)
#define LAPACK_ssbgv LAPACK_GLOBAL(ssbgv,SSBGV)
#define LAPACK_dsbgv LAPACK_GLOBAL(dsbgv,DSBGV)
#define LAPACK_chbgv LAPACK_GLOBAL(chbgv,CHBGV)
#define LAPACK_zhbgv LAPACK_GLOBAL(zhbgv,ZHBGV)
#define LAPACK_ssbgvd LAPACK_GLOBAL(ssbgvd,SSBGVD)
#define LAPACK_dsbgvd LAPACK_GLOBAL(dsbgvd,DSBGVD)
#define LAPACK_chbgvd LAPACK_GLOBAL(chbgvd,CHBGVD)
#define LAPACK_zhbgvd LAPACK_GLOBAL(zhbgvd,ZHBGVD)
#define LAPACK_ssbgvx LAPACK_GLOBAL(ssbgvx,SSBGVX)
#define LAPACK_dsbgvx LAPACK_GLOBAL(dsbgvx,DSBGVX)
#define LAPACK_chbgvx LAPACK_GLOBAL(chbgvx,CHBGVX)
#define LAPACK_zhbgvx LAPACK_GLOBAL(zhbgvx,ZHBGVX)
#define LAPACK_sgges LAPACK_GLOBAL(sgges,SGGES)
#define LAPACK_dgges LAPACK_GLOBAL(dgges,DGGES)
#define LAPACK_cgges LAPACK_GLOBAL(cgges,CGGES)
#define LAPACK_zgges LAPACK_GLOBAL(zgges,ZGGES)
#define LAPACK_sggesx LAPACK_GLOBAL(sggesx,SGGESX)
#define LAPACK_dggesx LAPACK_GLOBAL(dggesx,DGGESX)
#define LAPACK_cggesx LAPACK_GLOBAL(cggesx,CGGESX)
#define LAPACK_zggesx LAPACK_GLOBAL(zggesx,ZGGESX)
#define LAPACK_sggev LAPACK_GLOBAL(sggev,SGGEV)
#define LAPACK_dggev LAPACK_GLOBAL(dggev,DGGEV)
#define LAPACK_cggev LAPACK_GLOBAL(cggev,CGGEV)
#define LAPACK_zggev LAPACK_GLOBAL(zggev,ZGGEV)
#define LAPACK_sggevx LAPACK_GLOBAL(sggevx,SGGEVX)
#define LAPACK_dggevx LAPACK_GLOBAL(dggevx,DGGEVX)
#define LAPACK_cggevx LAPACK_GLOBAL(cggevx,CGGEVX)
#define LAPACK_zggevx LAPACK_GLOBAL(zggevx,ZGGEVX)
#define LAPACK_dsfrk LAPACK_GLOBAL(dsfrk,DSFRK)
#define LAPACK_ssfrk LAPACK_GLOBAL(ssfrk,SSFRK)
#define LAPACK_zhfrk LAPACK_GLOBAL(zhfrk,ZHFRK)
#define LAPACK_chfrk LAPACK_GLOBAL(chfrk,CHFRK)
#define LAPACK_dtfsm LAPACK_GLOBAL(dtfsm,DTFSM)
#define LAPACK_stfsm LAPACK_GLOBAL(stfsm,STFSM)
#define LAPACK_ztfsm LAPACK_GLOBAL(ztfsm,ZTFSM)
#define LAPACK_ctfsm LAPACK_GLOBAL(ctfsm,CTFSM)
#define LAPACK_dtfttp LAPACK_GLOBAL(dtfttp,DTFTTP)
#define LAPACK_stfttp LAPACK_GLOBAL(stfttp,STFTTP)
#define LAPACK_ztfttp LAPACK_GLOBAL(ztfttp,ZTFTTP)
#define LAPACK_ctfttp LAPACK_GLOBAL(ctfttp,CTFTTP)
#define LAPACK_dtfttr LAPACK_GLOBAL(dtfttr,DTFTTR)
#define LAPACK_stfttr LAPACK_GLOBAL(stfttr,STFTTR)
#define LAPACK_ztfttr LAPACK_GLOBAL(ztfttr,ZTFTTR)
#define LAPACK_ctfttr LAPACK_GLOBAL(ctfttr,CTFTTR)
#define LAPACK_dtpttf LAPACK_GLOBAL(dtpttf,DTPTTF)
#define LAPACK_stpttf LAPACK_GLOBAL(stpttf,STPTTF)
#define LAPACK_ztpttf LAPACK_GLOBAL(ztpttf,ZTPTTF)
#define LAPACK_ctpttf LAPACK_GLOBAL(ctpttf,CTPTTF)
#define LAPACK_dtpttr LAPACK_GLOBAL(dtpttr,DTPTTR)
#define LAPACK_stpttr LAPACK_GLOBAL(stpttr,STPTTR)
#define LAPACK_ztpttr LAPACK_GLOBAL(ztpttr,ZTPTTR)
#define LAPACK_ctpttr LAPACK_GLOBAL(ctpttr,CTPTTR)
#define LAPACK_dtrttf LAPACK_GLOBAL(dtrttf,DTRTTF)
#define LAPACK_strttf LAPACK_GLOBAL(strttf,STRTTF)
#define LAPACK_ztrttf LAPACK_GLOBAL(ztrttf,ZTRTTF)
#define LAPACK_ctrttf LAPACK_GLOBAL(ctrttf,CTRTTF)
#define LAPACK_dtrttp LAPACK_GLOBAL(dtrttp,DTRTTP)
#define LAPACK_strttp LAPACK_GLOBAL(strttp,STRTTP)
#define LAPACK_ztrttp LAPACK_GLOBAL(ztrttp,ZTRTTP)
#define LAPACK_ctrttp LAPACK_GLOBAL(ctrttp,CTRTTP)
#define LAPACK_sgeqrfp LAPACK_GLOBAL(sgeqrfp,SGEQRFP)
#define LAPACK_dgeqrfp LAPACK_GLOBAL(dgeqrfp,DGEQRFP)
#define LAPACK_cgeqrfp LAPACK_GLOBAL(cgeqrfp,CGEQRFP)
#define LAPACK_zgeqrfp LAPACK_GLOBAL(zgeqrfp,ZGEQRFP)
#define LAPACK_clacgv LAPACK_GLOBAL(clacgv,CLACGV)
#define LAPACK_zlacgv LAPACK_GLOBAL(zlacgv,ZLACGV)
#define LAPACK_slarnv LAPACK_GLOBAL(slarnv,SLARNV)
#define LAPACK_dlarnv LAPACK_GLOBAL(dlarnv,DLARNV)
#define LAPACK_clarnv LAPACK_GLOBAL(clarnv,CLARNV)
#define LAPACK_zlarnv LAPACK_GLOBAL(zlarnv,ZLARNV)
#define LAPACK_sgeqr2 LAPACK_GLOBAL(sgeqr2,SGEQR2)
#define LAPACK_dgeqr2 LAPACK_GLOBAL(dgeqr2,DGEQR2)
#define LAPACK_cgeqr2 LAPACK_GLOBAL(cgeqr2,CGEQR2)
#define LAPACK_zgeqr2 LAPACK_GLOBAL(zgeqr2,ZGEQR2)
#define LAPACK_slacn2 LAPACK_GLOBAL(slacn2,SLACN2)
#define LAPACK_dlacn2 LAPACK_GLOBAL(dlacn2,DLACN2)
#define LAPACK_clacn2 LAPACK_GLOBAL(clacn2,CLACN2)
#define LAPACK_zlacn2 LAPACK_GLOBAL(zlacn2,ZLACN2)
#define LAPACK_slacpy LAPACK_GLOBAL(slacpy,SLACPY)
#define LAPACK_dlacpy LAPACK_GLOBAL(dlacpy,DLACPY)
#define LAPACK_clacpy LAPACK_GLOBAL(clacpy,CLACPY)
#define LAPACK_zlacpy LAPACK_GLOBAL(zlacpy,ZLACPY)
#define LAPACK_clacp2 LAPACK_GLOBAL(clacp2,CLACP2)
#define LAPACK_zlacp2 LAPACK_GLOBAL(zlacp2,ZLACP2)
#define LAPACK_sgetf2 LAPACK_GLOBAL(sgetf2,SGETF2)
#define LAPACK_dgetf2 LAPACK_GLOBAL(dgetf2,DGETF2)
#define LAPACK_cgetf2 LAPACK_GLOBAL(cgetf2,CGETF2)
#define LAPACK_zgetf2 LAPACK_GLOBAL(zgetf2,ZGETF2)
#define LAPACK_slaswp LAPACK_GLOBAL(slaswp,SLASWP)
#define LAPACK_dlaswp LAPACK_GLOBAL(dlaswp,DLASWP)
#define LAPACK_claswp LAPACK_GLOBAL(claswp,CLASWP)
#define LAPACK_zlaswp LAPACK_GLOBAL(zlaswp,ZLASWP)
#define LAPACK_slange LAPACK_GLOBAL(slange,SLANGE)
#define LAPACK_dlange LAPACK_GLOBAL(dlange,DLANGE)
#define LAPACK_clange LAPACK_GLOBAL(clange,CLANGE)
#define LAPACK_zlange LAPACK_GLOBAL(zlange,ZLANGE)
#define LAPACK_clanhe LAPACK_GLOBAL(clanhe,CLANHE)
#define LAPACK_zlanhe LAPACK_GLOBAL(zlanhe,ZLANHE)
#define LAPACK_slansy LAPACK_GLOBAL(slansy,SLANSY)
#define LAPACK_dlansy LAPACK_GLOBAL(dlansy,DLANSY)
#define LAPACK_clansy LAPACK_GLOBAL(clansy,CLANSY)
#define LAPACK_zlansy LAPACK_GLOBAL(zlansy,ZLANSY)
#define LAPACK_slantr LAPACK_GLOBAL(slantr,SLANTR)
#define LAPACK_dlantr LAPACK_GLOBAL(dlantr,DLANTR)
#define LAPACK_clantr LAPACK_GLOBAL(clantr,CLANTR)
#define LAPACK_zlantr LAPACK_GLOBAL(zlantr,ZLANTR)
#define LAPACK_slamch LAPACK_GLOBAL(slamch,SLAMCH)
#define LAPACK_dlamch LAPACK_GLOBAL(dlamch,DLAMCH)
#define LAPACK_sgelq2 LAPACK_GLOBAL(sgelq2,SGELQ2)
#define LAPACK_dgelq2 LAPACK_GLOBAL(dgelq2,DGELQ2)
#define LAPACK_cgelq2 LAPACK_GLOBAL(cgelq2,CGELQ2)
#define LAPACK_zgelq2 LAPACK_GLOBAL(zgelq2,ZGELQ2)
#define LAPACK_slarfb LAPACK_GLOBAL(slarfb,SLARFB)
#define LAPACK_dlarfb LAPACK_GLOBAL(dlarfb,DLARFB)
#define LAPACK_clarfb LAPACK_GLOBAL(clarfb,CLARFB)
#define LAPACK_zlarfb LAPACK_GLOBAL(zlarfb,ZLARFB)
#define LAPACK_slarfg LAPACK_GLOBAL(slarfg,SLARFG)
#define LAPACK_dlarfg LAPACK_GLOBAL(dlarfg,DLARFG)
#define LAPACK_clarfg LAPACK_GLOBAL(clarfg,CLARFG)
#define LAPACK_zlarfg LAPACK_GLOBAL(zlarfg,ZLARFG)
#define LAPACK_slarft LAPACK_GLOBAL(slarft,SLARFT)
#define LAPACK_dlarft LAPACK_GLOBAL(dlarft,DLARFT)
#define LAPACK_clarft LAPACK_GLOBAL(clarft,CLARFT)
#define LAPACK_zlarft LAPACK_GLOBAL(zlarft,ZLARFT)
#define LAPACK_slarfx LAPACK_GLOBAL(slarfx,SLARFX)
#define LAPACK_dlarfx LAPACK_GLOBAL(dlarfx,DLARFX)
#define LAPACK_clarfx LAPACK_GLOBAL(clarfx,CLARFX)
#define LAPACK_zlarfx LAPACK_GLOBAL(zlarfx,ZLARFX)
#define LAPACK_slatms LAPACK_GLOBAL(slatms,SLATMS)
#define LAPACK_dlatms LAPACK_GLOBAL(dlatms,DLATMS)
#define LAPACK_clatms LAPACK_GLOBAL(clatms,CLATMS)
#define LAPACK_zlatms LAPACK_GLOBAL(zlatms,ZLATMS)
#define LAPACK_slag2d LAPACK_GLOBAL(slag2d,SLAG2D)
#define LAPACK_dlag2s LAPACK_GLOBAL(dlag2s,DLAG2S)
#define LAPACK_clag2z LAPACK_GLOBAL(clag2z,CLAG2Z)
#define LAPACK_zlag2c LAPACK_GLOBAL(zlag2c,ZLAG2C)
#define LAPACK_slauum LAPACK_GLOBAL(slauum,SLAUUM)
#define LAPACK_dlauum LAPACK_GLOBAL(dlauum,DLAUUM)
#define LAPACK_clauum LAPACK_GLOBAL(clauum,CLAUUM)
#define LAPACK_zlauum LAPACK_GLOBAL(zlauum,ZLAUUM)
#define LAPACK_slagge LAPACK_GLOBAL(slagge,SLAGGE)
#define LAPACK_dlagge LAPACK_GLOBAL(dlagge,DLAGGE)
#define LAPACK_clagge LAPACK_GLOBAL(clagge,CLAGGE)
#define LAPACK_zlagge LAPACK_GLOBAL(zlagge,ZLAGGE)
#define LAPACK_slaset LAPACK_GLOBAL(slaset,SLASET)
#define LAPACK_dlaset LAPACK_GLOBAL(dlaset,DLASET)
#define LAPACK_claset LAPACK_GLOBAL(claset,CLASET)
#define LAPACK_zlaset LAPACK_GLOBAL(zlaset,ZLASET)
#define LAPACK_slasrt LAPACK_GLOBAL(slasrt,SLASRT)
#define LAPACK_dlasrt LAPACK_GLOBAL(dlasrt,DLASRT)
#define LAPACK_slagsy LAPACK_GLOBAL(slagsy,SLAGSY)
#define LAPACK_dlagsy LAPACK_GLOBAL(dlagsy,DLAGSY)
#define LAPACK_clagsy LAPACK_GLOBAL(clagsy,CLAGSY)
#define LAPACK_zlagsy LAPACK_GLOBAL(zlagsy,ZLAGSY)
#define LAPACK_claghe LAPACK_GLOBAL(claghe,CLAGHE)
#define LAPACK_zlaghe LAPACK_GLOBAL(zlaghe,ZLAGHE)
#define LAPACK_slapmr LAPACK_GLOBAL(slapmr,SLAPMR)
#define LAPACK_dlapmr LAPACK_GLOBAL(dlapmr,DLAPMR)
#define LAPACK_clapmr LAPACK_GLOBAL(clapmr,CLAPMR)
#define LAPACK_zlapmr LAPACK_GLOBAL(zlapmr,ZLAPMR)
#define LAPACK_slapy2 LAPACK_GLOBAL(slapy2,SLAPY2)
#define LAPACK_dlapy2 LAPACK_GLOBAL(dlapy2,DLAPY2)
#define LAPACK_slapy3 LAPACK_GLOBAL(slapy3,SLAPY3)
#define LAPACK_dlapy3 LAPACK_GLOBAL(dlapy3,DLAPY3)
#define LAPACK_slartgp LAPACK_GLOBAL(slartgp,SLARTGP)
#define LAPACK_dlartgp LAPACK_GLOBAL(dlartgp,DLARTGP)
#define LAPACK_slartgs LAPACK_GLOBAL(slartgs,SLARTGS)
#define LAPACK_dlartgs LAPACK_GLOBAL(dlartgs,DLARTGS)
// LAPACK 3.3.0
#define LAPACK_cbbcsd LAPACK_GLOBAL(cbbcsd,CBBCSD)
#define LAPACK_cheswapr LAPACK_GLOBAL(cheswapr,CHESWAPR)
#define LAPACK_chetri2 LAPACK_GLOBAL(chetri2,CHETRI2)
#define LAPACK_chetri2x LAPACK_GLOBAL(chetri2x,CHETRI2X)
#define LAPACK_chetrs2 LAPACK_GLOBAL(chetrs2,CHETRS2)
#define LAPACK_csyconv LAPACK_GLOBAL(csyconv,CSYCONV)
#define LAPACK_csyswapr LAPACK_GLOBAL(csyswapr,CSYSWAPR)
#define LAPACK_csytri2 LAPACK_GLOBAL(csytri2,CSYTRI2)
#define LAPACK_csytri2x LAPACK_GLOBAL(csytri2x,CSYTRI2X)
#define LAPACK_csytrs2 LAPACK_GLOBAL(csytrs2,CSYTRS2)
#define LAPACK_cunbdb LAPACK_GLOBAL(cunbdb,CUNBDB)
#define LAPACK_cuncsd LAPACK_GLOBAL(cuncsd,CUNCSD)
#define LAPACK_dbbcsd LAPACK_GLOBAL(dbbcsd,DBBCSD)
#define LAPACK_dorbdb LAPACK_GLOBAL(dorbdb,DORBDB)
#define LAPACK_dorcsd LAPACK_GLOBAL(dorcsd,DORCSD)
#define LAPACK_dsyconv LAPACK_GLOBAL(dsyconv,DSYCONV)
#define LAPACK_dsyswapr LAPACK_GLOBAL(dsyswapr,DSYSWAPR)
#define LAPACK_dsytri2 LAPACK_GLOBAL(dsytri2,DSYTRI2)
#define LAPACK_dsytri2x LAPACK_GLOBAL(dsytri2x,DSYTRI2X)
#define LAPACK_dsytrs2 LAPACK_GLOBAL(dsytrs2,DSYTRS2)
#define LAPACK_sbbcsd LAPACK_GLOBAL(sbbcsd,SBBCSD)
#define LAPACK_sorbdb LAPACK_GLOBAL(sorbdb,SORBDB)
#define LAPACK_sorcsd LAPACK_GLOBAL(sorcsd,SORCSD)
#define LAPACK_ssyconv LAPACK_GLOBAL(ssyconv,SSYCONV)
#define LAPACK_ssyswapr LAPACK_GLOBAL(ssyswapr,SSYSWAPR)
#define LAPACK_ssytri2 LAPACK_GLOBAL(ssytri2,SSYTRI2)
#define LAPACK_ssytri2x LAPACK_GLOBAL(ssytri2x,SSYTRI2X)
#define LAPACK_ssytrs2 LAPACK_GLOBAL(ssytrs2,SSYTRS2)
#define LAPACK_zbbcsd LAPACK_GLOBAL(zbbcsd,ZBBCSD)
#define LAPACK_zheswapr LAPACK_GLOBAL(zheswapr,ZHESWAPR)
#define LAPACK_zhetri2 LAPACK_GLOBAL(zhetri2,ZHETRI2)
#define LAPACK_zhetri2x LAPACK_GLOBAL(zhetri2x,ZHETRI2X)
#define LAPACK_zhetrs2 LAPACK_GLOBAL(zhetrs2,ZHETRS2)
#define LAPACK_zsyconv LAPACK_GLOBAL(zsyconv,ZSYCONV)
#define LAPACK_zsyswapr LAPACK_GLOBAL(zsyswapr,ZSYSWAPR)
#define LAPACK_zsytri2 LAPACK_GLOBAL(zsytri2,ZSYTRI2)
#define LAPACK_zsytri2x LAPACK_GLOBAL(zsytri2x,ZSYTRI2X)
#define LAPACK_zsytrs2 LAPACK_GLOBAL(zsytrs2,ZSYTRS2)
#define LAPACK_zunbdb LAPACK_GLOBAL(zunbdb,ZUNBDB)
#define LAPACK_zuncsd LAPACK_GLOBAL(zuncsd,ZUNCSD)
// LAPACK 3.4.0
#define LAPACK_sgemqrt LAPACK_GLOBAL(sgemqrt,SGEMQRT)
#define LAPACK_dgemqrt LAPACK_GLOBAL(dgemqrt,DGEMQRT)
#define LAPACK_cgemqrt LAPACK_GLOBAL(cgemqrt,CGEMQRT)
#define LAPACK_zgemqrt LAPACK_GLOBAL(zgemqrt,ZGEMQRT)
#define LAPACK_sgeqrt LAPACK_GLOBAL(sgeqrt,SGEQRT)
#define LAPACK_dgeqrt LAPACK_GLOBAL(dgeqrt,DGEQRT)
#define LAPACK_cgeqrt LAPACK_GLOBAL(cgeqrt,CGEQRT)
#define LAPACK_zgeqrt LAPACK_GLOBAL(zgeqrt,ZGEQRT)
#define LAPACK_sgeqrt2 LAPACK_GLOBAL(sgeqrt2,SGEQRT2)
#define LAPACK_dgeqrt2 LAPACK_GLOBAL(dgeqrt2,DGEQRT2)
#define LAPACK_cgeqrt2 LAPACK_GLOBAL(cgeqrt2,CGEQRT2)
#define LAPACK_zgeqrt2 LAPACK_GLOBAL(zgeqrt2,ZGEQRT2)
#define LAPACK_sgeqrt3 LAPACK_GLOBAL(sgeqrt3,SGEQRT3)
#define LAPACK_dgeqrt3 LAPACK_GLOBAL(dgeqrt3,DGEQRT3)
#define LAPACK_cgeqrt3 LAPACK_GLOBAL(cgeqrt3,CGEQRT3)
#define LAPACK_zgeqrt3 LAPACK_GLOBAL(zgeqrt3,ZGEQRT3)
#define LAPACK_stpmqrt LAPACK_GLOBAL(stpmqrt,STPMQRT)
#define LAPACK_dtpmqrt LAPACK_GLOBAL(dtpmqrt,DTPMQRT)
#define LAPACK_ctpmqrt LAPACK_GLOBAL(ctpmqrt,CTPMQRT)
#define LAPACK_ztpmqrt LAPACK_GLOBAL(ztpmqrt,ZTPMQRT)
#define LAPACK_dtpqrt LAPACK_GLOBAL(dtpqrt,DTPQRT)
#define LAPACK_ctpqrt LAPACK_GLOBAL(ctpqrt,CTPQRT)
#define LAPACK_ztpqrt LAPACK_GLOBAL(ztpqrt,ZTPQRT)
#define LAPACK_stpqrt2 LAPACK_GLOBAL(stpqrt2,STPQRT2)
#define LAPACK_dtpqrt2 LAPACK_GLOBAL(dtpqrt2,DTPQRT2)
#define LAPACK_ctpqrt2 LAPACK_GLOBAL(ctpqrt2,CTPQRT2)
#define LAPACK_ztpqrt2 LAPACK_GLOBAL(ztpqrt2,ZTPQRT2)
#define LAPACK_stprfb LAPACK_GLOBAL(stprfb,STPRFB)
#define LAPACK_dtprfb LAPACK_GLOBAL(dtprfb,DTPRFB)
#define LAPACK_ctprfb LAPACK_GLOBAL(ctprfb,CTPRFB)
#define LAPACK_ztprfb LAPACK_GLOBAL(ztprfb,ZTPRFB)
// LAPACK 3.X.X
#define LAPACK_ssysv_rook LAPACK_GLOBAL(ssysv_rook,SSYSV_ROOK)
#define LAPACK_dsysv_rook LAPACK_GLOBAL(dsysv_rook,DSYSV_ROOK)
#define LAPACK_csysv_rook LAPACK_GLOBAL(csysv_rook,CSYSV_ROOK)
#define LAPACK_zsysv_rook LAPACK_GLOBAL(zsysv_rook,ZSYSV_ROOK)
#define LAPACK_csyr LAPACK_GLOBAL(csyr,CSYR)
#define LAPACK_zsyr LAPACK_GLOBAL(zsyr,ZSYR)
#define LAPACK_ilaver LAPACK_GLOBAL(ilaver,ILAVER)

void LAPACK_sgetrf( lapack_int* m, lapack_int* n, float* a, lapack_int* lda,
					lapack_int* ipiv, lapack_int *info );
void LAPACK_dgetrf( lapack_int* m, lapack_int* n, double* a, lapack_int* lda,
					lapack_int* ipiv, lapack_int *info );
void LAPACK_cgetrf( lapack_int* m, lapack_int* n, lapack_complex_float* a,
					lapack_int* lda, lapack_int* ipiv, lapack_int *info );
void LAPACK_zgetrf( lapack_int* m, lapack_int* n, lapack_complex_double* a,
					lapack_int* lda, lapack_int* ipiv, lapack_int *info );
void LAPACK_sgbtrf( lapack_int* m, lapack_int* n, lapack_int* kl,
					lapack_int* ku, float* ab, lapack_int* ldab,
					lapack_int* ipiv, lapack_int *info );
void LAPACK_dgbtrf( lapack_int* m, lapack_int* n, lapack_int* kl,
					lapack_int* ku, double* ab, lapack_int* ldab,
					lapack_int* ipiv, lapack_int *info );
void LAPACK_cgbtrf( lapack_int* m, lapack_int* n, lapack_int* kl,
					lapack_int* ku, lapack_complex_float* ab, lapack_int* ldab,
					lapack_int* ipiv, lapack_int *info );
void LAPACK_zgbtrf( lapack_int* m, lapack_int* n, lapack_int* kl,
					lapack_int* ku, lapack_complex_double* ab, lapack_int* ldab,
					lapack_int* ipiv, lapack_int *info );
void LAPACK_sgttrf( lapack_int* n, float* dl, float* d, float* du, float* du2,
					lapack_int* ipiv, lapack_int *info );
void LAPACK_dgttrf( lapack_int* n, double* dl, double* d, double* du,
					double* du2, lapack_int* ipiv, lapack_int *info );
void LAPACK_cgttrf( lapack_int* n, lapack_complex_float* dl,
					lapack_complex_float* d, lapack_complex_float* du,
					lapack_complex_float* du2, lapack_int* ipiv,
					lapack_int *info );
void LAPACK_zgttrf( lapack_int* n, lapack_complex_double* dl,
					lapack_complex_double* d, lapack_complex_double* du,
					lapack_complex_double* du2, lapack_int* ipiv,
					lapack_int *info );
void LAPACK_spotrf( char* uplo, lapack_int* n, float* a, lapack_int* lda,
					lapack_int *info );
void LAPACK_dpotrf( char* uplo, lapack_int* n, double* a, lapack_int* lda,
					lapack_int *info );
void LAPACK_cpotrf( char* uplo, lapack_int* n, lapack_complex_float* a,
					lapack_int* lda, lapack_int *info );
void LAPACK_zpotrf( char* uplo, lapack_int* n, lapack_complex_double* a,
					lapack_int* lda, lapack_int *info );
void LAPACK_dpstrf( char* uplo, lapack_int* n, double* a, lapack_int* lda,
					lapack_int* piv, lapack_int* rank, double* tol,
					double* work, lapack_int *info );
void LAPACK_spstrf( char* uplo, lapack_int* n, float* a, lapack_int* lda,
					lapack_int* piv, lapack_int* rank, float* tol, float* work,
					lapack_int *info );
void LAPACK_zpstrf( char* uplo, lapack_int* n, lapack_complex_double* a,
					lapack_int* lda, lapack_int* piv, lapack_int* rank,
					double* tol, double* work, lapack_int *info );
void LAPACK_cpstrf( char* uplo, lapack_int* n, lapack_complex_float* a,
					lapack_int* lda, lapack_int* piv, lapack_int* rank,
					float* tol, float* work, lapack_int *info );
void LAPACK_dpftrf( char* transr, char* uplo, lapack_int* n, double* a,
					lapack_int *info );
void LAPACK_spftrf( char* transr, char* uplo, lapack_int* n, float* a,
					lapack_int *info );
void LAPACK_zpftrf( char* transr, char* uplo, lapack_int* n,
					lapack_complex_double* a, lapack_int *info );
void LAPACK_cpftrf( char* transr, char* uplo, lapack_int* n,
					lapack_complex_float* a, lapack_int *info );
void LAPACK_spptrf( char* uplo, lapack_int* n, float* ap, lapack_int *info );
void LAPACK_dpptrf( char* uplo, lapack_int* n, double* ap, lapack_int *info );
void LAPACK_cpptrf( char* uplo, lapack_int* n, lapack_complex_float* ap,
					lapack_int *info );
void LAPACK_zpptrf( char* uplo, lapack_int* n, lapack_complex_double* ap,
					lapack_int *info );
void LAPACK_spbtrf( char* uplo, lapack_int* n, lapack_int* kd, float* ab,
					lapack_int* ldab, lapack_int *info );
void LAPACK_dpbtrf( char* uplo, lapack_int* n, lapack_int* kd, double* ab,
					lapack_int* ldab, lapack_int *info );
void LAPACK_cpbtrf( char* uplo, lapack_int* n, lapack_int* kd,
					lapack_complex_float* ab, lapack_int* ldab,
					lapack_int *info );
void LAPACK_zpbtrf( char* uplo, lapack_int* n, lapack_int* kd,
					lapack_complex_double* ab, lapack_int* ldab,
					lapack_int *info );
void LAPACK_spttrf( lapack_int* n, float* d, float* e, lapack_int *info );
void LAPACK_dpttrf( lapack_int* n, double* d, double* e, lapack_int *info );
void LAPACK_cpttrf( lapack_int* n, float* d, lapack_complex_float* e,
					lapack_int *info );
void LAPACK_zpttrf( lapack_int* n, double* d, lapack_complex_double* e,
					lapack_int *info );
void LAPACK_ssytrf( char* uplo, lapack_int* n, float* a, lapack_int* lda,
					lapack_int* ipiv, float* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_dsytrf( char* uplo, lapack_int* n, double* a, lapack_int* lda,
					lapack_int* ipiv, double* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_csytrf( char* uplo, lapack_int* n, lapack_complex_float* a,
					lapack_int* lda, lapack_int* ipiv,
					lapack_complex_float* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_zsytrf( char* uplo, lapack_int* n, lapack_complex_double* a,
					lapack_int* lda, lapack_int* ipiv,
					lapack_complex_double* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_chetrf( char* uplo, lapack_int* n, lapack_complex_float* a,
					lapack_int* lda, lapack_int* ipiv,
					lapack_complex_float* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_zhetrf( char* uplo, lapack_int* n, lapack_complex_double* a,
					lapack_int* lda, lapack_int* ipiv,
					lapack_complex_double* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_ssptrf( char* uplo, lapack_int* n, float* ap, lapack_int* ipiv,
					lapack_int *info );
void LAPACK_dsptrf( char* uplo, lapack_int* n, double* ap, lapack_int* ipiv,
					lapack_int *info );
void LAPACK_csptrf( char* uplo, lapack_int* n, lapack_complex_float* ap,
					lapack_int* ipiv, lapack_int *info );
void LAPACK_zsptrf( char* uplo, lapack_int* n, lapack_complex_double* ap,
					lapack_int* ipiv, lapack_int *info );
void LAPACK_chptrf( char* uplo, lapack_int* n, lapack_complex_float* ap,
					lapack_int* ipiv, lapack_int *info );
void LAPACK_zhptrf( char* uplo, lapack_int* n, lapack_complex_double* ap,
					lapack_int* ipiv, lapack_int *info );
void LAPACK_sgetrs( char* trans, lapack_int* n, lapack_int* nrhs,
					const float* a, lapack_int* lda, const lapack_int* ipiv,
					float* b, lapack_int* ldb, lapack_int *info );
void LAPACK_dgetrs( char* trans, lapack_int* n, lapack_int* nrhs,
					const double* a, lapack_int* lda, const lapack_int* ipiv,
					double* b, lapack_int* ldb, lapack_int *info );
void LAPACK_cgetrs( char* trans, lapack_int* n, lapack_int* nrhs,
					const lapack_complex_float* a, lapack_int* lda,
					const lapack_int* ipiv, lapack_complex_float* b,
					lapack_int* ldb, lapack_int *info );
void LAPACK_zgetrs( char* trans, lapack_int* n, lapack_int* nrhs,
					const lapack_complex_double* a, lapack_int* lda,
					const lapack_int* ipiv, lapack_complex_double* b,
					lapack_int* ldb, lapack_int *info );
void LAPACK_sgbtrs( char* trans, lapack_int* n, lapack_int* kl, lapack_int* ku,
					lapack_int* nrhs, const float* ab, lapack_int* ldab,
					const lapack_int* ipiv, float* b, lapack_int* ldb,
					lapack_int *info );
void LAPACK_dgbtrs( char* trans, lapack_int* n, lapack_int* kl, lapack_int* ku,
					lapack_int* nrhs, const double* ab, lapack_int* ldab,
					const lapack_int* ipiv, double* b, lapack_int* ldb,
					lapack_int *info );
void LAPACK_cgbtrs( char* trans, lapack_int* n, lapack_int* kl, lapack_int* ku,
					lapack_int* nrhs, const lapack_complex_float* ab,
					lapack_int* ldab, const lapack_int* ipiv,
					lapack_complex_float* b, lapack_int* ldb,
					lapack_int *info );
void LAPACK_zgbtrs( char* trans, lapack_int* n, lapack_int* kl, lapack_int* ku,
					lapack_int* nrhs, const lapack_complex_double* ab,
					lapack_int* ldab, const lapack_int* ipiv,
					lapack_complex_double* b, lapack_int* ldb,
					lapack_int *info );
void LAPACK_sgttrs( char* trans, lapack_int* n, lapack_int* nrhs,
					const float* dl, const float* d, const float* du,
					const float* du2, const lapack_int* ipiv, float* b,
					lapack_int* ldb, lapack_int *info );
void LAPACK_dgttrs( char* trans, lapack_int* n, lapack_int* nrhs,
					const double* dl, const double* d, const double* du,
					const double* du2, const lapack_int* ipiv, double* b,
					lapack_int* ldb, lapack_int *info );
void LAPACK_cgttrs( char* trans, lapack_int* n, lapack_int* nrhs,
					const lapack_complex_float* dl,
					const lapack_complex_float* d,
					const lapack_complex_float* du,
					const lapack_complex_float* du2, const lapack_int* ipiv,
					lapack_complex_float* b, lapack_int* ldb,
					lapack_int *info );
void LAPACK_zgttrs( char* trans, lapack_int* n, lapack_int* nrhs,
					const lapack_complex_double* dl,
					const lapack_complex_double* d,
					const lapack_complex_double* du,
					const lapack_complex_double* du2, const lapack_int* ipiv,
					lapack_complex_double* b, lapack_int* ldb,
					lapack_int *info );
void LAPACK_spotrs( char* uplo, lapack_int* n, lapack_int* nrhs, const float* a,
					lapack_int* lda, float* b, lapack_int* ldb,
					lapack_int *info );
void LAPACK_dpotrs( char* uplo, lapack_int* n, lapack_int* nrhs,
					const double* a, lapack_int* lda, double* b,
					lapack_int* ldb, lapack_int *info );
void LAPACK_cpotrs( char* uplo, lapack_int* n, lapack_int* nrhs,
					const lapack_complex_float* a, lapack_int* lda,
					lapack_complex_float* b, lapack_int* ldb,
					lapack_int *info );
void LAPACK_zpotrs( char* uplo, lapack_int* n, lapack_int* nrhs,
					const lapack_complex_double* a, lapack_int* lda,
					lapack_complex_double* b, lapack_int* ldb,
					lapack_int *info );
void LAPACK_dpftrs( char* transr, char* uplo, lapack_int* n, lapack_int* nrhs,
					const double* a, double* b, lapack_int* ldb,
					lapack_int *info );
void LAPACK_spftrs( char* transr, char* uplo, lapack_int* n, lapack_int* nrhs,
					const float* a, float* b, lapack_int* ldb,
					lapack_int *info );
void LAPACK_zpftrs( char* transr, char* uplo, lapack_int* n, lapack_int* nrhs,
					const lapack_complex_double* a, lapack_complex_double* b,
					lapack_int* ldb, lapack_int *info );
void LAPACK_cpftrs( char* transr, char* uplo, lapack_int* n, lapack_int* nrhs,
					const lapack_complex_float* a, lapack_complex_float* b,
					lapack_int* ldb, lapack_int *info );
void LAPACK_spptrs( char* uplo, lapack_int* n, lapack_int* nrhs,
					const float* ap, float* b, lapack_int* ldb,
					lapack_int *info );
void LAPACK_dpptrs( char* uplo, lapack_int* n, lapack_int* nrhs,
					const double* ap, double* b, lapack_int* ldb,
					lapack_int *info );
void LAPACK_cpptrs( char* uplo, lapack_int* n, lapack_int* nrhs,
					const lapack_complex_float* ap, lapack_complex_float* b,
					lapack_int* ldb, lapack_int *info );
void LAPACK_zpptrs( char* uplo, lapack_int* n, lapack_int* nrhs,
					const lapack_complex_double* ap, lapack_complex_double* b,
					lapack_int* ldb, lapack_int *info );
void LAPACK_spbtrs( char* uplo, lapack_int* n, lapack_int* kd, lapack_int* nrhs,
					const float* ab, lapack_int* ldab, float* b,
					lapack_int* ldb, lapack_int *info );
void LAPACK_dpbtrs( char* uplo, lapack_int* n, lapack_int* kd, lapack_int* nrhs,
					const double* ab, lapack_int* ldab, double* b,
					lapack_int* ldb, lapack_int *info );
void LAPACK_cpbtrs( char* uplo, lapack_int* n, lapack_int* kd, lapack_int* nrhs,
					const lapack_complex_float* ab, lapack_int* ldab,
					lapack_complex_float* b, lapack_int* ldb,
					lapack_int *info );
void LAPACK_zpbtrs( char* uplo, lapack_int* n, lapack_int* kd, lapack_int* nrhs,
					const lapack_complex_double* ab, lapack_int* ldab,
					lapack_complex_double* b, lapack_int* ldb,
					lapack_int *info );
void LAPACK_spttrs( lapack_int* n, lapack_int* nrhs, const float* d,
					const float* e, float* b, lapack_int* ldb,
					lapack_int *info );
void LAPACK_dpttrs( lapack_int* n, lapack_int* nrhs, const double* d,
					const double* e, double* b, lapack_int* ldb,
					lapack_int *info );
void LAPACK_cpttrs( char* uplo, lapack_int* n, lapack_int* nrhs, const float* d,
					const lapack_complex_float* e, lapack_complex_float* b,
					lapack_int* ldb, lapack_int *info );
void LAPACK_zpttrs( char* uplo, lapack_int* n, lapack_int* nrhs,
					const double* d, const lapack_complex_double* e,
					lapack_complex_double* b, lapack_int* ldb,
					lapack_int *info );
void LAPACK_ssytrs( char* uplo, lapack_int* n, lapack_int* nrhs, const float* a,
					lapack_int* lda, const lapack_int* ipiv, float* b,
					lapack_int* ldb, lapack_int *info );
void LAPACK_dsytrs( char* uplo, lapack_int* n, lapack_int* nrhs,
					const double* a, lapack_int* lda, const lapack_int* ipiv,
					double* b, lapack_int* ldb, lapack_int *info );
void LAPACK_csytrs( char* uplo, lapack_int* n, lapack_int* nrhs,
					const lapack_complex_float* a, lapack_int* lda,
					const lapack_int* ipiv, lapack_complex_float* b,
					lapack_int* ldb, lapack_int *info );
void LAPACK_zsytrs( char* uplo, lapack_int* n, lapack_int* nrhs,
					const lapack_complex_double* a, lapack_int* lda,
					const lapack_int* ipiv, lapack_complex_double* b,
					lapack_int* ldb, lapack_int *info );
void LAPACK_chetrs( char* uplo, lapack_int* n, lapack_int* nrhs,
					const lapack_complex_float* a, lapack_int* lda,
					const lapack_int* ipiv, lapack_complex_float* b,
					lapack_int* ldb, lapack_int *info );
void LAPACK_zhetrs( char* uplo, lapack_int* n, lapack_int* nrhs,
					const lapack_complex_double* a, lapack_int* lda,
					const lapack_int* ipiv, lapack_complex_double* b,
					lapack_int* ldb, lapack_int *info );
void LAPACK_ssptrs( char* uplo, lapack_int* n, lapack_int* nrhs,
					const float* ap, const lapack_int* ipiv, float* b,
					lapack_int* ldb, lapack_int *info );
void LAPACK_dsptrs( char* uplo, lapack_int* n, lapack_int* nrhs,
					const double* ap, const lapack_int* ipiv, double* b,
					lapack_int* ldb, lapack_int *info );
void LAPACK_csptrs( char* uplo, lapack_int* n, lapack_int* nrhs,
					const lapack_complex_float* ap, const lapack_int* ipiv,
					lapack_complex_float* b, lapack_int* ldb,
					lapack_int *info );
void LAPACK_zsptrs( char* uplo, lapack_int* n, lapack_int* nrhs,
					const lapack_complex_double* ap, const lapack_int* ipiv,
					lapack_complex_double* b, lapack_int* ldb,
					lapack_int *info );
void LAPACK_chptrs( char* uplo, lapack_int* n, lapack_int* nrhs,
					const lapack_complex_float* ap, const lapack_int* ipiv,
					lapack_complex_float* b, lapack_int* ldb,
					lapack_int *info );
void LAPACK_zhptrs( char* uplo, lapack_int* n, lapack_int* nrhs,
					const lapack_complex_double* ap, const lapack_int* ipiv,
					lapack_complex_double* b, lapack_int* ldb,
					lapack_int *info );
void LAPACK_strtrs( char* uplo, char* trans, char* diag, lapack_int* n,
					lapack_int* nrhs, const float* a, lapack_int* lda, float* b,
					lapack_int* ldb, lapack_int *info );
void LAPACK_dtrtrs( char* uplo, char* trans, char* diag, lapack_int* n,
					lapack_int* nrhs, const double* a, lapack_int* lda,
					double* b, lapack_int* ldb, lapack_int *info );
void LAPACK_ctrtrs( char* uplo, char* trans, char* diag, lapack_int* n,
					lapack_int* nrhs, const lapack_complex_float* a,
					lapack_int* lda, lapack_complex_float* b, lapack_int* ldb,
					lapack_int *info );
void LAPACK_ztrtrs( char* uplo, char* trans, char* diag, lapack_int* n,
					lapack_int* nrhs, const lapack_complex_double* a,
					lapack_int* lda, lapack_complex_double* b, lapack_int* ldb,
					lapack_int *info );
void LAPACK_stptrs( char* uplo, char* trans, char* diag, lapack_int* n,
					lapack_int* nrhs, const float* ap, float* b,
					lapack_int* ldb, lapack_int *info );
void LAPACK_dtptrs( char* uplo, char* trans, char* diag, lapack_int* n,
					lapack_int* nrhs, const double* ap, double* b,
					lapack_int* ldb, lapack_int *info );
void LAPACK_ctptrs( char* uplo, char* trans, char* diag, lapack_int* n,
					lapack_int* nrhs, const lapack_complex_float* ap,
					lapack_complex_float* b, lapack_int* ldb,
					lapack_int *info );
void LAPACK_ztptrs( char* uplo, char* trans, char* diag, lapack_int* n,
					lapack_int* nrhs, const lapack_complex_double* ap,
					lapack_complex_double* b, lapack_int* ldb,
					lapack_int *info );
void LAPACK_stbtrs( char* uplo, char* trans, char* diag, lapack_int* n,
					lapack_int* kd, lapack_int* nrhs, const float* ab,
					lapack_int* ldab, float* b, lapack_int* ldb,
					lapack_int *info );
void LAPACK_dtbtrs( char* uplo, char* trans, char* diag, lapack_int* n,
					lapack_int* kd, lapack_int* nrhs, const double* ab,
					lapack_int* ldab, double* b, lapack_int* ldb,
					lapack_int *info );
void LAPACK_ctbtrs( char* uplo, char* trans, char* diag, lapack_int* n,
					lapack_int* kd, lapack_int* nrhs,
					const lapack_complex_float* ab, lapack_int* ldab,
					lapack_complex_float* b, lapack_int* ldb,
					lapack_int *info );
void LAPACK_ztbtrs( char* uplo, char* trans, char* diag, lapack_int* n,
					lapack_int* kd, lapack_int* nrhs,
					const lapack_complex_double* ab, lapack_int* ldab,
					lapack_complex_double* b, lapack_int* ldb,
					lapack_int *info );
void LAPACK_sgecon( char* norm, lapack_int* n, const float* a, lapack_int* lda,
					float* anorm, float* rcond, float* work, lapack_int* iwork,
					lapack_int *info );
void LAPACK_dgecon( char* norm, lapack_int* n, const double* a, lapack_int* lda,
					double* anorm, double* rcond, double* work,
					lapack_int* iwork, lapack_int *info );
void LAPACK_cgecon( char* norm, lapack_int* n, const lapack_complex_float* a,
					lapack_int* lda, float* anorm, float* rcond,
					lapack_complex_float* work, float* rwork,
					lapack_int *info );
void LAPACK_zgecon( char* norm, lapack_int* n, const lapack_complex_double* a,
					lapack_int* lda, double* anorm, double* rcond,
					lapack_complex_double* work, double* rwork,
					lapack_int *info );
void LAPACK_sgbcon( char* norm, lapack_int* n, lapack_int* kl, lapack_int* ku,
					const float* ab, lapack_int* ldab, const lapack_int* ipiv,
					float* anorm, float* rcond, float* work, lapack_int* iwork,
					lapack_int *info );
void LAPACK_dgbcon( char* norm, lapack_int* n, lapack_int* kl, lapack_int* ku,
					const double* ab, lapack_int* ldab, const lapack_int* ipiv,
					double* anorm, double* rcond, double* work,
					lapack_int* iwork, lapack_int *info );
void LAPACK_cgbcon( char* norm, lapack_int* n, lapack_int* kl, lapack_int* ku,
					const lapack_complex_float* ab, lapack_int* ldab,
					const lapack_int* ipiv, float* anorm, float* rcond,
					lapack_complex_float* work, float* rwork,
					lapack_int *info );
void LAPACK_zgbcon( char* norm, lapack_int* n, lapack_int* kl, lapack_int* ku,
					const lapack_complex_double* ab, lapack_int* ldab,
					const lapack_int* ipiv, double* anorm, double* rcond,
					lapack_complex_double* work, double* rwork,
					lapack_int *info );
void LAPACK_sgtcon( char* norm, lapack_int* n, const float* dl, const float* d,
					const float* du, const float* du2, const lapack_int* ipiv,
					float* anorm, float* rcond, float* work, lapack_int* iwork,
					lapack_int *info );
void LAPACK_dgtcon( char* norm, lapack_int* n, const double* dl,
					const double* d, const double* du, const double* du2,
					const lapack_int* ipiv, double* anorm, double* rcond,
					double* work, lapack_int* iwork, lapack_int *info );
void LAPACK_cgtcon( char* norm, lapack_int* n, const lapack_complex_float* dl,
					const lapack_complex_float* d,
					const lapack_complex_float* du,
					const lapack_complex_float* du2, const lapack_int* ipiv,
					float* anorm, float* rcond, lapack_complex_float* work,
					lapack_int *info );
void LAPACK_zgtcon( char* norm, lapack_int* n, const lapack_complex_double* dl,
					const lapack_complex_double* d,
					const lapack_complex_double* du,
					const lapack_complex_double* du2, const lapack_int* ipiv,
					double* anorm, double* rcond, lapack_complex_double* work,
					lapack_int *info );
void LAPACK_spocon( char* uplo, lapack_int* n, const float* a, lapack_int* lda,
					float* anorm, float* rcond, float* work, lapack_int* iwork,
					lapack_int *info );
void LAPACK_dpocon( char* uplo, lapack_int* n, const double* a, lapack_int* lda,
					double* anorm, double* rcond, double* work,
					lapack_int* iwork, lapack_int *info );
void LAPACK_cpocon( char* uplo, lapack_int* n, const lapack_complex_float* a,
					lapack_int* lda, float* anorm, float* rcond,
					lapack_complex_float* work, float* rwork,
					lapack_int *info );
void LAPACK_zpocon( char* uplo, lapack_int* n, const lapack_complex_double* a,
					lapack_int* lda, double* anorm, double* rcond,
					lapack_complex_double* work, double* rwork,
					lapack_int *info );
void LAPACK_sppcon( char* uplo, lapack_int* n, const float* ap, float* anorm,
					float* rcond, float* work, lapack_int* iwork,
					lapack_int *info );
void LAPACK_dppcon( char* uplo, lapack_int* n, const double* ap, double* anorm,
					double* rcond, double* work, lapack_int* iwork,
					lapack_int *info );
void LAPACK_cppcon( char* uplo, lapack_int* n, const lapack_complex_float* ap,
					float* anorm, float* rcond, lapack_complex_float* work,
					float* rwork, lapack_int *info );
void LAPACK_zppcon( char* uplo, lapack_int* n, const lapack_complex_double* ap,
					double* anorm, double* rcond, lapack_complex_double* work,
					double* rwork, lapack_int *info );
void LAPACK_spbcon( char* uplo, lapack_int* n, lapack_int* kd, const float* ab,
					lapack_int* ldab, float* anorm, float* rcond, float* work,
					lapack_int* iwork, lapack_int *info );
void LAPACK_dpbcon( char* uplo, lapack_int* n, lapack_int* kd, const double* ab,
					lapack_int* ldab, double* anorm, double* rcond,
					double* work, lapack_int* iwork, lapack_int *info );
void LAPACK_cpbcon( char* uplo, lapack_int* n, lapack_int* kd,
					const lapack_complex_float* ab, lapack_int* ldab,
					float* anorm, float* rcond, lapack_complex_float* work,
					float* rwork, lapack_int *info );
void LAPACK_zpbcon( char* uplo, lapack_int* n, lapack_int* kd,
					const lapack_complex_double* ab, lapack_int* ldab,
					double* anorm, double* rcond, lapack_complex_double* work,
					double* rwork, lapack_int *info );
void LAPACK_sptcon( lapack_int* n, const float* d, const float* e, float* anorm,
					float* rcond, float* work, lapack_int *info );
void LAPACK_dptcon( lapack_int* n, const double* d, const double* e,
					double* anorm, double* rcond, double* work,
					lapack_int *info );
void LAPACK_cptcon( lapack_int* n, const float* d,
					const lapack_complex_float* e, float* anorm, float* rcond,
					float* work, lapack_int *info );
void LAPACK_zptcon( lapack_int* n, const double* d,
					const lapack_complex_double* e, double* anorm,
					double* rcond, double* work, lapack_int *info );
void LAPACK_ssycon( char* uplo, lapack_int* n, const float* a, lapack_int* lda,
					const lapack_int* ipiv, float* anorm, float* rcond,
					float* work, lapack_int* iwork, lapack_int *info );
void LAPACK_dsycon( char* uplo, lapack_int* n, const double* a, lapack_int* lda,
					const lapack_int* ipiv, double* anorm, double* rcond,
					double* work, lapack_int* iwork, lapack_int *info );
void LAPACK_csycon( char* uplo, lapack_int* n, const lapack_complex_float* a,
					lapack_int* lda, const lapack_int* ipiv, float* anorm,
					float* rcond, lapack_complex_float* work,
					lapack_int *info );
void LAPACK_zsycon( char* uplo, lapack_int* n, const lapack_complex_double* a,
					lapack_int* lda, const lapack_int* ipiv, double* anorm,
					double* rcond, lapack_complex_double* work,
					lapack_int *info );
void LAPACK_checon( char* uplo, lapack_int* n, const lapack_complex_float* a,
					lapack_int* lda, const lapack_int* ipiv, float* anorm,
					float* rcond, lapack_complex_float* work,
					lapack_int *info );
void LAPACK_zhecon( char* uplo, lapack_int* n, const lapack_complex_double* a,
					lapack_int* lda, const lapack_int* ipiv, double* anorm,
					double* rcond, lapack_complex_double* work,
					lapack_int *info );
void LAPACK_sspcon( char* uplo, lapack_int* n, const float* ap,
					const lapack_int* ipiv, float* anorm, float* rcond,
					float* work, lapack_int* iwork, lapack_int *info );
void LAPACK_dspcon( char* uplo, lapack_int* n, const double* ap,
					const lapack_int* ipiv, double* anorm, double* rcond,
					double* work, lapack_int* iwork, lapack_int *info );
void LAPACK_cspcon( char* uplo, lapack_int* n, const lapack_complex_float* ap,
					const lapack_int* ipiv, float* anorm, float* rcond,
					lapack_complex_float* work, lapack_int *info );
void LAPACK_zspcon( char* uplo, lapack_int* n, const lapack_complex_double* ap,
					const lapack_int* ipiv, double* anorm, double* rcond,
					lapack_complex_double* work, lapack_int *info );
void LAPACK_chpcon( char* uplo, lapack_int* n, const lapack_complex_float* ap,
					const lapack_int* ipiv, float* anorm, float* rcond,
					lapack_complex_float* work, lapack_int *info );
void LAPACK_zhpcon( char* uplo, lapack_int* n, const lapack_complex_double* ap,
					const lapack_int* ipiv, double* anorm, double* rcond,
					lapack_complex_double* work, lapack_int *info );
void LAPACK_strcon( char* norm, char* uplo, char* diag, lapack_int* n,
					const float* a, lapack_int* lda, float* rcond, float* work,
					lapack_int* iwork, lapack_int *info );
void LAPACK_dtrcon( char* norm, char* uplo, char* diag, lapack_int* n,
					const double* a, lapack_int* lda, double* rcond,
					double* work, lapack_int* iwork, lapack_int *info );
void LAPACK_ctrcon( char* norm, char* uplo, char* diag, lapack_int* n,
					const lapack_complex_float* a, lapack_int* lda,
					float* rcond, lapack_complex_float* work, float* rwork,
					lapack_int *info );
void LAPACK_ztrcon( char* norm, char* uplo, char* diag, lapack_int* n,
					const lapack_complex_double* a, lapack_int* lda,
					double* rcond, lapack_complex_double* work, double* rwork,
					lapack_int *info );
void LAPACK_stpcon( char* norm, char* uplo, char* diag, lapack_int* n,
					const float* ap, float* rcond, float* work,
					lapack_int* iwork, lapack_int *info );
void LAPACK_dtpcon( char* norm, char* uplo, char* diag, lapack_int* n,
					const double* ap, double* rcond, double* work,
					lapack_int* iwork, lapack_int *info );
void LAPACK_ctpcon( char* norm, char* uplo, char* diag, lapack_int* n,
					const lapack_complex_float* ap, float* rcond,
					lapack_complex_float* work, float* rwork,
					lapack_int *info );
void LAPACK_ztpcon( char* norm, char* uplo, char* diag, lapack_int* n,
					const lapack_complex_double* ap, double* rcond,
					lapack_complex_double* work, double* rwork,
					lapack_int *info );
void LAPACK_stbcon( char* norm, char* uplo, char* diag, lapack_int* n,
					lapack_int* kd, const float* ab, lapack_int* ldab,
					float* rcond, float* work, lapack_int* iwork,
					lapack_int *info );
void LAPACK_dtbcon( char* norm, char* uplo, char* diag, lapack_int* n,
					lapack_int* kd, const double* ab, lapack_int* ldab,
					double* rcond, double* work, lapack_int* iwork,
					lapack_int *info );
void LAPACK_ctbcon( char* norm, char* uplo, char* diag, lapack_int* n,
					lapack_int* kd, const lapack_complex_float* ab,
					lapack_int* ldab, float* rcond, lapack_complex_float* work,
					float* rwork, lapack_int *info );
void LAPACK_ztbcon( char* norm, char* uplo, char* diag, lapack_int* n,
					lapack_int* kd, const lapack_complex_double* ab,
					lapack_int* ldab, double* rcond,
					lapack_complex_double* work, double* rwork,
					lapack_int *info );
void LAPACK_sgerfs( char* trans, lapack_int* n, lapack_int* nrhs,
					const float* a, lapack_int* lda, const float* af,
					lapack_int* ldaf, const lapack_int* ipiv, const float* b,
					lapack_int* ldb, float* x, lapack_int* ldx, float* ferr,
					float* berr, float* work, lapack_int* iwork,
					lapack_int *info );
void LAPACK_dgerfs( char* trans, lapack_int* n, lapack_int* nrhs,
					const double* a, lapack_int* lda, const double* af,
					lapack_int* ldaf, const lapack_int* ipiv, const double* b,
					lapack_int* ldb, double* x, lapack_int* ldx, double* ferr,
					double* berr, double* work, lapack_int* iwork,
					lapack_int *info );
void LAPACK_cgerfs( char* trans, lapack_int* n, lapack_int* nrhs,
					const lapack_complex_float* a, lapack_int* lda,
					const lapack_complex_float* af, lapack_int* ldaf,
					const lapack_int* ipiv, const lapack_complex_float* b,
					lapack_int* ldb, lapack_complex_float* x, lapack_int* ldx,
					float* ferr, float* berr, lapack_complex_float* work,
					float* rwork, lapack_int *info );
void LAPACK_zgerfs( char* trans, lapack_int* n, lapack_int* nrhs,
					const lapack_complex_double* a, lapack_int* lda,
					const lapack_complex_double* af, lapack_int* ldaf,
					const lapack_int* ipiv, const lapack_complex_double* b,
					lapack_int* ldb, lapack_complex_double* x, lapack_int* ldx,
					double* ferr, double* berr, lapack_complex_double* work,
					double* rwork, lapack_int *info );
void LAPACK_dgerfsx( char* trans, char* equed, lapack_int* n, lapack_int* nrhs,
					const double* a, lapack_int* lda, const double* af,
					lapack_int* ldaf, const lapack_int* ipiv, const double* r,
					const double* c, const double* b, lapack_int* ldb,
					double* x, lapack_int* ldx, double* rcond, double* berr,
					lapack_int* n_err_bnds, double* err_bnds_norm,
					double* err_bnds_comp, lapack_int* nparams, double* params,
					double* work, lapack_int* iwork, lapack_int *info );
void LAPACK_sgerfsx( char* trans, char* equed, lapack_int* n, lapack_int* nrhs,
					const float* a, lapack_int* lda, const float* af,
					lapack_int* ldaf, const lapack_int* ipiv, const float* r,
					const float* c, const float* b, lapack_int* ldb, float* x,
					lapack_int* ldx, float* rcond, float* berr,
					lapack_int* n_err_bnds, float* err_bnds_norm,
					float* err_bnds_comp, lapack_int* nparams, float* params,
					float* work, lapack_int* iwork, lapack_int *info );
void LAPACK_zgerfsx( char* trans, char* equed, lapack_int* n, lapack_int* nrhs,
					const lapack_complex_double* a, lapack_int* lda,
					const lapack_complex_double* af, lapack_int* ldaf,
					const lapack_int* ipiv, const double* r, const double* c,
					const lapack_complex_double* b, lapack_int* ldb,
					lapack_complex_double* x, lapack_int* ldx, double* rcond,
					double* berr, lapack_int* n_err_bnds,
					double* err_bnds_norm, double* err_bnds_comp,
					lapack_int* nparams, double* params,
					lapack_complex_double* work, double* rwork,
					lapack_int *info );
void LAPACK_cgerfsx( char* trans, char* equed, lapack_int* n, lapack_int* nrhs,
					const lapack_complex_float* a, lapack_int* lda,
					const lapack_complex_float* af, lapack_int* ldaf,
					const lapack_int* ipiv, const float* r, const float* c,
					const lapack_complex_float* b, lapack_int* ldb,
					lapack_complex_float* x, lapack_int* ldx, float* rcond,
					float* berr, lapack_int* n_err_bnds, float* err_bnds_norm,
					float* err_bnds_comp, lapack_int* nparams, float* params,
					lapack_complex_float* work, float* rwork,
					lapack_int *info );
void LAPACK_sgbrfs( char* trans, lapack_int* n, lapack_int* kl, lapack_int* ku,
					lapack_int* nrhs, const float* ab, lapack_int* ldab,
					const float* afb, lapack_int* ldafb, const lapack_int* ipiv,
					const float* b, lapack_int* ldb, float* x, lapack_int* ldx,
					float* ferr, float* berr, float* work, lapack_int* iwork,
					lapack_int *info );
void LAPACK_dgbrfs( char* trans, lapack_int* n, lapack_int* kl, lapack_int* ku,
					lapack_int* nrhs, const double* ab, lapack_int* ldab,
					const double* afb, lapack_int* ldafb,
					const lapack_int* ipiv, const double* b, lapack_int* ldb,
					double* x, lapack_int* ldx, double* ferr, double* berr,
					double* work, lapack_int* iwork, lapack_int *info );
void LAPACK_cgbrfs( char* trans, lapack_int* n, lapack_int* kl, lapack_int* ku,
					lapack_int* nrhs, const lapack_complex_float* ab,
					lapack_int* ldab, const lapack_complex_float* afb,
					lapack_int* ldafb, const lapack_int* ipiv,
					const lapack_complex_float* b, lapack_int* ldb,
					lapack_complex_float* x, lapack_int* ldx, float* ferr,
					float* berr, lapack_complex_float* work, float* rwork,
					lapack_int *info );
void LAPACK_zgbrfs( char* trans, lapack_int* n, lapack_int* kl, lapack_int* ku,
					lapack_int* nrhs, const lapack_complex_double* ab,
					lapack_int* ldab, const lapack_complex_double* afb,
					lapack_int* ldafb, const lapack_int* ipiv,
					const lapack_complex_double* b, lapack_int* ldb,
					lapack_complex_double* x, lapack_int* ldx, double* ferr,
					double* berr, lapack_complex_double* work, double* rwork,
					lapack_int *info );
void LAPACK_dgbrfsx( char* trans, char* equed, lapack_int* n, lapack_int* kl,
					lapack_int* ku, lapack_int* nrhs, const double* ab,
					lapack_int* ldab, const double* afb, lapack_int* ldafb,
					const lapack_int* ipiv, const double* r, const double* c,
					const double* b, lapack_int* ldb, double* x,
					lapack_int* ldx, double* rcond, double* berr,
					lapack_int* n_err_bnds, double* err_bnds_norm,
					double* err_bnds_comp, lapack_int* nparams, double* params,
					double* work, lapack_int* iwork, lapack_int *info );
void LAPACK_sgbrfsx( char* trans, char* equed, lapack_int* n, lapack_int* kl,
					lapack_int* ku, lapack_int* nrhs, const float* ab,
					lapack_int* ldab, const float* afb, lapack_int* ldafb,
					const lapack_int* ipiv, const float* r, const float* c,
					const float* b, lapack_int* ldb, float* x, lapack_int* ldx,
					float* rcond, float* berr, lapack_int* n_err_bnds,
					float* err_bnds_norm, float* err_bnds_comp,
					lapack_int* nparams, float* params, float* work,
					lapack_int* iwork, lapack_int *info );
void LAPACK_zgbrfsx( char* trans, char* equed, lapack_int* n, lapack_int* kl,
					lapack_int* ku, lapack_int* nrhs,
					const lapack_complex_double* ab, lapack_int* ldab,
					const lapack_complex_double* afb, lapack_int* ldafb,
					const lapack_int* ipiv, const double* r, const double* c,
					const lapack_complex_double* b, lapack_int* ldb,
					lapack_complex_double* x, lapack_int* ldx, double* rcond,
					double* berr, lapack_int* n_err_bnds,
					double* err_bnds_norm, double* err_bnds_comp,
					lapack_int* nparams, double* params,
					lapack_complex_double* work, double* rwork,
					lapack_int *info );
void LAPACK_cgbrfsx( char* trans, char* equed, lapack_int* n, lapack_int* kl,
					lapack_int* ku, lapack_int* nrhs,
					const lapack_complex_float* ab, lapack_int* ldab,
					const lapack_complex_float* afb, lapack_int* ldafb,
					const lapack_int* ipiv, const float* r, const float* c,
					const lapack_complex_float* b, lapack_int* ldb,
					lapack_complex_float* x, lapack_int* ldx, float* rcond,
					float* berr, lapack_int* n_err_bnds, float* err_bnds_norm,
					float* err_bnds_comp, lapack_int* nparams, float* params,
					lapack_complex_float* work, float* rwork,
					lapack_int *info );
void LAPACK_sgtrfs( char* trans, lapack_int* n, lapack_int* nrhs,
					const float* dl, const float* d, const float* du,
					const float* dlf, const float* df, const float* duf,
					const float* du2, const lapack_int* ipiv, const float* b,
					lapack_int* ldb, float* x, lapack_int* ldx, float* ferr,
					float* berr, float* work, lapack_int* iwork,
					lapack_int *info );
void LAPACK_dgtrfs( char* trans, lapack_int* n, lapack_int* nrhs,
					const double* dl, const double* d, const double* du,
					const double* dlf, const double* df, const double* duf,
					const double* du2, const lapack_int* ipiv, const double* b,
					lapack_int* ldb, double* x, lapack_int* ldx, double* ferr,
					double* berr, double* work, lapack_int* iwork,
					lapack_int *info );
void LAPACK_cgtrfs( char* trans, lapack_int* n, lapack_int* nrhs,
					const lapack_complex_float* dl,
					const lapack_complex_float* d,
					const lapack_complex_float* du,
					const lapack_complex_float* dlf,
					const lapack_complex_float* df,
					const lapack_complex_float* duf,
					const lapack_complex_float* du2, const lapack_int* ipiv,
					const lapack_complex_float* b, lapack_int* ldb,
					lapack_complex_float* x, lapack_int* ldx, float* ferr,
					float* berr, lapack_complex_float* work, float* rwork,
					lapack_int *info );
void LAPACK_zgtrfs( char* trans, lapack_int* n, lapack_int* nrhs,
					const lapack_complex_double* dl,
					const lapack_complex_double* d,
					const lapack_complex_double* du,
					const lapack_complex_double* dlf,
					const lapack_complex_double* df,
					const lapack_complex_double* duf,
					const lapack_complex_double* du2, const lapack_int* ipiv,
					const lapack_complex_double* b, lapack_int* ldb,
					lapack_complex_double* x, lapack_int* ldx, double* ferr,
					double* berr, lapack_complex_double* work, double* rwork,
					lapack_int *info );
void LAPACK_sporfs( char* uplo, lapack_int* n, lapack_int* nrhs, const float* a,
					lapack_int* lda, const float* af, lapack_int* ldaf,
					const float* b, lapack_int* ldb, float* x, lapack_int* ldx,
					float* ferr, float* berr, float* work, lapack_int* iwork,
					lapack_int *info );
void LAPACK_dporfs( char* uplo, lapack_int* n, lapack_int* nrhs,
					const double* a, lapack_int* lda, const double* af,
					lapack_int* ldaf, const double* b, lapack_int* ldb,
					double* x, lapack_int* ldx, double* ferr, double* berr,
					double* work, lapack_int* iwork, lapack_int *info );
void LAPACK_cporfs( char* uplo, lapack_int* n, lapack_int* nrhs,
					const lapack_complex_float* a, lapack_int* lda,
					const lapack_complex_float* af, lapack_int* ldaf,
					const lapack_complex_float* b, lapack_int* ldb,
					lapack_complex_float* x, lapack_int* ldx, float* ferr,
					float* berr, lapack_complex_float* work, float* rwork,
					lapack_int *info );
void LAPACK_zporfs( char* uplo, lapack_int* n, lapack_int* nrhs,
					const lapack_complex_double* a, lapack_int* lda,
					const lapack_complex_double* af, lapack_int* ldaf,
					const lapack_complex_double* b, lapack_int* ldb,
					lapack_complex_double* x, lapack_int* ldx, double* ferr,
					double* berr, lapack_complex_double* work, double* rwork,
					lapack_int *info );
void LAPACK_dporfsx( char* uplo, char* equed, lapack_int* n, lapack_int* nrhs,
					const double* a, lapack_int* lda, const double* af,
					lapack_int* ldaf, const double* s, const double* b,
					lapack_int* ldb, double* x, lapack_int* ldx, double* rcond,
					double* berr, lapack_int* n_err_bnds,
					double* err_bnds_norm, double* err_bnds_comp,
					lapack_int* nparams, double* params, double* work,
					lapack_int* iwork, lapack_int *info );
void LAPACK_sporfsx( char* uplo, char* equed, lapack_int* n, lapack_int* nrhs,
					const float* a, lapack_int* lda, const float* af,
					lapack_int* ldaf, const float* s, const float* b,
					lapack_int* ldb, float* x, lapack_int* ldx, float* rcond,
					float* berr, lapack_int* n_err_bnds, float* err_bnds_norm,
					float* err_bnds_comp, lapack_int* nparams, float* params,
					float* work, lapack_int* iwork, lapack_int *info );
void LAPACK_zporfsx( char* uplo, char* equed, lapack_int* n, lapack_int* nrhs,
					const lapack_complex_double* a, lapack_int* lda,
					const lapack_complex_double* af, lapack_int* ldaf,
					const double* s, const lapack_complex_double* b,
					lapack_int* ldb, lapack_complex_double* x, lapack_int* ldx,
					double* rcond, double* berr, lapack_int* n_err_bnds,
					double* err_bnds_norm, double* err_bnds_comp,
					lapack_int* nparams, double* params,
					lapack_complex_double* work, double* rwork,
					lapack_int *info );
void LAPACK_cporfsx( char* uplo, char* equed, lapack_int* n, lapack_int* nrhs,
					const lapack_complex_float* a, lapack_int* lda,
					const lapack_complex_float* af, lapack_int* ldaf,
					const float* s, const lapack_complex_float* b,
					lapack_int* ldb, lapack_complex_float* x, lapack_int* ldx,
					float* rcond, float* berr, lapack_int* n_err_bnds,
					float* err_bnds_norm, float* err_bnds_comp,
					lapack_int* nparams, float* params,
					lapack_complex_float* work, float* rwork,
					lapack_int *info );
void LAPACK_spprfs( char* uplo, lapack_int* n, lapack_int* nrhs,
					const float* ap, const float* afp, const float* b,
					lapack_int* ldb, float* x, lapack_int* ldx, float* ferr,
					float* berr, float* work, lapack_int* iwork,
					lapack_int *info );
void LAPACK_dpprfs( char* uplo, lapack_int* n, lapack_int* nrhs,
					const double* ap, const double* afp, const double* b,
					lapack_int* ldb, double* x, lapack_int* ldx, double* ferr,
					double* berr, double* work, lapack_int* iwork,
					lapack_int *info );
void LAPACK_cpprfs( char* uplo, lapack_int* n, lapack_int* nrhs,
					const lapack_complex_float* ap,
					const lapack_complex_float* afp,
					const lapack_complex_float* b, lapack_int* ldb,
					lapack_complex_float* x, lapack_int* ldx, float* ferr,
					float* berr, lapack_complex_float* work, float* rwork,
					lapack_int *info );
void LAPACK_zpprfs( char* uplo, lapack_int* n, lapack_int* nrhs,
					const lapack_complex_double* ap,
					const lapack_complex_double* afp,
					const lapack_complex_double* b, lapack_int* ldb,
					lapack_complex_double* x, lapack_int* ldx, double* ferr,
					double* berr, lapack_complex_double* work, double* rwork,
					lapack_int *info );
void LAPACK_spbrfs( char* uplo, lapack_int* n, lapack_int* kd, lapack_int* nrhs,
					const float* ab, lapack_int* ldab, const float* afb,
					lapack_int* ldafb, const float* b, lapack_int* ldb,
					float* x, lapack_int* ldx, float* ferr, float* berr,
					float* work, lapack_int* iwork, lapack_int *info );
void LAPACK_dpbrfs( char* uplo, lapack_int* n, lapack_int* kd, lapack_int* nrhs,
					const double* ab, lapack_int* ldab, const double* afb,
					lapack_int* ldafb, const double* b, lapack_int* ldb,
					double* x, lapack_int* ldx, double* ferr, double* berr,
					double* work, lapack_int* iwork, lapack_int *info );
void LAPACK_cpbrfs( char* uplo, lapack_int* n, lapack_int* kd, lapack_int* nrhs,
					const lapack_complex_float* ab, lapack_int* ldab,
					const lapack_complex_float* afb, lapack_int* ldafb,
					const lapack_complex_float* b, lapack_int* ldb,
					lapack_complex_float* x, lapack_int* ldx, float* ferr,
					float* berr, lapack_complex_float* work, float* rwork,
					lapack_int *info );
void LAPACK_zpbrfs( char* uplo, lapack_int* n, lapack_int* kd, lapack_int* nrhs,
					const lapack_complex_double* ab, lapack_int* ldab,
					const lapack_complex_double* afb, lapack_int* ldafb,
					const lapack_complex_double* b, lapack_int* ldb,
					lapack_complex_double* x, lapack_int* ldx, double* ferr,
					double* berr, lapack_complex_double* work, double* rwork,
					lapack_int *info );
void LAPACK_sptrfs( lapack_int* n, lapack_int* nrhs, const float* d,
					const float* e, const float* df, const float* ef,
					const float* b, lapack_int* ldb, float* x, lapack_int* ldx,
					float* ferr, float* berr, float* work, lapack_int *info );
void LAPACK_dptrfs( lapack_int* n, lapack_int* nrhs, const double* d,
					const double* e, const double* df, const double* ef,
					const double* b, lapack_int* ldb, double* x,
					lapack_int* ldx, double* ferr, double* berr, double* work,
					lapack_int *info );
void LAPACK_cptrfs( char* uplo, lapack_int* n, lapack_int* nrhs, const float* d,
					const lapack_complex_float* e, const float* df,
					const lapack_complex_float* ef,
					const lapack_complex_float* b, lapack_int* ldb,
					lapack_complex_float* x, lapack_int* ldx, float* ferr,
					float* berr, lapack_complex_float* work, float* rwork,
					lapack_int *info );
void LAPACK_zptrfs( char* uplo, lapack_int* n, lapack_int* nrhs,
					const double* d, const lapack_complex_double* e,
					const double* df, const lapack_complex_double* ef,
					const lapack_complex_double* b, lapack_int* ldb,
					lapack_complex_double* x, lapack_int* ldx, double* ferr,
					double* berr, lapack_complex_double* work, double* rwork,
					lapack_int *info );
void LAPACK_ssyrfs( char* uplo, lapack_int* n, lapack_int* nrhs, const float* a,
					lapack_int* lda, const float* af, lapack_int* ldaf,
					const lapack_int* ipiv, const float* b, lapack_int* ldb,
					float* x, lapack_int* ldx, float* ferr, float* berr,
					float* work, lapack_int* iwork, lapack_int *info );
void LAPACK_dsyrfs( char* uplo, lapack_int* n, lapack_int* nrhs,
					const double* a, lapack_int* lda, const double* af,
					lapack_int* ldaf, const lapack_int* ipiv, const double* b,
					lapack_int* ldb, double* x, lapack_int* ldx, double* ferr,
					double* berr, double* work, lapack_int* iwork,
					lapack_int *info );
void LAPACK_csyrfs( char* uplo, lapack_int* n, lapack_int* nrhs,
					const lapack_complex_float* a, lapack_int* lda,
					const lapack_complex_float* af, lapack_int* ldaf,
					const lapack_int* ipiv, const lapack_complex_float* b,
					lapack_int* ldb, lapack_complex_float* x, lapack_int* ldx,
					float* ferr, float* berr, lapack_complex_float* work,
					float* rwork, lapack_int *info );
void LAPACK_zsyrfs( char* uplo, lapack_int* n, lapack_int* nrhs,
					const lapack_complex_double* a, lapack_int* lda,
					const lapack_complex_double* af, lapack_int* ldaf,
					const lapack_int* ipiv, const lapack_complex_double* b,
					lapack_int* ldb, lapack_complex_double* x, lapack_int* ldx,
					double* ferr, double* berr, lapack_complex_double* work,
					double* rwork, lapack_int *info );
void LAPACK_dsyrfsx( char* uplo, char* equed, lapack_int* n, lapack_int* nrhs,
					const double* a, lapack_int* lda, const double* af,
					lapack_int* ldaf, const lapack_int* ipiv, const double* s,
					const double* b, lapack_int* ldb, double* x,
					lapack_int* ldx, double* rcond, double* berr,
					lapack_int* n_err_bnds, double* err_bnds_norm,
					double* err_bnds_comp, lapack_int* nparams, double* params,
					double* work, lapack_int* iwork, lapack_int *info );
void LAPACK_ssyrfsx( char* uplo, char* equed, lapack_int* n, lapack_int* nrhs,
					const float* a, lapack_int* lda, const float* af,
					lapack_int* ldaf, const lapack_int* ipiv, const float* s,
					const float* b, lapack_int* ldb, float* x, lapack_int* ldx,
					float* rcond, float* berr, lapack_int* n_err_bnds,
					float* err_bnds_norm, float* err_bnds_comp,
					lapack_int* nparams, float* params, float* work,
					lapack_int* iwork, lapack_int *info );
void LAPACK_zsyrfsx( char* uplo, char* equed, lapack_int* n, lapack_int* nrhs,
					const lapack_complex_double* a, lapack_int* lda,
					const lapack_complex_double* af, lapack_int* ldaf,
					const lapack_int* ipiv, const double* s,
					const lapack_complex_double* b, lapack_int* ldb,
					lapack_complex_double* x, lapack_int* ldx, double* rcond,
					double* berr, lapack_int* n_err_bnds,
					double* err_bnds_norm, double* err_bnds_comp,
					lapack_int* nparams, double* params,
					lapack_complex_double* work, double* rwork,
					lapack_int *info );
void LAPACK_csyrfsx( char* uplo, char* equed, lapack_int* n, lapack_int* nrhs,
					const lapack_complex_float* a, lapack_int* lda,
					const lapack_complex_float* af, lapack_int* ldaf,
					const lapack_int* ipiv, const float* s,
					const lapack_complex_float* b, lapack_int* ldb,
					lapack_complex_float* x, lapack_int* ldx, float* rcond,
					float* berr, lapack_int* n_err_bnds, float* err_bnds_norm,
					float* err_bnds_comp, lapack_int* nparams, float* params,
					lapack_complex_float* work, float* rwork,
					lapack_int *info );
void LAPACK_cherfs( char* uplo, lapack_int* n, lapack_int* nrhs,
					const lapack_complex_float* a, lapack_int* lda,
					const lapack_complex_float* af, lapack_int* ldaf,
					const lapack_int* ipiv, const lapack_complex_float* b,
					lapack_int* ldb, lapack_complex_float* x, lapack_int* ldx,
					float* ferr, float* berr, lapack_complex_float* work,
					float* rwork, lapack_int *info );
void LAPACK_zherfs( char* uplo, lapack_int* n, lapack_int* nrhs,
					const lapack_complex_double* a, lapack_int* lda,
					const lapack_complex_double* af, lapack_int* ldaf,
					const lapack_int* ipiv, const lapack_complex_double* b,
					lapack_int* ldb, lapack_complex_double* x, lapack_int* ldx,
					double* ferr, double* berr, lapack_complex_double* work,
					double* rwork, lapack_int *info );
void LAPACK_zherfsx( char* uplo, char* equed, lapack_int* n, lapack_int* nrhs,
					const lapack_complex_double* a, lapack_int* lda,
					const lapack_complex_double* af, lapack_int* ldaf,
					const lapack_int* ipiv, const double* s,
					const lapack_complex_double* b, lapack_int* ldb,
					lapack_complex_double* x, lapack_int* ldx, double* rcond,
					double* berr, lapack_int* n_err_bnds,
					double* err_bnds_norm, double* err_bnds_comp,
					lapack_int* nparams, double* params,
					lapack_complex_double* work, double* rwork,
					lapack_int *info );
void LAPACK_cherfsx( char* uplo, char* equed, lapack_int* n, lapack_int* nrhs,
					const lapack_complex_float* a, lapack_int* lda,
					const lapack_complex_float* af, lapack_int* ldaf,
					const lapack_int* ipiv, const float* s,
					const lapack_complex_float* b, lapack_int* ldb,
					lapack_complex_float* x, lapack_int* ldx, float* rcond,
					float* berr, lapack_int* n_err_bnds, float* err_bnds_norm,
					float* err_bnds_comp, lapack_int* nparams, float* params,
					lapack_complex_float* work, float* rwork,
					lapack_int *info );
void LAPACK_ssprfs( char* uplo, lapack_int* n, lapack_int* nrhs,
					const float* ap, const float* afp, const lapack_int* ipiv,
					const float* b, lapack_int* ldb, float* x, lapack_int* ldx,
					float* ferr, float* berr, float* work, lapack_int* iwork,
					lapack_int *info );
void LAPACK_dsprfs( char* uplo, lapack_int* n, lapack_int* nrhs,
					const double* ap, const double* afp, const lapack_int* ipiv,
					const double* b, lapack_int* ldb, double* x,
					lapack_int* ldx, double* ferr, double* berr, double* work,
					lapack_int* iwork, lapack_int *info );
void LAPACK_csprfs( char* uplo, lapack_int* n, lapack_int* nrhs,
					const lapack_complex_float* ap,
					const lapack_complex_float* afp, const lapack_int* ipiv,
					const lapack_complex_float* b, lapack_int* ldb,
					lapack_complex_float* x, lapack_int* ldx, float* ferr,
					float* berr, lapack_complex_float* work, float* rwork,
					lapack_int *info );
void LAPACK_zsprfs( char* uplo, lapack_int* n, lapack_int* nrhs,
					const lapack_complex_double* ap,
					const lapack_complex_double* afp, const lapack_int* ipiv,
					const lapack_complex_double* b, lapack_int* ldb,
					lapack_complex_double* x, lapack_int* ldx, double* ferr,
					double* berr, lapack_complex_double* work, double* rwork,
					lapack_int *info );
void LAPACK_chprfs( char* uplo, lapack_int* n, lapack_int* nrhs,
					const lapack_complex_float* ap,
					const lapack_complex_float* afp, const lapack_int* ipiv,
					const lapack_complex_float* b, lapack_int* ldb,
					lapack_complex_float* x, lapack_int* ldx, float* ferr,
					float* berr, lapack_complex_float* work, float* rwork,
					lapack_int *info );
void LAPACK_zhprfs( char* uplo, lapack_int* n, lapack_int* nrhs,
					const lapack_complex_double* ap,
					const lapack_complex_double* afp, const lapack_int* ipiv,
					const lapack_complex_double* b, lapack_int* ldb,
					lapack_complex_double* x, lapack_int* ldx, double* ferr,
					double* berr, lapack_complex_double* work, double* rwork,
					lapack_int *info );
void LAPACK_strrfs( char* uplo, char* trans, char* diag, lapack_int* n,
					lapack_int* nrhs, const float* a, lapack_int* lda,
					const float* b, lapack_int* ldb, const float* x,
					lapack_int* ldx, float* ferr, float* berr, float* work,
					lapack_int* iwork, lapack_int *info );
void LAPACK_dtrrfs( char* uplo, char* trans, char* diag, lapack_int* n,
					lapack_int* nrhs, const double* a, lapack_int* lda,
					const double* b, lapack_int* ldb, const double* x,
					lapack_int* ldx, double* ferr, double* berr, double* work,
					lapack_int* iwork, lapack_int *info );
void LAPACK_ctrrfs( char* uplo, char* trans, char* diag, lapack_int* n,
					lapack_int* nrhs, const lapack_complex_float* a,
					lapack_int* lda, const lapack_complex_float* b,
					lapack_int* ldb, const lapack_complex_float* x,
					lapack_int* ldx, float* ferr, float* berr,
					lapack_complex_float* work, float* rwork,
					lapack_int *info );
void LAPACK_ztrrfs( char* uplo, char* trans, char* diag, lapack_int* n,
					lapack_int* nrhs, const lapack_complex_double* a,
					lapack_int* lda, const lapack_complex_double* b,
					lapack_int* ldb, const lapack_complex_double* x,
					lapack_int* ldx, double* ferr, double* berr,
					lapack_complex_double* work, double* rwork,
					lapack_int *info );
void LAPACK_stprfs( char* uplo, char* trans, char* diag, lapack_int* n,
					lapack_int* nrhs, const float* ap, const float* b,
					lapack_int* ldb, const float* x, lapack_int* ldx,
					float* ferr, float* berr, float* work, lapack_int* iwork,
					lapack_int *info );
void LAPACK_dtprfs( char* uplo, char* trans, char* diag, lapack_int* n,
					lapack_int* nrhs, const double* ap, const double* b,
					lapack_int* ldb, const double* x, lapack_int* ldx,
					double* ferr, double* berr, double* work, lapack_int* iwork,
					lapack_int *info );
void LAPACK_ctprfs( char* uplo, char* trans, char* diag, lapack_int* n,
					lapack_int* nrhs, const lapack_complex_float* ap,
					const lapack_complex_float* b, lapack_int* ldb,
					const lapack_complex_float* x, lapack_int* ldx, float* ferr,
					float* berr, lapack_complex_float* work, float* rwork,
					lapack_int *info );
void LAPACK_ztprfs( char* uplo, char* trans, char* diag, lapack_int* n,
					lapack_int* nrhs, const lapack_complex_double* ap,
					const lapack_complex_double* b, lapack_int* ldb,
					const lapack_complex_double* x, lapack_int* ldx,
					double* ferr, double* berr, lapack_complex_double* work,
					double* rwork, lapack_int *info );
void LAPACK_stbrfs( char* uplo, char* trans, char* diag, lapack_int* n,
					lapack_int* kd, lapack_int* nrhs, const float* ab,
					lapack_int* ldab, const float* b, lapack_int* ldb,
					const float* x, lapack_int* ldx, float* ferr, float* berr,
					float* work, lapack_int* iwork, lapack_int *info );
void LAPACK_dtbrfs( char* uplo, char* trans, char* diag, lapack_int* n,
					lapack_int* kd, lapack_int* nrhs, const double* ab,
					lapack_int* ldab, const double* b, lapack_int* ldb,
					const double* x, lapack_int* ldx, double* ferr,
					double* berr, double* work, lapack_int* iwork,
					lapack_int *info );
void LAPACK_ctbrfs( char* uplo, char* trans, char* diag, lapack_int* n,
					lapack_int* kd, lapack_int* nrhs,
					const lapack_complex_float* ab, lapack_int* ldab,
					const lapack_complex_float* b, lapack_int* ldb,
					const lapack_complex_float* x, lapack_int* ldx, float* ferr,
					float* berr, lapack_complex_float* work, float* rwork,
					lapack_int *info );
void LAPACK_ztbrfs( char* uplo, char* trans, char* diag, lapack_int* n,
					lapack_int* kd, lapack_int* nrhs,
					const lapack_complex_double* ab, lapack_int* ldab,
					const lapack_complex_double* b, lapack_int* ldb,
					const lapack_complex_double* x, lapack_int* ldx,
					double* ferr, double* berr, lapack_complex_double* work,
					double* rwork, lapack_int *info );
void LAPACK_sgetri( lapack_int* n, float* a, lapack_int* lda,
					const lapack_int* ipiv, float* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_dgetri( lapack_int* n, double* a, lapack_int* lda,
					const lapack_int* ipiv, double* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_cgetri( lapack_int* n, lapack_complex_float* a, lapack_int* lda,
					const lapack_int* ipiv, lapack_complex_float* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_zgetri( lapack_int* n, lapack_complex_double* a, lapack_int* lda,
					const lapack_int* ipiv, lapack_complex_double* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_spotri( char* uplo, lapack_int* n, float* a, lapack_int* lda,
					lapack_int *info );
void LAPACK_dpotri( char* uplo, lapack_int* n, double* a, lapack_int* lda,
					lapack_int *info );
void LAPACK_cpotri( char* uplo, lapack_int* n, lapack_complex_float* a,
					lapack_int* lda, lapack_int *info );
void LAPACK_zpotri( char* uplo, lapack_int* n, lapack_complex_double* a,
					lapack_int* lda, lapack_int *info );
void LAPACK_dpftri( char* transr, char* uplo, lapack_int* n, double* a,
					lapack_int *info );
void LAPACK_spftri( char* transr, char* uplo, lapack_int* n, float* a,
					lapack_int *info );
void LAPACK_zpftri( char* transr, char* uplo, lapack_int* n,
					lapack_complex_double* a, lapack_int *info );
void LAPACK_cpftri( char* transr, char* uplo, lapack_int* n,
					lapack_complex_float* a, lapack_int *info );
void LAPACK_spptri( char* uplo, lapack_int* n, float* ap, lapack_int *info );
void LAPACK_dpptri( char* uplo, lapack_int* n, double* ap, lapack_int *info );
void LAPACK_cpptri( char* uplo, lapack_int* n, lapack_complex_float* ap,
					lapack_int *info );
void LAPACK_zpptri( char* uplo, lapack_int* n, lapack_complex_double* ap,
					lapack_int *info );
void LAPACK_ssytri( char* uplo, lapack_int* n, float* a, lapack_int* lda,
					const lapack_int* ipiv, float* work, lapack_int *info );
void LAPACK_dsytri( char* uplo, lapack_int* n, double* a, lapack_int* lda,
					const lapack_int* ipiv, double* work, lapack_int *info );
void LAPACK_csytri( char* uplo, lapack_int* n, lapack_complex_float* a,
					lapack_int* lda, const lapack_int* ipiv,
					lapack_complex_float* work, lapack_int *info );
void LAPACK_zsytri( char* uplo, lapack_int* n, lapack_complex_double* a,
					lapack_int* lda, const lapack_int* ipiv,
					lapack_complex_double* work, lapack_int *info );
void LAPACK_chetri( char* uplo, lapack_int* n, lapack_complex_float* a,
					lapack_int* lda, const lapack_int* ipiv,
					lapack_complex_float* work, lapack_int *info );
void LAPACK_zhetri( char* uplo, lapack_int* n, lapack_complex_double* a,
					lapack_int* lda, const lapack_int* ipiv,
					lapack_complex_double* work, lapack_int *info );
void LAPACK_ssptri( char* uplo, lapack_int* n, float* ap,
					const lapack_int* ipiv, float* work, lapack_int *info );
void LAPACK_dsptri( char* uplo, lapack_int* n, double* ap,
					const lapack_int* ipiv, double* work, lapack_int *info );
void LAPACK_csptri( char* uplo, lapack_int* n, lapack_complex_float* ap,
					const lapack_int* ipiv, lapack_complex_float* work,
					lapack_int *info );
void LAPACK_zsptri( char* uplo, lapack_int* n, lapack_complex_double* ap,
					const lapack_int* ipiv, lapack_complex_double* work,
					lapack_int *info );
void LAPACK_chptri( char* uplo, lapack_int* n, lapack_complex_float* ap,
					const lapack_int* ipiv, lapack_complex_float* work,
					lapack_int *info );
void LAPACK_zhptri( char* uplo, lapack_int* n, lapack_complex_double* ap,
					const lapack_int* ipiv, lapack_complex_double* work,
					lapack_int *info );
void LAPACK_strtri( char* uplo, char* diag, lapack_int* n, float* a,
					lapack_int* lda, lapack_int *info );
void LAPACK_dtrtri( char* uplo, char* diag, lapack_int* n, double* a,
					lapack_int* lda, lapack_int *info );
void LAPACK_ctrtri( char* uplo, char* diag, lapack_int* n,
					lapack_complex_float* a, lapack_int* lda,
					lapack_int *info );
void LAPACK_ztrtri( char* uplo, char* diag, lapack_int* n,
					lapack_complex_double* a, lapack_int* lda,
					lapack_int *info );
void LAPACK_dtftri( char* transr, char* uplo, char* diag, lapack_int* n,
					double* a, lapack_int *info );
void LAPACK_stftri( char* transr, char* uplo, char* diag, lapack_int* n,
					float* a, lapack_int *info );
void LAPACK_ztftri( char* transr, char* uplo, char* diag, lapack_int* n,
					lapack_complex_double* a, lapack_int *info );
void LAPACK_ctftri( char* transr, char* uplo, char* diag, lapack_int* n,
					lapack_complex_float* a, lapack_int *info );
void LAPACK_stptri( char* uplo, char* diag, lapack_int* n, float* ap,
					lapack_int *info );
void LAPACK_dtptri( char* uplo, char* diag, lapack_int* n, double* ap,
					lapack_int *info );
void LAPACK_ctptri( char* uplo, char* diag, lapack_int* n,
					lapack_complex_float* ap, lapack_int *info );
void LAPACK_ztptri( char* uplo, char* diag, lapack_int* n,
					lapack_complex_double* ap, lapack_int *info );
void LAPACK_sgeequ( lapack_int* m, lapack_int* n, const float* a,
					lapack_int* lda, float* r, float* c, float* rowcnd,
					float* colcnd, float* amax, lapack_int *info );
void LAPACK_dgeequ( lapack_int* m, lapack_int* n, const double* a,
					lapack_int* lda, double* r, double* c, double* rowcnd,
					double* colcnd, double* amax, lapack_int *info );
void LAPACK_cgeequ( lapack_int* m, lapack_int* n, const lapack_complex_float* a,
					lapack_int* lda, float* r, float* c, float* rowcnd,
					float* colcnd, float* amax, lapack_int *info );
void LAPACK_zgeequ( lapack_int* m, lapack_int* n,
					const lapack_complex_double* a, lapack_int* lda, double* r,
					double* c, double* rowcnd, double* colcnd, double* amax,
					lapack_int *info );
void LAPACK_dgeequb( lapack_int* m, lapack_int* n, const double* a,
					lapack_int* lda, double* r, double* c, double* rowcnd,
					double* colcnd, double* amax, lapack_int *info );
void LAPACK_sgeequb( lapack_int* m, lapack_int* n, const float* a,
					lapack_int* lda, float* r, float* c, float* rowcnd,
					float* colcnd, float* amax, lapack_int *info );
void LAPACK_zgeequb( lapack_int* m, lapack_int* n,
					const lapack_complex_double* a, lapack_int* lda, double* r,
					double* c, double* rowcnd, double* colcnd, double* amax,
					lapack_int *info );
void LAPACK_cgeequb( lapack_int* m, lapack_int* n,
					const lapack_complex_float* a, lapack_int* lda, float* r,
					float* c, float* rowcnd, float* colcnd, float* amax,
					lapack_int *info );
void LAPACK_sgbequ( lapack_int* m, lapack_int* n, lapack_int* kl,
					lapack_int* ku, const float* ab, lapack_int* ldab, float* r,
					float* c, float* rowcnd, float* colcnd, float* amax,
					lapack_int *info );
void LAPACK_dgbequ( lapack_int* m, lapack_int* n, lapack_int* kl,
					lapack_int* ku, const double* ab, lapack_int* ldab,
					double* r, double* c, double* rowcnd, double* colcnd,
					double* amax, lapack_int *info );
void LAPACK_cgbequ( lapack_int* m, lapack_int* n, lapack_int* kl,
					lapack_int* ku, const lapack_complex_float* ab,
					lapack_int* ldab, float* r, float* c, float* rowcnd,
					float* colcnd, float* amax, lapack_int *info );
void LAPACK_zgbequ( lapack_int* m, lapack_int* n, lapack_int* kl,
					lapack_int* ku, const lapack_complex_double* ab,
					lapack_int* ldab, double* r, double* c, double* rowcnd,
					double* colcnd, double* amax, lapack_int *info );
void LAPACK_dgbequb( lapack_int* m, lapack_int* n, lapack_int* kl,
					lapack_int* ku, const double* ab, lapack_int* ldab,
					double* r, double* c, double* rowcnd, double* colcnd,
					double* amax, lapack_int *info );
void LAPACK_sgbequb( lapack_int* m, lapack_int* n, lapack_int* kl,
					lapack_int* ku, const float* ab, lapack_int* ldab,
					float* r, float* c, float* rowcnd, float* colcnd,
					float* amax, lapack_int *info );
void LAPACK_zgbequb( lapack_int* m, lapack_int* n, lapack_int* kl,
					lapack_int* ku, const lapack_complex_double* ab,
					lapack_int* ldab, double* r, double* c, double* rowcnd,
					double* colcnd, double* amax, lapack_int *info );
void LAPACK_cgbequb( lapack_int* m, lapack_int* n, lapack_int* kl,
					lapack_int* ku, const lapack_complex_float* ab,
					lapack_int* ldab, float* r, float* c, float* rowcnd,
					float* colcnd, float* amax, lapack_int *info );
void LAPACK_spoequ( lapack_int* n, const float* a, lapack_int* lda, float* s,
					float* scond, float* amax, lapack_int *info );
void LAPACK_dpoequ( lapack_int* n, const double* a, lapack_int* lda, double* s,
					double* scond, double* amax, lapack_int *info );
void LAPACK_cpoequ( lapack_int* n, const lapack_complex_float* a,
					lapack_int* lda, float* s, float* scond, float* amax,
					lapack_int *info );
void LAPACK_zpoequ( lapack_int* n, const lapack_complex_double* a,
					lapack_int* lda, double* s, double* scond, double* amax,
					lapack_int *info );
void LAPACK_dpoequb( lapack_int* n, const double* a, lapack_int* lda, double* s,
					double* scond, double* amax, lapack_int *info );
void LAPACK_spoequb( lapack_int* n, const float* a, lapack_int* lda, float* s,
					float* scond, float* amax, lapack_int *info );
void LAPACK_zpoequb( lapack_int* n, const lapack_complex_double* a,
					lapack_int* lda, double* s, double* scond, double* amax,
					lapack_int *info );
void LAPACK_cpoequb( lapack_int* n, const lapack_complex_float* a,
					lapack_int* lda, float* s, float* scond, float* amax,
					lapack_int *info );
void LAPACK_sppequ( char* uplo, lapack_int* n, const float* ap, float* s,
					float* scond, float* amax, lapack_int *info );
void LAPACK_dppequ( char* uplo, lapack_int* n, const double* ap, double* s,
					double* scond, double* amax, lapack_int *info );
void LAPACK_cppequ( char* uplo, lapack_int* n, const lapack_complex_float* ap,
					float* s, float* scond, float* amax, lapack_int *info );
void LAPACK_zppequ( char* uplo, lapack_int* n, const lapack_complex_double* ap,
					double* s, double* scond, double* amax, lapack_int *info );
void LAPACK_spbequ( char* uplo, lapack_int* n, lapack_int* kd, const float* ab,
					lapack_int* ldab, float* s, float* scond, float* amax,
					lapack_int *info );
void LAPACK_dpbequ( char* uplo, lapack_int* n, lapack_int* kd, const double* ab,
					lapack_int* ldab, double* s, double* scond, double* amax,
					lapack_int *info );
void LAPACK_cpbequ( char* uplo, lapack_int* n, lapack_int* kd,
					const lapack_complex_float* ab, lapack_int* ldab, float* s,
					float* scond, float* amax, lapack_int *info );
void LAPACK_zpbequ( char* uplo, lapack_int* n, lapack_int* kd,
					const lapack_complex_double* ab, lapack_int* ldab,
					double* s, double* scond, double* amax, lapack_int *info );
void LAPACK_dsyequb( char* uplo, lapack_int* n, const double* a,
					lapack_int* lda, double* s, double* scond, double* amax,
					double* work, lapack_int *info );
void LAPACK_ssyequb( char* uplo, lapack_int* n, const float* a, lapack_int* lda,
					float* s, float* scond, float* amax, float* work,
					lapack_int *info );
void LAPACK_zsyequb( char* uplo, lapack_int* n, const lapack_complex_double* a,
					lapack_int* lda, double* s, double* scond, double* amax,
					lapack_complex_double* work, lapack_int *info );
void LAPACK_csyequb( char* uplo, lapack_int* n, const lapack_complex_float* a,
					lapack_int* lda, float* s, float* scond, float* amax,
					lapack_complex_float* work, lapack_int *info );
void LAPACK_zheequb( char* uplo, lapack_int* n, const lapack_complex_double* a,
					lapack_int* lda, double* s, double* scond, double* amax,
					lapack_complex_double* work, lapack_int *info );
void LAPACK_cheequb( char* uplo, lapack_int* n, const lapack_complex_float* a,
					lapack_int* lda, float* s, float* scond, float* amax,
					lapack_complex_float* work, lapack_int *info );
void LAPACK_sgesv( lapack_int* n, lapack_int* nrhs, float* a, lapack_int* lda,
				lapack_int* ipiv, float* b, lapack_int* ldb,
				lapack_int *info );
void LAPACK_dgesv( lapack_int* n, lapack_int* nrhs, double* a, lapack_int* lda,
				lapack_int* ipiv, double* b, lapack_int* ldb,
				lapack_int *info );
void LAPACK_cgesv( lapack_int* n, lapack_int* nrhs, lapack_complex_float* a,
				lapack_int* lda, lapack_int* ipiv, lapack_complex_float* b,
				lapack_int* ldb, lapack_int *info );
void LAPACK_zgesv( lapack_int* n, lapack_int* nrhs, lapack_complex_double* a,
				lapack_int* lda, lapack_int* ipiv, lapack_complex_double* b,
				lapack_int* ldb, lapack_int *info );
void LAPACK_dsgesv( lapack_int* n, lapack_int* nrhs, double* a, lapack_int* lda,
					lapack_int* ipiv, double* b, lapack_int* ldb, double* x,
					lapack_int* ldx, double* work, float* swork,
					lapack_int* iter, lapack_int *info );
void LAPACK_zcgesv( lapack_int* n, lapack_int* nrhs, lapack_complex_double* a,
					lapack_int* lda, lapack_int* ipiv, lapack_complex_double* b,
					lapack_int* ldb, lapack_complex_double* x, lapack_int* ldx,
					lapack_complex_double* work, lapack_complex_float* swork,
					double* rwork, lapack_int* iter, lapack_int *info );
void LAPACK_sgesvx( char* fact, char* trans, lapack_int* n, lapack_int* nrhs,
					float* a, lapack_int* lda, float* af, lapack_int* ldaf,
					lapack_int* ipiv, char* equed, float* r, float* c, float* b,
					lapack_int* ldb, float* x, lapack_int* ldx, float* rcond,
					float* ferr, float* berr, float* work, lapack_int* iwork,
					lapack_int *info );
void LAPACK_dgesvx( char* fact, char* trans, lapack_int* n, lapack_int* nrhs,
					double* a, lapack_int* lda, double* af, lapack_int* ldaf,
					lapack_int* ipiv, char* equed, double* r, double* c,
					double* b, lapack_int* ldb, double* x, lapack_int* ldx,
					double* rcond, double* ferr, double* berr, double* work,
					lapack_int* iwork, lapack_int *info );
void LAPACK_cgesvx( char* fact, char* trans, lapack_int* n, lapack_int* nrhs,
					lapack_complex_float* a, lapack_int* lda,
					lapack_complex_float* af, lapack_int* ldaf,
					lapack_int* ipiv, char* equed, float* r, float* c,
					lapack_complex_float* b, lapack_int* ldb,
					lapack_complex_float* x, lapack_int* ldx, float* rcond,
					float* ferr, float* berr, lapack_complex_float* work,
					float* rwork, lapack_int *info );
void LAPACK_zgesvx( char* fact, char* trans, lapack_int* n, lapack_int* nrhs,
					lapack_complex_double* a, lapack_int* lda,
					lapack_complex_double* af, lapack_int* ldaf,
					lapack_int* ipiv, char* equed, double* r, double* c,
					lapack_complex_double* b, lapack_int* ldb,
					lapack_complex_double* x, lapack_int* ldx, double* rcond,
					double* ferr, double* berr, lapack_complex_double* work,
					double* rwork, lapack_int *info );
void LAPACK_dgesvxx( char* fact, char* trans, lapack_int* n, lapack_int* nrhs,
					double* a, lapack_int* lda, double* af, lapack_int* ldaf,
					lapack_int* ipiv, char* equed, double* r, double* c,
					double* b, lapack_int* ldb, double* x, lapack_int* ldx,
					double* rcond, double* rpvgrw, double* berr,
					lapack_int* n_err_bnds, double* err_bnds_norm,
					double* err_bnds_comp, lapack_int* nparams, double* params,
					double* work, lapack_int* iwork, lapack_int *info );
void LAPACK_sgesvxx( char* fact, char* trans, lapack_int* n, lapack_int* nrhs,
					float* a, lapack_int* lda, float* af, lapack_int* ldaf,
					lapack_int* ipiv, char* equed, float* r, float* c,
					float* b, lapack_int* ldb, float* x, lapack_int* ldx,
					float* rcond, float* rpvgrw, float* berr,
					lapack_int* n_err_bnds, float* err_bnds_norm,
					float* err_bnds_comp, lapack_int* nparams, float* params,
					float* work, lapack_int* iwork, lapack_int *info );
void LAPACK_zgesvxx( char* fact, char* trans, lapack_int* n, lapack_int* nrhs,
					lapack_complex_double* a, lapack_int* lda,
					lapack_complex_double* af, lapack_int* ldaf,
					lapack_int* ipiv, char* equed, double* r, double* c,
					lapack_complex_double* b, lapack_int* ldb,
					lapack_complex_double* x, lapack_int* ldx, double* rcond,
					double* rpvgrw, double* berr, lapack_int* n_err_bnds,
					double* err_bnds_norm, double* err_bnds_comp,
					lapack_int* nparams, double* params,
					lapack_complex_double* work, double* rwork,
					lapack_int *info );
void LAPACK_cgesvxx( char* fact, char* trans, lapack_int* n, lapack_int* nrhs,
					lapack_complex_float* a, lapack_int* lda,
					lapack_complex_float* af, lapack_int* ldaf,
					lapack_int* ipiv, char* equed, float* r, float* c,
					lapack_complex_float* b, lapack_int* ldb,
					lapack_complex_float* x, lapack_int* ldx, float* rcond,
					float* rpvgrw, float* berr, lapack_int* n_err_bnds,
					float* err_bnds_norm, float* err_bnds_comp,
					lapack_int* nparams, float* params,
					lapack_complex_float* work, float* rwork,
					lapack_int *info );
void LAPACK_sgbsv( lapack_int* n, lapack_int* kl, lapack_int* ku,
				lapack_int* nrhs, float* ab, lapack_int* ldab,
				lapack_int* ipiv, float* b, lapack_int* ldb,
				lapack_int *info );
void LAPACK_dgbsv( lapack_int* n, lapack_int* kl, lapack_int* ku,
				lapack_int* nrhs, double* ab, lapack_int* ldab,
				lapack_int* ipiv, double* b, lapack_int* ldb,
				lapack_int *info );
void LAPACK_cgbsv( lapack_int* n, lapack_int* kl, lapack_int* ku,
				lapack_int* nrhs, lapack_complex_float* ab, lapack_int* ldab,
				lapack_int* ipiv, lapack_complex_float* b, lapack_int* ldb,
				lapack_int *info );
void LAPACK_zgbsv( lapack_int* n, lapack_int* kl, lapack_int* ku,
				lapack_int* nrhs, lapack_complex_double* ab,
				lapack_int* ldab, lapack_int* ipiv, lapack_complex_double* b,
				lapack_int* ldb, lapack_int *info );
void LAPACK_sgbsvx( char* fact, char* trans, lapack_int* n, lapack_int* kl,
					lapack_int* ku, lapack_int* nrhs, float* ab,
					lapack_int* ldab, float* afb, lapack_int* ldafb,
					lapack_int* ipiv, char* equed, float* r, float* c, float* b,
					lapack_int* ldb, float* x, lapack_int* ldx, float* rcond,
					float* ferr, float* berr, float* work, lapack_int* iwork,
					lapack_int *info );
void LAPACK_dgbsvx( char* fact, char* trans, lapack_int* n, lapack_int* kl,
					lapack_int* ku, lapack_int* nrhs, double* ab,
					lapack_int* ldab, double* afb, lapack_int* ldafb,
					lapack_int* ipiv, char* equed, double* r, double* c,
					double* b, lapack_int* ldb, double* x, lapack_int* ldx,
					double* rcond, double* ferr, double* berr, double* work,
					lapack_int* iwork, lapack_int *info );
void LAPACK_cgbsvx( char* fact, char* trans, lapack_int* n, lapack_int* kl,
					lapack_int* ku, lapack_int* nrhs, lapack_complex_float* ab,
					lapack_int* ldab, lapack_complex_float* afb,
					lapack_int* ldafb, lapack_int* ipiv, char* equed, float* r,
					float* c, lapack_complex_float* b, lapack_int* ldb,
					lapack_complex_float* x, lapack_int* ldx, float* rcond,
					float* ferr, float* berr, lapack_complex_float* work,
					float* rwork, lapack_int *info );
void LAPACK_zgbsvx( char* fact, char* trans, lapack_int* n, lapack_int* kl,
					lapack_int* ku, lapack_int* nrhs, lapack_complex_double* ab,
					lapack_int* ldab, lapack_complex_double* afb,
					lapack_int* ldafb, lapack_int* ipiv, char* equed, double* r,
					double* c, lapack_complex_double* b, lapack_int* ldb,
					lapack_complex_double* x, lapack_int* ldx, double* rcond,
					double* ferr, double* berr, lapack_complex_double* work,
					double* rwork, lapack_int *info );
void LAPACK_dgbsvxx( char* fact, char* trans, lapack_int* n, lapack_int* kl,
					lapack_int* ku, lapack_int* nrhs, double* ab,
					lapack_int* ldab, double* afb, lapack_int* ldafb,
					lapack_int* ipiv, char* equed, double* r, double* c,
					double* b, lapack_int* ldb, double* x, lapack_int* ldx,
					double* rcond, double* rpvgrw, double* berr,
					lapack_int* n_err_bnds, double* err_bnds_norm,
					double* err_bnds_comp, lapack_int* nparams, double* params,
					double* work, lapack_int* iwork, lapack_int *info );
void LAPACK_sgbsvxx( char* fact, char* trans, lapack_int* n, lapack_int* kl,
					lapack_int* ku, lapack_int* nrhs, float* ab,
					lapack_int* ldab, float* afb, lapack_int* ldafb,
					lapack_int* ipiv, char* equed, float* r, float* c,
					float* b, lapack_int* ldb, float* x, lapack_int* ldx,
					float* rcond, float* rpvgrw, float* berr,
					lapack_int* n_err_bnds, float* err_bnds_norm,
					float* err_bnds_comp, lapack_int* nparams, float* params,
					float* work, lapack_int* iwork, lapack_int *info );
void LAPACK_zgbsvxx( char* fact, char* trans, lapack_int* n, lapack_int* kl,
					lapack_int* ku, lapack_int* nrhs,
					lapack_complex_double* ab, lapack_int* ldab,
					lapack_complex_double* afb, lapack_int* ldafb,
					lapack_int* ipiv, char* equed, double* r, double* c,
					lapack_complex_double* b, lapack_int* ldb,
					lapack_complex_double* x, lapack_int* ldx, double* rcond,
					double* rpvgrw, double* berr, lapack_int* n_err_bnds,
					double* err_bnds_norm, double* err_bnds_comp,
					lapack_int* nparams, double* params,
					lapack_complex_double* work, double* rwork,
					lapack_int *info );
void LAPACK_cgbsvxx( char* fact, char* trans, lapack_int* n, lapack_int* kl,
					lapack_int* ku, lapack_int* nrhs, lapack_complex_float* ab,
					lapack_int* ldab, lapack_complex_float* afb,
					lapack_int* ldafb, lapack_int* ipiv, char* equed, float* r,
					float* c, lapack_complex_float* b, lapack_int* ldb,
					lapack_complex_float* x, lapack_int* ldx, float* rcond,
					float* rpvgrw, float* berr, lapack_int* n_err_bnds,
					float* err_bnds_norm, float* err_bnds_comp,
					lapack_int* nparams, float* params,
					lapack_complex_float* work, float* rwork,
					lapack_int *info );
void LAPACK_sgtsv( lapack_int* n, lapack_int* nrhs, float* dl, float* d,
				float* du, float* b, lapack_int* ldb, lapack_int *info );
void LAPACK_dgtsv( lapack_int* n, lapack_int* nrhs, double* dl, double* d,
				double* du, double* b, lapack_int* ldb, lapack_int *info );
void LAPACK_cgtsv( lapack_int* n, lapack_int* nrhs, lapack_complex_float* dl,
				lapack_complex_float* d, lapack_complex_float* du,
				lapack_complex_float* b, lapack_int* ldb, lapack_int *info );
void LAPACK_zgtsv( lapack_int* n, lapack_int* nrhs, lapack_complex_double* dl,
				lapack_complex_double* d, lapack_complex_double* du,
				lapack_complex_double* b, lapack_int* ldb,
				lapack_int *info );
void LAPACK_sgtsvx( char* fact, char* trans, lapack_int* n, lapack_int* nrhs,
					const float* dl, const float* d, const float* du,
					float* dlf, float* df, float* duf, float* du2,
					lapack_int* ipiv, const float* b, lapack_int* ldb, float* x,
					lapack_int* ldx, float* rcond, float* ferr, float* berr,
					float* work, lapack_int* iwork, lapack_int *info );
void LAPACK_dgtsvx( char* fact, char* trans, lapack_int* n, lapack_int* nrhs,
					const double* dl, const double* d, const double* du,
					double* dlf, double* df, double* duf, double* du2,
					lapack_int* ipiv, const double* b, lapack_int* ldb,
					double* x, lapack_int* ldx, double* rcond, double* ferr,
					double* berr, double* work, lapack_int* iwork,
					lapack_int *info );
void LAPACK_cgtsvx( char* fact, char* trans, lapack_int* n, lapack_int* nrhs,
					const lapack_complex_float* dl,
					const lapack_complex_float* d,
					const lapack_complex_float* du, lapack_complex_float* dlf,
					lapack_complex_float* df, lapack_complex_float* duf,
					lapack_complex_float* du2, lapack_int* ipiv,
					const lapack_complex_float* b, lapack_int* ldb,
					lapack_complex_float* x, lapack_int* ldx, float* rcond,
					float* ferr, float* berr, lapack_complex_float* work,
					float* rwork, lapack_int *info );
void LAPACK_zgtsvx( char* fact, char* trans, lapack_int* n, lapack_int* nrhs,
					const lapack_complex_double* dl,
					const lapack_complex_double* d,
					const lapack_complex_double* du, lapack_complex_double* dlf,
					lapack_complex_double* df, lapack_complex_double* duf,
					lapack_complex_double* du2, lapack_int* ipiv,
					const lapack_complex_double* b, lapack_int* ldb,
					lapack_complex_double* x, lapack_int* ldx, double* rcond,
					double* ferr, double* berr, lapack_complex_double* work,
					double* rwork, lapack_int *info );
void LAPACK_sposv( char* uplo, lapack_int* n, lapack_int* nrhs, float* a,
				lapack_int* lda, float* b, lapack_int* ldb,
				lapack_int *info );
void LAPACK_dposv( char* uplo, lapack_int* n, lapack_int* nrhs, double* a,
				lapack_int* lda, double* b, lapack_int* ldb,
				lapack_int *info );
void LAPACK_cposv( char* uplo, lapack_int* n, lapack_int* nrhs,
				lapack_complex_float* a, lapack_int* lda,
				lapack_complex_float* b, lapack_int* ldb, lapack_int *info );
void LAPACK_zposv( char* uplo, lapack_int* n, lapack_int* nrhs,
				lapack_complex_double* a, lapack_int* lda,
				lapack_complex_double* b, lapack_int* ldb,
				lapack_int *info );
void LAPACK_dsposv( char* uplo, lapack_int* n, lapack_int* nrhs, double* a,
					lapack_int* lda, double* b, lapack_int* ldb, double* x,
					lapack_int* ldx, double* work, float* swork,
					lapack_int* iter, lapack_int *info );
void LAPACK_zcposv( char* uplo, lapack_int* n, lapack_int* nrhs,
					lapack_complex_double* a, lapack_int* lda,
					lapack_complex_double* b, lapack_int* ldb,
					lapack_complex_double* x, lapack_int* ldx,
					lapack_complex_double* work, lapack_complex_float* swork,
					double* rwork, lapack_int* iter, lapack_int *info );
void LAPACK_sposvx( char* fact, char* uplo, lapack_int* n, lapack_int* nrhs,
					float* a, lapack_int* lda, float* af, lapack_int* ldaf,
					char* equed, float* s, float* b, lapack_int* ldb, float* x,
					lapack_int* ldx, float* rcond, float* ferr, float* berr,
					float* work, lapack_int* iwork, lapack_int *info );
void LAPACK_dposvx( char* fact, char* uplo, lapack_int* n, lapack_int* nrhs,
					double* a, lapack_int* lda, double* af, lapack_int* ldaf,
					char* equed, double* s, double* b, lapack_int* ldb,
					double* x, lapack_int* ldx, double* rcond, double* ferr,
					double* berr, double* work, lapack_int* iwork,
					lapack_int *info );
void LAPACK_cposvx( char* fact, char* uplo, lapack_int* n, lapack_int* nrhs,
					lapack_complex_float* a, lapack_int* lda,
					lapack_complex_float* af, lapack_int* ldaf, char* equed,
					float* s, lapack_complex_float* b, lapack_int* ldb,
					lapack_complex_float* x, lapack_int* ldx, float* rcond,
					float* ferr, float* berr, lapack_complex_float* work,
					float* rwork, lapack_int *info );
void LAPACK_zposvx( char* fact, char* uplo, lapack_int* n, lapack_int* nrhs,
					lapack_complex_double* a, lapack_int* lda,
					lapack_complex_double* af, lapack_int* ldaf, char* equed,
					double* s, lapack_complex_double* b, lapack_int* ldb,
					lapack_complex_double* x, lapack_int* ldx, double* rcond,
					double* ferr, double* berr, lapack_complex_double* work,
					double* rwork, lapack_int *info );
void LAPACK_dposvxx( char* fact, char* uplo, lapack_int* n, lapack_int* nrhs,
					double* a, lapack_int* lda, double* af, lapack_int* ldaf,
					char* equed, double* s, double* b, lapack_int* ldb,
					double* x, lapack_int* ldx, double* rcond, double* rpvgrw,
					double* berr, lapack_int* n_err_bnds,
					double* err_bnds_norm, double* err_bnds_comp,
					lapack_int* nparams, double* params, double* work,
					lapack_int* iwork, lapack_int *info );
void LAPACK_sposvxx( char* fact, char* uplo, lapack_int* n, lapack_int* nrhs,
					float* a, lapack_int* lda, float* af, lapack_int* ldaf,
					char* equed, float* s, float* b, lapack_int* ldb, float* x,
					lapack_int* ldx, float* rcond, float* rpvgrw, float* berr,
					lapack_int* n_err_bnds, float* err_bnds_norm,
					float* err_bnds_comp, lapack_int* nparams, float* params,
					float* work, lapack_int* iwork, lapack_int *info );
void LAPACK_zposvxx( char* fact, char* uplo, lapack_int* n, lapack_int* nrhs,
					lapack_complex_double* a, lapack_int* lda,
					lapack_complex_double* af, lapack_int* ldaf, char* equed,
					double* s, lapack_complex_double* b, lapack_int* ldb,
					lapack_complex_double* x, lapack_int* ldx, double* rcond,
					double* rpvgrw, double* berr, lapack_int* n_err_bnds,
					double* err_bnds_norm, double* err_bnds_comp,
					lapack_int* nparams, double* params,
					lapack_complex_double* work, double* rwork,
					lapack_int *info );
void LAPACK_cposvxx( char* fact, char* uplo, lapack_int* n, lapack_int* nrhs,
					lapack_complex_float* a, lapack_int* lda,
					lapack_complex_float* af, lapack_int* ldaf, char* equed,
					float* s, lapack_complex_float* b, lapack_int* ldb,
					lapack_complex_float* x, lapack_int* ldx, float* rcond,
					float* rpvgrw, float* berr, lapack_int* n_err_bnds,
					float* err_bnds_norm, float* err_bnds_comp,
					lapack_int* nparams, float* params,
					lapack_complex_float* work, float* rwork,
					lapack_int *info );
void LAPACK_sppsv( char* uplo, lapack_int* n, lapack_int* nrhs, float* ap,
				float* b, lapack_int* ldb, lapack_int *info );
void LAPACK_dppsv( char* uplo, lapack_int* n, lapack_int* nrhs, double* ap,
				double* b, lapack_int* ldb, lapack_int *info );
void LAPACK_cppsv( char* uplo, lapack_int* n, lapack_int* nrhs,
				lapack_complex_float* ap, lapack_complex_float* b,
				lapack_int* ldb, lapack_int *info );
void LAPACK_zppsv( char* uplo, lapack_int* n, lapack_int* nrhs,
				lapack_complex_double* ap, lapack_complex_double* b,
				lapack_int* ldb, lapack_int *info );
void LAPACK_sppsvx( char* fact, char* uplo, lapack_int* n, lapack_int* nrhs,
					float* ap, float* afp, char* equed, float* s, float* b,
					lapack_int* ldb, float* x, lapack_int* ldx, float* rcond,
					float* ferr, float* berr, float* work, lapack_int* iwork,
					lapack_int *info );
void LAPACK_dppsvx( char* fact, char* uplo, lapack_int* n, lapack_int* nrhs,
					double* ap, double* afp, char* equed, double* s, double* b,
					lapack_int* ldb, double* x, lapack_int* ldx, double* rcond,
					double* ferr, double* berr, double* work, lapack_int* iwork,
					lapack_int *info );
void LAPACK_cppsvx( char* fact, char* uplo, lapack_int* n, lapack_int* nrhs,
					lapack_complex_float* ap, lapack_complex_float* afp,
					char* equed, float* s, lapack_complex_float* b,
					lapack_int* ldb, lapack_complex_float* x, lapack_int* ldx,
					float* rcond, float* ferr, float* berr,
					lapack_complex_float* work, float* rwork,
					lapack_int *info );
void LAPACK_zppsvx( char* fact, char* uplo, lapack_int* n, lapack_int* nrhs,
					lapack_complex_double* ap, lapack_complex_double* afp,
					char* equed, double* s, lapack_complex_double* b,
					lapack_int* ldb, lapack_complex_double* x, lapack_int* ldx,
					double* rcond, double* ferr, double* berr,
					lapack_complex_double* work, double* rwork,
					lapack_int *info );
void LAPACK_spbsv( char* uplo, lapack_int* n, lapack_int* kd, lapack_int* nrhs,
				float* ab, lapack_int* ldab, float* b, lapack_int* ldb,
				lapack_int *info );
void LAPACK_dpbsv( char* uplo, lapack_int* n, lapack_int* kd, lapack_int* nrhs,
				double* ab, lapack_int* ldab, double* b, lapack_int* ldb,
				lapack_int *info );
void LAPACK_cpbsv( char* uplo, lapack_int* n, lapack_int* kd, lapack_int* nrhs,
				lapack_complex_float* ab, lapack_int* ldab,
				lapack_complex_float* b, lapack_int* ldb, lapack_int *info );
void LAPACK_zpbsv( char* uplo, lapack_int* n, lapack_int* kd, lapack_int* nrhs,
				lapack_complex_double* ab, lapack_int* ldab,
				lapack_complex_double* b, lapack_int* ldb,
				lapack_int *info );
void LAPACK_spbsvx( char* fact, char* uplo, lapack_int* n, lapack_int* kd,
					lapack_int* nrhs, float* ab, lapack_int* ldab, float* afb,
					lapack_int* ldafb, char* equed, float* s, float* b,
					lapack_int* ldb, float* x, lapack_int* ldx, float* rcond,
					float* ferr, float* berr, float* work, lapack_int* iwork,
					lapack_int *info );
void LAPACK_dpbsvx( char* fact, char* uplo, lapack_int* n, lapack_int* kd,
					lapack_int* nrhs, double* ab, lapack_int* ldab, double* afb,
					lapack_int* ldafb, char* equed, double* s, double* b,
					lapack_int* ldb, double* x, lapack_int* ldx, double* rcond,
					double* ferr, double* berr, double* work, lapack_int* iwork,
					lapack_int *info );
void LAPACK_cpbsvx( char* fact, char* uplo, lapack_int* n, lapack_int* kd,
					lapack_int* nrhs, lapack_complex_float* ab,
					lapack_int* ldab, lapack_complex_float* afb,
					lapack_int* ldafb, char* equed, float* s,
					lapack_complex_float* b, lapack_int* ldb,
					lapack_complex_float* x, lapack_int* ldx, float* rcond,
					float* ferr, float* berr, lapack_complex_float* work,
					float* rwork, lapack_int *info );
void LAPACK_zpbsvx( char* fact, char* uplo, lapack_int* n, lapack_int* kd,
					lapack_int* nrhs, lapack_complex_double* ab,
					lapack_int* ldab, lapack_complex_double* afb,
					lapack_int* ldafb, char* equed, double* s,
					lapack_complex_double* b, lapack_int* ldb,
					lapack_complex_double* x, lapack_int* ldx, double* rcond,
					double* ferr, double* berr, lapack_complex_double* work,
					double* rwork, lapack_int *info );
void LAPACK_sptsv( lapack_int* n, lapack_int* nrhs, float* d, float* e,
				float* b, lapack_int* ldb, lapack_int *info );
void LAPACK_dptsv( lapack_int* n, lapack_int* nrhs, double* d, double* e,
				double* b, lapack_int* ldb, lapack_int *info );
void LAPACK_cptsv( lapack_int* n, lapack_int* nrhs, float* d,
				lapack_complex_float* e, lapack_complex_float* b,
				lapack_int* ldb, lapack_int *info );
void LAPACK_zptsv( lapack_int* n, lapack_int* nrhs, double* d,
				lapack_complex_double* e, lapack_complex_double* b,
				lapack_int* ldb, lapack_int *info );
void LAPACK_sptsvx( char* fact, lapack_int* n, lapack_int* nrhs, const float* d,
					const float* e, float* df, float* ef, const float* b,
					lapack_int* ldb, float* x, lapack_int* ldx, float* rcond,
					float* ferr, float* berr, float* work, lapack_int *info );
void LAPACK_dptsvx( char* fact, lapack_int* n, lapack_int* nrhs,
					const double* d, const double* e, double* df, double* ef,
					const double* b, lapack_int* ldb, double* x,
					lapack_int* ldx, double* rcond, double* ferr, double* berr,
					double* work, lapack_int *info );
void LAPACK_cptsvx( char* fact, lapack_int* n, lapack_int* nrhs, const float* d,
					const lapack_complex_float* e, float* df,
					lapack_complex_float* ef, const lapack_complex_float* b,
					lapack_int* ldb, lapack_complex_float* x, lapack_int* ldx,
					float* rcond, float* ferr, float* berr,
					lapack_complex_float* work, float* rwork,
					lapack_int *info );
void LAPACK_zptsvx( char* fact, lapack_int* n, lapack_int* nrhs,
					const double* d, const lapack_complex_double* e, double* df,
					lapack_complex_double* ef, const lapack_complex_double* b,
					lapack_int* ldb, lapack_complex_double* x, lapack_int* ldx,
					double* rcond, double* ferr, double* berr,
					lapack_complex_double* work, double* rwork,
					lapack_int *info );
void LAPACK_ssysv( char* uplo, lapack_int* n, lapack_int* nrhs, float* a,
				lapack_int* lda, lapack_int* ipiv, float* b, lapack_int* ldb,
				float* work, lapack_int* lwork, lapack_int *info );
void LAPACK_dsysv( char* uplo, lapack_int* n, lapack_int* nrhs, double* a,
				lapack_int* lda, lapack_int* ipiv, double* b,
				lapack_int* ldb, double* work, lapack_int* lwork,
				lapack_int *info );
void LAPACK_csysv( char* uplo, lapack_int* n, lapack_int* nrhs,
				lapack_complex_float* a, lapack_int* lda, lapack_int* ipiv,
				lapack_complex_float* b, lapack_int* ldb,
				lapack_complex_float* work, lapack_int* lwork,
				lapack_int *info );
void LAPACK_zsysv( char* uplo, lapack_int* n, lapack_int* nrhs,
				lapack_complex_double* a, lapack_int* lda, lapack_int* ipiv,
				lapack_complex_double* b, lapack_int* ldb,
				lapack_complex_double* work, lapack_int* lwork,
				lapack_int *info );
void LAPACK_ssysvx( char* fact, char* uplo, lapack_int* n, lapack_int* nrhs,
					const float* a, lapack_int* lda, float* af,
					lapack_int* ldaf, lapack_int* ipiv, const float* b,
					lapack_int* ldb, float* x, lapack_int* ldx, float* rcond,
					float* ferr, float* berr, float* work, lapack_int* lwork,
					lapack_int* iwork, lapack_int *info );
void LAPACK_dsysvx( char* fact, char* uplo, lapack_int* n, lapack_int* nrhs,
					const double* a, lapack_int* lda, double* af,
					lapack_int* ldaf, lapack_int* ipiv, const double* b,
					lapack_int* ldb, double* x, lapack_int* ldx, double* rcond,
					double* ferr, double* berr, double* work, lapack_int* lwork,
					lapack_int* iwork, lapack_int *info );
void LAPACK_csysvx( char* fact, char* uplo, lapack_int* n, lapack_int* nrhs,
					const lapack_complex_float* a, lapack_int* lda,
					lapack_complex_float* af, lapack_int* ldaf,
					lapack_int* ipiv, const lapack_complex_float* b,
					lapack_int* ldb, lapack_complex_float* x, lapack_int* ldx,
					float* rcond, float* ferr, float* berr,
					lapack_complex_float* work, lapack_int* lwork, float* rwork,
					lapack_int *info );
void LAPACK_zsysvx( char* fact, char* uplo, lapack_int* n, lapack_int* nrhs,
					const lapack_complex_double* a, lapack_int* lda,
					lapack_complex_double* af, lapack_int* ldaf,
					lapack_int* ipiv, const lapack_complex_double* b,
					lapack_int* ldb, lapack_complex_double* x, lapack_int* ldx,
					double* rcond, double* ferr, double* berr,
					lapack_complex_double* work, lapack_int* lwork,
					double* rwork, lapack_int *info );
void LAPACK_dsysvxx( char* fact, char* uplo, lapack_int* n, lapack_int* nrhs,
					double* a, lapack_int* lda, double* af, lapack_int* ldaf,
					lapack_int* ipiv, char* equed, double* s, double* b,
					lapack_int* ldb, double* x, lapack_int* ldx, double* rcond,
					double* rpvgrw, double* berr, lapack_int* n_err_bnds,
					double* err_bnds_norm, double* err_bnds_comp,
					lapack_int* nparams, double* params, double* work,
					lapack_int* iwork, lapack_int *info );
void LAPACK_ssysvxx( char* fact, char* uplo, lapack_int* n, lapack_int* nrhs,
					float* a, lapack_int* lda, float* af, lapack_int* ldaf,
					lapack_int* ipiv, char* equed, float* s, float* b,
					lapack_int* ldb, float* x, lapack_int* ldx, float* rcond,
					float* rpvgrw, float* berr, lapack_int* n_err_bnds,
					float* err_bnds_norm, float* err_bnds_comp,
					lapack_int* nparams, float* params, float* work,
					lapack_int* iwork, lapack_int *info );
void LAPACK_zsysvxx( char* fact, char* uplo, lapack_int* n, lapack_int* nrhs,
					lapack_complex_double* a, lapack_int* lda,
					lapack_complex_double* af, lapack_int* ldaf,
					lapack_int* ipiv, char* equed, double* s,
					lapack_complex_double* b, lapack_int* ldb,
					lapack_complex_double* x, lapack_int* ldx, double* rcond,
					double* rpvgrw, double* berr, lapack_int* n_err_bnds,
					double* err_bnds_norm, double* err_bnds_comp,
					lapack_int* nparams, double* params,
					lapack_complex_double* work, double* rwork,
					lapack_int *info );
void LAPACK_csysvxx( char* fact, char* uplo, lapack_int* n, lapack_int* nrhs,
					lapack_complex_float* a, lapack_int* lda,
					lapack_complex_float* af, lapack_int* ldaf,
					lapack_int* ipiv, char* equed, float* s,
					lapack_complex_float* b, lapack_int* ldb,
					lapack_complex_float* x, lapack_int* ldx, float* rcond,
					float* rpvgrw, float* berr, lapack_int* n_err_bnds,
					float* err_bnds_norm, float* err_bnds_comp,
					lapack_int* nparams, float* params,
					lapack_complex_float* work, float* rwork,
					lapack_int *info );
void LAPACK_chesv( char* uplo, lapack_int* n, lapack_int* nrhs,
				lapack_complex_float* a, lapack_int* lda, lapack_int* ipiv,
				lapack_complex_float* b, lapack_int* ldb,
				lapack_complex_float* work, lapack_int* lwork,
				lapack_int *info );
void LAPACK_zhesv( char* uplo, lapack_int* n, lapack_int* nrhs,
				lapack_complex_double* a, lapack_int* lda, lapack_int* ipiv,
				lapack_complex_double* b, lapack_int* ldb,
				lapack_complex_double* work, lapack_int* lwork,
				lapack_int *info );
void LAPACK_chesvx( char* fact, char* uplo, lapack_int* n, lapack_int* nrhs,
					const lapack_complex_float* a, lapack_int* lda,
					lapack_complex_float* af, lapack_int* ldaf,
					lapack_int* ipiv, const lapack_complex_float* b,
					lapack_int* ldb, lapack_complex_float* x, lapack_int* ldx,
					float* rcond, float* ferr, float* berr,
					lapack_complex_float* work, lapack_int* lwork, float* rwork,
					lapack_int *info );
void LAPACK_zhesvx( char* fact, char* uplo, lapack_int* n, lapack_int* nrhs,
					const lapack_complex_double* a, lapack_int* lda,
					lapack_complex_double* af, lapack_int* ldaf,
					lapack_int* ipiv, const lapack_complex_double* b,
					lapack_int* ldb, lapack_complex_double* x, lapack_int* ldx,
					double* rcond, double* ferr, double* berr,
					lapack_complex_double* work, lapack_int* lwork,
					double* rwork, lapack_int *info );
void LAPACK_zhesvxx( char* fact, char* uplo, lapack_int* n, lapack_int* nrhs,
					lapack_complex_double* a, lapack_int* lda,
					lapack_complex_double* af, lapack_int* ldaf,
					lapack_int* ipiv, char* equed, double* s,
					lapack_complex_double* b, lapack_int* ldb,
					lapack_complex_double* x, lapack_int* ldx, double* rcond,
					double* rpvgrw, double* berr, lapack_int* n_err_bnds,
					double* err_bnds_norm, double* err_bnds_comp,
					lapack_int* nparams, double* params,
					lapack_complex_double* work, double* rwork,
					lapack_int *info );
void LAPACK_chesvxx( char* fact, char* uplo, lapack_int* n, lapack_int* nrhs,
					lapack_complex_float* a, lapack_int* lda,
					lapack_complex_float* af, lapack_int* ldaf,
					lapack_int* ipiv, char* equed, float* s,
					lapack_complex_float* b, lapack_int* ldb,
					lapack_complex_float* x, lapack_int* ldx, float* rcond,
					float* rpvgrw, float* berr, lapack_int* n_err_bnds,
					float* err_bnds_norm, float* err_bnds_comp,
					lapack_int* nparams, float* params,
					lapack_complex_float* work, float* rwork,
					lapack_int *info );
void LAPACK_sspsv( char* uplo, lapack_int* n, lapack_int* nrhs, float* ap,
				lapack_int* ipiv, float* b, lapack_int* ldb,
				lapack_int *info );
void LAPACK_dspsv( char* uplo, lapack_int* n, lapack_int* nrhs, double* ap,
				lapack_int* ipiv, double* b, lapack_int* ldb,
				lapack_int *info );
void LAPACK_cspsv( char* uplo, lapack_int* n, lapack_int* nrhs,
				lapack_complex_float* ap, lapack_int* ipiv,
				lapack_complex_float* b, lapack_int* ldb, lapack_int *info );
void LAPACK_zspsv( char* uplo, lapack_int* n, lapack_int* nrhs,
				lapack_complex_double* ap, lapack_int* ipiv,
				lapack_complex_double* b, lapack_int* ldb,
				lapack_int *info );
void LAPACK_sspsvx( char* fact, char* uplo, lapack_int* n, lapack_int* nrhs,
					const float* ap, float* afp, lapack_int* ipiv,
					const float* b, lapack_int* ldb, float* x, lapack_int* ldx,
					float* rcond, float* ferr, float* berr, float* work,
					lapack_int* iwork, lapack_int *info );
void LAPACK_dspsvx( char* fact, char* uplo, lapack_int* n, lapack_int* nrhs,
					const double* ap, double* afp, lapack_int* ipiv,
					const double* b, lapack_int* ldb, double* x,
					lapack_int* ldx, double* rcond, double* ferr, double* berr,
					double* work, lapack_int* iwork, lapack_int *info );
void LAPACK_cspsvx( char* fact, char* uplo, lapack_int* n, lapack_int* nrhs,
					const lapack_complex_float* ap, lapack_complex_float* afp,
					lapack_int* ipiv, const lapack_complex_float* b,
					lapack_int* ldb, lapack_complex_float* x, lapack_int* ldx,
					float* rcond, float* ferr, float* berr,
					lapack_complex_float* work, float* rwork,
					lapack_int *info );
void LAPACK_zspsvx( char* fact, char* uplo, lapack_int* n, lapack_int* nrhs,
					const lapack_complex_double* ap, lapack_complex_double* afp,
					lapack_int* ipiv, const lapack_complex_double* b,
					lapack_int* ldb, lapack_complex_double* x, lapack_int* ldx,
					double* rcond, double* ferr, double* berr,
					lapack_complex_double* work, double* rwork,
					lapack_int *info );
void LAPACK_chpsv( char* uplo, lapack_int* n, lapack_int* nrhs,
				lapack_complex_float* ap, lapack_int* ipiv,
				lapack_complex_float* b, lapack_int* ldb, lapack_int *info );
void LAPACK_zhpsv( char* uplo, lapack_int* n, lapack_int* nrhs,
				lapack_complex_double* ap, lapack_int* ipiv,
				lapack_complex_double* b, lapack_int* ldb,
				lapack_int *info );
void LAPACK_chpsvx( char* fact, char* uplo, lapack_int* n, lapack_int* nrhs,
					const lapack_complex_float* ap, lapack_complex_float* afp,
					lapack_int* ipiv, const lapack_complex_float* b,
					lapack_int* ldb, lapack_complex_float* x, lapack_int* ldx,
					float* rcond, float* ferr, float* berr,
					lapack_complex_float* work, float* rwork,
					lapack_int *info );
void LAPACK_zhpsvx( char* fact, char* uplo, lapack_int* n, lapack_int* nrhs,
					const lapack_complex_double* ap, lapack_complex_double* afp,
					lapack_int* ipiv, const lapack_complex_double* b,
					lapack_int* ldb, lapack_complex_double* x, lapack_int* ldx,
					double* rcond, double* ferr, double* berr,
					lapack_complex_double* work, double* rwork,
					lapack_int *info );
void LAPACK_sgeqrf( lapack_int* m, lapack_int* n, float* a, lapack_int* lda,
					float* tau, float* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_dgeqrf( lapack_int* m, lapack_int* n, double* a, lapack_int* lda,
					double* tau, double* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_cgeqrf( lapack_int* m, lapack_int* n, lapack_complex_float* a,
					lapack_int* lda, lapack_complex_float* tau,
					lapack_complex_float* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_zgeqrf( lapack_int* m, lapack_int* n, lapack_complex_double* a,
					lapack_int* lda, lapack_complex_double* tau,
					lapack_complex_double* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_sgeqpf( lapack_int* m, lapack_int* n, float* a, lapack_int* lda,
					lapack_int* jpvt, float* tau, float* work,
					lapack_int *info );
void LAPACK_dgeqpf( lapack_int* m, lapack_int* n, double* a, lapack_int* lda,
					lapack_int* jpvt, double* tau, double* work,
					lapack_int *info );
void LAPACK_cgeqpf( lapack_int* m, lapack_int* n, lapack_complex_float* a,
					lapack_int* lda, lapack_int* jpvt,
					lapack_complex_float* tau, lapack_complex_float* work,
					float* rwork, lapack_int *info );
void LAPACK_zgeqpf( lapack_int* m, lapack_int* n, lapack_complex_double* a,
					lapack_int* lda, lapack_int* jpvt,
					lapack_complex_double* tau, lapack_complex_double* work,
					double* rwork, lapack_int *info );
void LAPACK_sgeqp3( lapack_int* m, lapack_int* n, float* a, lapack_int* lda,
					lapack_int* jpvt, float* tau, float* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_dgeqp3( lapack_int* m, lapack_int* n, double* a, lapack_int* lda,
					lapack_int* jpvt, double* tau, double* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_cgeqp3( lapack_int* m, lapack_int* n, lapack_complex_float* a,
					lapack_int* lda, lapack_int* jpvt,
					lapack_complex_float* tau, lapack_complex_float* work,
					lapack_int* lwork, float* rwork, lapack_int *info );
void LAPACK_zgeqp3( lapack_int* m, lapack_int* n, lapack_complex_double* a,
					lapack_int* lda, lapack_int* jpvt,
					lapack_complex_double* tau, lapack_complex_double* work,
					lapack_int* lwork, double* rwork, lapack_int *info );
void LAPACK_sorgqr( lapack_int* m, lapack_int* n, lapack_int* k, float* a,
					lapack_int* lda, const float* tau, float* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_dorgqr( lapack_int* m, lapack_int* n, lapack_int* k, double* a,
					lapack_int* lda, const double* tau, double* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_sormqr( char* side, char* trans, lapack_int* m, lapack_int* n,
					lapack_int* k, const float* a, lapack_int* lda,
					const float* tau, float* c, lapack_int* ldc, float* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_dormqr( char* side, char* trans, lapack_int* m, lapack_int* n,
					lapack_int* k, const double* a, lapack_int* lda,
					const double* tau, double* c, lapack_int* ldc, double* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_cungqr( lapack_int* m, lapack_int* n, lapack_int* k,
					lapack_complex_float* a, lapack_int* lda,
					const lapack_complex_float* tau, lapack_complex_float* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_zungqr( lapack_int* m, lapack_int* n, lapack_int* k,
					lapack_complex_double* a, lapack_int* lda,
					const lapack_complex_double* tau,
					lapack_complex_double* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_cunmqr( char* side, char* trans, lapack_int* m, lapack_int* n,
					lapack_int* k, const lapack_complex_float* a,
					lapack_int* lda, const lapack_complex_float* tau,
					lapack_complex_float* c, lapack_int* ldc,
					lapack_complex_float* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_zunmqr( char* side, char* trans, lapack_int* m, lapack_int* n,
					lapack_int* k, const lapack_complex_double* a,
					lapack_int* lda, const lapack_complex_double* tau,
					lapack_complex_double* c, lapack_int* ldc,
					lapack_complex_double* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_sgelqf( lapack_int* m, lapack_int* n, float* a, lapack_int* lda,
					float* tau, float* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_dgelqf( lapack_int* m, lapack_int* n, double* a, lapack_int* lda,
					double* tau, double* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_cgelqf( lapack_int* m, lapack_int* n, lapack_complex_float* a,
					lapack_int* lda, lapack_complex_float* tau,
					lapack_complex_float* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_zgelqf( lapack_int* m, lapack_int* n, lapack_complex_double* a,
					lapack_int* lda, lapack_complex_double* tau,
					lapack_complex_double* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_sorglq( lapack_int* m, lapack_int* n, lapack_int* k, float* a,
					lapack_int* lda, const float* tau, float* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_dorglq( lapack_int* m, lapack_int* n, lapack_int* k, double* a,
					lapack_int* lda, const double* tau, double* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_sormlq( char* side, char* trans, lapack_int* m, lapack_int* n,
					lapack_int* k, const float* a, lapack_int* lda,
					const float* tau, float* c, lapack_int* ldc, float* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_dormlq( char* side, char* trans, lapack_int* m, lapack_int* n,
					lapack_int* k, const double* a, lapack_int* lda,
					const double* tau, double* c, lapack_int* ldc, double* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_cunglq( lapack_int* m, lapack_int* n, lapack_int* k,
					lapack_complex_float* a, lapack_int* lda,
					const lapack_complex_float* tau, lapack_complex_float* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_zunglq( lapack_int* m, lapack_int* n, lapack_int* k,
					lapack_complex_double* a, lapack_int* lda,
					const lapack_complex_double* tau,
					lapack_complex_double* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_cunmlq( char* side, char* trans, lapack_int* m, lapack_int* n,
					lapack_int* k, const lapack_complex_float* a,
					lapack_int* lda, const lapack_complex_float* tau,
					lapack_complex_float* c, lapack_int* ldc,
					lapack_complex_float* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_zunmlq( char* side, char* trans, lapack_int* m, lapack_int* n,
					lapack_int* k, const lapack_complex_double* a,
					lapack_int* lda, const lapack_complex_double* tau,
					lapack_complex_double* c, lapack_int* ldc,
					lapack_complex_double* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_sgeqlf( lapack_int* m, lapack_int* n, float* a, lapack_int* lda,
					float* tau, float* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_dgeqlf( lapack_int* m, lapack_int* n, double* a, lapack_int* lda,
					double* tau, double* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_cgeqlf( lapack_int* m, lapack_int* n, lapack_complex_float* a,
					lapack_int* lda, lapack_complex_float* tau,
					lapack_complex_float* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_zgeqlf( lapack_int* m, lapack_int* n, lapack_complex_double* a,
					lapack_int* lda, lapack_complex_double* tau,
					lapack_complex_double* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_sorgql( lapack_int* m, lapack_int* n, lapack_int* k, float* a,
					lapack_int* lda, const float* tau, float* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_dorgql( lapack_int* m, lapack_int* n, lapack_int* k, double* a,
					lapack_int* lda, const double* tau, double* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_cungql( lapack_int* m, lapack_int* n, lapack_int* k,
					lapack_complex_float* a, lapack_int* lda,
					const lapack_complex_float* tau, lapack_complex_float* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_zungql( lapack_int* m, lapack_int* n, lapack_int* k,
					lapack_complex_double* a, lapack_int* lda,
					const lapack_complex_double* tau,
					lapack_complex_double* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_sormql( char* side, char* trans, lapack_int* m, lapack_int* n,
					lapack_int* k, const float* a, lapack_int* lda,
					const float* tau, float* c, lapack_int* ldc, float* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_dormql( char* side, char* trans, lapack_int* m, lapack_int* n,
					lapack_int* k, const double* a, lapack_int* lda,
					const double* tau, double* c, lapack_int* ldc, double* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_cunmql( char* side, char* trans, lapack_int* m, lapack_int* n,
					lapack_int* k, const lapack_complex_float* a,
					lapack_int* lda, const lapack_complex_float* tau,
					lapack_complex_float* c, lapack_int* ldc,
					lapack_complex_float* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_zunmql( char* side, char* trans, lapack_int* m, lapack_int* n,
					lapack_int* k, const lapack_complex_double* a,
					lapack_int* lda, const lapack_complex_double* tau,
					lapack_complex_double* c, lapack_int* ldc,
					lapack_complex_double* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_sgerqf( lapack_int* m, lapack_int* n, float* a, lapack_int* lda,
					float* tau, float* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_dgerqf( lapack_int* m, lapack_int* n, double* a, lapack_int* lda,
					double* tau, double* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_cgerqf( lapack_int* m, lapack_int* n, lapack_complex_float* a,
					lapack_int* lda, lapack_complex_float* tau,
					lapack_complex_float* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_zgerqf( lapack_int* m, lapack_int* n, lapack_complex_double* a,
					lapack_int* lda, lapack_complex_double* tau,
					lapack_complex_double* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_sorgrq( lapack_int* m, lapack_int* n, lapack_int* k, float* a,
					lapack_int* lda, const float* tau, float* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_dorgrq( lapack_int* m, lapack_int* n, lapack_int* k, double* a,
					lapack_int* lda, const double* tau, double* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_cungrq( lapack_int* m, lapack_int* n, lapack_int* k,
					lapack_complex_float* a, lapack_int* lda,
					const lapack_complex_float* tau, lapack_complex_float* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_zungrq( lapack_int* m, lapack_int* n, lapack_int* k,
					lapack_complex_double* a, lapack_int* lda,
					const lapack_complex_double* tau,
					lapack_complex_double* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_sormrq( char* side, char* trans, lapack_int* m, lapack_int* n,
					lapack_int* k, const float* a, lapack_int* lda,
					const float* tau, float* c, lapack_int* ldc, float* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_dormrq( char* side, char* trans, lapack_int* m, lapack_int* n,
					lapack_int* k, const double* a, lapack_int* lda,
					const double* tau, double* c, lapack_int* ldc, double* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_cunmrq( char* side, char* trans, lapack_int* m, lapack_int* n,
					lapack_int* k, const lapack_complex_float* a,
					lapack_int* lda, const lapack_complex_float* tau,
					lapack_complex_float* c, lapack_int* ldc,
					lapack_complex_float* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_zunmrq( char* side, char* trans, lapack_int* m, lapack_int* n,
					lapack_int* k, const lapack_complex_double* a,
					lapack_int* lda, const lapack_complex_double* tau,
					lapack_complex_double* c, lapack_int* ldc,
					lapack_complex_double* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_stzrzf( lapack_int* m, lapack_int* n, float* a, lapack_int* lda,
					float* tau, float* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_dtzrzf( lapack_int* m, lapack_int* n, double* a, lapack_int* lda,
					double* tau, double* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_ctzrzf( lapack_int* m, lapack_int* n, lapack_complex_float* a,
					lapack_int* lda, lapack_complex_float* tau,
					lapack_complex_float* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_ztzrzf( lapack_int* m, lapack_int* n, lapack_complex_double* a,
					lapack_int* lda, lapack_complex_double* tau,
					lapack_complex_double* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_sormrz( char* side, char* trans, lapack_int* m, lapack_int* n,
					lapack_int* k, lapack_int* l, const float* a,
					lapack_int* lda, const float* tau, float* c,
					lapack_int* ldc, float* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_dormrz( char* side, char* trans, lapack_int* m, lapack_int* n,
					lapack_int* k, lapack_int* l, const double* a,
					lapack_int* lda, const double* tau, double* c,
					lapack_int* ldc, double* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_cunmrz( char* side, char* trans, lapack_int* m, lapack_int* n,
					lapack_int* k, lapack_int* l, const lapack_complex_float* a,
					lapack_int* lda, const lapack_complex_float* tau,
					lapack_complex_float* c, lapack_int* ldc,
					lapack_complex_float* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_zunmrz( char* side, char* trans, lapack_int* m, lapack_int* n,
					lapack_int* k, lapack_int* l,
					const lapack_complex_double* a, lapack_int* lda,
					const lapack_complex_double* tau, lapack_complex_double* c,
					lapack_int* ldc, lapack_complex_double* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_sggqrf( lapack_int* n, lapack_int* m, lapack_int* p, float* a,
					lapack_int* lda, float* taua, float* b, lapack_int* ldb,
					float* taub, float* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_dggqrf( lapack_int* n, lapack_int* m, lapack_int* p, double* a,
					lapack_int* lda, double* taua, double* b, lapack_int* ldb,
					double* taub, double* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_cggqrf( lapack_int* n, lapack_int* m, lapack_int* p,
					lapack_complex_float* a, lapack_int* lda,
					lapack_complex_float* taua, lapack_complex_float* b,
					lapack_int* ldb, lapack_complex_float* taub,
					lapack_complex_float* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_zggqrf( lapack_int* n, lapack_int* m, lapack_int* p,
					lapack_complex_double* a, lapack_int* lda,
					lapack_complex_double* taua, lapack_complex_double* b,
					lapack_int* ldb, lapack_complex_double* taub,
					lapack_complex_double* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_sggrqf( lapack_int* m, lapack_int* p, lapack_int* n, float* a,
					lapack_int* lda, float* taua, float* b, lapack_int* ldb,
					float* taub, float* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_dggrqf( lapack_int* m, lapack_int* p, lapack_int* n, double* a,
					lapack_int* lda, double* taua, double* b, lapack_int* ldb,
					double* taub, double* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_cggrqf( lapack_int* m, lapack_int* p, lapack_int* n,
					lapack_complex_float* a, lapack_int* lda,
					lapack_complex_float* taua, lapack_complex_float* b,
					lapack_int* ldb, lapack_complex_float* taub,
					lapack_complex_float* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_zggrqf( lapack_int* m, lapack_int* p, lapack_int* n,
					lapack_complex_double* a, lapack_int* lda,
					lapack_complex_double* taua, lapack_complex_double* b,
					lapack_int* ldb, lapack_complex_double* taub,
					lapack_complex_double* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_sgebrd( lapack_int* m, lapack_int* n, float* a, lapack_int* lda,
					float* d, float* e, float* tauq, float* taup, float* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_dgebrd( lapack_int* m, lapack_int* n, double* a, lapack_int* lda,
					double* d, double* e, double* tauq, double* taup,
					double* work, lapack_int* lwork, lapack_int *info );
void LAPACK_cgebrd( lapack_int* m, lapack_int* n, lapack_complex_float* a,
					lapack_int* lda, float* d, float* e,
					lapack_complex_float* tauq, lapack_complex_float* taup,
					lapack_complex_float* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_zgebrd( lapack_int* m, lapack_int* n, lapack_complex_double* a,
					lapack_int* lda, double* d, double* e,
					lapack_complex_double* tauq, lapack_complex_double* taup,
					lapack_complex_double* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_sgbbrd( char* vect, lapack_int* m, lapack_int* n, lapack_int* ncc,
					lapack_int* kl, lapack_int* ku, float* ab, lapack_int* ldab,
					float* d, float* e, float* q, lapack_int* ldq, float* pt,
					lapack_int* ldpt, float* c, lapack_int* ldc, float* work,
					lapack_int *info );
void LAPACK_dgbbrd( char* vect, lapack_int* m, lapack_int* n, lapack_int* ncc,
					lapack_int* kl, lapack_int* ku, double* ab,
					lapack_int* ldab, double* d, double* e, double* q,
					lapack_int* ldq, double* pt, lapack_int* ldpt, double* c,
					lapack_int* ldc, double* work, lapack_int *info );
void LAPACK_cgbbrd( char* vect, lapack_int* m, lapack_int* n, lapack_int* ncc,
					lapack_int* kl, lapack_int* ku, lapack_complex_float* ab,
					lapack_int* ldab, float* d, float* e,
					lapack_complex_float* q, lapack_int* ldq,
					lapack_complex_float* pt, lapack_int* ldpt,
					lapack_complex_float* c, lapack_int* ldc,
					lapack_complex_float* work, float* rwork,
					lapack_int *info );
void LAPACK_zgbbrd( char* vect, lapack_int* m, lapack_int* n, lapack_int* ncc,
					lapack_int* kl, lapack_int* ku, lapack_complex_double* ab,
					lapack_int* ldab, double* d, double* e,
					lapack_complex_double* q, lapack_int* ldq,
					lapack_complex_double* pt, lapack_int* ldpt,
					lapack_complex_double* c, lapack_int* ldc,
					lapack_complex_double* work, double* rwork,
					lapack_int *info );
void LAPACK_sorgbr( char* vect, lapack_int* m, lapack_int* n, lapack_int* k,
					float* a, lapack_int* lda, const float* tau, float* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_dorgbr( char* vect, lapack_int* m, lapack_int* n, lapack_int* k,
					double* a, lapack_int* lda, const double* tau, double* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_sormbr( char* vect, char* side, char* trans, lapack_int* m,
					lapack_int* n, lapack_int* k, const float* a,
					lapack_int* lda, const float* tau, float* c,
					lapack_int* ldc, float* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_dormbr( char* vect, char* side, char* trans, lapack_int* m,
					lapack_int* n, lapack_int* k, const double* a,
					lapack_int* lda, const double* tau, double* c,
					lapack_int* ldc, double* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_cungbr( char* vect, lapack_int* m, lapack_int* n, lapack_int* k,
					lapack_complex_float* a, lapack_int* lda,
					const lapack_complex_float* tau, lapack_complex_float* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_zungbr( char* vect, lapack_int* m, lapack_int* n, lapack_int* k,
					lapack_complex_double* a, lapack_int* lda,
					const lapack_complex_double* tau,
					lapack_complex_double* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_cunmbr( char* vect, char* side, char* trans, lapack_int* m,
					lapack_int* n, lapack_int* k, const lapack_complex_float* a,
					lapack_int* lda, const lapack_complex_float* tau,
					lapack_complex_float* c, lapack_int* ldc,
					lapack_complex_float* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_zunmbr( char* vect, char* side, char* trans, lapack_int* m,
					lapack_int* n, lapack_int* k,
					const lapack_complex_double* a, lapack_int* lda,
					const lapack_complex_double* tau, lapack_complex_double* c,
					lapack_int* ldc, lapack_complex_double* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_sbdsqr( char* uplo, lapack_int* n, lapack_int* ncvt,
					lapack_int* nru, lapack_int* ncc, float* d, float* e,
					float* vt, lapack_int* ldvt, float* u, lapack_int* ldu,
					float* c, lapack_int* ldc, float* work, lapack_int *info );
void LAPACK_dbdsqr( char* uplo, lapack_int* n, lapack_int* ncvt,
					lapack_int* nru, lapack_int* ncc, double* d, double* e,
					double* vt, lapack_int* ldvt, double* u, lapack_int* ldu,
					double* c, lapack_int* ldc, double* work,
					lapack_int *info );
void LAPACK_cbdsqr( char* uplo, lapack_int* n, lapack_int* ncvt,
					lapack_int* nru, lapack_int* ncc, float* d, float* e,
					lapack_complex_float* vt, lapack_int* ldvt,
					lapack_complex_float* u, lapack_int* ldu,
					lapack_complex_float* c, lapack_int* ldc, float* work,
					lapack_int *info );
void LAPACK_zbdsqr( char* uplo, lapack_int* n, lapack_int* ncvt,
					lapack_int* nru, lapack_int* ncc, double* d, double* e,
					lapack_complex_double* vt, lapack_int* ldvt,
					lapack_complex_double* u, lapack_int* ldu,
					lapack_complex_double* c, lapack_int* ldc, double* work,
					lapack_int *info );
void LAPACK_sbdsdc( char* uplo, char* compq, lapack_int* n, float* d, float* e,
					float* u, lapack_int* ldu, float* vt, lapack_int* ldvt,
					float* q, lapack_int* iq, float* work, lapack_int* iwork,
					lapack_int *info );
void LAPACK_dbdsdc( char* uplo, char* compq, lapack_int* n, double* d,
					double* e, double* u, lapack_int* ldu, double* vt,
					lapack_int* ldvt, double* q, lapack_int* iq, double* work,
					lapack_int* iwork, lapack_int *info );
void LAPACK_ssytrd( char* uplo, lapack_int* n, float* a, lapack_int* lda,
					float* d, float* e, float* tau, float* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_dsytrd( char* uplo, lapack_int* n, double* a, lapack_int* lda,
					double* d, double* e, double* tau, double* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_sorgtr( char* uplo, lapack_int* n, float* a, lapack_int* lda,
					const float* tau, float* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_dorgtr( char* uplo, lapack_int* n, double* a, lapack_int* lda,
					const double* tau, double* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_sormtr( char* side, char* uplo, char* trans, lapack_int* m,
					lapack_int* n, const float* a, lapack_int* lda,
					const float* tau, float* c, lapack_int* ldc, float* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_dormtr( char* side, char* uplo, char* trans, lapack_int* m,
					lapack_int* n, const double* a, lapack_int* lda,
					const double* tau, double* c, lapack_int* ldc, double* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_chetrd( char* uplo, lapack_int* n, lapack_complex_float* a,
					lapack_int* lda, float* d, float* e,
					lapack_complex_float* tau, lapack_complex_float* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_zhetrd( char* uplo, lapack_int* n, lapack_complex_double* a,
					lapack_int* lda, double* d, double* e,
					lapack_complex_double* tau, lapack_complex_double* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_cungtr( char* uplo, lapack_int* n, lapack_complex_float* a,
					lapack_int* lda, const lapack_complex_float* tau,
					lapack_complex_float* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_zungtr( char* uplo, lapack_int* n, lapack_complex_double* a,
					lapack_int* lda, const lapack_complex_double* tau,
					lapack_complex_double* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_cunmtr( char* side, char* uplo, char* trans, lapack_int* m,
					lapack_int* n, const lapack_complex_float* a,
					lapack_int* lda, const lapack_complex_float* tau,
					lapack_complex_float* c, lapack_int* ldc,
					lapack_complex_float* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_zunmtr( char* side, char* uplo, char* trans, lapack_int* m,
					lapack_int* n, const lapack_complex_double* a,
					lapack_int* lda, const lapack_complex_double* tau,
					lapack_complex_double* c, lapack_int* ldc,
					lapack_complex_double* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_ssptrd( char* uplo, lapack_int* n, float* ap, float* d, float* e,
					float* tau, lapack_int *info );
void LAPACK_dsptrd( char* uplo, lapack_int* n, double* ap, double* d, double* e,
					double* tau, lapack_int *info );
void LAPACK_sopgtr( char* uplo, lapack_int* n, const float* ap,
					const float* tau, float* q, lapack_int* ldq, float* work,
					lapack_int *info );
void LAPACK_dopgtr( char* uplo, lapack_int* n, const double* ap,
					const double* tau, double* q, lapack_int* ldq, double* work,
					lapack_int *info );
void LAPACK_sopmtr( char* side, char* uplo, char* trans, lapack_int* m,
					lapack_int* n, const float* ap, const float* tau, float* c,
					lapack_int* ldc, float* work, lapack_int *info );
void LAPACK_dopmtr( char* side, char* uplo, char* trans, lapack_int* m,
					lapack_int* n, const double* ap, const double* tau,
					double* c, lapack_int* ldc, double* work,
					lapack_int *info );
void LAPACK_chptrd( char* uplo, lapack_int* n, lapack_complex_float* ap,
					float* d, float* e, lapack_complex_float* tau,
					lapack_int *info );
void LAPACK_zhptrd( char* uplo, lapack_int* n, lapack_complex_double* ap,
					double* d, double* e, lapack_complex_double* tau,
					lapack_int *info );
void LAPACK_cupgtr( char* uplo, lapack_int* n, const lapack_complex_float* ap,
					const lapack_complex_float* tau, lapack_complex_float* q,
					lapack_int* ldq, lapack_complex_float* work,
					lapack_int *info );
void LAPACK_zupgtr( char* uplo, lapack_int* n, const lapack_complex_double* ap,
					const lapack_complex_double* tau, lapack_complex_double* q,
					lapack_int* ldq, lapack_complex_double* work,
					lapack_int *info );
void LAPACK_cupmtr( char* side, char* uplo, char* trans, lapack_int* m,
					lapack_int* n, const lapack_complex_float* ap,
					const lapack_complex_float* tau, lapack_complex_float* c,
					lapack_int* ldc, lapack_complex_float* work,
					lapack_int *info );
void LAPACK_zupmtr( char* side, char* uplo, char* trans, lapack_int* m,
					lapack_int* n, const lapack_complex_double* ap,
					const lapack_complex_double* tau, lapack_complex_double* c,
					lapack_int* ldc, lapack_complex_double* work,
					lapack_int *info );
void LAPACK_ssbtrd( char* vect, char* uplo, lapack_int* n, lapack_int* kd,
					float* ab, lapack_int* ldab, float* d, float* e, float* q,
					lapack_int* ldq, float* work, lapack_int *info );
void LAPACK_dsbtrd( char* vect, char* uplo, lapack_int* n, lapack_int* kd,
					double* ab, lapack_int* ldab, double* d, double* e,
					double* q, lapack_int* ldq, double* work,
					lapack_int *info );
void LAPACK_chbtrd( char* vect, char* uplo, lapack_int* n, lapack_int* kd,
					lapack_complex_float* ab, lapack_int* ldab, float* d,
					float* e, lapack_complex_float* q, lapack_int* ldq,
					lapack_complex_float* work, lapack_int *info );
void LAPACK_zhbtrd( char* vect, char* uplo, lapack_int* n, lapack_int* kd,
					lapack_complex_double* ab, lapack_int* ldab, double* d,
					double* e, lapack_complex_double* q, lapack_int* ldq,
					lapack_complex_double* work, lapack_int *info );
void LAPACK_ssterf( lapack_int* n, float* d, float* e, lapack_int *info );
void LAPACK_dsterf( lapack_int* n, double* d, double* e, lapack_int *info );
void LAPACK_ssteqr( char* compz, lapack_int* n, float* d, float* e, float* z,
					lapack_int* ldz, float* work, lapack_int *info );
void LAPACK_dsteqr( char* compz, lapack_int* n, double* d, double* e, double* z,
					lapack_int* ldz, double* work, lapack_int *info );
void LAPACK_csteqr( char* compz, lapack_int* n, float* d, float* e,
					lapack_complex_float* z, lapack_int* ldz, float* work,
					lapack_int *info );
void LAPACK_zsteqr( char* compz, lapack_int* n, double* d, double* e,
					lapack_complex_double* z, lapack_int* ldz, double* work,
					lapack_int *info );
void LAPACK_sstemr( char* jobz, char* range, lapack_int* n, float* d, float* e,
					float* vl, float* vu, lapack_int* il, lapack_int* iu,
					lapack_int* m, float* w, float* z, lapack_int* ldz,
					lapack_int* nzc, lapack_int* isuppz, lapack_logical* tryrac,
					float* work, lapack_int* lwork, lapack_int* iwork,
					lapack_int* liwork, lapack_int *info );
void LAPACK_dstemr( char* jobz, char* range, lapack_int* n, double* d,
					double* e, double* vl, double* vu, lapack_int* il,
					lapack_int* iu, lapack_int* m, double* w, double* z,
					lapack_int* ldz, lapack_int* nzc, lapack_int* isuppz,
					lapack_logical* tryrac, double* work, lapack_int* lwork,
					lapack_int* iwork, lapack_int* liwork, lapack_int *info );
void LAPACK_cstemr( char* jobz, char* range, lapack_int* n, float* d, float* e,
					float* vl, float* vu, lapack_int* il, lapack_int* iu,
					lapack_int* m, float* w, lapack_complex_float* z,
					lapack_int* ldz, lapack_int* nzc, lapack_int* isuppz,
					lapack_logical* tryrac, float* work, lapack_int* lwork,
					lapack_int* iwork, lapack_int* liwork, lapack_int *info );
void LAPACK_zstemr( char* jobz, char* range, lapack_int* n, double* d,
					double* e, double* vl, double* vu, lapack_int* il,
					lapack_int* iu, lapack_int* m, double* w,
					lapack_complex_double* z, lapack_int* ldz, lapack_int* nzc,
					lapack_int* isuppz, lapack_logical* tryrac, double* work,
					lapack_int* lwork, lapack_int* iwork, lapack_int* liwork,
					lapack_int *info );
void LAPACK_sstedc( char* compz, lapack_int* n, float* d, float* e, float* z,
					lapack_int* ldz, float* work, lapack_int* lwork,
					lapack_int* iwork, lapack_int* liwork, lapack_int *info );
void LAPACK_dstedc( char* compz, lapack_int* n, double* d, double* e, double* z,
					lapack_int* ldz, double* work, lapack_int* lwork,
					lapack_int* iwork, lapack_int* liwork, lapack_int *info );
void LAPACK_cstedc( char* compz, lapack_int* n, float* d, float* e,
					lapack_complex_float* z, lapack_int* ldz,
					lapack_complex_float* work, lapack_int* lwork, float* rwork,
					lapack_int* lrwork, lapack_int* iwork, lapack_int* liwork,
					lapack_int *info );
void LAPACK_zstedc( char* compz, lapack_int* n, double* d, double* e,
					lapack_complex_double* z, lapack_int* ldz,
					lapack_complex_double* work, lapack_int* lwork,
					double* rwork, lapack_int* lrwork, lapack_int* iwork,
					lapack_int* liwork, lapack_int *info );
void LAPACK_sstegr( char* jobz, char* range, lapack_int* n, float* d, float* e,
					float* vl, float* vu, lapack_int* il, lapack_int* iu,
					float* abstol, lapack_int* m, float* w, float* z,
					lapack_int* ldz, lapack_int* isuppz, float* work,
					lapack_int* lwork, lapack_int* iwork, lapack_int* liwork,
					lapack_int *info );
void LAPACK_dstegr( char* jobz, char* range, lapack_int* n, double* d,
					double* e, double* vl, double* vu, lapack_int* il,
					lapack_int* iu, double* abstol, lapack_int* m, double* w,
					double* z, lapack_int* ldz, lapack_int* isuppz,
					double* work, lapack_int* lwork, lapack_int* iwork,
					lapack_int* liwork, lapack_int *info );
void LAPACK_cstegr( char* jobz, char* range, lapack_int* n, float* d, float* e,
					float* vl, float* vu, lapack_int* il, lapack_int* iu,
					float* abstol, lapack_int* m, float* w,
					lapack_complex_float* z, lapack_int* ldz,
					lapack_int* isuppz, float* work, lapack_int* lwork,
					lapack_int* iwork, lapack_int* liwork, lapack_int *info );
void LAPACK_zstegr( char* jobz, char* range, lapack_int* n, double* d,
					double* e, double* vl, double* vu, lapack_int* il,
					lapack_int* iu, double* abstol, lapack_int* m, double* w,
					lapack_complex_double* z, lapack_int* ldz,
					lapack_int* isuppz, double* work, lapack_int* lwork,
					lapack_int* iwork, lapack_int* liwork, lapack_int *info );
void LAPACK_spteqr( char* compz, lapack_int* n, float* d, float* e, float* z,
					lapack_int* ldz, float* work, lapack_int *info );
void LAPACK_dpteqr( char* compz, lapack_int* n, double* d, double* e, double* z,
					lapack_int* ldz, double* work, lapack_int *info );
void LAPACK_cpteqr( char* compz, lapack_int* n, float* d, float* e,
					lapack_complex_float* z, lapack_int* ldz, float* work,
					lapack_int *info );
void LAPACK_zpteqr( char* compz, lapack_int* n, double* d, double* e,
					lapack_complex_double* z, lapack_int* ldz, double* work,
					lapack_int *info );
void LAPACK_sstebz( char* range, char* order, lapack_int* n, float* vl,
					float* vu, lapack_int* il, lapack_int* iu, float* abstol,
					const float* d, const float* e, lapack_int* m,
					lapack_int* nsplit, float* w, lapack_int* iblock,
					lapack_int* isplit, float* work, lapack_int* iwork,
					lapack_int *info );
void LAPACK_dstebz( char* range, char* order, lapack_int* n, double* vl,
					double* vu, lapack_int* il, lapack_int* iu, double* abstol,
					const double* d, const double* e, lapack_int* m,
					lapack_int* nsplit, double* w, lapack_int* iblock,
					lapack_int* isplit, double* work, lapack_int* iwork,
					lapack_int *info );
void LAPACK_sstein( lapack_int* n, const float* d, const float* e,
					lapack_int* m, const float* w, const lapack_int* iblock,
					const lapack_int* isplit, float* z, lapack_int* ldz,
					float* work, lapack_int* iwork, lapack_int* ifailv,
					lapack_int *info );
void LAPACK_dstein( lapack_int* n, const double* d, const double* e,
					lapack_int* m, const double* w, const lapack_int* iblock,
					const lapack_int* isplit, double* z, lapack_int* ldz,
					double* work, lapack_int* iwork, lapack_int* ifailv,
					lapack_int *info );
void LAPACK_cstein( lapack_int* n, const float* d, const float* e,
					lapack_int* m, const float* w, const lapack_int* iblock,
					const lapack_int* isplit, lapack_complex_float* z,
					lapack_int* ldz, float* work, lapack_int* iwork,
					lapack_int* ifailv, lapack_int *info );
void LAPACK_zstein( lapack_int* n, const double* d, const double* e,
					lapack_int* m, const double* w, const lapack_int* iblock,
					const lapack_int* isplit, lapack_complex_double* z,
					lapack_int* ldz, double* work, lapack_int* iwork,
					lapack_int* ifailv, lapack_int *info );
void LAPACK_sdisna( char* job, lapack_int* m, lapack_int* n, const float* d,
					float* sep, lapack_int *info );
void LAPACK_ddisna( char* job, lapack_int* m, lapack_int* n, const double* d,
					double* sep, lapack_int *info );
void LAPACK_ssygst( lapack_int* itype, char* uplo, lapack_int* n, float* a,
					lapack_int* lda, const float* b, lapack_int* ldb,
					lapack_int *info );
void LAPACK_dsygst( lapack_int* itype, char* uplo, lapack_int* n, double* a,
					lapack_int* lda, const double* b, lapack_int* ldb,
					lapack_int *info );
void LAPACK_chegst( lapack_int* itype, char* uplo, lapack_int* n,
					lapack_complex_float* a, lapack_int* lda,
					const lapack_complex_float* b, lapack_int* ldb,
					lapack_int *info );
void LAPACK_zhegst( lapack_int* itype, char* uplo, lapack_int* n,
					lapack_complex_double* a, lapack_int* lda,
					const lapack_complex_double* b, lapack_int* ldb,
					lapack_int *info );
void LAPACK_sspgst( lapack_int* itype, char* uplo, lapack_int* n, float* ap,
					const float* bp, lapack_int *info );
void LAPACK_dspgst( lapack_int* itype, char* uplo, lapack_int* n, double* ap,
					const double* bp, lapack_int *info );
void LAPACK_chpgst( lapack_int* itype, char* uplo, lapack_int* n,
					lapack_complex_float* ap, const lapack_complex_float* bp,
					lapack_int *info );
void LAPACK_zhpgst( lapack_int* itype, char* uplo, lapack_int* n,
					lapack_complex_double* ap, const lapack_complex_double* bp,
					lapack_int *info );
void LAPACK_ssbgst( char* vect, char* uplo, lapack_int* n, lapack_int* ka,
					lapack_int* kb, float* ab, lapack_int* ldab,
					const float* bb, lapack_int* ldbb, float* x,
					lapack_int* ldx, float* work, lapack_int *info );
void LAPACK_dsbgst( char* vect, char* uplo, lapack_int* n, lapack_int* ka,
					lapack_int* kb, double* ab, lapack_int* ldab,
					const double* bb, lapack_int* ldbb, double* x,
					lapack_int* ldx, double* work, lapack_int *info );
void LAPACK_chbgst( char* vect, char* uplo, lapack_int* n, lapack_int* ka,
					lapack_int* kb, lapack_complex_float* ab, lapack_int* ldab,
					const lapack_complex_float* bb, lapack_int* ldbb,
					lapack_complex_float* x, lapack_int* ldx,
					lapack_complex_float* work, float* rwork,
					lapack_int *info );
void LAPACK_zhbgst( char* vect, char* uplo, lapack_int* n, lapack_int* ka,
					lapack_int* kb, lapack_complex_double* ab, lapack_int* ldab,
					const lapack_complex_double* bb, lapack_int* ldbb,
					lapack_complex_double* x, lapack_int* ldx,
					lapack_complex_double* work, double* rwork,
					lapack_int *info );
void LAPACK_spbstf( char* uplo, lapack_int* n, lapack_int* kb, float* bb,
					lapack_int* ldbb, lapack_int *info );
void LAPACK_dpbstf( char* uplo, lapack_int* n, lapack_int* kb, double* bb,
					lapack_int* ldbb, lapack_int *info );
void LAPACK_cpbstf( char* uplo, lapack_int* n, lapack_int* kb,
					lapack_complex_float* bb, lapack_int* ldbb,
					lapack_int *info );
void LAPACK_zpbstf( char* uplo, lapack_int* n, lapack_int* kb,
					lapack_complex_double* bb, lapack_int* ldbb,
					lapack_int *info );
void LAPACK_sgehrd( lapack_int* n, lapack_int* ilo, lapack_int* ihi, float* a,
					lapack_int* lda, float* tau, float* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_dgehrd( lapack_int* n, lapack_int* ilo, lapack_int* ihi, double* a,
					lapack_int* lda, double* tau, double* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_cgehrd( lapack_int* n, lapack_int* ilo, lapack_int* ihi,
					lapack_complex_float* a, lapack_int* lda,
					lapack_complex_float* tau, lapack_complex_float* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_zgehrd( lapack_int* n, lapack_int* ilo, lapack_int* ihi,
					lapack_complex_double* a, lapack_int* lda,
					lapack_complex_double* tau, lapack_complex_double* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_sorghr( lapack_int* n, lapack_int* ilo, lapack_int* ihi, float* a,
					lapack_int* lda, const float* tau, float* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_dorghr( lapack_int* n, lapack_int* ilo, lapack_int* ihi, double* a,
					lapack_int* lda, const double* tau, double* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_sormhr( char* side, char* trans, lapack_int* m, lapack_int* n,
					lapack_int* ilo, lapack_int* ihi, const float* a,
					lapack_int* lda, const float* tau, float* c,
					lapack_int* ldc, float* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_dormhr( char* side, char* trans, lapack_int* m, lapack_int* n,
					lapack_int* ilo, lapack_int* ihi, const double* a,
					lapack_int* lda, const double* tau, double* c,
					lapack_int* ldc, double* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_cunghr( lapack_int* n, lapack_int* ilo, lapack_int* ihi,
					lapack_complex_float* a, lapack_int* lda,
					const lapack_complex_float* tau, lapack_complex_float* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_zunghr( lapack_int* n, lapack_int* ilo, lapack_int* ihi,
					lapack_complex_double* a, lapack_int* lda,
					const lapack_complex_double* tau,
					lapack_complex_double* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_cunmhr( char* side, char* trans, lapack_int* m, lapack_int* n,
					lapack_int* ilo, lapack_int* ihi,
					const lapack_complex_float* a, lapack_int* lda,
					const lapack_complex_float* tau, lapack_complex_float* c,
					lapack_int* ldc, lapack_complex_float* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_zunmhr( char* side, char* trans, lapack_int* m, lapack_int* n,
					lapack_int* ilo, lapack_int* ihi,
					const lapack_complex_double* a, lapack_int* lda,
					const lapack_complex_double* tau, lapack_complex_double* c,
					lapack_int* ldc, lapack_complex_double* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_sgebal( char* job, lapack_int* n, float* a, lapack_int* lda,
					lapack_int* ilo, lapack_int* ihi, float* scale,
					lapack_int *info );
void LAPACK_dgebal( char* job, lapack_int* n, double* a, lapack_int* lda,
					lapack_int* ilo, lapack_int* ihi, double* scale,
					lapack_int *info );
void LAPACK_cgebal( char* job, lapack_int* n, lapack_complex_float* a,
					lapack_int* lda, lapack_int* ilo, lapack_int* ihi,
					float* scale, lapack_int *info );
void LAPACK_zgebal( char* job, lapack_int* n, lapack_complex_double* a,
					lapack_int* lda, lapack_int* ilo, lapack_int* ihi,
					double* scale, lapack_int *info );
void LAPACK_sgebak( char* job, char* side, lapack_int* n, lapack_int* ilo,
					lapack_int* ihi, const float* scale, lapack_int* m,
					float* v, lapack_int* ldv, lapack_int *info );
void LAPACK_dgebak( char* job, char* side, lapack_int* n, lapack_int* ilo,
					lapack_int* ihi, const double* scale, lapack_int* m,
					double* v, lapack_int* ldv, lapack_int *info );
void LAPACK_cgebak( char* job, char* side, lapack_int* n, lapack_int* ilo,
					lapack_int* ihi, const float* scale, lapack_int* m,
					lapack_complex_float* v, lapack_int* ldv,
					lapack_int *info );
void LAPACK_zgebak( char* job, char* side, lapack_int* n, lapack_int* ilo,
					lapack_int* ihi, const double* scale, lapack_int* m,
					lapack_complex_double* v, lapack_int* ldv,
					lapack_int *info );
void LAPACK_shseqr( char* job, char* compz, lapack_int* n, lapack_int* ilo,
					lapack_int* ihi, float* h, lapack_int* ldh, float* wr,
					float* wi, float* z, lapack_int* ldz, float* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_dhseqr( char* job, char* compz, lapack_int* n, lapack_int* ilo,
					lapack_int* ihi, double* h, lapack_int* ldh, double* wr,
					double* wi, double* z, lapack_int* ldz, double* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_chseqr( char* job, char* compz, lapack_int* n, lapack_int* ilo,
					lapack_int* ihi, lapack_complex_float* h, lapack_int* ldh,
					lapack_complex_float* w, lapack_complex_float* z,
					lapack_int* ldz, lapack_complex_float* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_zhseqr( char* job, char* compz, lapack_int* n, lapack_int* ilo,
					lapack_int* ihi, lapack_complex_double* h, lapack_int* ldh,
					lapack_complex_double* w, lapack_complex_double* z,
					lapack_int* ldz, lapack_complex_double* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_shsein( char* job, char* eigsrc, char* initv,
					lapack_logical* select, lapack_int* n, const float* h,
					lapack_int* ldh, float* wr, const float* wi, float* vl,
					lapack_int* ldvl, float* vr, lapack_int* ldvr,
					lapack_int* mm, lapack_int* m, float* work,
					lapack_int* ifaill, lapack_int* ifailr, lapack_int *info );
void LAPACK_dhsein( char* job, char* eigsrc, char* initv,
					lapack_logical* select, lapack_int* n, const double* h,
					lapack_int* ldh, double* wr, const double* wi, double* vl,
					lapack_int* ldvl, double* vr, lapack_int* ldvr,
					lapack_int* mm, lapack_int* m, double* work,
					lapack_int* ifaill, lapack_int* ifailr, lapack_int *info );
void LAPACK_chsein( char* job, char* eigsrc, char* initv,
					const lapack_logical* select, lapack_int* n,
					const lapack_complex_float* h, lapack_int* ldh,
					lapack_complex_float* w, lapack_complex_float* vl,
					lapack_int* ldvl, lapack_complex_float* vr,
					lapack_int* ldvr, lapack_int* mm, lapack_int* m,
					lapack_complex_float* work, float* rwork,
					lapack_int* ifaill, lapack_int* ifailr, lapack_int *info );
void LAPACK_zhsein( char* job, char* eigsrc, char* initv,
					const lapack_logical* select, lapack_int* n,
					const lapack_complex_double* h, lapack_int* ldh,
					lapack_complex_double* w, lapack_complex_double* vl,
					lapack_int* ldvl, lapack_complex_double* vr,
					lapack_int* ldvr, lapack_int* mm, lapack_int* m,
					lapack_complex_double* work, double* rwork,
					lapack_int* ifaill, lapack_int* ifailr, lapack_int *info );
void LAPACK_strevc( char* side, char* howmny, lapack_logical* select,
					lapack_int* n, const float* t, lapack_int* ldt, float* vl,
					lapack_int* ldvl, float* vr, lapack_int* ldvr,
					lapack_int* mm, lapack_int* m, float* work,
					lapack_int *info );
void LAPACK_dtrevc( char* side, char* howmny, lapack_logical* select,
					lapack_int* n, const double* t, lapack_int* ldt, double* vl,
					lapack_int* ldvl, double* vr, lapack_int* ldvr,
					lapack_int* mm, lapack_int* m, double* work,
					lapack_int *info );
void LAPACK_ctrevc( char* side, char* howmny, const lapack_logical* select,
					lapack_int* n, lapack_complex_float* t, lapack_int* ldt,
					lapack_complex_float* vl, lapack_int* ldvl,
					lapack_complex_float* vr, lapack_int* ldvr, lapack_int* mm,
					lapack_int* m, lapack_complex_float* work, float* rwork,
					lapack_int *info );
void LAPACK_ztrevc( char* side, char* howmny, const lapack_logical* select,
					lapack_int* n, lapack_complex_double* t, lapack_int* ldt,
					lapack_complex_double* vl, lapack_int* ldvl,
					lapack_complex_double* vr, lapack_int* ldvr, lapack_int* mm,
					lapack_int* m, lapack_complex_double* work, double* rwork,
					lapack_int *info );
void LAPACK_strsna( char* job, char* howmny, const lapack_logical* select,
					lapack_int* n, const float* t, lapack_int* ldt,
					const float* vl, lapack_int* ldvl, const float* vr,
					lapack_int* ldvr, float* s, float* sep, lapack_int* mm,
					lapack_int* m, float* work, lapack_int* ldwork,
					lapack_int* iwork, lapack_int *info );
void LAPACK_dtrsna( char* job, char* howmny, const lapack_logical* select,
					lapack_int* n, const double* t, lapack_int* ldt,
					const double* vl, lapack_int* ldvl, const double* vr,
					lapack_int* ldvr, double* s, double* sep, lapack_int* mm,
					lapack_int* m, double* work, lapack_int* ldwork,
					lapack_int* iwork, lapack_int *info );
void LAPACK_ctrsna( char* job, char* howmny, const lapack_logical* select,
					lapack_int* n, const lapack_complex_float* t,
					lapack_int* ldt, const lapack_complex_float* vl,
					lapack_int* ldvl, const lapack_complex_float* vr,
					lapack_int* ldvr, float* s, float* sep, lapack_int* mm,
					lapack_int* m, lapack_complex_float* work,
					lapack_int* ldwork, float* rwork, lapack_int *info );
void LAPACK_ztrsna( char* job, char* howmny, const lapack_logical* select,
					lapack_int* n, const lapack_complex_double* t,
					lapack_int* ldt, const lapack_complex_double* vl,
					lapack_int* ldvl, const lapack_complex_double* vr,
					lapack_int* ldvr, double* s, double* sep, lapack_int* mm,
					lapack_int* m, lapack_complex_double* work,
					lapack_int* ldwork, double* rwork, lapack_int *info );
void LAPACK_strexc( char* compq, lapack_int* n, float* t, lapack_int* ldt,
					float* q, lapack_int* ldq, lapack_int* ifst,
					lapack_int* ilst, float* work, lapack_int *info );
void LAPACK_dtrexc( char* compq, lapack_int* n, double* t, lapack_int* ldt,
					double* q, lapack_int* ldq, lapack_int* ifst,
					lapack_int* ilst, double* work, lapack_int *info );
void LAPACK_ctrexc( char* compq, lapack_int* n, lapack_complex_float* t,
					lapack_int* ldt, lapack_complex_float* q, lapack_int* ldq,
					lapack_int* ifst, lapack_int* ilst, lapack_int *info );
void LAPACK_ztrexc( char* compq, lapack_int* n, lapack_complex_double* t,
					lapack_int* ldt, lapack_complex_double* q, lapack_int* ldq,
					lapack_int* ifst, lapack_int* ilst, lapack_int *info );
void LAPACK_strsen( char* job, char* compq, const lapack_logical* select,
					lapack_int* n, float* t, lapack_int* ldt, float* q,
					lapack_int* ldq, float* wr, float* wi, lapack_int* m,
					float* s, float* sep, float* work, lapack_int* lwork,
					lapack_int* iwork, lapack_int* liwork, lapack_int *info );
void LAPACK_dtrsen( char* job, char* compq, const lapack_logical* select,
					lapack_int* n, double* t, lapack_int* ldt, double* q,
					lapack_int* ldq, double* wr, double* wi, lapack_int* m,
					double* s, double* sep, double* work, lapack_int* lwork,
					lapack_int* iwork, lapack_int* liwork, lapack_int *info );
void LAPACK_ctrsen( char* job, char* compq, const lapack_logical* select,
					lapack_int* n, lapack_complex_float* t, lapack_int* ldt,
					lapack_complex_float* q, lapack_int* ldq,
					lapack_complex_float* w, lapack_int* m, float* s,
					float* sep, lapack_complex_float* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_ztrsen( char* job, char* compq, const lapack_logical* select,
					lapack_int* n, lapack_complex_double* t, lapack_int* ldt,
					lapack_complex_double* q, lapack_int* ldq,
					lapack_complex_double* w, lapack_int* m, double* s,
					double* sep, lapack_complex_double* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_strsyl( char* trana, char* tranb, lapack_int* isgn, lapack_int* m,
					lapack_int* n, const float* a, lapack_int* lda,
					const float* b, lapack_int* ldb, float* c, lapack_int* ldc,
					float* scale, lapack_int *info );
void LAPACK_dtrsyl( char* trana, char* tranb, lapack_int* isgn, lapack_int* m,
					lapack_int* n, const double* a, lapack_int* lda,
					const double* b, lapack_int* ldb, double* c,
					lapack_int* ldc, double* scale, lapack_int *info );
void LAPACK_ctrsyl( char* trana, char* tranb, lapack_int* isgn, lapack_int* m,
					lapack_int* n, const lapack_complex_float* a,
					lapack_int* lda, const lapack_complex_float* b,
					lapack_int* ldb, lapack_complex_float* c, lapack_int* ldc,
					float* scale, lapack_int *info );
void LAPACK_ztrsyl( char* trana, char* tranb, lapack_int* isgn, lapack_int* m,
					lapack_int* n, const lapack_complex_double* a,
					lapack_int* lda, const lapack_complex_double* b,
					lapack_int* ldb, lapack_complex_double* c, lapack_int* ldc,
					double* scale, lapack_int *info );
void LAPACK_sgghrd( char* compq, char* compz, lapack_int* n, lapack_int* ilo,
					lapack_int* ihi, float* a, lapack_int* lda, float* b,
					lapack_int* ldb, float* q, lapack_int* ldq, float* z,
					lapack_int* ldz, lapack_int *info );
void LAPACK_dgghrd( char* compq, char* compz, lapack_int* n, lapack_int* ilo,
					lapack_int* ihi, double* a, lapack_int* lda, double* b,
					lapack_int* ldb, double* q, lapack_int* ldq, double* z,
					lapack_int* ldz, lapack_int *info );
void LAPACK_cgghrd( char* compq, char* compz, lapack_int* n, lapack_int* ilo,
					lapack_int* ihi, lapack_complex_float* a, lapack_int* lda,
					lapack_complex_float* b, lapack_int* ldb,
					lapack_complex_float* q, lapack_int* ldq,
					lapack_complex_float* z, lapack_int* ldz,
					lapack_int *info );
void LAPACK_zgghrd( char* compq, char* compz, lapack_int* n, lapack_int* ilo,
					lapack_int* ihi, lapack_complex_double* a, lapack_int* lda,
					lapack_complex_double* b, lapack_int* ldb,
					lapack_complex_double* q, lapack_int* ldq,
					lapack_complex_double* z, lapack_int* ldz,
					lapack_int *info );
void LAPACK_sggbal( char* job, lapack_int* n, float* a, lapack_int* lda,
					float* b, lapack_int* ldb, lapack_int* ilo, lapack_int* ihi,
					float* lscale, float* rscale, float* work,
					lapack_int *info );
void LAPACK_dggbal( char* job, lapack_int* n, double* a, lapack_int* lda,
					double* b, lapack_int* ldb, lapack_int* ilo,
					lapack_int* ihi, double* lscale, double* rscale,
					double* work, lapack_int *info );
void LAPACK_cggbal( char* job, lapack_int* n, lapack_complex_float* a,
					lapack_int* lda, lapack_complex_float* b, lapack_int* ldb,
					lapack_int* ilo, lapack_int* ihi, float* lscale,
					float* rscale, float* work, lapack_int *info );
void LAPACK_zggbal( char* job, lapack_int* n, lapack_complex_double* a,
					lapack_int* lda, lapack_complex_double* b, lapack_int* ldb,
					lapack_int* ilo, lapack_int* ihi, double* lscale,
					double* rscale, double* work, lapack_int *info );
void LAPACK_sggbak( char* job, char* side, lapack_int* n, lapack_int* ilo,
					lapack_int* ihi, const float* lscale, const float* rscale,
					lapack_int* m, float* v, lapack_int* ldv,
					lapack_int *info );
void LAPACK_dggbak( char* job, char* side, lapack_int* n, lapack_int* ilo,
					lapack_int* ihi, const double* lscale, const double* rscale,
					lapack_int* m, double* v, lapack_int* ldv,
					lapack_int *info );
void LAPACK_cggbak( char* job, char* side, lapack_int* n, lapack_int* ilo,
					lapack_int* ihi, const float* lscale, const float* rscale,
					lapack_int* m, lapack_complex_float* v, lapack_int* ldv,
					lapack_int *info );
void LAPACK_zggbak( char* job, char* side, lapack_int* n, lapack_int* ilo,
					lapack_int* ihi, const double* lscale, const double* rscale,
					lapack_int* m, lapack_complex_double* v, lapack_int* ldv,
					lapack_int *info );
void LAPACK_shgeqz( char* job, char* compq, char* compz, lapack_int* n,
					lapack_int* ilo, lapack_int* ihi, float* h, lapack_int* ldh,
					float* t, lapack_int* ldt, float* alphar, float* alphai,
					float* beta, float* q, lapack_int* ldq, float* z,
					lapack_int* ldz, float* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_dhgeqz( char* job, char* compq, char* compz, lapack_int* n,
					lapack_int* ilo, lapack_int* ihi, double* h,
					lapack_int* ldh, double* t, lapack_int* ldt, double* alphar,
					double* alphai, double* beta, double* q, lapack_int* ldq,
					double* z, lapack_int* ldz, double* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_chgeqz( char* job, char* compq, char* compz, lapack_int* n,
					lapack_int* ilo, lapack_int* ihi, lapack_complex_float* h,
					lapack_int* ldh, lapack_complex_float* t, lapack_int* ldt,
					lapack_complex_float* alpha, lapack_complex_float* beta,
					lapack_complex_float* q, lapack_int* ldq,
					lapack_complex_float* z, lapack_int* ldz,
					lapack_complex_float* work, lapack_int* lwork, float* rwork,
					lapack_int *info );
void LAPACK_zhgeqz( char* job, char* compq, char* compz, lapack_int* n,
					lapack_int* ilo, lapack_int* ihi, lapack_complex_double* h,
					lapack_int* ldh, lapack_complex_double* t, lapack_int* ldt,
					lapack_complex_double* alpha, lapack_complex_double* beta,
					lapack_complex_double* q, lapack_int* ldq,
					lapack_complex_double* z, lapack_int* ldz,
					lapack_complex_double* work, lapack_int* lwork,
					double* rwork, lapack_int *info );
void LAPACK_stgevc( char* side, char* howmny, const lapack_logical* select,
					lapack_int* n, const float* s, lapack_int* lds,
					const float* p, lapack_int* ldp, float* vl,
					lapack_int* ldvl, float* vr, lapack_int* ldvr,
					lapack_int* mm, lapack_int* m, float* work,
					lapack_int *info );
void LAPACK_dtgevc( char* side, char* howmny, const lapack_logical* select,
					lapack_int* n, const double* s, lapack_int* lds,
					const double* p, lapack_int* ldp, double* vl,
					lapack_int* ldvl, double* vr, lapack_int* ldvr,
					lapack_int* mm, lapack_int* m, double* work,
					lapack_int *info );
void LAPACK_ctgevc( char* side, char* howmny, const lapack_logical* select,
					lapack_int* n, const lapack_complex_float* s,
					lapack_int* lds, const lapack_complex_float* p,
					lapack_int* ldp, lapack_complex_float* vl, lapack_int* ldvl,
					lapack_complex_float* vr, lapack_int* ldvr, lapack_int* mm,
					lapack_int* m, lapack_complex_float* work, float* rwork,
					lapack_int *info );
void LAPACK_ztgevc( char* side, char* howmny, const lapack_logical* select,
					lapack_int* n, const lapack_complex_double* s,
					lapack_int* lds, const lapack_complex_double* p,
					lapack_int* ldp, lapack_complex_double* vl,
					lapack_int* ldvl, lapack_complex_double* vr,
					lapack_int* ldvr, lapack_int* mm, lapack_int* m,
					lapack_complex_double* work, double* rwork,
					lapack_int *info );
void LAPACK_stgexc( lapack_logical* wantq, lapack_logical* wantz, lapack_int* n,
					float* a, lapack_int* lda, float* b, lapack_int* ldb,
					float* q, lapack_int* ldq, float* z, lapack_int* ldz,
					lapack_int* ifst, lapack_int* ilst, float* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_dtgexc( lapack_logical* wantq, lapack_logical* wantz, lapack_int* n,
					double* a, lapack_int* lda, double* b, lapack_int* ldb,
					double* q, lapack_int* ldq, double* z, lapack_int* ldz,
					lapack_int* ifst, lapack_int* ilst, double* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_ctgexc( lapack_logical* wantq, lapack_logical* wantz, lapack_int* n,
					lapack_complex_float* a, lapack_int* lda,
					lapack_complex_float* b, lapack_int* ldb,
					lapack_complex_float* q, lapack_int* ldq,
					lapack_complex_float* z, lapack_int* ldz, lapack_int* ifst,
					lapack_int* ilst, lapack_int *info );
void LAPACK_ztgexc( lapack_logical* wantq, lapack_logical* wantz, lapack_int* n,
					lapack_complex_double* a, lapack_int* lda,
					lapack_complex_double* b, lapack_int* ldb,
					lapack_complex_double* q, lapack_int* ldq,
					lapack_complex_double* z, lapack_int* ldz, lapack_int* ifst,
					lapack_int* ilst, lapack_int *info );
void LAPACK_stgsen( lapack_int* ijob, lapack_logical* wantq,
					lapack_logical* wantz, const lapack_logical* select,
					lapack_int* n, float* a, lapack_int* lda, float* b,
					lapack_int* ldb, float* alphar, float* alphai, float* beta,
					float* q, lapack_int* ldq, float* z, lapack_int* ldz,
					lapack_int* m, float* pl, float* pr, float* dif,
					float* work, lapack_int* lwork, lapack_int* iwork,
					lapack_int* liwork, lapack_int *info );
void LAPACK_dtgsen( lapack_int* ijob, lapack_logical* wantq,
					lapack_logical* wantz, const lapack_logical* select,
					lapack_int* n, double* a, lapack_int* lda, double* b,
					lapack_int* ldb, double* alphar, double* alphai,
					double* beta, double* q, lapack_int* ldq, double* z,
					lapack_int* ldz, lapack_int* m, double* pl, double* pr,
					double* dif, double* work, lapack_int* lwork,
					lapack_int* iwork, lapack_int* liwork, lapack_int *info );
void LAPACK_ctgsen( lapack_int* ijob, lapack_logical* wantq,
					lapack_logical* wantz, const lapack_logical* select,
					lapack_int* n, lapack_complex_float* a, lapack_int* lda,
					lapack_complex_float* b, lapack_int* ldb,
					lapack_complex_float* alpha, lapack_complex_float* beta,
					lapack_complex_float* q, lapack_int* ldq,
					lapack_complex_float* z, lapack_int* ldz, lapack_int* m,
					float* pl, float* pr, float* dif,
					lapack_complex_float* work, lapack_int* lwork,
					lapack_int* iwork, lapack_int* liwork, lapack_int *info );
void LAPACK_ztgsen( lapack_int* ijob, lapack_logical* wantq,
					lapack_logical* wantz, const lapack_logical* select,
					lapack_int* n, lapack_complex_double* a, lapack_int* lda,
					lapack_complex_double* b, lapack_int* ldb,
					lapack_complex_double* alpha, lapack_complex_double* beta,
					lapack_complex_double* q, lapack_int* ldq,
					lapack_complex_double* z, lapack_int* ldz, lapack_int* m,
					double* pl, double* pr, double* dif,
					lapack_complex_double* work, lapack_int* lwork,
					lapack_int* iwork, lapack_int* liwork, lapack_int *info );
void LAPACK_stgsyl( char* trans, lapack_int* ijob, lapack_int* m, lapack_int* n,
					const float* a, lapack_int* lda, const float* b,
					lapack_int* ldb, float* c, lapack_int* ldc, const float* d,
					lapack_int* ldd, const float* e, lapack_int* lde, float* f,
					lapack_int* ldf, float* scale, float* dif, float* work,
					lapack_int* lwork, lapack_int* iwork, lapack_int *info );
void LAPACK_dtgsyl( char* trans, lapack_int* ijob, lapack_int* m, lapack_int* n,
					const double* a, lapack_int* lda, const double* b,
					lapack_int* ldb, double* c, lapack_int* ldc,
					const double* d, lapack_int* ldd, const double* e,
					lapack_int* lde, double* f, lapack_int* ldf, double* scale,
					double* dif, double* work, lapack_int* lwork,
					lapack_int* iwork, lapack_int *info );
void LAPACK_ctgsyl( char* trans, lapack_int* ijob, lapack_int* m, lapack_int* n,
					const lapack_complex_float* a, lapack_int* lda,
					const lapack_complex_float* b, lapack_int* ldb,
					lapack_complex_float* c, lapack_int* ldc,
					const lapack_complex_float* d, lapack_int* ldd,
					const lapack_complex_float* e, lapack_int* lde,
					lapack_complex_float* f, lapack_int* ldf, float* scale,
					float* dif, lapack_complex_float* work, lapack_int* lwork,
					lapack_int* iwork, lapack_int *info );
void LAPACK_ztgsyl( char* trans, lapack_int* ijob, lapack_int* m, lapack_int* n,
					const lapack_complex_double* a, lapack_int* lda,
					const lapack_complex_double* b, lapack_int* ldb,
					lapack_complex_double* c, lapack_int* ldc,
					const lapack_complex_double* d, lapack_int* ldd,
					const lapack_complex_double* e, lapack_int* lde,
					lapack_complex_double* f, lapack_int* ldf, double* scale,
					double* dif, lapack_complex_double* work, lapack_int* lwork,
					lapack_int* iwork, lapack_int *info );
void LAPACK_stgsna( char* job, char* howmny, const lapack_logical* select,
					lapack_int* n, const float* a, lapack_int* lda,
					const float* b, lapack_int* ldb, const float* vl,
					lapack_int* ldvl, const float* vr, lapack_int* ldvr,
					float* s, float* dif, lapack_int* mm, lapack_int* m,
					float* work, lapack_int* lwork, lapack_int* iwork,
					lapack_int *info );
void LAPACK_dtgsna( char* job, char* howmny, const lapack_logical* select,
					lapack_int* n, const double* a, lapack_int* lda,
					const double* b, lapack_int* ldb, const double* vl,
					lapack_int* ldvl, const double* vr, lapack_int* ldvr,
					double* s, double* dif, lapack_int* mm, lapack_int* m,
					double* work, lapack_int* lwork, lapack_int* iwork,
					lapack_int *info );
void LAPACK_ctgsna( char* job, char* howmny, const lapack_logical* select,
					lapack_int* n, const lapack_complex_float* a,
					lapack_int* lda, const lapack_complex_float* b,
					lapack_int* ldb, const lapack_complex_float* vl,
					lapack_int* ldvl, const lapack_complex_float* vr,
					lapack_int* ldvr, float* s, float* dif, lapack_int* mm,
					lapack_int* m, lapack_complex_float* work,
					lapack_int* lwork, lapack_int* iwork, lapack_int *info );
void LAPACK_ztgsna( char* job, char* howmny, const lapack_logical* select,
					lapack_int* n, const lapack_complex_double* a,
					lapack_int* lda, const lapack_complex_double* b,
					lapack_int* ldb, const lapack_complex_double* vl,
					lapack_int* ldvl, const lapack_complex_double* vr,
					lapack_int* ldvr, double* s, double* dif, lapack_int* mm,
					lapack_int* m, lapack_complex_double* work,
					lapack_int* lwork, lapack_int* iwork, lapack_int *info );
void LAPACK_sggsvp( char* jobu, char* jobv, char* jobq, lapack_int* m,
					lapack_int* p, lapack_int* n, float* a, lapack_int* lda,
					float* b, lapack_int* ldb, float* tola, float* tolb,
					lapack_int* k, lapack_int* l, float* u, lapack_int* ldu,
					float* v, lapack_int* ldv, float* q, lapack_int* ldq,
					lapack_int* iwork, float* tau, float* work,
					lapack_int *info );
void LAPACK_dggsvp( char* jobu, char* jobv, char* jobq, lapack_int* m,
					lapack_int* p, lapack_int* n, double* a, lapack_int* lda,
					double* b, lapack_int* ldb, double* tola, double* tolb,
					lapack_int* k, lapack_int* l, double* u, lapack_int* ldu,
					double* v, lapack_int* ldv, double* q, lapack_int* ldq,
					lapack_int* iwork, double* tau, double* work,
					lapack_int *info );
void LAPACK_cggsvp( char* jobu, char* jobv, char* jobq, lapack_int* m,
					lapack_int* p, lapack_int* n, lapack_complex_float* a,
					lapack_int* lda, lapack_complex_float* b, lapack_int* ldb,
					float* tola, float* tolb, lapack_int* k, lapack_int* l,
					lapack_complex_float* u, lapack_int* ldu,
					lapack_complex_float* v, lapack_int* ldv,
					lapack_complex_float* q, lapack_int* ldq, lapack_int* iwork,
					float* rwork, lapack_complex_float* tau,
					lapack_complex_float* work, lapack_int *info );
void LAPACK_zggsvp( char* jobu, char* jobv, char* jobq, lapack_int* m,
					lapack_int* p, lapack_int* n, lapack_complex_double* a,
					lapack_int* lda, lapack_complex_double* b, lapack_int* ldb,
					double* tola, double* tolb, lapack_int* k, lapack_int* l,
					lapack_complex_double* u, lapack_int* ldu,
					lapack_complex_double* v, lapack_int* ldv,
					lapack_complex_double* q, lapack_int* ldq,
					lapack_int* iwork, double* rwork,
					lapack_complex_double* tau, lapack_complex_double* work,
					lapack_int *info );
void LAPACK_stgsja( char* jobu, char* jobv, char* jobq, lapack_int* m,
					lapack_int* p, lapack_int* n, lapack_int* k, lapack_int* l,
					float* a, lapack_int* lda, float* b, lapack_int* ldb,
					float* tola, float* tolb, float* alpha, float* beta,
					float* u, lapack_int* ldu, float* v, lapack_int* ldv,
					float* q, lapack_int* ldq, float* work, lapack_int* ncycle,
					lapack_int *info );
void LAPACK_dtgsja( char* jobu, char* jobv, char* jobq, lapack_int* m,
					lapack_int* p, lapack_int* n, lapack_int* k, lapack_int* l,
					double* a, lapack_int* lda, double* b, lapack_int* ldb,
					double* tola, double* tolb, double* alpha, double* beta,
					double* u, lapack_int* ldu, double* v, lapack_int* ldv,
					double* q, lapack_int* ldq, double* work,
					lapack_int* ncycle, lapack_int *info );
void LAPACK_ctgsja( char* jobu, char* jobv, char* jobq, lapack_int* m,
					lapack_int* p, lapack_int* n, lapack_int* k, lapack_int* l,
					lapack_complex_float* a, lapack_int* lda,
					lapack_complex_float* b, lapack_int* ldb, float* tola,
					float* tolb, float* alpha, float* beta,
					lapack_complex_float* u, lapack_int* ldu,
					lapack_complex_float* v, lapack_int* ldv,
					lapack_complex_float* q, lapack_int* ldq,
					lapack_complex_float* work, lapack_int* ncycle,
					lapack_int *info );
void LAPACK_ztgsja( char* jobu, char* jobv, char* jobq, lapack_int* m,
					lapack_int* p, lapack_int* n, lapack_int* k, lapack_int* l,
					lapack_complex_double* a, lapack_int* lda,
					lapack_complex_double* b, lapack_int* ldb, double* tola,
					double* tolb, double* alpha, double* beta,
					lapack_complex_double* u, lapack_int* ldu,
					lapack_complex_double* v, lapack_int* ldv,
					lapack_complex_double* q, lapack_int* ldq,
					lapack_complex_double* work, lapack_int* ncycle,
					lapack_int *info );
void LAPACK_sgels( char* trans, lapack_int* m, lapack_int* n, lapack_int* nrhs,
				float* a, lapack_int* lda, float* b, lapack_int* ldb,
				float* work, lapack_int* lwork, lapack_int *info );
void LAPACK_dgels( char* trans, lapack_int* m, lapack_int* n, lapack_int* nrhs,
				double* a, lapack_int* lda, double* b, lapack_int* ldb,
				double* work, lapack_int* lwork, lapack_int *info );
void LAPACK_cgels( char* trans, lapack_int* m, lapack_int* n, lapack_int* nrhs,
				lapack_complex_float* a, lapack_int* lda,
				lapack_complex_float* b, lapack_int* ldb,
				lapack_complex_float* work, lapack_int* lwork,
				lapack_int *info );
void LAPACK_zgels( char* trans, lapack_int* m, lapack_int* n, lapack_int* nrhs,
				lapack_complex_double* a, lapack_int* lda,
				lapack_complex_double* b, lapack_int* ldb,
				lapack_complex_double* work, lapack_int* lwork,
				lapack_int *info );
void LAPACK_sgelsy( lapack_int* m, lapack_int* n, lapack_int* nrhs, float* a,
					lapack_int* lda, float* b, lapack_int* ldb,
					lapack_int* jpvt, float* rcond, lapack_int* rank,
					float* work, lapack_int* lwork, lapack_int *info );
void LAPACK_dgelsy( lapack_int* m, lapack_int* n, lapack_int* nrhs, double* a,
					lapack_int* lda, double* b, lapack_int* ldb,
					lapack_int* jpvt, double* rcond, lapack_int* rank,
					double* work, lapack_int* lwork, lapack_int *info );
void LAPACK_cgelsy( lapack_int* m, lapack_int* n, lapack_int* nrhs,
					lapack_complex_float* a, lapack_int* lda,
					lapack_complex_float* b, lapack_int* ldb, lapack_int* jpvt,
					float* rcond, lapack_int* rank, lapack_complex_float* work,
					lapack_int* lwork, float* rwork, lapack_int *info );
void LAPACK_zgelsy( lapack_int* m, lapack_int* n, lapack_int* nrhs,
					lapack_complex_double* a, lapack_int* lda,
					lapack_complex_double* b, lapack_int* ldb, lapack_int* jpvt,
					double* rcond, lapack_int* rank,
					lapack_complex_double* work, lapack_int* lwork,
					double* rwork, lapack_int *info );
void LAPACK_sgelss( lapack_int* m, lapack_int* n, lapack_int* nrhs, float* a,
					lapack_int* lda, float* b, lapack_int* ldb, float* s,
					float* rcond, lapack_int* rank, float* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_dgelss( lapack_int* m, lapack_int* n, lapack_int* nrhs, double* a,
					lapack_int* lda, double* b, lapack_int* ldb, double* s,
					double* rcond, lapack_int* rank, double* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_cgelss( lapack_int* m, lapack_int* n, lapack_int* nrhs,
					lapack_complex_float* a, lapack_int* lda,
					lapack_complex_float* b, lapack_int* ldb, float* s,
					float* rcond, lapack_int* rank, lapack_complex_float* work,
					lapack_int* lwork, float* rwork, lapack_int *info );
void LAPACK_zgelss( lapack_int* m, lapack_int* n, lapack_int* nrhs,
					lapack_complex_double* a, lapack_int* lda,
					lapack_complex_double* b, lapack_int* ldb, double* s,
					double* rcond, lapack_int* rank,
					lapack_complex_double* work, lapack_int* lwork,
					double* rwork, lapack_int *info );
void LAPACK_sgelsd( lapack_int* m, lapack_int* n, lapack_int* nrhs, float* a,
					lapack_int* lda, float* b, lapack_int* ldb, float* s,
					float* rcond, lapack_int* rank, float* work,
					lapack_int* lwork, lapack_int* iwork, lapack_int *info );
void LAPACK_dgelsd( lapack_int* m, lapack_int* n, lapack_int* nrhs, double* a,
					lapack_int* lda, double* b, lapack_int* ldb, double* s,
					double* rcond, lapack_int* rank, double* work,
					lapack_int* lwork, lapack_int* iwork, lapack_int *info );
void LAPACK_cgelsd( lapack_int* m, lapack_int* n, lapack_int* nrhs,
					lapack_complex_float* a, lapack_int* lda,
					lapack_complex_float* b, lapack_int* ldb, float* s,
					float* rcond, lapack_int* rank, lapack_complex_float* work,
					lapack_int* lwork, float* rwork, lapack_int* iwork,
					lapack_int *info );
void LAPACK_zgelsd( lapack_int* m, lapack_int* n, lapack_int* nrhs,
					lapack_complex_double* a, lapack_int* lda,
					lapack_complex_double* b, lapack_int* ldb, double* s,
					double* rcond, lapack_int* rank,
					lapack_complex_double* work, lapack_int* lwork,
					double* rwork, lapack_int* iwork, lapack_int *info );
void LAPACK_sgglse( lapack_int* m, lapack_int* n, lapack_int* p, float* a,
					lapack_int* lda, float* b, lapack_int* ldb, float* c,
					float* d, float* x, float* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_dgglse( lapack_int* m, lapack_int* n, lapack_int* p, double* a,
					lapack_int* lda, double* b, lapack_int* ldb, double* c,
					double* d, double* x, double* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_cgglse( lapack_int* m, lapack_int* n, lapack_int* p,
					lapack_complex_float* a, lapack_int* lda,
					lapack_complex_float* b, lapack_int* ldb,
					lapack_complex_float* c, lapack_complex_float* d,
					lapack_complex_float* x, lapack_complex_float* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_zgglse( lapack_int* m, lapack_int* n, lapack_int* p,
					lapack_complex_double* a, lapack_int* lda,
					lapack_complex_double* b, lapack_int* ldb,
					lapack_complex_double* c, lapack_complex_double* d,
					lapack_complex_double* x, lapack_complex_double* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_sggglm( lapack_int* n, lapack_int* m, lapack_int* p, float* a,
					lapack_int* lda, float* b, lapack_int* ldb, float* d,
					float* x, float* y, float* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_dggglm( lapack_int* n, lapack_int* m, lapack_int* p, double* a,
					lapack_int* lda, double* b, lapack_int* ldb, double* d,
					double* x, double* y, double* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_cggglm( lapack_int* n, lapack_int* m, lapack_int* p,
					lapack_complex_float* a, lapack_int* lda,
					lapack_complex_float* b, lapack_int* ldb,
					lapack_complex_float* d, lapack_complex_float* x,
					lapack_complex_float* y, lapack_complex_float* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_zggglm( lapack_int* n, lapack_int* m, lapack_int* p,
					lapack_complex_double* a, lapack_int* lda,
					lapack_complex_double* b, lapack_int* ldb,
					lapack_complex_double* d, lapack_complex_double* x,
					lapack_complex_double* y, lapack_complex_double* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_ssyev( char* jobz, char* uplo, lapack_int* n, float* a,
				lapack_int* lda, float* w, float* work, lapack_int* lwork,
				lapack_int *info );
void LAPACK_dsyev( char* jobz, char* uplo, lapack_int* n, double* a,
				lapack_int* lda, double* w, double* work, lapack_int* lwork,
				lapack_int *info );
void LAPACK_cheev( char* jobz, char* uplo, lapack_int* n,
				lapack_complex_float* a, lapack_int* lda, float* w,
				lapack_complex_float* work, lapack_int* lwork, float* rwork,
				lapack_int *info );
void LAPACK_zheev( char* jobz, char* uplo, lapack_int* n,
				lapack_complex_double* a, lapack_int* lda, double* w,
				lapack_complex_double* work, lapack_int* lwork,
				double* rwork, lapack_int *info );
void LAPACK_ssyevd( char* jobz, char* uplo, lapack_int* n, float* a,
					lapack_int* lda, float* w, float* work, lapack_int* lwork,
					lapack_int* iwork, lapack_int* liwork, lapack_int *info );
void LAPACK_dsyevd( char* jobz, char* uplo, lapack_int* n, double* a,
					lapack_int* lda, double* w, double* work, lapack_int* lwork,
					lapack_int* iwork, lapack_int* liwork, lapack_int *info );
void LAPACK_cheevd( char* jobz, char* uplo, lapack_int* n,
					lapack_complex_float* a, lapack_int* lda, float* w,
					lapack_complex_float* work, lapack_int* lwork, float* rwork,
					lapack_int* lrwork, lapack_int* iwork, lapack_int* liwork,
					lapack_int *info );
void LAPACK_zheevd( char* jobz, char* uplo, lapack_int* n,
					lapack_complex_double* a, lapack_int* lda, double* w,
					lapack_complex_double* work, lapack_int* lwork,
					double* rwork, lapack_int* lrwork, lapack_int* iwork,
					lapack_int* liwork, lapack_int *info );
void LAPACK_ssyevx( char* jobz, char* range, char* uplo, lapack_int* n,
					float* a, lapack_int* lda, float* vl, float* vu,
					lapack_int* il, lapack_int* iu, float* abstol,
					lapack_int* m, float* w, float* z, lapack_int* ldz,
					float* work, lapack_int* lwork, lapack_int* iwork,
					lapack_int* ifail, lapack_int *info );
void LAPACK_dsyevx( char* jobz, char* range, char* uplo, lapack_int* n,
					double* a, lapack_int* lda, double* vl, double* vu,
					lapack_int* il, lapack_int* iu, double* abstol,
					lapack_int* m, double* w, double* z, lapack_int* ldz,
					double* work, lapack_int* lwork, lapack_int* iwork,
					lapack_int* ifail, lapack_int *info );
void LAPACK_cheevx( char* jobz, char* range, char* uplo, lapack_int* n,
					lapack_complex_float* a, lapack_int* lda, float* vl,
					float* vu, lapack_int* il, lapack_int* iu, float* abstol,
					lapack_int* m, float* w, lapack_complex_float* z,
					lapack_int* ldz, lapack_complex_float* work,
					lapack_int* lwork, float* rwork, lapack_int* iwork,
					lapack_int* ifail, lapack_int *info );
void LAPACK_zheevx( char* jobz, char* range, char* uplo, lapack_int* n,
					lapack_complex_double* a, lapack_int* lda, double* vl,
					double* vu, lapack_int* il, lapack_int* iu, double* abstol,
					lapack_int* m, double* w, lapack_complex_double* z,
					lapack_int* ldz, lapack_complex_double* work,
					lapack_int* lwork, double* rwork, lapack_int* iwork,
					lapack_int* ifail, lapack_int *info );
void LAPACK_ssyevr( char* jobz, char* range, char* uplo, lapack_int* n,
					float* a, lapack_int* lda, float* vl, float* vu,
					lapack_int* il, lapack_int* iu, float* abstol,
					lapack_int* m, float* w, float* z, lapack_int* ldz,
					lapack_int* isuppz, float* work, lapack_int* lwork,
					lapack_int* iwork, lapack_int* liwork, lapack_int *info );
void LAPACK_dsyevr( char* jobz, char* range, char* uplo, lapack_int* n,
					double* a, lapack_int* lda, double* vl, double* vu,
					lapack_int* il, lapack_int* iu, double* abstol,
					lapack_int* m, double* w, double* z, lapack_int* ldz,
					lapack_int* isuppz, double* work, lapack_int* lwork,
					lapack_int* iwork, lapack_int* liwork, lapack_int *info );
void LAPACK_cheevr( char* jobz, char* range, char* uplo, lapack_int* n,
					lapack_complex_float* a, lapack_int* lda, float* vl,
					float* vu, lapack_int* il, lapack_int* iu, float* abstol,
					lapack_int* m, float* w, lapack_complex_float* z,
					lapack_int* ldz, lapack_int* isuppz,
					lapack_complex_float* work, lapack_int* lwork, float* rwork,
					lapack_int* lrwork, lapack_int* iwork, lapack_int* liwork,
					lapack_int *info );
void LAPACK_zheevr( char* jobz, char* range, char* uplo, lapack_int* n,
					lapack_complex_double* a, lapack_int* lda, double* vl,
					double* vu, lapack_int* il, lapack_int* iu, double* abstol,
					lapack_int* m, double* w, lapack_complex_double* z,
					lapack_int* ldz, lapack_int* isuppz,
					lapack_complex_double* work, lapack_int* lwork,
					double* rwork, lapack_int* lrwork, lapack_int* iwork,
					lapack_int* liwork, lapack_int *info );
void LAPACK_sspev( char* jobz, char* uplo, lapack_int* n, float* ap, float* w,
				float* z, lapack_int* ldz, float* work, lapack_int *info );
void LAPACK_dspev( char* jobz, char* uplo, lapack_int* n, double* ap, double* w,
				double* z, lapack_int* ldz, double* work, lapack_int *info );
void LAPACK_chpev( char* jobz, char* uplo, lapack_int* n,
				lapack_complex_float* ap, float* w, lapack_complex_float* z,
				lapack_int* ldz, lapack_complex_float* work, float* rwork,
				lapack_int *info );
void LAPACK_zhpev( char* jobz, char* uplo, lapack_int* n,
				lapack_complex_double* ap, double* w,
				lapack_complex_double* z, lapack_int* ldz,
				lapack_complex_double* work, double* rwork,
				lapack_int *info );
void LAPACK_sspevd( char* jobz, char* uplo, lapack_int* n, float* ap, float* w,
					float* z, lapack_int* ldz, float* work, lapack_int* lwork,
					lapack_int* iwork, lapack_int* liwork, lapack_int *info );
void LAPACK_dspevd( char* jobz, char* uplo, lapack_int* n, double* ap,
					double* w, double* z, lapack_int* ldz, double* work,
					lapack_int* lwork, lapack_int* iwork, lapack_int* liwork,
					lapack_int *info );
void LAPACK_chpevd( char* jobz, char* uplo, lapack_int* n,
					lapack_complex_float* ap, float* w, lapack_complex_float* z,
					lapack_int* ldz, lapack_complex_float* work,
					lapack_int* lwork, float* rwork, lapack_int* lrwork,
					lapack_int* iwork, lapack_int* liwork, lapack_int *info );
void LAPACK_zhpevd( char* jobz, char* uplo, lapack_int* n,
					lapack_complex_double* ap, double* w,
					lapack_complex_double* z, lapack_int* ldz,
					lapack_complex_double* work, lapack_int* lwork,
					double* rwork, lapack_int* lrwork, lapack_int* iwork,
					lapack_int* liwork, lapack_int *info );
void LAPACK_sspevx( char* jobz, char* range, char* uplo, lapack_int* n,
					float* ap, float* vl, float* vu, lapack_int* il,
					lapack_int* iu, float* abstol, lapack_int* m, float* w,
					float* z, lapack_int* ldz, float* work, lapack_int* iwork,
					lapack_int* ifail, lapack_int *info );
void LAPACK_dspevx( char* jobz, char* range, char* uplo, lapack_int* n,
					double* ap, double* vl, double* vu, lapack_int* il,
					lapack_int* iu, double* abstol, lapack_int* m, double* w,
					double* z, lapack_int* ldz, double* work, lapack_int* iwork,
					lapack_int* ifail, lapack_int *info );
void LAPACK_chpevx( char* jobz, char* range, char* uplo, lapack_int* n,
					lapack_complex_float* ap, float* vl, float* vu,
					lapack_int* il, lapack_int* iu, float* abstol,
					lapack_int* m, float* w, lapack_complex_float* z,
					lapack_int* ldz, lapack_complex_float* work, float* rwork,
					lapack_int* iwork, lapack_int* ifail, lapack_int *info );
void LAPACK_zhpevx( char* jobz, char* range, char* uplo, lapack_int* n,
					lapack_complex_double* ap, double* vl, double* vu,
					lapack_int* il, lapack_int* iu, double* abstol,
					lapack_int* m, double* w, lapack_complex_double* z,
					lapack_int* ldz, lapack_complex_double* work, double* rwork,
					lapack_int* iwork, lapack_int* ifail, lapack_int *info );
void LAPACK_ssbev( char* jobz, char* uplo, lapack_int* n, lapack_int* kd,
				float* ab, lapack_int* ldab, float* w, float* z,
				lapack_int* ldz, float* work, lapack_int *info );
void LAPACK_dsbev( char* jobz, char* uplo, lapack_int* n, lapack_int* kd,
				double* ab, lapack_int* ldab, double* w, double* z,
				lapack_int* ldz, double* work, lapack_int *info );
void LAPACK_chbev( char* jobz, char* uplo, lapack_int* n, lapack_int* kd,
				lapack_complex_float* ab, lapack_int* ldab, float* w,
				lapack_complex_float* z, lapack_int* ldz,
				lapack_complex_float* work, float* rwork, lapack_int *info );
void LAPACK_zhbev( char* jobz, char* uplo, lapack_int* n, lapack_int* kd,
				lapack_complex_double* ab, lapack_int* ldab, double* w,
				lapack_complex_double* z, lapack_int* ldz,
				lapack_complex_double* work, double* rwork,
				lapack_int *info );
void LAPACK_ssbevd( char* jobz, char* uplo, lapack_int* n, lapack_int* kd,
					float* ab, lapack_int* ldab, float* w, float* z,
					lapack_int* ldz, float* work, lapack_int* lwork,
					lapack_int* iwork, lapack_int* liwork, lapack_int *info );
void LAPACK_dsbevd( char* jobz, char* uplo, lapack_int* n, lapack_int* kd,
					double* ab, lapack_int* ldab, double* w, double* z,
					lapack_int* ldz, double* work, lapack_int* lwork,
					lapack_int* iwork, lapack_int* liwork, lapack_int *info );
void LAPACK_chbevd( char* jobz, char* uplo, lapack_int* n, lapack_int* kd,
					lapack_complex_float* ab, lapack_int* ldab, float* w,
					lapack_complex_float* z, lapack_int* ldz,
					lapack_complex_float* work, lapack_int* lwork, float* rwork,
					lapack_int* lrwork, lapack_int* iwork, lapack_int* liwork,
					lapack_int *info );
void LAPACK_zhbevd( char* jobz, char* uplo, lapack_int* n, lapack_int* kd,
					lapack_complex_double* ab, lapack_int* ldab, double* w,
					lapack_complex_double* z, lapack_int* ldz,
					lapack_complex_double* work, lapack_int* lwork,
					double* rwork, lapack_int* lrwork, lapack_int* iwork,
					lapack_int* liwork, lapack_int *info );
void LAPACK_ssbevx( char* jobz, char* range, char* uplo, lapack_int* n,
					lapack_int* kd, float* ab, lapack_int* ldab, float* q,
					lapack_int* ldq, float* vl, float* vu, lapack_int* il,
					lapack_int* iu, float* abstol, lapack_int* m, float* w,
					float* z, lapack_int* ldz, float* work, lapack_int* iwork,
					lapack_int* ifail, lapack_int *info );
void LAPACK_dsbevx( char* jobz, char* range, char* uplo, lapack_int* n,
					lapack_int* kd, double* ab, lapack_int* ldab, double* q,
					lapack_int* ldq, double* vl, double* vu, lapack_int* il,
					lapack_int* iu, double* abstol, lapack_int* m, double* w,
					double* z, lapack_int* ldz, double* work, lapack_int* iwork,
					lapack_int* ifail, lapack_int *info );
void LAPACK_chbevx( char* jobz, char* range, char* uplo, lapack_int* n,
					lapack_int* kd, lapack_complex_float* ab, lapack_int* ldab,
					lapack_complex_float* q, lapack_int* ldq, float* vl,
					float* vu, lapack_int* il, lapack_int* iu, float* abstol,
					lapack_int* m, float* w, lapack_complex_float* z,
					lapack_int* ldz, lapack_complex_float* work, float* rwork,
					lapack_int* iwork, lapack_int* ifail, lapack_int *info );
void LAPACK_zhbevx( char* jobz, char* range, char* uplo, lapack_int* n,
					lapack_int* kd, lapack_complex_double* ab, lapack_int* ldab,
					lapack_complex_double* q, lapack_int* ldq, double* vl,
					double* vu, lapack_int* il, lapack_int* iu, double* abstol,
					lapack_int* m, double* w, lapack_complex_double* z,
					lapack_int* ldz, lapack_complex_double* work, double* rwork,
					lapack_int* iwork, lapack_int* ifail, lapack_int *info );
void LAPACK_sstev( char* jobz, lapack_int* n, float* d, float* e, float* z,
				lapack_int* ldz, float* work, lapack_int *info );
void LAPACK_dstev( char* jobz, lapack_int* n, double* d, double* e, double* z,
				lapack_int* ldz, double* work, lapack_int *info );
void LAPACK_sstevd( char* jobz, lapack_int* n, float* d, float* e, float* z,
					lapack_int* ldz, float* work, lapack_int* lwork,
					lapack_int* iwork, lapack_int* liwork, lapack_int *info );
void LAPACK_dstevd( char* jobz, lapack_int* n, double* d, double* e, double* z,
					lapack_int* ldz, double* work, lapack_int* lwork,
					lapack_int* iwork, lapack_int* liwork, lapack_int *info );
void LAPACK_sstevx( char* jobz, char* range, lapack_int* n, float* d, float* e,
					float* vl, float* vu, lapack_int* il, lapack_int* iu,
					float* abstol, lapack_int* m, float* w, float* z,
					lapack_int* ldz, float* work, lapack_int* iwork,
					lapack_int* ifail, lapack_int *info );
void LAPACK_dstevx( char* jobz, char* range, lapack_int* n, double* d,
					double* e, double* vl, double* vu, lapack_int* il,
					lapack_int* iu, double* abstol, lapack_int* m, double* w,
					double* z, lapack_int* ldz, double* work, lapack_int* iwork,
					lapack_int* ifail, lapack_int *info );
void LAPACK_sstevr( char* jobz, char* range, lapack_int* n, float* d, float* e,
					float* vl, float* vu, lapack_int* il, lapack_int* iu,
					float* abstol, lapack_int* m, float* w, float* z,
					lapack_int* ldz, lapack_int* isuppz, float* work,
					lapack_int* lwork, lapack_int* iwork, lapack_int* liwork,
					lapack_int *info );
void LAPACK_dstevr( char* jobz, char* range, lapack_int* n, double* d,
					double* e, double* vl, double* vu, lapack_int* il,
					lapack_int* iu, double* abstol, lapack_int* m, double* w,
					double* z, lapack_int* ldz, lapack_int* isuppz,
					double* work, lapack_int* lwork, lapack_int* iwork,
					lapack_int* liwork, lapack_int *info );
void LAPACK_sgees( char* jobvs, char* sort, LAPACK_S_SELECT2 select,
				lapack_int* n, float* a, lapack_int* lda, lapack_int* sdim,
				float* wr, float* wi, float* vs, lapack_int* ldvs,
				float* work, lapack_int* lwork, lapack_logical* bwork,
				lapack_int *info );
void LAPACK_dgees( char* jobvs, char* sort, LAPACK_D_SELECT2 select,
				lapack_int* n, double* a, lapack_int* lda, lapack_int* sdim,
				double* wr, double* wi, double* vs, lapack_int* ldvs,
				double* work, lapack_int* lwork, lapack_logical* bwork,
				lapack_int *info );
void LAPACK_cgees( char* jobvs, char* sort, LAPACK_C_SELECT1 select,
				lapack_int* n, lapack_complex_float* a, lapack_int* lda,
				lapack_int* sdim, lapack_complex_float* w,
				lapack_complex_float* vs, lapack_int* ldvs,
				lapack_complex_float* work, lapack_int* lwork, float* rwork,
				lapack_logical* bwork, lapack_int *info );
void LAPACK_zgees( char* jobvs, char* sort, LAPACK_Z_SELECT1 select,
				lapack_int* n, lapack_complex_double* a, lapack_int* lda,
				lapack_int* sdim, lapack_complex_double* w,
				lapack_complex_double* vs, lapack_int* ldvs,
				lapack_complex_double* work, lapack_int* lwork,
				double* rwork, lapack_logical* bwork, lapack_int *info );
void LAPACK_sgeesx( char* jobvs, char* sort, LAPACK_S_SELECT2 select,
					char* sense, lapack_int* n, float* a, lapack_int* lda,
					lapack_int* sdim, float* wr, float* wi, float* vs,
					lapack_int* ldvs, float* rconde, float* rcondv, float* work,
					lapack_int* lwork, lapack_int* iwork, lapack_int* liwork,
					lapack_logical* bwork, lapack_int *info );
void LAPACK_dgeesx( char* jobvs, char* sort, LAPACK_D_SELECT2 select,
					char* sense, lapack_int* n, double* a, lapack_int* lda,
					lapack_int* sdim, double* wr, double* wi, double* vs,
					lapack_int* ldvs, double* rconde, double* rcondv,
					double* work, lapack_int* lwork, lapack_int* iwork,
					lapack_int* liwork, lapack_logical* bwork,
					lapack_int *info );
void LAPACK_cgeesx( char* jobvs, char* sort, LAPACK_C_SELECT1 select,
					char* sense, lapack_int* n, lapack_complex_float* a,
					lapack_int* lda, lapack_int* sdim, lapack_complex_float* w,
					lapack_complex_float* vs, lapack_int* ldvs, float* rconde,
					float* rcondv, lapack_complex_float* work,
					lapack_int* lwork, float* rwork, lapack_logical* bwork,
					lapack_int *info );
void LAPACK_zgeesx( char* jobvs, char* sort, LAPACK_Z_SELECT1 select,
					char* sense, lapack_int* n, lapack_complex_double* a,
					lapack_int* lda, lapack_int* sdim, lapack_complex_double* w,
					lapack_complex_double* vs, lapack_int* ldvs, double* rconde,
					double* rcondv, lapack_complex_double* work,
					lapack_int* lwork, double* rwork, lapack_logical* bwork,
					lapack_int *info );
void LAPACK_sgeev( char* jobvl, char* jobvr, lapack_int* n, float* a,
				lapack_int* lda, float* wr, float* wi, float* vl,
				lapack_int* ldvl, float* vr, lapack_int* ldvr, float* work,
				lapack_int* lwork, lapack_int *info );
void LAPACK_dgeev( char* jobvl, char* jobvr, lapack_int* n, double* a,
				lapack_int* lda, double* wr, double* wi, double* vl,
				lapack_int* ldvl, double* vr, lapack_int* ldvr, double* work,
				lapack_int* lwork, lapack_int *info );
void LAPACK_cgeev( char* jobvl, char* jobvr, lapack_int* n,
				lapack_complex_float* a, lapack_int* lda,
				lapack_complex_float* w, lapack_complex_float* vl,
				lapack_int* ldvl, lapack_complex_float* vr, lapack_int* ldvr,
				lapack_complex_float* work, lapack_int* lwork, float* rwork,
				lapack_int *info );
void LAPACK_zgeev( char* jobvl, char* jobvr, lapack_int* n,
				lapack_complex_double* a, lapack_int* lda,
				lapack_complex_double* w, lapack_complex_double* vl,
				lapack_int* ldvl, lapack_complex_double* vr,
				lapack_int* ldvr, lapack_complex_double* work,
				lapack_int* lwork, double* rwork, lapack_int *info );
void LAPACK_sgeevx( char* balanc, char* jobvl, char* jobvr, char* sense,
					lapack_int* n, float* a, lapack_int* lda, float* wr,
					float* wi, float* vl, lapack_int* ldvl, float* vr,
					lapack_int* ldvr, lapack_int* ilo, lapack_int* ihi,
					float* scale, float* abnrm, float* rconde, float* rcondv,
					float* work, lapack_int* lwork, lapack_int* iwork,
					lapack_int *info );
void LAPACK_dgeevx( char* balanc, char* jobvl, char* jobvr, char* sense,
					lapack_int* n, double* a, lapack_int* lda, double* wr,
					double* wi, double* vl, lapack_int* ldvl, double* vr,
					lapack_int* ldvr, lapack_int* ilo, lapack_int* ihi,
					double* scale, double* abnrm, double* rconde,
					double* rcondv, double* work, lapack_int* lwork,
					lapack_int* iwork, lapack_int *info );
void LAPACK_cgeevx( char* balanc, char* jobvl, char* jobvr, char* sense,
					lapack_int* n, lapack_complex_float* a, lapack_int* lda,
					lapack_complex_float* w, lapack_complex_float* vl,
					lapack_int* ldvl, lapack_complex_float* vr,
					lapack_int* ldvr, lapack_int* ilo, lapack_int* ihi,
					float* scale, float* abnrm, float* rconde, float* rcondv,
					lapack_complex_float* work, lapack_int* lwork, float* rwork,
					lapack_int *info );
void LAPACK_zgeevx( char* balanc, char* jobvl, char* jobvr, char* sense,
					lapack_int* n, lapack_complex_double* a, lapack_int* lda,
					lapack_complex_double* w, lapack_complex_double* vl,
					lapack_int* ldvl, lapack_complex_double* vr,
					lapack_int* ldvr, lapack_int* ilo, lapack_int* ihi,
					double* scale, double* abnrm, double* rconde,
					double* rcondv, lapack_complex_double* work,
					lapack_int* lwork, double* rwork, lapack_int *info );
void LAPACK_sgesvd( char* jobu, char* jobvt, lapack_int* m, lapack_int* n,
					float* a, lapack_int* lda, float* s, float* u,
					lapack_int* ldu, float* vt, lapack_int* ldvt, float* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_dgesvd( char* jobu, char* jobvt, lapack_int* m, lapack_int* n,
					double* a, lapack_int* lda, double* s, double* u,
					lapack_int* ldu, double* vt, lapack_int* ldvt, double* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_cgesvd( char* jobu, char* jobvt, lapack_int* m, lapack_int* n,
					lapack_complex_float* a, lapack_int* lda, float* s,
					lapack_complex_float* u, lapack_int* ldu,
					lapack_complex_float* vt, lapack_int* ldvt,
					lapack_complex_float* work, lapack_int* lwork, float* rwork,
					lapack_int *info );
void LAPACK_zgesvd( char* jobu, char* jobvt, lapack_int* m, lapack_int* n,
					lapack_complex_double* a, lapack_int* lda, double* s,
					lapack_complex_double* u, lapack_int* ldu,
					lapack_complex_double* vt, lapack_int* ldvt,
					lapack_complex_double* work, lapack_int* lwork,
					double* rwork, lapack_int *info );
void LAPACK_sgesdd( char* jobz, lapack_int* m, lapack_int* n, float* a,
					lapack_int* lda, float* s, float* u, lapack_int* ldu,
					float* vt, lapack_int* ldvt, float* work, lapack_int* lwork,
					lapack_int* iwork, lapack_int *info );
void LAPACK_dgesdd( char* jobz, lapack_int* m, lapack_int* n, double* a,
					lapack_int* lda, double* s, double* u, lapack_int* ldu,
					double* vt, lapack_int* ldvt, double* work,
					lapack_int* lwork, lapack_int* iwork, lapack_int *info );
void LAPACK_cgesdd( char* jobz, lapack_int* m, lapack_int* n,
					lapack_complex_float* a, lapack_int* lda, float* s,
					lapack_complex_float* u, lapack_int* ldu,
					lapack_complex_float* vt, lapack_int* ldvt,
					lapack_complex_float* work, lapack_int* lwork, float* rwork,
					lapack_int* iwork, lapack_int *info );
void LAPACK_zgesdd( char* jobz, lapack_int* m, lapack_int* n,
					lapack_complex_double* a, lapack_int* lda, double* s,
					lapack_complex_double* u, lapack_int* ldu,
					lapack_complex_double* vt, lapack_int* ldvt,
					lapack_complex_double* work, lapack_int* lwork,
					double* rwork, lapack_int* iwork, lapack_int *info );
void LAPACK_dgejsv( char* joba, char* jobu, char* jobv, char* jobr, char* jobt,
					char* jobp, lapack_int* m, lapack_int* n, double* a,
					lapack_int* lda, double* sva, double* u, lapack_int* ldu,
					double* v, lapack_int* ldv, double* work, lapack_int* lwork,
					lapack_int* iwork, lapack_int *info );
void LAPACK_sgejsv( char* joba, char* jobu, char* jobv, char* jobr, char* jobt,
					char* jobp, lapack_int* m, lapack_int* n, float* a,
					lapack_int* lda, float* sva, float* u, lapack_int* ldu,
					float* v, lapack_int* ldv, float* work, lapack_int* lwork,
					lapack_int* iwork, lapack_int *info );
void LAPACK_dgesvj( char* joba, char* jobu, char* jobv, lapack_int* m,
					lapack_int* n, double* a, lapack_int* lda, double* sva,
					lapack_int* mv, double* v, lapack_int* ldv, double* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_sgesvj( char* joba, char* jobu, char* jobv, lapack_int* m,
					lapack_int* n, float* a, lapack_int* lda, float* sva,
					lapack_int* mv, float* v, lapack_int* ldv, float* work,
					lapack_int* lwork, lapack_int *info );
void LAPACK_sggsvd( char* jobu, char* jobv, char* jobq, lapack_int* m,
					lapack_int* n, lapack_int* p, lapack_int* k, lapack_int* l,
					float* a, lapack_int* lda, float* b, lapack_int* ldb,
					float* alpha, float* beta, float* u, lapack_int* ldu,
					float* v, lapack_int* ldv, float* q, lapack_int* ldq,
					float* work, lapack_int* iwork, lapack_int *info );
void LAPACK_dggsvd( char* jobu, char* jobv, char* jobq, lapack_int* m,
					lapack_int* n, lapack_int* p, lapack_int* k, lapack_int* l,
					double* a, lapack_int* lda, double* b, lapack_int* ldb,
					double* alpha, double* beta, double* u, lapack_int* ldu,
					double* v, lapack_int* ldv, double* q, lapack_int* ldq,
					double* work, lapack_int* iwork, lapack_int *info );
void LAPACK_cggsvd( char* jobu, char* jobv, char* jobq, lapack_int* m,
					lapack_int* n, lapack_int* p, lapack_int* k, lapack_int* l,
					lapack_complex_float* a, lapack_int* lda,
					lapack_complex_float* b, lapack_int* ldb, float* alpha,
					float* beta, lapack_complex_float* u, lapack_int* ldu,
					lapack_complex_float* v, lapack_int* ldv,
					lapack_complex_float* q, lapack_int* ldq,
					lapack_complex_float* work, float* rwork, lapack_int* iwork,
					lapack_int *info );
void LAPACK_zggsvd( char* jobu, char* jobv, char* jobq, lapack_int* m,
					lapack_int* n, lapack_int* p, lapack_int* k, lapack_int* l,
					lapack_complex_double* a, lapack_int* lda,
					lapack_complex_double* b, lapack_int* ldb, double* alpha,
					double* beta, lapack_complex_double* u, lapack_int* ldu,
					lapack_complex_double* v, lapack_int* ldv,
					lapack_complex_double* q, lapack_int* ldq,
					lapack_complex_double* work, double* rwork,
					lapack_int* iwork, lapack_int *info );
void LAPACK_ssygv( lapack_int* itype, char* jobz, char* uplo, lapack_int* n,
				float* a, lapack_int* lda, float* b, lapack_int* ldb,
				float* w, float* work, lapack_int* lwork, lapack_int *info );
void LAPACK_dsygv( lapack_int* itype, char* jobz, char* uplo, lapack_int* n,
				double* a, lapack_int* lda, double* b, lapack_int* ldb,
				double* w, double* work, lapack_int* lwork,
				lapack_int *info );
void LAPACK_chegv( lapack_int* itype, char* jobz, char* uplo, lapack_int* n,
				lapack_complex_float* a, lapack_int* lda,
				lapack_complex_float* b, lapack_int* ldb, float* w,
				lapack_complex_float* work, lapack_int* lwork, float* rwork,
				lapack_int *info );
void LAPACK_zhegv( lapack_int* itype, char* jobz, char* uplo, lapack_int* n,
				lapack_complex_double* a, lapack_int* lda,
				lapack_complex_double* b, lapack_int* ldb, double* w,
				lapack_complex_double* work, lapack_int* lwork,
				double* rwork, lapack_int *info );
void LAPACK_ssygvd( lapack_int* itype, char* jobz, char* uplo, lapack_int* n,
					float* a, lapack_int* lda, float* b, lapack_int* ldb,
					float* w, float* work, lapack_int* lwork, lapack_int* iwork,
					lapack_int* liwork, lapack_int *info );
void LAPACK_dsygvd( lapack_int* itype, char* jobz, char* uplo, lapack_int* n,
					double* a, lapack_int* lda, double* b, lapack_int* ldb,
					double* w, double* work, lapack_int* lwork,
					lapack_int* iwork, lapack_int* liwork, lapack_int *info );
void LAPACK_chegvd( lapack_int* itype, char* jobz, char* uplo, lapack_int* n,
					lapack_complex_float* a, lapack_int* lda,
					lapack_complex_float* b, lapack_int* ldb, float* w,
					lapack_complex_float* work, lapack_int* lwork, float* rwork,
					lapack_int* lrwork, lapack_int* iwork, lapack_int* liwork,
					lapack_int *info );
void LAPACK_zhegvd( lapack_int* itype, char* jobz, char* uplo, lapack_int* n,
					lapack_complex_double* a, lapack_int* lda,
					lapack_complex_double* b, lapack_int* ldb, double* w,
					lapack_complex_double* work, lapack_int* lwork,
					double* rwork, lapack_int* lrwork, lapack_int* iwork,
					lapack_int* liwork, lapack_int *info );
void LAPACK_ssygvx( lapack_int* itype, char* jobz, char* range, char* uplo,
					lapack_int* n, float* a, lapack_int* lda, float* b,
					lapack_int* ldb, float* vl, float* vu, lapack_int* il,
					lapack_int* iu, float* abstol, lapack_int* m, float* w,
					float* z, lapack_int* ldz, float* work, lapack_int* lwork,
					lapack_int* iwork, lapack_int* ifail, lapack_int *info );
void LAPACK_dsygvx( lapack_int* itype, char* jobz, char* range, char* uplo,
					lapack_int* n, double* a, lapack_int* lda, double* b,
					lapack_int* ldb, double* vl, double* vu, lapack_int* il,
					lapack_int* iu, double* abstol, lapack_int* m, double* w,
					double* z, lapack_int* ldz, double* work, lapack_int* lwork,
					lapack_int* iwork, lapack_int* ifail, lapack_int *info );
void LAPACK_chegvx( lapack_int* itype, char* jobz, char* range, char* uplo,
					lapack_int* n, lapack_complex_float* a, lapack_int* lda,
					lapack_complex_float* b, lapack_int* ldb, float* vl,
					float* vu, lapack_int* il, lapack_int* iu, float* abstol,
					lapack_int* m, float* w, lapack_complex_float* z,
					lapack_int* ldz, lapack_complex_float* work,
					lapack_int* lwork, float* rwork, lapack_int* iwork,
					lapack_int* ifail, lapack_int *info );
void LAPACK_zhegvx( lapack_int* itype, char* jobz, char* range, char* uplo,
					lapack_int* n, lapack_complex_double* a, lapack_int* lda,
					lapack_complex_double* b, lapack_int* ldb, double* vl,
					double* vu, lapack_int* il, lapack_int* iu, double* abstol,
					lapack_int* m, double* w, lapack_complex_double* z,
					lapack_int* ldz, lapack_complex_double* work,
					lapack_int* lwork, double* rwork, lapack_int* iwork,
					lapack_int* ifail, lapack_int *info );
void LAPACK_sspgv( lapack_int* itype, char* jobz, char* uplo, lapack_int* n,
				float* ap, float* bp, float* w, float* z, lapack_int* ldz,
				float* work, lapack_int *info );
void LAPACK_dspgv( lapack_int* itype, char* jobz, char* uplo, lapack_int* n,
				double* ap, double* bp, double* w, double* z,
				lapack_int* ldz, double* work, lapack_int *info );
void LAPACK_chpgv( lapack_int* itype, char* jobz, char* uplo, lapack_int* n,
				lapack_complex_float* ap, lapack_complex_float* bp, float* w,
				lapack_complex_float* z, lapack_int* ldz,
				lapack_complex_float* work, float* rwork, lapack_int *info );
void LAPACK_zhpgv( lapack_int* itype, char* jobz, char* uplo, lapack_int* n,
				lapack_complex_double* ap, lapack_complex_double* bp,
				double* w, lapack_complex_double* z, lapack_int* ldz,
				lapack_complex_double* work, double* rwork,
				lapack_int *info );
void LAPACK_sspgvd( lapack_int* itype, char* jobz, char* uplo, lapack_int* n,
					float* ap, float* bp, float* w, float* z, lapack_int* ldz,
					float* work, lapack_int* lwork, lapack_int* iwork,
					lapack_int* liwork, lapack_int *info );
void LAPACK_dspgvd( lapack_int* itype, char* jobz, char* uplo, lapack_int* n,
					double* ap, double* bp, double* w, double* z,
					lapack_int* ldz, double* work, lapack_int* lwork,
					lapack_int* iwork, lapack_int* liwork, lapack_int *info );
void LAPACK_chpgvd( lapack_int* itype, char* jobz, char* uplo, lapack_int* n,
					lapack_complex_float* ap, lapack_complex_float* bp,
					float* w, lapack_complex_float* z, lapack_int* ldz,
					lapack_complex_float* work, lapack_int* lwork, float* rwork,
					lapack_int* lrwork, lapack_int* iwork, lapack_int* liwork,
					lapack_int *info );
void LAPACK_zhpgvd( lapack_int* itype, char* jobz, char* uplo, lapack_int* n,
					lapack_complex_double* ap, lapack_complex_double* bp,
					double* w, lapack_complex_double* z, lapack_int* ldz,
					lapack_complex_double* work, lapack_int* lwork,
					double* rwork, lapack_int* lrwork, lapack_int* iwork,
					lapack_int* liwork, lapack_int *info );
void LAPACK_sspgvx( lapack_int* itype, char* jobz, char* range, char* uplo,
					lapack_int* n, float* ap, float* bp, float* vl, float* vu,
					lapack_int* il, lapack_int* iu, float* abstol,
					lapack_int* m, float* w, float* z, lapack_int* ldz,
					float* work, lapack_int* iwork, lapack_int* ifail,
					lapack_int *info );
void LAPACK_dspgvx( lapack_int* itype, char* jobz, char* range, char* uplo,
					lapack_int* n, double* ap, double* bp, double* vl,
					double* vu, lapack_int* il, lapack_int* iu, double* abstol,
					lapack_int* m, double* w, double* z, lapack_int* ldz,
					double* work, lapack_int* iwork, lapack_int* ifail,
					lapack_int *info );
void LAPACK_chpgvx( lapack_int* itype, char* jobz, char* range, char* uplo,
					lapack_int* n, lapack_complex_float* ap,
					lapack_complex_float* bp, float* vl, float* vu,
					lapack_int* il, lapack_int* iu, float* abstol,
					lapack_int* m, float* w, lapack_complex_float* z,
					lapack_int* ldz, lapack_complex_float* work, float* rwork,
					lapack_int* iwork, lapack_int* ifail, lapack_int *info );
void LAPACK_zhpgvx( lapack_int* itype, char* jobz, char* range, char* uplo,
					lapack_int* n, lapack_complex_double* ap,
					lapack_complex_double* bp, double* vl, double* vu,
					lapack_int* il, lapack_int* iu, double* abstol,
					lapack_int* m, double* w, lapack_complex_double* z,
					lapack_int* ldz, lapack_complex_double* work, double* rwork,
					lapack_int* iwork, lapack_int* ifail, lapack_int *info );
void LAPACK_ssbgv( char* jobz, char* uplo, lapack_int* n, lapack_int* ka,
				lapack_int* kb, float* ab, lapack_int* ldab, float* bb,
				lapack_int* ldbb, float* w, float* z, lapack_int* ldz,
				float* work, lapack_int *info );
void LAPACK_dsbgv( char* jobz, char* uplo, lapack_int* n, lapack_int* ka,
				lapack_int* kb, double* ab, lapack_int* ldab, double* bb,
				lapack_int* ldbb, double* w, double* z, lapack_int* ldz,
				double* work, lapack_int *info );
void LAPACK_chbgv( char* jobz, char* uplo, lapack_int* n, lapack_int* ka,
				lapack_int* kb, lapack_complex_float* ab, lapack_int* ldab,
				lapack_complex_float* bb, lapack_int* ldbb, float* w,
				lapack_complex_float* z, lapack_int* ldz,
				lapack_complex_float* work, float* rwork, lapack_int *info );
void LAPACK_zhbgv( char* jobz, char* uplo, lapack_int* n, lapack_int* ka,
				lapack_int* kb, lapack_complex_double* ab, lapack_int* ldab,
				lapack_complex_double* bb, lapack_int* ldbb, double* w,
				lapack_complex_double* z, lapack_int* ldz,
				lapack_complex_double* work, double* rwork,
				lapack_int *info );
void LAPACK_ssbgvd( char* jobz, char* uplo, lapack_int* n, lapack_int* ka,
					lapack_int* kb, float* ab, lapack_int* ldab, float* bb,
					lapack_int* ldbb, float* w, float* z, lapack_int* ldz,
					float* work, lapack_int* lwork, lapack_int* iwork,
					lapack_int* liwork, lapack_int *info );
void LAPACK_dsbgvd( char* jobz, char* uplo, lapack_int* n, lapack_int* ka,
					lapack_int* kb, double* ab, lapack_int* ldab, double* bb,
					lapack_int* ldbb, double* w, double* z, lapack_int* ldz,
					double* work, lapack_int* lwork, lapack_int* iwork,
					lapack_int* liwork, lapack_int *info );
void LAPACK_chbgvd( char* jobz, char* uplo, lapack_int* n, lapack_int* ka,
					lapack_int* kb, lapack_complex_float* ab, lapack_int* ldab,
					lapack_complex_float* bb, lapack_int* ldbb, float* w,
					lapack_complex_float* z, lapack_int* ldz,
					lapack_complex_float* work, lapack_int* lwork, float* rwork,
					lapack_int* lrwork, lapack_int* iwork, lapack_int* liwork,
					lapack_int *info );
void LAPACK_zhbgvd( char* jobz, char* uplo, lapack_int* n, lapack_int* ka,
					lapack_int* kb, lapack_complex_double* ab, lapack_int* ldab,
					lapack_complex_double* bb, lapack_int* ldbb, double* w,
					lapack_complex_double* z, lapack_int* ldz,
					lapack_complex_double* work, lapack_int* lwork,
					double* rwork, lapack_int* lrwork, lapack_int* iwork,
					lapack_int* liwork, lapack_int *info );
void LAPACK_ssbgvx( char* jobz, char* range, char* uplo, lapack_int* n,
					lapack_int* ka, lapack_int* kb, float* ab, lapack_int* ldab,
					float* bb, lapack_int* ldbb, float* q, lapack_int* ldq,
					float* vl, float* vu, lapack_int* il, lapack_int* iu,
					float* abstol, lapack_int* m, float* w, float* z,
					lapack_int* ldz, float* work, lapack_int* iwork,
					lapack_int* ifail, lapack_int *info );
void LAPACK_dsbgvx( char* jobz, char* range, char* uplo, lapack_int* n,
					lapack_int* ka, lapack_int* kb, double* ab,
					lapack_int* ldab, double* bb, lapack_int* ldbb, double* q,
					lapack_int* ldq, double* vl, double* vu, lapack_int* il,
					lapack_int* iu, double* abstol, lapack_int* m, double* w,
					double* z, lapack_int* ldz, double* work, lapack_int* iwork,
					lapack_int* ifail, lapack_int *info );
void LAPACK_chbgvx( char* jobz, char* range, char* uplo, lapack_int* n,
					lapack_int* ka, lapack_int* kb, lapack_complex_float* ab,
					lapack_int* ldab, lapack_complex_float* bb,
					lapack_int* ldbb, lapack_complex_float* q, lapack_int* ldq,
					float* vl, float* vu, lapack_int* il, lapack_int* iu,
					float* abstol, lapack_int* m, float* w,
					lapack_complex_float* z, lapack_int* ldz,
					lapack_complex_float* work, float* rwork, lapack_int* iwork,
					lapack_int* ifail, lapack_int *info );
void LAPACK_zhbgvx( char* jobz, char* range, char* uplo, lapack_int* n,
					lapack_int* ka, lapack_int* kb, lapack_complex_double* ab,
					lapack_int* ldab, lapack_complex_double* bb,
					lapack_int* ldbb, lapack_complex_double* q, lapack_int* ldq,
					double* vl, double* vu, lapack_int* il, lapack_int* iu,
					double* abstol, lapack_int* m, double* w,
					lapack_complex_double* z, lapack_int* ldz,
					lapack_complex_double* work, double* rwork,
					lapack_int* iwork, lapack_int* ifail, lapack_int *info );
void LAPACK_sgges( char* jobvsl, char* jobvsr, char* sort,
				LAPACK_S_SELECT3 selctg, lapack_int* n, float* a,
				lapack_int* lda, float* b, lapack_int* ldb, lapack_int* sdim,
				float* alphar, float* alphai, float* beta, float* vsl,
				lapack_int* ldvsl, float* vsr, lapack_int* ldvsr,
				float* work, lapack_int* lwork, lapack_logical* bwork,
				lapack_int *info );
void LAPACK_dgges( char* jobvsl, char* jobvsr, char* sort,
				LAPACK_D_SELECT3 selctg, lapack_int* n, double* a,
				lapack_int* lda, double* b, lapack_int* ldb,
				lapack_int* sdim, double* alphar, double* alphai,
				double* beta, double* vsl, lapack_int* ldvsl, double* vsr,
				lapack_int* ldvsr, double* work, lapack_int* lwork,
				lapack_logical* bwork, lapack_int *info );
void LAPACK_cgges( char* jobvsl, char* jobvsr, char* sort,
				LAPACK_C_SELECT2 selctg, lapack_int* n,
				lapack_complex_float* a, lapack_int* lda,
				lapack_complex_float* b, lapack_int* ldb, lapack_int* sdim,
				lapack_complex_float* alpha, lapack_complex_float* beta,
				lapack_complex_float* vsl, lapack_int* ldvsl,
				lapack_complex_float* vsr, lapack_int* ldvsr,
				lapack_complex_float* work, lapack_int* lwork, float* rwork,
				lapack_logical* bwork, lapack_int *info );
void LAPACK_zgges( char* jobvsl, char* jobvsr, char* sort,
				LAPACK_Z_SELECT2 selctg, lapack_int* n,
				lapack_complex_double* a, lapack_int* lda,
				lapack_complex_double* b, lapack_int* ldb, lapack_int* sdim,
				lapack_complex_double* alpha, lapack_complex_double* beta,
				lapack_complex_double* vsl, lapack_int* ldvsl,
				lapack_complex_double* vsr, lapack_int* ldvsr,
				lapack_complex_double* work, lapack_int* lwork,
				double* rwork, lapack_logical* bwork, lapack_int *info );
void LAPACK_sggesx( char* jobvsl, char* jobvsr, char* sort,
					LAPACK_S_SELECT3 selctg, char* sense, lapack_int* n,
					float* a, lapack_int* lda, float* b, lapack_int* ldb,
					lapack_int* sdim, float* alphar, float* alphai, float* beta,
					float* vsl, lapack_int* ldvsl, float* vsr,
					lapack_int* ldvsr, float* rconde, float* rcondv,
					float* work, lapack_int* lwork, lapack_int* iwork,
					lapack_int* liwork, lapack_logical* bwork,
					lapack_int *info );
void LAPACK_dggesx( char* jobvsl, char* jobvsr, char* sort,
					LAPACK_D_SELECT3 selctg, char* sense, lapack_int* n,
					double* a, lapack_int* lda, double* b, lapack_int* ldb,
					lapack_int* sdim, double* alphar, double* alphai,
					double* beta, double* vsl, lapack_int* ldvsl, double* vsr,
					lapack_int* ldvsr, double* rconde, double* rcondv,
					double* work, lapack_int* lwork, lapack_int* iwork,
					lapack_int* liwork, lapack_logical* bwork,
					lapack_int *info );
void LAPACK_cggesx( char* jobvsl, char* jobvsr, char* sort,
					LAPACK_C_SELECT2 selctg, char* sense, lapack_int* n,
					lapack_complex_float* a, lapack_int* lda,
					lapack_complex_float* b, lapack_int* ldb, lapack_int* sdim,
					lapack_complex_float* alpha, lapack_complex_float* beta,
					lapack_complex_float* vsl, lapack_int* ldvsl,
					lapack_complex_float* vsr, lapack_int* ldvsr, float* rconde,
					float* rcondv, lapack_complex_float* work,
					lapack_int* lwork, float* rwork, lapack_int* iwork,
					lapack_int* liwork, lapack_logical* bwork,
					lapack_int *info );
void LAPACK_zggesx( char* jobvsl, char* jobvsr, char* sort,
					LAPACK_Z_SELECT2 selctg, char* sense, lapack_int* n,
					lapack_complex_double* a, lapack_int* lda,
					lapack_complex_double* b, lapack_int* ldb, lapack_int* sdim,
					lapack_complex_double* alpha, lapack_complex_double* beta,
					lapack_complex_double* vsl, lapack_int* ldvsl,
					lapack_complex_double* vsr, lapack_int* ldvsr,
					double* rconde, double* rcondv, lapack_complex_double* work,
					lapack_int* lwork, double* rwork, lapack_int* iwork,
					lapack_int* liwork, lapack_logical* bwork,
					lapack_int *info );
void LAPACK_sggev( char* jobvl, char* jobvr, lapack_int* n, float* a,
				lapack_int* lda, float* b, lapack_int* ldb, float* alphar,
				float* alphai, float* beta, float* vl, lapack_int* ldvl,
				float* vr, lapack_int* ldvr, float* work, lapack_int* lwork,
				lapack_int *info );
void LAPACK_dggev( char* jobvl, char* jobvr, lapack_int* n, double* a,
				lapack_int* lda, double* b, lapack_int* ldb, double* alphar,
				double* alphai, double* beta, double* vl, lapack_int* ldvl,
				double* vr, lapack_int* ldvr, double* work,
				lapack_int* lwork, lapack_int *info );
void LAPACK_cggev( char* jobvl, char* jobvr, lapack_int* n,
				lapack_complex_float* a, lapack_int* lda,
				lapack_complex_float* b, lapack_int* ldb,
				lapack_complex_float* alpha, lapack_complex_float* beta,
				lapack_complex_float* vl, lapack_int* ldvl,
				lapack_complex_float* vr, lapack_int* ldvr,
				lapack_complex_float* work, lapack_int* lwork, float* rwork,
				lapack_int *info );
void LAPACK_zggev( char* jobvl, char* jobvr, lapack_int* n,
				lapack_complex_double* a, lapack_int* lda,
				lapack_complex_double* b, lapack_int* ldb,
				lapack_complex_double* alpha, lapack_complex_double* beta,
				lapack_complex_double* vl, lapack_int* ldvl,
				lapack_complex_double* vr, lapack_int* ldvr,
				lapack_complex_double* work, lapack_int* lwork,
				double* rwork, lapack_int *info );
void LAPACK_sggevx( char* balanc, char* jobvl, char* jobvr, char* sense,
					lapack_int* n, float* a, lapack_int* lda, float* b,
					lapack_int* ldb, float* alphar, float* alphai, float* beta,
					float* vl, lapack_int* ldvl, float* vr, lapack_int* ldvr,
					lapack_int* ilo, lapack_int* ihi, float* lscale,
					float* rscale, float* abnrm, float* bbnrm, float* rconde,
					float* rcondv, float* work, lapack_int* lwork,
					lapack_int* iwork, lapack_logical* bwork,
					lapack_int *info );
void LAPACK_dggevx( char* balanc, char* jobvl, char* jobvr, char* sense,
					lapack_int* n, double* a, lapack_int* lda, double* b,
					lapack_int* ldb, double* alphar, double* alphai,
					double* beta, double* vl, lapack_int* ldvl, double* vr,
					lapack_int* ldvr, lapack_int* ilo, lapack_int* ihi,
					double* lscale, double* rscale, double* abnrm,
					double* bbnrm, double* rconde, double* rcondv, double* work,
					lapack_int* lwork, lapack_int* iwork, lapack_logical* bwork,
					lapack_int *info );
void LAPACK_cggevx( char* balanc, char* jobvl, char* jobvr, char* sense,
					lapack_int* n, lapack_complex_float* a, lapack_int* lda,
					lapack_complex_float* b, lapack_int* ldb,
					lapack_complex_float* alpha, lapack_complex_float* beta,
					lapack_complex_float* vl, lapack_int* ldvl,
					lapack_complex_float* vr, lapack_int* ldvr, lapack_int* ilo,
					lapack_int* ihi, float* lscale, float* rscale, float* abnrm,
					float* bbnrm, float* rconde, float* rcondv,
					lapack_complex_float* work, lapack_int* lwork, float* rwork,
					lapack_int* iwork, lapack_logical* bwork,
					lapack_int *info );
void LAPACK_zggevx( char* balanc, char* jobvl, char* jobvr, char* sense,
					lapack_int* n, lapack_complex_double* a, lapack_int* lda,
					lapack_complex_double* b, lapack_int* ldb,
					lapack_complex_double* alpha, lapack_complex_double* beta,
					lapack_complex_double* vl, lapack_int* ldvl,
					lapack_complex_double* vr, lapack_int* ldvr,
					lapack_int* ilo, lapack_int* ihi, double* lscale,
					double* rscale, double* abnrm, double* bbnrm,
					double* rconde, double* rcondv, lapack_complex_double* work,
					lapack_int* lwork, double* rwork, lapack_int* iwork,
					lapack_logical* bwork, lapack_int *info );
void LAPACK_dsfrk( char* transr, char* uplo, char* trans, lapack_int* n,
				lapack_int* k, double* alpha, const double* a,
				lapack_int* lda, double* beta, double* c );
void LAPACK_ssfrk( char* transr, char* uplo, char* trans, lapack_int* n,
				lapack_int* k, float* alpha, const float* a, lapack_int* lda,
				float* beta, float* c );
void LAPACK_zhfrk( char* transr, char* uplo, char* trans, lapack_int* n,
				lapack_int* k, double* alpha, const lapack_complex_double* a,
				lapack_int* lda, double* beta, lapack_complex_double* c );
void LAPACK_chfrk( char* transr, char* uplo, char* trans, lapack_int* n,
				lapack_int* k, float* alpha, const lapack_complex_float* a,
				lapack_int* lda, float* beta, lapack_complex_float* c );
void LAPACK_dtfsm( char* transr, char* side, char* uplo, char* trans,
				char* diag, lapack_int* m, lapack_int* n, double* alpha,
				const double* a, double* b, lapack_int* ldb );
void LAPACK_stfsm( char* transr, char* side, char* uplo, char* trans,
				char* diag, lapack_int* m, lapack_int* n, float* alpha,
				const float* a, float* b, lapack_int* ldb );
void LAPACK_ztfsm( char* transr, char* side, char* uplo, char* trans,
				char* diag, lapack_int* m, lapack_int* n,
				lapack_complex_double* alpha, const lapack_complex_double* a,
				lapack_complex_double* b, lapack_int* ldb );
void LAPACK_ctfsm( char* transr, char* side, char* uplo, char* trans,
				char* diag, lapack_int* m, lapack_int* n,
				lapack_complex_float* alpha, const lapack_complex_float* a,
				lapack_complex_float* b, lapack_int* ldb );
void LAPACK_dtfttp( char* transr, char* uplo, lapack_int* n, const double* arf,
					double* ap, lapack_int *info );
void LAPACK_stfttp( char* transr, char* uplo, lapack_int* n, const float* arf,
					float* ap, lapack_int *info );
void LAPACK_ztfttp( char* transr, char* uplo, lapack_int* n,
					const lapack_complex_double* arf, lapack_complex_double* ap,
					lapack_int *info );
void LAPACK_ctfttp( char* transr, char* uplo, lapack_int* n,
					const lapack_complex_float* arf, lapack_complex_float* ap,
					lapack_int *info );
void LAPACK_dtfttr( char* transr, char* uplo, lapack_int* n, const double* arf,
					double* a, lapack_int* lda, lapack_int *info );
void LAPACK_stfttr( char* transr, char* uplo, lapack_int* n, const float* arf,
					float* a, lapack_int* lda, lapack_int *info );
void LAPACK_ztfttr( char* transr, char* uplo, lapack_int* n,
					const lapack_complex_double* arf, lapack_complex_double* a,
					lapack_int* lda, lapack_int *info );
void LAPACK_ctfttr( char* transr, char* uplo, lapack_int* n,
					const lapack_complex_float* arf, lapack_complex_float* a,
					lapack_int* lda, lapack_int *info );
void LAPACK_dtpttf( char* transr, char* uplo, lapack_int* n, const double* ap,
					double* arf, lapack_int *info );
void LAPACK_stpttf( char* transr, char* uplo, lapack_int* n, const float* ap,
					float* arf, lapack_int *info );
void LAPACK_ztpttf( char* transr, char* uplo, lapack_int* n,
					const lapack_complex_double* ap, lapack_complex_double* arf,
					lapack_int *info );
void LAPACK_ctpttf( char* transr, char* uplo, lapack_int* n,
					const lapack_complex_float* ap, lapack_complex_float* arf,
					lapack_int *info );
void LAPACK_dtpttr( char* uplo, lapack_int* n, const double* ap, double* a,
					lapack_int* lda, lapack_int *info );
void LAPACK_stpttr( char* uplo, lapack_int* n, const float* ap, float* a,
					lapack_int* lda, lapack_int *info );
void LAPACK_ztpttr( char* uplo, lapack_int* n, const lapack_complex_double* ap,
					lapack_complex_double* a, lapack_int* lda,
					lapack_int *info );
void LAPACK_ctpttr( char* uplo, lapack_int* n, const lapack_complex_float* ap,
					lapack_complex_float* a, lapack_int* lda,
					lapack_int *info );
void LAPACK_dtrttf( char* transr, char* uplo, lapack_int* n, const double* a,
					lapack_int* lda, double* arf, lapack_int *info );
void LAPACK_strttf( char* transr, char* uplo, lapack_int* n, const float* a,
					lapack_int* lda, float* arf, lapack_int *info );
void LAPACK_ztrttf( char* transr, char* uplo, lapack_int* n,
					const lapack_complex_double* a, lapack_int* lda,
					lapack_complex_double* arf, lapack_int *info );
void LAPACK_ctrttf( char* transr, char* uplo, lapack_int* n,
					const lapack_complex_float* a, lapack_int* lda,
					lapack_complex_float* arf, lapack_int *info );
void LAPACK_dtrttp( char* uplo, lapack_int* n, const double* a, lapack_int* lda,
					double* ap, lapack_int *info );
void LAPACK_strttp( char* uplo, lapack_int* n, const float* a, lapack_int* lda,
					float* ap, lapack_int *info );
void LAPACK_ztrttp( char* uplo, lapack_int* n, const lapack_complex_double* a,
					lapack_int* lda, lapack_complex_double* ap,
					lapack_int *info );
void LAPACK_ctrttp( char* uplo, lapack_int* n, const lapack_complex_float* a,
					lapack_int* lda, lapack_complex_float* ap,
					lapack_int *info );
void LAPACK_sgeqrfp( lapack_int* m, lapack_int* n, float* a, lapack_int* lda,
					float* tau, float* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_dgeqrfp( lapack_int* m, lapack_int* n, double* a, lapack_int* lda,
					double* tau, double* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_cgeqrfp( lapack_int* m, lapack_int* n, lapack_complex_float* a,
					lapack_int* lda, lapack_complex_float* tau,
					lapack_complex_float* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_zgeqrfp( lapack_int* m, lapack_int* n, lapack_complex_double* a,
					lapack_int* lda, lapack_complex_double* tau,
					lapack_complex_double* work, lapack_int* lwork,
					lapack_int *info );
void LAPACK_clacgv( lapack_int* n, lapack_complex_float* x, lapack_int* incx );
void LAPACK_zlacgv( lapack_int* n, lapack_complex_double* x, lapack_int* incx );
void LAPACK_slarnv( lapack_int* idist, lapack_int* iseed, lapack_int* n,
					float* x );
void LAPACK_dlarnv( lapack_int* idist, lapack_int* iseed, lapack_int* n,
					double* x );
void LAPACK_clarnv( lapack_int* idist, lapack_int* iseed, lapack_int* n,
					lapack_complex_float* x );
void LAPACK_zlarnv( lapack_int* idist, lapack_int* iseed, lapack_int* n,
					lapack_complex_double* x );
void LAPACK_sgeqr2( lapack_int* m, lapack_int* n, float* a, lapack_int* lda,
					float* tau, float* work, lapack_int *info );
void LAPACK_dgeqr2( lapack_int* m, lapack_int* n, double* a, lapack_int* lda,
					double* tau, double* work, lapack_int *info );
void LAPACK_cgeqr2( lapack_int* m, lapack_int* n, lapack_complex_float* a,
					lapack_int* lda, lapack_complex_float* tau,
					lapack_complex_float* work, lapack_int *info );
void LAPACK_zgeqr2( lapack_int* m, lapack_int* n, lapack_complex_double* a,
					lapack_int* lda, lapack_complex_double* tau,
					lapack_complex_double* work, lapack_int *info );
void LAPACK_slacn2( lapack_int* n, float* v, float* x, lapack_int* isgn,
					float* est, lapack_int* kase, lapack_int* isave );
void LAPACK_dlacn2( lapack_int* n, double* v, double* x, lapack_int* isgn,
					double* est, lapack_int* kase, lapack_int* isave );
void LAPACK_clacn2( lapack_int* n, lapack_complex_float* v,
					lapack_complex_float* x, float* est,
					lapack_int* kase, lapack_int* isave );
void LAPACK_zlacn2( lapack_int* n, lapack_complex_double* v,
					lapack_complex_double* x, double* est,
					lapack_int* kase, lapack_int* isave );
void LAPACK_slacpy( char* uplo, lapack_int* m, lapack_int* n, const float* a,
					lapack_int* lda, float* b, lapack_int* ldb );
void LAPACK_dlacpy( char* uplo, lapack_int* m, lapack_int* n, const double* a,
					lapack_int* lda, double* b, lapack_int* ldb );
void LAPACK_clacpy( char* uplo, lapack_int* m, lapack_int* n,
					const lapack_complex_float* a, lapack_int* lda,
					lapack_complex_float* b, lapack_int* ldb );
void LAPACK_zlacpy( char* uplo, lapack_int* m, lapack_int* n,
					const lapack_complex_double* a, lapack_int* lda,
					lapack_complex_double* b, lapack_int* ldb );

void LAPACK_clacp2( char* uplo, lapack_int* m, lapack_int* n, const float* a,
					lapack_int* lda, lapack_complex_float* b, lapack_int* ldb );
void LAPACK_zlacp2( char* uplo, lapack_int* m, lapack_int* n, const double* a,
					lapack_int* lda, lapack_complex_double* b,
					lapack_int* ldb );

void LAPACK_sgetf2( lapack_int* m, lapack_int* n, float* a, lapack_int* lda,
					lapack_int* ipiv, lapack_int *info );
void LAPACK_dgetf2( lapack_int* m, lapack_int* n, double* a, lapack_int* lda,
					lapack_int* ipiv, lapack_int *info );
void LAPACK_cgetf2( lapack_int* m, lapack_int* n, lapack_complex_float* a,
					lapack_int* lda, lapack_int* ipiv, lapack_int *info );
void LAPACK_zgetf2( lapack_int* m, lapack_int* n, lapack_complex_double* a,
					lapack_int* lda, lapack_int* ipiv, lapack_int *info );
void LAPACK_slaswp( lapack_int* n, float* a, lapack_int* lda, lapack_int* k1,
					lapack_int* k2, const lapack_int* ipiv, lapack_int* incx );
void LAPACK_dlaswp( lapack_int* n, double* a, lapack_int* lda, lapack_int* k1,
					lapack_int* k2, const lapack_int* ipiv, lapack_int* incx );
void LAPACK_claswp( lapack_int* n, lapack_complex_float* a, lapack_int* lda,
					lapack_int* k1, lapack_int* k2, const lapack_int* ipiv,
					lapack_int* incx );
void LAPACK_zlaswp( lapack_int* n, lapack_complex_double* a, lapack_int* lda,
					lapack_int* k1, lapack_int* k2, const lapack_int* ipiv,
					lapack_int* incx );
float LAPACK_slange( char* norm, lapack_int* m, lapack_int* n, const float* a,
					lapack_int* lda, float* work );
double LAPACK_dlange( char* norm, lapack_int* m, lapack_int* n, const double* a,
					lapack_int* lda, double* work );
float LAPACK_clange( char* norm, lapack_int* m, lapack_int* n,
					const lapack_complex_float* a, lapack_int* lda, float* work );
double LAPACK_zlange( char* norm, lapack_int* m, lapack_int* n,
					const lapack_complex_double* a, lapack_int* lda, double* work );
float LAPACK_clanhe( char* norm, char* uplo, lapack_int* n,
					const lapack_complex_float* a, lapack_int* lda, float* work );
double LAPACK_zlanhe( char* norm, char* uplo, lapack_int* n,
					const lapack_complex_double* a, lapack_int* lda, double* work );
float LAPACK_slansy( char* norm, char* uplo, lapack_int* n, const float* a,
					lapack_int* lda, float* work );
double LAPACK_dlansy( char* norm, char* uplo, lapack_int* n, const double* a,
					lapack_int* lda, double* work );
float LAPACK_clansy( char* norm, char* uplo, lapack_int* n,
					const lapack_complex_float* a, lapack_int* lda, float* work );
double LAPACK_zlansy( char* norm, char* uplo, lapack_int* n,
					const lapack_complex_double* a, lapack_int* lda, double* work );
float LAPACK_slantr( char* norm, char* uplo, char* diag, lapack_int* m,
					lapack_int* n, const float* a, lapack_int* lda, float* work );
double LAPACK_dlantr( char* norm, char* uplo, char* diag, lapack_int* m,
					lapack_int* n, const double* a, lapack_int* lda, double* work );
float LAPACK_clantr( char* norm, char* uplo, char* diag, lapack_int* m,
					lapack_int* n, const lapack_complex_float* a, lapack_int* lda,
					float* work );
double LAPACK_zlantr( char* norm, char* uplo, char* diag, lapack_int* m,
					lapack_int* n, const lapack_complex_double* a, lapack_int* lda,
					double* work );
float LAPACK_slamch( char* cmach );
double LAPACK_dlamch( char* cmach );
void LAPACK_sgelq2( lapack_int* m, lapack_int* n, float* a, lapack_int* lda,
					float* tau, float* work, lapack_int *info );
void LAPACK_dgelq2( lapack_int* m, lapack_int* n, double* a, lapack_int* lda,
					double* tau, double* work, lapack_int *info );
void LAPACK_cgelq2( lapack_int* m, lapack_int* n, lapack_complex_float* a,
					lapack_int* lda, lapack_complex_float* tau,
					lapack_complex_float* work, lapack_int *info );
void LAPACK_zgelq2( lapack_int* m, lapack_int* n, lapack_complex_double* a,
					lapack_int* lda, lapack_complex_double* tau,
					lapack_complex_double* work, lapack_int *info );
void LAPACK_slarfb( char* side, char* trans, char* direct, char* storev,
					lapack_int* m, lapack_int* n, lapack_int* k, const float* v,
					lapack_int* ldv, const float* t, lapack_int* ldt, float* c,
					lapack_int* ldc, float* work, lapack_int* ldwork );
void LAPACK_dlarfb( char* side, char* trans, char* direct, char* storev,
					lapack_int* m, lapack_int* n, lapack_int* k,
					const double* v, lapack_int* ldv, const double* t,
					lapack_int* ldt, double* c, lapack_int* ldc, double* work,
					lapack_int* ldwork );
void LAPACK_clarfb( char* side, char* trans, char* direct, char* storev,
					lapack_int* m, lapack_int* n, lapack_int* k,
					const lapack_complex_float* v, lapack_int* ldv,
					const lapack_complex_float* t, lapack_int* ldt,
					lapack_complex_float* c, lapack_int* ldc,
					lapack_complex_float* work, lapack_int* ldwork );
void LAPACK_zlarfb( char* side, char* trans, char* direct, char* storev,
					lapack_int* m, lapack_int* n, lapack_int* k,
					const lapack_complex_double* v, lapack_int* ldv,
					const lapack_complex_double* t, lapack_int* ldt,
					lapack_complex_double* c, lapack_int* ldc,
					lapack_complex_double* work, lapack_int* ldwork );
void LAPACK_slarfg( lapack_int* n, float* alpha, float* x, lapack_int* incx,
					float* tau );
void LAPACK_dlarfg( lapack_int* n, double* alpha, double* x, lapack_int* incx,
					double* tau );
void LAPACK_clarfg( lapack_int* n, lapack_complex_float* alpha,
					lapack_complex_float* x, lapack_int* incx,
					lapack_complex_float* tau );
void LAPACK_zlarfg( lapack_int* n, lapack_complex_double* alpha,
					lapack_complex_double* x, lapack_int* incx,
					lapack_complex_double* tau );
void LAPACK_slarft( char* direct, char* storev, lapack_int* n, lapack_int* k,
					const float* v, lapack_int* ldv, const float* tau, float* t,
					lapack_int* ldt );
void LAPACK_dlarft( char* direct, char* storev, lapack_int* n, lapack_int* k,
					const double* v, lapack_int* ldv, const double* tau,
					double* t, lapack_int* ldt );
void LAPACK_clarft( char* direct, char* storev, lapack_int* n, lapack_int* k,
					const lapack_complex_float* v, lapack_int* ldv,
					const lapack_complex_float* tau, lapack_complex_float* t,
					lapack_int* ldt );
void LAPACK_zlarft( char* direct, char* storev, lapack_int* n, lapack_int* k,
					const lapack_complex_double* v, lapack_int* ldv,
					const lapack_complex_double* tau, lapack_complex_double* t,
					lapack_int* ldt );
void LAPACK_slarfx( char* side, lapack_int* m, lapack_int* n, const float* v,
					float* tau, float* c, lapack_int* ldc, float* work );
void LAPACK_dlarfx( char* side, lapack_int* m, lapack_int* n, const double* v,
					double* tau, double* c, lapack_int* ldc, double* work );
void LAPACK_clarfx( char* side, lapack_int* m, lapack_int* n,
					const lapack_complex_float* v, lapack_complex_float* tau,
					lapack_complex_float* c, lapack_int* ldc,
					lapack_complex_float* work );
void LAPACK_zlarfx( char* side, lapack_int* m, lapack_int* n,
					const lapack_complex_double* v, lapack_complex_double* tau,
					lapack_complex_double* c, lapack_int* ldc,
					lapack_complex_double* work );
void LAPACK_slatms( lapack_int* m, lapack_int* n, char* dist, lapack_int* iseed,
					char* sym, float* d, lapack_int* mode, float* cond,
					float* dmax, lapack_int* kl, lapack_int* ku, char* pack,
					float* a, lapack_int* lda, float* work, lapack_int *info );
void LAPACK_dlatms( lapack_int* m, lapack_int* n, char* dist, lapack_int* iseed,
					char* sym, double* d, lapack_int* mode, double* cond,
					double* dmax, lapack_int* kl, lapack_int* ku, char* pack,
					double* a, lapack_int* lda, double* work,
					lapack_int *info );
void LAPACK_clatms( lapack_int* m, lapack_int* n, char* dist, lapack_int* iseed,
					char* sym, float* d, lapack_int* mode, float* cond,
					float* dmax, lapack_int* kl, lapack_int* ku, char* pack,
					lapack_complex_float* a, lapack_int* lda,
					lapack_complex_float* work, lapack_int *info );
void LAPACK_zlatms( lapack_int* m, lapack_int* n, char* dist, lapack_int* iseed,
					char* sym, double* d, lapack_int* mode, double* cond,
					double* dmax, lapack_int* kl, lapack_int* ku, char* pack,
					lapack_complex_double* a, lapack_int* lda,
					lapack_complex_double* work, lapack_int *info );
void LAPACK_slag2d( lapack_int* m, lapack_int* n, const float* sa,
					lapack_int* ldsa, double* a, lapack_int* lda,
					lapack_int *info );
void LAPACK_dlag2s( lapack_int* m, lapack_int* n, const double* a,
					lapack_int* lda, float* sa, lapack_int* ldsa,
					lapack_int *info );
void LAPACK_clag2z( lapack_int* m, lapack_int* n,
					const lapack_complex_float* sa, lapack_int* ldsa,
					lapack_complex_double* a, lapack_int* lda,
					lapack_int *info );
void LAPACK_zlag2c( lapack_int* m, lapack_int* n,
					const lapack_complex_double* a, lapack_int* lda,
					lapack_complex_float* sa, lapack_int* ldsa,
					lapack_int *info );
void LAPACK_slauum( char* uplo, lapack_int* n, float* a, lapack_int* lda,
					lapack_int *info );
void LAPACK_dlauum( char* uplo, lapack_int* n, double* a, lapack_int* lda,
					lapack_int *info );
void LAPACK_clauum( char* uplo, lapack_int* n, lapack_complex_float* a,
					lapack_int* lda, lapack_int *info );
void LAPACK_zlauum( char* uplo, lapack_int* n, lapack_complex_double* a,
					lapack_int* lda, lapack_int *info );
void LAPACK_slagge( lapack_int* m, lapack_int* n, lapack_int* kl,
					lapack_int* ku, const float* d, float* a, lapack_int* lda,
					lapack_int* iseed, float* work, lapack_int *info );
void LAPACK_dlagge( lapack_int* m, lapack_int* n, lapack_int* kl,
					lapack_int* ku, const double* d, double* a, lapack_int* lda,
					lapack_int* iseed, double* work, lapack_int *info );
void LAPACK_clagge( lapack_int* m, lapack_int* n, lapack_int* kl,
					lapack_int* ku, const float* d, lapack_complex_float* a,
					lapack_int* lda, lapack_int* iseed,
					lapack_complex_float* work, lapack_int *info );
void LAPACK_zlagge( lapack_int* m, lapack_int* n, lapack_int* kl,
					lapack_int* ku, const double* d, lapack_complex_double* a,
					lapack_int* lda, lapack_int* iseed,
					lapack_complex_double* work, lapack_int *info );
void LAPACK_slaset( char* uplo, lapack_int* m, lapack_int* n, float* alpha,
					float* beta, float* a, lapack_int* lda );
void LAPACK_dlaset( char* uplo, lapack_int* m, lapack_int* n, double* alpha,
					double* beta, double* a, lapack_int* lda );
void LAPACK_claset( char* uplo, lapack_int* m, lapack_int* n,
					lapack_complex_float* alpha, lapack_complex_float* beta,
					lapack_complex_float* a, lapack_int* lda );
void LAPACK_zlaset( char* uplo, lapack_int* m, lapack_int* n,
					lapack_complex_double* alpha, lapack_complex_double* beta,
					lapack_complex_double* a, lapack_int* lda );
void LAPACK_slasrt( char* id, lapack_int* n, float* d, lapack_int *info );
void LAPACK_dlasrt( char* id, lapack_int* n, double* d, lapack_int *info );
void LAPACK_claghe( lapack_int* n, lapack_int* k, const float* d,
					lapack_complex_float* a, lapack_int* lda, lapack_int* iseed,
					lapack_complex_float* work, lapack_int *info );
void LAPACK_zlaghe( lapack_int* n, lapack_int* k, const double* d,
					lapack_complex_double* a, lapack_int* lda,
					lapack_int* iseed, lapack_complex_double* work,
					lapack_int *info );
void LAPACK_slagsy( lapack_int* n, lapack_int* k, const float* d, float* a,
					lapack_int* lda, lapack_int* iseed, float* work,
					lapack_int *info );
void LAPACK_dlagsy( lapack_int* n, lapack_int* k, const double* d, double* a,
					lapack_int* lda, lapack_int* iseed, double* work,
					lapack_int *info );
void LAPACK_clagsy( lapack_int* n, lapack_int* k, const float* d,
					lapack_complex_float* a, lapack_int* lda, lapack_int* iseed,
					lapack_complex_float* work, lapack_int *info );
void LAPACK_zlagsy( lapack_int* n, lapack_int* k, const double* d,
					lapack_complex_double* a, lapack_int* lda,
					lapack_int* iseed, lapack_complex_double* work,
					lapack_int *info );
void LAPACK_slapmr( lapack_logical* forwrd, lapack_int* m, lapack_int* n,
					float* x, lapack_int* ldx, lapack_int* k );
void LAPACK_dlapmr( lapack_logical* forwrd, lapack_int* m, lapack_int* n,
					double* x, lapack_int* ldx, lapack_int* k );
void LAPACK_clapmr( lapack_logical* forwrd, lapack_int* m, lapack_int* n,
					lapack_complex_float* x, lapack_int* ldx, lapack_int* k );
void LAPACK_zlapmr( lapack_logical* forwrd, lapack_int* m, lapack_int* n,
					lapack_complex_double* x, lapack_int* ldx, lapack_int* k );
float LAPACK_slapy2( float* x, float* y );
double LAPACK_dlapy2( double* x, double* y );
float LAPACK_slapy3( float* x, float* y, float* z );
double LAPACK_dlapy3( double* x, double* y, double* z );
void LAPACK_slartgp( float* f, float* g, float* cs, float* sn, float* r );
void LAPACK_dlartgp( double* f, double* g, double* cs, double* sn, double* r );
void LAPACK_slartgs( float* x, float* y, float* sigma, float* cs, float* sn );
void LAPACK_dlartgs( double* x, double* y, double* sigma, double* cs,
					double* sn );
// LAPACK 3.3.0
void LAPACK_cbbcsd( char* jobu1, char* jobu2,
					char* jobv1t, char* jobv2t, char* trans,
					lapack_int* m, lapack_int* p, lapack_int* q,
					float* theta, float* phi,
					lapack_complex_float* u1, lapack_int* ldu1,
					lapack_complex_float* u2, lapack_int* ldu2,
					lapack_complex_float* v1t, lapack_int* ldv1t,
					lapack_complex_float* v2t, lapack_int* ldv2t,
					float* b11d, float* b11e, float* b12d,
					float* b12e, float* b21d, float* b21e,
					float* b22d, float* b22e, float* rwork,
					lapack_int* lrwork , lapack_int *info );
void LAPACK_cheswapr( char* uplo, lapack_int* n,
					lapack_complex_float* a, lapack_int* i1,
					lapack_int* i2 );
void LAPACK_chetri2( char* uplo, lapack_int* n,
					lapack_complex_float* a, lapack_int* lda,
					const lapack_int* ipiv,
					lapack_complex_float* work, lapack_int* lwork , lapack_int *info );
void LAPACK_chetri2x( char* uplo, lapack_int* n,
					lapack_complex_float* a, lapack_int* lda,
					const lapack_int* ipiv,
					lapack_complex_float* work, lapack_int* nb , lapack_int *info );
void LAPACK_chetrs2( char* uplo, lapack_int* n,
					lapack_int* nrhs, const lapack_complex_float* a,
					lapack_int* lda, const lapack_int* ipiv,
					lapack_complex_float* b, lapack_int* ldb,
					lapack_complex_float* work , lapack_int *info );
void LAPACK_csyconv( char* uplo, char* way,
					lapack_int* n, lapack_complex_float* a,
					lapack_int* lda, const lapack_int* ipiv,
					lapack_complex_float* work , lapack_int *info );
void LAPACK_csyswapr( char* uplo, lapack_int* n,
					lapack_complex_float* a, lapack_int* i1,
					lapack_int* i2 );
void LAPACK_csytri2( char* uplo, lapack_int* n,
					lapack_complex_float* a, lapack_int* lda,
					const lapack_int* ipiv,
					lapack_complex_float* work, lapack_int* lwork , lapack_int *info );
void LAPACK_csytri2x( char* uplo, lapack_int* n,
					lapack_complex_float* a, lapack_int* lda,
					const lapack_int* ipiv,
					lapack_complex_float* work, lapack_int* nb , lapack_int *info );
void LAPACK_csytrs2( char* uplo, lapack_int* n,
					lapack_int* nrhs, const lapack_complex_float* a,
					lapack_int* lda, const lapack_int* ipiv,
					lapack_complex_float* b, lapack_int* ldb,
					lapack_complex_float* work , lapack_int *info );
void LAPACK_cunbdb( char* trans, char* signs,
					lapack_int* m, lapack_int* p, lapack_int* q,
					lapack_complex_float* x11, lapack_int* ldx11,
					lapack_complex_float* x12, lapack_int* ldx12,
					lapack_complex_float* x21, lapack_int* ldx21,
					lapack_complex_float* x22, lapack_int* ldx22,
					float* theta, float* phi,
					lapack_complex_float* taup1,
					lapack_complex_float* taup2,
					lapack_complex_float* tauq1,
					lapack_complex_float* tauq2,
					lapack_complex_float* work, lapack_int* lwork , lapack_int *info );
void LAPACK_cuncsd( char* jobu1, char* jobu2,
					char* jobv1t, char* jobv2t, char* trans,
					char* signs, lapack_int* m, lapack_int* p,
					lapack_int* q, lapack_complex_float* x11,
					lapack_int* ldx11, lapack_complex_float* x12,
					lapack_int* ldx12, lapack_complex_float* x21,
					lapack_int* ldx21, lapack_complex_float* x22,
					lapack_int* ldx22, float* theta,
					lapack_complex_float* u1, lapack_int* ldu1,
					lapack_complex_float* u2, lapack_int* ldu2,
					lapack_complex_float* v1t, lapack_int* ldv1t,
					lapack_complex_float* v2t, lapack_int* ldv2t,
					lapack_complex_float* work, lapack_int* lwork,
					float* rwork, lapack_int* lrwork,
					lapack_int* iwork , lapack_int *info );
void LAPACK_dbbcsd( char* jobu1, char* jobu2,
					char* jobv1t, char* jobv2t, char* trans,
					lapack_int* m, lapack_int* p, lapack_int* q,
					double* theta, double* phi, double* u1,
					lapack_int* ldu1, double* u2, lapack_int* ldu2,
					double* v1t, lapack_int* ldv1t, double* v2t,
					lapack_int* ldv2t, double* b11d, double* b11e,
					double* b12d, double* b12e, double* b21d,
					double* b21e, double* b22d, double* b22e,
					double* work, lapack_int* lwork , lapack_int *info );
void LAPACK_dorbdb( char* trans, char* signs,
					lapack_int* m, lapack_int* p, lapack_int* q,
					double* x11, lapack_int* ldx11, double* x12,
					lapack_int* ldx12, double* x21, lapack_int* ldx21,
					double* x22, lapack_int* ldx22, double* theta,
					double* phi, double* taup1, double* taup2,
					double* tauq1, double* tauq2, double* work,
					lapack_int* lwork , lapack_int *info );
void LAPACK_dorcsd( char* jobu1, char* jobu2,
					char* jobv1t, char* jobv2t, char* trans,
					char* signs, lapack_int* m, lapack_int* p,
					lapack_int* q, double* x11, lapack_int* ldx11,
					double* x12, lapack_int* ldx12, double* x21,
					lapack_int* ldx21, double* x22, lapack_int* ldx22,
					double* theta, double* u1, lapack_int* ldu1,
					double* u2, lapack_int* ldu2, double* v1t,
					lapack_int* ldv1t, double* v2t, lapack_int* ldv2t,
					double* work, lapack_int* lwork,
					lapack_int* iwork , lapack_int *info );
void LAPACK_dsyconv( char* uplo, char* way,
					lapack_int* n, double* a, lapack_int* lda,
					const lapack_int* ipiv, double* work , lapack_int *info );
void LAPACK_dsyswapr( char* uplo, lapack_int* n,
					double* a, lapack_int* i1, lapack_int* i2 );
void LAPACK_dsytri2( char* uplo, lapack_int* n,
					double* a, lapack_int* lda,
					const lapack_int* ipiv,
					lapack_complex_double* work, lapack_int* lwork , lapack_int *info );
void LAPACK_dsytri2x( char* uplo, lapack_int* n,
					double* a, lapack_int* lda,
					const lapack_int* ipiv, double* work,
					lapack_int* nb , lapack_int *info );
void LAPACK_dsytrs2( char* uplo, lapack_int* n,
					lapack_int* nrhs, const double* a,
					lapack_int* lda, const lapack_int* ipiv,
					double* b, lapack_int* ldb, double* work , lapack_int *info );
void LAPACK_sbbcsd( char* jobu1, char* jobu2,
					char* jobv1t, char* jobv2t, char* trans,
					lapack_int* m, lapack_int* p, lapack_int* q,
					float* theta, float* phi, float* u1,
					lapack_int* ldu1, float* u2, lapack_int* ldu2,
					float* v1t, lapack_int* ldv1t, float* v2t,
					lapack_int* ldv2t, float* b11d, float* b11e,
					float* b12d, float* b12e, float* b21d,
					float* b21e, float* b22d, float* b22e,
					float* work, lapack_int* lwork , lapack_int *info );
void LAPACK_sorbdb( char* trans, char* signs,
					lapack_int* m, lapack_int* p, lapack_int* q,
					float* x11, lapack_int* ldx11, float* x12,
					lapack_int* ldx12, float* x21, lapack_int* ldx21,
					float* x22, lapack_int* ldx22, float* theta,
					float* phi, float* taup1, float* taup2,
					float* tauq1, float* tauq2, float* work,
					lapack_int* lwork , lapack_int *info );
void LAPACK_sorcsd( char* jobu1, char* jobu2,
					char* jobv1t, char* jobv2t, char* trans,
					char* signs, lapack_int* m, lapack_int* p,
					lapack_int* q, float* x11, lapack_int* ldx11,
					float* x12, lapack_int* ldx12, float* x21,
					lapack_int* ldx21, float* x22, lapack_int* ldx22,
					float* theta, float* u1, lapack_int* ldu1,
					float* u2, lapack_int* ldu2, float* v1t,
					lapack_int* ldv1t, float* v2t, lapack_int* ldv2t,
					float* work, lapack_int* lwork,
					lapack_int* iwork , lapack_int *info );
void LAPACK_ssyconv( char* uplo, char* way,
					lapack_int* n, float* a, lapack_int* lda,
					const lapack_int* ipiv, float* work , lapack_int *info );
void LAPACK_ssyswapr( char* uplo, lapack_int* n,
					float* a, lapack_int* i1, lapack_int* i2 );
void LAPACK_ssytri2( char* uplo, lapack_int* n,
					float* a, lapack_int* lda,
					const lapack_int* ipiv,
					lapack_complex_float* work, lapack_int* lwork , lapack_int *info );
void LAPACK_ssytri2x( char* uplo, lapack_int* n,
					float* a, lapack_int* lda,
					const lapack_int* ipiv, float* work,
					lapack_int* nb , lapack_int *info );
void LAPACK_ssytrs2( char* uplo, lapack_int* n,
					lapack_int* nrhs, const float* a,
					lapack_int* lda, const lapack_int* ipiv,
					float* b, lapack_int* ldb, float* work , lapack_int *info );
void LAPACK_zbbcsd( char* jobu1, char* jobu2,
					char* jobv1t, char* jobv2t, char* trans,
					lapack_int* m, lapack_int* p, lapack_int* q,
					double* theta, double* phi,
					lapack_complex_double* u1, lapack_int* ldu1,
					lapack_complex_double* u2, lapack_int* ldu2,
					lapack_complex_double* v1t, lapack_int* ldv1t,
					lapack_complex_double* v2t, lapack_int* ldv2t,
					double* b11d, double* b11e, double* b12d,
					double* b12e, double* b21d, double* b21e,
					double* b22d, double* b22e, double* rwork,
					lapack_int* lrwork , lapack_int *info );
void LAPACK_zheswapr( char* uplo, lapack_int* n,
					lapack_complex_double* a, lapack_int* i1,
					lapack_int* i2 );
void LAPACK_zhetri2( char* uplo, lapack_int* n,
					lapack_complex_double* a, lapack_int* lda,
					const lapack_int* ipiv,
					lapack_complex_double* work, lapack_int* lwork , lapack_int *info );
void LAPACK_zhetri2x( char* uplo, lapack_int* n,
					lapack_complex_double* a, lapack_int* lda,
					const lapack_int* ipiv,
					lapack_complex_double* work, lapack_int* nb , lapack_int *info );
void LAPACK_zhetrs2( char* uplo, lapack_int* n,
					lapack_int* nrhs,
					const lapack_complex_double* a, lapack_int* lda,
					const lapack_int* ipiv,
					lapack_complex_double* b, lapack_int* ldb,
					lapack_complex_double* work , lapack_int *info );
void LAPACK_zsyconv( char* uplo, char* way,
					lapack_int* n, lapack_complex_double* a,
					lapack_int* lda, const lapack_int* ipiv,
					lapack_complex_double* work , lapack_int *info );
void LAPACK_zsyswapr( char* uplo, lapack_int* n,
					lapack_complex_double* a, lapack_int* i1,
					lapack_int* i2 );
void LAPACK_zsytri2( char* uplo, lapack_int* n,
					lapack_complex_double* a, lapack_int* lda,
					const lapack_int* ipiv,
					lapack_complex_double* work, lapack_int* lwork , lapack_int *info );
void LAPACK_zsytri2x( char* uplo, lapack_int* n,
					lapack_complex_double* a, lapack_int* lda,
					const lapack_int* ipiv,
					lapack_complex_double* work, lapack_int* nb , lapack_int *info );
void LAPACK_zsytrs2( char* uplo, lapack_int* n,
					lapack_int* nrhs,
					const lapack_complex_double* a, lapack_int* lda,
					const lapack_int* ipiv,
					lapack_complex_double* b, lapack_int* ldb,
					lapack_complex_double* work , lapack_int *info );
void LAPACK_zunbdb( char* trans, char* signs,
					lapack_int* m, lapack_int* p, lapack_int* q,
					lapack_complex_double* x11, lapack_int* ldx11,
					lapack_complex_double* x12, lapack_int* ldx12,
					lapack_complex_double* x21, lapack_int* ldx21,
					lapack_complex_double* x22, lapack_int* ldx22,
					double* theta, double* phi,
					lapack_complex_double* taup1,
					lapack_complex_double* taup2,
					lapack_complex_double* tauq1,
					lapack_complex_double* tauq2,
					lapack_complex_double* work, lapack_int* lwork , lapack_int *info );
void LAPACK_zuncsd( char* jobu1, char* jobu2,
					char* jobv1t, char* jobv2t, char* trans,
					char* signs, lapack_int* m, lapack_int* p,
					lapack_int* q, lapack_complex_double* x11,
					lapack_int* ldx11, lapack_complex_double* x12,
					lapack_int* ldx12, lapack_complex_double* x21,
					lapack_int* ldx21, lapack_complex_double* x22,
					lapack_int* ldx22, double* theta,
					lapack_complex_double* u1, lapack_int* ldu1,
					lapack_complex_double* u2, lapack_int* ldu2,
					lapack_complex_double* v1t, lapack_int* ldv1t,
					lapack_complex_double* v2t, lapack_int* ldv2t,
					lapack_complex_double* work, lapack_int* lwork,
					double* rwork, lapack_int* lrwork,
					lapack_int* iwork , lapack_int *info );
// LAPACK 3.4.0
void LAPACK_sgemqrt( char* side, char* trans, lapack_int* m, lapack_int* n,
					lapack_int* k, lapack_int* nb, const float* v,
					lapack_int* ldv, const float* t, lapack_int* ldt, float* c,
					lapack_int* ldc, float* work, lapack_int *info );
void LAPACK_dgemqrt( char* side, char* trans, lapack_int* m, lapack_int* n,
					lapack_int* k, lapack_int* nb, const double* v,
					lapack_int* ldv, const double* t, lapack_int* ldt,
					double* c, lapack_int* ldc, double* work,
					lapack_int *info );
void LAPACK_cgemqrt( char* side, char* trans, lapack_int* m, lapack_int* n,
					lapack_int* k, lapack_int* nb,
					const lapack_complex_float* v, lapack_int* ldv,
					const lapack_complex_float* t, lapack_int* ldt,
					lapack_complex_float* c, lapack_int* ldc,
					lapack_complex_float* work, lapack_int *info );
void LAPACK_zgemqrt( char* side, char* trans, lapack_int* m, lapack_int* n,
					lapack_int* k, lapack_int* nb,
					const lapack_complex_double* v, lapack_int* ldv,
					const lapack_complex_double* t, lapack_int* ldt,
					lapack_complex_double* c, lapack_int* ldc,
					lapack_complex_double* work, lapack_int *info );
void LAPACK_sgeqrt( lapack_int* m, lapack_int* n, lapack_int* nb, float* a,
					lapack_int* lda, float* t, lapack_int* ldt, float* work,
					lapack_int *info );
void LAPACK_dgeqrt( lapack_int* m, lapack_int* n, lapack_int* nb, double* a,
					lapack_int* lda, double* t, lapack_int* ldt, double* work,
					lapack_int *info );
void LAPACK_cgeqrt( lapack_int* m, lapack_int* n, lapack_int* nb,
					lapack_complex_float* a, lapack_int* lda,
					lapack_complex_float* t, lapack_int* ldt,
					lapack_complex_float* work, lapack_int *info );
void LAPACK_zgeqrt( lapack_int* m, lapack_int* n, lapack_int* nb,
					lapack_complex_double* a, lapack_int* lda,
					lapack_complex_double* t, lapack_int* ldt,
					lapack_complex_double* work, lapack_int *info );
void LAPACK_sgeqrt2( lapack_int* m, lapack_int* n, float* a, lapack_int* lda,
					float* t, lapack_int* ldt, lapack_int *info );
void LAPACK_dgeqrt2( lapack_int* m, lapack_int* n, double* a, lapack_int* lda,
					double* t, lapack_int* ldt, lapack_int *info );
void LAPACK_cgeqrt2( lapack_int* m, lapack_int* n, lapack_complex_float* a,
					lapack_int* lda, lapack_complex_float* t, lapack_int* ldt,
					lapack_int *info );
void LAPACK_zgeqrt2( lapack_int* m, lapack_int* n, lapack_complex_double* a,
					lapack_int* lda, lapack_complex_double* t, lapack_int* ldt,
					lapack_int *info );
void LAPACK_sgeqrt3( lapack_int* m, lapack_int* n, float* a, lapack_int* lda,
					float* t, lapack_int* ldt, lapack_int *info );
void LAPACK_dgeqrt3( lapack_int* m, lapack_int* n, double* a, lapack_int* lda,
					double* t, lapack_int* ldt, lapack_int *info );
void LAPACK_cgeqrt3( lapack_int* m, lapack_int* n, lapack_complex_float* a,
					lapack_int* lda, lapack_complex_float* t, lapack_int* ldt,
					lapack_int *info );
void LAPACK_zgeqrt3( lapack_int* m, lapack_int* n, lapack_complex_double* a,
					lapack_int* lda, lapack_complex_double* t, lapack_int* ldt,
					lapack_int *info );
void LAPACK_stpmqrt( char* side, char* trans, lapack_int* m, lapack_int* n,
					lapack_int* k, lapack_int* l, lapack_int* nb,
					const float* v, lapack_int* ldv, const float* t,
					lapack_int* ldt, float* a, lapack_int* lda, float* b,
					lapack_int* ldb, float* work, lapack_int *info );
void LAPACK_dtpmqrt( char* side, char* trans, lapack_int* m, lapack_int* n,
					lapack_int* k, lapack_int* l, lapack_int* nb,
					const double* v, lapack_int* ldv, const double* t,
					lapack_int* ldt, double* a, lapack_int* lda, double* b,
					lapack_int* ldb, double* work, lapack_int *info );
void LAPACK_ctpmqrt( char* side, char* trans, lapack_int* m, lapack_int* n,
					lapack_int* k, lapack_int* l, lapack_int* nb,
					const lapack_complex_float* v, lapack_int* ldv,
					const lapack_complex_float* t, lapack_int* ldt,
					lapack_complex_float* a, lapack_int* lda,
					lapack_complex_float* b, lapack_int* ldb,
					lapack_complex_float* work, lapack_int *info );
void LAPACK_ztpmqrt( char* side, char* trans, lapack_int* m, lapack_int* n,
					lapack_int* k, lapack_int* l, lapack_int* nb,
					const lapack_complex_double* v, lapack_int* ldv,
					const lapack_complex_double* t, lapack_int* ldt,
					lapack_complex_double* a, lapack_int* lda,
					lapack_complex_double* b, lapack_int* ldb,
					lapack_complex_double* work, lapack_int *info );
void LAPACK_dtpqrt( lapack_int* m, lapack_int* n, lapack_int* l, lapack_int* nb,
					double* a, lapack_int* lda, double* b, lapack_int* ldb,
					double* t, lapack_int* ldt, double* work,
					lapack_int *info );
void LAPACK_ctpqrt( lapack_int* m, lapack_int* n, lapack_int* l, lapack_int* nb,
					lapack_complex_float* a, lapack_int* lda,
					lapack_complex_float* b, lapack_int* ldb,
					lapack_complex_float* t, lapack_int* ldt,
					lapack_complex_float* work, lapack_int *info );
void LAPACK_ztpqrt( lapack_int* m, lapack_int* n, lapack_int* l, lapack_int* nb,
					lapack_complex_double* a, lapack_int* lda,
					lapack_complex_double* b, lapack_int* ldb,
					lapack_complex_double* t, lapack_int* ldt,
					lapack_complex_double* work, lapack_int *info );
void LAPACK_stpqrt2( lapack_int* m, lapack_int* n, lapack_int* l,
					float* a, lapack_int* lda,
					float* b, lapack_int* ldb,
					float* t, lapack_int* ldt,
					lapack_int *info );
void LAPACK_dtpqrt2( lapack_int* m, lapack_int* n, lapack_int* l,
					double* a, lapack_int* lda,
					double* b, lapack_int* ldb,
					double* t, lapack_int* ldt,
					lapack_int *info );
void LAPACK_ctpqrt2( lapack_int* m, lapack_int* n, lapack_int* l,
					lapack_complex_float* a, lapack_int* lda,
					lapack_complex_float* b, lapack_int* ldb,
					lapack_complex_float* t, lapack_int* ldt,
					lapack_int *info );
void LAPACK_ztpqrt2( lapack_int* m, lapack_int* n, lapack_int* l,
					lapack_complex_double* a, lapack_int* lda,
					lapack_complex_double* b, lapack_int* ldb,
					lapack_complex_double* t, lapack_int* ldt,
					lapack_int *info );
void LAPACK_stprfb( char* side, char* trans, char* direct, char* storev,
					lapack_int* m, lapack_int* n, lapack_int* k, lapack_int* l,
					const float* v, lapack_int* ldv, const float* t,
					lapack_int* ldt, float* a, lapack_int* lda, float* b,
					lapack_int* ldb, const float* work,
					lapack_int* ldwork );
void LAPACK_dtprfb( char* side, char* trans, char* direct, char* storev,
					lapack_int* m, lapack_int* n, lapack_int* k, lapack_int* l,
					const double* v, lapack_int* ldv, const double* t,
					lapack_int* ldt, double* a, lapack_int* lda, double* b,
					lapack_int* ldb, const double* work,
					lapack_int* ldwork );
void LAPACK_ctprfb( char* side, char* trans, char* direct, char* storev,
					lapack_int* m, lapack_int* n, lapack_int* k, lapack_int* l,
					const lapack_complex_float* v, lapack_int* ldv,
					const lapack_complex_float* t, lapack_int* ldt,
					lapack_complex_float* a, lapack_int* lda,
					lapack_complex_float* b, lapack_int* ldb,
					const float* work, lapack_int* ldwork );
void LAPACK_ztprfb( char* side, char* trans, char* direct, char* storev,
					lapack_int* m, lapack_int* n, lapack_int* k, lapack_int* l,
					const lapack_complex_double* v, lapack_int* ldv,
					const lapack_complex_double* t, lapack_int* ldt,
					lapack_complex_double* a, lapack_int* lda,
					lapack_complex_double* b, lapack_int* ldb,
					const double* work, lapack_int* ldwork );
// LAPACK 3.5.0
void LAPACK_ssysv_rook( char* uplo, lapack_int* n, lapack_int* nrhs, float* a,
						lapack_int* lda, lapack_int* ipiv, float* b,
						lapack_int* ldb, float* work, lapack_int* lwork,
						lapack_int *info );
void LAPACK_dsysv_rook( char* uplo, lapack_int* n, lapack_int* nrhs, double* a,
						lapack_int* lda, lapack_int* ipiv, double* b,
						lapack_int* ldb, double* work, lapack_int* lwork,
						lapack_int *info );
void LAPACK_csysv_rook( char* uplo, lapack_int* n, lapack_int* nrhs,
						lapack_complex_float* a, lapack_int* lda,
						lapack_int* ipiv, lapack_complex_float* b,
						lapack_int* ldb, lapack_complex_float* work,
						lapack_int* lwork, lapack_int *info );
void LAPACK_zsysv_rook( char* uplo, lapack_int* n, lapack_int* nrhs,
						lapack_complex_double* a, lapack_int* lda,
						lapack_int* ipiv, lapack_complex_double* b,
						lapack_int* ldb, lapack_complex_double* work,
						lapack_int* lwork, lapack_int *info );
void LAPACK_csyr( char* uplo, lapack_int* n, lapack_complex_float* alpha,
					const lapack_complex_float* x, lapack_int* incx,
					lapack_complex_float* a, lapack_int* lda );
void LAPACK_zsyr( char* uplo, lapack_int* n, lapack_complex_double* alpha,
					const lapack_complex_double* x, lapack_int* incx,
					lapack_complex_double* a, lapack_int* lda );
void LAPACK_ilaver( const lapack_int* vers_major, const lapack_int* vers_minor,
					const lapack_int* vers_patch );

}
