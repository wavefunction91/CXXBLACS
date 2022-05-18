/*
 *  A simple C++ Wrapper for BLACS along with minimal extra functionality to 
 *  aid the the high-level development of distributed memory linear algebra.
 *  Copyright (C) 2016-2018 David Williams-Young

 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.

 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.

 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef __INCLUDED_CXXBLACS_LAPACK_HPP__
#define __INCLUDED_CXXBLACS_LAPACK_HPP__

#include <cxxblacs/config.hpp>
#include <cxxblacs/proto.hpp>

namespace CXXBLACS {

  template <typename Field>
  inline void LACOPY(const char UPLO, const CB_INT M, const CB_INT N,
    Field *A, const CB_INT LDA, Field *B, const CB_INT LDB);

#ifdef LAPACK_FORTRAN_STRLEN_END
  #define LACOPY_IMPL(F,FUNC)\
  template <>\
  inline void LACOPY(const char UPLO, const CB_INT M, const CB_INT N,\
    F *A, const CB_INT LDA, F *B, const CB_INT LDB) {\
    FUNC(cc(&UPLO),cc(&M),cc(&N),ToLapackType(cc(A)),cc(&LDA),\
      ToLapackType(cc(B)),cc(&LDB),1);\
  }
#else
#define LACOPY_IMPL(F,FUNC)\
  template <>\
  inline void LACOPY(const char UPLO, const CB_INT M, const CB_INT N,\
    F *A, const CB_INT LDA, F *B, const CB_INT LDB) {\
    FUNC(cc(&UPLO),cc(&M),cc(&N),ToLapackType(cc(A)),cc(&LDA),\
      ToLapackType(cc(B)),cc(&LDB));\
  }
#endif
  LACOPY_IMPL(float                    ,slacpy_);
  LACOPY_IMPL(double                   ,dlacpy_);
  LACOPY_IMPL(std::complex<float> ,clacpy_);
  LACOPY_IMPL(std::complex<double>,zlacpy_);




  template <typename Field>
  inline void GEMM(const char TRANSA, const char TRANSB, const CB_INT M,
    const CB_INT N, const CB_INT K, const Field ALPHA, const Field *A,
    const CB_INT LDA, const Field *B, const CB_INT LDB, const Field BETA,
    Field *C, const CB_INT LDC);

  #define GEMM_IMPL(F,FUNC)\
  template <>\
  inline void GEMM(const char TRANSA, const char TRANSB, const CB_INT M,\
    const CB_INT N, const CB_INT K, const F ALPHA, const F *A,\
    const CB_INT LDA, const F *B, const CB_INT LDB, const F BETA,\
    F *C, const CB_INT LDC) {\
    FUNC(cc(&TRANSA),cc(&TRANSB),cc(&M),cc(&N),cc(&K),ToBlasType(cc(&ALPHA)),\
        ToBlasType(cc(A)),cc(&LDA),ToBlasType(cc(B)),cc(&LDB),ToBlasType(cc(&BETA)),\
        ToBlasType(cc(C)),cc(&LDC));\
  };

  GEMM_IMPL(float                  ,sgemm_);
  GEMM_IMPL(double                 ,dgemm_);
  GEMM_IMPL(std::complex<float> ,cgemm_);
  GEMM_IMPL(std::complex<double>,zgemm_);





  template <typename Field>
  inline void TRMM(const char SIDE, const char UPLO, const char TRANSA, 
    const char DIAG, const CB_INT M, const CB_INT N, const Field ALPHA,
    const Field *A, const CB_INT LDA, Field *B, const CB_INT LDB);


  #define TRMM_IMPL(F,FUNC)\
  template <>\
  inline void TRMM(const char SIDE, const char UPLO, const char TRANSA, \
    const char DIAG, const CB_INT M, const CB_INT N, const F ALPHA,\
    const F *A, const CB_INT LDA, F *B, const CB_INT LDB){\
    FUNC(cc(&SIDE),cc(&UPLO),cc(&TRANSA),cc(&DIAG),cc(&M),cc(&N),\
      ToBlasType(cc(&ALPHA)),ToBlasType(cc(A)),cc(&LDA),ToBlasType(cc(B)),\
      cc(&LDB));\
  }

  TRMM_IMPL(float                  ,strmm_);
  TRMM_IMPL(double                 ,dtrmm_);
  TRMM_IMPL(std::complex<float> ,ctrmm_);
  TRMM_IMPL(std::complex<double>,ztrmm_);

};

#endif
