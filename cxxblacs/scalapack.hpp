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
#ifndef __INCLUDED_CXXBLACS_SCALAPACK_HPP__
#define __INCLUDED_CXXBLACS_SCALAPACK_HPP__

#include <cxxblacs/config.hpp>
#include <cxxblacs/proto.hpp>

namespace CXXBLACS {

  template <typename Field>
  inline void PGEMM(const char TRANSA, const char TRANSB, const CB_INT M,
    const CB_INT N, const CB_INT K, const Field ALPHA, const Field* A,
    const CB_INT IA, const CB_INT JA, const CB_INT* DESCA,
    const Field* B, const CB_INT IB, const CB_INT JB, 
    const CB_INT* DESCB, const Field BETA, Field* C,
    const CB_INT IC, const CB_INT JC, const CB_INT* DESCC);

  #define PGEMM_IMPL(F,FUNC)\
  template <>\
  inline void PGEMM(const char TRANSA, const char TRANSB, const CB_INT M,\
    const CB_INT N, const CB_INT K, const F ALPHA, const F* A,\
    const CB_INT IA, const CB_INT JA, const CB_INT* DESCA,\
    const F* B, const CB_INT IB, const CB_INT JB, \
    const CB_INT* DESCB, const F BETA, F* C,\
    const CB_INT IC, const CB_INT JC, const CB_INT* DESCC){\
    \
    FUNC(&TRANSA,&TRANSB,&M,&N,&K,ToPblasType(&ALPHA),ToPblasType(A),&IA,&JA,\
      DESCA,ToPblasType(B),&IB,&JB,DESCB,ToPblasType(&BETA),ToPblasType(C),\
      &IC,&JC,DESCC);\
  }

  PGEMM_IMPL(float                   ,psgemm_);
  PGEMM_IMPL(double                  ,pdgemm_);
  PGEMM_IMPL(std::complex<float> ,pcgemm_);
  PGEMM_IMPL(std::complex<double>,pzgemm_);

  template <typename Field>
  inline void PGEMM(const char TRANSA, const char TRANSB, const CB_INT M,
    const CB_INT N, const CB_INT K, const Field ALPHA, const Field* A,
    const CB_INT IA, const CB_INT JA, const ScaLAPACK_Desc_t DESCA,
    const Field* B, const CB_INT IB, const CB_INT JB, 
    const ScaLAPACK_Desc_t DESCB, const Field BETA, Field* C,
    const CB_INT IC, const CB_INT JC, const ScaLAPACK_Desc_t DESCC) {

    PGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,IA,JA,&DESCA[0],B,IB,JB,&DESCB[0],
      BETA,C,IC,JC,&DESCC[0]);

  }







  template <typename Field>
  inline void PTRMM(const char SIDE, const char UPLO, const char TRANSA,
    const char DIAG, const CB_INT M, const CB_INT N, const Field ALPHA,
    const Field *A, const CB_INT IA, const CB_INT JA, const CB_INT *DESCA,
    Field *B, const CB_INT IB, const CB_INT JB, const CB_INT *DESCB);


  #define PTRMM_IMPL(F,FUNC)\
  template <>\
  inline void PTRMM(const char SIDE, const char UPLO, const char TRANSA,\
    const char DIAG, const CB_INT M, const CB_INT N, const F ALPHA,\
    const F *A, const CB_INT IA, const CB_INT JA, const CB_INT *DESCA,\
    F *B, const CB_INT IB, const CB_INT JB, const CB_INT *DESCB) {\
    FUNC(&SIDE,&UPLO,&TRANSA,&DIAG,&M,&N,ToPblasType(&ALPHA),ToPblasType(A),\
      &IA,&JA,DESCA,ToPblasType(B),&IB,&JB,DESCB);\
  }

  PTRMM_IMPL(float                   ,pstrmm_);
  PTRMM_IMPL(double                  ,pdtrmm_);
  PTRMM_IMPL(std::complex<float> ,pctrmm_);
  PTRMM_IMPL(std::complex<double>,pztrmm_);

  template <typename Field>
  inline void PTRMM(const char SIDE, const char UPLO, const char TRANSA,
    const char DIAG, const CB_INT M, const CB_INT N, const Field ALPHA,
    const Field *A, const CB_INT IA, const CB_INT JA, 
    const ScaLAPACK_Desc_t DESCA, Field *B, const CB_INT IB, const CB_INT JB, 
    const ScaLAPACK_Desc_t DESCB) {

    PTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,IA,JA,&DESCA[0],B,IB,JB,&DESCB[0]);

  }












  template <typename Field>
  inline void PGEMR2D(const CB_INT M, const CB_INT N, const Field *A,
    const CB_INT IA, const CB_INT JA, const CB_INT *DESCA, Field *B,
    const CB_INT IB, const CB_INT JB, const CB_INT *DESCB, const CB_INT ICTXT);

  #define PGEMR2D_IMPL(F,FUNC)\
  template <>\
  inline void PGEMR2D(const CB_INT M, const CB_INT N, const F *A,\
    const CB_INT IA, const CB_INT JA, const CB_INT *DESCA, F *B,\
    const CB_INT IB, const CB_INT JB, const CB_INT *DESCB, const CB_INT ICTXT){\
    \
    FUNC(&M,&N,ToScalapackType(A),&IA,&JA,DESCA,ToScalapackType(B),&IB,&JB,\
      DESCB,&ICTXT);\
  }

  PGEMR2D_IMPL(float                       ,psgemr2d_);
  PGEMR2D_IMPL(double                      ,pdgemr2d_);
  PGEMR2D_IMPL(std::complex<float> ,pcgemr2d_);
  PGEMR2D_IMPL(std::complex<double>,pzgemr2d_);

  template <typename Field>
  inline void PGEMR2D(const CB_INT M, const CB_INT N, const Field *A,
    const CB_INT IA, const CB_INT JA, const ScaLAPACK_Desc_t DESCA, Field *B,
    const CB_INT IB, const CB_INT JB, const ScaLAPACK_Desc_t DESCB, 
    const CB_INT ICTXT) {

    PGEMR2D(M,N,A,IA,JA,&DESCA[0],B,IB,JB,&DESCB[0],ICTXT);

  }




  template <typename Field>
  inline CB_INT PSYEV(const char JOBZ, const char UPLO, const CB_INT N,
    Field *A, const CB_INT IA, const CB_INT JA, const CB_INT *DESCA,
    Field *W, Field *Z, const CB_INT IZ, const CB_INT JZ, const CB_INT *DESCZ,
    Field *WORK, const CB_INT LWORK);

  template <typename Field, typename RealField>
  inline CB_INT PHEEV(const char JOBZ, const char UPLO, const CB_INT N,
    Field *A, const CB_INT IA, const CB_INT JA, const CB_INT *DESCA,
    RealField *W, Field *Z, const CB_INT IZ, const CB_INT JZ, 
    const CB_INT *DESCZ, Field *WORK, const CB_INT LWORK, RealField *RWORK,
    const CB_INT LRWORK);

  #define PSYEV_IMPL(F,FUNC)\
  template <>\
  inline CB_INT PSYEV(const char JOBZ, const char UPLO, const CB_INT N,\
    F *A, const CB_INT IA, const CB_INT JA, const CB_INT *DESCA,\
    F *W, F *Z, const CB_INT IZ, const CB_INT JZ, const CB_INT *DESCZ,\
    F *WORK, const CB_INT LWORK) {\
    \
    if( DESCA[4] != DESCA[5] ) {\
      std::runtime_error err("MB must be the same as NB in P?SYEV");\
      throw err;\
    }\
    CB_INT INFO;\
    FUNC(&JOBZ,&UPLO,&N,A,&IA,&JA,DESCA,W,Z,&IZ,&JZ,DESCZ,WORK,&LWORK,&INFO);\
    return INFO;\
    \
  }

  #define PHEEV_IMPL(F,RF,FUNC)\
  template <>\
  inline CB_INT PHEEV(const char JOBZ, const char UPLO, const CB_INT N,\
    F *A, const CB_INT IA, const CB_INT JA, const CB_INT *DESCA,\
    RF *W, F *Z, const CB_INT IZ, const CB_INT JZ, const CB_INT *DESCZ,\
    F *WORK, const CB_INT LWORK, RF* RWORK, const CB_INT LRWORK) {\
    \
    if( DESCA[4] != DESCA[5] ) {\
      std::runtime_error err("MB must be the same as NB in P?HEEV");\
      throw err;\
    }\
    CB_INT INFO;\
    FUNC(&JOBZ,&UPLO,&N,ToScalapackType(A),&IA,&JA,DESCA,W,ToScalapackType(Z),\
      &IZ,&JZ,DESCZ,ToScalapackType(WORK),&LWORK,RWORK,&LRWORK,&INFO);\
    return INFO;\
    \
  }

  PSYEV_IMPL(float ,pssyev_);
  PSYEV_IMPL(double,pdsyev_);

  PHEEV_IMPL(std::complex<float> ,float ,pcheev_);
  PHEEV_IMPL(std::complex<double>,double,pzheev_);

  template <typename Field>
  inline CB_INT PSYEV(const char JOBZ, const char UPLO, const CB_INT N,
    Field *A, const CB_INT IA, const CB_INT JA, const ScaLAPACK_Desc_t DESCA,
    Field *W, Field *Z, const CB_INT IZ, const CB_INT JZ, 
    const ScaLAPACK_Desc_t DESCZ, Field *WORK, const CB_INT LWORK) {

    return PSYEV(JOBZ,UPLO,N,A,IA,JA,&DESCA[0],W,Z,IZ,JZ,&DESCZ[0],WORK,LWORK);

  }

  template <typename Field, typename RealField>
  inline CB_INT PHEEV(const char JOBZ, const char UPLO, const CB_INT N,
    Field *A, const CB_INT IA, const CB_INT JA, const ScaLAPACK_Desc_t DESCA,
    RealField *W, Field *Z, const CB_INT IZ, const CB_INT JZ, 
    const ScaLAPACK_Desc_t DESCZ, Field *WORK, const CB_INT LWORK, 
    RealField *RWORK, const CB_INT LRWORK) {

    return PHEEV(JOBZ,UPLO,N,A,IA,JA,&DESCA[0],W,Z,IZ,JZ,&DESCZ[0],WORK,LWORK,
        RWORK,LRWORK);

  }





  template <typename Field>
  inline CB_INT PGESV(const CB_INT N, const CB_INT NRHS, Field *A, 
    const CB_INT IA, const CB_INT JA, const CB_INT *DESCA, 
    CB_INT *IPIV, Field *B, const CB_INT IB, const CB_INT JB,
    const CB_INT *DESCB);

  #define PGESV_IMPL(F,FUNC)\
  template <>\
  inline CB_INT PGESV(const CB_INT N, const CB_INT NRHS, F *A, \
    const CB_INT IA, const CB_INT JA, const CB_INT *DESCA, \
    CB_INT *IPIV, F *B, const CB_INT IB, const CB_INT JB, \
    const CB_INT *DESCB) {\
    \
    CB_INT INFO;\
    FUNC(&N,&NRHS,ToScalapackType(A),&IA,&JA,DESCA,IPIV,ToScalapackType(B),\
      &IB,&JB,DESCB,&INFO);\
    return INFO;\
    \
  }

  PGESV_IMPL(float                       ,psgesv_);
  PGESV_IMPL(double                      ,pdgesv_);
  PGESV_IMPL(std::complex<float> ,pcgesv_);
  PGESV_IMPL(std::complex<double>,pzgesv_);

  template <typename Field>
  inline CB_INT PGESV(const CB_INT N, const CB_INT NRHS, Field *A, 
    const CB_INT IA, const CB_INT JA, const ScaLAPACK_Desc_t DESCA, 
    CB_INT *IPIV, Field *B, const CB_INT IB, const CB_INT JB,
    const ScaLAPACK_Desc_t DESCB) {

    return PGESV(N,NRHS,A,IA,JA,&DESCA[0],IPIV,B,IB,JB,&DESCB[0]);

  }


  template <typename Field>
  inline CB_INT PPOTRF(const char UPLO, const CB_INT N, Field *A,
    const CB_INT IA, const CB_INT JA, const CB_INT *DESCA);

  #define PPOTRF_IMPL(F,FUNC)\
  template <>\
  inline CB_INT PPOTRF(const char UPLO, const CB_INT N, F *A,\
    const CB_INT IA, const CB_INT JA, const CB_INT *DESCA) {\
    \
    CB_INT INFO;\
    FUNC(&UPLO,&N,ToScalapackType(A),&IA,&JA,DESCA,&INFO);\
    return INFO;\
    \
  }

  PPOTRF_IMPL(float                       ,pspotrf_);
  PPOTRF_IMPL(double                      ,pdpotrf_);
  PPOTRF_IMPL(std::complex<float> ,pcpotrf_);
  PPOTRF_IMPL(std::complex<double>,pzpotrf_);


  template <typename Field>
  inline CB_INT PPOTRF(const char UPLO, const CB_INT N, Field *A,
    const CB_INT IA, const CB_INT JA, const ScaLAPACK_Desc_t DESCA) {

    return PPOTRF(UPLO,N,A,IA,JA,&DESCA[0]);

  }



};



#endif
