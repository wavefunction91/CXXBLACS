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
    FUNC(&TRANSA,&TRANSB,&M,&N,&K,&ALPHA,A,&IA,&JA,DESCA,B,&IB,&JB,\
      DESCB,&BETA,C,&IC,&JC,DESCC);\
  }

  PGEMM_IMPL(float               ,psgemm_);
  PGEMM_IMPL(double              ,pdgemm_);
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

};



#endif
