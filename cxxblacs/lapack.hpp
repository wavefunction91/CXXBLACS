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

  #define LACOPY_IMPL(F,FUNC)\
  template <>\
  inline void LACOPY(const char UPLO, const CB_INT M, const CB_INT N,\
    F *A, const CB_INT LDA, F *B, const CB_INT LDB) {\
    FUNC(&UPLO,&M,&N,A,&LDA,B,&LDB);\
  }

  LACOPY_IMPL(float               ,slacpy_);
  LACOPY_IMPL(double              ,dlacpy_);
  LACOPY_IMPL(std::complex<float> ,clacpy_);
  LACOPY_IMPL(std::complex<double>,zlacpy_);


};

#endif
