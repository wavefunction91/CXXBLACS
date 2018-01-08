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
#ifndef __INCLUDED_CXXBLACS_BLACS_COLLECTIVE_HPP__
#define __INCLUDED_CXXBLACS_BLACS_COLLECTIVE_HPP__

#include <cxxblacs/proto.hpp>

namespace CXXBLACS {


  /**
   * \brief C++ Wrapper for ?GERV2D
   *
   * A templated function that encompasses the functionaliy of all of the BLACS
   * collective matrix sum routines: [I/S/D/C/Z]GSUM2D. Template deduction is based
   * on inut parameters
   *
   * See BLACS Documentaion for specifics.
   */
  template<typename Field>
  inline void GSUM2D(const int ICONTXT, const char SCOPE[], 
    const char TOP[], const int M, const int N, Field *A, const int LDA,
    const int RDest, const int CDest);
  
  // Specializations

  #define GSUM2D_IMPL(FIELD,FUNC)\
    template<>\
    inline void GSUM2D(const int ICONTXT, const char SCOPE[], \
      const char TOP[], const int M, const int N, FIELD *A, const int LDA,\
      const int RDest, const int CDest) {\
        FUNC(&ICONTXT,SCOPE,TOP,&M,&N,A,&LDA,&RDest,&CDest);\
    }

  GSUM2D_IMPL(int                 ,igsum2d_);
  GSUM2D_IMPL(float               ,sgsum2d_);
  GSUM2D_IMPL(double              ,dgsum2d_);
  GSUM2D_IMPL(std::complex<float> ,cgsum2d_);
  GSUM2D_IMPL(std::complex<double>,zgsum2d_);

};

#endif
