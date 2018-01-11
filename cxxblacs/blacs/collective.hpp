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
  inline void GSUM2D(const CB_INT ICONTXT, const char SCOPE[], 
    const char TOP[], const CB_INT M, const CB_INT N, Field *A, 
    const CB_INT LDA, const CB_INT RDest, const CB_INT CDest);
  
  #define GSUM2D_IMPL(FIELD,FUNC)\
    template<>\
    inline void GSUM2D(const CB_INT ICONTXT, const char SCOPE[], \
      const char TOP[], const CB_INT M, const CB_INT N, FIELD *A,\
      const CB_INT LDA, const CB_INT RDest, const CB_INT CDest) {\
        FUNC(&ICONTXT,SCOPE,TOP,&M,&N,ToBlacsType(A),&LDA,&RDest,&CDest);\
    }

  GSUM2D_IMPL(CB_INT              ,igsum2d_);
  GSUM2D_IMPL(float               ,sgsum2d_);
  GSUM2D_IMPL(double              ,dgsum2d_);
  GSUM2D_IMPL(std::complex<float> ,cgsum2d_);
  GSUM2D_IMPL(std::complex<double>,zgsum2d_);


  // Broadcast Send
  template<typename Field>
  inline void GEBS2D(const CB_INT ICONTXT, const char SCOPE[], 
    const char TOP[], const CB_INT M, const CB_INT N, Field *A, 
    const CB_INT LDA);

  #define GEBS2D_IMPL(FIELD,FUNC)\
    template<>\
    inline void GEBS2D(const CB_INT ICONTXT, const char SCOPE[], \
      const char TOP[], const CB_INT M, const CB_INT N, FIELD *A,\
      const CB_INT LDA) {\
        FUNC(&ICONTXT,SCOPE,TOP,&M,&N,ToBlacsType(A),&LDA);\
    }

  GEBS2D_IMPL(float               ,sgebs2d_);
  GEBS2D_IMPL(double              ,dgebs2d_);
  GEBS2D_IMPL(std::complex<float> ,cgebs2d_);
  GEBS2D_IMPL(std::complex<double>,zgebs2d_);


  // Broadcast Send
  template<typename Field>
  inline void GEBR2D(const CB_INT ICONTXT, const char SCOPE[], 
    const char TOP[], const CB_INT M, const CB_INT N, Field *A, 
    const CB_INT LDA, const CB_INT RSrc, const CB_INT CSrc);
  
  #define GEBR2D_IMPL(FIELD,FUNC)\
    template<>\
    inline void GEBR2D(const CB_INT ICONTXT, const char SCOPE[], \
      const char TOP[], const CB_INT M, const CB_INT N, FIELD *A,\
      const CB_INT LDA, const CB_INT RSrc, const CB_INT CSrc) {\
        FUNC(&ICONTXT,SCOPE,TOP,&M,&N,ToBlacsType(A),&LDA,&RSrc,&CSrc);\
    }

  GEBR2D_IMPL(float               ,sgebr2d_);
  GEBR2D_IMPL(double              ,dgebr2d_);
  GEBR2D_IMPL(std::complex<float> ,cgebr2d_);
  GEBR2D_IMPL(std::complex<double>,zgebr2d_);

};

#endif
