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
#ifndef __INCLUDED_CXXBLACS_BLACS_POINTTOPOINT_HPP__
#define __INCLUDED_CXXBLACS_BLACS_POINTTOPOINT_HPP__

#include <cxxblacs/proto.hpp>

namespace CXXBLACS {



  /**
   * \brief C++ Wrapper for ?GESD2D
   *
   * A templated function that encompasses the functionaliy of all of the BLACS
   * point-to-point send routines: [I/S/D/C/Z]GESD2D. Template deduction is based
   * on inut parameters
   *
   * See BLACS Documentaion for specifics.
   */
  template<typename Field>
  inline void GESD2D(const int ICONTXT, const int M, const int N,
    const Field *A, const int LDA, const int RDest, const int CDest);




  #define GESD2D_IMPL(FIELD,FUNC)\
  template<>\
  inline void GESD2D(const int ICONTXT, const int M, const int N,\
    const FIELD *A, const int LDA, const int RDest, const int CDest){\
      FUNC(&ICONTXT,&M,&N,A,&LDA,&RDest,&CDest);\
  }


  GESD2D_IMPL(int                 ,igesd2d_);
  GESD2D_IMPL(float               ,sgesd2d_);
  GESD2D_IMPL(double              ,dgesd2d_);
  GESD2D_IMPL(std::complex<float> ,cgesd2d_);
  GESD2D_IMPL(std::complex<double>,zgesd2d_);






  /**
   * \brief C++ Wrapper for ?GERV2D
   *
   * A templated function that encompasses the functionaliy of all of the BLACS
   * point-to-point recieve routines: [I/S/D/C/Z]GERV2D. Template deduction is based
   * on inut parameters
   *
   * See BLACS Documentaion for specifics.
   */
  template<typename Field>
  inline void GERV2D(const int ICONTXT, const int M, const int N,
    Field *A, const int LDA, const int RDest, const int CDest);



  #define GERV2D_IMPL(FIELD,FUNC)\
  template<>\
  inline void GERV2D(const int ICONTXT, const int M, const int N,\
    FIELD *A, const int LDA, const int RDest, const int CDest){\
      FUNC(&ICONTXT,&M,&N,A,&LDA,&RDest,&CDest);\
  }

  

  GERV2D_IMPL(int                 ,igerv2d_);
  GERV2D_IMPL(float               ,sgerv2d_);
  GERV2D_IMPL(double              ,dgerv2d_);
  GERV2D_IMPL(std::complex<float> ,cgerv2d_);
  GERV2D_IMPL(std::complex<double>,zgerv2d_);
  
  

}; // namespace CXXBLACS

#endif
