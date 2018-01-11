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
#ifndef __INCLUDED_CXXBLACS_CONFIG_HPP__
#define __INCLUDED_CXXBLACS_CONFIG_HPP__

#include <mpi.h>

#include <iostream>
#include <iomanip>
#include <sstream>
#include <cassert>
#include <complex>
#include <tuple>

#ifdef MKL_ILP64
  #define CB_INT int64_t
#elif !defined(CB_INT)
  #define CB_INT int32_t
#endif


// BLACS types

#ifndef CXXBLACS_BLACS_Complex16
  #define CXXBLACS_BLACS_Complex16 std::complex<double>
#endif

#ifndef CXXBLACS_BLACS_Complex8
  #define CXXBLACS_BLACS_Complex8 std::complex<float>
#endif

namespace CXXBLACS {

  typedef std::pair<CB_INT,CB_INT> INDX;

  static constexpr CB_INT DESCINIT_LEN_MAX = 9;
  typedef std::array<CB_INT,DESCINIT_LEN_MAX> ScaLAPACK_Desc_t;


  // Type conversions for BLACS

  template <typename T>
  struct CXXBLACS_BLACS_TYPE { typedef T type; };

  template <>
  struct CXXBLACS_BLACS_TYPE<std::complex<double>> {
    typedef CXXBLACS_BLACS_Complex16 type;
  };

  template <>
  struct CXXBLACS_BLACS_TYPE<std::complex<float>> {
    typedef CXXBLACS_BLACS_Complex8 type;
  };

  template <typename T, 
            typename BLACS_TYPE = typename CXXBLACS_BLACS_TYPE<T>::type >
  BLACS_TYPE* ToBlacsType(T* ptr){ return reinterpret_cast<BLACS_TYPE*>(ptr); }

};

#endif
