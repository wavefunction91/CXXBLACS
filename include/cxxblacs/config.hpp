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


template <typename T>
inline T* cc(const T *x){ return const_cast<T*>(x); };

// BLACS types

#ifndef CXXBLACS_BLACS_Complex16
  #define CXXBLACS_BLACS_Complex16 std::complex<double>
#endif

#ifndef CXXBLACS_BLACS_Complex8
  #define CXXBLACS_BLACS_Complex8 std::complex<float>
#endif

// PBLAS types

#ifndef CXXBLACS_PBLAS_Complex16
  #define CXXBLACS_PBLAS_Complex16 std::complex<double>
#endif

#ifndef CXXBLACS_PBLAS_Complex8
  #define CXXBLACS_PBLAS_Complex8 std::complex<float>
#endif


// SCALAPCK types

#ifndef CXXBLACS_SCALAPACK_Complex16
  #define CXXBLACS_SCALAPACK_Complex16 std::complex<double>
#endif

#ifndef CXXBLACS_SCALAPACK_Complex8
  #define CXXBLACS_SCALAPACK_Complex8 std::complex<float>
#endif


// BLAS types


#ifndef CXXBLACS_BLAS_Complex16
  #define CXXBLACS_BLAS_Complex16 std::complex<double>
#endif

#ifndef CXXBLACS_BLAS_Complex8
  #define CXXBLACS_BLAS_Complex8 std::complex<float>
#endif


// LAPACK types

#ifndef CXXBLACS_LAPACK_Complex16
  #define CXXBLACS_LAPACK_Complex16 std::complex<double>
#endif

#ifndef CXXBLACS_LAPACK_Complex8
  #define CXXBLACS_LAPACK_Complex8 std::complex<float>
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




  // Type conversions for PBLAS

  template <typename T>
  struct CXXBLACS_PBLAS_TYPE { typedef T type; };

  template <>
  struct CXXBLACS_PBLAS_TYPE<std::complex<double>> {
    typedef CXXBLACS_PBLAS_Complex16 type;
  };

  template <>
  struct CXXBLACS_PBLAS_TYPE<std::complex<float>> {
    typedef CXXBLACS_PBLAS_Complex8 type;
  };

  template <typename T, 
            typename PBLAS_TYPE = typename CXXBLACS_PBLAS_TYPE<T>::type >
  PBLAS_TYPE* ToPblasType(T* ptr){ return reinterpret_cast<PBLAS_TYPE*>(ptr); }



  // Type conversions for SCALAPACK

  template <typename T>
  struct CXXBLACS_SCALAPACK_TYPE { typedef T type; };

  template <>
  struct CXXBLACS_SCALAPACK_TYPE<std::complex<double>> {
    typedef CXXBLACS_SCALAPACK_Complex16 type;
  };

  template <>
  struct CXXBLACS_SCALAPACK_TYPE<std::complex<float>> {
    typedef CXXBLACS_SCALAPACK_Complex8 type;
  };

  template <typename T, 
            typename SCALAPACK_TYPE = 
              typename CXXBLACS_SCALAPACK_TYPE<T>::type >
  SCALAPACK_TYPE* ToScalapackType(T* ptr){ 
    return reinterpret_cast<SCALAPACK_TYPE*>(ptr); }



  // Type conversions for BLAS

  template <typename T>
  struct CXXBLACS_BLAS_TYPE { typedef T type; };

  template <>
  struct CXXBLACS_BLAS_TYPE<std::complex<double>> {
    typedef CXXBLACS_BLAS_Complex16 type;
  };

  template <>
  struct CXXBLACS_BLAS_TYPE<std::complex<float>> {
    typedef CXXBLACS_BLAS_Complex8 type;
  };

  template <typename T, 
            typename BLAS_TYPE = typename CXXBLACS_BLAS_TYPE<T>::type >
  BLAS_TYPE* ToBlasType(T* ptr){ return reinterpret_cast<BLAS_TYPE*>(ptr); }



  // Type conversions for LAPACK

  template <typename T>
  struct CXXBLACS_LAPACK_TYPE { typedef T type; };

  template <>
  struct CXXBLACS_LAPACK_TYPE<std::complex<double>> {
    typedef CXXBLACS_LAPACK_Complex16 type;
  };

  template <>
  struct CXXBLACS_LAPACK_TYPE<std::complex<float>> {
    typedef CXXBLACS_LAPACK_Complex8 type;
  };

  template <typename T, 
            typename LAPACK_TYPE = 
              typename CXXBLACS_LAPACK_TYPE<T>::type >
  LAPACK_TYPE* ToLapackType(T* ptr){ 
    return reinterpret_cast<LAPACK_TYPE*>(ptr); }



};

#endif
