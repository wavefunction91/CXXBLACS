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
#ifndef __INCLUDED_CXXBLACS_PROTO_HPP__
#define __INCLUDED_CXXBLACS_PROTO_HPP__

#include <cxxblacs/config.hpp>


// External function prototypes

extern "C" {

  // BLACS grid info
  void blacs_pinfo_(CB_INT*,CB_INT*);
  void blacs_get_(const CB_INT*,const CB_INT*,CB_INT*);
  void blacs_gridinit_(CB_INT*,const char*,const CB_INT*,const CB_INT*);
  void blacs_gridinfo_(const CB_INT*,const CB_INT*,const CB_INT*,CB_INT*,CB_INT*);
  void blacs_barrier_(const CB_INT*,const char*);
  void blacs_gridexit_(const CB_INT*);
  void blacs_exit_(const CB_INT*);
  CB_INT  blacs_pnum_(CB_INT*,CB_INT*,CB_INT*);


  // BLACS Point-to-point communication
  #define gesd_rv2d(F,FUNC) \
  void FUNC(const CB_INT*, const CB_INT *, const CB_INT *, const F *, \
    const CB_INT*, const CB_INT *, const CB_INT *);

  gesd_rv2d(CB_INT              ,igesd2d_);
  gesd_rv2d(float               ,sgesd2d_);
  gesd_rv2d(double              ,dgesd2d_);
  gesd_rv2d(std::complex<float> ,cgesd2d_);
  gesd_rv2d(std::complex<double>,zgesd2d_);

  gesd_rv2d(CB_INT              ,igerv2d_);
  gesd_rv2d(float               ,sgerv2d_);
  gesd_rv2d(double              ,dgerv2d_);
  gesd_rv2d(std::complex<float> ,cgerv2d_);
  gesd_rv2d(std::complex<double>,zgerv2d_);



  // BLACS collectives

  // GSUM
  #define gesum(F,FUNC) \
  void FUNC(const CB_INT*, const char*, const char*, const CB_INT*, \
    const CB_INT*, F*, const CB_INT*, const CB_INT*, const CB_INT*);

  gesum(CB_INT              ,igsum2d_);
  gesum(float               ,sgsum2d_);
  gesum(double              ,dgsum2d_);
  gesum(std::complex<float> ,cgsum2d_);
  gesum(std::complex<double>,zgsum2d_);


  // Broadcast send
  #define gebs2d(F,FUNC) \
  void FUNC(const CB_INT *, const char*, const char*, const CB_INT*,\
    const CB_INT*, const F *, const CB_INT*);

  gebs2d(float               ,sgebs2d_);
  gebs2d(double              ,dgebs2d_);
  gebs2d(std::complex<float> ,cgebs2d_);
  gebs2d(std::complex<double>,zgebs2d_);


  // Broadcast recv
  #define gebr2d(F,FUNC) \
  void FUNC(const CB_INT *, const char*, const char*, const CB_INT*,\
    const CB_INT*, const F *, const CB_INT*, const CB_INT *, const CB_INT*);

  gebr2d(float               ,sgebr2d_);
  gebr2d(double              ,dgebr2d_);
  gebr2d(std::complex<float> ,cgebr2d_);
  gebr2d(std::complex<double>,zgebr2d_);
    
    
    
    

  // ScaLAPACK + PBLAS
      
  void descinit_(CB_INT*, const CB_INT*, const CB_INT*, const CB_INT*,
    const CB_INT*, const CB_INT*, const CB_INT*, const CB_INT*, 
    const CB_INT*, const CB_INT*);

  #define pgemm(F,FUNC)\
  void FUNC(const char*, const char*,\
    const CB_INT*, const CB_INT*,const CB_INT*, const F*,\
    const F*, const CB_INT*, const CB_INT*, const CB_INT*,\
    const F*, const CB_INT*, const CB_INT*, const CB_INT*,\
    const F*, F*, const CB_INT*, const CB_INT*, const CB_INT*);

  pgemm(float               ,psgemm_);
  pgemm(double              ,pdgemm_);
  pgemm(std::complex<float> ,pcgemm_);
  pgemm(std::complex<double>,pzgemm_);



  #define pgemr2d(F,FUNC)\
  void FUNC(const CB_INT*, const CB_INT*,\
    const F*, const CB_INT*, const CB_INT*, const CB_INT*,\
    F*, const CB_INT*, const CB_INT*, const CB_INT*,\
    const CB_INT*);

  pgemr2d(float               ,psgemr2d_);
  pgemr2d(double              ,pdgemr2d_);
  pgemr2d(std::complex<float> ,pcgemr2d_);
  pgemr2d(std::complex<double>,pzgemr2d_);


  // LAPACK + BLAS
    
  #define lacpy(F,FUNC)\
  void FUNC(const char*, const CB_INT*, const CB_INT*, F *, const CB_INT *,\
    F *, const CB_INT*);

  lacpy(float               ,slacpy_);
  lacpy(double              ,dlacpy_);
  lacpy(std::complex<float> ,clacpy_);
  lacpy(std::complex<double>,zlacpy_);


  #define gemm(F,FUNC)\
  void FUNC(const char*, const char*, const CB_INT*, const CB_INT*,\
    const CB_INT*, const F*, const F*, const CB_INT*, const F*,\
    const CB_INT*,const F*, F*, const CB_INT*);

  gemm(float               ,sgemm_);
  gemm(double              ,dgemm_);
  gemm(std::complex<float> ,cgemm_);
  gemm(std::complex<double>,zgemm_);

};

#endif

