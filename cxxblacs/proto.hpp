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


#ifndef CXXBLACS_HAS_BLACS
#define CXXBLACS_HAS_BLACS


extern "C" {

  // BLACS grid info
  void blacs_pinfo_(CB_INT*,CB_INT*);
  void blacs_get_(const CB_INT*,const CB_INT*,CB_INT*);
  void blacs_gridinit_(CB_INT*,const char*,const CB_INT*,const CB_INT*);
  void blacs_gridinfo_(const CB_INT*,CB_INT*,CB_INT*,CB_INT*,CB_INT*);
  void blacs_barrier_(const CB_INT*,const char*);
  void blacs_gridexit_(const CB_INT*);
  void blacs_exit_(const CB_INT*);
  CB_INT  blacs_pnum_(CB_INT*,CB_INT*,CB_INT*);
  CB_INT  Csys2blacs_handle( MPI_Comm );
  void Cfree_blacs_system_handle( CB_INT );


  // BLACS Point-to-point communication
  #define gesd_rv2d(F,FUNC) \
  void FUNC(const CB_INT*, const CB_INT *, const CB_INT *, F *, \
    const CB_INT*, const CB_INT *, const CB_INT *);

  gesd_rv2d(CB_INT                  ,igesd2d_);
  gesd_rv2d(float                   ,sgesd2d_);
  gesd_rv2d(double                  ,dgesd2d_);
  gesd_rv2d(CXXBLACS_BLACS_Complex8 ,cgesd2d_);
  gesd_rv2d(CXXBLACS_BLACS_Complex16,zgesd2d_);

  gesd_rv2d(CB_INT                  ,igerv2d_);
  gesd_rv2d(float                   ,sgerv2d_);
  gesd_rv2d(double                  ,dgerv2d_);
  gesd_rv2d(CXXBLACS_BLACS_Complex8 ,cgerv2d_);
  gesd_rv2d(CXXBLACS_BLACS_Complex16,zgerv2d_);



  // BLACS collectives

  // GSUM
  #define gesum(F,FUNC) \
  void FUNC(const CB_INT*, const char*, const char*, const CB_INT*, \
    const CB_INT*, F*, const CB_INT*, const CB_INT*, const CB_INT*);

  gesum(CB_INT                  ,igsum2d_);
  gesum(float                   ,sgsum2d_);
  gesum(double                  ,dgsum2d_);
  gesum(CXXBLACS_BLACS_Complex8 ,cgsum2d_);
  gesum(CXXBLACS_BLACS_Complex16,zgsum2d_);


  // Broadcast send
  #define gebs2d(F,FUNC) \
  void FUNC(const CB_INT *, const char*, const char*, const CB_INT*,\
    const CB_INT*, const F *, const CB_INT*);

  gebs2d(float                   ,sgebs2d_);
  gebs2d(double                  ,dgebs2d_);
  gebs2d(CXXBLACS_BLACS_Complex8 ,cgebs2d_);
  gebs2d(CXXBLACS_BLACS_Complex16,zgebs2d_);


  // Broadcast recv
  #define gebr2d(F,FUNC) \
  void FUNC(const CB_INT *, const char*, const char*, const CB_INT*,\
    const CB_INT*, const F *, const CB_INT*, const CB_INT *, const CB_INT*);

  gebr2d(float                   ,sgebr2d_);
  gebr2d(double                  ,dgebr2d_);
  gebr2d(CXXBLACS_BLACS_Complex8 ,cgebr2d_);
  gebr2d(CXXBLACS_BLACS_Complex16,zgebr2d_);

};
    
#endif
    






#ifndef CXXBLACS_HAS_PBLAS
#define CXXBLACS_HAS_PBLAS

extern "C" {

  // PBLAS
    
  #define pgemm(F,FUNC)\
  void FUNC(const char*, const char*,\
    const CB_INT*, const CB_INT*,const CB_INT*, const F*,\
    const F*, const CB_INT*, const CB_INT*, const CB_INT*,\
    const F*, const CB_INT*, const CB_INT*, const CB_INT*,\
    const F*, F*, const CB_INT*, const CB_INT*, const CB_INT*);

  pgemm(float                   ,psgemm_);
  pgemm(double                  ,pdgemm_);
  pgemm(CXXBLACS_PBLAS_Complex8 ,pcgemm_);
  pgemm(CXXBLACS_PBLAS_Complex16,pzgemm_);



  #define ptrmm(F,FUNC)\
  void FUNC(const char*, const char*, const char*, const char*,\
    const CB_INT*, const CB_INT*, const F*, const F*, const CB_INT*,\
    const CB_INT*, const CB_INT*, F*, const CB_INT*, const CB_INT*,\
    const CB_INT*);

  ptrmm(float                   ,pstrmm_);
  ptrmm(double                  ,pdtrmm_);
  ptrmm(CXXBLACS_PBLAS_Complex8 ,pctrmm_);
  ptrmm(CXXBLACS_PBLAS_Complex16,pztrmm_);

}


#endif

    
#ifndef CXXBLACS_HAS_SCALAPACK
#define CXXBLACS_HAS_SCALAPACK

extern "C" {

  // ScaLAPACK
      
  void descinit_(CB_INT*, const CB_INT*, const CB_INT*, const CB_INT*,
    const CB_INT*, const CB_INT*, const CB_INT*, const CB_INT*, 
    const CB_INT*, const CB_INT*);

  #define pgemr2d(F,FUNC)\
  void FUNC(const CB_INT*, const CB_INT*,\
    const F*, const CB_INT*, const CB_INT*, const CB_INT*,\
    F*, const CB_INT*, const CB_INT*, const CB_INT*,\
    const CB_INT*);

  pgemr2d(float                       ,psgemr2d_);
  pgemr2d(double                      ,pdgemr2d_);
  pgemr2d(CXXBLACS_SCALAPACK_Complex8 ,pcgemr2d_);
  pgemr2d(CXXBLACS_SCALAPACK_Complex16,pzgemr2d_);

  #define plascl(F,RF,FUNC)\
  void FUNC(const char*, const RF*, const RF*, const CB_INT*, const CB_INT*,\
    F*, const CB_INT*, const CB_INT*, const CB_INT*, CB_INT*);

  plascl(float                       ,float ,pslascl_);
  plascl(double                      ,double,pdlascl_);
  plascl(CXXBLACS_SCALAPACK_Complex8 ,float ,pclascl_);
  plascl(CXXBLACS_SCALAPACK_Complex16,double,pzlascl_);

  #define psyev(F,FUNC)\
  void FUNC(const char*, const char*, const CB_INT*, F*, const CB_INT*, \
    const CB_INT*, const CB_INT*, F*, F*, const CB_INT*, const CB_INT*,\
    const CB_INT*, F*, const CB_INT*, CB_INT*);

  #define psyevd(F,FUNC)\
  void FUNC(const char*, const char*, const CB_INT*, F*, const CB_INT*, \
    const CB_INT*, const CB_INT*, F*, F*, const CB_INT*, const CB_INT*,\
    const CB_INT*, F*, const CB_INT*, CB_INT*, const CB_INT*, CB_INT*);

  #define pheev(F,RF,FUNC)\
  void FUNC(const char*, const char*, const CB_INT*, F*, const CB_INT*, \
    const CB_INT*, const CB_INT*, RF*, F*, const CB_INT*, const CB_INT*,\
    const CB_INT*, F*, const CB_INT*, RF*, const CB_INT*, CB_INT*);

  #define pheevd(F,RF,FUNC)\
  void FUNC(const char*, const char*, const CB_INT*, F*, const CB_INT*, \
    const CB_INT*, const CB_INT*, RF*, F*, const CB_INT*, const CB_INT*,\
    const CB_INT*, F*, const CB_INT*, RF*, const CB_INT*, \
    CB_INT*, const CB_INT*, CB_INT*);

  psyev(float ,pssyev_);
  psyev(double,pdsyev_);
  psyevd(float ,pssyevd_);
  psyevd(double,pdsyevd_);

  pheev(CXXBLACS_SCALAPACK_Complex8 ,float ,pcheev_);
  pheev(CXXBLACS_SCALAPACK_Complex16,double,pzheev_);
  pheevd(CXXBLACS_SCALAPACK_Complex8 ,float ,pcheevd_);
  pheevd(CXXBLACS_SCALAPACK_Complex16,double,pzheevd_);




  #define pgesv(F,FUNC)\
  void FUNC(const CB_INT*, const CB_INT*, F*, const CB_INT*, const CB_INT*, \
    const CB_INT*, CB_INT*, F*, const CB_INT*, const CB_INT*,\
    const CB_INT*, CB_INT*);

  pgesv(float                       ,psgesv_);
  pgesv(double                      ,pdgesv_);
  pgesv(CXXBLACS_SCALAPACK_Complex8 ,pcgesv_);
  pgesv(CXXBLACS_SCALAPACK_Complex16,pzgesv_);




  #define ppotrf(F,FUNC)\
  void FUNC(const char*, const CB_INT*, F*, const CB_INT*, const CB_INT*,\
    const CB_INT*, CB_INT*);

  ppotrf(float                       ,pspotrf_);
  ppotrf(double                      ,pdpotrf_);
  ppotrf(CXXBLACS_SCALAPACK_Complex8 ,pcpotrf_);
  ppotrf(CXXBLACS_SCALAPACK_Complex16,pzpotrf_);



}

#endif




#ifndef CXXBLACS_HAS_LAPACK
#define CXXBLACS_HAS_LAPACK

extern "C" {

  // LAPACK
    
  #define lacpy(F,FUNC)\
  void FUNC(const char*, const CB_INT*, const CB_INT*, F *, const CB_INT *,\
    F *, const CB_INT*);

  lacpy(float                    ,slacpy_);
  lacpy(double                   ,dlacpy_);
  lacpy(CXXBLACS_LAPACK_Complex8 ,clacpy_);
  lacpy(CXXBLACS_LAPACK_Complex16,zlacpy_);

}

#endif

#ifndef CXXBLACS_HAS_BLAS
#define CXXBLACS_HAS_BLAS

extern "C" {

  #define gemm(F,FUNC)\
  void FUNC(const char*, const char*, const CB_INT*, const CB_INT*,\
    const CB_INT*, const F*, const F*, const CB_INT*, const F*,\
    const CB_INT*,const F*, F*, const CB_INT*);

  gemm(float                  ,sgemm_);
  gemm(double                 ,dgemm_);
  gemm(CXXBLACS_BLAS_Complex8 ,cgemm_);
  gemm(CXXBLACS_BLAS_Complex16,zgemm_);

  #define trmm(F,FUNC)\
  void FUNC(const char*, const char*, const char*, const char*,\
    const CB_INT*, const CB_INT*, const F*, const F*, const CB_INT*,\
    F*, const CB_INT*);
      
  trmm(float                  ,strmm_);
  trmm(double                 ,dtrmm_);
  trmm(CXXBLACS_BLAS_Complex8 ,ctrmm_);
  trmm(CXXBLACS_BLAS_Complex16,ztrmm_);

}

#endif

#endif

