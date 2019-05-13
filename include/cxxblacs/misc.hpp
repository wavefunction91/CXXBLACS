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
#ifndef __INCLUDED_CXXBLACS_MISC_HPP__
#define __INCLUDED_CXXBLACS_MISC_HPP__

#include <cxxblacs/config.hpp>

namespace CXXBLACS {

  /**
   * \brief Implementation of ScaLAPACK's NUMROC
   *
   * See ScaLAPACK Documentation for details.
   *
   */
  inline CB_INT NumRoc(const CB_INT N, const CB_INT NB, const CB_INT IPROC,
    const CB_INT ISRCPROC, const CB_INT NPROCS) {

    CB_INT DIST    = (NPROCS + IPROC - ISRCPROC) % NPROCS;
    CB_INT NBLOCKS = N / NB;


    CB_INT DIM = ( NBLOCKS / NPROCS ) * NB;

    CB_INT EXTRA = NBLOCKS % NPROCS;

    if(      DIST <  EXTRA ) DIM += NB;     // Extra block
    else if( DIST == EXTRA ) DIM += N % NB; // Last block

    return DIM;

  };


  inline INDX GetLocalDims(const CB_INT M, const CB_INT N, const CB_INT MB,
    const CB_INT NB,   const CB_INT iProc, const CB_INT jProc,
    const CB_INT iSrc, const CB_INT jSrc,  const CB_INT nProcRow, 
    const CB_INT nProcCol) {

    return  { NumRoc(M,MB,iProc,iSrc,nProcRow),
              NumRoc(N,NB,jProc,jSrc,nProcCol) };
  }

  /**
   * \brief Calculate local 2D coordinates from global (distributed) 2D
   *        coordinates using general BLACS grid and coordinates
   *
   * @param[in]  ICONTXT BLACS context
   * @param[in]  MB Distributed row block size
   * @param[in]  NB Distributed column block size
   * @param[in]  NPROW Number of rows in process grid
   * @param[in]  NPCOL Number of columnss in process grid
   * @param[in]  I Distributed row coordinate
   * @param[in]  J Distributed column coordinate
   * @param[out] L Local row block of local buffer
   * @param[out] M Local column block of local buffer
   * @param[out] Pr Process row the contains distributed element (I,J)
   * @param[out] Pc Process column the contains distributed element (I,J)
   * @param[out] iX Local row coordinate
   * @param[out] iY Local column coordinate
   */
  inline void LocalFromGlobal(const CB_INT ICONTXT, const CB_INT MB, 
    const CB_INT NB, const CB_INT NPROW, const CB_INT NPCOL, const CB_INT I, const CB_INT J, 
    CB_INT &L, CB_INT &M, CB_INT &Pr, CB_INT &Pc, CB_INT &iX, CB_INT &iY) {
    
    L  = I / (NPROW * MB);
    M  = J / (NPCOL * NB);
    Pr = (I / MB) % NPROW;
    Pc = (J / NB) % NPCOL;
    iX = I % MB;
    iY = J % NB;

  }
  
  /**
   * \brief Calculate global (distributed) 2D coordinates from local
   *        coordinates using general BLACS grid and coordinates
   *
   * @param[in]  ICONTXT BLACS context
   * @param[in]  MB Distributed row block size
   * @param[in]  NB Distributed column block size
   * @param[in]  NPROW Number of rows in process grid
   * @param[in]  NPCOL Number of columnss in process grid
   * @param[out] I Distributed row coordinate
   * @param[out] J Distributed column coordinate
   * @param[in]  Pr Process row the contains distributed element (I,J)
   * @param[in]  Pc Process column the contains distributed element (I,J)
   * @param[in]  iX Local row coordinate
   * @param[in]  iY Local column coordinate
   */
  inline void GlobalFromLocal(const CB_INT ICONTXT, const CB_INT MB,
    const CB_INT NB, const CB_INT NPROW, const CB_INT NPCOL, CB_INT &I, 
    CB_INT &J, const CB_INT Pr, const CB_INT Pc, const CB_INT iX, 
    const CB_INT iY) {
  
    CB_INT L = iX / MB;
    CB_INT M = iY / NB;
  
    I = (L * (NPROW - 1) + Pr) * MB + iX;
    J = (M * (NPCOL - 1) + Pc) * NB + iY;

  }




  inline ScaLAPACK_Desc_t DescInit(const CB_INT M,
    const CB_INT N, const CB_INT MB, const CB_INT NB, const CB_INT ISRC,
    const CB_INT JSRC, const CB_INT ICTXT, const CB_INT LDD) {


    ScaLAPACK_Desc_t desc;

    CB_INT INFO;
    descinit_(&desc[0],&M,&N,&MB,&NB,&ISRC,&JSRC,&ICTXT,&LDD,&INFO);

    // Note that this is not always fatal, and useful information can
    // be obtained, such as smallest LDD
    if( INFO != 0 ) {

      std::stringstream ss;
      ss << "DESCINIT RECIEVED ILLEGAL ARG(" << -INFO << ")";
      std::runtime_error err(ss.str());

      throw err;

    };

    return desc;

  };



}; // namespace CXXBLACS

#endif

