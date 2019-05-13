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
#ifndef __INCLUDED_CXXBLACS_BLACS_GRIDMANIP_HPP__
#define __INCLUDED_CXXBLACS_BLACS_GRIDMANIP_HPP__

#include <cxxblacs/proto.hpp>

namespace CXXBLACS {

  /**
   * \brief C++ Wrapper for BLACS_PINFO
   *
   * See BLACS Documentaion
   */
  inline void BlacsPINFO(CB_INT &IAM, CB_INT &NPROCS){
    Cblacs_pinfo(&IAM,&NPROCS);
  }



  /**
   * \brief C++ Wrapper for BLACS_GET
   *
   * See BLACS Documentation.
   */
  inline CB_INT BlacsGet(const CB_INT ICONTXT, const CB_INT WHAT){
    CB_INT VAL;
    Cblacs_get(ICONTXT,WHAT,&VAL);
    return VAL;
  }



  
  /**
   * \brief C++ Wrapper for BLACS_GRIDINIT
   *
   * See BLACS Documentation.
   */
  inline void BlacsGridInit(CB_INT &ICONTXT, const char ORDER[], 
    const CB_INT NPROW, const CB_INT NPCOL){
    Cblacs_gridinit(&ICONTXT,ORDER,NPROW,NPCOL);
  }



  
  /**
   * \brief C++ Wrapper for BLACS_GRIDINFO
   *
   * See BLACS Documentation.
   */
  
  inline void BlacsGridInfo(const CB_INT ICONTXT, CB_INT NPROW, 
    CB_INT NPCOL, CB_INT &MYROW, CB_INT &MYCOL) {
    Cblacs_gridinfo(ICONTXT,&NPROW,&NPCOL,&MYROW,&MYCOL);
  }



  
  /**
   * \brief C++ Wrapper for BLACS_BARRIER
   *
   * See BLACS Documentaion
   */
  inline void BlacsBarrier(const CB_INT ICONTXT, const char SCOPE[]){
    Cblacs_barrier(ICONTXT,SCOPE);
  }



  
  /**
   * \brief C++ Wrapper for BLACS_GRIDEXIT
   *
   * See BLACS Documentaion
   */
  inline void BlacsGridExit(const CB_INT ICONTXT){
    Cblacs_gridexit(ICONTXT);
  }




  
  /**
   * \brief C++ Wrapper for BLACS_EXIT
   *
   * See BLACS Documentaion
   */
  inline void BlacsExit(const CB_INT CONTINUE){
    Cblacs_exit(CONTINUE);
  }



}; // namespace CXXBLACS

#endif

