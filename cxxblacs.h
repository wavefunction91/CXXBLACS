/*
 *  A simple C++ Wrapper for BLACS along with minimal extra functionality to 
 *  aid the the high-level development of distributed memory linear algebra.
 *  Copyright (C) 2016 David Williams-Young

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
#ifndef __INCLUDED_CXXBLACS_H__
#define __INCLUDED_CXXBLACS_H__
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <mpi.h>
#include <cmath>
#include <functional>

extern "C" {
  // BLACS
  void blacs_pinfo_(int*,int*);
  void blacs_get_(const int*,const int*,int*);
  void blacs_gridinit_(int*,const char*,const int*,const int*);
  void blacs_gridinfo_(const int*,const int*,const int*,int*,int*);
  void blacs_barrier_(const int*,const char*);
  void blacs_gridexit_(const int*);
  void blacs_exit_(const int*);
  int  blacs_pnum_(int*,int*,int*);
  void dgesd2d_(const int*, const int *, const int *, const double *,
    const int*, const int *, const int *);
  void dgerv2d_(const int*, const int *, const int *, const double *,
    const int*, const int *, const int *);

  // SCALAPACK
  int  numroc_(const int*,const int*,const int*,const int*,const int*);

};


/**
 * \brief A class to wrap the FORTRAN functionality of the BLACS (Basic Linear Algebra Communication
 *        Subroutines) and provide minimal extra functionality to aid in development.
 *
 * BlacsGrid is a class that is designed to wrap the FORTRAN routines of the BLACS (Basic Linear
 * Algebra Communication Subroutines) that are used in conjunction with ScaLAPACK (Scalable Linear
 * Albegra PACKage) and PBLAS (Parallel Basic Linear Algebra Subroutines). In addition to providing
 * simple C++ wrappers to the BLACS routines, BlacsGrid also provides extra (but minimal) functionality
 * to aid the in the development of distributed memory linear algebra programs, such as:
 *   - Basic MPI execution design (such as ring execution)
 *   - Gathering a distributed matrix from the BLACS grid to root process
 *   - Distributing a matrix from a root process to the BLACS grid
 *   - Conversion of local and distrubuted memory element coordinates
 *
 * 
 */
class BlacsGrid {
  // Global MPI information
  int iProc_; ///< Current MPI process (rank)
  int nProc_; ///< Size of MPI universe (size)

  // BLACS grid information
  int IContxt_;  ///< BLACS context
  int nProcRow_; ///< Number of rows on BLACS grid
  int nProcCol_; ///< Number of columns on BLACS grid
  int iProcRow_; ///< Current row coordinate of BLACS grid
  int iProcCol_; ///< Current column coordinate of BLACS grid
  int mb_;       ///< Row block size (default = 2)
  int nb_;       ///< Column block size (default = 1)

public:
  /**
   * \brief Constructor
   *
   * Initialize the BLACS grid. Attempts to create as close to a square
   * BLACS grid as possible. Obtains MPI and BLACS information about the
   * current proces, such as iProc, nProc, iProcRow, etc... The default
   * BLACS process grid is taken to be row-major
   *
   *   @param[in] MB Block size for row distribution
   *   @param[in] NB Block size for column distribution
   */
  BlacsGrid(int mb = 2, int nb = 2) : IContxt_(0), mb_(mb), nb_(nb){
    BlacsGrid::PINFO(this->iProc_,this->nProc_);
    BlacsGrid::SquareGrid(this->nProc_,this->nProcRow_,this->nProcCol_);
    BlacsGrid::Get(0,0,this->IContxt_);
    BlacsGrid::GridInit(this->IContxt_,"row-major",this->nProcRow_,
      this->nProcCol_);
    BlacsGrid::GridInfo(this->IContxt_,this->nProcRow_,this->nProcCol_,
      this->iProcRow_,this->iProcCol_);
  };

  /**
   * \brief Destructor
   *
   * Properly destroy the BLACS grid and cleanup the MPI environment
   */
  ~BlacsGrid() {
    BlacsGrid::GridExit(this->IContxt_);
    BlacsGrid::Exit(0);
  };

  // Getters 
  /** Access to private #iProc_ variable */
  inline int iProc()    { return this->iProc_;    };
  /** Access to private #iProcRow_ variable */
  inline int iProcRow() { return this->iProcRow_; };
  /** Access to private #iProcCol_ variable */
  inline int iProcCol() { return this->iProcCol_; };
  /** Access to private #nProc_ variable */
  inline int nProc()    { return this->nProc_;    };
  /** Access to private #nProcRow_ variable */
  inline int nProcRow() { return this->nProcRow_; };
  /** Access to private #nProcCol_ variable */
  inline int nProcCol() { return this->nProcCol_; };
  /** Access to private #IContxt_ variable */
  inline int iContxt()  { return this->IContxt_;  };
  /** Access to private #nb_ variable */
  inline int NB()       { return this->nb_;       };
  /** Access to private #mb_ variable */
  inline int MB()       { return this->mb_;       };

  // Setters
  /**
   * \brief Set the block size for distribution
   *   @param[in] NB The blocksize for both row and columns
   */
  inline void setBlockSize(int NB) { 
    this->nb_ = NB; this->mb_ = NB;
  };
  /**
   * \brief Set the block size for distribution
   *   @param[in] MB The blocksize for rows
   *   @param[in] NB The blocksize for columns
   */
  inline void setBlockSize(int MB, int NB) {
    this->nb_ = NB; this->mb_ = MB;
  }

  #include "cxxblacs/blacs_print.h"
  #include "cxxblacs/blacs_wrappers.h"
  #include "cxxblacs/mpi_exec.h"
  #include "cxxblacs/global_local_conv.h"
  #include "cxxblacs/blacsgrid_misc.h"
  #include "cxxblacs/blacs_gather_distribute.h"

};
#endif
