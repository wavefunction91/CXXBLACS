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
/**
 * \brief C++ Wrapper for BLACS_PINFO
 *
 * See BLACS Documentaion
 */
static inline void PINFO(int &IAM, int &NPROCS){
  blacs_pinfo_(&IAM,&NPROCS);
}

/**
 * \brief C++ Wrapper for BLACS_GET
 *
 * See BLACS Documentation.
 */
static inline void Get(const int ICONTXT, const int WHAT, int &VAL){
  blacs_get_(&ICONTXT,&WHAT,&VAL);
}

/**
 * \brief C++ Wrapper for BLACS_GRIDINIT
 *
 * See BLACS Documentation.
 */
static inline void GridInit(int &ICONTXT, const char ORDER[], 
  const int NPROW, const int NPCOL){
  blacs_gridinit_(&ICONTXT,ORDER,&NPROW,&NPCOL);
}

/**
 * \brief C++ Wrapper for BLACS_GRIDINFO
 *
 * See BLACS Documentation.
 */

static inline void GridInfo(const int ICONTXT, const int NPROW, 
  const int NPCOL, int &MYROW, int &MYCOL) {
  blacs_gridinfo_(&ICONTXT,&NPROW,&NPCOL,&MYROW,&MYCOL);
}

/**
 * \brief C++ Wrapper for BLACS_BARRIER
 *
 * See BLACS Documentaion
 */
static inline void Barrier(const int ICONTXT, const char SCOPE[]){
  blacs_barrier_(&ICONTXT,SCOPE);
}

/**
 * \brief C++ Wrapper for BLACS_GRIDEXIT
 *
 * See BLACS Documentaion
 */
static inline void GridExit(const int ICONTXT){
  blacs_gridexit_(&ICONTXT);
}

/**
 * \brief C++ Wrapper for BLACS_EXIT
 *
 * See BLACS Documentaion
 */
static inline void Exit(const int CONTINUE){
  blacs_exit_(&CONTINUE);
}

/**
 * \brief C++ Wrapper for NUMROC
 *
 * See BLACS/ScaLAPACK Documentaion
 */
static inline int NumRoc(const int N, const int NB, const int P, 
  const int SRC, const int NP) {
  numroc_(&N,&NB,&P,&SRC,&NP);
}

/**
 * \brief C++ Wrapper for ?GESD2D
 *
 * A templated function that encompasses the functionaliy of all of the BLACS
 * point-to-point send routines: [I/S/D/C/Z]GESD2D. Template deduction is based
 * on inut parameters
 *
 * TODO: Add implementation for Integer, Float, Complex and Complex Double
 *
 * See BLACS Documentaion for specifics.
 */
template<typename Field>
static inline void GESD2D(const int ICONTXT, const int M, const int N,
  Field *A, const int LDA, const int RDest, const int CDest){
  if(typeid(Field).hash_code() == typeid(double).hash_code())
    dgesd2d_(&ICONTXT,&M,&N,A,&LDA,&RDest,&CDest);
}

/**
 * \brief C++ Wrapper for ?GERV2D
 *
 * A templated function that encompasses the functionaliy of all of the BLACS
 * point-to-point receive routines: [I/S/D/C/Z]GERV2D. Template deduction is based
 * on inut parameters
 *
 * TODO: Add implementation for Integer, Float, Complex and Complex Double
 *
 * See BLACS Documentaion for specifics.
 */
template<typename Field>
static inline void GERV2D(const int ICONTXT, const int M, const int N,
  Field *A, const int LDA, const int RDest, const int CDest){
  if(typeid(Field).hash_code() == typeid(double).hash_code())
    dgerv2d_(&ICONTXT,&M,&N,A,&LDA,&RDest,&CDest);
}


template<typename Field>
static inline void GSUM2D(const int ICONTXT, const char SCOPE[], 
  const char TOP[], const int M, const int N, Field &A, const int LDA,
  const int RDest, const int CDest) {
  
  if(typeid(Field).hash_code() == typeid(double).hash_code())
    dgsum2d_(&ICONTXT,SCOPE,TOP,&M,&N,&A,&LDA,&RDest,&CDest);
}
