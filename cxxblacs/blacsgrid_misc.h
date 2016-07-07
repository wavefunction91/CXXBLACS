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
 * \brief Compute the number of rows and columns closest to a square process
 *        grid
 *
 *  @param[in]  NPROCS Number of MPI processes
 *  @param[out] NPROW Number of BLACS process rows
 *  @param[out] NPCOL Number of BLACS process columns
 */
static inline void SquareGrid(const int NPROCS, int &NPROW, int &NPCOL){
  NPROW = int(std::sqrt(NPROCS));
  NPCOL = NPROCS / NPROW;
}

/**
 * \brief Obtains the number of rows and columns of the local buffer of
 *        a distributed matrix. Wrapper for ScaLAPACK #NumRoc.
 *
 * @param [in] N Number of columns of distributed matrix
 * @param [in] M Number of rows of distributed matrix
 * @param [in] MB Distributed row block size
 * @param [in] NB Distributed column block size
 * @param [in] NPROW Number of rows in process grid
 * @param [in] NPCOL Number of columnss in process grid
 * @param [in] Pr Process row 
 * @param [in] Pc Process column 
 * @param [out] NLocR Number of rows of local buffer
 * @param [out] NLocC Number of columns of local buffer
 */
static inline void GetLocalDims(const int N, const int M, const int MB, 
  const int NB, const int NPROW, const int NPCOL, const int Pr,
  const int Pc, int &NLocR, int &NLocC) {

  NLocR = NumRoc(M,MB,Pr,0,NPROW);
  NLocC = NumRoc(N,NB,Pc,0,NPCOL);
}

/**
 * \brief Wrapper for GESD2D that uses the current grid
 *
 * @param [in] N Number of columns of buffer 
 * @param [in] M Number of rows of buffer
 * @param [in] A Buffer to send
 * @param [in] LDA Leading dimension of buffer
 * @param [in] RDest Destination process row
 * @param [in] CDest Destination process column
 */
template<typename Field>
inline void Send(const int M, const int N, Field *A, const int LDA,
  const int RDest, const int CDest) {

  BlacsGrid::GESD2D(this->IContxt_,M,N,A,LDA,RDest,CDest);
}

/**
 * \brief Wrapper for GERV2D that uses the current grid
 */
template<typename Field>
inline void Recv(const int M, const int N, Field *A, const int LDA,
  const int RDest, const int CDest) {

  BlacsGrid::GERV2D(this->IContxt_,M,N,A,LDA,RDest,CDest);
}

template<typename Field>
inline void PDOT(const int N, const Field *ALoc, const int AStride,
  const Field * BLoc, const int BStride, Field &LDDOT){

  int NLocR,NLocC;
  this->getLocalDims(N,N,NLocR,NLocC); 
  int NLOC = NLocR*NLocC;
  LDDOT = ddot_(&NLOC,ALoc,&AStride,BLoc,&BStride);
  BlacsGrid::GSUM2D(this->IContxt_,"All","1-tree",1,1,LDDOT,1,0,0);
}
