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
 * \brief Gather a distributed matrix to the root process
 *
 * Note: The gathered buffer "A" need only be allocated on root process
 *
 * @param [in]  N Number of rows / columns of distributed matrix
 * @param [out] A Buffer for gathered matrix on root process
 * @param [in]  ALoc Local buffer of distributed matrix A
 */
inline void gather(const int N, std::vector<double> &A, 
  std::vector<double> &ALoc){

  int NLocR,NLocC;
  this->getLocalDims(N,N,NLocR,NLocC);

  BlacsGrid::RootExecute([&] {
    int I,J;
    for(auto iLocR = 0; iLocR < NLocR; iLocR++)
    for(auto iLocC = 0; iLocC < NLocC; iLocC++) {
      this->globalFromLocal(I,J,iLocR,iLocC);
      A[J*N + I] = ALoc[iLocR + iLocC*NLocR];
    }

    for(auto iPc = 0; iPc < this->nProcCol_; iPc++)
    for(auto iPr = 0; iPr < this->nProcRow_; iPr++){
      if(iPr == 0 && iPc == 0) continue;
      BlacsGrid::GetLocalDims(N,N,this->mb_,this->nb_,this->nProcRow_,
        this->nProcCol_,iPr,iPc,NLocR,NLocC);

      this->Recv(NLocR,NLocC,ALoc.data(),NLocR,iPr,iPc);
      for(auto iLocR = 0; iLocR < NLocR; iLocR++) 
      for(auto iLocC = 0; iLocC < NLocC; iLocC++) {
        BlacsGrid::GlobalFromLocal(this->IContxt_,this->mb_,this->nb_,
          this->nProcRow_,this->nProcCol_,I,J,iPr,iPc,iLocR,iLocC);
        A[J*N + I] = ALoc[iLocR + iLocC*NLocR];
      }
    }
  });

  BlacsGrid::NotRootExecute( [&] {
    this->Send(NLocR,NLocC,ALoc.data(),NLocR,0,0);
  });

  BlacsGrid::Barrier(this->IContxt_,"All");

};

/**
 * \brief Distribute a matrix from root process to the BLACS grid
 *
 * Note: The gathered buffer "A" need only be allocated on root process
 *
 * @param [in]   N Number of rows / columns of distributed matrix
 * @param [in]   A Buffer for gathered matrix on root process
 * @param [out]  ALoc Local buffer of distributed matrix A
 */
inline void distribute(const int N, std::vector<double> &A,
  std::vector<double> &ALoc){

  int NLocR,NLocC;
  this->getLocalDims(N,N,NLocR,NLocC);

  BlacsGrid::RootExecute([&] {
    int I,J;
    
    // Send local parts of matricies out to BLACS grid
    for(auto iPc = 0; iPc < this->nProcCol_; iPc++)
    for(auto iPr = 0; iPr < this->nProcRow_; iPr++){
      if(iPr == 0 && iPc == 0) continue;
      // Get local variables for BLACS coordinate (iPr,iPc)
      BlacsGrid::GetLocalDims(N,N,this->mb_,this->nb_,this->nProcRow_,
        this->nProcCol_,iPr,iPc,NLocR,NLocC);

      // Copy the local parts of the matrix to ALoc (tmp buffer)
      for(auto iLocR = 0; iLocR < NLocR; iLocR++)
      for(auto iLocC = 0; iLocC < NLocC; iLocC++) {
        BlacsGrid::GlobalFromLocal(this->IContxt_,this->mb_,this->nb_,
          this->nProcRow_,this->nProcCol_,I,J,iPr,iPc,iLocR,iLocC);
        ALoc[iLocR + iLocC*NLocR] = A[J*N + I];
      }

      // Send buffer
      this->Send(NLocR,NLocC,ALoc.data(),NLocR,iPr,iPc);
    }

    // Take care of local buffer
    this->getLocalDims(N,N,NLocR,NLocC);
    for(auto iLocR = 0; iLocR < NLocR; iLocR++)
    for(auto iLocC = 0; iLocC < NLocC; iLocC++) {
      this->globalFromLocal(I,J,iLocR,iLocC);
      ALoc[iLocR + iLocC*NLocR] = A[J*N + I];
    }
  });

  // Recieve local buffer
  BlacsGrid::NotRootExecute( [&] {
    this->Recv(NLocR,NLocC,ALoc.data(),NLocR,0,0);
  });
  
};
