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

#include "scalapack_ut.hpp"

BOOST_AUTO_TEST_SUITE(PGESV)

template <typename T, typename U>
using Diag_t = 
  std::function<void(const CB_INT,T*,const CB_INT,const CB_INT,
      const ScaLAPACK_Desc_t,U*,T*,const CB_INT,const CB_INT,
      const ScaLAPACK_Desc_t)>;


template <typename Field, CB_INT MB, typename RealType>
void pgesv_test(const CB_INT N, const CB_INT NRHS, 
  const Diag_t<Field,RealType> &diag) {

  // Synchronize processes
  MPI_Barrier(MPI_COMM_WORLD);

  BlacsGrid grid(MPI_COMM_WORLD,MB);

  std::vector<Field> A, ALoc, ZLoc, B, BLoc;
  std::vector<RealType> W;


  // Allocate local buffers
  CB_INT NLoc,MLoc;
  CB_INT BMLoc, NRHSLoc;
  std::tie(MLoc,NLoc)     = grid.getLocalDims(N,N);
  std::tie(BMLoc,NRHSLoc) = grid.getLocalDims(N,NRHS);

  ALoc.resize(MLoc  * NLoc);
  ZLoc.resize(MLoc  * NLoc);
  BLoc.resize(BMLoc * NRHSLoc);
  W.resize(N);

  // Get DESCA and DESCB
  auto DescA = grid.descInit(N,N   ,0,0,MLoc );
  auto DescB = grid.descInit(N,NRHS,0,0,BMLoc);

  // Form Random hermetian matrix A and RHS B on root process
  RootExecute(MPI_COMM_WORLD,[&]() {

    A.resize(N*N);
    B.resize(N*NRHS);
    for(auto i = 0; i < N; i++)
    for(auto j = 0; j <= i; j++) {

      A[i + j*N] = generate<Field>();
      A[j + i*N] = std::conj(A[i + j*N]);
      if( i == j ) A[i + j*N] = std::real(A[i + j*N]);

    }

    for(auto &X : B ) X = generate<Field>();

  });


  // Distribute A and B to Grid
  grid.Scatter(N,N   ,A.data(),N,ALoc.data(),MLoc ,0,0);
  grid.Scatter(N,NRHS,B.data(),N,BLoc.data(),BMLoc,0,0);

  // Diagonalize the matrix AZ = ZW
  diag(N,ALoc.data(),1,1,DescA,W.data(),ZLoc.data(),1,1,DescA);

  // Make a copy of B and Z for linear system
  std::vector<Field> ZLocCpy(ZLoc);
  std::vector<Field> BLocCpy(BLoc);

  // Allocate IPIV
  std::vector<CB_INT> IPIV(MLoc + MB,0);


  // Solve ZX = B
  CXXBLACS::PGESV(N,NRHS,ZLocCpy.data(),1,1,DescA,&IPIV[0],BLocCpy.data(),1,1,
      DescB);

  // Get analytic solution X = Z**H * B
  PGEMM('C','N',N,NRHS,N,Field(1.),ZLoc.data(),1,1,DescA,BLoc.data(),1,1,DescB,
    Field(0.), ZLocCpy.data(),1,1,DescA);


  // Compute differences on each MPI process
  for(auto j = 0; j < NRHSLoc; j++)
  for(auto i = 0; i < BMLoc;   i++)
    BLocCpy[i + j*BMLoc] -= ZLocCpy[i + j*MLoc];


  // Get max difference on each MPI process
  RealType maxDiff = std::abs(*std::max_element(BLocCpy.begin(), BLocCpy.end(),
      [&](Field x, Field y){ return std::abs(x) < std::abs(y); }));

  BOOST_CHECK_MESSAGE( maxDiff < 1e-10, "MAX DIFF " << maxDiff );

};

template <typename T>
void DIAG(const CB_INT N, T* A, const CB_INT IA, const CB_INT JA,
  const ScaLAPACK_Desc_t DescA, T* W, T* Z, const CB_INT IZ, const CB_INT JZ,
  const ScaLAPACK_Desc_t DescZ) {

  // Get LWORK
  CB_INT LWORK = -1;
  std::vector<T> WORK(5);

  PSYEV('V','U',N,A,1,1,DescA,W,Z,1,1,DescZ,WORK.data(),LWORK);

  // Resize WORK array
  LWORK = CB_INT(WORK[0]);
  WORK.resize(LWORK);

  // Diagonalize the matrix
  PSYEV('V','U',N,A,1,1,DescA,W,Z,1,1,DescZ,WORK.data(),LWORK);

};

template <typename T, typename RT = typename T::value_type>
void DIAGC(const CB_INT N, T* A, const CB_INT IA, const CB_INT JA,
  const ScaLAPACK_Desc_t DescA, RT* W, T* Z, const CB_INT IZ, const CB_INT JZ, 
  const ScaLAPACK_Desc_t DescZ) {

  // Real workspace
  CB_INT LRWORK = 4*N-2;
  std::vector<RT> RWORK(LRWORK);

  // Get LWORK
  CB_INT LWORK = -1;
  std::vector<T> WORK(5);

  PHEEV('V','U',N,A,1,1,DescA,W,Z,1,1,DescZ,WORK.data(),LWORK,
    RWORK.data(),LRWORK);

  // Resize WORK array
  LWORK = CB_INT(std::real(WORK[0]));
  WORK.resize(LWORK);

  // Diagonalize the matrix
  PHEEV('V','U',N,A,1,1,DescA,W,Z,1,1,DescZ,WORK.data(),LWORK,
    RWORK.data(),LRWORK);

};

#define TEST_IMPL_F(NAME,F,RF,DIAGF,MB,N,NRHS) \
  BOOST_AUTO_TEST_CASE(NAME) { pgesv_test<F,MB,RF>(N,NRHS,DIAGF<F>); };

#define TEST_IMPL(NAME,MB,N,NRHS) \
  TEST_IMPL_F(NAME##_Double,double,double,DIAG,MB,N,NRHS);\
  TEST_IMPL_F(NAME##_CDouble,std::complex<double>,double,DIAGC,MB,N,NRHS);
//FIXME: Need single precision tests

TEST_IMPL(PGESV_2x2,2,CXXBLACS_N,CXXBLACS_NRHS);


BOOST_AUTO_TEST_SUITE_END()

