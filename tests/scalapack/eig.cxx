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

BOOST_AUTO_TEST_SUITE(PSYEV)

template <typename Field, CB_INT MB>
void psyev_test( CB_INT N ) {

  // Synchronize processes
  MPI_Barrier(MPI_COMM_WORLD);

  BlacsGrid grid(MPI_COMM_WORLD,MB,MB);

  std::vector<Field> A, ALoc, Z, ZLoc, W;

  // Allocate local buffers
  CB_INT NLoc,MLoc;
  std::tie(MLoc,NLoc) = grid.getLocalDims(N,N);

  ALoc.resize(MLoc * NLoc);
  ZLoc.resize(MLoc * NLoc);
  W.resize(N);


  // Get DESC
  auto DescA = grid.descInit(N,N,0,0,MLoc);

  // Form Random symmetic matrix on root process
  RootExecute(MPI_COMM_WORLD,[&]() {

    A.resize(N*N);
    Z.resize(N*N);
    for(auto i = 0; i < N; i++)
    for(auto j = 0; j <= i; j++) {

      A[i + j*N] = generate<Field>();
      A[j + i*N] = A[i + j*N];

    }

  });

  // Distribute to Grid 
  grid.Scatter(N,N,A.data(),N,ALoc.data(),MLoc,0,0);


  // Get LWORK
  CB_INT LWORK = -1;
  std::vector<Field> WORK(5);

  CXXBLACS::PSYEV('V','U',N,ALoc.data(),1,1,&DescA[0],W.data(),
    ZLoc.data(),1,1,&DescA[0],WORK.data(),LWORK);

  LWORK = CB_INT(WORK[0]);
  WORK.resize(LWORK);

  // Diagonalize the matrix
  CXXBLACS::PSYEV('V','U',N,ALoc.data(),1,1,&DescA[0],W.data(),
    ZLoc.data(),1,1,&DescA[0],WORK.data(),LWORK);



  // Gather the eigenvectors to root process
  grid.Gather(N,N,Z.data(),N,ZLoc.data(),MLoc,0,0);


  // Check results on root process
  RootExecute(MPI_COMM_WORLD,[&](){

    std::vector<Field> TMP(Z);
    GEMM('N','N',N,N,N,1.,A.data(),N,Z.data(),N,0.,TMP.data(),N);
    GEMM('C','N',N,N,N,1.,Z.data(),N,TMP.data(),N,0.,A.data(),N);
    
    for(auto k = 0; k < N; k++) A[k*(N+1)] -= W[k];

    std::vector<Field> DIFF(N*N,0.);
    for(auto k = 0; k < N*N; k++) DIFF[k] = std::abs(A[k]);

    Field maxDiff = *std::max_element(DIFF.begin(),DIFF.end());

    BOOST_CHECK_MESSAGE( maxDiff < 1e-10, 
      "MAX DIFF " << maxDiff
      << " " << std::numeric_limits<Field>::epsilon() 
      << " " << std::numeric_limits<Field>::epsilon() * N*N*N 
    );

  });

  NotRootExecute(MPI_COMM_WORLD,[&](){ BOOST_CHECK(true); });


  // Synchronize processes
  MPI_Barrier(MPI_COMM_WORLD);

};


#define PSYEV_TEST_IMPL_F(NAME,F,MB,N)\
  BOOST_AUTO_TEST_CASE(NAME) { psyev_test<F,MB>(N); };

#define PSYEV_TEST_IMPL(NAME,MB,N)\
  PSYEV_TEST_IMPL_F(NAME##_Double, double, MB, N)\
//FIXME: Need to determine proper error criteria for single precision
//PSYEV_TEST_IMPL_F(NAME##_Float , float , MB, N)\


PSYEV_TEST_IMPL(PSYEV_2x2,2,CXXBLACS_N);

BOOST_AUTO_TEST_SUITE_END()



BOOST_AUTO_TEST_SUITE(PHEEV)

template <typename Field, CB_INT MB, 
  typename RealType = typename Field::value_type>
void pheev_test( CB_INT N ) {

  // Synchronize processes
  MPI_Barrier(MPI_COMM_WORLD);

  BlacsGrid grid(MPI_COMM_WORLD,MB);

  std::vector<Field> A, ALoc, Z, ZLoc;
  std::vector<RealType> W;

  // Allocate local buffers
  CB_INT NLoc,MLoc;
  std::tie(MLoc,NLoc) = grid.getLocalDims(N,N);

  ALoc.resize(MLoc * NLoc);
  ZLoc.resize(MLoc * NLoc);
  W.resize(N);


  // Get DESC
  auto DescA = grid.descInit(N,N,0,0,MLoc);

  // Form Random hermetian matrix on root process
  RootExecute(MPI_COMM_WORLD,[&]() {

    A.resize(N*N);
    Z.resize(N*N);
    for(auto i = 0; i < N; i++)
    for(auto j = 0; j <= i; j++) {

      A[i + j*N] = generate<Field>();
      A[j + i*N] = std::conj(A[i + j*N]);
      if( i == j ) A[i + j*N] = A[i + j*N].real();

    }

  });

  // Distribute to Grid 
  grid.Scatter(N,N,A.data(),N,ALoc.data(),MLoc,0,0);


  // Real workspace
  CB_INT LRWORK = 4*N-2;
  std::vector<RealType> RWORK(LRWORK);

  // Get LWORK and allocate WORK
  CB_INT LWORK  = -1;
  std::vector<Field>    WORK(5);

  CXXBLACS::PHEEV('V','U',N,ALoc.data(),1,1,&DescA[0],W.data(),
    ZLoc.data(),1,1,&DescA[0],WORK.data(),LWORK,RWORK.data(),
    LRWORK);

  LWORK  = CB_INT(std::real(WORK[0]));
  WORK .resize(LWORK );

  // Diagonalize the matrix
  CXXBLACS::PHEEV('V','U',N,ALoc.data(),1,1,&DescA[0],W.data(),
    ZLoc.data(),1,1,&DescA[0],WORK.data(),LWORK,RWORK.data(),
    LRWORK);



  // Gather the eigenvectors to root process
  grid.Gather(N,N,Z.data(),N,ZLoc.data(),MLoc,0,0);


  // Check results on root process
  RootExecute(MPI_COMM_WORLD,[&](){

    std::vector<Field> TMP(Z);
    GEMM('N','N',N,N,N,Field(1.),A.data(),N,Z.data(),N,Field(0.),TMP.data(),N);
    GEMM('C','N',N,N,N,Field(1.),Z.data(),N,TMP.data(),N,Field(0.),A.data(),N);
    
    for(auto k = 0; k < N; k++) A[k*(N+1)] -= W[k];

    std::vector<RealType> DIFF(N*N,0.);
    for(auto k = 0; k < N*N; k++) DIFF[k] = std::abs(A[k]);

    RealType maxDiff = *std::max_element(DIFF.begin(),DIFF.end());

    BOOST_CHECK_MESSAGE( maxDiff < 1e-10, 
      "MAX DIFF " << maxDiff
      << " " << std::numeric_limits<Field>::epsilon() 
      << " " << std::numeric_limits<Field>::epsilon() * Field(N*N*N)
    );

  });

  NotRootExecute(MPI_COMM_WORLD,[&](){ BOOST_CHECK(true); });


  // Synchronize processes
  MPI_Barrier(MPI_COMM_WORLD);

};


#define PHEEV_TEST_IMPL_F(NAME,F,MB,N)\
  BOOST_AUTO_TEST_CASE(NAME) { pheev_test<F,MB>(N); };

#define PHEEV_TEST_IMPL(NAME,MB,N)\
  PHEEV_TEST_IMPL_F(NAME##_CDouble, std::complex<double>, MB, N)\
//FIXME: Need to determine proper error criteria for single precision
//PHEEV_TEST_IMPL_F(NAME##_CFloat , std::complex<float> , MB, N)\


PHEEV_TEST_IMPL(PHEEV_2x2,2,CXXBLACS_N);


BOOST_AUTO_TEST_SUITE_END()
