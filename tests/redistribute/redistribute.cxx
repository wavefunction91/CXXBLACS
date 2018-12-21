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

#include "redistribute.hpp"

template <typename Field, size_t M, size_t N>
void redistribute_test() {

  // Synchronize processes
  MPI_Barrier(MPI_COMM_WORLD);

  BlacsGrid grid_0(MPI_COMM_WORLD,2,2);
  BlacsGrid grid_1(MPI_COMM_WORLD,2,1);
  BlacsGrid grid_2(MPI_COMM_WORLD,1,2);
  BlacsGrid grid_3(MPI_COMM_WORLD,1,1);

  std::vector<Field> A, ALoc0, ALoc1, ALoc2, ALoc3;

  // Allocate local buffers
  CB_INT N0,M0,N1,M1,N2,M2,N3,M3;
  std::tie(M0,N0) = grid_0.getLocalDims(M,N);
  std::tie(M1,N1) = grid_1.getLocalDims(M,N);
  std::tie(M2,N2) = grid_2.getLocalDims(M,N);
  std::tie(M3,N3) = grid_3.getLocalDims(M,N);

  ALoc0.resize(M0 * N0);
  ALoc1.resize(M1 * N1);
  ALoc2.resize(M2 * N2);
  ALoc3.resize(M3 * N3);


  // Get DESC
  auto DescA0 = grid_0.descInit(M,N,0,0,M0);
  auto DescA1 = grid_1.descInit(M,N,0,0,M1);
  auto DescA2 = grid_2.descInit(M,N,0,0,M2);
  auto DescA3 = grid_3.descInit(M,N,0,0,M3);

  // Form full matrix on root process
  RootExecute(MPI_COMM_WORLD,[&]() {
    A.resize(M*N);
    for(auto k = 0ul; k < M*N; k++) A[k] = generate(Field(k));
  });

  // Distribute to Grid 0
  grid_0.Scatter(M,N,A.data(),M,ALoc0.data(),M0,0,0);

  // Redistribute Grid 0 -> Grid 1
  PGEMR2D(M,N,ALoc0.data(),1,1,DescA0,ALoc1.data(),1,1,DescA1,grid_0.iContxt());

  // Redistribute Grid 0 -> Grid 2
  PGEMR2D(M,N,ALoc0.data(),1,1,DescA0,ALoc2.data(),1,1,DescA2,grid_0.iContxt());

  // Redistribute Grid 0 -> Grid 3
  PGEMR2D(M,N,ALoc0.data(),1,1,DescA0,ALoc3.data(),1,1,DescA3,grid_0.iContxt());

  auto check = [&](BlacsGrid & grid, Field* ALoc, CB_INT NLocR, CB_INT NLocC) {
    for(auto iLocR = 0; iLocR < NLocR; iLocR++)
    for(auto iLocC = 0; iLocC < NLocC; iLocC++) {

      CB_INT I,J;
      std::tie(I,J) = grid.globalFromLocal(iLocR,iLocC);
      CB_INT K = I + J*M;

      EXPECT_NEAR( std::abs(ALoc[iLocR + iLocC*NLocR]),  std::abs(generate(Field(K))),  1e-16 ) <<
        "Redistributed Buffer Not Correct! " << 
        "(" << iLocR << ", " << iLocC << ") -> " <<
        "(" << I     << ", " << J     << "): " <<
        ALoc[iLocR + iLocC*NLocR] << ", " << K << std::endl; 

    }

  };


  // Test the local buffers
  check(grid_0,ALoc0.data(),M0,N0);
  check(grid_1,ALoc1.data(),M1,N1);
  check(grid_2,ALoc2.data(),M2,N2);
  check(grid_3,ALoc3.data(),M3,N3);

  // Synchronize processes
  MPI_Barrier(MPI_COMM_WORLD);

};

#define TEST_IMPL_F(NAME,F,M,N)\
  TEST(REDIST,NAME) { redistribute_test<F,M,N>(); };

#define TEST_IMPL(NAME,M,N) \
  TEST_IMPL_F(NAME##_Float,        float,               M,N)\
  TEST_IMPL_F(NAME##_ComplexFloat, std::complex<float>, M,N)\
  TEST_IMPL_F(NAME##_Double,       double,              M,N)\
  TEST_IMPL_F(NAME##_ComplexDouble,std::complex<double>,M,N)


TEST_IMPL(Redistribute_SquareMatrix,CXXBLACS_N,CXXBLACS_N);

