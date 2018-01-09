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

#include <ut.hpp>
#include <cxxblacs.hpp>

using namespace CXXBLACS;

constexpr CB_INT CXXBLACS_M = 20;
constexpr CB_INT CXXBLACS_N = 15;

namespace CXXBLACS {
  template <typename Field>
  inline Field generate(Field x);
  
  template <>
  inline double generate(double x){ return x; }
  
  template <>
  inline std::complex<double> generate(std::complex<double> x){ 
    return std::complex<double>(std::real(x),-std::real(x)); 
  }

  template <>
  inline float generate(float x){ return x; }
  
  template <>
  inline std::complex<float> generate(std::complex<float> x){ 
    return std::complex<float>(std::real(x),-std::real(x)); 
  }
}



BOOST_AUTO_TEST_SUITE(SCATTER)



template <typename Field, size_t MB, size_t NB, size_t M, size_t N>
void scatter_test() {

  // Synchronize processes
  MPI_Barrier(MPI_COMM_WORLD);

  BlacsGrid grid(MPI_COMM_WORLD,MB,NB);

  std::vector<Field> A, ALoc;

  // Form full matrix on root process
  RootExecute(MPI_COMM_WORLD,[&]() {
    A.resize(M*N);
    for(auto k = 0ul; k < M*N; k++) A[k] = generate(Field(k));
  });

  // Allocate Local Buffers
  CB_INT NLocR, NLocC;
  std::tie(NLocR, NLocC) = grid.getLocalDims(M,N);

  ALoc.resize(NLocR * NLocC);

  // Scatter the matrix onto the process grid
  grid.Scatter(M,N,A.data(),M,ALoc.data(),NLocR,0,0);

  // Test the local buffers
  for(auto iLocR = 0; iLocR < NLocR; iLocR++)
  for(auto iLocC = 0; iLocC < NLocC; iLocC++) {

    CB_INT I,J;
    std::tie(I,J) = grid.globalFromLocal(iLocR,iLocC);
    CB_INT K = I + J*M;

    BOOST_CHECK_MESSAGE(
      (std::abs(ALoc[iLocR + iLocC*NLocR] -  generate(Field(K)))) < 1e-16, 
      "Scattered Buffer Not Correct! " << 
      "(" << iLocR << ", " << iLocC << ") -> " <<
      "(" << I     << ", " << J     << "): " <<
      ALoc[iLocR + iLocC*NLocR] << ", " << K 
    );

  }

  // Synchronize processes
  MPI_Barrier(MPI_COMM_WORLD);

};

#define TEST_IMPL_F(NAME,F,MB,NB,M,N)\
  BOOST_AUTO_TEST_CASE(NAME) { scatter_test<F,MB,NB,M,N>(); };

#define TEST_IMPL(NAME,MB,NB,M,N) \
  TEST_IMPL_F(NAME##_Float,        float,               MB,NB,M,N)\
  TEST_IMPL_F(NAME##_ComplexFloat, std::complex<float>, MB,NB,M,N)\
  TEST_IMPL_F(NAME##_Double,       double,              MB,NB,M,N)\
  TEST_IMPL_F(NAME##_ComplexDouble,std::complex<double>,MB,NB,M,N)


TEST_IMPL(Scatter_2x2_SquareMatrix,2,2,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(Scatter_1x2_SquareMatrix,1,2,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(Scatter_2x1_SquareMatrix,2,1,CXXBLACS_N,CXXBLACS_N);

TEST_IMPL(Scatter_2x2_RectangularMatrix_BS,2,2,CXXBLACS_M,CXXBLACS_N);
TEST_IMPL(Scatter_1x2_RectangularMatrix_BS,1,2,CXXBLACS_M,CXXBLACS_N);
TEST_IMPL(Scatter_2x1_RectangularMatrix_BS,2,1,CXXBLACS_M,CXXBLACS_N);

TEST_IMPL(Scatter_2x2_RectangularMatrix_SB,2,2,CXXBLACS_N,CXXBLACS_M);
TEST_IMPL(Scatter_1x2_RectangularMatrix_SB,1,2,CXXBLACS_N,CXXBLACS_M);
TEST_IMPL(Scatter_2x1_RectangularMatrix_SB,2,1,CXXBLACS_N,CXXBLACS_M);


BOOST_AUTO_TEST_SUITE_END()
