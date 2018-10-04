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

template <typename T> T SmartConj(const T);
template<> inline double SmartConj(const double x){ return x; }
template<> inline std::complex<double> SmartConj( const std::complex<double>  x ){ return std::conj(x); }

BOOST_AUTO_TEST_SUITE(PPOTRF)

template <typename Field, typename RealType, CB_INT MB>
void ppotrf_test( CB_INT N ) {

  // Synchronize processes
  MPI_Barrier(MPI_COMM_WORLD);

  BlacsGrid grid(MPI_COMM_WORLD,MB,MB);

  std::vector<Field> A, ALoc, L;

  // Allocate local buffers
  CB_INT NLoc,MLoc;
  std::tie(MLoc,NLoc) = grid.getLocalDims(N,N);

  ALoc.resize(MLoc * NLoc);

  // Get DESC
  auto DescA = grid.descInit(N,N,0,0,MLoc);

  // Form Random symmetic matrix on root process
  RootExecute(MPI_COMM_WORLD,[&]() {

    A.resize(N*N);
    L.resize(N*N);

    for(auto i = 0; i < N; i++)
    for(auto j = 0; j <= i; j++) {

      A[i + j*N] = generate<Field>();
      A[j + i*N] = SmartConj(A[i + j*N]);

      if(i == j) A[j + i*N] = std::real(A[j+i*N]) + N;

    }

  });

  // Distribute to Grid 
  grid.Scatter(N,N,A.data(),N,ALoc.data(),MLoc,0,0);


  // Get LWORK
  CXXBLACS::PPOTRF('L',N,ALoc.data(),1,1,DescA);
  grid.printLocalBuffer(std::cout,N,N,ALoc.data(),MLoc);

  // Gather the eigenvectors to root process
  grid.Gather(N,N,L.data(),N,ALoc.data(),MLoc,0,0);


  // Check results on root process
  RootExecute(MPI_COMM_WORLD,[&](){

    std::vector<Field> LCpy(L);
    for(auto j = 0; j < N; j++)
    for(auto i = 0; i < j; i++)
      LCpy[i + j*N] = 0.;

    TRMM('R','L','C','N',N,N,Field(1.),L.data(),N,LCpy.data(),N);
    
    for(auto k = 0; k < N*N; k++) A[k] -= LCpy[k];

    std::vector<RealType> DIFF(N*N,0.);
    for(auto k = 0; k < N*N; k++) DIFF[k] = std::abs(A[k]);

    RealType maxDiff = *std::max_element(DIFF.begin(),DIFF.end());

    BOOST_CHECK_MESSAGE( maxDiff < 1e-10, 
      "MAX DIFF " << maxDiff
      << " " << std::numeric_limits<Field>::epsilon() 
    );

  });

  NotRootExecute(MPI_COMM_WORLD,[&](){ BOOST_CHECK(true); });


  // Synchronize processes
  MPI_Barrier(MPI_COMM_WORLD);

};


#define TEST_IMPL_F(NAME,F,RF,MB,N)\
  BOOST_AUTO_TEST_CASE(NAME) { ppotrf_test<F,RF,MB>(N); };

#define TEST_IMPL(NAME,MB,N)\
  TEST_IMPL_F(NAME##_Double, double, double, MB, N)\
  TEST_IMPL_F(NAME##_CDouble, std::complex<double>, double, MB, N)\
//FIXME: Need to determine proper error criteria for single precision
//TEST_IMPL_F(NAME##_Float , float , float, MB, N)\
//TEST_IMPL_F(NAME##_CFloat , std::complex<float> , float, MB, N)\


TEST_IMPL(PPOTRF_2x2,2,CXXBLACS_N);

BOOST_AUTO_TEST_SUITE_END()



