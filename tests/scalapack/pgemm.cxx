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




template <typename Field, CB_INT MB, CB_INT NB, typename RealType = Field>
void pgemm_test(char TRANSA, char TRANSB, Field ALPHA, Field BETA,
    CB_INT M, CB_INT N, CB_INT K) {

  // Synchronize processes
  MPI_Barrier(MPI_COMM_WORLD);

  BlacsGrid grid(MPI_COMM_WORLD,MB,NB);


  // Make sure that all of the ALPHA's and BETA's are the same
  grid.Broadcast("All","I",1,1,&ALPHA,1,0,0);
  grid.Broadcast("All","I",1,1,&BETA ,1,0,0);

  std::vector<Field> A, ALoc, B, BLoc, C, CLoc, TrueAns;

  // Get local dims and allocate local buffers
  CB_INT ANLocR, ANLocC;
  CB_INT BNLocR, BNLocC;
  CB_INT CNLocR, CNLocC;
  std::tie(ANLocR,ANLocC) = grid.getLocalDims(M,K);
  std::tie(BNLocR,BNLocC) = grid.getLocalDims(K,N);
  std::tie(CNLocR,CNLocC) = grid.getLocalDims(M,N);

  ALoc.resize(ANLocR*ANLocC);
  BLoc.resize(BNLocR*BNLocC);
  CLoc.resize(CNLocR*CNLocC);





  // Generate random matricies on root process and get serial product
  RootExecute(MPI_COMM_WORLD,[&](){

    A.resize(M*K);
    B.resize(K*N);
    C.resize(M*N);
    TrueAns.resize(M*N);

    for(auto &X : A ) X = generate<Field>();
    for(auto &X : B ) X = generate<Field>();
    for(auto &X : C ) X = generate<Field>();

    std::copy(C.begin(),C.end(),TrueAns.begin());

    GEMM(TRANSA,TRANSB,M,N,K,ALPHA,&A[0],M,&B[0],K,BETA,&TrueAns[0],M);

  });




  // Scatter the matrices onto the process grid
  grid.Scatter(M,K,A.data(),M,ALoc.data(),ANLocR,0,0);
  grid.Scatter(K,N,B.data(),K,BLoc.data(),BNLocR,0,0);
  grid.Scatter(M,N,C.data(),M,CLoc.data(),CNLocR,0,0);




  // Get ScaLAPACK descriptiors
  auto DESCA = grid.descInit(M,K,0,0,ANLocR);
  auto DESCB = grid.descInit(K,N,0,0,BNLocR);
  auto DESCC = grid.descInit(M,N,0,0,CNLocR);


  // Get Product
  CXXBLACS::PGEMM(TRANSA,TRANSB,M,N,K,ALPHA,
    ALoc.data(),1,1,DESCA,
    BLoc.data(),1,1,DESCB,BETA,
    CLoc.data(),1,1,DESCC);


  // Gather product to root process
  grid.Gather(M,N,C.data(),M,CLoc.data(),CNLocR,0,0);

  // Test result on root process
  RootExecute(MPI_COMM_WORLD,[&](){

    std::vector<RealType> DIFF(TrueAns.size(),0.);
    for(auto k = 0; k < DIFF.size(); k++) 
      DIFF[k] = std::abs(TrueAns[k] - C[k]);

    RealType maxDiff = *std::max_element(DIFF.begin(),DIFF.end());

    std::cout << "MAX DIFF " << maxDiff << std::endl;


    EXPECT_NEAR( maxDiff, 0., 1.e-10 ) <<
      "MAX DIFF " << maxDiff
      << " " << std::numeric_limits<RealType>::epsilon() 
      << " " << std::numeric_limits<RealType>::epsilon() * M*N*K ;

  });

  NotRootExecute(MPI_COMM_WORLD,[&](){ EXPECT_TRUE(true); });

  // Synchronize processes
  MPI_Barrier(MPI_COMM_WORLD);

};

#define TEST_IMPL_F(NAME,F,RF,TA,TB,MB,NB,M,N,K)\
  TEST(PGEMM,NAME) { \
    pgemm_test<F,MB,NB,RF>(TA,TB,generate<F>(),generate<F>(),M,N,K); \
  };

#define TEST_IMPL(NAME,TA,TB,MB,NB,M,N,K) \
  TEST_IMPL_F(NAME##_Double, double,              double,TA,TB,MB,NB,M,N,K)\
  TEST_IMPL_F(NAME##_CDouble,std::complex<double>,double,TA,TB,MB,NB,M,N,K)\
//FIXME: Need to determine proper error criteria for single precision
//TEST_IMPL_F(NAME##_Float,  float,               float, TA,TB,MB,NB,M,N,K)\
//TEST_IMPL_F(NAME##_CFloat, std::complex<float>, float, TA,TB,MB,NB,M,N,K)




TEST_IMPL(PGEMM_2x2_Square_NN,'N','N',2,2,CXXBLACS_N,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PGEMM_2x1_Square_NN,'N','N',2,1,CXXBLACS_N,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PGEMM_1x2_Square_NN,'N','N',1,2,CXXBLACS_N,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PGEMM_2x2_Square_TN,'T','N',2,2,CXXBLACS_N,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PGEMM_2x1_Square_TN,'T','N',2,1,CXXBLACS_N,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PGEMM_1x2_Square_TN,'T','N',1,2,CXXBLACS_N,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PGEMM_2x2_Square_CN,'C','N',2,2,CXXBLACS_N,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PGEMM_2x1_Square_CN,'C','N',2,1,CXXBLACS_N,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PGEMM_1x2_Square_CN,'C','N',1,2,CXXBLACS_N,CXXBLACS_N,CXXBLACS_N);

TEST_IMPL(PGEMM_2x2_Square_NT,'N','T',2,2,CXXBLACS_N,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PGEMM_2x1_Square_NT,'N','T',2,1,CXXBLACS_N,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PGEMM_1x2_Square_NT,'N','T',1,2,CXXBLACS_N,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PGEMM_2x2_Square_TT,'T','T',2,2,CXXBLACS_N,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PGEMM_2x1_Square_TT,'T','T',2,1,CXXBLACS_N,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PGEMM_1x2_Square_TT,'T','T',1,2,CXXBLACS_N,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PGEMM_2x2_Square_CT,'C','T',2,2,CXXBLACS_N,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PGEMM_2x1_Square_CT,'C','T',2,1,CXXBLACS_N,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PGEMM_1x2_Square_CT,'C','T',1,2,CXXBLACS_N,CXXBLACS_N,CXXBLACS_N);

TEST_IMPL(PGEMM_2x2_Square_NC,'N','C',2,2,CXXBLACS_N,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PGEMM_2x1_Square_NC,'N','C',2,1,CXXBLACS_N,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PGEMM_1x2_Square_NC,'N','C',1,2,CXXBLACS_N,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PGEMM_2x2_Square_TC,'T','C',2,2,CXXBLACS_N,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PGEMM_2x1_Square_TC,'T','C',2,1,CXXBLACS_N,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PGEMM_1x2_Square_TC,'T','C',1,2,CXXBLACS_N,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PGEMM_2x2_Square_CC,'C','C',2,2,CXXBLACS_N,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PGEMM_2x1_Square_CC,'C','C',2,1,CXXBLACS_N,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PGEMM_1x2_Square_CC,'C','C',1,2,CXXBLACS_N,CXXBLACS_N,CXXBLACS_N);


TEST_IMPL(PGEMM_2x2_Rect_NN,'N','N',2,2,CXXBLACS_M,CXXBLACS_N,CXXBLACS_K);
TEST_IMPL(PGEMM_2x1_Rect_NN,'N','N',2,1,CXXBLACS_M,CXXBLACS_N,CXXBLACS_K);
TEST_IMPL(PGEMM_1x2_Rect_NN,'N','N',1,2,CXXBLACS_M,CXXBLACS_N,CXXBLACS_K);


