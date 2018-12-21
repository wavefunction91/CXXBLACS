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
void ptrmm_test(char SIDE, char UPLO, char TRANSA, char DIAG,  Field ALPHA,
    CB_INT M, CB_INT N) {

  // Synchronize processes
  MPI_Barrier(MPI_COMM_WORLD);

  BlacsGrid grid(MPI_COMM_WORLD,MB,NB);


  // Make sure that all of the ALPHA's and BETA's are the same
  grid.Broadcast("All","I",1,1,&ALPHA,1,0,0);

  std::vector<Field> A, ALoc, B, BLoc, TrueAns;

  CB_INT AM = (SIDE == 'L') ? M : N;


  // Get local dims and allocate local buffers
  CB_INT ANLocR, ANLocC;
  CB_INT BNLocR, BNLocC;
  std::tie(ANLocR,ANLocC) = grid.getLocalDims(AM,AM);
  std::tie(BNLocR,BNLocC) = grid.getLocalDims(M,N);

  ALoc.resize(ANLocR*ANLocC);
  BLoc.resize(BNLocR*BNLocC);





  // Generate random matricies on root process and get serial product
  RootExecute(MPI_COMM_WORLD,[&](){

    A.resize(AM*AM);
    B.resize(M*N);
    TrueAns.resize(M*N);

    for(auto &X : A ) X = generate<Field>();
    for(auto &X : B ) X = generate<Field>();

    std::copy(B.begin(),B.end(),TrueAns.begin());

    TRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,&A[0],AM,&TrueAns[0],M);

  });




  // Scatter the matrices onto the process grid
  grid.Scatter(AM,AM,A.data(),AM,ALoc.data(),ANLocR,0,0);
  grid.Scatter(M ,N ,B.data(),M ,BLoc.data(),BNLocR,0,0);




  // Get ScaLAPACK descriptiors
  auto DESCA = grid.descInit(AM,AM,0,0,ANLocR);
  auto DESCB = grid.descInit(M ,N ,0,0,BNLocR);


  // Get Product
  CXXBLACS::PTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,
    ALoc.data(),1,1,DESCA,
    BLoc.data(),1,1,DESCB);


  // Gather product to root process
  grid.Gather(M,N,B.data(),M,BLoc.data(),BNLocR,0,0);

  // Test result on root process
  RootExecute(MPI_COMM_WORLD,[&](){

    std::vector<RealType> DIFF(TrueAns.size(),0.);
    for(auto k = 0; k < DIFF.size(); k++) 
      DIFF[k] = std::abs(TrueAns[k] - B[k]);

    RealType maxDiff = *std::max_element(DIFF.begin(),DIFF.end());

    std::cout << "MAX DIFF " << maxDiff << std::endl;


    EXPECT_NEAR( maxDiff, 0., 1e-10 ) << 
      "MAX DIFF " << maxDiff
      << " " << std::numeric_limits<RealType>::epsilon() 
      << " " << std::numeric_limits<RealType>::epsilon() * M*N*M ;

  });

  NotRootExecute(MPI_COMM_WORLD,[&](){ EXPECT_TRUE(true); });

  // Synchronize processes
  MPI_Barrier(MPI_COMM_WORLD);

};

#define TEST_IMPL_F(NAME,F,RF,SIDE,UPLO,TA,DIAG,MB,NB,M,N)\
  TEST(PTRMM,NAME) { \
    ptrmm_test<F,MB,NB,RF>(SIDE,UPLO,TA,DIAG,generate<F>(),M,N); \
  };

#define TEST_IMPL(NAME,SIDE,UPLO,TA,DIAG,MB,NB,M,N)\
  TEST_IMPL_F(NAME##_Double, double,              double,SIDE,UPLO,TA,DIAG,MB,NB,M,N)\
  TEST_IMPL_F(NAME##_CDouble,std::complex<double>,double,SIDE,UPLO,TA,DIAG,MB,NB,M,N)\
//FIXME: Need to determine proper error criteria for single precision
//TEST_IMPL_F(NAME##_float, float,              float,SIDE,UPLO,TA,DIAG,MB,NB,M,N)\
//TEST_IMPL_F(NAME##_Cfloat,std::complex<float>,float,SIDE,UPLO,TATA,TB,MB,NB,M,N)\




TEST_IMPL(PTRMM_2x2_Square_LUNN,'L','U','N','N',2,2,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_2x1_Square_LUNN,'L','U','N','N',2,1,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_1x2_Square_LUNN,'L','U','N','N',1,2,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_2x2_Square_LUNU,'L','U','N','U',2,2,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_2x1_Square_LUNU,'L','U','N','U',2,1,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_1x2_Square_LUNU,'L','U','N','U',1,2,CXXBLACS_N,CXXBLACS_N);

TEST_IMPL(PTRMM_2x2_Square_LUTN,'L','U','T','N',2,2,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_2x1_Square_LUTN,'L','U','T','N',2,1,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_1x2_Square_LUTN,'L','U','T','N',1,2,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_2x2_Square_LUTU,'L','U','T','U',2,2,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_2x1_Square_LUTU,'L','U','T','U',2,1,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_1x2_Square_LUTU,'L','U','T','U',1,2,CXXBLACS_N,CXXBLACS_N);

TEST_IMPL(PTRMM_2x2_Square_LUCN,'L','U','C','N',2,2,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_2x1_Square_LUCN,'L','U','C','N',2,1,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_1x2_Square_LUCN,'L','U','C','N',1,2,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_2x2_Square_LUCU,'L','U','C','U',2,2,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_2x1_Square_LUCU,'L','U','C','U',2,1,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_1x2_Square_LUCU,'L','U','C','U',1,2,CXXBLACS_N,CXXBLACS_N);

TEST_IMPL(PTRMM_2x2_Square_LLNN,'L','L','N','N',2,2,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_2x1_Square_LLNN,'L','L','N','N',2,1,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_1x2_Square_LLNN,'L','L','N','N',1,2,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_2x2_Square_LLNU,'L','L','N','U',2,2,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_2x1_Square_LLNU,'L','L','N','U',2,1,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_1x2_Square_LLNU,'L','L','N','U',1,2,CXXBLACS_N,CXXBLACS_N);

TEST_IMPL(PTRMM_2x2_Square_LLTN,'L','L','T','N',2,2,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_2x1_Square_LLTN,'L','L','T','N',2,1,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_1x2_Square_LLTN,'L','L','T','N',1,2,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_2x2_Square_LLTU,'L','L','T','U',2,2,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_2x1_Square_LLTU,'L','L','T','U',2,1,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_1x2_Square_LLTU,'L','L','T','U',1,2,CXXBLACS_N,CXXBLACS_N);

TEST_IMPL(PTRMM_2x2_Square_LLCN,'L','L','C','N',2,2,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_2x1_Square_LLCN,'L','L','C','N',2,1,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_1x2_Square_LLCN,'L','L','C','N',1,2,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_2x2_Square_LLCU,'L','L','C','U',2,2,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_2x1_Square_LLCU,'L','L','C','U',2,1,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_1x2_Square_LLCU,'L','L','C','U',1,2,CXXBLACS_N,CXXBLACS_N);

TEST_IMPL(PTRMM_2x2_Square_RUNN,'R','U','N','N',2,2,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_2x1_Square_RUNN,'R','U','N','N',2,1,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_1x2_Square_RUNN,'R','U','N','N',1,2,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_2x2_Square_RUNU,'R','U','N','U',2,2,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_2x1_Square_RUNU,'R','U','N','U',2,1,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_1x2_Square_RUNU,'R','U','N','U',1,2,CXXBLACS_N,CXXBLACS_N);

TEST_IMPL(PTRMM_2x2_Square_RUTN,'R','U','T','N',2,2,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_2x1_Square_RUTN,'R','U','T','N',2,1,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_1x2_Square_RUTN,'R','U','T','N',1,2,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_2x2_Square_RUTU,'R','U','T','U',2,2,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_2x1_Square_RUTU,'R','U','T','U',2,1,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_1x2_Square_RUTU,'R','U','T','U',1,2,CXXBLACS_N,CXXBLACS_N);

TEST_IMPL(PTRMM_2x2_Square_RUCN,'R','U','C','N',2,2,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_2x1_Square_RUCN,'R','U','C','N',2,1,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_1x2_Square_RUCN,'R','U','C','N',1,2,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_2x2_Square_RUCU,'R','U','C','U',2,2,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_2x1_Square_RUCU,'R','U','C','U',2,1,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_1x2_Square_RUCU,'R','U','C','U',1,2,CXXBLACS_N,CXXBLACS_N);

TEST_IMPL(PTRMM_2x2_Square_RLNN,'R','L','N','N',2,2,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_2x1_Square_RLNN,'R','L','N','N',2,1,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_1x2_Square_RLNN,'R','L','N','N',1,2,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_2x2_Square_RLNU,'R','L','N','U',2,2,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_2x1_Square_RLNU,'R','L','N','U',2,1,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_1x2_Square_RLNU,'R','L','N','U',1,2,CXXBLACS_N,CXXBLACS_N);

TEST_IMPL(PTRMM_2x2_Square_RLTN,'R','L','T','N',2,2,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_2x1_Square_RLTN,'R','L','T','N',2,1,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_1x2_Square_RLTN,'R','L','T','N',1,2,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_2x2_Square_RLTU,'R','L','T','U',2,2,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_2x1_Square_RLTU,'R','L','T','U',2,1,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_1x2_Square_RLTU,'R','L','T','U',1,2,CXXBLACS_N,CXXBLACS_N);

TEST_IMPL(PTRMM_2x2_Square_RLCN,'R','L','C','N',2,2,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_2x1_Square_RLCN,'R','L','C','N',2,1,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_1x2_Square_RLCN,'R','L','C','N',1,2,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_2x2_Square_RLCU,'R','L','C','U',2,2,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_2x1_Square_RLCU,'R','L','C','U',2,1,CXXBLACS_N,CXXBLACS_N);
TEST_IMPL(PTRMM_1x2_Square_RLCU,'R','L','C','U',1,2,CXXBLACS_N,CXXBLACS_N);




TEST_IMPL(PTRMM_2x2_Rect_LUNN,'L','U','N','N',2,2,CXXBLACS_M,CXXBLACS_N);
TEST_IMPL(PTRMM_2x1_Rect_LUNN,'L','U','N','N',2,1,CXXBLACS_M,CXXBLACS_N);
TEST_IMPL(PTRMM_1x2_Rect_LUNN,'L','U','N','N',1,2,CXXBLACS_M,CXXBLACS_N);
TEST_IMPL(PTRMM_2x2_Rect_LUNU,'L','U','N','U',2,2,CXXBLACS_M,CXXBLACS_N);
TEST_IMPL(PTRMM_2x1_Rect_LUNU,'L','U','N','U',2,1,CXXBLACS_M,CXXBLACS_N);
TEST_IMPL(PTRMM_1x2_Rect_LUNU,'L','U','N','U',1,2,CXXBLACS_M,CXXBLACS_N);


