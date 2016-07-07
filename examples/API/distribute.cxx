#include "cxxblacs.h"

void LocalIdentity(const int N, double *A);

int main() {
  // Initialize the BLACS Grid (NB = MB = 2)
  BlacsGrid Grid;

  // Leading dimension of distributed Matrix
  const int N = 100;

  // Create Identity on Root process
  std::vector<double> A;
  BlacsGrid::RootExecute([&] {
    A.resize(N*N);
    LocalIdentity(N,A.data());
  });


  // Allocate local buffers of distributed matrix
  int NLocR,NLocC;
  Grid.getLocalDims(N,N,NLocR,NLocC);
  
  
  std::vector<double> ALoc(NLocR*NLocC,0.0);

  // Distribute the matrix to ALoc
  Grid.distribute(N,A,ALoc);

  // Compute dot product in parallel to verify
  double LDDOT;
  Grid.PDOT(N,ALoc.data(),1,ALoc.data(),1,LDDOT);
  BlacsGrid::RootExecute([&] { 
    std::cout << "Root Process has total DDOT = " << LDDOT << std::endl; 
  });

};

void LocalIdentity(const int N, double *A){
  // Build identity on root process
  for(auto J = 0ul; J < N; ++J) 
  for(auto I = 0ul; I < N; ++I) {
    double el = 0.0;
    if( I == J ) el = 1.0;
    A[I + J*N]  = el;
  }
};
