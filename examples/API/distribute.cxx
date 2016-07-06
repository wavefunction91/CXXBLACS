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

  // Compute local scalar product
  int Stride = 1;
  int NLOC = NLocR*NLocC;
  double LDDOT = ddot_(&NLOC,ALoc.data(),&Stride,ALoc.data(),&Stride);

  // Output local scalar products in a Ring
  BlacsGrid::RingExecute([&] { 
    std::cout << "Process " << Grid.iProc() << " has local scalar product "
	     <<  "LDDOT = " << LDDOT << std::endl;
  });

  // Communicate and sum the local doc product and send to BLACS coordinate
  // (0,0)
  BlacsGrid::GSUM2D(Grid.iContxt(),"All","1-tree",1,1,LDDOT,1,0,0);

  // Ouput final result on root process
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
