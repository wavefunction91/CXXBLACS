/**
 * \brief Execute a command in a ring around the MPI processes
 *
 * Pass a token from 0 -> 1 -> 2 -> ... -> N -> 0 and execute a command
 * in the same order 
 *
 * Example usage:
 * @code
 *   int iProc, nProc;
 *   BlacsGrid::PINFO(iProc,nProc);
 *   BlacsGrid::RingExecute([&] {
 *     std::cout << "Hello from process" << iProc << std::endl;
 *   });
 * @endcode
 */
template<typename T>
static void RingExecute(const T& operation){
  MPI_Barrier(MPI_COMM_WORLD);
  int iProc,nProc;
  BlacsGrid::PINFO(iProc,nProc);
  int iToken = 0;

  if(nProc == 1) operation();
  else {
    BlacsGrid::NotRootExecute([&]{
      MPI_Recv(&iToken,1,MPI_INT,iProc-1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      operation();
    });
 
    MPI_Send(&iToken,1,MPI_INT,(iProc+1)%nProc,0,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
 
    BlacsGrid::RootExecute([&]{
      MPI_Recv(&iToken,1,MPI_INT,nProc-1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      operation();
    });
  }
};

/**
 * \brief Execute command only on root process
 *
 * Example usage:
 * @code
 *   BlacsGrid::RootExecute([&] {
 *     std::cout << "Hello from root process" << std::endl;
 *   });
 * @endcode
 */

template<typename T>
static void RootExecute(const T& operation){
  int iProc,nProc;
  BlacsGrid::PINFO(iProc,nProc);
  if(iProc == 0) operation();
}

/**
 * \brief Execute command on every process but root
 *
 * Example usage:
 * @code
 *   BlacsGrid::NotRootExecute([&] {
 *     std::cout << "Hello from not root process" << std::endl;
 *   });
 * @endcode
 */
template<typename T>
static void NotRootExecute(const T& operation){
  int iProc,nProc;
  BlacsGrid::PINFO(iProc,nProc);
  if(iProc != 0) operation();
}

