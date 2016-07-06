/**
 * \brief Print a greeting from each MPI process
 */
void printMPIInfo(){
  std::stringstream ss;
  ss << "Greetings from process " << this->iProc_ << " out of "
     << this->nProc_ << std::endl;
  std::cout << ss.str();
  BlacsGrid::Barrier(this->IContxt_,"A");
}

/**
 * \brief Print a greeting from each BLACS coordinate
 */
void printBLACSInfo(){
  std::stringstream ss;
  ss << "Greetings from process coordinate (" << this->iProcRow_ << "," 
      << this->iProcCol_ << ") out of (" << this->nProcRow_ << "," 
      << this->nProcCol_ << ")" << std::endl;
  std::cout << ss.str();
  BlacsGrid::Barrier(this->IContxt_,"A");
}

