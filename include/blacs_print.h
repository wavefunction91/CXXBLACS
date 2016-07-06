/*
 *  A simple C++ Wrapper for BLACS along with minimal extra functionality to 
 *  aid the the high-level development of distributed memory linear algebra.
 *  Copyright (C) 2016 David Williams-Young

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

