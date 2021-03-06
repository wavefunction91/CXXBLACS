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


#include <mpi.h>
#include <gtest/gtest.h>

int main(int argc, char **argv) {

  MPI_Init(NULL,NULL);

  // Only get print from root
  ::testing::TestEventListeners& listeners =
      ::testing::UnitTest::GetInstance()->listeners();

  int rank; MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  if( rank != 0 ) {
      delete listeners.Release(listeners.default_result_printer());
  }

  // Init GT and run tests
  ::testing::InitGoogleTest(&argc, argv);
  auto gt_return = RUN_ALL_TESTS();

  MPI_Finalize();

  return gt_return; // return GT result

}
