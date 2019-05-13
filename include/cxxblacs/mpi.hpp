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
#ifndef __INCLUDED_CXXBLACS_MPI_HPP__
#define __INCLUDED_CXXBLACS_MPI_HPP__

#include <cxxblacs/config.hpp>

#define CXXBLACS_MPI_ROOT 0
#define CXXBLACS_MPI_DEFAULT_TAG 0


namespace CXXBLACS {

  template <typename Func>
  inline void RootExecute(const MPI_Comm c, const Func& op) {

    int iProc; MPI_Comm_rank(c, &iProc);

    if(iProc == CXXBLACS_MPI_ROOT) op();

  }

  template <typename Func>
  inline void NotRootExecute(const MPI_Comm c, const Func& op) {

    int iProc; MPI_Comm_rank(c, &iProc);

    if(iProc != CXXBLACS_MPI_ROOT) op();

  }


  template <typename Func, typename... Args>
  inline void RingExecute(const MPI_Comm c, const Func& op, Args... args) {

    MPI_Barrier(c);

    int iProc, nProc;
    MPI_Comm_rank(c,&iProc);
    MPI_Comm_size(c,&nProc);

    if( nProc == 1 ) op(args...);
    else {

      int iToken = 0;
      int source = (iProc == 0) ? nProc - 1 : iProc - 1;
      int dest = (iProc+1) % nProc;

      auto run = [&]() {
        MPI_Recv(&iToken,1,MPI_INT,source,CXXBLACS_MPI_DEFAULT_TAG,c,
          MPI_STATUS_IGNORE);
        op(args...);
      };

      NotRootExecute(c, run);

      MPI_Send(&iToken,1,MPI_INT,dest,CXXBLACS_MPI_DEFAULT_TAG,c);
      MPI_Barrier(c);

      RootExecute(c, run);

    }
  }

};

#endif
