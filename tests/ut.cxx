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

/**
 *
 *  This file acts as a generic template for CXXBLACS UTs through
 *  the Boost.Test framework. When compiling, one must add a definition
 *  of BOOST_TEST_MODULE in the compiler invocation
 *
 *  CXX -DBOOST_TEST_MODULE=ModName 
 *
 */


#include <mpi.h>

#include <boost/test/impl/unit_test_main.ipp>
#include <boost/test/impl/framework.ipp>
#include <boost/test/unit_test.hpp>

using namespace boost::unit_test;


/**
 *  \brief A global fixture for CXXBLACS UTs through the
 *  Boost.Test framework.
 *
 *  Ensures that MPI_Init() and MPI_Finalize()
 *  only get called once throughout the test module
 */
struct CXXBLACSConfig {

  CXXBLACSConfig() { MPI_Init(NULL,NULL); }
  ~CXXBLACSConfig() { MPI_Finalize(); }

};

BOOST_GLOBAL_FIXTURE( CXXBLACSConfig );
