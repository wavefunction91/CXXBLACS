#
# A simple C++ Wrapper for BLACS along with minimal extra functionality to 
# aid the the high-level development of distributed memory linear algebra.
# Copyright (C) 2016-2018 David Williams-Young
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#
cmake_minimum_required(VERSION 3.0)
project(CXXBLACS CXX)
set(CXXBLACS_VERSION_MAJOR 0)
set(CXXBLACS_VERSION_MINOR 1)
set(CXXBLACS_VERSION_PATCH 0)



# C++11
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 11)
if(${CMAKE_CXX_COMPILER} MATCHES "icpc.*$")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -mkl=parallel")
endif()



# MPI
find_package(MPI REQUIRED)
if( NOT ${MPI_CXX_FOUND} )
  message( FATAL_ERROR "No Suitable C++ MPI Compiler Found" )
endif()


# BLACS


# Boost
set(Boost_USE_STATIC_LIBS    ON)
set(Boost_USE_MULTITHREADED  ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost REQUIRED unit_test_framework)


## Setup compile env

# Flags
if( MPI_CXX_COMPILE_FLAGS )
  add_definitions( ${MPI_CXX_COMPILE_FLAGS} )
endif()

# Includes
include_directories( ${PROJECT_SOURCE_DIR}   ) # Project source
include_directories( ${MPI_CXX_INCLUDE_PATH} ) # MPI Headers  
include_directories( ${Boost_INCLUDE_DIRS}   ) # Boost Headers

# Link
set( CXXBLACS_LINKER ${MPI_CXX_LIBRARIES}      ) # MPI Libs
list(APPEND CXXBLACS_LINKER ${Boost_LIBRARIES} ) # Boost Libs
if( CXXBLACS_BLACS_LIBRARIES )
  list(APPEND CXXBLACS_LINKER ${CXXBLACS_BLACS_LIBRARIES} ) # BLACS Libs
else()
  message( FATAL_ERROR 
    "User Must Specify A Valid BLACS Library -DCXXBLACS_BLACS_LIBRARIES=...")
endif()

enable_testing()
add_subdirectory( tests )
