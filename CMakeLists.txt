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
cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(CXXBLACS VERSION 0.2.0 LANGUAGES C CXX)

# Append local cmake directory to find CMAKE Modules
if( CMAKE_MODULE_PATH )
  list( APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
else()
  set( CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
endif()
list( APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules")

# MPI / ScaLAPACK / BLACS
find_package( ScaLAPACK REQUIRED )


add_library( CXXBLACS::CXXBLACS INTERFACE IMPORTED )
target_compile_features( CXXBLACS::CXXBLACS INTERFACE cxx_std_11 )
target_link_libraries( CXXBLACS::CXXBLACS INTERFACE ScaLAPACK::scalapack )
target_include_directories( CXXBLACS::CXXBLACS INTERFACE ${PROJECT_SOURCE_DIR} )


enable_testing()
add_subdirectory( tests )
