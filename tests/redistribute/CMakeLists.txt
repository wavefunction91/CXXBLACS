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

add_executable( redistribute_test ../ut.cxx redistribute.cxx )

target_compile_definitions(redistribute_test PUBLIC 
  BOOST_TEST_MODULE=REDISTRIBUTE)
target_include_directories(redistribute_test PUBLIC 
  ${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/tests)
target_link_libraries( redistribute_test PUBLIC ${CXXBLACS_LINKER} )



add_test( NAME REDIST_SQP COMMAND ${MPIEXEC} -np 4 "./redistribute_test" )
add_test( NAME REDIST_RTP COMMAND ${MPIEXEC} -np 2 "./redistribute_test" )
add_test( NAME REDIST_SER COMMAND ${MPIEXEC} -np 1 "./redistribute_test" )

