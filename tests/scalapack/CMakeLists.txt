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

add_executable( scalapack_test ../ut.cxx pgemm.cxx eig.cxx )

target_compile_definitions(scalapack_test PUBLIC 
  BOOST_TEST_MODULE=SCALAPACK)
target_include_directories(scalapack_test PUBLIC 
  ${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/tests)
target_link_libraries( scalapack_test PUBLIC ${CXXBLACS_LINKER} )



add_test( NAME PGEMM_SQP COMMAND ${MPIEXEC} -np 4 "./scalapack_test" "--run_test=PGEMM" )
add_test( NAME PGEMM_RTP COMMAND ${MPIEXEC} -np 2 "./scalapack_test" "--run_test=PGEMM" )
add_test( NAME PGEMM_SER COMMAND ${MPIEXEC} -np 1 "./scalapack_test" "--run_test=PGEMM" )

add_test( NAME PSYEV_SQP COMMAND ${MPIEXEC} -np 4 "./scalapack_test" "--run_test=PSYEV" )
add_test( NAME PSYEV_RTP COMMAND ${MPIEXEC} -np 2 "./scalapack_test" "--run_test=PSYEV" )
add_test( NAME PSYEV_SER COMMAND ${MPIEXEC} -np 1 "./scalapack_test" "--run_test=PSYEV" )

add_test( NAME PHEEV_SQP COMMAND ${MPIEXEC} -np 4 "./scalapack_test" "--run_test=PHEEV" )
add_test( NAME PHEEV_RTP COMMAND ${MPIEXEC} -np 2 "./scalapack_test" "--run_test=PHEEV" )
add_test( NAME PHEEV_SER COMMAND ${MPIEXEC} -np 1 "./scalapack_test" "--run_test=PHEEV" )
