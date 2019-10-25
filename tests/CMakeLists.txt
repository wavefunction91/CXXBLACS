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

include( HandleGTest )



add_library( ut_framework INTERFACE IMPORTED )
target_link_libraries( ut_framework INTERFACE cxxblacs )
target_link_libraries( ut_framework INTERFACE CXXBLACS::gtest )
target_include_directories( ut_framework INTERFACE ${PROJECT_SOURCE_DIR}/tests )

add_subdirectory(scatter_gather)
add_subdirectory(redistribute)
add_subdirectory(scalapack)


