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
project(CXXBLACS VERSION 0.3.0 LANGUAGES C CXX)

# Append local cmake directory to find CMAKE Modules
if( CMAKE_MODULE_PATH )
  list( APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
else()
  set( CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
endif()
list( APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules")

# MPI / ScaLAPACK / BLACS
find_package( ScaLAPACK REQUIRED )


add_library( cxxblacs INTERFACE )
add_library( CXXBLACS::cxxblacs ALIAS cxxblacs )

target_compile_features( cxxblacs INTERFACE cxx_std_11 )



target_link_libraries( cxxblacs INTERFACE ScaLAPACK::scalapack )

target_include_directories( cxxblacs 
  INTERFACE 
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include> 
    $<INSTALL_INTERFACE:include>
)


# INSTALL rules
include( GNUInstallDirs )
set( INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/CXXBLACS )

# Targets
install(TARGETS cxxblacs
    EXPORT cxxblacs-targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

set_target_properties( cxxblacs PROPERTIES EXPORT_NAME CXXBLACS )

# Install Headers
install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} )

# Export target to script
install(EXPORT cxxblacs-targets
  FILE
    CXXBLACSTargets.cmake
  NAMESPACE
    CXXBLACS::
  DESTINATION
    ${INSTALL_CONFIGDIR}
)


# Create a ConfigVersion.cmake file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/CXXBLACSConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

# Setup CXXBLACSConfig.cmake
configure_package_config_file(${PROJECT_SOURCE_DIR}/cmake/CXXBLACSConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/CXXBLACSConfig.cmake
    INSTALL_DESTINATION ${INSTALL_CONFIGDIR}
)

#Install the config, configversion and custom find modules
install(DIRECTORY 
  ${PROJECT_SOURCE_DIR}/cmake/ 
  DESTINATION ${INSTALL_CONFIGDIR} 
  FILES_MATCHING PATTERN "*.cmake"
)


install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/CXXBLACSConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/CXXBLACSConfigVersion.cmake
    DESTINATION ${INSTALL_CONFIGDIR}
)









option( ENABLE_CXXBLACS_TESTS "Enable Testing" OFF )
if( ENABLE_CXXBLACS_TESTS )
  enable_testing()
  add_subdirectory( tests )
endif()
