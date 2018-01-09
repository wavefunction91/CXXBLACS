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

#ifndef __INCLUDED_TESTS_UT_HPP__
#define __INCLUDED_TESTS_UT_HPP__


/**
 *  This header contains the configuration information common to
 *  all CXXBLACS UTs through the Boost.Test framework. Should be
 *  included in all UTs
 */

// BOOST_NO_MAIN has strange behaviour if BOOST_TEST_MODULE is defined.
// This is already taken care of in ut.cxx
#ifdef BOOST_TEST_MODULE
  #undef BOOST_TEST_MODULE
#endif



#define BOOST_NO_MAIN
#include <boost/test/unit_test.hpp>

#include <cxxblacs.hpp>


#endif
