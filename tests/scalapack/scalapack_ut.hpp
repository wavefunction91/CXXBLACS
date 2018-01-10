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

#include <ut.hpp>
#include <cxxblacs.hpp>

#include <algorithm>
#include <random>

using namespace CXXBLACS;

constexpr CB_INT CXXBLACS_M = 200;
constexpr CB_INT CXXBLACS_N = 150;
constexpr CB_INT CXXBLACS_K = 7;

constexpr CB_INT CXXBLACS_NRHS = 15;

static std::default_random_engine gen;

namespace CXXBLACS {
  template <typename Field>
  inline Field generate();
  
  template <>
  inline double generate(){ 
    std::uniform_real_distribution<double> dis(-10,10);
    return dis(gen); 
  }
  
  template <>
  inline std::complex<double> generate(){ 
    std::uniform_real_distribution<double> dis(-10,10);
    return std::complex<double>(dis(gen),dis(gen)); 
  }

  template <>
  inline float generate(){ 
    std::uniform_real_distribution<float> dis(-10,10);
    return dis(gen); 
  }
  
  template <>
  inline std::complex<float> generate(){ 
    std::uniform_real_distribution<float> dis(-10,10);
    return std::complex<float>(dis(gen),dis(gen)); 
  }
}




