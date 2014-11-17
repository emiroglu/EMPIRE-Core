/*  Copyright &copy; 2013, TU Muenchen, Chair of Structural Analysis,
 *  Stefan Sicklinger, Tianyang Wang, Munich
 *
 *  All rights reserved.
 *
 *  This file is part of EMPIRE.
 *
 *  EMPIRE is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  EMPIRE is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with EMPIRE.  If not, see http://www.gnu.org/licenses/.
 */
#include "AuxiliaryParameters.h"
#include <limits>
#define GIT_SHA1 "ee80451f92adffe8e0f54e38f0dccaaaf0626aed" //Cmake will set this value
#define GIT_TAG  "v1.0.0-162-gee80451"  //Cmake will set this value

namespace EMPIRE {

const int AuxiliaryParameters::mklSetNumThreads = 1;
const int AuxiliaryParameters::mapperSetNumThreads= 1;
const double AuxiliaryParameters::machineEpsilon= std::numeric_limits<double>::epsilon();
const std::string AuxiliaryParameters::gitSHA1(GIT_SHA1);
const std::string AuxiliaryParameters::gitTAG(GIT_TAG);
} /* namespace EMPIRE */
