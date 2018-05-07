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
#include <math.h>
#include "cppunit/TestFixture.h"
#include "cppunit/TestAssert.h"
#include "cppunit/extensions/HelperMacros.h"


#include "AuxiliaryParameters.h"
#include "GeometryMath.h"

namespace EMPIRE {
using namespace std;
using namespace MathLibrary;

/********//**
 * \brief This class manages tests the MathLibrary of EMPIRE
 **************************************************************************************************/
class TestGeometryMath: public CppUnit::TestFixture {
private:
    int numVertices;
    std::vector<double> polygon;
    double pointOut[2];
    double pointIn[2];
public:
    /***********************************************************************************************
     * \brief Set up some test vectors and matrices
     * \author Andreas Apostolatos
     ***********/
    void setUp() {
        numVertices = 15;
        static double vertices[30] = {0.5, 0.1, 0, 0.1, 0, 0, 0.5, 0, 0.483654, 0.00274729,
                           0.468042, 0.0115466, 0.464645, 0.0146447, 0.456387, 0.0255479,
                           0.450655, 0.0419322, 0.45, 0.05, 0.450655, 0.0580678, 0.456387, 0.0744521,
                           0.464645, 0.0853553, 0.468042, 0.0884534, 0.483654, 0.0972527};
        for (int i = 0; i < 2*numVertices; i++) {
            polygon.push_back(vertices[i]);
        }
        pointOut[0] = (0.45 + 0.5)/2.0;
        pointOut[1] = 0.075;
        pointIn[0] = (0.5 + 0)/2.0;
        pointIn[1] = 0.025;
    }

    /***********************************************************************************************
     * \brief Delete test vectors and matrices
     * \author Andreas Apostolatos
     ***********/
    void tearDown() {
    }

    /***********************************************************************************************
     * \brief Test dense dot product function
     * \author Andreas Apostolatos
     ***********/
    void testFindIfPointIsInside2DPolygon() {

        bool isInsideNodeInside = findIfPointIsInside2DPolygon(numVertices, polygon, pointIn);
        CPPUNIT_ASSERT(isInsideNodeInside);
        bool isInsideNodeOutside = findIfPointIsInside2DPolygon(numVertices, polygon, pointOut);
        CPPUNIT_ASSERT(!isInsideNodeOutside);
    }

    CPPUNIT_TEST_SUITE( TestGeometryMath );
    CPPUNIT_TEST(testFindIfPointIsInside2DPolygon);
    CPPUNIT_TEST_SUITE_END();
};

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION( EMPIRE::TestGeometryMath);
