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
#include "cppunit/TestFixture.h"
#include "cppunit/TestAssert.h"
#include "cppunit/extensions/HelperMacros.h"

#include "CurveSurfaceCorotate2DMapper.h"
#include "KinematicMotion.h"
#include <iostream>
#include <string>
#include <math.h>

using namespace std;

namespace EMPIRE {
/********//**
 * \brief Test the class CurveSurfaceCorotate2DMapper.
 ***********/
class TestCurveSurfaceCorotate2DMapper: public CppUnit::TestFixture {
private:

public:
    void setUp() {
    }
    void tearDown() {
    }


    /***********************************************************************************************
     * \brief Test consistent and conservative mapping
     ***********/
    void testBending() {
        // define curve and surface in Q
        int curveNumNodes = 3;
        int curveNumElements = 2;
        double curveNodeCoors[] = { 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 1.0, 0.0, 0.0 };
        int curveNodeIDs[] = { 1, 2, 3 };
        int curveElems[] = { 1, 2, 2, 3 };
        int surfaceNumNodes = 5;
        double surfaceNodeCoors[] = { 0.0, -1.0, 0.0, 0.25, -1.0, 0.0, 0.5, -1.0, 0.0,
                0.75, -1.0, 0.0, 1.0, -1.0, 0.0 };

        int surfaceNumSections = 5;
        int surfaceNumRootSectionNodes = 1;
        int surfaceNumNormalSectionNodes = 1;
        int surfaceNumTipSectionNodes = 1;

        KinematicMotion KM_O_Q;
        double axis[] = { 1.0, 1.0, 1.0 };
        //KM_O_Q.addRotation(axis, false, M_PI / 2.0);
        double translation_O_Q[] = { 2.0, 3.0, 4.0 };
        //KM_O_Q.addTranslation(translation_O_Q);

        // transform the coordinates to O
        for (int i = 0; i < curveNumNodes; i++) {
            KM_O_Q.move(&curveNodeCoors[i * 3]);
        }
        for (int i = 0; i < surfaceNumNodes; i++) {
            KM_O_Q.move(&surfaceNodeCoors[i * 3]);
        }

        // construct the mapper
        CurveSurfaceCorotate2DMapper *mapper = new CurveSurfaceCorotate2DMapper(curveNumNodes,
                curveNumElements, curveNodeCoors, curveNodeIDs, curveElems, surfaceNumNodes,
                surfaceNodeCoors, surfaceNumSections, surfaceNumRootSectionNodes,
                surfaceNumNormalSectionNodes, surfaceNumTipSectionNodes, KM_O_Q.getRotationMatrix(),
                KM_O_Q.getTranslationVector());
        { // test consistent mapping
          // define deformation of curve/beam
            double curveDispRot[18];
            // curve node 1
            curveDispRot[0] = 0.0;
            curveDispRot[1] = 0.0;
            curveDispRot[2] = 0.0;
            curveDispRot[3] = 0.0;
            curveDispRot[4] = 0.0;
            curveDispRot[5] = 0.0;
            // curve node 2
            curveDispRot[6] = 1.0 - 0.5;
            curveDispRot[7] = 1.0 - 0.0;
            curveDispRot[8] = 0.0;
            curveDispRot[9] = 0.0;
            curveDispRot[10] = 0.0;
            curveDispRot[11] = 1.5707963;
            // curve node 2
            curveDispRot[12] = 0.0 - 1.0;
            curveDispRot[13] = 2.0 - 0.0;
            curveDispRot[14] = 0.0;
            curveDispRot[15] = 0.0;
            curveDispRot[16] = 0.0;
            curveDispRot[17] = 3.1415926;

            // rotate the DOFs form Q to O
            KinematicMotion ROT_O_Q;
            ROT_O_Q.addRotation(KM_O_Q.getRotationMatrix());
            for (int i = 0; i < curveNumNodes * 2; i++) {
                //ROT_O_Q.move(&curveDispRot[i * 3]);
            }

            // consistent mapping
            double surfaceDisp[15];
            mapper->consistentMapping(curveDispRot, surfaceDisp);

            /*cout << "surfaceDisp" << endl;
             for (int i = 0; i < 5; i++) {
             for (int j = 0; j < 3; j++)
             cout << " " << surfaceDisp[i * 3 + j];
             cout << endl;
             }
             cout << "surfaceNodeCoors" << endl;
             for (int i = 0; i < 5; i++) {
             for (int j = 0; j < 3; j++)
             cout << " " << surfaceNodeCoors[i * 3 + j] + surfaceDisp[i * 3 + j];
             cout << endl;
             }*/
            for (int i = 0; i < 5; i++) {
                for (int j = 0; j < 3; j++)
                    surfaceNodeCoors[i * 3 + j] += surfaceDisp[i * 3 + j];
            }
            const double TOL = 2E-2; // cannot be very exact to the circle due to cubic interpolation
            CPPUNIT_ASSERT(fabs(surfaceNodeCoors[0 * 3 + 0] - 0.0) < TOL);
            CPPUNIT_ASSERT(fabs(surfaceNodeCoors[0 * 3 + 1] - (-1.0)) < TOL);
            CPPUNIT_ASSERT(fabs(surfaceNodeCoors[0 * 3 + 2] - 0.0) < TOL);
            CPPUNIT_ASSERT(fabs(surfaceNodeCoors[1 * 3 + 0] - sqrt(2.0)) < TOL);
            CPPUNIT_ASSERT(fabs(surfaceNodeCoors[1 * 3 + 1] - (1.0 - sqrt(2.0))) < TOL);
            CPPUNIT_ASSERT(fabs(surfaceNodeCoors[1 * 3 + 2] - 0.0) < TOL);
            CPPUNIT_ASSERT(fabs(surfaceNodeCoors[2 * 3 + 0] - 2.0) < TOL);
            CPPUNIT_ASSERT(fabs(surfaceNodeCoors[2 * 3 + 1] - 1.0) < TOL);
            CPPUNIT_ASSERT(fabs(surfaceNodeCoors[2 * 3 + 2] - 0.0) < TOL);
            CPPUNIT_ASSERT(fabs(surfaceNodeCoors[3 * 3 + 0] - sqrt(2.0)) < TOL);
            CPPUNIT_ASSERT(fabs(surfaceNodeCoors[3 * 3 + 1] - (1.0 + sqrt(2.0))) < TOL);
            CPPUNIT_ASSERT(fabs(surfaceNodeCoors[3 * 3 + 2] - 0.0) < TOL);
            CPPUNIT_ASSERT(fabs(surfaceNodeCoors[4 * 3 + 0] - 0.0) < TOL);
            CPPUNIT_ASSERT(fabs(surfaceNodeCoors[4 * 3 + 1] - 3.0) < TOL);
            CPPUNIT_ASSERT(fabs(surfaceNodeCoors[4 * 3 + 2] - 0.0) < TOL);
        }

        delete mapper;
    }

    CPPUNIT_TEST_SUITE (TestCurveSurfaceCorotate2DMapper);
    CPPUNIT_TEST (testBending);
    CPPUNIT_TEST_SUITE_END();
};

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION (EMPIRE::TestCurveSurfaceCorotate2DMapper);
