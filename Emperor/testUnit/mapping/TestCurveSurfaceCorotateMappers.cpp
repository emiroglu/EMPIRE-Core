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
#include "CurveSurfaceCorotate3DMapper.h"
#include "KinematicMotion.h"
#include <iostream>
#include <string>
#include <math.h>

using namespace std;

namespace EMPIRE {
/********//**
 * \brief Test the class CurveSurfaceCorotate2DMapper.
 ***********/
class TestCurveSurfaceCorotateMappers: public CppUnit::TestFixture {
private:

public:
    void setUp() {
    }
    void tearDown() {
    }

    /***********************************************************************************************
     * \brief Bend a beam into half circle
     ***********/
    void testBending() {
        // length of the beam is 10.0, the beam is bended and prolonged and twisted
        // counter = 0 for corotate2D; counter = 1 for corotate3D
        for (int counter = 0; counter < 2; counter++) {
            // two parameters to control: curveNumNodes and RADIUS
            // curve mesh in Q
            int curveNumNodes = 21;
            int curveNumElements = curveNumNodes - 1;
            double *curveNodeCoors = new double[curveNumNodes * 3];
            for (int i = 0; i < curveNumNodes; i++) {
                curveNodeCoors[i * 3 + 0] = (double) (i) / (double) (curveNumElements) * 10.0;
                curveNodeCoors[i * 3 + 1] = 0.0;
                curveNodeCoors[i * 3 + 2] = 0.0;
            }
            int *curveNodeIDs = new int[curveNumNodes];
            for (int i = 0; i < curveNumNodes; i++) {
                curveNodeIDs[i] = i + 1;
            }
            int *curveElems = new int[curveNumElements * 2];
            for (int i = 0; i < curveNumElements; i++) {
                curveElems[i * 2 + 0] = i + 1;
                curveElems[i * 2 + 1] = i + 2;
            }
            double RADIUS = 5.0;
            //double RADIUS = 10.0 / M_PI;
            // define deformation of curve/beam in Q
            double *curveDispRot = new double[6 * curveNumNodes];
            double axisX[] = { 1.0, 0.0, 0.0 };
            double axisY[] = { 0.0, 1.0, 0.0 };
            double axisZ[] = { 0.0, 0.0, 1.0 };
            for (int i = 0; i < curveNumNodes; i++) {
                double angle = 0.0 + (double) (i) * (M_PI / curveNumElements);
                curveDispRot[i * 6 + 0] = RADIUS * sin(angle) - curveNodeCoors[i * 3 + 0];
                curveDispRot[i * 6 + 1] = (RADIUS - RADIUS * cos(angle))
                        - curveNodeCoors[i * 3 + 1];
                curveDispRot[i * 6 + 2] = 0.0;

                curveDispRot[i * 6 + 3] = angle / 2.0;
                //curveDispRot[i * 6 + 3] = 0.0;
                curveDispRot[i * 6 + 4] = 0.0;
                curveDispRot[i * 6 + 5] = angle;

                if (counter == 1) {
                    KinematicMotion rot;
                    rot.addRotation(axisX, true, curveDispRot[i * 6 + 3]);
                    rot.addRotation(axisZ, true, curveDispRot[i * 6 + 5]);
                    rot.getRotationVector(&curveDispRot[i * 6 + 3]);
                }
            }

            // surface mesh in Q
            int surfaceNumNodes = 11;
            double *surfaceNodeCoors = new double[surfaceNumNodes * 3];
            for (int i = 0; i < surfaceNumNodes; i++) {
                surfaceNodeCoors[i * 3 + 0] = (double) (i) / (double) (surfaceNumNodes - 1) * 10.0;
                surfaceNodeCoors[i * 3 + 1] = 0.0;
                surfaceNodeCoors[i * 3 + 2] = -1.0;
            }

            // section information
            int numRootSectionNodes = 1;
            int numTipSectionNodes = 1;
            int numNormalSectionNodes = 1;
            int numSections = surfaceNumNodes;

            // compute km
            KinematicMotion *km_Q_O = new KinematicMotion();
            double newX[] = { 0.0, 0.0, -1.0 };
            double newY[] = { 0.0, 1.0, 0.0 };
            double newZ[] = { 1.0, 0.0, 0.0 };
            //km_Q_O->addRotation(newX, newY, newZ, true);
            double translate_Q_O[] = { -1.0, 0.0, 0.0 };
            //km_Q_O->addTranslation(translate_Q_O);

            KinematicMotion *km_O_Q = km_Q_O->newInverse();
            // transform coordinates from Q to O with km_O_Q
            for (int i = 0; i < surfaceNumNodes; i++) {
                km_O_Q->move(&surfaceNodeCoors[i * 3]);
            }
            for (int i = 0; i < curveNumNodes; i++) {
                km_O_Q->move(&curveNodeCoors[i * 3]);
            }
            KinematicMotion rot_O_Q;
            rot_O_Q.addRotation(km_O_Q->getRotationMatrix());
            // transform coordinates from Q to O with rot_O_Q
            for (int i = 0; i < curveNumNodes; i++) {
                rot_O_Q.move(&curveDispRot[i * 6]);
                rot_O_Q.move(&curveDispRot[i * 6 + 3]); // Can be done since O and Q are parallel
            }

            // disp of surface
            double *displacement = new double[surfaceNumNodes * 3];

            if (counter == 0) {
                CurveSurfaceCorotate2DMapper *mapper = new CurveSurfaceCorotate2DMapper(
                        curveNumNodes, curveNumElements, curveNodeCoors, curveNodeIDs, curveElems,
                        surfaceNumNodes, surfaceNodeCoors, numSections, numRootSectionNodes,
                        numNormalSectionNodes, numTipSectionNodes, km_O_Q->getRotationMatrix(),
                        km_O_Q->getTranslationVector());

                mapper->consistentMapping(curveDispRot, displacement);
                delete mapper;
            }
            if (counter == 1) {
                CurveSurfaceCorotate3DMapper *mapper = new CurveSurfaceCorotate3DMapper(
                        curveNumNodes, curveNumElements, curveNodeCoors, curveNodeIDs, curveElems,
                        surfaceNumNodes, surfaceNodeCoors, numSections, numRootSectionNodes,
                        numNormalSectionNodes, numTipSectionNodes, km_O_Q->getRotationMatrix(),
                        km_O_Q->getTranslationVector());

                mapper->consistentMapping(curveDispRot, displacement);
                delete mapper;
            }

            for (int i = 0; i < surfaceNumNodes * 3; i++) {
                surfaceNodeCoors[i] += displacement[i];
            }
            /*cout << "surfaceNodeCoors: " << endl;
            for (int i = 0; i < surfaceNumNodes; i++) {
                cout << "  " << surfaceNodeCoors[i * 3 + 0] << "  " << surfaceNodeCoors[i * 3 + 1]
                        << "  " << surfaceNodeCoors[i * 3 + 2] << endl;
            }*/
            double tipNode[3];
            tipNode[0] = surfaceNodeCoors[(surfaceNumNodes - 1) * 3 + 0];
            tipNode[1] = surfaceNodeCoors[(surfaceNumNodes - 1) * 3 + 1];
            tipNode[2] = surfaceNodeCoors[(surfaceNumNodes - 1) * 3 + 2];

            /*cout << "tipNode: " << "  " << tipNode[0] << "  " << tipNode[1] << "  " << tipNode[2]
                    << endl;*/

            const double EPS = 1E-3; // corotate3D has bigger error due to big torsion
            CPPUNIT_ASSERT(fabs(tipNode[0] - 0.0) < EPS);
            CPPUNIT_ASSERT(fabs(tipNode[1] - (2.0 * RADIUS - 1.0)) < EPS);
            CPPUNIT_ASSERT(fabs(tipNode[2] - 0.0) < EPS);

            /*// without torsion
            const double EPS = 1E-6;
            CPPUNIT_ASSERT(fabs(tipNode[0] - 0.0) < EPS);
            CPPUNIT_ASSERT(fabs(tipNode[1] - (2.0 * RADIUS)) < EPS);
            CPPUNIT_ASSERT(fabs(tipNode[2] - (-1.0)) < EPS);*/

            delete[] surfaceNodeCoors;
            delete[] displacement;

            delete[] curveNodeCoors;
            delete[] curveNodeIDs;
            delete[] curveElems;
            delete[] curveDispRot;

            delete km_Q_O;
            delete km_O_Q;
        }
    }

    CPPUNIT_TEST_SUITE (TestCurveSurfaceCorotateMappers);
    CPPUNIT_TEST (testBending);CPPUNIT_TEST_SUITE_END();
};

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION (EMPIRE::TestCurveSurfaceCorotateMappers);
