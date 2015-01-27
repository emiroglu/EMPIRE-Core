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

#include "CurveSurfaceLinearMapper.h"
#include "KinematicMotion.h"
#include <iostream>
#include <string>
#include <math.h>

using namespace std;

namespace EMPIRE {
/********//**
 * \brief Test the class CurveSurfaceLinearMapper.
 ***********/
class TestCurveSurfaceLinearMapper: public CppUnit::TestFixture {
private:

public:
    void setUp() {
    }
    void tearDown() {
    }

    /***********************************************************************************************
     * \brief Test constructor
     ***********/
    void testConstructor() {
        const double TOL = 1E-6;
        const double DISTURB = 1E-8;

        /*
         * curve:
         *
         *                              ---(3,1)
         *                          ----
         * (0,0)-----------(1,0)----
         *
         * surface: (no elements, just 4 sections)
         *
         *                    /(5/3,1.5)---(3,1.5)
         *                   /
         * (0,.5)--(2/3,.5) /
         * (0,.1)                          (3,.1)
         *(0,-.5)--(2/3,-.5)---(5/3,-.5)---(3,-.5)
         *
         *
         *Above are the coordinates in Q, they will be transformed to O later
         *
         */

        /*
         * Ordering:
         *            -1000
         *           -
         * 10<-100<--
         *
         *
         *
         *        /8---9
         *       /
         * 6---7/
         * 5           0
         * 4---3---2---1
         */
        // curve and surface
        int curveNumNodes = 3;
        int curveNumElements = 2;
        double curveNodeCoors[] = { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 1.0, 0.0 };
        curveNodeCoors[1 * 3 + 0] += DISTURB;
        curveNodeCoors[2 * 3 + 0] -= DISTURB;
        int curveNodeIDs[] = { 100, 10, 1000 };
        int curveElems[] = { 100, 10, 1000, 100 };
        int surfaceNumNodes = 10;
        double surfaceNodeCoors[] = { 3.0, 0.1, 0.0, 3.0, -0.5, 0.0, 5.0 / 3.0, -0.5, 0.0, 2.0
                / 3.0, -0.5, 0.0, 0.0, -0.5, 0.0, 0.0, 0.1, 0.0, 0.0, 0.5, 0.0, 2.0 / 3.0, 0.5, 0.0,
                5.0 / 3.0, 1.5, 0.0, 3.0, 1.5, 0.0, }; // clockwise from (3,.1), corresponding to below

        int surfaceNumSections = 4;
        int surfaceNumRootSectionNodes = 3;
        int surfaceNumNormalSectionNodes = 2;
        int surfaceNumTipSectionNodes = 3;

        KinematicMotion KM_O_Q;
        double axis[] = { 1.0, 1.0, 1.0 };
        KM_O_Q.addRotation(axis, false, M_PI / 2.0);
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
        CurveSurfaceLinearMapper *mapper = new CurveSurfaceLinearMapper(curveNumNodes, curveNumElements,
                curveNodeCoors, curveNodeIDs, curveElems, surfaceNumNodes, surfaceNodeCoors,
                surfaceNumSections, surfaceNumRootSectionNodes, surfaceNumNormalSectionNodes,
                surfaceNumTipSectionNodes, KM_O_Q.getRotationMatrix(),
                KM_O_Q.getTranslationVector());

        /*cout << "============sortedPosToUnsortedPos============" << endl;
         for (int i = 0; i < surfaceNumNodes; i++) {
         cout << (mapper->sortedPosToUnsortedPos)[i] << endl;
         }*/

        // test sorting of surface nodes according to x in Q
        for (int i = 0; i < surfaceNumSections; i++) {
            // compute the x of a section in the Q system
            set<int> tmp;
            if (i == 0) {
                for (int j = 0; j < surfaceNumRootSectionNodes; j++) { // root
                    tmp.insert((mapper->sortedPosToUnsortedPos)[j]);
                }
                set<int>::iterator it = tmp.begin();
                CPPUNIT_ASSERT(*it == 4);
                it++;
                CPPUNIT_ASSERT(*it == 5);
                it++;
                CPPUNIT_ASSERT(*it == 6);
            } else if (i == surfaceNumSections - 1) { // tip
                for (int j = 0; j < surfaceNumTipSectionNodes; j++) {
                    tmp.insert((mapper->sortedPosToUnsortedPos)[surfaceNumNodes - j - 1]);
                }
                set<int>::iterator it = tmp.begin();
                CPPUNIT_ASSERT(*it == 0);
                it++;
                CPPUNIT_ASSERT(*it == 1);
                it++;
                CPPUNIT_ASSERT(*it == 9);
            } else { // normal
                for (int j = 0; j < surfaceNumNormalSectionNodes; j++) {
                    tmp.insert(
                            (mapper->sortedPosToUnsortedPos)[surfaceNumRootSectionNodes
                                    + (i - 1) * surfaceNumNormalSectionNodes + j]);
                }
                if (i == 1) {
                    set<int>::iterator it = tmp.begin();
                    CPPUNIT_ASSERT(*it == 3);
                    it++;
                    CPPUNIT_ASSERT(*it == 7);
                }
                if (i == 2) {
                    set<int>::iterator it = tmp.begin();
                    CPPUNIT_ASSERT(*it == 2);
                    it++;
                    CPPUNIT_ASSERT(*it == 8);
                }
            }
        }

        // test relation between section and curve/beam element
        /*cout << "============sectionToCurveElem============" << endl;
         for (int i = 0; i < surfaceNumSections; i++) {
         cout << (mapper->sectionToCurveElem)[i] << endl;
         }*/
        CPPUNIT_ASSERT(mapper->sectionToCurveElem[0] == 0);
        CPPUNIT_ASSERT(mapper->sectionToCurveElem[1] == 0);
        CPPUNIT_ASSERT(mapper->sectionToCurveElem[2] == 1);
        CPPUNIT_ASSERT(mapper->sectionToCurveElem[3] == 1);

        // test section P
        double sectionP[] = { 0.0, 0.0, 0.0, 2.0 / 3.0, 0.0, 0.0, 5.0 / 3.0, 1.0 / 3.0, 0.0, 3.0,
                1.0, 0.0 };
        for (int i = 0; i < surfaceNumSections; i++) {
            KM_O_Q.move(&sectionP[i * 3]);
        }

        for (int i = 0; i < surfaceNumSections * 3; i++) {
            //cout << "sectionP: " << mapper->sectionP[i] << "  "<< sectionP[i] << endl;
            CPPUNIT_ASSERT(fabs(mapper->sectionP[i] - sectionP[i]) < TOL);
        }

        // test shape functions and their derivatives
        const double SQRT5 = sqrt(5.0);
        CPPUNIT_ASSERT(fabs(mapper->shapeFuncOfSection[0 * 10 + 0] - 0.0) < TOL);
        CPPUNIT_ASSERT(fabs(mapper->shapeFuncOfSection[0 * 10 + 1] - 0.0) < TOL);
        CPPUNIT_ASSERT(fabs(mapper->shapeFuncOfSection[0 * 10 + 2] - 0.0) < TOL);
        CPPUNIT_ASSERT(fabs(mapper->shapeFuncOfSection[0 * 10 + 3] - 0.0) < TOL);
        CPPUNIT_ASSERT(fabs(mapper->shapeFuncOfSection[0 * 10 + 4] - 0.0) < TOL);
        CPPUNIT_ASSERT(fabs(mapper->shapeFuncOfSection[0 * 10 + 5 + 0] - 1.0) < TOL);
        CPPUNIT_ASSERT(fabs(mapper->shapeFuncOfSection[0 * 10 + 5 + 1] - 1.0) < TOL);
        CPPUNIT_ASSERT(fabs(mapper->shapeFuncOfSection[0 * 10 + 5 + 2] - 0.0) < TOL);
        CPPUNIT_ASSERT(fabs(mapper->shapeFuncOfSection[0 * 10 + 5 + 3] - 0.0) < TOL);
        CPPUNIT_ASSERT(fabs(mapper->shapeFuncOfSection[0 * 10 + 5 + 4] - 1.0) < TOL);

        CPPUNIT_ASSERT(fabs(mapper->shapeFuncOfSection[1 * 10 + 0] - 2.0 / 3.0) < TOL);
        CPPUNIT_ASSERT(fabs(mapper->shapeFuncOfSection[1 * 10 + 1] - 20.0 / 27.0) < TOL);
        CPPUNIT_ASSERT(fabs(mapper->shapeFuncOfSection[1 * 10 + 2] - (-4.0) / 27.0) < TOL);
        CPPUNIT_ASSERT(fabs(mapper->shapeFuncOfSection[1 * 10 + 3] - 4.0 / 3.0) < TOL);
        CPPUNIT_ASSERT(fabs(mapper->shapeFuncOfSection[1 * 10 + 4] - 0.0) < TOL);
        CPPUNIT_ASSERT(fabs(mapper->shapeFuncOfSection[1 * 10 + 5 + 0] - 1.0 / 3.0) < TOL);
        CPPUNIT_ASSERT(fabs(mapper->shapeFuncOfSection[1 * 10 + 5 + 1] - 7.0 / 27.0) < TOL);
        CPPUNIT_ASSERT(fabs(mapper->shapeFuncOfSection[1 * 10 + 5 + 2] - 2.0 / 27.0) < TOL);
        CPPUNIT_ASSERT(fabs(mapper->shapeFuncOfSection[1 * 10 + 5 + 3] - (-4.0) / 3.0) < TOL);
        CPPUNIT_ASSERT(fabs(mapper->shapeFuncOfSection[1 * 10 + 5 + 4] - (-1.0) / 3.0) < TOL);

        CPPUNIT_ASSERT(fabs(mapper->shapeFuncOfSection[2 * 10 + 0] - 1.0 / 3.0) < TOL);
        CPPUNIT_ASSERT(fabs(mapper->shapeFuncOfSection[2 * 10 + 1] - 7.0 / 27.0) < TOL);
        CPPUNIT_ASSERT(fabs(mapper->shapeFuncOfSection[2 * 10 + 2] - (-2.0) * SQRT5 / 27.0) < TOL);
        CPPUNIT_ASSERT(fabs(mapper->shapeFuncOfSection[2 * 10 + 3] - 4.0 / 3.0 / SQRT5) < TOL);
        CPPUNIT_ASSERT(fabs(mapper->shapeFuncOfSection[2 * 10 + 4] - (-1.0) / 3.0) < TOL);
        CPPUNIT_ASSERT(fabs(mapper->shapeFuncOfSection[2 * 10 + 5 + 0] - 2.0 / 3.0) < TOL);
        CPPUNIT_ASSERT(fabs(mapper->shapeFuncOfSection[2 * 10 + 5 + 1] - 20.0 / 27.0) < TOL);
        CPPUNIT_ASSERT(fabs(mapper->shapeFuncOfSection[2 * 10 + 5 + 2] - 4.0 * SQRT5 / 27.0) < TOL);
        CPPUNIT_ASSERT(
                fabs(mapper->shapeFuncOfSection[2 * 10 + 5 + 3] - (-4.0) / 3.0 / SQRT5) < TOL);
        CPPUNIT_ASSERT(fabs(mapper->shapeFuncOfSection[2 * 10 + 5 + 4] - 0.0) < TOL);

        CPPUNIT_ASSERT(fabs(mapper->shapeFuncOfSection[3 * 10 + 0] - 1.0) < TOL);
        CPPUNIT_ASSERT(fabs(mapper->shapeFuncOfSection[3 * 10 + 1] - 1.0) < TOL);
        CPPUNIT_ASSERT(fabs(mapper->shapeFuncOfSection[3 * 10 + 2] - 0.0) < TOL);
        CPPUNIT_ASSERT(fabs(mapper->shapeFuncOfSection[3 * 10 + 3] - 0.0) < TOL);
        CPPUNIT_ASSERT(fabs(mapper->shapeFuncOfSection[3 * 10 + 4] - 1.0) < TOL);
        CPPUNIT_ASSERT(fabs(mapper->shapeFuncOfSection[3 * 10 + 5 + 0] - 0.0) < TOL);
        CPPUNIT_ASSERT(fabs(mapper->shapeFuncOfSection[3 * 10 + 5 + 1] - 0.0) < TOL);
        CPPUNIT_ASSERT(fabs(mapper->shapeFuncOfSection[3 * 10 + 5 + 2] - 0.0) < TOL);
        CPPUNIT_ASSERT(fabs(mapper->shapeFuncOfSection[3 * 10 + 5 + 3] - 0.0) < TOL);
        CPPUNIT_ASSERT(fabs(mapper->shapeFuncOfSection[3 * 10 + 5 + 4] - 0.0) < TOL);

        // test rotation from O to element local orientation
        KinematicMotion ROT_O_Q;
        ROT_O_Q.addRotation(KM_O_Q.getRotationMatrix());
        // check x axis direction
        {
            double elemXAxisInQ[] = { -1.0, 0.0, 0.0 };
            double elemXAxisInLocal[] = { 1.0, 0.0, 0.0 };

            double *elemXAxisInO = elemXAxisInQ;
            ROT_O_Q.move(elemXAxisInO);

            double *elemXAxisInO_ = elemXAxisInLocal;
            mapper->ROT_O_ELEM[0]->move(elemXAxisInO_);

            for (int i = 0; i < 3; i++) {
                CPPUNIT_ASSERT(fabs(elemXAxisInO[i] - elemXAxisInO_[i]) < TOL);
            }
        }
        {
            double elemXAxisInQ[] = { -2.0, -1.0, 0.0 };
            for (int i = 0; i < 3; i++)
                elemXAxisInQ[i] /= SQRT5;
            double elemXAxisInLocal[] = { 1.0, 0.0, 0.0 };

            double *elemXAxisInO = elemXAxisInQ;
            ROT_O_Q.move(elemXAxisInO);

            double *elemXAxisInO_ = elemXAxisInLocal;
            mapper->ROT_O_ELEM[1]->move(elemXAxisInO_);

            //cout << elemXAxisInO[0] << '\t' << elemXAxisInO[1] << '\t' << elemXAxisInO[2] << endl;
            //cout << elemXAxisInO_[0] << '\t' << elemXAxisInO_[1] << '\t' << elemXAxisInO_[2] << endl;

            for (int i = 0; i < 3; i++) {
                CPPUNIT_ASSERT(fabs(elemXAxisInO[i] - elemXAxisInO_[i]) < TOL);
            }
        }
        // check whether local y is in global XY plane
        {
            double elemYAxisInLocal[] = { 0.0, 1.0, 0.0 };
            double *elemYAxisInO = elemYAxisInLocal;
            mapper->ROT_O_ELEM[0]->move(elemYAxisInO);
            CPPUNIT_ASSERT(fabs(elemYAxisInO[2] - 0.0) < TOL); // z coordinate is 0.0
        }
        {
            double elemYAxisInLocal[] = { 0.0, 1.0, 0.0 };
            double *elemYAxisInO = elemYAxisInLocal;
            mapper->ROT_O_ELEM[1]->move(elemYAxisInO);
            CPPUNIT_ASSERT(fabs(elemYAxisInO[2] - 0.0) < TOL); // z coordinate is 0.0
        }
        // check the correctness of the rotation matrces
        mapper->ROT_O_ELEM[0]->checkRotationCorrectness();
        mapper->ROT_O_ELEM[1]->checkRotationCorrectness();

        delete mapper;
    }

    /***********************************************************************************************
     * \brief Test consistent and conservative mapping
     ***********/
    void testConsistentConservativeMapping() {
        // define curve and surface in Q
        int curveNumNodes = 2;
        int curveNumElements = 1;
        double curveNodeCoors[] = { 0.0, 0.0, 0.0, 10.0, 0.0, 0.0 };
        int curveNodeIDs[] = { 1, 2 };
        int curveElems[] = { 1, 2 };
        int surfaceNumNodes = 3;
        double surfaceNodeCoors[] = { 0.0, 1.0, 1.0, 5.0, 1.0, 1.0, 10.0, 1.0, 1.0 }; // clockwise from (3,.1), corresponding to below

        int surfaceNumSections = 3;
        int surfaceNumRootSectionNodes = 1;
        int surfaceNumNormalSectionNodes = 1;
        int surfaceNumTipSectionNodes = 1;

        KinematicMotion KM_O_Q;
        double axis[] = { 1.0, 1.0, 1.0 };
        KM_O_Q.addRotation(axis, false, M_PI / 2.0);
        double translation_O_Q[] = { 2.0, 3.0, 4.0 };
        KM_O_Q.addTranslation(translation_O_Q);

        // transform the coordinates to O
        for (int i = 0; i < curveNumNodes; i++) {
            KM_O_Q.move(&curveNodeCoors[i * 3]);
        }
        for (int i = 0; i < surfaceNumNodes; i++) {
            KM_O_Q.move(&surfaceNodeCoors[i * 3]);
        }

        // construct the mapper
        CurveSurfaceLinearMapper *mapper = new CurveSurfaceLinearMapper(curveNumNodes, curveNumElements,
                curveNodeCoors, curveNodeIDs, curveElems, surfaceNumNodes, surfaceNodeCoors,
                surfaceNumSections, surfaceNumRootSectionNodes, surfaceNumNormalSectionNodes,
                surfaceNumTipSectionNodes, KM_O_Q.getRotationMatrix(),
                KM_O_Q.getTranslationVector());
        { // test consistent mapping
          // define deformation of curve/beam
          // analytical disp and rot at x=10, attention to the right hand rule: disp_y' = rot_z, disp_z' = -rot_y!
          // pure bending
            const double disp_x = 0.0; // 0.2; // disp_x
            const double disp_y = -10.0 / 180.0 * M_PI / 20.0 * 10.0 * 10.0; // disp_y = b*x^2, b = rot_z / 2x
            const double disp_z = -20.0 / 180.0 * M_PI / 20.0 * 10.0 * 10.0; // disp_z = a*x^2, a = - rot_y / 2x
            const double rot_x = 0.0; // 5.0 / 180.0 * M_PI; // rot_x
            const double rot_y = 20.0 / 180.0 * M_PI; // rot_y
            const double rot_z = -10.0 / 180.0 * M_PI; // rot_z

            /*// pure torsion
             const double disp_x = 0.2; // disp_x
             const double disp_y = 0.0; // disp_y
             const double disp_z = 0.0; // disp_z
             const double rot_x = 1.0 / 180.0 * M_PI; // rot_x
             const double rot_y = 0.0; // rot_y
             const double rot_z = 0.0; // rot_z*/

            double curveDispRot[12];
            // curve node 1
            curveDispRot[0] = 0.0;
            curveDispRot[1] = 0.0;
            curveDispRot[2] = 0.0;
            curveDispRot[3] = 0.0;
            curveDispRot[4] = 0.0;
            curveDispRot[5] = 0.0;
            // curve node 2
            curveDispRot[6] = disp_x;
            curveDispRot[7] = disp_y;
            curveDispRot[8] = disp_z;
            curveDispRot[9] = rot_x;
            curveDispRot[10] = rot_y;
            curveDispRot[11] = rot_z;

            // rotate the DOFs form Q to O
            KinematicMotion ROT_O_Q;
            ROT_O_Q.addRotation(KM_O_Q.getRotationMatrix());
            for (int i = 0; i < curveNumNodes * 2; i++) {
                ROT_O_Q.move(&curveDispRot[i * 3]);
            }

            // consistent mapping
            double surfaceDisp[9];
            mapper->consistentMapping(curveDispRot, surfaceDisp);

            /*cout << "surfaceDisp" << endl;
             for (int i=0; i<3; i++) {
             for (int j=0; j<3; j++)
             cout << " " << surfaceDisp[i*3+j];
             cout << endl;
             }*/

            // calculate the displacement of beam according to kinematic formula of linear beam (book of Wunderlich)
            // u_x = u - y*phi_z + z*phi_y + ...
            // u_y = v - z*phi_x
            // u_z = w + y*phi_x
            double surfaceDispRef[9];
            // surface node 1
            surfaceDispRef[0] = 0.0;
            surfaceDispRef[1] = 0.0;
            surfaceDispRef[2] = 0.0;
            // surface node 3
            surfaceDispRef[2 * 3 + 0] = disp_x - 1.0 * rot_z + 1.0 * rot_y;
            surfaceDispRef[2 * 3 + 1] = disp_y - 1.0 * rot_x;
            surfaceDispRef[2 * 3 + 2] = disp_z + 1.0 * rot_x;
            // interpolate DOFs on surface node 2
            double node2DispRot[6];
            node2DispRot[0] = disp_x / 2.0; //disp_x
            node2DispRot[1] = disp_y / 4.0; //disp_y
            node2DispRot[2] = disp_z / 4.0; //disp_z
            node2DispRot[3] = rot_x / 2.0; //rot_x
            node2DispRot[4] = rot_y / 2.0; //rot_y
            node2DispRot[5] = rot_z / 2.0; //rot_z
            // surface node 2
            surfaceDispRef[1 * 3 + 0] = node2DispRot[0] - 1.0 * node2DispRot[5]
                    + 1.0 * node2DispRot[4];
            surfaceDispRef[1 * 3 + 1] = node2DispRot[1] - 1.0 * node2DispRot[3];
            surfaceDispRef[1 * 3 + 2] = node2DispRot[2] + 1.0 * node2DispRot[3];
            for (int i = 0; i < surfaceNumNodes; i++) {
                ROT_O_Q.move(&surfaceDispRef[i * 3]);
            }

            /*cout << "surfaceDispRef" << endl;
             for (int i=0; i<3; i++) {
             for (int j=0; j<3; j++)
             cout << " " << surfaceDispRef[i*3+j];
             cout << endl;
             }*/

            // 10 meter long, rotation >20deg, two different linear kinematics give a difference of 0.15 meter, around 10%
            const double TOL = 0.15;
            for (int i = 0; i < surfaceNumNodes; i++) {
                for (int j = 0; j < 3; j++) {
                    CPPUNIT_ASSERT(fabs(surfaceDisp[i * 3 + j] - surfaceDispRef[i * 3 + j]) < TOL);
                }
            }
        }

        { // test conservative mapping
            double surfaceForce[] = { 1.0, 0.0, 0.0, 0.0, -0.5, 0.5, -1.0, 0.0, 0.0 };
            double curveForceMoment[12];
            // rotate the DOFs form Q to O
            KinematicMotion ROT_O_Q;
            ROT_O_Q.addRotation(KM_O_Q.getRotationMatrix());
            for (int i = 0; i < surfaceNumNodes; i++) {
                ROT_O_Q.move(&surfaceForce[i * 3]);
            }
            mapper->conservativeMapping(surfaceForce, curveForceMoment);
            /*cout << "curveForceMoment" << endl;
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 3; j++)
                    cout << " " << curveForceMoment[i * 3 + j];
                cout << endl;
            }*/

            double curveForceMomentRef[12];
            // node1 force
            curveForceMomentRef[0 * 6 + 0] = 1.0;
            curveForceMomentRef[0 * 6 + 1] = -0.25;
            curveForceMomentRef[0 * 6 + 2] = 0.25;
            // node1 moment
            curveForceMomentRef[0 * 6 + 3 + 0] = 0.5;
            curveForceMomentRef[0 * 6 + 3 + 1] = 1.0;
            curveForceMomentRef[0 * 6 + 3 + 2] = -1.0;
            // node2 force
            curveForceMomentRef[1 * 6 + 0] = -1.0;
            curveForceMomentRef[1 * 6 + 1] = -0.25;
            curveForceMomentRef[1 * 6 + 2] = 0.25;
            // node2 moment
            curveForceMomentRef[1 * 6 + 3 + 0] = 0.5;
            curveForceMomentRef[1 * 6 + 3 + 1] = -1.0;
            curveForceMomentRef[1 * 6 + 3 + 2] = 1.0;
            for (int i = 0; i < curveNumNodes * 2; i++) {
                ROT_O_Q.move(&curveForceMomentRef[i * 3]);
            }
            /*cout << "curveForceMomentRef" << endl;
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 3; j++)
                    cout << " " << curveForceMomentRef[i * 3 + j];
                cout << endl;
            }*/

            const double TOL = 1E-10;
            for (int i = 0; i < curveNumNodes * 2; i++) {
                for (int j = 0; j < 3; j++) {
                    CPPUNIT_ASSERT(
                            fabs(curveForceMoment[i * 3 + j] - curveForceMomentRef[i * 3 + j])
                                    < TOL);
                }
            }
        }

        delete mapper;
    }

    CPPUNIT_TEST_SUITE (TestCurveSurfaceLinearMapper);
    CPPUNIT_TEST (testConstructor);
    CPPUNIT_TEST (testConsistentConservativeMapping);CPPUNIT_TEST_SUITE_END();
};

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION (EMPIRE::TestCurveSurfaceLinearMapper);
