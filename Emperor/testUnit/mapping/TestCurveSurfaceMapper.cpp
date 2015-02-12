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

#include "CurveSurfaceMapper.h"
#include "KinematicMotion.h"
#include "FEMesh.h"
#include "SectionMesh.h"
#include "MapperAdapter.h"
#include "MappingFilter.h"
#include "DataField.h"
#include "ConnectionIOSetup.h"

#include <iostream>
#include <string>
#include <math.h>

using namespace std;

namespace EMPIRE {
/********//**
 * \brief Test the class CurveSurfaceCorotate2DMapper.
 ***********/
class TestCurveSurfaceMapper: public CppUnit::TestFixture {
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
        CurveSurfaceMapper *mapper = new CurveSurfaceMapper(EMPIRE_CurveSurfaceMapper_linear,
                curveNumNodes, curveNumElements, curveNodeCoors, curveNodeIDs, curveElems,
                surfaceNumNodes, surfaceNodeCoors, surfaceNumSections, surfaceNumRootSectionNodes,
                surfaceNumNormalSectionNodes, surfaceNumTipSectionNodes, KM_O_Q.getRotationMatrix(),
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
        double shapeFuncOfSectionRef[4 * 10];
        shapeFuncOfSectionRef[0 * 10 + 0] = 0.0;
        shapeFuncOfSectionRef[0 * 10 + 1] = 0.0;
        shapeFuncOfSectionRef[0 * 10 + 2] = 0.0;
        shapeFuncOfSectionRef[0 * 10 + 3] = 0.0;
        shapeFuncOfSectionRef[0 * 10 + 4] = 0.0;
        shapeFuncOfSectionRef[0 * 10 + 5 + 0] = 1.0;
        shapeFuncOfSectionRef[0 * 10 + 5 + 1] = 1.0;
        shapeFuncOfSectionRef[0 * 10 + 5 + 2] = 0.0;
        shapeFuncOfSectionRef[0 * 10 + 5 + 3] = 0.0;
        shapeFuncOfSectionRef[0 * 10 + 5 + 4] = 1.0;

        shapeFuncOfSectionRef[1 * 10 + 0] = 2.0 / 3.0;
        shapeFuncOfSectionRef[1 * 10 + 1] = 20.0 / 27.0;
        shapeFuncOfSectionRef[1 * 10 + 2] = (-4.0) / 27.0;
        shapeFuncOfSectionRef[1 * 10 + 3] = 4.0 / 3.0;
        shapeFuncOfSectionRef[1 * 10 + 4] = 0.0;
        shapeFuncOfSectionRef[1 * 10 + 5 + 0] = 1.0 / 3.0;
        shapeFuncOfSectionRef[1 * 10 + 5 + 1] = 7.0 / 27.0;
        shapeFuncOfSectionRef[1 * 10 + 5 + 2] = 2.0 / 27.0;
        shapeFuncOfSectionRef[1 * 10 + 5 + 3] = (-4.0) / 3.0;
        shapeFuncOfSectionRef[1 * 10 + 5 + 4] = (-1.0) / 3.0;

        shapeFuncOfSectionRef[2 * 10 + 0] = 1.0 / 3.0;
        shapeFuncOfSectionRef[2 * 10 + 1] = 7.0 / 27.0;
        shapeFuncOfSectionRef[2 * 10 + 2] = (-2.0) / 27.0;
        shapeFuncOfSectionRef[2 * 10 + 3] = 4.0 / 3.0;
        shapeFuncOfSectionRef[2 * 10 + 4] = (-1.0) / 3.0;
        shapeFuncOfSectionRef[2 * 10 + 5 + 0] = 2.0 / 3.0;
        shapeFuncOfSectionRef[2 * 10 + 5 + 1] = 20.0 / 27.0;
        shapeFuncOfSectionRef[2 * 10 + 5 + 2] = 4.0 / 27.0;
        shapeFuncOfSectionRef[2 * 10 + 5 + 3] = (-4.0) / 3.0;
        shapeFuncOfSectionRef[2 * 10 + 5 + 4] = 0.0;

        shapeFuncOfSectionRef[3 * 10 + 0] = 1.0;
        shapeFuncOfSectionRef[3 * 10 + 1] = 1.0;
        shapeFuncOfSectionRef[3 * 10 + 2] = 0.0;
        shapeFuncOfSectionRef[3 * 10 + 3] = 0.0;
        shapeFuncOfSectionRef[3 * 10 + 4] = 1.0;
        shapeFuncOfSectionRef[3 * 10 + 5 + 0] = 0.0;
        shapeFuncOfSectionRef[3 * 10 + 5 + 1] = 0.0;
        shapeFuncOfSectionRef[3 * 10 + 5 + 2] = 0.0;
        shapeFuncOfSectionRef[3 * 10 + 5 + 3] = 0.0;
        shapeFuncOfSectionRef[3 * 10 + 5 + 4] = 0.0;

        for (int i = 0; i < 40; i++) {
            CPPUNIT_ASSERT(fabs(mapper->shapeFuncOfSection[i] - shapeFuncOfSectionRef[i]) < TOL);
        }

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
     * \brief Test bending with specific parabolic deformation
     ***********/
    void testParabolicMapping() {
        // counter=0 : linear;
        // counter=1 : corotate3D;
        for (int counter = 0; counter < 2; counter++) {
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

            // test consistent mapping
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

            // construct the mapper
            CurveSurfaceMapper *mapper;
            if (counter == 0) {
                mapper = new CurveSurfaceMapper(EMPIRE_CurveSurfaceMapper_linear, curveNumNodes,
                        curveNumElements, curveNodeCoors, curveNodeIDs, curveElems, surfaceNumNodes,
                        surfaceNodeCoors, surfaceNumSections, surfaceNumRootSectionNodes,
                        surfaceNumNormalSectionNodes, surfaceNumTipSectionNodes,
                        KM_O_Q.getRotationMatrix(), KM_O_Q.getTranslationVector());
            }
            if (counter == 1) {
                mapper = new CurveSurfaceMapper(EMPIRE_CurveSurfaceMapper_corotate3D, curveNumNodes,
                        curveNumElements, curveNodeCoors, curveNodeIDs, curveElems, surfaceNumNodes,
                        surfaceNodeCoors, surfaceNumSections, surfaceNumRootSectionNodes,
                        surfaceNumNormalSectionNodes, surfaceNumTipSectionNodes,
                        KM_O_Q.getRotationMatrix(), KM_O_Q.getTranslationVector());
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

            delete mapper;
        }
    }

    /***********************************************************************************************
     * \brief Bend a beam into half circle
     ***********/
    void testCircleConsistent() {
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

                curveDispRot[i * 6 + 3] = angle / 2.0; // with torsion
                //curveDispRot[i * 6 + 3] = 0.0; // without torsion
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
                CurveSurfaceMapper *mapper = new CurveSurfaceMapper(
                        EMPIRE_CurveSurfaceMapper_corotate2D, curveNumNodes, curveNumElements,
                        curveNodeCoors, curveNodeIDs, curveElems, surfaceNumNodes, surfaceNodeCoors,
                        numSections, numRootSectionNodes, numNormalSectionNodes, numTipSectionNodes,
                        km_O_Q->getRotationMatrix(), km_O_Q->getTranslationVector());

                mapper->consistentMapping(curveDispRot, displacement);
                delete mapper;
            }
            if (counter == 1) {
                CurveSurfaceMapper *mapper = new CurveSurfaceMapper(
                        EMPIRE_CurveSurfaceMapper_corotate3D, curveNumNodes, curveNumElements,
                        curveNodeCoors, curveNodeIDs, curveElems, surfaceNumNodes, surfaceNodeCoors,
                        numSections, numRootSectionNodes, numNormalSectionNodes, numTipSectionNodes,
                        km_O_Q->getRotationMatrix(), km_O_Q->getTranslationVector());

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
            // with torsion
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

    /***********************************************************************************************
     * \brief Bend a beam into half circle, test the conservative mapping of forces/moments
     ***********/
    void testCircleConservative() {
        // length of the beam is 10.0, the beam is bended and prolonged and twisted
        // counter = 0 for corotate2D; counter = 1 for corotate3D
        for (int counter = 0; counter < 2; counter++) {
            // two parameters to control: curveNumNodes and RADIUS
            // curve mesh in Q
            int curveNumNodes = 2; // fix this to 2 so that the forces/moments can be calculated by hand which are the reference values
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

                curveDispRot[i * 6 + 3] = angle / 2.0; // with torsion
                //curveDispRot[i * 6 + 3] = 0.0; // without torsion
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

            // force/moments on the beam/curve nodes
            double *curveForcesMoments = new double[curveNumNodes * 6];

            // forces on the surface nodes
            double *surfaceForces = new double[surfaceNumNodes * 3];
            for (int i = 0; i < surfaceNumNodes * 3; i++) {
                surfaceForces[i] = 0.0;
            }
            surfaceForces[5 * 3 + 0] = 1.0;

            if (counter == 0) {
                CurveSurfaceMapper *mapper = new CurveSurfaceMapper(
                        EMPIRE_CurveSurfaceMapper_corotate2D, curveNumNodes, curveNumElements,
                        curveNodeCoors, curveNodeIDs, curveElems, surfaceNumNodes, surfaceNodeCoors,
                        numSections, numRootSectionNodes, numNormalSectionNodes, numTipSectionNodes,
                        km_O_Q->getRotationMatrix(), km_O_Q->getTranslationVector());

                mapper->consistentMapping(curveDispRot, displacement);
                mapper->conservativeMapping(surfaceForces, curveForcesMoments);
                delete mapper;
            }
            if (counter == 1) {
                CurveSurfaceMapper *mapper = new CurveSurfaceMapper(
                        EMPIRE_CurveSurfaceMapper_corotate3D, curveNumNodes, curveNumElements,
                        curveNodeCoors, curveNodeIDs, curveElems, surfaceNumNodes, surfaceNodeCoors,
                        numSections, numRootSectionNodes, numNormalSectionNodes, numTipSectionNodes,
                        km_O_Q->getRotationMatrix(), km_O_Q->getTranslationVector());

                mapper->consistentMapping(curveDispRot, displacement);
                mapper->conservativeMapping(surfaceForces, curveForcesMoments);
                delete mapper;
            }

            /*for (int i = 0; i < 2; i++) {
             cout << "curveForcesMoments" << endl;
             for (int j = 0; j < 6; j++) {
             cout << curveForcesMoments[i*6 + j] << "  ";
             }
             cout << endl;
             }*/

            // with torsion
            double curveForcesMomentsRef[12];
            curveForcesMomentsRef[0] = 1.0 / 2.0;
            curveForcesMomentsRef[1] = 0.0;
            curveForcesMomentsRef[2] = 0.0;
            curveForcesMomentsRef[3] = 0.0;
            curveForcesMomentsRef[4] = -1.0 / sqrt(2.0) / 2.0;
            curveForcesMomentsRef[5] = 0.0;
            curveForcesMomentsRef[6] = 1.0 / 2.0;
            curveForcesMomentsRef[7] = 0.0;
            curveForcesMomentsRef[8] = 0.0;
            curveForcesMomentsRef[9] = 0.0;
            curveForcesMomentsRef[10] = -1.0 / sqrt(2.0) / 2.0;
            curveForcesMomentsRef[11] = 0.0;
            const double EPS = 2E-1; // corotate3D has bigger error due to big torsion
            for (int i = 0; i < 12; i++)
                CPPUNIT_ASSERT(fabs(curveForcesMoments[i] - curveForcesMomentsRef[i]) < EPS);

            /*// without torsion
             double curveForcesMomentsRef[12];
             curveForcesMomentsRef[0] = 1.0 / 2.0;
             curveForcesMomentsRef[1] = 0.0;
             curveForcesMomentsRef[2] = 0.0;
             curveForcesMomentsRef[3] = 0.0;
             curveForcesMomentsRef[4] = -1.0 / 2.0;
             curveForcesMomentsRef[5] = 0.0;
             curveForcesMomentsRef[6] = 1.0 / 2.0;
             curveForcesMomentsRef[7] = 0.0;
             curveForcesMomentsRef[8] = 0.0;
             curveForcesMomentsRef[9] = 0.0;
             curveForcesMomentsRef[10] = -1.0 / 2.0;
             curveForcesMomentsRef[11] = 0.0;
             const double EPS = 1E-8; // corotate3D has bigger error due to big torsion
             for (int i = 0; i < 12; i++)
             CPPUNIT_ASSERT(fabs(curveForcesMoments[i] - curveForcesMomentsRef[i]) < EPS);*/

            delete[] surfaceNodeCoors;
            delete[] displacement;

            delete[] curveNodeCoors;
            delete[] curveNodeIDs;
            delete[] curveElems;
            delete[] curveDispRot;

            delete[] curveForcesMoments;
            delete[] surfaceForces;

            delete km_Q_O;
            delete km_O_Q;
        }
    }

    void testMapperAdapter() { // no assertion, just test if calls between classes work or not.
        /*
         * 4-----5-----6
         * |     |     | surface mesh B
         * 1-----2-----3
         */
        int numNodesMeshB = 6;
        int numElemsMeshB = 2;
        SectionMesh *surfaceMeshB = new SectionMesh("", numNodesMeshB, numElemsMeshB);
        for (int i = 0; i < numElemsMeshB; i++)
            surfaceMeshB->numNodesPerElem[i] = 4;
        surfaceMeshB->initElems();

        for (int i = 0; i < numNodesMeshB; i++)
            surfaceMeshB->nodeIDs[i] = i + 1;

        surfaceMeshB->nodes[0 * 3 + 0] = 0;
        surfaceMeshB->nodes[0 * 3 + 1] = 0;
        surfaceMeshB->nodes[0 * 3 + 2] = 0;

        surfaceMeshB->nodes[1 * 3 + 0] = 1.5;
        surfaceMeshB->nodes[1 * 3 + 1] = 0;
        surfaceMeshB->nodes[1 * 3 + 2] = 0;

        surfaceMeshB->nodes[2 * 3 + 0] = 3;
        surfaceMeshB->nodes[2 * 3 + 1] = 0;
        surfaceMeshB->nodes[2 * 3 + 2] = 0;

        surfaceMeshB->nodes[3 * 3 + 0] = 0;
        surfaceMeshB->nodes[3 * 3 + 1] = 1;
        surfaceMeshB->nodes[3 * 3 + 2] = 0;

        surfaceMeshB->nodes[4 * 3 + 0] = 1.5;
        surfaceMeshB->nodes[4 * 3 + 1] = 1;
        surfaceMeshB->nodes[4 * 3 + 2] = 0;

        surfaceMeshB->nodes[5 * 3 + 0] = 3;
        surfaceMeshB->nodes[5 * 3 + 1] = 1;
        surfaceMeshB->nodes[5 * 3 + 2] = 0;

        surfaceMeshB->elems[0 * 4 + 0] = 1;
        surfaceMeshB->elems[0 * 4 + 1] = 2;
        surfaceMeshB->elems[0 * 4 + 2] = 5;
        surfaceMeshB->elems[0 * 4 + 3] = 4;

        surfaceMeshB->elems[1 * 4 + 0] = 2;
        surfaceMeshB->elems[1 * 4 + 1] = 3;
        surfaceMeshB->elems[1 * 4 + 2] = 6;
        surfaceMeshB->elems[1 * 4 + 3] = 5;

        surfaceMeshB->setNumSections(3);
        surfaceMeshB->setNumRootSectionNodes(2);
        surfaceMeshB->setNumNormalSectionNodes(2);
        surfaceMeshB->setNumTipSectionNodes(2);
        double rotation[] = { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 };
        double translation[] = { 0.0, 0.5, 0.0 };
        surfaceMeshB->setRotationGlobal2Root(rotation);
        surfaceMeshB->setTranslationGlobal2Root(translation);
        /*
         * 1-----------2 beam mesh A
         */
        int numNodesMeshA = 2;
        int numElemsMeshA = 1;
        FEMesh *beamMeshA = new FEMesh("", numNodesMeshA, numElemsMeshA);
        for (int i = 0; i < numElemsMeshA; i++)
            beamMeshA->numNodesPerElem[i] = 2;
        beamMeshA->initElems();

        for (int i = 0; i < numNodesMeshA; i++)
            beamMeshA->nodeIDs[i] = i + 1;

        beamMeshA->nodes[0 * 3 + 0] = 0;
        beamMeshA->nodes[0 * 3 + 1] = 0.5;
        beamMeshA->nodes[0 * 3 + 2] = 0;

        beamMeshA->nodes[1 * 3 + 0] = 3.0;
        beamMeshA->nodes[1 * 3 + 1] = 0.5;
        beamMeshA->nodes[1 * 3 + 2] = 0;

        beamMeshA->elems[0] = 1;
        beamMeshA->elems[1] = 2;

        DataField *d_A, *d_B, *f_A, *f_B;
        d_A = new DataField("d_A", EMPIRE_DataField_atNode, beamMeshA->numNodes,
                EMPIRE_DataField_doubleVector, EMPIRE_DataField_field);
        d_B = new DataField("d_B", EMPIRE_DataField_atNode, surfaceMeshB->numNodes,
                EMPIRE_DataField_vector, EMPIRE_DataField_field);
        f_A = new DataField("f_A", EMPIRE_DataField_atNode, beamMeshA->numNodes,
                EMPIRE_DataField_doubleVector, EMPIRE_DataField_fieldIntegral);
        f_B = new DataField("f_B", EMPIRE_DataField_atNode, surfaceMeshB->numNodes,
                EMPIRE_DataField_vector, EMPIRE_DataField_fieldIntegral);
        {
            MapperAdapter *mapper = new MapperAdapter("", beamMeshA, surfaceMeshB);
            mapper->initCurveSurfaceMapper(EMPIRE_CurveSurfaceMapper_linear);
            AbstractFilter *filterConsistent = new MappingFilter(mapper);
            ConnectionIOSetup::setupIOForFilter(filterConsistent, beamMeshA, d_A, surfaceMeshB, d_B);
            filterConsistent->filtering();
            AbstractFilter *filterConservative = new MappingFilter(mapper);
            ConnectionIOSetup::setupIOForFilter(filterConservative, surfaceMeshB, f_B, beamMeshA, f_A);
            filterConservative->filtering();
            delete mapper;
            delete filterConsistent;
            delete filterConservative;
        }
        {
            MapperAdapter *mapper = new MapperAdapter("", beamMeshA, surfaceMeshB);
            mapper->initCurveSurfaceMapper(EMPIRE_CurveSurfaceMapper_corotate2D);
            AbstractFilter *filterConsistent = new MappingFilter(mapper);
            ConnectionIOSetup::setupIOForFilter(filterConsistent, beamMeshA, d_A, surfaceMeshB, d_B);
            filterConsistent->filtering();
            AbstractFilter *filterConservative = new MappingFilter(mapper);
            ConnectionIOSetup::setupIOForFilter(filterConservative, surfaceMeshB, f_B, beamMeshA, f_A);
            filterConservative->filtering();
            delete mapper;
            delete filterConsistent;
            delete filterConservative;
        }
        {
            MapperAdapter *mapper = new MapperAdapter("", beamMeshA, surfaceMeshB);
            mapper->initCurveSurfaceMapper(EMPIRE_CurveSurfaceMapper_corotate3D);
            AbstractFilter *filterConsistent = new MappingFilter(mapper);
            ConnectionIOSetup::setupIOForFilter(filterConsistent, beamMeshA, d_A, surfaceMeshB, d_B);
            filterConsistent->filtering();
            AbstractFilter *filterConservative = new MappingFilter(mapper);
            ConnectionIOSetup::setupIOForFilter(filterConservative, surfaceMeshB, f_B, beamMeshA, f_A);
            filterConservative->filtering();
            delete mapper;
            delete filterConsistent;
            delete filterConservative;
        }

        delete surfaceMeshB;
        delete beamMeshA;
        delete d_A;
        delete d_B;
        delete f_A;
        delete f_B;
    }
    void testMemoryLeak() {
        for (int i = 0; i < 1000000000; i++) {
            testCircleConsistent();
            testCircleConservative();
            testMapperAdapter();
        }
    }

    CPPUNIT_TEST_SUITE (TestCurveSurfaceMapper);
    CPPUNIT_TEST (testConstructor);
    CPPUNIT_TEST (testParabolicMapping);
    CPPUNIT_TEST (testCircleConsistent);
    CPPUNIT_TEST (testCircleConservative);
    CPPUNIT_TEST (testMapperAdapter);
    //CPPUNIT_TEST (testMemoryLeak);
    CPPUNIT_TEST_SUITE_END();
};

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION (EMPIRE::TestCurveSurfaceMapper);
