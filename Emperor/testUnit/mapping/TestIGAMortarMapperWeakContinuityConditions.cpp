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
/*
 * TestIGAMortarMapperWeakContinuityConditions.cpp
 *
 *  Created on: May 8, 2013
 *      Author: Andreas Apostolatos
 */

#include "cppunit/TestFixture.h"
#include "cppunit/TestAssert.h"
#include "cppunit/extensions/HelperMacros.h"
#include "WeakIGAPatchContinuityCondition.h"
#include "IGAMortarMapper.h"
#include "FEMesh.h"
#include "IGAPatchSurface.h"
#include "IGAMesh.h"
#include "DataField.h"
#include "MathLibrary.h"
#include <iostream>
#include <math.h>

using namespace std;

namespace EMPIRE {
class TestIGAMortarMapperWeakContinuityConditions: public CppUnit::TestFixture {

private:
    IGAMortarMapper* theMapper;
    IGAMesh* theIGAMesh;
    FEMesh* theFEMesh;
    bool isMappingIGA2FEM;
    bool isWeakConditions;
    double Tol;
    double TolRel10;
    double TolRel100;
    double TolRel1000;
    double TolRel1E8;
    double TolRel1E9;
    double TolRel1E11;

public:
    void setUp() {
        isMappingIGA2FEM = false;
        isWeakConditions = true;

        // Assign a tolerance value (corresponding to maximum accuracy provided by MATLAB)
        Tol = 1e-15;
        TolRel10 = Tol*1e1;
        TolRel100 = Tol*1e2;
        TolRel1000 = Tol*1e3;
        TolRel1E8 = Tol*1e8;
        TolRel1E9 = Tol*1e9;
        TolRel1E11 = Tol*1e11;

        // Patch 1

        // Polynomial orders
        int p1 = 2;
        int q1 = 2;

        // Number of knots at each parametric direction
        int noUKnots1 = 6;
        int noVKnots1 = 8;

        // Knot vectors
        double knotVectorU1[noUKnots1];
        knotVectorU1[0] = 0.0000000000e+00;
        knotVectorU1[1] = 0.0000000000e+00;
        knotVectorU1[2] = 0.0000000000e+00;
        knotVectorU1[3] = 5.0000000000e+00;
        knotVectorU1[4] = 5.0000000000e+00;
        knotVectorU1[5] = 5.0000000000e+00;

        double knotVectorV1[noVKnots1];
        knotVectorV1[0] = 0.00000000;
        knotVectorV1[1] = 0.00000000;
        knotVectorV1[2] = 0.00000000;
        knotVectorV1[3] = 0.3333333333333333;
        knotVectorV1[4] = 0.6666666666666666;
        knotVectorV1[5] = 1.00000000;
        knotVectorV1[6] = 1.00000000;
        knotVectorV1[7] = 1.00000000;

        // Number of Control Points at each parametric direction
        int noUCP1 = noUKnots1 - p1 - 1;
        int noVCP1 = noVKnots1 - q1 - 1;

        // Contol Point net
        double CP1[4*noUCP1*noVCP1];
        CP1[0 * 4 + 0] = -0.0039080073992120;
        CP1[0 * 4 + 1] = -0.0052340654608283;
        CP1[0 * 4 + 2] = -0.0000000001332505;
        CP1[0 * 4 + 3] = 1.0000000;

        CP1[1 * 4 + 0] = 0.1140281817914990;
        CP1[1 * 4 + 1] = 0.2038696439224300;
        CP1[1 * 4 + 2] = 0.0000000004657780;
        CP1[1 * 4 + 3] = 1.0000000;

        CP1[2 * 4 + 0] = 0.3160230427898070;
        CP1[2 * 4 + 1] = 0.2566812221037400;
        CP1[2 * 4 + 2] = 0.0000000002586590;
        CP1[2 * 4 + 3] = 1.0000000;

        CP1[3 * 4 + 0] = -0.0038308271877550;
        CP1[3 * 4 + 1] = -0.0045726953562137;
        CP1[3 * 4 + 2] = 0.0130601936149476;
        CP1[3 * 4 + 3] = 0.9023689270621825;

        CP1[4 * 4 + 0] = 0.1203657359237186;
        CP1[4 * 4 + 1] = 0.2022538651308512;
        CP1[4 * 4 + 2] = 0.0130601942127190;
        CP1[4 * 4 + 3] = 0.9023689270621825;

        CP1[5 * 4 + 0] = 0.3292518318252257;
        CP1[5 * 4 + 1] = 0.2573832786712888;
        CP1[5 * 4 + 2] = 0.0130601940087691;
        CP1[5 * 4 + 3] = 0.9023689270621825;

        CP1[6 * 4 + 0] = -0.0037332280765124;
        CP1[6 * 4 + 1] = -0.0041310995735450;
        CP1[6 * 4 + 2] = 0.0499999998652280;
        CP1[6 * 4 + 3] = 0.8372815451036375;

        CP1[7 * 4 + 0] = 0.1253036273101167;
        CP1[7 * 4 + 1] = 0.2010548510603844;
        CP1[7 * 4 + 2] = 0.0500000004655298;
        CP1[7 * 4 + 3] = 0.8372815451036375;

        CP1[8 * 4 + 0] = 0.3398741993916353;
        CP1[8 * 4 + 1] = 0.2580315193260972;
        CP1[8 * 4 + 2] = 0.0500000002615814;
        CP1[8 * 4 + 3] = 0.8372815451036375;

        CP1[9 * 4 + 0] = -0.0037302104905827;
        CP1[9 * 4 + 1] = -0.0048093172135863;
        CP1[9 * 4 + 2] = 0.0869398061142833;
        CP1[9 * 4 + 3] = 0.9023689270621824;

        CP1[10 * 4 + 0] = 0.1200644264840723;
        CP1[10 * 4 + 1] = 0.2024974689562943;
        CP1[10 * 4 + 2] = 0.0869398067218844;
        CP1[10 * 4 + 3] = 0.9023689270621824;

        CP1[11 * 4 + 0] = 0.3295002434320827;
        CP1[11 * 4 + 1] = 0.2576316988787566;
        CP1[11 * 4 + 2] = 0.0869398065109145;
        CP1[11 * 4 + 3] = 0.9023689270621824;

        CP1[12 * 4 + 0] = -0.0037718173276041;
        CP1[12 * 4 + 1] = -0.0055543457781634;
        CP1[12 * 4 + 2] = 0.0999999998609426;
        CP1[12 * 4 + 3] = 1.0000000;

        CP1[13 * 4 + 0] = 0.1136203433778480;
        CP1[13 * 4 + 1] = 0.2041993747063200;
        CP1[13 * 4 + 2] = 0.1000000004732761;
        CP1[13 * 4 + 3] = 1.0000000;

        CP1[14 * 4 + 0] = 0.3163592811625310;
        CP1[14 * 4 + 1] = 0.2570174721178500;
        CP1[14 * 4 + 2] = 0.1000000002566551;
        CP1[14 * 4 + 3] = 1.0000000;

        int dofIndex1[noUCP1*noVCP1];
        for(int i = 0; i < noUCP1*noVCP1; i++)
            dofIndex1[i] = i;

        // Patch 1 Trimming Loop 1

        // Number of trimming curves
        int noTrCurves1_1 = 4;

        // Patch 1 Trimming Curve 1

        // Direction
        int d1_1 = 1;

        // Polynomial orders
        int p1_1 = 1;

        // Number of knots at each parametric direction
        int noUKnots1_1 = 4;

        // Knot vectors
        double knotVectorU1_1[noUKnots1_1];
        knotVectorU1_1[0] = 0.0000000000e+00;
        knotVectorU1_1[1] = 0.0000000000e+00;
        knotVectorU1_1[2] = 1.0000000000e+00;
        knotVectorU1_1[3] = 1.0000000000e+00;

        // Number of Control Points at each parametric direction
        int noUCP1_1 = noUKnots1_1 - p1_1 - 1;

        // Contol Point net
        double CP1_1[4*noUCP1_1];
        CP1_1[0 * 4 + 0] = 0.0000000;
        CP1_1[0 * 4 + 1] = 0.0000000;
        CP1_1[0 * 4 + 2] = 0.0000000;
        CP1_1[0 * 4 + 3] = 1.0000000;

        CP1_1[1 * 4 + 0] = 5.0000000;
        CP1_1[1 * 4 + 1] = 0.0000000;
        CP1_1[1 * 4 + 2] = 0.0000000;
        CP1_1[1 * 4 + 3] = 1.0000000;

        // Patch 1 Trimming Curve 2

        // Direction
        int d1_2 = 1;

        // Polynomial orders
        int p1_2 = 1;

        // Number of knots at each parametric direction
        int noUKnots1_2 = 4;

        // Knot vectors
        double knotVectorU1_2[noUKnots1_2];
        knotVectorU1_2[0] = 0.0000000000e+00;
        knotVectorU1_2[1] = 0.0000000000e+00;
        knotVectorU1_2[2] = 5.0000000000e+00;
        knotVectorU1_2[3] = 5.0000000000e+00;

        // Number of Control Points at each parametric direction
        int noUCP1_2 = noUKnots1_2 - p1_2 - 1;

        // Contol Point net
        double CP1_2[4*noUCP1_2];
        CP1_2[0 * 4 + 0] = 5.0000000;
        CP1_2[0 * 4 + 1] = 0.0000000;
        CP1_2[0 * 4 + 2] = 0.0000000;
        CP1_2[0 * 4 + 3] = 1.0000000;

        CP1_2[1 * 4 + 0] = 5.0000000;
        CP1_2[1 * 4 + 1] = 1.0000000;
        CP1_2[1 * 4 + 2] = 0.0000000;
        CP1_2[1 * 4 + 3] = 1.0000000;

        // Patch 1 Trimming Curve 3

        // Direction
        int d1_3 = 1;

        // Polynomial orders
        int p1_3 = 1;

        // Number of knots at each parametric direction
        int noUKnots1_3 = 4;

        // Knot vectors
        double knotVectorU1_3[noUKnots1_3];
        knotVectorU1_3[0] = 0.0000000000e+00;
        knotVectorU1_3[1] = 0.0000000000e+00;
        knotVectorU1_3[2] = 1.0000000000e+00;
        knotVectorU1_3[3] = 1.0000000000e+00;

        // Number of Control Points at each parametric direction
        int noUCP1_3 = noUKnots1_3 - p1_3 - 1;

        // Contol Point net
        double CP1_3[4*noUCP1_3];
        CP1_3[0 * 4 + 0] = 5.0000000;
        CP1_3[0 * 4 + 1] = 1.0000000;
        CP1_3[0 * 4 + 2] = 0.0000000;
        CP1_3[0 * 4 + 3] = 1.0000000;

        CP1_3[1 * 4 + 0] = 0.0000000;
        CP1_3[1 * 4 + 1] = 1.0000000;
        CP1_3[1 * 4 + 2] = 0.0000000;
        CP1_3[1 * 4 + 3] = 1.0000000;

        // Patch 1 Trimming Curve 4

        // Direction
        int d1_4 = 1;

        // Polynomial orders
        int p1_4 = 1;

        // Number of knots at each parametric direction
        int noUKnots1_4 = 4;

        // Knot vectors
        double knotVectorU1_4[noUKnots1_4];
        knotVectorU1_4[0] = 0.0000000000e+00;
        knotVectorU1_4[1] = 0.0000000000e+00;
        knotVectorU1_4[2] = 5.0000000000e+00;
        knotVectorU1_4[3] = 5.0000000000e+00;

        // Number of Control Points at each parametric direction
        int noUCP1_4 = noUKnots1_4 - p1_4 - 1;

        // Contol Point net
        double CP1_4[4*noUCP1_3];
        CP1_4[0 * 4 + 0] = 0.0000000;
        CP1_4[0 * 4 + 1] = 1.0000000;
        CP1_4[0 * 4 + 2] = 0.0000000;
        CP1_4[0 * 4 + 3] = 1.0000000;

        CP1_4[1 * 4 + 0] = 0.0000000;
        CP1_4[1 * 4 + 1] = 0.0000000;
        CP1_4[1 * 4 + 2] = 0.0000000;
        CP1_4[1 * 4 + 3] = 1.0000000;

        // Patch 2

        // Polynomial orders
        int p2 = 2;
        int q2 = 2;

        // Number of knots at each parametric direction
        int noUKnots2 = 6;
        int noVKnots2 = 7;

        // Knot vectors
        double knotVectorU2[noUKnots2];
        knotVectorU2[0] = 0.0000000000e+00;
        knotVectorU2[1] = 0.0000000000e+00;
        knotVectorU2[2] = 0.0000000000e+00;
        knotVectorU2[3] = 5.0000000000e+00;
        knotVectorU2[4] = 5.0000000000e+00;
        knotVectorU2[5] = 5.0000000000e+00;

        double knotVectorV2[noVKnots2];
        knotVectorV2[0] = 0.00000000;
        knotVectorV2[1] = 0.00000000;
        knotVectorV2[2] = 0.00000000;
        knotVectorV2[3] = 0.5000000000000000;
        knotVectorV2[4] = 1.00000000;
        knotVectorV2[5] = 1.00000000;
        knotVectorV2[6] = 1.00000000;

        // Number of Control Points at each parametric direction
        int noUCP2 = noUKnots2 - p2 - 1;
        int noVCP2 = noVKnots2 - q2 - 1;

        // Contol Point net
        double CP2[4*noUCP2*noVCP2];
        CP2[0 * 4 + 0] = 0.3160038569049120;
        CP2[0 * 4 + 1] = 0.2566887936212000;
        CP2[0 * 4 + 2] = 0.0000000002582877;
        CP2[0 * 4 + 3] = 1.0000000;

        CP2[1 * 4 + 0] = 0.7364859004510570;
        CP2[1 * 4 + 1] = 0.3526685822793200;
        CP2[1 * 4 + 2] = 0.0000000005144890;
        CP2[1 * 4 + 3] = 1.0000000;

        CP2[2 * 4 + 0] = 1.0093446589356891;
        CP2[2 * 4 + 1] = -0.0074272016071648;
        CP2[2 * 4 + 2] = -0.0000000000605252;
        CP2[2 * 4 + 3] = 1.0000000;

        CP2[3 * 4 + 0] = 0.3369879570374723;
        CP2[3 * 4 + 1] = 0.2578016940013632;
        CP2[3 * 4 + 2] = 0.0207106783801090;
        CP2[3 * 4 + 3] = 0.8535533905932737;

        CP2[4 * 4 + 0] = 0.7466464001136912;
        CP2[4 * 4 + 1] = 0.3516596333694540;
        CP2[4 * 4 + 2] = 0.0207106786301504;
        CP2[4 * 4 + 3] = 0.8535533905932737;

        CP2[5 * 4 + 0] = 1.0093699098678452;
        CP2[5 * 4 + 1] = -0.0076537680707467;
        CP2[5 * 4 + 2] = 0.0207106780584029;
        CP2[5 * 4 + 3] = 0.8535533905932737;

        CP2[6 * 4 + 0] = 0.3371855088540514;
        CP2[6 * 4 + 1] = 0.2579987704666362;
        CP2[6 * 4 + 2] = 0.0792893221416341;
        CP2[6 * 4 + 3] = 0.8535533905932737;

        CP2[7 * 4 + 0] = 0.7463812657955617;
        CP2[7 * 4 + 1] = 0.3518599426152602;
        CP2[7 * 4 + 2] = 0.0792893223966293;
        CP2[7 * 4 + 3] = 0.8535533905932737;

        CP2[8 * 4 + 0] = 1.0096523717583676;
        CP2[8 * 4 + 1] = -0.0082828461960626;
        CP2[8 * 4 + 2] = 0.0792893218178884;
        CP2[8 * 4 + 3] = 0.8535533905932737;

        CP2[9 * 4 + 0] = 0.3163410989506300;
        CP2[9 * 4 + 1] = 0.2570252241914800;
        CP2[9 * 4 + 2] = 0.1000000002562983;
        CP2[9 * 4 + 3] = 1.0000000;

        CP2[10 * 4 + 0] = 0.7360332878586530;
        CP2[10 * 4 + 1] = 0.3530105315511700;
        CP2[10 * 4 + 2] = 0.1000000005209562;
        CP2[10 * 4 + 3] = 1.0000000;

        CP2[11 * 4 + 0] = 1.0098268515444271;
        CP2[11 * 4 + 1] = -0.0085011051407877;
        CP2[11 * 4 + 2] = 0.0999999999340034;
        CP2[11 * 4 + 3] = 1.0000000;

        int dofIndex2[noUCP2*noVCP2];
        for(int i = 0; i < noUCP2*noVCP2; i++)
            dofIndex2[i] = i;

        // Patch 2 Trimming Loop 1

        // Number of trimming curves
        int noTrCurves2_1 = 4;

        // Patch 2 Trimming Curve 1

        // Direction
        int d2_1 = 1;

        // Polynomial orders
        int p2_1 = 1;

        // Number of knots at each parametric direction
        int noUKnots2_1 = 4;

        // Knot vectors
        double knotVectorU2_1[noUKnots2_1];
        knotVectorU2_1[0] = 0.0000000000e+00;
        knotVectorU2_1[1] = 0.0000000000e+00;
        knotVectorU2_1[2] = 1.0000000000e+00;
        knotVectorU2_1[3] = 1.0000000000e+00;

        // Number of Control Points at each parametric direction
        int noUCP2_1 = noUKnots2_1 - p2_1 - 1;

        // Contol Point net
        double CP2_1[4*noUCP2_1];
        CP2_1[0 * 4 + 0] = 0.0000000;
        CP2_1[0 * 4 + 1] = 0.0000000;
        CP2_1[0 * 4 + 2] = 0.0000000;
        CP2_1[0 * 4 + 3] = 1.0000000;

        CP2_1[1 * 4 + 0] = 5.0000000;
        CP2_1[1 * 4 + 1] = 0.0000000;
        CP2_1[1 * 4 + 2] = 0.0000000;
        CP2_1[1 * 4 + 3] = 1.0000000;

        // Patch 2 Trimming Curve 2

        // Direction
        int d2_2 = 1;

        // Polynomial orders
        int p2_2 = 1;

        // Number of knots at each parametric direction
        int noUKnots2_2 = 4;

        // Knot vectors
        double knotVectorU2_2[noUKnots2_2];
        knotVectorU2_2[0] = 0.0000000000e+00;
        knotVectorU2_2[1] = 0.0000000000e+00;
        knotVectorU2_2[2] = 5.0000000000e+00;
        knotVectorU2_2[3] = 5.0000000000e+00;

        // Number of Control Points at each parametric direction
        int noUCP2_2 = noUKnots2_2 - p2_2 - 1;

        // Contol Point net
        double CP2_2[4*noUCP2_2];
        CP2_2[0 * 4 + 0] = 5.0000000;
        CP2_2[0 * 4 + 1] = 0.0000000;
        CP2_2[0 * 4 + 2] = 0.0000000;
        CP2_2[0 * 4 + 3] = 1.0000000;

        CP2_2[1 * 4 + 0] = 5.0000000;
        CP2_2[1 * 4 + 1] = 1.0000000;
        CP2_2[1 * 4 + 2] = 0.0000000;
        CP2_2[1 * 4 + 3] = 1.0000000;

        // Patch 1 Trimming Curve 3

        // Direction
        int d2_3 = 1;

        // Polynomial orders
        int p2_3 = 1;

        // Number of knots at each parametric direction
        int noUKnots2_3 = 4;

        // Knot vectors
        double knotVectorU2_3[noUKnots2_3];
        knotVectorU2_3[0] = 0.0000000000e+00;
        knotVectorU2_3[1] = 0.0000000000e+00;
        knotVectorU2_3[2] = 1.0000000000e+00;
        knotVectorU2_3[3] = 1.0000000000e+00;

        // Number of Control Points at each parametric direction
        int noUCP2_3 = noUKnots2_3 - p2_3 - 1;

        // Contol Point net
        double CP2_3[4*noUCP2_3];
        CP2_3[0 * 4 + 0] = 5.0000000;
        CP2_3[0 * 4 + 1] = 1.0000000;
        CP2_3[0 * 4 + 2] = 0.0000000;
        CP2_3[0 * 4 + 3] = 1.0000000;

        CP2_3[1 * 4 + 0] = 0.0000000;
        CP2_3[1 * 4 + 1] = 1.0000000;
        CP2_3[1 * 4 + 2] = 0.0000000;
        CP2_3[1 * 4 + 3] = 1.0000000;

        // Patch 1 Trimming Curve 4

        // Direction
        int d2_4 = 1;

        // Polynomial orders
        int p2_4 = 1;

        // Number of knots at each parametric direction
        int noUKnots2_4 = 4;

        // Knot vectors
        double knotVectorU2_4[noUKnots2_4];
        knotVectorU2_4[0] = 0.0000000000e+00;
        knotVectorU2_4[1] = 0.0000000000e+00;
        knotVectorU2_4[2] = 5.0000000000e+00;
        knotVectorU2_4[3] = 5.0000000000e+00;

        // Number of Control Points at each parametric direction
        int noUCP2_4 = noUKnots2_4 - p2_4 - 1;

        // Contol Point net
        double CP2_4[4*noUCP2_3];
        CP2_4[0 * 4 + 0] = 0.0000000;
        CP2_4[0 * 4 + 1] = 1.0000000;
        CP2_4[0 * 4 + 2] = 0.0000000;
        CP2_4[0 * 4 + 3] = 1.0000000;

        CP2_4[1 * 4 + 0] = 0.0000000;
        CP2_4[1 * 4 + 1] = 0.0000000;
        CP2_4[1 * 4 + 2] = 0.0000000;
        CP2_4[1 * 4 + 3] = 1.0000000;

        // Total number of Gauss Points
        int noCPs = noUCP1*noVCP1 + noUCP2*noVCP2;

        // Construct the NURBS multipatch geometry
        theIGAMesh = new IGAMesh("modifiedCavity2Patches", noCPs);

        // Add the underlying patches of the multipatch geometry
        theIGAMesh->addPatch(p1, noUKnots1, knotVectorU1, q1, noVKnots1, knotVectorV1,
                             noUCP1, noVCP1, CP1, dofIndex1);
        theIGAMesh->getSurfacePatch(0)->addTrimLoop(0, noTrCurves1_1);
        theIGAMesh->getSurfacePatch(0)->addTrimCurve(d1_1, p1_1, noUKnots1_1, knotVectorU1_1, noUCP1_1, CP1_1);
        theIGAMesh->getSurfacePatch(0)->addTrimCurve(d1_2, p1_2, noUKnots1_2, knotVectorU1_2, noUCP1_2, CP1_2);
        theIGAMesh->getSurfacePatch(0)->addTrimCurve(d1_3, p1_3, noUKnots1_3, knotVectorU1_3, noUCP1_3, CP1_3);
        theIGAMesh->getSurfacePatch(0)->addTrimCurve(d1_4, p1_4, noUKnots1_4, knotVectorU1_4, noUCP1_4, CP1_4);
        theIGAMesh->getSurfacePatch(0)->linearizeTrimming();

        theIGAMesh->addPatch(p2, noUKnots2, knotVectorU2, q2, noVKnots2, knotVectorV2,
                             noUCP2, noVCP2, CP2, dofIndex2);
        theIGAMesh->getSurfacePatch(1)->addTrimLoop(0, noTrCurves2_1);
        theIGAMesh->getSurfacePatch(1)->addTrimCurve(d2_1, p2_1, noUKnots2_1, knotVectorU2_1, noUCP2_1, CP2_1);
        theIGAMesh->getSurfacePatch(1)->addTrimCurve(d2_2, p2_2, noUKnots2_2, knotVectorU2_2, noUCP2_2, CP2_2);
        theIGAMesh->getSurfacePatch(1)->addTrimCurve(d2_3, p2_3, noUKnots2_3, knotVectorU2_3, noUCP2_3, CP2_3);
        theIGAMesh->getSurfacePatch(1)->addTrimCurve(d2_4, p2_4, noUKnots2_4, knotVectorU2_4, noUCP2_4, CP2_4);
        theIGAMesh->getSurfacePatch(1)->linearizeTrimming();

        // Initialize auxiliary variables
        int connectionCtr = 0;
        int masterPatchIndex = 0;
        int masterPatchBLIndex = 0;
        int masterPatchBLTrCurveIndex = 1;
        int slavePatchIndex = 1;
        int slavePatchBLIndex = 0;
        int slavePatchBLTrCurveIndex = 3;

        // Add the coupling data corresponding to the coupling between the patches
        theIGAMesh->addWeakContinuityCondition(connectionCtr,
                                               masterPatchIndex, masterPatchBLIndex, masterPatchBLTrCurveIndex,
                                               slavePatchIndex, slavePatchBLIndex, slavePatchBLTrCurveIndex);

        // Compute a bounding box for the multipatch geometry
        theIGAMesh->computeBoundingBox();

        // FEM mesh

        // Number of nodes
        int noNodes = 62;

        // Number of elements
        int noElements = 30;

        // Number of nodes per element
        int noNodesPerElement = 4;

        // Initialize the Finite Element mesh
        theFEMesh = new FEMesh("FEMMesh", noNodes, noElements);

        // Nodal ids and coordinates
        for (int i = 0; i < noNodes; i++)
            theFEMesh->nodeIDs[i] = i + 1;

        theFEMesh->nodes[0 * 3 + 0] = -0.000202020000;
        theFEMesh->nodes[0 * 3 + 1] = 0.000136472000;
        theFEMesh->nodes[0 * 3 + 2] = -0.000000000002;

        theFEMesh->nodes[1 * 3 + 0] = 0.019840400000;
        theFEMesh->nodes[1 * 3 + 1] = 0.034537500000;
        theFEMesh->nodes[1 * 3 + 2] = 0.000000000002;

        theFEMesh->nodes[2 * 3 + 0] = 0.042295900000;
        theFEMesh->nodes[2 * 3 + 1] = 0.067276700000;
        theFEMesh->nodes[2 * 3 + 2] = 0.000000000027;

        theFEMesh->nodes[3 * 3 + 0] = 0.067189800000;
        theFEMesh->nodes[3 * 3 + 1] = 0.098408700000;
        theFEMesh->nodes[3 * 3 + 2] = 0.000000000073;

        theFEMesh->nodes[4 * 3 + 0] = 0.094350500000;
        theFEMesh->nodes[4 * 3 + 1] = 0.127573000000;
        theFEMesh->nodes[4 * 3 + 2] = 0.000000000135;

        theFEMesh->nodes[5 * 3 + 0] = 0.123590900000;
        theFEMesh->nodes[5 * 3 + 1] = 0.154396000000;
        theFEMesh->nodes[5 * 3 + 2] = 0.000000000207;

        theFEMesh->nodes[6 * 3 + 0] = 0.154985100000;
        theFEMesh->nodes[6 * 3 + 1] = 0.178969000000;
        theFEMesh->nodes[6 * 3 + 2] = 0.000000000292;

        theFEMesh->nodes[7 * 3 + 0] = 0.188190800000;
        theFEMesh->nodes[7 * 3 + 1] = 0.201040000000;
        theFEMesh->nodes[7 * 3 + 2] = 0.000000000363;

        theFEMesh->nodes[8 * 3 + 0] = 0.222854800000;
        theFEMesh->nodes[8 * 3 + 1] = 0.220377000000;
        theFEMesh->nodes[8 * 3 + 2] = 0.000000000392;

        theFEMesh->nodes[9 * 3 + 0] = 0.259083900000;
        theFEMesh->nodes[9 * 3 + 1] = 0.237019000000;
        theFEMesh->nodes[9 * 3 + 2] = 0.000000000380;

        theFEMesh->nodes[10 * 3 + 0] = 0.296473600000;
        theFEMesh->nodes[10 * 3 + 1] = 0.250920000000;
        theFEMesh->nodes[10 * 3 + 2] = 0.000000000346;

        theFEMesh->nodes[11 * 3 + 0] = 0.334619300000;
        theFEMesh->nodes[11 * 3 + 1] = 0.262049000000;
        theFEMesh->nodes[11 * 3 + 2] = 0.000000000310;

        theFEMesh->nodes[12 * 3 + 0] = 0.373627500000;
        theFEMesh->nodes[12 * 3 + 1] = 0.270402000000;
        theFEMesh->nodes[12 * 3 + 2] = 0.000000000266;

        theFEMesh->nodes[13 * 3 + 0] = 0.413146500000;
        theFEMesh->nodes[13 * 3 + 1] = 0.276031000000;
        theFEMesh->nodes[13 * 3 + 2] = 0.000000000232;

        theFEMesh->nodes[14 * 3 + 0] = 0.452829900000;
        theFEMesh->nodes[14 * 3 + 1] = 0.278993000000;
        theFEMesh->nodes[14 * 3 + 2] = 0.000000000221;

        theFEMesh->nodes[15 * 3 + 0] = 0.492767920000;
        theFEMesh->nodes[15 * 3 + 1] = 0.279270000000;
        theFEMesh->nodes[15 * 3 + 2] = 0.000000000233;

        theFEMesh->nodes[16 * 3 + 0] = 0.532663480000;
        theFEMesh->nodes[16 * 3 + 1] = 0.276930000000;
        theFEMesh->nodes[16 * 3 + 2] = 0.000000000258;

        theFEMesh->nodes[17 * 3 + 0] = 0.572223120000;
        theFEMesh->nodes[17 * 3 + 1] = 0.272040000000;
        theFEMesh->nodes[17 * 3 + 2] = 0.000000000288;

        theFEMesh->nodes[18 * 3 + 0] = 0.611526300000;
        theFEMesh->nodes[18 * 3 + 1] = 0.264582000000;
        theFEMesh->nodes[18 * 3 + 2] = 0.000000000323;

        theFEMesh->nodes[19 * 3 + 0] = 0.650297900000;
        theFEMesh->nodes[19 * 3 + 1] = 0.254628000000;
        theFEMesh->nodes[19 * 3 + 2] = 0.000000000356;

        theFEMesh->nodes[20 * 3 + 0] = 0.688264200000;
        theFEMesh->nodes[20 * 3 + 1] = 0.242256000000;
        theFEMesh->nodes[20 * 3 + 2] = 0.000000000378;

        theFEMesh->nodes[21 * 3 + 0] = 0.725500300000;
        theFEMesh->nodes[21 * 3 + 1] = 0.227439000000;
        theFEMesh->nodes[21 * 3 + 2] = 0.000000000392;

        theFEMesh->nodes[22 * 3 + 0] = 0.761743200000;
        theFEMesh->nodes[22 * 3 + 1] = 0.210287000000;
        theFEMesh->nodes[22 * 3 + 2] = 0.000000000381;

        theFEMesh->nodes[23 * 3 + 0] = 0.796732300000;
        theFEMesh->nodes[23 * 3 + 1] = 0.190915000000;
        theFEMesh->nodes[23 * 3 + 2] = 0.000000000332;

        theFEMesh->nodes[24 * 3 + 0] = 0.830535100000;
        theFEMesh->nodes[24 * 3 + 1] = 0.169280000000;
        theFEMesh->nodes[24 * 3 + 2] = 0.000000000241;

        theFEMesh->nodes[25 * 3 + 0] = 0.862930000000;
        theFEMesh->nodes[25 * 3 + 1] = 0.145561000000;
        theFEMesh->nodes[25 * 3 + 2] = 0.000000000143;

        theFEMesh->nodes[26 * 3 + 0] = 0.893703800000;
        theFEMesh->nodes[26 * 3 + 1] = 0.119942000000;
        theFEMesh->nodes[26 * 3 + 2] = 0.000000000076;

        theFEMesh->nodes[27 * 3 + 0] = 0.922899800000;
        theFEMesh->nodes[27 * 3 + 1] = 0.092359500000;
        theFEMesh->nodes[27 * 3 + 2] = 0.000000000034;

        theFEMesh->nodes[28 * 3 + 0] = 0.950416000000;
        theFEMesh->nodes[28 * 3 + 1] = 0.063074900000;
        theFEMesh->nodes[28 * 3 + 2] = 0.000000000011;

        theFEMesh->nodes[29 * 3 + 0] = 0.976161180000;
        theFEMesh->nodes[29 * 3 + 1] = 0.032342000000;
        theFEMesh->nodes[29 * 3 + 2] = -0.000000000001;

        theFEMesh->nodes[30 * 3 + 0] = 1.000147000000;
        theFEMesh->nodes[30 * 3 + 1] = 0.000122479000;
        theFEMesh->nodes[30 * 3 + 2] = -0.000000000001;

        theFEMesh->nodes[31 * 3 + 0] = -0.000202020000;
        theFEMesh->nodes[31 * 3 + 1] = 0.000136472000;
        theFEMesh->nodes[31 * 3 + 2] = 0.100000000000;

        theFEMesh->nodes[32 * 3 + 0] = 0.019840400000;
        theFEMesh->nodes[32 * 3 + 1] = 0.034537500000;
        theFEMesh->nodes[32 * 3 + 2] = 0.100000000000;

        theFEMesh->nodes[33 * 3 + 0] = 0.042295900000;
        theFEMesh->nodes[33 * 3 + 1] = 0.067276700000;
        theFEMesh->nodes[33 * 3 + 2] = 0.100000000000;

        theFEMesh->nodes[34 * 3 + 0] = 0.067189800000;
        theFEMesh->nodes[34 * 3 + 1] = 0.098408700000;
        theFEMesh->nodes[34 * 3 + 2] = 0.100000000000;

        theFEMesh->nodes[35 * 3 + 0] = 0.094350500000;
        theFEMesh->nodes[35 * 3 + 1] = 0.127573000000;
        theFEMesh->nodes[35 * 3 + 2] = 0.100000000000;

        theFEMesh->nodes[36 * 3 + 0] = 0.123590900000;
        theFEMesh->nodes[36 * 3 + 1] = 0.154396000000;
        theFEMesh->nodes[36 * 3 + 2] = 0.100000000000;

        theFEMesh->nodes[37 * 3 + 0] = 0.154985100000;
        theFEMesh->nodes[37 * 3 + 1] = 0.178969000000;
        theFEMesh->nodes[37 * 3 + 2] = 0.100000000000;

        theFEMesh->nodes[38 * 3 + 0] = 0.188190800000;
        theFEMesh->nodes[38 * 3 + 1] = 0.201040000000;
        theFEMesh->nodes[38 * 3 + 2] = 0.100000000000;

        theFEMesh->nodes[39 * 3 + 0] = 0.222854800000;
        theFEMesh->nodes[39 * 3 + 1] = 0.220377000000;
        theFEMesh->nodes[39 * 3 + 2] = 0.100000000000;

        theFEMesh->nodes[40 * 3 + 0] = 0.259083900000;
        theFEMesh->nodes[40 * 3 + 1] = 0.237019000000;
        theFEMesh->nodes[40 * 3 + 2] = 0.100000000000;

        theFEMesh->nodes[41 * 3 + 0] = 0.296473600000;
        theFEMesh->nodes[41 * 3 + 1] = 0.250920000000;
        theFEMesh->nodes[41 * 3 + 2] = 0.100000000000;

        theFEMesh->nodes[42 * 3 + 0] = 0.334619300000;
        theFEMesh->nodes[42 * 3 + 1] = 0.262049000000;
        theFEMesh->nodes[42 * 3 + 2] = 0.100000000000;

        theFEMesh->nodes[43 * 3 + 0] = 0.373627500000;
        theFEMesh->nodes[43 * 3 + 1] = 0.270402000000;
        theFEMesh->nodes[43 * 3 + 2] = 0.100000000000;

        theFEMesh->nodes[44 * 3 + 0] = 0.413146500000;
        theFEMesh->nodes[44 * 3 + 1] = 0.276031000000;
        theFEMesh->nodes[44 * 3 + 2] = 0.100000000000;

        theFEMesh->nodes[45 * 3 + 0] = 0.452829900000;
        theFEMesh->nodes[45 * 3 + 1] = 0.278993000000;
        theFEMesh->nodes[45 * 3 + 2] = 0.100000000000;

        theFEMesh->nodes[46 * 3 + 0] = 0.492767920000;
        theFEMesh->nodes[46 * 3 + 1] = 0.279270000000;
        theFEMesh->nodes[46 * 3 + 2] = 0.100000000000;

        theFEMesh->nodes[47 * 3 + 0] = 0.532663480000;
        theFEMesh->nodes[47 * 3 + 1] = 0.276930000000;
        theFEMesh->nodes[47 * 3 + 2] = 0.100000000000;

        theFEMesh->nodes[48 * 3 + 0] = 0.572223120000;
        theFEMesh->nodes[48 * 3 + 1] = 0.272040000000;
        theFEMesh->nodes[48 * 3 + 2] = 0.100000000000;

        theFEMesh->nodes[49 * 3 + 0] = 0.611526300000;
        theFEMesh->nodes[49 * 3 + 1] = 0.264582000000;
        theFEMesh->nodes[49 * 3 + 2] = 0.100000000000;

        theFEMesh->nodes[50 * 3 + 0] = 0.650297900000;
        theFEMesh->nodes[50 * 3 + 1] = 0.254628000000;
        theFEMesh->nodes[50 * 3 + 2] = 0.100000000000;

        theFEMesh->nodes[51 * 3 + 0] = 0.688264200000;
        theFEMesh->nodes[51 * 3 + 1] = 0.242256000000;
        theFEMesh->nodes[51 * 3 + 2] = 0.100000000000;

        theFEMesh->nodes[52 * 3 + 0] = 0.725500300000;
        theFEMesh->nodes[52 * 3 + 1] = 0.227439000000;
        theFEMesh->nodes[52 * 3 + 2] = 0.100000000000;

        theFEMesh->nodes[53 * 3 + 0] = 0.761743200000;
        theFEMesh->nodes[53 * 3 + 1] = 0.210287000000;
        theFEMesh->nodes[53 * 3 + 2] = 0.100000000000;

        theFEMesh->nodes[54 * 3 + 0] = 0.796732300000;
        theFEMesh->nodes[54 * 3 + 1] = 0.190915000000;
        theFEMesh->nodes[54 * 3 + 2] = 0.100000000000;

        theFEMesh->nodes[55 * 3 + 0] = 0.830535100000;
        theFEMesh->nodes[55 * 3 + 1] = 0.169280000000;
        theFEMesh->nodes[55 * 3 + 2] = 0.100000000000;

        theFEMesh->nodes[56 * 3 + 0] = 0.862930000000;
        theFEMesh->nodes[56 * 3 + 1] = 0.145561000000;
        theFEMesh->nodes[56 * 3 + 2] = 0.100000000000;

        theFEMesh->nodes[57 * 3 + 0] = 0.893703800000;
        theFEMesh->nodes[57 * 3 + 1] = 0.119942000000;
        theFEMesh->nodes[57 * 3 + 2] = 0.100000000000;

        theFEMesh->nodes[58 * 3 + 0] = 0.922899800000;
        theFEMesh->nodes[58 * 3 + 1] = 0.092359500000;
        theFEMesh->nodes[58 * 3 + 2] = 0.100000000000;

        theFEMesh->nodes[59 * 3 + 0] = 0.950416000000;
        theFEMesh->nodes[59 * 3 + 1] = 0.063074900000;
        theFEMesh->nodes[59 * 3 + 2] = 0.100000000000;

        theFEMesh->nodes[60 * 3 + 0] = 0.976161180000;
        theFEMesh->nodes[60 * 3 + 1] = 0.032342000000;
        theFEMesh->nodes[60 * 3 + 2] = 0.100000000000;

        theFEMesh->nodes[61 * 3 + 0] = 1.000147000000;
        theFEMesh->nodes[61 * 3 + 1] = 0.000122479000;
        theFEMesh->nodes[61 * 3 + 2] = 0.100000000000;

        // Assign the nodal ids
        for (int i = 0; i < noNodes; i++)
            theFEMesh->nodeIDs[i] = i + 1;

        // Get the number of nodes per element
        for (int i = 0; i < noElements; i++)
            theFEMesh->numNodesPerElem[i] = noNodesPerElement;

        // Initialize the elements of the mesh
        theFEMesh->initElems();

        // Assign the elements in the mesh
        theFEMesh->elems[0 * 4 + 0] = 1;
        theFEMesh->elems[0 * 4 + 1] = 2;
        theFEMesh->elems[0 * 4 + 2] = 33;
        theFEMesh->elems[0 * 4 + 3] = 32;
        theFEMesh->elems[1 * 4 + 0] = 2;
        theFEMesh->elems[1 * 4 + 1] = 3;
        theFEMesh->elems[1 * 4 + 2] = 34;
        theFEMesh->elems[1 * 4 + 3] = 33;
        theFEMesh->elems[2 * 4 + 0] = 3;
        theFEMesh->elems[2 * 4 + 1] = 4;
        theFEMesh->elems[2 * 4 + 2] = 35;
        theFEMesh->elems[2 * 4 + 3] = 34;
        theFEMesh->elems[3 * 4 + 0] = 4;
        theFEMesh->elems[3 * 4 + 1] = 5;
        theFEMesh->elems[3 * 4 + 2] = 36;
        theFEMesh->elems[3 * 4 + 3] = 35;
        theFEMesh->elems[4 * 4 + 0] = 5;
        theFEMesh->elems[4 * 4 + 1] = 6;
        theFEMesh->elems[4 * 4 + 2] = 37;
        theFEMesh->elems[4 * 4 + 3] = 36;
        theFEMesh->elems[5 * 4 + 0] = 6;
        theFEMesh->elems[5 * 4 + 1] = 7;
        theFEMesh->elems[5 * 4 + 2] = 38;
        theFEMesh->elems[5 * 4 + 3] = 37;
        theFEMesh->elems[6 * 4 + 0] = 7;
        theFEMesh->elems[6 * 4 + 1] = 8;
        theFEMesh->elems[6 * 4 + 2] = 39;
        theFEMesh->elems[6 * 4 + 3] = 38;
        theFEMesh->elems[7 * 4 + 0] = 8;
        theFEMesh->elems[7 * 4 + 1] = 9;
        theFEMesh->elems[7 * 4 + 2] = 40;
        theFEMesh->elems[7 * 4 + 3] = 39;
        theFEMesh->elems[8 * 4 + 0] = 9;
        theFEMesh->elems[8 * 4 + 1] = 10;
        theFEMesh->elems[8 * 4 + 2] = 41;
        theFEMesh->elems[8 * 4 + 3] = 40;
        theFEMesh->elems[9 * 4 + 0] = 10;
        theFEMesh->elems[9 * 4 + 1] = 11;
        theFEMesh->elems[9 * 4 + 2] = 42;
        theFEMesh->elems[9 * 4 + 3] = 41;
        theFEMesh->elems[10 * 4 + 0] = 11;
        theFEMesh->elems[10 * 4 + 1] = 12;
        theFEMesh->elems[10 * 4 + 2] = 43;
        theFEMesh->elems[10 * 4 + 3] = 42;
        theFEMesh->elems[11 * 4 + 0] = 12;
        theFEMesh->elems[11 * 4 + 1] = 13;
        theFEMesh->elems[11 * 4 + 2] = 44;
        theFEMesh->elems[11 * 4 + 3] = 43;
        theFEMesh->elems[12 * 4 + 0] = 13;
        theFEMesh->elems[12 * 4 + 1] = 14;
        theFEMesh->elems[12 * 4 + 2] = 45;
        theFEMesh->elems[12 * 4 + 3] = 44;
        theFEMesh->elems[13 * 4 + 0] = 14;
        theFEMesh->elems[13 * 4 + 1] = 15;
        theFEMesh->elems[13 * 4 + 2] = 46;
        theFEMesh->elems[13 * 4 + 3] = 45;
        theFEMesh->elems[14 * 4 + 0] = 15;
        theFEMesh->elems[14 * 4 + 1] = 16;
        theFEMesh->elems[14 * 4 + 2] = 47;
        theFEMesh->elems[14 * 4 + 3] = 46;
        theFEMesh->elems[15 * 4 + 0] = 16;
        theFEMesh->elems[15 * 4 + 1] = 17;
        theFEMesh->elems[15 * 4 + 2] = 48;
        theFEMesh->elems[15 * 4 + 3] = 47;
        theFEMesh->elems[16 * 4 + 0] = 17;
        theFEMesh->elems[16 * 4 + 1] = 18;
        theFEMesh->elems[16 * 4 + 2] = 49;
        theFEMesh->elems[16 * 4 + 3] = 48;
        theFEMesh->elems[17 * 4 + 0] = 18;
        theFEMesh->elems[17 * 4 + 1] = 19;
        theFEMesh->elems[17 * 4 + 2] = 50;
        theFEMesh->elems[17 * 4 + 3] = 49;
        theFEMesh->elems[18 * 4 + 0] = 19;
        theFEMesh->elems[18 * 4 + 1] = 20;
        theFEMesh->elems[18 * 4 + 2] = 51;
        theFEMesh->elems[18 * 4 + 3] = 50;
        theFEMesh->elems[19 * 4 + 0] = 20;
        theFEMesh->elems[19 * 4 + 1] = 21;
        theFEMesh->elems[19 * 4 + 2] = 52;
        theFEMesh->elems[19 * 4 + 3] = 51;
        theFEMesh->elems[20 * 4 + 0] = 21;
        theFEMesh->elems[20 * 4 + 1] = 22;
        theFEMesh->elems[20 * 4 + 2] = 53;
        theFEMesh->elems[20 * 4 + 3] = 52;
        theFEMesh->elems[21 * 4 + 0] = 22;
        theFEMesh->elems[21 * 4 + 1] = 23;
        theFEMesh->elems[21 * 4 + 2] = 54;
        theFEMesh->elems[21 * 4 + 3] = 53;
        theFEMesh->elems[22 * 4 + 0] = 23;
        theFEMesh->elems[22 * 4 + 1] = 24;
        theFEMesh->elems[22 * 4 + 2] = 55;
        theFEMesh->elems[22 * 4 + 3] = 54;
        theFEMesh->elems[23 * 4 + 0] = 24;
        theFEMesh->elems[23 * 4 + 1] = 25;
        theFEMesh->elems[23 * 4 + 2] = 56;
        theFEMesh->elems[23 * 4 + 3] = 55;
        theFEMesh->elems[24 * 4 + 0] = 25;
        theFEMesh->elems[24 * 4 + 1] = 26;
        theFEMesh->elems[24 * 4 + 2] = 57;
        theFEMesh->elems[24 * 4 + 3] = 56;
        theFEMesh->elems[25 * 4 + 0] = 26;
        theFEMesh->elems[25 * 4 + 1] = 27;
        theFEMesh->elems[25 * 4 + 2] = 58;
        theFEMesh->elems[25 * 4 + 3] = 57;
        theFEMesh->elems[26 * 4 + 0] = 27;
        theFEMesh->elems[26 * 4 + 1] = 28;
        theFEMesh->elems[26 * 4 + 2] = 59;
        theFEMesh->elems[26 * 4 + 3] = 58;
        theFEMesh->elems[27 * 4 + 0] = 28;
        theFEMesh->elems[27 * 4 + 1] = 29;
        theFEMesh->elems[27 * 4 + 2] = 60;
        theFEMesh->elems[27 * 4 + 3] = 59;
        theFEMesh->elems[28 * 4 + 0] = 29;
        theFEMesh->elems[28 * 4 + 1] = 30;
        theFEMesh->elems[28 * 4 + 2] = 61;
        theFEMesh->elems[28 * 4 + 3] = 60;
        theFEMesh->elems[29 * 4 + 0] = 30;
        theFEMesh->elems[29 * 4 + 1] = 31;
        theFEMesh->elems[29 * 4 + 2] = 62;
        theFEMesh->elems[29 * 4 + 3] = 61;

        // Initialize the isogeometric mortar-based mapper
        theMapper = new IGAMortarMapper("Test IGA Mortar Mapper", theIGAMesh, theFEMesh, isWeakConditions, isMappingIGA2FEM);

        // Assign the mapper's properties
        theMapper->setParametersProjection(0.05, 2, 1e-3);
        theMapper->setParametersNewtonRaphson(10, 1e-6);
        theMapper->setParametersNewtonRaphsonBoundary(0, 1e-6);
        theMapper->setParametersBisection(100, 1e-6);

    }

    void tearDown() {

        delete theMapper;
        delete theFEMesh;
        delete theIGAMesh;

    }

    void addPrecomputedGPData() {

        // Initialize auxiliary variables
        int noGPsConnection = 12;

        // Initialize the necessary for the coupling constituents for each Gauss Point
        double* trCurveMasterGPs = new double[noGPsConnection*2]; // The parametric coordinates of the Gauss Points on the master curve
        double* trCurveSlaveGPs = new double[noGPsConnection*2]; // The parametric coordinates of the Gauss Points on the slave
        double* trCurveGPWeights = new double[noGPsConnection]; // The Gauss Point weights for each Gauss Point
        double* trCurveMasterGPTangents = new double[noGPsConnection*3]; // The parametric components of the tangent to the trimming curve vectors on the master curve
        double* trCurveSlaveGPTangents = new double[noGPsConnection*3]; // The parametric components of the tangent to the trimming curve vectors on the slave curve
        double* trCurveGPJacobianProducts = new double[noGPsConnection]; // The products of the Jacobian transformations with the Gauss weights at each Gauss Point

        // Assign the necessary for the coupling constituents for each Gauss Point

        // The parametric coordinates of the Gauss Points on the master curve
        trCurveMasterGPs[2*0 + 0] = 5;
        trCurveMasterGPs[2*0 + 1] = 0.037567221793086;

        trCurveMasterGPs[2*1 + 0] = 5;
        trCurveMasterGPs[2*1 + 1] = 0.166666666666667;

        trCurveMasterGPs[2*2 + 0] = 5;
        trCurveMasterGPs[2*2 + 1] = 0.295766111540247;

        trCurveMasterGPs[2*3 + 0] = 5;
        trCurveMasterGPs[2*3 + 1] = 0.352116944229876;

        trCurveMasterGPs[2*4 + 0] = 5;
        trCurveMasterGPs[2*4 + 1] = 0.416666666666667;

        trCurveMasterGPs[2*5 + 0] = 5;
        trCurveMasterGPs[2*5 + 1] = 0.481216389103457;

        trCurveMasterGPs[2*6 + 0] = 5;
        trCurveMasterGPs[2*6 + 1] = 0.518783610896543;

        trCurveMasterGPs[2*7 + 0] = 5;
        trCurveMasterGPs[2*7 + 1] = 0.583333333333333;

        trCurveMasterGPs[2*8 + 0] = 5;
        trCurveMasterGPs[2*8 + 1] = 0.647883055770124;

        trCurveMasterGPs[2*9 + 0] = 5;
        trCurveMasterGPs[2*9 + 1] = 0.704233888459753;

        trCurveMasterGPs[2*10 + 0] = 5;
        trCurveMasterGPs[2*10 + 1] = 0.833333333333333;

        trCurveMasterGPs[2*11 + 0] = 5;
        trCurveMasterGPs[2*11 + 1] = 0.962432778206914;

        // The parametric coordinates of the Gauss Points on the slave
        trCurveSlaveGPs[2*0 + 0] = 0;
        trCurveSlaveGPs[2*0 + 1] = 0.037567221793086;

        trCurveSlaveGPs[2*1 + 0] = 0;
        trCurveSlaveGPs[2*1 + 1] = 0.166666666666667;

        trCurveSlaveGPs[2*2 + 0] = 0;
        trCurveSlaveGPs[2*2 + 1] = 0.295766111540247;

        trCurveSlaveGPs[2*3 + 0] = 0;
        trCurveSlaveGPs[2*3 + 1] = 0.352116944229876;

        trCurveSlaveGPs[2*4 + 0] = 0;
        trCurveSlaveGPs[2*4 + 1] = 0.416666666666667;

        trCurveSlaveGPs[2*5 + 0] = 0;
        trCurveSlaveGPs[2*5 + 1] = 0.481216389103457;

        trCurveSlaveGPs[2*6 + 0] = 0;
        trCurveSlaveGPs[2*6 + 1] = 0.518783610896543;

        trCurveSlaveGPs[2*7 + 0] = 0;
        trCurveSlaveGPs[2*7 + 1] = 0.583333333333333;

        trCurveSlaveGPs[2*8 + 0] = 0;
        trCurveSlaveGPs[2*8 + 1] = 0.647883055770124;

        trCurveSlaveGPs[2*9 + 0] = 0;
        trCurveSlaveGPs[2*9 + 1] = 0.704233888459753;

        trCurveSlaveGPs[2*10 + 0] = 0;
        trCurveSlaveGPs[2*10 + 1] = 0.833333333333333;

        trCurveSlaveGPs[2*11 + 0] = 0;
        trCurveSlaveGPs[2*11 + 1] = 0.962432778206914;

        // The Gauss Point weights for each Gauss Point
        trCurveGPWeights[0] = 0.555555555555556;

        trCurveGPWeights[1] = 0.888888888888889;

        trCurveGPWeights[2] = 0.555555555555556;

        trCurveGPWeights[3] = 0.555555555555556;

        trCurveGPWeights[4] = 0.888888888888889;

        trCurveGPWeights[5] = 0.555555555555556;

        trCurveGPWeights[6] = 0.555555555555556;

        trCurveGPWeights[7] = 0.888888888888889;

        trCurveGPWeights[8] = 0.555555555555556;

        trCurveGPWeights[9] = 0.555555555555556;

        trCurveGPWeights[10] = 0.888888888888889;

        trCurveGPWeights[11] = 0.555555555555556;

        // The parametric components of the tangent to the trimming curve vectors on the master curve
        trCurveMasterGPTangents[3*0 + 0] = 0.672567413298545;
        trCurveMasterGPTangents[3*0 + 1] = 0.035933526804658;
        trCurveMasterGPTangents[3*0 + 2] = 0.739162942943085;

        trCurveMasterGPTangents[3*1 + 0] = 0.518642424700041;
        trCurveMasterGPTangents[3*1 + 1] = 0.028618645620947;
        trCurveMasterGPTangents[3*1 + 2] = 0.854512146446196;

        trCurveMasterGPTangents[3*2 + 0] = 0.334518371397290;
        trCurveMasterGPTangents[3*2 + 1] = 0.019707912983376;
        trCurveMasterGPTangents[3*2 + 2] = 0.942183133665396;

        trCurveMasterGPTangents[3*3 + 0] = 0.246600757308154;
        trCurveMasterGPTangents[3*3 + 1] = 0.015405035206640;
        trCurveMasterGPTangents[3*3 + 2] = 0.968994711742705;

        trCurveMasterGPTangents[3*4 + 0] = 0.141959112266735;
        trCurveMasterGPTangents[3*4 + 1] = 0.010248227973116;
        trCurveMasterGPTangents[3*4 + 2] = 0.989819470543923;

        trCurveMasterGPTangents[3*5 + 0] = 0.034767353927938;
        trCurveMasterGPTangents[3*5 + 1] = 0.004928153317023;
        trCurveMasterGPTangents[3*5 + 2] = 0.999383282032341;

        trCurveMasterGPTangents[3*6 + 0] = -0.028052951297859;
        trCurveMasterGPTangents[3*6 + 1] = 0.001793096951330;
        trCurveMasterGPTangents[3*6 + 2] = 0.999604830283849;

        trCurveMasterGPTangents[3*7 + 0] = -0.135434083146395;
        trCurveMasterGPTangents[3*7 + 1] = -0.003595176569249;
        trCurveMasterGPTangents[3*7 + 2] = 0.990779836203650;

        trCurveMasterGPTangents[3*8 + 0] = -0.240485196250373;
        trCurveMasterGPTangents[3*8 + 1] = -0.008903380554228;
        trCurveMasterGPTangents[3*8 + 2] = 0.970611972004841;

        trCurveMasterGPTangents[3*9 + 0] = -0.328908196739009;
        trCurveMasterGPTangents[3*9 + 1] = -0.013401303233212;
        trCurveMasterGPTangents[3*9 + 2] = 0.944266807205222;

        trCurveMasterGPTangents[3*10 + 0] = -0.514499528145185;
        trCurveMasterGPTangents[3*10 + 1] = -0.022943557518762;
        trCurveMasterGPTangents[3*10 + 2] = 0.857183661012484;

        trCurveMasterGPTangents[3*11 + 0] = -0.669943646656621;
        trCurveMasterGPTangents[3*11 + 1] = -0.031071504836861;
        trCurveMasterGPTangents[3*11 + 2] = 0.741761465628677;

        // The parametric components of the tangent to the trimming curve vectors on the slave curve
        trCurveSlaveGPTangents[3*0 + 0] = -0.672674587247831;
        trCurveSlaveGPTangents[3*0 + 1] = -0.035915780365672;
        trCurveSlaveGPTangents[3*0 + 2] = -0.739066273342036;

        trCurveSlaveGPTangents[3*1 + 0] = -0.518754457663932;
        trCurveSlaveGPTangents[3*1 + 1] = -0.028606853848593;
        trCurveSlaveGPTangents[3*1 + 2] = -0.854444533347065;

        trCurveSlaveGPTangents[3*2 + 0] = -0.334609248379411;
        trCurveSlaveGPTangents[3*2 + 1] = -0.019701829572375;
        trCurveSlaveGPTangents[3*2 + 2] = -0.942150990452415;

        trCurveSlaveGPTangents[3*3 + 0] = -0.246673935138192;
        trCurveSlaveGPTangents[3*3 + 1] = -0.015401266353883;
        trCurveSlaveGPTangents[3*3 + 2] = -0.968976145587772;

        trCurveSlaveGPTangents[3*4 + 0] = -0.142006906143282;
        trCurveSlaveGPTangents[3*4 + 1] = -0.010246968609278;
        trCurveSlaveGPTangents[3*4 + 2] = -0.989812627845257;

        trCurveSlaveGPTangents[3*5 + 0] = -0.034786160987690;
        trCurveSlaveGPTangents[3*5 + 1] = -0.004929280359464;
        trCurveSlaveGPTangents[3*5 + 2] = -0.999382622021654;

        trCurveSlaveGPTangents[3*6 + 0] = 0.028051720822084;
        trCurveSlaveGPTangents[3*6 + 1] = -0.001795576153355;
        trCurveSlaveGPTangents[3*6 + 2] = -0.999604860364933;

        trCurveSlaveGPTangents[3*7 + 0] = 0.135462489104364;
        trCurveSlaveGPTangents[3*7 + 1] = 0.003590395052853;
        trCurveSlaveGPTangents[3*7 + 2] = -0.990775970191554;

        trCurveSlaveGPTangents[3*8 + 0] = 0.240540374596422;
        trCurveSlaveGPTangents[3*8 + 1] = 0.008896271275466;
        trCurveSlaveGPTangents[3*8 + 2] = -0.970598364178719;

        trCurveSlaveGPTangents[3*9 + 0] = 0.328982756872293;
        trCurveSlaveGPTangents[3*9 + 1] = 0.013392103616478;
        trCurveSlaveGPTangents[3*9 + 2] = -0.944240963547669;

        trCurveSlaveGPTangents[3*10 + 0] = 0.514599927063366;
        trCurveSlaveGPTangents[3*10 + 1] = 0.022929302160749;
        trCurveSlaveGPTangents[3*10 + 2] = -0.857123772957441;

        trCurveSlaveGPTangents[3*11 + 0] = 0.670043695210792;
        trCurveSlaveGPTangents[3*11 + 1] = 0.031051982197132;
        trCurveSlaveGPTangents[3*11 + 2] = -0.741671909209117;

        // The products of the Jacobian transformations with the Gauss weights at each Gauss Point
        trCurveGPJacobianProducts[0] = 0.009522090368294;

        trCurveGPJacobianProducts[1] = 0.016198303357023;

        trCurveGPJacobianProducts[2] = 0.010570156158216;

        trCurveGPJacobianProducts[3] = 0.005351640216454;

        trCurveGPJacobianProducts[4] = 0.008643208768130;

        trCurveGPJacobianProducts[5] = 0.005423354699369;

        trCurveGPJacobianProducts[6] = 0.005422152689563;

        trCurveGPJacobianProducts[7] = 0.008634830881354;

        trCurveGPJacobianProducts[8] = 0.005342723166447;

        trCurveGPJacobianProducts[9] = 0.010546831441419;

        trCurveGPJacobianProducts[10] = 0.016147819419794;

        trCurveGPJacobianProducts[11] = 0.009488732784896;

        // Add the Gauss Point data of the weak continuity condition
        theIGAMesh->getWeakIGAPatchContinuityConditions()[0]->addWeakContinuityConditionGPData(noGPsConnection,
                                                                                               trCurveMasterGPs, trCurveSlaveGPs, trCurveGPWeights,
                                                                                               trCurveMasterGPTangents, trCurveSlaveGPTangents,
                                                                                               trCurveGPJacobianProducts);

        // delete the pointers after object creation
        delete trCurveMasterGPs;
        delete trCurveSlaveGPs;
        delete trCurveGPWeights;
        delete trCurveMasterGPTangents;
        delete trCurveSlaveGPTangents;
        delete trCurveGPJacobianProducts;

    }

    /***********************************************************************************************
     * \brief Test case: testComputeBOperatorMatricesIGAPatchContinuityConditions with the precomputed GP data
     ***********/

    void testComputeBOperatorMatricesIGAPatchContinuityConditions() {

        // Add the precomputed GP data to the continuity condition
        addPrecomputedGPData();

        // Get the weak patch continuity conditions
        std::vector<WeakIGAPatchContinuityCondition*> weakIGAPatchContinuityConditions = theIGAMesh->getWeakIGAPatchContinuityConditions();

        // Initialize constant array sizes
        const int noCoord = 3;

        // Initialize varying array sizes
        int indexMaster;
        int indexSlave;
        int pMaster;
        int qMaster;
        int pSlave;
        int qSlave;
        int noLocalBasisFctsMaster;
        int noLocalBasisFctsSlave;
        int noDOFsLocMaster;
        int noDOFsLocSlave;
        int noGPsOnContCond;
        int uKnotSpanMaster;
        int uKnotSpanSlave;
        int vKnotSpanMaster;
        int vKnotSpanSlave;
        double uGPMaster;
        double vGPMaster;
        double uGPSlave;
        double vGPSlave;
        double tangentTrCurveVctMaster[noCoord];
        double tangentTrCurveVctSlave[noCoord];
        double normalTrCurveVctMaster[noCoord];
        double normalTrCurveVctSlave[noCoord];

        // Initialize pointers
        double* trCurveMasterGPs;
        double* trCurveSlaveGPs;
        double* trCurveGPWeights;
        double* trCurveMasterGPTangents;
        double* trCurveSlaveGPTangents;
        double* trCurveGPJacobianProducts;
        IGAPatchSurface* patchMaster;
        IGAPatchSurface* patchSlave;

        // Get the index of the master and slave patches
        indexMaster = weakIGAPatchContinuityConditions[0]->getMasterPatchIndex();
        indexSlave = weakIGAPatchContinuityConditions[0]->getSlavePatchIndex();

        // Get the number of Gauss Points for the given condition
        noGPsOnContCond = weakIGAPatchContinuityConditions[0]->getTrCurveNumGP();

        // Get the parametric coordinates of the Gauss Points
        trCurveMasterGPs = weakIGAPatchContinuityConditions[0]->getTrCurveMasterGPs();
        trCurveSlaveGPs = weakIGAPatchContinuityConditions[0]->getTrCurveSlaveGPs();

        // Get the corresponding Gauss weights
        trCurveGPWeights = weakIGAPatchContinuityConditions[0]->getTrCurveGPWeights();

        // Get the tangent vectors at the trimming curve of the given condition in the Cartesian space
        trCurveMasterGPTangents = weakIGAPatchContinuityConditions[0]->getTrCurveMasterGPTangents();
        trCurveSlaveGPTangents = weakIGAPatchContinuityConditions[0]->getTrCurveSlaveGPTangents();

        // Get the product of the Jacobian transformations
        trCurveGPJacobianProducts = weakIGAPatchContinuityConditions[0]->getTrCurveGPJacobianProducts();

        // Get the master and the slave patch
        patchMaster = theIGAMesh->getSurfacePatch(indexMaster);
        patchSlave = theIGAMesh->getSurfacePatch(indexSlave);

        // Get the polynomial orders of the master and the slave patch
        pMaster = patchMaster->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree();
        qMaster = patchMaster->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree();
        pSlave = patchSlave->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree();
        qSlave = patchSlave->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree();

        // get the number of local basis functions for master and slave patch
        noLocalBasisFctsMaster = (pMaster + 1)*(qMaster + 1);
        noLocalBasisFctsSlave = (pSlave + 1)*(qSlave + 1);

        // get the number of the local DOFs for the master and slave patch
        noDOFsLocMaster = noCoord*noLocalBasisFctsMaster;
        noDOFsLocSlave = noCoord*noLocalBasisFctsSlave;

        // Initialize pointers
        double* BDisplacementsGCMaster = new double[noCoord*noDOFsLocMaster];
        double* BDisplacementsGCSlave = new double[noCoord*noDOFsLocSlave];
        double* BOperatorOmegaTMaster = new double[noDOFsLocMaster];
        double* BOperatorOmegaTSlave = new double[noDOFsLocSlave];
        double* BOperatorOmegaNMaster = new double[noDOFsLocMaster];
        double* BOperatorOmegaNSlave = new double[noDOFsLocSlave];

        // Expected values for the displacement B-operator matrix of the master patch by MATLAB
        static double expectedBOperatorDisplacementMatrixMaster[][81] = {{0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 8.0433384464093327e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.9023366204009751e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 5.4324933189692538e-03, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 8.0433384464093327e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.9023366204009751e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 5.4324933189692538e-03, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 8.0433384464093327e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.9023366204009751e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 5.4324933189692538e-03},
                                                                         {0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 2.7214119759298838e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 6.1392940120350581e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.1392940120350581e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 2.7214119759298838e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 6.1392940120350581e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.1392940120350581e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 2.7214119759298838e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 6.1392940120350581e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.1392940120350581e-01},
                                                                         {0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.4466796706116473e-02, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 6.1013467358138296e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 3.7539852971250054e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.4466796706116473e-02, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 6.1013467358138296e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 3.7539852971250054e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.4466796706116473e-02, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 6.1013467358138296e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 3.7539852971250054e-01},
                                                                         {0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 4.6374022315160546e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 5.3460608600939241e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.6536908390020913e-03, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 4.6374022315160546e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 5.3460608600939241e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.6536908390020913e-03, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 4.6374022315160546e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 5.3460608600939241e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.6536908390020913e-03},
                                                                         {0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 2.9592460606833648e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 6.7119488214629286e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 3.2880511785370700e-02, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 2.9592460606833648e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 6.7119488214629286e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 3.2880511785370700e-02, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 2.9592460606833648e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 6.7119488214629286e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 3.2880511785370700e-02},
                                                                         {0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.6357456796751663e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 7.3241003321827625e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.0401539881420721e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.6357456796751663e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 7.3241003321827625e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.0401539881420721e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.6357456796751663e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 7.3241003321827625e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.0401539881420721e-01},
                                                                         {0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.0401539881420724e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 7.3241003321827614e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.6357456796751657e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.0401539881420724e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 7.3241003321827614e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.6357456796751657e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.0401539881420724e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 7.3241003321827614e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.6357456796751657e-01},
                                                                         {0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 3.2880511785370749e-02, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 6.7119488214629286e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 2.9592460606833626e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 3.2880511785370749e-02, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 6.7119488214629286e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 2.9592460606833626e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 3.2880511785370749e-02, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 6.7119488214629286e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 2.9592460606833626e-01},
                                                                         {0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.6536908390020822e-03, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 5.3460608600939230e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 4.6374022315160557e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.6536908390020822e-03, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 5.3460608600939230e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 4.6374022315160557e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.6536908390020822e-03, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 5.3460608600939230e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 4.6374022315160557e-01},
                                                                         {0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 3.7539852971250048e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 6.1013467358138296e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.4466796706116515e-02, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 3.7539852971250048e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 6.1013467358138296e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.4466796706116515e-02, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 3.7539852971250048e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 6.1013467358138296e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.4466796706116515e-02},
                                                                         {0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.1392940120350591e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 6.1392940120350581e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 2.7214119759298822e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.1392940120350591e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 6.1392940120350581e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 2.7214119759298822e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.1392940120350591e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 6.1392940120350581e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 2.7214119759298822e-01},
                                                                         {0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 5.4324933189692651e-03, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.9023366204009765e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 8.0433384464093305e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 5.4324933189692651e-03, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.9023366204009765e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 8.0433384464093305e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 5.4324933189692651e-03, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.9023366204009765e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 8.0433384464093305e-01}};
        static double expectedBOperatorOmegaTMaster[][27] = {{0000000000000000, 0000000000000000, 0000000000000000, 1.2688646455931618e+00, -4.8438025276657193e+00, -9.1906948398051924e-01, -1.1893664609088134e+01, 4.5403237371924675e+01, 8.6148701777422314e+00, 0000000000000000, 0000000000000000, 0000000000000000, 3.0010022551289456e-01, -1.1456117371864183e+00, -2.1736988288111531e-01, 9.7029687343121243e+00, -3.7040408245555263e+01, -7.0280959428533514e+00, 0000000000000000, 0000000000000000, 0000000000000000, 8.5699473617678289e-03, -3.2715177964067060e-02, -6.2074210412239568e-03, 6.1316105630818563e-01, -2.3406996835532117e+00, -4.4412744698602180e-01},
                                                             {0000000000000000, 0000000000000000, 0000000000000000, 3.7226397778604853e-01, -1.4129958263095190e+00, -1.7862105984706447e-01, -4.6783240159795509e+00, 1.7757432099707270e+01, 2.2447704959584200e+00, 0000000000000000, 0000000000000000, 0000000000000000, 8.3979861554674196e-01, -3.1876088193256567e+00, -4.0295523531226368e-01, 1.5311541472261383e+00, -5.8117748387422035e+00, -7.3468396860022667e-01, 0000000000000000, 0000000000000000, 0000000000000000, 1.5584487925356447e-01, -5.9153766433871258e-01, -7.4778058488399718e-02, 1.7792623961670582e+00, -6.7535149509911774e+00, -8.5373217371046473e-01},
                                                             {0000000000000000, 0000000000000000, 0000000000000000, 1.7907107964798094e-02, -6.7700630630254027e-02, -4.9417340313932219e-03, -6.1324871863390684e-01, 2.3184829770572515e+00, 1.6923514777142482e-01, 0000000000000000, 0000000000000000, 0000000000000000, 7.5522921174867497e-01, -2.8552625028165144e+00, -2.0841678648151182e-01, -2.2116895111800479e+00, 8.3616391300903050e+00, 6.1034877020674250e-01, 0000000000000000, 0000000000000000, 0000000000000000, 4.6467107666938218e-01, -1.7567618952204711e+00, -1.2823292725409702e-01, 1.5871308334310996e+00, -6.0003970784803169e+00, -4.3799247021116555e-01},
                                                             {0000000000000000, 0000000000000000, 0000000000000000, 5.5740850754883642e-01, -2.1050552552485393e+00, -1.0838955930278911e-01, -2.1350145528295359e+00, 8.0628902207276791e+00, 4.1515922945607331e-01, 0000000000000000, 0000000000000000, 0000000000000000, 6.4258816822883280e-01, -2.4267365534820504e+00, -1.2495296972374867e-01, 8.3756783057771333e-01, -3.1630779572022591e+00, -1.6286728102110576e-01, 0000000000000000, 0000000000000000, 0000000000000000, 1.9877105683239193e-03, -7.5065961877101980e-03, -3.8651557987433484e-04, 9.5462335905829518e-02, -3.6051385860712126e-01, -1.8562903828555462e-02},
                                                             {0000000000000000, 0000000000000000, 0000000000000000, 3.4775224449698239e-01, -1.3123015855458826e+00, -3.6287257594218789e-02, -1.0710668021091185e+00, 4.0418507281426068e+00, 1.1176369833347541e-01, 0000000000000000, 0000000000000000, 0000000000000000, 7.8874659955570126e-01, -2.9764679583537581e+00, -8.2304144653450914e-02, -3.1117746680939845e-01, 1.1742805101684064e+00, 3.2470752020486347e-02, 0000000000000000, 0000000000000000, 0000000000000000, 3.8639138277442461e-02, -1.4581128728287573e-01, -4.0319175104687518e-03, 2.0710628658839081e-01, -7.8155040712849710e-01, -2.1611130595823529e-02},
                                                             {0000000000000000, 0000000000000000, 0000000000000000, 1.9020160421184584e-01, -7.1761988658063092e-01, -3.0781640259431799e-03, -3.2449089674596487e-01, 1.2242857860437739e+00, 5.2514604661113799e-03, 0000000000000000, 0000000000000000, 0000000000000000, 8.5163338647259279e-01, -3.2131645615777518e+00, -1.3782571731689074e-02, -8.2489113970860339e-01, 3.1122675782467870e+00, 1.3349783468399143e-02, 0000000000000000, 0000000000000000, 0000000000000000, 1.2094725948550827e-01, -4.5632716397889872e-01, -1.9573730974952511e-03, -1.3400213715378848e-02, 5.0558247846720639e-02, 2.1686492061700861e-04},
                                                             {0000000000000000, 0000000000000000, 0000000000000000, 1.2087873116307407e-01, -4.5613636261462437e-01, 4.2105657675387620e-03, -3.8676087282129637e-02, 1.4594436592188070e-01, -1.3472031644076615e-03, 0000000000000000, 0000000000000000, 0000000000000000, 8.5115085377567801e-01, -3.2118210601814230e+00, 2.9648116036926478e-02, -8.3071065133613253e-01, 3.1346898767038018e+00, -2.8936122985333261e-02, 0000000000000000, 0000000000000000, 0000000000000000, 1.9009383660374632e-01, -7.1731983244360276e-01, 6.6215337745731346e-03, -2.9273668292423688e-01, 1.1046430126139672e+00, -1.0196889429297466e-02},
                                                             {0000000000000000, 0000000000000000, 0000000000000000, 3.8542781958847930e-02, -1.4554327064864342e-01, 4.7404604014314646e-03, 1.9292719476252929e-01, -7.2852174896940536e-01, 2.3728534388293619e-02, 0000000000000000, 0000000000000000, 0000000000000000, 7.8677966338617633e-01, -2.9709968940831901e+00, 9.6767738325581651e-02, -3.3695270892716001e-01, 1.2723834870959867e+00, -4.1442544949913347e-02, 0000000000000000, 0000000000000000, 0000000000000000, 3.4688503762963080e-01, -1.3098894358377886e+00, 4.2664143612883114e-02, -1.0281819688100242e+00, 3.8825678624430409e+00, -1.2645833177827653e-01},
                                                             {0000000000000000, 0000000000000000, 0000000000000000, 1.9790592447205151e-03, -7.4825860674192337e-03, 4.2170728506910266e-04, 9.2282126917440963e-02, -3.4890767367692993e-01, 1.9663911177278649e-02, 0000000000000000, 0000000000000000, 0000000000000000, 6.3979136356539168e-01, -2.4189745485591154e+00, 1.3632976357810986e-01, 7.9192945341286791e-01, -2.9941904520007854e+00, 1.6874806585800883e-01, 0000000000000000, 0000000000000000, 0000000000000000, 5.5498243936016878e-01, -2.0983221596305404e+00, 1.1825827770843518e-01, -2.0809644425005898e+00, 7.8678774199347910e+00, -4.4342172560690180e-01},
                                                             {0000000000000000, 0000000000000000, 0000000000000000, 4.6193687683245760e-01, -1.7491960093759249e+00, 1.3607734387121850e-01, 1.5417570020611673e+00, -5.8381032792298875e+00, 4.5417070655614550e-01, 0000000000000000, 0000000000000000, 0000000000000000, 7.5078532080886262e-01, -2.8429656797744753e+00, 2.2116630517512340e-01, -2.1730963989183962e+00, 8.2287683439395956e+00, -6.4015063696284558e-01, 0000000000000000, 0000000000000000, 0000000000000000, 1.7801739642697958e-02, -6.7409062806324388e-02, 5.2440356428696874e-03, -5.9918454042678937e-01, 2.2689056872470159e+00, -1.7650775428251159e-01},
                                                             {0000000000000000, 0000000000000000, 0000000000000000, 1.5445027340042236e-01, -5.8771779756869202e-01, 7.6973300698662539e-02, 1.7536160571327863e+00, -6.6729008902897569e+00, 8.7394870273707037e-01, 0000000000000000, 0000000000000000, 0000000000000000, 8.3228352701568620e-01, -3.1670247690803017e+00, 4.1478469918556105e-01, 1.5055378410261855e+00, -5.7289078523690922e+00, 7.5031409397431892e-01, 0000000000000000, 0000000000000000, 0000000000000000, 3.6893270681441931e-01, -1.4038713763742263e+00, 1.8386479708957348e-01, -4.6148204053894997e+00, 1.7560422685682074e+01, -2.2998855936851861e+00},
                                                             {0000000000000000, 0000000000000000, 0000000000000000, 8.4733769733794650e-03, -3.2454596155753544e-02, 6.2934948011267406e-03, 6.0698309717273435e-01, -2.3248571796107438e+00, 4.5082910608482585e-01, 0000000000000000, 0000000000000000, 0000000000000000, 2.9671854834384825e-01, -1.1364867500495868e+00, 2.2038399179767670e-01, 9.6053991139795976e+00, -3.6790449680029845e+01, 7.1342900919546723e+00, 0000000000000000, 0000000000000000, 0000000000000000, 1.2545664537298287e+00, -4.8052208386662141e+00, 9.3181354718692022e-01, -1.1772140590199387e+01, 4.5089469044512136e+01, -8.7436102318252242e+00}};
        static double expectedBOperatorOmegaNMaster[][27] = {{0000000000000000, 0000000000000000, 0000000000000000, -1.5579964493267173e-16, 5.9475430776267812e-16, 1.1284946725399597e-16, -1.2103084210714748e+01, 4.6202682134787885e+01, 8.7665578820776311e+00, 0000000000000000, 0000000000000000, 0000000000000000, -3.6848302725990667e-17, 1.4066583264356587e-16, 2.6690120722925999e-17, 1.1394848505601823e+01, -4.3499041592415146e+01, -8.2535655575650129e+00, 0000000000000000, 0000000000000000, 0000000000000000, -1.0522751663798979e-18, 4.0169872558219883e-18, 7.6218846315020998e-19, 7.0823570511292555e-01, -2.7036405423727485e+00, -5.1299232451261834e-01},
                                                             {0000000000000000, 0000000000000000, 0000000000000000, -2.3685548623303610e-17, 8.9902819895762288e-17, 1.1364886345745770e-17, -7.2854278117928599e+00, 2.7653170076152357e+01, 3.4957205500276340e+00, 0000000000000000, 0000000000000000, 0000000000000000, -5.3432757745224083e-17, 2.0281377782300913e-16, 2.5638300744985648e-17, 4.0114176405211657e+00, -1.5226067312100833e+01, -1.9247730459993873e+00, 0000000000000000, 0000000000000000, 0000000000000000, -9.9157363739734910e-18, 3.7636985975748944e-17, 4.7578047997466255e-18, 3.2740101712716951e+00, -1.2427102764051519e+01, -1.5709475040282452e+00},
                                                             {0000000000000000, 0000000000000000, 0000000000000000, 2.8678105249163121e-19, -1.0842207543874786e-18, -7.9141516845857649e-20, -1.7120997542490752e+00, 6.4728616866789652e+00, 4.7247951133956445e-01, 0000000000000000, 0000000000000000, 0000000000000000, 1.2094941776387064e-17, -4.5726824638980257e-17, -3.3377799196602686e-18, -4.1885285289019309e+00, 1.5835389130221742e+01, 1.1558870373387753e+00, 0000000000000000, 0000000000000000, 0000000000000000, 7.4416740375735711e-18, -2.8134416025128377e-17, -2.0536411526810919e-18, 5.9006282831510068e+00, -2.2308250816900706e+01, -1.6283665486783405e+00},
                                                             {0000000000000000, 0000000000000000, 0000000000000000, 3.2248411285720333e-17, -1.2178624246145775e-16, -6.2707889099953762e-18, -6.3161617606246248e+00, 2.3853007851764502e+01, 1.2281943681300036e+00, 0000000000000000, 0000000000000000, 0000000000000000, 3.7176410578135074e-17, -1.4039684970528393e-16, -7.2290514128375885e-18, 5.9260071596376411e+00, -2.2379587582707909e+01, -1.1523277735441120e+00, 0000000000000000, 0000000000000000, 0000000000000000, 1.1499736199965803e-19, -4.3428795567116404e-19, -2.2361541345929639e-20, 3.9015460098698096e-01, -1.4734202690565938e+00, -7.5866594585891442e-02},
                                                             {0000000000000000, 0000000000000000, 0000000000000000, -1.1041673667446503e-17, 4.1667612762154299e-17, 1.1521767665985850e-18, -5.1213235367315999e+00, 1.9326175757881106e+01, 5.3439996244892007e-01, 0000000000000000, 0000000000000000, 0000000000000000, -2.5043929108781727e-17, 9.4507478809495154e-17, 2.6132843745012915e-18, 3.3813587936407128e+00, -1.2760126142716311e+01, -3.5283808948753581e-01, 0000000000000000, 0000000000000000, 0000000000000000, -1.2268526297162772e-18, 4.6297347513504741e-18, 1.2801964073317602e-19, 1.7399647430908869e+00, -6.5660496151647960e+00, -1.8156187296138537e-01},
                                                             {0000000000000000, 0000000000000000, 0000000000000000, 6.6881722687122747e-19, -2.5234095394691235e-18, -1.0823931460604746e-20, -3.8488490625884428e+00, 1.4521489654127043e+01, 6.2288584656464259e-02, 0000000000000000, 0000000000000000, 0000000000000000, 2.9946491892736501e-18, -1.1298641882966271e-17, -4.8464477694284113e-20, 7.6645627843291442e-01, -2.8917961542817774e+00, -1.2404091718925198e-02, 0000000000000000, 0000000000000000, 0000000000000000, 4.2529405060472377e-19, -1.6046103797234480e-18, -6.8828275788937107e-21, 3.0823927841555281e+00, -1.1629693499845263e+01, -4.9884492937538882e-02},
                                                             {0000000000000000, 0000000000000000, 0000000000000000, 4.3446570403439590e-19, -1.6394580255123475e-18, 1.5133732816147247e-20, -3.0825519978432738e+00, 1.1632021963978039e+01, -1.0737445532305522e-01, 0000000000000000, 0000000000000000, 0000000000000000, 3.0592301173830657e-18, -1.1544016757271320e-17, 1.0656208485427608e-19, -7.6649586791391922e-01, 2.8923751414773142e+00, -2.6699331067963027e-02, 0000000000000000, 0000000000000000, 0000000000000000, 6.8324056480397063e-19, -2.5782109310860106e-18, 2.3799301212687278e-20, 3.8490478657571927e+00, -1.4524397105455353e+01, 1.3407378639101877e-01},
                                                             {0000000000000000, 0000000000000000, 0000000000000000, -2.0432166764180464e-18, 7.7154897134109783e-18, -2.5129965336818855e-19, -1.7403431656485449e+00, 6.5717943414138071e+00, -2.1404858296080834e-01, 0000000000000000, 0000000000000000, 0000000000000000, -4.1708492412758465e-17, 1.5749746362518257e-16, -5.1298180006100436e-18, -3.3820942007506085e+00, 1.2771290150893057e+01, -4.1597110581398439e-01, 0000000000000000, 0000000000000000, 0000000000000000, -1.8388950087762390e-17, 6.9439407420698681e-17, -2.2616968803136933e-18, 5.1224373663991534e+00, -1.9343084492306861e+01, 6.3001968877479264e-01},
                                                             {0000000000000000, 0000000000000000, 0000000000000000, -1.6992346639543373e-20, 6.4246028286920262e-20, -3.6208094262114350e-21, -3.9028739415651276e-01, 1.4756299113305709e+00, -8.3164280112130201e-02, 0000000000000000, 0000000000000000, 0000000000000000, -5.4932951884543235e-18, 2.0769491439430952e-17, -1.1705372672324342e-18, -5.9280241377058553e+00, 2.2413149549945324e+01, -1.2631713636693487e+00, 0000000000000000, 0000000000000000, 0000000000000000, -4.7651195959012963e-18, 1.8016346702603487e-17, -1.0153741749660945e-18, 6.3183115318623688e+00, -2.3888779461275895e+01, 1.3463356437814791e+00},
                                                             {0000000000000000, 0000000000000000, 0000000000000000, -1.7509209836793847e-17, 6.6301353084988214e-17, -5.1578622261389156e-18, -5.9029302982988838e+00, 2.2352365960065178e+01, -1.7388849350098299e+00, 0000000000000000, 0000000000000000, 0000000000000000, -2.8457692779515483e-17, 1.0775949083629630e-16, -8.3830658264248929e-18, 4.1901626016918367e+00, -1.5866704021930568e+01, 1.2343379059419486e+00, 0000000000000000, 0000000000000000, 0000000000000000, -6.7475538433148349e-19, 2.5550664707066265e-18, -1.9876940999437416e-19, 1.7127676966070493e+00, -6.4856619381346086e+00, 5.0454702906788129e-01},
                                                             {0000000000000000, 0000000000000000, 0000000000000000, -2.6153017352416666e-17, 9.9518074133082592e-17, -1.3033865363422543e-17, -3.2747665062032163e+00, 1.2461218204436092e+01, -1.6320436438877775e+00, 0000000000000000, 0000000000000000, 0000000000000000, -1.4093031397710857e-16, 5.3627133133364853e-16, -7.0235365703711584e-17, -4.0123443252679198e+00, 1.5267866595613631e+01, -1.9996299097167893e+00, 0000000000000000, 0000000000000000, 0000000000000000, -6.2471261919858614e-17, 2.3771710893440089e-16, -3.1133769613443975e-17, 7.2871108314711366e+00, -2.7729084800049730e+01, 3.6316735536045655e+00},
                                                             {0000000000000000, 0000000000000000, 0000000000000000, -9.2843029786041825e-19, 3.5560592276833717e-18, -6.8957999521915299e-19, -7.0799180976672949e-01, 2.7117391731476412e+00, -5.2585206441371468e-01, 0000000000000000, 0000000000000000, 0000000000000000, -3.2511534785371210e-17, 1.2452517280629470e-16, -2.4147535957765746e-17, -1.1390924458139832e+01, 4.3629340968903001e+01, -8.4604667162280869e+00, 0000000000000000, 0000000000000000, 0000000000000000, -1.3746319914497103e-16, 5.2650939861921905e-16, -1.0209907241033380e-16, 1.2098916267906560e+01, -4.6341080142050643e+01, 8.9863187806418043e+00}};

        // Expected values for the displacement B-operator matrix of the slave patch by MATLAB
        static double expectedBOperatorDisplacementMatrixSlave[][81] = {{8.7388487398058579e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.2365376195435124e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 2.4613640650629183e-03, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 8.7388487398058579e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.2365376195435124e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 2.4613640650629183e-03, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 8.7388487398058579e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.2365376195435124e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 2.4613640650629183e-03, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000},
                                                                        {4.8380657349864603e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 4.6457408385121851e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 5.1619342650135389e-02, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 4.8380657349864603e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 4.6457408385121851e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 5.1619342650135389e-02, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 4.8380657349864603e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 4.6457408385121851e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 5.1619342650135389e-02, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000},
                                                                        {1.9003225191245110e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 6.3988149176664488e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.7008625632090402e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.9003225191245110e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 6.3988149176664488e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.7008625632090402e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.9003225191245110e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 6.3988149176664488e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.7008625632090402e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000},
                                                                        {1.0097092424061095e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 6.5472312570339131e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 2.4430595005599781e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.0097092424061095e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 6.5472312570339131e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 2.4430595005599781e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.0097092424061095e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 6.5472312570339131e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 2.4430595005599781e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000},
                                                                        {3.2389326261654590e-02, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 6.2203543311750786e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 3.4557524062083761e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 3.2389326261654590e-02, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 6.2203543311750786e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 3.4557524062083761e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 3.2389326261654590e-02, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 6.2203543311750786e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 3.4557524062083761e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000},
                                                                        {1.6530360257260390e-03, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 5.3532065494549308e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 4.6302630902878095e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.6530360257260390e-03, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 5.3532065494549308e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 4.6302630902878095e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.6530360257260390e-03, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 5.3532065494549308e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 4.6302630902878095e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000},
                                                                        {4.6302630902878100e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 5.3532065494549308e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.6530360257260293e-03, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 4.6302630902878100e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 5.3532065494549308e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.6530360257260293e-03, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 4.6302630902878100e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 5.3532065494549308e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.6530360257260293e-03, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000},
                                                                        {3.4557524062083772e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 6.2203543311750786e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 3.2389326261654500e-02, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 3.4557524062083772e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 6.2203543311750786e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 3.2389326261654500e-02, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 3.4557524062083772e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 6.2203543311750786e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 3.2389326261654500e-02, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000},
                                                                        {2.4430595005599781e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 6.5472312570339131e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.0097092424061095e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 2.4430595005599781e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 6.5472312570339131e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.0097092424061095e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 2.4430595005599781e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 6.5472312570339131e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.0097092424061095e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000},
                                                                        {1.7008625632090402e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 6.3988149176664488e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.9003225191245110e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.7008625632090402e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 6.3988149176664488e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.9003225191245110e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.7008625632090402e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 6.3988149176664488e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.9003225191245110e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000},
                                                                        {5.1619342650135444e-02, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 4.6457408385121879e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 4.8380657349864581e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 5.1619342650135444e-02, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 4.6457408385121879e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 4.8380657349864581e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 5.1619342650135444e-02, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 4.6457408385121879e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 4.8380657349864581e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000},
                                                                        {2.4613640650629235e-03, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.2365376195435138e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 8.7388487398058567e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 2.4613640650629235e-03, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.2365376195435138e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 8.7388487398058567e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 2.4613640650629235e-03, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.2365376195435138e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 8.7388487398058567e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000}};
        static double expectedBOperatorOmegaTSlave[][27] = {{5.6587802797774680e+00, -2.4779335555319317e+01, -3.9462611416764992e+00, 5.9792537652015554e-01, -2.6182662710516524e+00, -4.1697495967745368e-01, 0000000000000000, 0000000000000000, 0000000000000000, -6.0899659623403490e+00, 2.6667462357672417e+01, 4.2469569135243637e+00, 8.4605792337277261e-02, -3.7048183788671518e-01, -5.9001504591812547e-02, 0000000000000000, 0000000000000000, 0000000000000000, -2.5302958917717866e-01, 1.1079958552291531e+00, 1.7645513451593769e-01, 1.6841028826282457e-03, -7.3745486438937888e-03, -1.1744420945360642e-03, 0000000000000000, 0000000000000000, 0000000000000000},
                                                            {2.9678591760692901e+00, -1.2976584683297375e+01, -1.3674040504876139e+00, 2.9353252483382736e-01, -1.2834334245107459e+00, -1.3524144495940224e-01, 0000000000000000, 0000000000000000, 0000000000000000, -2.7677093141050286e+00, 1.2101455009331712e+01, 1.2751875012115323e+00, 2.8186389204898232e-01, -1.2324138199783690e+00, -1.2986526813056032e-01, 0000000000000000, 0000000000000000, 0000000000000000, -8.0686448907473585e-01, 3.5279117873412638e+00, 3.7175273660277419e-01, 3.1318210227664706e-02, -1.3693486888648546e-01, -1.4429474236728922e-02, 0000000000000000, 0000000000000000, 0000000000000000},
                                                            {1.1305757139261579e+00, -4.9371725246398830e+00, -2.9828526538833527e-01, 1.0624639588219539e-01, -4.6397315998408928e-01, -2.8031501120980819e-02, 0000000000000000, 0000000000000000, 0000000000000000, -7.7659222503510761e-01, 3.3913427902827840e+00, 2.0489208735846293e-01, 3.5775559994547562e-01, -1.5623023711105686e+00, -9.4388392355710157e-02, 0000000000000000, 0000000000000000, 0000000000000000, -9.1308014123795633e-01, 3.9873792887866508e+00, 2.4090235522940953e-01, 9.5094656519235266e-02, -4.1527402333489372e-01, -2.5089283722846113e-02, 0000000000000000, 0000000000000000, 0000000000000000},
                                                            {5.9786328029123359e-01, -2.6097351347443949e+00, -1.1071899198891479e-01, 5.5140726738686793e-02, -2.4069498273115583e-01, -1.0211574925741630e-02, 0000000000000000, 0000000000000000, 0000000000000000, -3.1331114503552371e-01, 1.3676355953954922e+00, 5.8022453127289784e-02, 3.5754757357553463e-01, -1.5607321872120024e+00, -6.6214648464581416e-02, 0000000000000000, 0000000000000000, 0000000000000000, -8.3065713680336273e-01, 3.6259044271558740e+00, 1.5383035537901879e-01, 1.3341670123343150e-01, -5.8237771786381354e-01, -2.4707593127070888e-02, 0000000000000000, 0000000000000000, 0000000000000000},
                                                            {1.9307226370034433e-01, -8.4246104526887777e-01, -1.8978261558098803e-02, 1.7376188036377679e-02, -7.5820116547839517e-02, -1.7080125084610734e-03, 0000000000000000, 0000000000000000, 0000000000000000, -8.2243380151857726e-02, 3.5886482439919898e-01, 8.0842082131828072e-03, 3.3370884481581969e-01, -1.4561216449781744e+00, -3.2802297024885081e-02, 0000000000000000, 0000000000000000, 0000000000000000, -6.4730771907613927e-01, 2.8244944518280120e+00, 6.3627861225420634e-02, 1.8539380267545538e-01, -8.0895646943231891e-01, -1.8223498347158376e-02, 0000000000000000, 0000000000000000, 0000000000000000},
                                                            {1.0877270901224780e-02, -4.7450356851592532e-02, -1.4457163994537997e-04, 8.7977436213344199e-04, -3.8378751261415643e-03, -1.1693229254886969e-05, 0000000000000000, 0000000000000000, 0000000000000000, -1.6728892052695293e-01, 7.2977119430093507e-01, 2.2234652244020906e-03, 2.8490691092753045e-01, -1.2428608899934612e+00, -3.7867457488741456e-03, 0000000000000000, 0000000000000000, 0000000000000000, -3.7580564536458011e-01, 1.6393920994819080e+00, 4.9948961411792900e-03, 2.4643060970064429e-01, -1.0750141718116479e+00, -3.2753507475069350e-03, 0000000000000000, 0000000000000000, 0000000000000000},
                                                            {-3.4596828728649398e-01, 1.5090993827419739e+00, -1.2419616155544624e-02, 2.4640245482486692e-01, -1.0747973330121754e+00, 8.8453884970521070e-03, 0000000000000000, 0000000000000000, 0000000000000000, -1.9435664108438921e-01, 8.4777564265412664e-01, -6.9770408683446378e-03, 2.8487436010644929e-01, -1.2426101951929980e+00, 1.0226458132414537e-02, 0000000000000000, 0000000000000000, 0000000000000000, 8.1684395922918766e-03, -3.5630396194333509e-02, 2.9323174422056037e-04, 8.7967384727488044e-04, -3.8371009965935793e-03, 3.1578650202057710e-05, 0000000000000000, 0000000000000000, 0000000000000000},
                                                            {-6.2143722800588552e-01, 2.7105633696728857e+00, -7.5142557606116628e-02, 1.8530139575822607e-01, -8.0824120772943031e-01, 2.2406158140761488e-02, 0000000000000000, 0000000000000000, 0000000000000000, -9.6110610404905233e-02, 4.1921193044129912e-01, -1.1621442606655600e-02, 3.3354251236480686e-01, -1.4548341739129742e+00, 4.0331084653370666e-02, 0000000000000000, 0000000000000000, 0000000000000000, 1.8133640316639665e-01, -7.9094684042072561e-01, 2.1926721649331109e-02, 1.7367527121360693e-02, -7.5753078051054301e-02, 2.1000357693090757e-03, 0000000000000000, 0000000000000000, 0000000000000000},
                                                            {-8.0941253728408791e-01, 3.5307787744397636e+00, -1.6823192285971544e-01, 1.3330302107661704e-01, -5.8148775279079012e-01, 2.7706296265160808e-02, 0000000000000000, 0000000000000000, 0000000000000000, -3.1434714149703163e-01, 1.3712293347069597e+00, -6.5335316199753551e-02, 3.5724291858214230e-01, -1.5583471420903297e+00, 7.4250966414172004e-02, 0000000000000000, 0000000000000000, 0000000000000000, 5.7811999605516717e-01, -2.5218460514583785e+00, 1.2015904634533250e-01, 5.5093743067193141e-02, -2.4032716280722671e-01, 1.1450930034803682e-02, 0000000000000000, 0000000000000000, 0000000000000000},
                                                            {-8.9608500070421171e-01, 3.9095625423127767e+00, -2.5675569756150424e-01, 9.4988081199038663e-02, -4.1442702860786851e-01, 2.7216983912376019e-02, 0000000000000000, 0000000000000000, 0000000000000000, -7.6754375336386715e-01, 3.3487451585277124e+00, -2.1992470764385239e-01, 3.5735465294159635e-01, -1.5591158923134525e+00, 1.0239301307438747e-01, 0000000000000000, 0000000000000000, 0000000000000000, 1.1051586974409082e+00, -4.8217379416357051e+00, 3.1666169175313391e-01, 1.0612732248653525e-01, -4.6302683828346403e-01, 3.0408716465459090e-02, 0000000000000000, 0000000000000000, 0000000000000000},
                                                            {-7.9892287727339262e-01, 3.4882145217382416e+00, -3.8634248640789354e-01, 3.1269493402318049e-02, -1.3652719689993109e-01, 1.5121276625594400e-02, 0000000000000000, 0000000000000000, 0000000000000000, -2.7419718993809723e+00, 1.1971851688940786e+01, -1.3259605794276330e+00, 2.8142544062086228e-01, -1.2287447720993794e+00, 1.3609148963034956e-01, 0000000000000000, 0000000000000000, 0000000000000000, 2.9351239202865154e+00, -1.2815181756700563e+01, 1.4193648793095366e+00, 2.9307592234466989e-01, -1.2796124849791608e+00, 1.4172542027004673e-01, 0000000000000000, 0000000000000000, 0000000000000000},
                                                            {-2.5164668515255922e-01, 1.0999180925088297e+00, -1.8129261215199699e-01, 1.6811560922099387e-03, -7.3481357444948178e-03, 1.2111472050872103e-03, 0000000000000000, 0000000000000000, 0000000000000000, -6.0568468096533614e+00, 2.6473765730127930e+01, -4.3635050422414832e+00, 8.4457751774693682e-02, -3.6915490928595901e-01, 6.0845492268005175e-02, 0000000000000000, 0000000000000000, 0000000000000000, 5.6254754414945083e+00, -2.4588292165711366e+01, 4.0527342403384576e+00, 5.9687914544450715e-01, -2.6088886118949373e+00, 4.3000677458192993e-01, 0000000000000000, 0000000000000000, 0000000000000000}};
        static double expectedBOperatorOmegaNSlave[][27] = {{7.0497249908647284e+00, -3.0870168567178773e+01, -4.9162636496727901e+00, 4.2689018465969178e-17, -1.8693171687116133e-16, -2.9770022234400921e-17, 0000000000000000, 0000000000000000, 0000000000000000, -6.7665221915947571e+00, 2.9630046695263236e+01, 4.7187666367622558e+00, 6.0404498173899364e-18, -2.6450635212869339e-17, -4.2124258610638641e-18, 0000000000000000, 0000000000000000, 0000000000000000, -2.8320279926996966e-01, 1.2401218719155296e+00, 1.9749701291053492e-01, 1.2023690894926529e-19, -5.2650758037655870e-19, -8.3849560881066792e-20, 0000000000000000, 0000000000000000, 0000000000000000},
                                                            {5.4705692817934786e+00, -2.3919364545139388e+01, -2.5204964759496642e+00, 5.0033799425791854e-17, -2.1876638908986604e-16, -2.3052448225230534e-17, 0000000000000000, 0000000000000000, 0000000000000000, -4.1696891549153801e+00, 1.8231432561921135e+01, 1.9211322038726191e+00, 4.8044834037166976e-17, -2.1006989229994031e-16, -2.2136057262136594e-17, 0000000000000000, 0000000000000000, 0000000000000000, -1.3008801268780992e+00, 5.6879319832182524e+00, 5.9936427207704701e-01, 5.3383148930185535e-18, -2.3341099144437814e-17, -2.4595619180151768e-18, 0000000000000000, 0000000000000000, 0000000000000000},
                                                            {3.5312079003022601e+00, -1.5420623678196490e+01, -9.3165568011824562e-01, 5.4063821456368297e-18, -2.3609423993769913e-17, -1.4263919817468919e-18, 0000000000000000, 0000000000000000, 0000000000000000, -1.1958059188330781e+00, 5.2220298512886938e+00, 3.1549526622448193e-01, 1.8204509169340548e-17, -7.9498260389213804e-17, -4.8029838090045790e-18, 0000000000000000, 0000000000000000, 0000000000000000, -2.3354019814691815e+00, 1.0198593826907796e+01, 6.1616041389376397e-01, 4.8389222889132905e-18, -2.1131352707661844e-17, -1.2766762998380699e-18, 0000000000000000, 0000000000000000, 0000000000000000},
                                                            {2.5966665750211138e+00, -1.1334718517497750e+01, -4.8087968469580072e-01, 1.6882524615488437e-18, -7.3693968344677620e-18, -3.1264942492276502e-19, 0000000000000000, 0000000000000000, 0000000000000000, 1.7590589667376003e-01, -7.6784745625221540e-01, -3.2576216346884394e-02, 1.0947091322722898e-17, -4.7785187332877832e-17, -2.0273044965599090e-18, 0000000000000000, 0000000000000000, 0000000000000000, -2.7725724716948736e+00, 1.2102565973749964e+01, 5.1345590104268424e-01, 4.0848405088400480e-18, -1.7830751857772789e-17, -7.5647633578356664e-19, 0000000000000000, 0000000000000000, 0000000000000000},
                                                            {1.4807655901155761e+00, -6.4612456649034993e+00, -1.4555356702639527e-01, -1.1293394635856504e-18, 4.9278155583880694e-18, 1.1100972929532837e-19, 0000000000000000, 0000000000000000, 0000000000000000, 1.7693710935096614e+00, -7.7205611636697684e+00, -1.7392237891861290e-01, -2.1688909386172196e-17, 9.4638457756809417e-17, 2.1319364436495502e-18, 0000000000000000, 0000000000000000, 0000000000000000, -3.2501366836252390e+00, 1.4181806828573270e+01, 3.1947594594500855e-01, -1.2049394103428997e-17, 5.2576920976005221e-17, 1.1844091353608612e-18, 0000000000000000, 0000000000000000, 0000000000000000},
                                                            {3.3563211523848591e-01, -1.4641414913301214e+00, -4.4609429845956098e-03, 1.2814452145145871e-21, -5.5901000597159035e-21, -1.7031904219805061e-23, 0000000000000000, 0000000000000000, 0000000000000000, 3.3576782281579787e+00, -1.4647334939592177e+01, -4.4627466968670579e-02, 4.1498435656259928e-19, -1.8103029689650029e-18, -5.5156269918015802e-21, 0000000000000000, 0000000000000000, 0000000000000000, -3.6933103433964640e+00, 1.6111476430922295e+01, 4.9088409953267305e-02, 3.5894126846913590e-19, -1.5658239490666718e-18, -4.7707488668701621e-21, 0000000000000000, 0000000000000000, 0000000000000000},
                                                            {3.6943360860065435e+00, -1.6114541453382277e+01, 1.3261977419272422e-01, -3.2306501967185842e-19, 1.4091962751736232e-18, -1.1597431571193001e-20, 0000000000000000, 0000000000000000, 0000000000000000, -3.3586107557035505e+00, 1.4650121426030044e+01, -1.2056786108600687e-01, -3.7350659033495213e-19, 1.6292203234736769e-18, -1.3408189865925468e-20, 0000000000000000, 0000000000000000, 0000000000000000, -3.3572533030299384e-01, 1.4644200273522352e+00, -1.2051913106717256e-02, -1.1533645189398414e-21, 5.0309284046236240e-21, -4.1403634781113095e-23, 0000000000000000, 0000000000000000, 0000000000000000},
                                                            {3.2540899684145908e+00, -1.4193576877117911e+01, 3.9347601316341446e-01, -4.3120523756362058e-18, 1.8808160648849560e-17, -5.2140204904777149e-19, 0000000000000000, 0000000000000000, 0000000000000000, -1.7715232577143063e+00, 7.7269687660850925e+00, -2.1420793998861035e-01, -7.7616942761451686e-18, 3.3854689167929199e-17, -9.3852368828598833e-19, 0000000000000000, 0000000000000000, 0000000000000000, -1.4825667107002847e+00, 6.4666081110328184e+00, -1.7926807317480500e-01, -4.0415068799752980e-19, 1.7628104679686560e-18, -4.8868839821289951e-20, 0000000000000000, 0000000000000000, 0000000000000000},
                                                            {2.7783813599920064e+00, -1.2119715820161304e+01, 5.7747121164864512e-01, 1.2713733584061103e-17, -5.5459210989129559e-17, 2.6424792661965459e-18, 0000000000000000, 0000000000000000, 0000000000000000, -1.7627444166763048e-01, 7.6893552848174451e-01, -3.6637668564256480e-02, 3.4071930665662096e-17, -1.4862686708827186e-16, 7.0816624989037241e-18, 0000000000000000, 0000000000000000, 0000000000000000, -2.6021069183243766e+00, 1.1350780291679561e+01, -5.4083354308438902e-01, 5.2545483654299102e-18, -2.2921127340303937e-17, 1.0921288398148297e-18, 0000000000000000, 0000000000000000, 0000000000000000},
                                                            {2.3418957136285532e+00, -1.0217543818733411e+01, 6.7102458706088453e-01, -4.8416370151374567e-18, 2.1123758017354210e-17, -1.3872767518530950e-18, 0000000000000000, 0000000000000000, 0000000000000000, 1.1991309324337804e+00, -5.2317328970883965e+00, 3.4358760472796052e-01, -1.8214722240658789e-17, 7.9469688405394108e-17, -5.2190737609870277e-18, 0000000000000000, 0000000000000000, 0000000000000000, -3.5410266460623343e+00, 1.5449276715821808e+01, -1.0146121917888444e+00, -5.4094152274910795e-18, 2.3600938675008976e-17, -1.5499625359678979e-18, 0000000000000000, 0000000000000000, 0000000000000000},
                                                            {1.3060141958000924e+00, -5.7022496325727579e+00, 6.3156130089982265e-01, -5.5970415029053147e-18, 2.4437504550924159e-17, -2.7066128562252256e-18, 0000000000000000, 0000000000000000, 0000000000000000, 4.1861452995380155e+00, -1.8277324682189445e+01, 2.0243327979388877e+00, -5.0373373526147803e-17, 2.1993754095831733e-16, -2.4359515706027020e-17, 0000000000000000, 0000000000000000, 0000000000000000, -5.4921594953381083e+00, 2.3979574314762210e+01, -2.6558940988387105e+00, -5.2458736051789403e-17, 2.2904253976603392e-16, -2.5367953649333683e-17, 0000000000000000, 0000000000000000, 0000000000000000},
                                                            {2.8446035737997255e-01, -1.2433427982334677e+00, 2.0493240835595147e-01, -1.9255444505224971e-19, 8.4163285432338545e-19, -1.3872107357121489e-19, 0000000000000000, 0000000000000000, 0000000000000000, 6.7965688397228377e+00, -2.9707003806788936e+01, 4.8964194297939061e+00, -9.6735309699641800e-18, 4.2281867237224285e-17, -6.9690554326787709e-18, 0000000000000000, 0000000000000000, 0000000000000000, -7.0810291971028105e+00, 3.0950346605022400e+01, -5.1013518381498573e+00, -6.8364700426625057e-17, 2.9881407276477399e-16, -4.9251652616913288e-17, 0000000000000000, 0000000000000000, 0000000000000000}};

        // Loop over all Gauss Points on the interface
        for(int iGP = 0; iGP < noGPsOnContCond; iGP++){

            // Get the parametric coordinates of the Gauss Point on the master patch
            uGPMaster = trCurveMasterGPs[2*iGP];
            vGPMaster = trCurveMasterGPs[2*iGP + 1];

            // Get the parametric coordinates of the Gauss Point on the slave patch
            uGPSlave = trCurveSlaveGPs[2*iGP];
            vGPSlave = trCurveSlaveGPs[2*iGP + 1];

            // Find the knot span indices of the Gauss point locations in the parameter space of the master patch
            uKnotSpanMaster = patchMaster->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(uGPMaster);
            vKnotSpanMaster = patchMaster->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(vGPMaster);

            // Find the knot span indices of the Gauss point locations in the parameter space of the slave patch
            uKnotSpanSlave = patchSlave->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(uGPSlave);
            vKnotSpanSlave = patchSlave->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(vGPSlave);

            // Get the tangent to the boundary vector on the master and the slave patch
            for(int iCoord = 0; iCoord < noCoord; iCoord++){
                tangentTrCurveVctMaster[iCoord] = trCurveMasterGPTangents[3*iGP + iCoord];
                tangentTrCurveVctSlave[iCoord] = trCurveSlaveGPTangents[3*iGP + iCoord];
            }

            // Compute the B-operator matrices needed for the computation of the patch weak continuity contributions at the master patch
            theMapper->computeDisplacementAndRotationBOperatorMatrices(BDisplacementsGCMaster, BOperatorOmegaTMaster, BOperatorOmegaNMaster, normalTrCurveVctMaster,
                                                                       patchMaster, tangentTrCurveVctMaster, uGPMaster, vGPMaster, uKnotSpanMaster, vKnotSpanMaster);

            // Compare the computed values against the expected ones for the master patch
            for(int i = 0; i < noCoord*noDOFsLocMaster; i++){
                CPPUNIT_ASSERT(fabs( BDisplacementsGCMaster[i] - expectedBOperatorDisplacementMatrixMaster[iGP][i] ) <= TolRel10);
            }
            for(int i = 0; i < noDOFsLocMaster; i++){
                CPPUNIT_ASSERT(fabs( BOperatorOmegaTMaster[i] + expectedBOperatorOmegaTMaster[iGP][i] ) <= TolRel1000);
                CPPUNIT_ASSERT(fabs( BOperatorOmegaNMaster[i] - expectedBOperatorOmegaNMaster[iGP][i] ) <= TolRel1000);
            }

            // Compute the B-operator matrices needed for the computation of the patch weak continuity contributions at the master patch
            theMapper->computeDisplacementAndRotationBOperatorMatrices(BDisplacementsGCSlave, BOperatorOmegaTSlave, BOperatorOmegaNSlave, normalTrCurveVctSlave,
                                                                       patchSlave, tangentTrCurveVctSlave, uGPSlave, vGPSlave, uKnotSpanSlave, vKnotSpanSlave);

            // Compare the computed values against the expected ones for the slave patch
            for(int i = 0; i < noCoord*noDOFsLocSlave; i++){
                CPPUNIT_ASSERT(fabs( BDisplacementsGCSlave[i] - expectedBOperatorDisplacementMatrixSlave[iGP][i] ) <= TolRel10);
            }
            for(int i = 0; i < noDOFsLocSlave; i++){
                CPPUNIT_ASSERT(fabs( BOperatorOmegaTSlave[i] + expectedBOperatorOmegaTSlave[iGP][i] ) <= TolRel1000);
                CPPUNIT_ASSERT(fabs( BOperatorOmegaNSlave[i] - expectedBOperatorOmegaNSlave[iGP][i] ) <= TolRel1000);
            }
        }

        // Delete pointers
        delete[] BDisplacementsGCMaster;
        delete[] BDisplacementsGCSlave;
        delete[] BOperatorOmegaTMaster;
        delete[] BOperatorOmegaTSlave;
        delete[] BOperatorOmegaNMaster;
        delete[] BOperatorOmegaNSlave;
    }

    /***********************************************************************************************
     * \brief Test case: testComputePenaltyFactorsIGAPatchContinuityConditions with the precomputed GP data
     ***********/

    void testComputePenaltyFactorsIGAPatchContinuityConditions() {

        // Initialize empty string
        std::string emptyString = "";

        // Add the precomputed GP data to the continuity condition
        addPrecomputedGPData();

        // Set the computation of the Penalty parameters to automatic
        theMapper->setParametersWeakPatchContinuityConditions(true, true, 0.0, 0.0 ,0.0);

        // Compute the Penalty parameters
        theMapper->computePenaltyParametersForPatchContinuityConditions(emptyString);
        int noWeakPatchContinuityConditions = theMapper->getNoWeakIGAPatchContinuityConditions();
        double* weakPatchContinuityAlphaPrimaryIJ = new double[noWeakPatchContinuityConditions];
        double* weakPatchContinuityAlphaSecondaryBendingIJ = new double[noWeakPatchContinuityConditions];
        double* weakPatchContinuityAlphaSecondaryTwistingIJ = new double[noWeakPatchContinuityConditions];
        theMapper->getPenaltyParameterForPatchContinuityPrimaryField(weakPatchContinuityAlphaPrimaryIJ);
        theMapper->getPenaltyParameterForPatchContinuitySecondaryFieldBending(weakPatchContinuityAlphaSecondaryBendingIJ);
        theMapper->getPenaltyParameterForPatchContinuitySecondaryFieldTwisting(weakPatchContinuityAlphaSecondaryTwistingIJ);

        // Define the expected solution in terms of the Penalty parameters
        double expectedAlphaPrimaryIJ = 5.527399039185964e+01;
        double expectedAlphaSecondaryBendingIJ = 1.051417998627184e+01;
        double expectedAlphaSecondaryTwistingIJ = 1.051417998627184e+01;

        // Check if the computed Penalty parameters match the expected ones
        CPPUNIT_ASSERT(fabs( weakPatchContinuityAlphaPrimaryIJ[0] - expectedAlphaPrimaryIJ ) <= Tol);
        CPPUNIT_ASSERT(fabs( weakPatchContinuityAlphaSecondaryBendingIJ[0] - expectedAlphaSecondaryBendingIJ ) <= Tol);
        CPPUNIT_ASSERT(fabs( weakPatchContinuityAlphaSecondaryTwistingIJ[0] - expectedAlphaSecondaryTwistingIJ ) <= Tol);

        // Delete pointers
        delete[] weakPatchContinuityAlphaPrimaryIJ;
        delete[] weakPatchContinuityAlphaSecondaryBendingIJ;
        delete[] weakPatchContinuityAlphaSecondaryTwistingIJ;
    }

    /***********************************************************************************************
     * \brief Test case: testCreateWeakContinuityConditionGPData
     ***********/

    void testCreateWeakContinuityConditionGPData() {

        // Initialize auxiliary variables
        int noGPsConnection = 12;

        // Initialize the necessary for the coupling constituents for each Gauss Point
        double* referenceTrCurveMasterGPs = new double[noGPsConnection*2]; // The parametric coordinates of the Gauss Points on the master curve
        double* referenceTrCurveSlaveGPs = new double[noGPsConnection*2]; // The parametric coordinates of the Gauss Points on the slave
        double* referenceTrCurveGPWeights = new double[noGPsConnection]; // The Gauss Point weights for each Gauss Point
        double* referenceTrCurveMasterGPTangents = new double[noGPsConnection*3]; // The parametric components of the tangent to the trimming curve vectors on the master curve
        double* referenceTrCurveSlaveGPTangents = new double[noGPsConnection*3]; // The parametric components of the tangent to the trimming curve vectors on the slave curve
        double* referenceTrCurveGPJacobianProducts = new double[noGPsConnection]; // The products of the Jacobian transformations with the Gauss weights at each Gauss Point

        // Assign the necessary for the coupling constituents for each Gauss Point

        // The parametric coordinates of the Gauss Points on the master curve
        referenceTrCurveMasterGPs[2*0 + 0] = 5;
        referenceTrCurveMasterGPs[2*0 + 1] = 0.037567221793086;

        referenceTrCurveMasterGPs[2*1 + 0] = 5;
        referenceTrCurveMasterGPs[2*1 + 1] = 0.166666666666667;

        referenceTrCurveMasterGPs[2*2 + 0] = 5;
        referenceTrCurveMasterGPs[2*2 + 1] = 0.295766111540247;

        referenceTrCurveMasterGPs[2*3 + 0] = 5;
        referenceTrCurveMasterGPs[2*3 + 1] = 0.352116944229876;

        referenceTrCurveMasterGPs[2*4 + 0] = 5;
        referenceTrCurveMasterGPs[2*4 + 1] = 0.416666666666667;

        referenceTrCurveMasterGPs[2*5 + 0] = 5;
        referenceTrCurveMasterGPs[2*5 + 1] = 0.481216389103457;

        referenceTrCurveMasterGPs[2*6 + 0] = 5;
        referenceTrCurveMasterGPs[2*6 + 1] = 0.518783610896543;

        referenceTrCurveMasterGPs[2*7 + 0] = 5;
        referenceTrCurveMasterGPs[2*7 + 1] = 0.583333333333333;

        referenceTrCurveMasterGPs[2*8 + 0] = 5;
        referenceTrCurveMasterGPs[2*8 + 1] = 0.647883055770124;

        referenceTrCurveMasterGPs[2*9 + 0] = 5;
        referenceTrCurveMasterGPs[2*9 + 1] = 0.704233888459753;

        referenceTrCurveMasterGPs[2*10 + 0] = 5;
        referenceTrCurveMasterGPs[2*10 + 1] = 0.833333333333333;

        referenceTrCurveMasterGPs[2*11 + 0] = 5;
        referenceTrCurveMasterGPs[2*11 + 1] = 0.962432778206914;

        // The parametric coordinates of the Gauss Points on the slave
        referenceTrCurveSlaveGPs[2*0 + 0] = 0;
        referenceTrCurveSlaveGPs[2*0 + 1] = 0.037567221793086;

        referenceTrCurveSlaveGPs[2*1 + 0] = 0;
        referenceTrCurveSlaveGPs[2*1 + 1] = 0.166666666666667;

        referenceTrCurveSlaveGPs[2*2 + 0] = 0;
        referenceTrCurveSlaveGPs[2*2 + 1] = 0.295766111540247;

        referenceTrCurveSlaveGPs[2*3 + 0] = 0;
        referenceTrCurveSlaveGPs[2*3 + 1] = 0.352116944229876;

        referenceTrCurveSlaveGPs[2*4 + 0] = 0;
        referenceTrCurveSlaveGPs[2*4 + 1] = 0.416666666666667;

        referenceTrCurveSlaveGPs[2*5 + 0] = 0;
        referenceTrCurveSlaveGPs[2*5 + 1] = 0.481216389103457;

        referenceTrCurveSlaveGPs[2*6 + 0] = 0;
        referenceTrCurveSlaveGPs[2*6 + 1] = 0.518783610896543;

        referenceTrCurveSlaveGPs[2*7 + 0] = 0;
        referenceTrCurveSlaveGPs[2*7 + 1] = 0.583333333333333;

        referenceTrCurveSlaveGPs[2*8 + 0] = 0;
        referenceTrCurveSlaveGPs[2*8 + 1] = 0.647883055770124;

        referenceTrCurveSlaveGPs[2*9 + 0] = 0;
        referenceTrCurveSlaveGPs[2*9 + 1] = 0.704233888459753;

        referenceTrCurveSlaveGPs[2*10 + 0] = 0;
        referenceTrCurveSlaveGPs[2*10 + 1] = 0.833333333333333;

        referenceTrCurveSlaveGPs[2*11 + 0] = 0;
        referenceTrCurveSlaveGPs[2*11 + 1] = 0.962432778206914;

        // The Gauss Point weights for each Gauss Point
        referenceTrCurveGPWeights[0] = 0.555555555555556;
        referenceTrCurveGPWeights[1] = 0.888888888888889;
        referenceTrCurveGPWeights[2] = 0.555555555555556;
        referenceTrCurveGPWeights[3] = 0.555555555555556;
        referenceTrCurveGPWeights[4] = 0.888888888888889;
        referenceTrCurveGPWeights[5] = 0.555555555555556;
        referenceTrCurveGPWeights[6] = 0.555555555555556;
        referenceTrCurveGPWeights[7] = 0.888888888888889;
        referenceTrCurveGPWeights[8] = 0.555555555555556;
        referenceTrCurveGPWeights[9] = 0.555555555555556;
        referenceTrCurveGPWeights[10] = 0.888888888888889;
        referenceTrCurveGPWeights[11] = 0.555555555555556;

        // The parametric components of the tangent to the trimming curve vectors on the master curve
        referenceTrCurveMasterGPTangents[3*0 + 0] = 0.672567413298545;
        referenceTrCurveMasterGPTangents[3*0 + 1] = 0.035933526804658;
        referenceTrCurveMasterGPTangents[3*0 + 2] = 0.739162942943085;

        referenceTrCurveMasterGPTangents[3*1 + 0] = 0.518642424700041;
        referenceTrCurveMasterGPTangents[3*1 + 1] = 0.028618645620947;
        referenceTrCurveMasterGPTangents[3*1 + 2] = 0.854512146446196;

        referenceTrCurveMasterGPTangents[3*2 + 0] = 0.334518371397290;
        referenceTrCurveMasterGPTangents[3*2 + 1] = 0.019707912983376;
        referenceTrCurveMasterGPTangents[3*2 + 2] = 0.942183133665396;

        referenceTrCurveMasterGPTangents[3*3 + 0] = 0.246600757308154;
        referenceTrCurveMasterGPTangents[3*3 + 1] = 0.015405035206640;
        referenceTrCurveMasterGPTangents[3*3 + 2] = 0.968994711742705;

        referenceTrCurveMasterGPTangents[3*4 + 0] = 0.141959112266735;
        referenceTrCurveMasterGPTangents[3*4 + 1] = 0.010248227973116;
        referenceTrCurveMasterGPTangents[3*4 + 2] = 0.989819470543923;

        referenceTrCurveMasterGPTangents[3*5 + 0] = 0.034767353927938;
        referenceTrCurveMasterGPTangents[3*5 + 1] = 0.004928153317023;
        referenceTrCurveMasterGPTangents[3*5 + 2] = 0.999383282032341;

        referenceTrCurveMasterGPTangents[3*6 + 0] = -0.028052951297859;
        referenceTrCurveMasterGPTangents[3*6 + 1] = 0.001793096951330;
        referenceTrCurveMasterGPTangents[3*6 + 2] = 0.999604830283849;

        referenceTrCurveMasterGPTangents[3*7 + 0] = -0.135434083146395;
        referenceTrCurveMasterGPTangents[3*7 + 1] = -0.003595176569249;
        referenceTrCurveMasterGPTangents[3*7 + 2] = 0.990779836203650;

        referenceTrCurveMasterGPTangents[3*8 + 0] = -0.240485196250373;
        referenceTrCurveMasterGPTangents[3*8 + 1] = -0.008903380554228;
        referenceTrCurveMasterGPTangents[3*8 + 2] = 0.970611972004841;

        referenceTrCurveMasterGPTangents[3*9 + 0] = -0.328908196739009;
        referenceTrCurveMasterGPTangents[3*9 + 1] = -0.013401303233212;
        referenceTrCurveMasterGPTangents[3*9 + 2] = 0.944266807205222;

        referenceTrCurveMasterGPTangents[3*10 + 0] = -0.514499528145185;
        referenceTrCurveMasterGPTangents[3*10 + 1] = -0.022943557518762;
        referenceTrCurveMasterGPTangents[3*10 + 2] = 0.857183661012484;

        referenceTrCurveMasterGPTangents[3*11 + 0] = -0.669943646656621;
        referenceTrCurveMasterGPTangents[3*11 + 1] = -0.031071504836861;
        referenceTrCurveMasterGPTangents[3*11 + 2] = 0.741761465628677;

        // The parametric components of the tangent to the trimming curve vectors on the slave curve
        referenceTrCurveSlaveGPTangents[3*0 + 0] = -0.672674587247831;
        referenceTrCurveSlaveGPTangents[3*0 + 1] = -0.035915780365672;
        referenceTrCurveSlaveGPTangents[3*0 + 2] = -0.739066273342036;

        referenceTrCurveSlaveGPTangents[3*1 + 0] = -0.518754457663932;
        referenceTrCurveSlaveGPTangents[3*1 + 1] = -0.028606853848593;
        referenceTrCurveSlaveGPTangents[3*1 + 2] = -0.854444533347065;

        referenceTrCurveSlaveGPTangents[3*2 + 0] = -0.334609248379411;
        referenceTrCurveSlaveGPTangents[3*2 + 1] = -0.019701829572375;
        referenceTrCurveSlaveGPTangents[3*2 + 2] = -0.942150990452415;

        referenceTrCurveSlaveGPTangents[3*3 + 0] = -0.246673935138192;
        referenceTrCurveSlaveGPTangents[3*3 + 1] = -0.015401266353883;
        referenceTrCurveSlaveGPTangents[3*3 + 2] = -0.968976145587772;

        referenceTrCurveSlaveGPTangents[3*4 + 0] = -0.142006906143282;
        referenceTrCurveSlaveGPTangents[3*4 + 1] = -0.010246968609278;
        referenceTrCurveSlaveGPTangents[3*4 + 2] = -0.989812627845257;

        referenceTrCurveSlaveGPTangents[3*5 + 0] = -0.034786160987690;
        referenceTrCurveSlaveGPTangents[3*5 + 1] = -0.004929280359464;
        referenceTrCurveSlaveGPTangents[3*5 + 2] = -0.999382622021654;

        referenceTrCurveSlaveGPTangents[3*6 + 0] = 0.028051720822084;
        referenceTrCurveSlaveGPTangents[3*6 + 1] = -0.001795576153355;
        referenceTrCurveSlaveGPTangents[3*6 + 2] = -0.999604860364933;

        referenceTrCurveSlaveGPTangents[3*7 + 0] = 0.135462489104364;
        referenceTrCurveSlaveGPTangents[3*7 + 1] = 0.003590395052853;
        referenceTrCurveSlaveGPTangents[3*7 + 2] = -0.990775970191554;

        referenceTrCurveSlaveGPTangents[3*8 + 0] = 0.240540374596422;
        referenceTrCurveSlaveGPTangents[3*8 + 1] = 0.008896271275466;
        referenceTrCurveSlaveGPTangents[3*8 + 2] = -0.970598364178719;

        referenceTrCurveSlaveGPTangents[3*9 + 0] = 0.328982756872293;
        referenceTrCurveSlaveGPTangents[3*9 + 1] = 0.013392103616478;
        referenceTrCurveSlaveGPTangents[3*9 + 2] = -0.944240963547669;

        referenceTrCurveSlaveGPTangents[3*10 + 0] = 0.514599927063366;
        referenceTrCurveSlaveGPTangents[3*10 + 1] = 0.022929302160749;
        referenceTrCurveSlaveGPTangents[3*10 + 2] = -0.857123772957441;

        referenceTrCurveSlaveGPTangents[3*11 + 0] = 0.670043695210792;
        referenceTrCurveSlaveGPTangents[3*11 + 1] = 0.031051982197132;
        referenceTrCurveSlaveGPTangents[3*11 + 2] = -0.741671909209117;

        // The products of the Jacobian transformations with the Gauss weights at each Gauss Point
        referenceTrCurveGPJacobianProducts[0] = 0.009522090368294;
        referenceTrCurveGPJacobianProducts[1] = 0.016198303357023;
        referenceTrCurveGPJacobianProducts[2] = 0.010570156158216;
        referenceTrCurveGPJacobianProducts[3] = 0.005351640216454;
        referenceTrCurveGPJacobianProducts[4] = 0.008643208768130;
        referenceTrCurveGPJacobianProducts[5] = 0.005423354699369;
        referenceTrCurveGPJacobianProducts[6] = 0.005422152689563;
        referenceTrCurveGPJacobianProducts[7] = 0.008634830881354;
        referenceTrCurveGPJacobianProducts[8] = 0.005342723166447;
        referenceTrCurveGPJacobianProducts[9] = 0.010546831441419;
        referenceTrCurveGPJacobianProducts[10] = 0.016147819419794;
        referenceTrCurveGPJacobianProducts[11] = 0.009488732784896;

        // Add the Gauss Point data of the weak continuity condition
        theIGAMesh->createWeakContinuityConditionGPData(0);

        int computedNoGPsConnection = theIGAMesh->getWeakIGAPatchContinuityConditions()[0]->getTrCurveNumGP();
        double* computedTrCurveMasterGPs = theIGAMesh->getWeakIGAPatchContinuityConditions()[0]->getTrCurveMasterGPs();
        double* computedTrCurveSlaveGPs = theIGAMesh->getWeakIGAPatchContinuityConditions()[0]->getTrCurveSlaveGPs();
        double* computedTrCurveGPWeights = theIGAMesh->getWeakIGAPatchContinuityConditions()[0]->getTrCurveGPWeights();
        double* computedTrCurveMasterGPTangents = theIGAMesh->getWeakIGAPatchContinuityConditions()[0]->getTrCurveMasterGPTangents();
        double* computedTrCurveSlaveGPTangents = theIGAMesh->getWeakIGAPatchContinuityConditions()[0]->getTrCurveSlaveGPTangents();
        double* computedTrCurveGPJacobianProducts = theIGAMesh->getWeakIGAPatchContinuityConditions()[0]->getTrCurveGPJacobianProducts();

        // Compare number of computed gauss points
        CPPUNIT_ASSERT(computedNoGPsConnection == noGPsConnection);

        double tmpVect2[2];
        double tmpVect3[3];

        // Compare Master GP coordinates
        for (int i = 0; i<computedNoGPsConnection; i++) {
            tmpVect2[0] = referenceTrCurveMasterGPs[i*2] - computedTrCurveMasterGPs[i*2];
            tmpVect2[1] = referenceTrCurveMasterGPs[i*2+1] - computedTrCurveMasterGPs[i*2+1];
            CPPUNIT_ASSERT(EMPIRE::MathLibrary::vector2norm(tmpVect2,2) < TolRel1E9);
        }

        // Compare Slave GP coordinates
        for (int i = 0; i<computedNoGPsConnection; i++) {
            tmpVect2[0] = referenceTrCurveSlaveGPs[i*2] - computedTrCurveSlaveGPs[i*2];
            tmpVect2[1] = referenceTrCurveSlaveGPs[i*2+1] - computedTrCurveSlaveGPs[i*2+1];
            CPPUNIT_ASSERT(EMPIRE::MathLibrary::vector2norm(tmpVect2,2) < TolRel1E11);
        }

        // Compare Master GP tangents
        for (int i = 0; i<computedNoGPsConnection; i++) {
            tmpVect3[0] = referenceTrCurveMasterGPTangents[i*3] - computedTrCurveMasterGPTangents[i*3];
            tmpVect3[1] = referenceTrCurveMasterGPTangents[i*3+1] - computedTrCurveMasterGPTangents[i*3+1];
            tmpVect3[2] = referenceTrCurveMasterGPTangents[i*3+2] - computedTrCurveMasterGPTangents[i*3+2];
            CPPUNIT_ASSERT(EMPIRE::MathLibrary::vector2norm(tmpVect3,3) < TolRel1E9);
        }

        // Compare Slave GP tangents
        for (int i = 0; i<computedNoGPsConnection; i++) {
            tmpVect3[0] = referenceTrCurveSlaveGPTangents[i*3] - computedTrCurveSlaveGPTangents[i*3];
            tmpVect3[1] = referenceTrCurveSlaveGPTangents[i*3+1] - computedTrCurveSlaveGPTangents[i*3+1];
            tmpVect3[2] = referenceTrCurveSlaveGPTangents[i*3+2] - computedTrCurveSlaveGPTangents[i*3+2];
            CPPUNIT_ASSERT(EMPIRE::MathLibrary::vector2norm(tmpVect3,3) < TolRel1E11);
        }

        // Compare GP Jacobian products
        for (int i = 0; i<computedNoGPsConnection; i++) {
            CPPUNIT_ASSERT(fabs(referenceTrCurveGPJacobianProducts[i] - computedTrCurveGPJacobianProducts[i]) < TolRel1E8);
        }

        // Compare GP Weights
        for (int i = 0; i<computedNoGPsConnection; i++) {
            CPPUNIT_ASSERT(fabs(referenceTrCurveGPWeights[i] - computedTrCurveGPWeights[i]) < Tol);
        }

        delete referenceTrCurveMasterGPs;
        delete referenceTrCurveSlaveGPs;
        delete referenceTrCurveGPWeights;
        delete referenceTrCurveMasterGPTangents;
        delete referenceTrCurveSlaveGPTangents;
        delete referenceTrCurveGPJacobianProducts;

    }

    /***********************************************************************************************
     * \brief Test case: testComputeBOperatorMatricesIGAPatchContinuityConditions4Leakage
     ***********/

    void testComputeBOperatorMatricesIGAPatchContinuityConditions4Leakage(){

        // Add the precomputed GP data to the continuity condition
        addPrecomputedGPData();

        // Get the weak patch continuity conditions
        std::vector<WeakIGAPatchContinuityCondition*> weakIGAPatchContinuityConditions = theIGAMesh->getWeakIGAPatchContinuityConditions();

        // Initialize constant array sizes
        const int noCoord = 3;

        // Initialize varying array sizes
        int indexMaster;
        int indexSlave;
        int pMaster;
        int qMaster;
        int pSlave;
        int qSlave;
        int noLocalBasisFctsMaster;
        int noLocalBasisFctsSlave;
        int noDOFsLocMaster;
        int noDOFsLocSlave;
        int noGPsOnContCond;
        int uKnotSpanMaster;
        int uKnotSpanSlave;
        int vKnotSpanMaster;
        int vKnotSpanSlave;
        double uGPMaster;
        double vGPMaster;
        double uGPSlave;
        double vGPSlave;
        double tangentTrCurveVctMaster[noCoord];
        double tangentTrCurveVctSlave[noCoord];
        double normalTrCurveVctMaster[noCoord];
        double normalTrCurveVctSlave[noCoord];

        // Initialize pointers
        double* trCurveMasterGPs;
        double* trCurveSlaveGPs;
        double* trCurveGPWeights;
        double* trCurveMasterGPTangents;
        double* trCurveSlaveGPTangents;
        double* trCurveGPJacobianProducts;
        IGAPatchSurface* patchMaster;
        IGAPatchSurface* patchSlave;

        // Get the index of the master and slave patches
        indexMaster = weakIGAPatchContinuityConditions[0]->getMasterPatchIndex();
        indexSlave = weakIGAPatchContinuityConditions[0]->getSlavePatchIndex();

        // Get the number of Gauss Points for the given condition
        noGPsOnContCond = weakIGAPatchContinuityConditions[0]->getTrCurveNumGP();

        // Get the parametric coordinates of the Gauss Points
        trCurveMasterGPs = weakIGAPatchContinuityConditions[0]->getTrCurveMasterGPs();
        trCurveSlaveGPs = weakIGAPatchContinuityConditions[0]->getTrCurveSlaveGPs();

        // Get the corresponding Gauss weights
        trCurveGPWeights = weakIGAPatchContinuityConditions[0]->getTrCurveGPWeights();

        // Get the tangent vectors at the trimming curve of the given condition in the Cartesian space
        trCurveMasterGPTangents = weakIGAPatchContinuityConditions[0]->getTrCurveMasterGPTangents();
        trCurveSlaveGPTangents = weakIGAPatchContinuityConditions[0]->getTrCurveSlaveGPTangents();

        // Get the product of the Jacobian transformations
        trCurveGPJacobianProducts = weakIGAPatchContinuityConditions[0]->getTrCurveGPJacobianProducts();

        // Get the master and the slave patch
        patchMaster = theIGAMesh->getSurfacePatch(indexMaster);
        patchSlave = theIGAMesh->getSurfacePatch(indexSlave);

        // Get the polynomial orders of the master and the slave patch
        pMaster = patchMaster->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree();
        qMaster = patchMaster->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree();
        pSlave = patchSlave->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree();
        qSlave = patchSlave->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree();

        // get the number of local basis functions for master and slave patch
        noLocalBasisFctsMaster = (pMaster + 1)*(qMaster + 1);
        noLocalBasisFctsSlave = (pSlave + 1)*(qSlave + 1);

        // get the number of the local DOFs for the master and slave patch
        noDOFsLocMaster = noCoord*noLocalBasisFctsMaster;
        noDOFsLocSlave = noCoord*noLocalBasisFctsSlave;

        // Initialize pointers
        double* BDisplacementsGCMaster = new double[noCoord*noDOFsLocMaster];
        double* BDisplacementsGCSlave = new double[noCoord*noDOFsLocSlave];
        double* BOperatorOmegaTMaster = new double[noDOFsLocMaster];
        double* BOperatorOmegaTSlave = new double[noDOFsLocSlave];
        double* BOperatorOmegaNMaster = new double[noDOFsLocMaster];
        double* BOperatorOmegaNSlave = new double[noDOFsLocSlave];

        // Create and destroy the same objects iteratively
        for(int i = 0; i < 1e10; i++)
            // Loop over all the Gauss points
            for(int iGP = 0; iGP < noGPsOnContCond; iGP++){

                // Get the parametric coordinates of the Gauss Point on the master patch
                uGPMaster = trCurveMasterGPs[2*iGP];
                vGPMaster = trCurveMasterGPs[2*iGP + 1];

                // Get the parametric coordinates of the Gauss Point on the slave patch
                uGPSlave = trCurveSlaveGPs[2*iGP];
                vGPSlave = trCurveSlaveGPs[2*iGP + 1];

                // Find the knot span indices of the Gauss point locations in the parameter space of the master patch
                uKnotSpanMaster = patchMaster->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(uGPMaster);
                vKnotSpanMaster = patchMaster->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(vGPMaster);

                // Find the knot span indices of the Gauss point locations in the parameter space of the slave patch
                uKnotSpanSlave = patchSlave->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(uGPSlave);
                vKnotSpanSlave = patchSlave->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(vGPSlave);

                // Get the tangent to the boundary vector on the master and the slave patch
                for(int iCoord = 0; iCoord < noCoord; iCoord++){
                    tangentTrCurveVctMaster[iCoord] = trCurveMasterGPTangents[3*iGP + iCoord];
                    tangentTrCurveVctSlave[iCoord] = trCurveSlaveGPTangents[3*iGP + iCoord];
                }

                // Compute the B-operator matrices needed for the computation of the patch weak continuity contributions at the master patch
                theMapper->computeDisplacementAndRotationBOperatorMatrices(BDisplacementsGCMaster, BOperatorOmegaTMaster, BOperatorOmegaNMaster, normalTrCurveVctMaster,
                                                                           patchMaster, tangentTrCurveVctMaster, uGPMaster, vGPMaster, uKnotSpanMaster, vKnotSpanMaster);

                // Compute the B-operator matrices needed for the computation of the patch weak continuity contributions at the slave patch
                theMapper->computeDisplacementAndRotationBOperatorMatrices(BDisplacementsGCSlave, BOperatorOmegaTSlave, BOperatorOmegaNSlave, normalTrCurveVctSlave,
                                                                           patchSlave, tangentTrCurveVctSlave, uGPSlave, vGPSlave, uKnotSpanSlave, vKnotSpanSlave);
            }

        // Delete pointers
        delete[] BDisplacementsGCMaster;
        delete[] BDisplacementsGCSlave;
        delete[] BOperatorOmegaTMaster;
        delete[] BOperatorOmegaTSlave;
        delete[] BOperatorOmegaNMaster;
        delete[] BOperatorOmegaNSlave;
    }

    /***********************************************************************************************
     * \brief Test case: testComputePenaltyFactorsIGAPatchContinuityConditions4Leakage
     ***********/

    void testComputePenaltyParametersIGAPatchContinuityConditions4Leakage(){

        // Initialize empty string
        std::string emptyString = "";

        // Add the precomputed GP data to the continuity condition
        addPrecomputedGPData();

        // Set the computation of the Penalty parameters to automatic
        theMapper->setParametersWeakPatchContinuityConditions(true, true, 0.0, 0.0, 0.0);

        // Create and destroy the same objects iteratively
        for(int i = 0; i < 1e10; i++)
            theMapper->computePenaltyParametersForPatchContinuityConditions(emptyString);
    }

    /***********************************************************************************************
     * \brief Test case: testCreateWeakContinuityConditionGPData4Leakage
     ***********/

    void testCreateWeakContinuityConditionGPData4Leakage() {
        for(int i = 0; i < 1e10; i++)
            testCreateWeakContinuityConditionGPData();
    }

    CPPUNIT_TEST_SUITE (TestIGAMortarMapperWeakContinuityConditions);

    // Make the tests
    CPPUNIT_TEST (testComputeBOperatorMatricesIGAPatchContinuityConditions);
    CPPUNIT_TEST (testComputePenaltyFactorsIGAPatchContinuityConditions);
    CPPUNIT_TEST (testCreateWeakContinuityConditionGPData);

    // Make the tests for memory leakage
    //    CPPUNIT_TEST (testComputeBOperatorMatricesIGAPatchContinuityConditions4Leakage);
    //    CPPUNIT_TEST (testComputePenaltyFactorsIGAPatchContinuityConditions4Leakage);
    //    CPPUNIT_TEST (testCreateWeakContinuityConditionGPData4Leakage);       // to use this leakage test comment out the assert in the function inside the class
    CPPUNIT_TEST_SUITE_END();
}
;

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION (EMPIRE::TestIGAMortarMapperWeakContinuityConditions);
