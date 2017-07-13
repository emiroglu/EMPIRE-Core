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
class TestIGAMortarMapperCreateInterfaceGPs: public CppUnit::TestFixture {

private:
    IGAMesh* theIGAMesh;
        
    double Tol;
    double TolRel10;
    double TolRel100;
    double TolRel1000;

public:
    void setUp() {

        // Assign a tolerance value (corresponding to maximum accuracy provided by MATLAB)
        Tol = 1e-15;
        TolRel10 = Tol*1e1;
        TolRel100 = Tol*1e2;
        TolRel1000 = Tol*1e3;

        // Patch 1

        // Polynomial orders
        int p1 = 1;
        int q1 = 2;

        // Number of knots at each parametric direction
        int noUKnots1 = 4;
        int noVKnots1 = 8;

        // Knot vectors
        double knotVectorU1[noUKnots1];
        knotVectorU1[0] = 0.0000000000e+00;
        knotVectorU1[1] = 0.0000000000e+00;
        knotVectorU1[2] = 1.0000000000e+00;
        knotVectorU1[3] = 1.0000000000e+00;
        
        double knotVectorV1[noVKnots1];
        knotVectorV1[0] = 0.00000000;
        knotVectorV1[1] = 0.00000000;
	knotVectorV1[2] = 0.00000000;
        knotVectorV1[3] = 0.50000000;
        knotVectorV1[4] = 0.50000000;
        knotVectorV1[5] = 1.00000000;
        knotVectorV1[6] = 1.00000000;
        knotVectorV1[7] = 1.00000000;

        // Number of Control Points at each parametric direction
        int noUCP1 = noUKnots1 - p1 - 1;
        int noVCP1 = noVKnots1 - q1 - 1;

        // Contol Point net
        double CP1[4*noUCP1*noVCP1];
        CP1[0 * 4 + 0] = 0.0000000;
        CP1[0 * 4 + 1] = 1.0000000;
        CP1[0 * 4 + 2] = 0.0000000;
        CP1[0 * 4 + 3] = 1.0000000;

        CP1[1 * 4 + 0] = 2.0000000;
        CP1[1 * 4 + 1] = 1.0000000;
        CP1[1 * 4 + 2] = 0.0000000;
        CP1[1 * 4 + 3] = 1.0000000;

        CP1[2 * 4 + 0] = 0.0000000;
        CP1[2 * 4 + 1] = 1.0000000;
        CP1[2 * 4 + 2] = 1.0000000;
        CP1[2 * 4 + 3] = 0.707106781;

        CP1[3 * 4 + 0] = 2.0000000;
        CP1[3 * 4 + 1] = 1.0000000;
        CP1[3 * 4 + 2] = 1.0000000;
        CP1[3 * 4 + 3] = 0.707106781;

        CP1[4 * 4 + 0] = 0.0000000;
        CP1[4 * 4 + 1] = 0.0000000;
        CP1[4 * 4 + 2] = 1.0000000;
        CP1[4 * 4 + 3] = 1.0000000;

        CP1[5 * 4 + 0] = 2.0000000;
        CP1[5 * 4 + 1] = 0.0000000;
        CP1[5 * 4 + 2] = 1.0000000;
        CP1[5 * 4 + 3] = 1.0000000;

        CP1[6 * 4 + 0] = 0.0000000;
        CP1[6 * 4 + 1] = -1.0000000;
        CP1[6 * 4 + 2] = 1.0000000;
        CP1[6 * 4 + 3] = 0.707106781;

        CP1[7 * 4 + 0] = 2.0000000;
        CP1[7 * 4 + 1] = -1.0000000;
        CP1[7 * 4 + 2] = 1.0000000;
        CP1[7 * 4 + 3] = 0.707106781;

        CP1[8 * 4 + 0] = 0.0000000;
        CP1[8 * 4 + 1] = -1.0000000;
        CP1[8 * 4 + 2] = 0.0000000;
        CP1[8 * 4 + 3] = 1.0000000;

        CP1[9 * 4 + 0] = 2.0000000;
        CP1[9 * 4 + 1] = -1.0000000;
        CP1[9 * 4 + 2] = 0.0000000;
        CP1[9 * 4 + 3] = 1.0000000;

        int dofIndex1[noUCP1*noVCP1];
        for(int i = 0; i < noUCP1*noVCP1; i++)
            dofIndex1[i] = i;
	
	// Trimming Loops
	int noTrimmingLoops1 = 1;
	int noTrimmingCurves1 = 5;
	
	// Trimming Curves
	const int direction1_1 = 1;
	const int p1_1 = 2;
	const int noUKnots1_1 = 6;
	const int noCP1_1 = noUKnots1_1 - p1_1 - 1;
	double knotVectorU1_1[noUKnots1_1] = {0.0, 0.0, 0.0, 0.5, 0.5, 0.5};
	double CP1_1[4*noCP1_1] = {1.0, 0.0,  0.0, 1.0,
				   1.0, 0.25, 0.0, 1.0,
				   1.0, 0.5,  0.0, 1.0};
	
	const int direction1_2 = 1;
	const int p1_2 = 2;
	const int noUKnots1_2 = 6;
	const int noCP1_2 = noUKnots1_2 - p1_2 - 1;
	double knotVectorU1_2[noUKnots1_2] = {0.5, 0.5, 0.5, 1.0, 1.0, 1.0};
	double CP1_2[4*noCP1_2] = {1.0, 0.5,  0.0, 1.0,
				   1.0, 0.75, 0.0, 1.0,
				   1.0, 1.0,  0.0, 1.0};
				   
	const int direction1_3 = 1;
	const int p1_3 = 2;
	const int noUKnots1_3 = 6;
	const int noCP1_3 = noUKnots1_3 - p1_3 - 1;
	double knotVectorU1_3[noUKnots1_3] = {0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
	double CP1_3[4*noCP1_3] = {1.0,  1.0, 0.0, 1.0,
				   0.75, 1.0, 0.0, 1.0,
				   0.5,  1.0, 0.0, 1.0};
				   
	const int direction1_4 = 0;
	const int p1_4 = 8;
	const int noUKnots1_4 = 40;
	const int noCP1_4 = noUKnots1_4 - p1_4 - 1;
	double knotVectorU1_4[noUKnots1_4] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.4375, 0.4375, 0.4375, 0.4375, 0.4375, 0.4375, 0.4375, 0.71875, 0.71875, 0.71875, 0.71875, 0.71875, 0.71875, 0.71875, 0.71875, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
	double CP1_4[4*noCP1_4] = {0.5, 		-3.44245186937e-17, 	0.0, 	1.0,
				   0.50030037542, 	0.0398228834965, 	0.0, 	1.0,
				   0.492191690635, 	0.0767294988321, 	0.0, 	1.0,
				   0.477227781906, 	0.111198936067, 	0.0, 	1.0,
				   0.456962485494, 	0.143710285262, 	0.0, 	1.0,
				   0.43294963766, 	0.174742636477, 	0.0, 	1.0,
				   0.406743074664, 	0.204775079773, 	0.0, 	1.0,
				   0.379896632766, 	0.23428670521, 		0.0, 	1.0,
				   0.334833104, 	0.285850231955, 	0.0, 	1.0,
				   0.314519261108, 	0.309305546752, 	0.0, 	1.0,
				   0.29936384996, 	0.331403598518, 	0.0, 	1.0,
				   0.28422414737, 	0.354951202152, 	0.0, 	1.0,
				   0.270503061637, 	0.38198057788, 		0.0, 	1.0,
				   0.269491278928, 	0.410972614516, 	0.0, 	1.0,
				   0.291021970955, 	0.437603754288, 	0.0, 	1.0,
				   0.357340844024, 	0.491403987965, 	0.0, 	1.0,
				   0.41157980266, 	0.519083903877, 	0.0, 	1.0,
				   0.462022701318, 	0.536681284459, 	0.0, 	1.0,
				   0.523119120189, 	0.570801866401, 	0.0, 	1.0,
				   0.551623424242, 	0.577303508464, 	0.0, 	1.0,
				   0.640275727184, 	0.616174436182, 	0.0, 	1.0,
				   0.632026292042, 	0.658214712613, 	0.0, 	1.0,
				   0.624010438262, 	0.695492901102, 	0.0, 	1.0,
				   0.614063834049, 	0.733227683504, 	0.0, 	1.0,
				   0.591880703549, 	0.7686881693, 		0.0, 	1.0,
				   0.566389565004, 	0.803274738077, 	0.0, 	1.0,
				   0.542569192962, 	0.838321121331, 	0.0, 	1.0,
				   0.522782014141, 	0.874961424859, 	0.0, 	1.0,
				   0.508107503301, 	0.91399715116, 		0.0, 	1.0,
				   0.499675579108, 	0.955764221839, 	0.0, 	1.0,
				   0.5,			1.0, 			0.0, 	1.0};
				   
	const int direction1_5 = 1;
	const int p1_5 = 2;
	const int noUKnots1_5 = 6;
	const int noCP1_5 = noUKnots1_5 - p1_5 - 1;
	double knotVectorU1_5[noUKnots1_5] = {0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
	double CP1_5[4*noCP1_5] = {0.5, -4.28187081941e-17, 0.0, 1.0,
				   0.75, 3.38443719739e-17, 0.0, 1.0,
				   1.0, 0.0, 0.0, 1.0};

        // Patch 2

        // Polynomial orders
        const int p2 = 3;
        const int q2 = 1;

        // Number of knots at each parametric direction
        const int noUKnots2 = 10;
        const int noVKnots2 = 4;

        // Knot vectors
        double knotVectorU2[noUKnots2];
        knotVectorU2[0] = 0.0000000000e+00;
        knotVectorU2[1] = 0.0000000000e+00;
        knotVectorU2[2] = 0.0000000000e+00;
	knotVectorU2[3] = 0.0000000000e+00;
        knotVectorU2[4] = 0.372578242;
        knotVectorU2[5] = 0.746413635;
        knotVectorU2[6] = 1.0000000000e+00;
	knotVectorU2[7] = 1.0000000000e+00;
	knotVectorU2[8] = 1.0000000000e+00;
	knotVectorU2[9] = 1.0000000000e+00;

        double knotVectorV2[noVKnots2];
        knotVectorV2[0] = 0.00000000;
        knotVectorV2[1] = 0.00000000;
        knotVectorV2[2] = 1.00000000;
        knotVectorV2[3] = 1.00000000;

        // Number of Control Points at each parametric direction
        int noUCP2 = noUKnots2 - p2 - 1;
        int noVCP2 = noVKnots2 - q2 - 1;

        // Contol Point net
        double CP2[4*noUCP2*noVCP2];
        CP2[0 * 4 + 0] = 1.0000000;
        CP2[0 * 4 + 1] = 1.0000000;
        CP2[0 * 4 + 2] = 0.0000000;
        CP2[0 * 4 + 3] = 1.0000000;

        CP2[1 * 4 + 0] = 0.736652716;
        CP2[1 * 4 + 1] = 0.751117883;
        CP2[1 * 4 + 2] = 0.0000000;
        CP2[1 * 4 + 3] = 1.0000000;

        CP2[2 * 4 + 0] = 0.219597039;
        CP2[2 * 4 + 1] = 0.105638419;
        CP2[2 * 4 + 2] = 0.0000000;
        CP2[2 * 4 + 3] = 1.0000000;

        CP2[3 * 4 + 0] = 1.4923261;
        CP2[3 * 4 + 1] = -0.288356341;
        CP2[3 * 4 + 2] = 0.0000000;
        CP2[3 * 4 + 3] = 1.0000000;

        CP2[4 * 4 + 0] = 1.15830572;
        CP2[4 * 4 + 1] = -0.810892925;
        CP2[4 * 4 + 2] = 0.0000000;
        CP2[4 * 4 + 3] = 1.0000000;

        CP2[5 * 4 + 0] = 1.0000000;
        CP2[5 * 4 + 1] = -1.0000000;
        CP2[5 * 4 + 2] = 0.0000000;
        CP2[5 * 4 + 3] = 1.0000000;

        CP2[6 * 4 + 0] = 1.0000000;
        CP2[6 * 4 + 1] = 1.0000000;
        CP2[6 * 4 + 2] = 2.0000000;
        CP2[6 * 4 + 3] = 1.0000000;

        CP2[7 * 4 + 0] = 0.736652716;
        CP2[7 * 4 + 1] = 0.751117883;
        CP2[7 * 4 + 2] = 2.0000000;
        CP2[7 * 4 + 3] = 1.0000000;

        CP2[8 * 4 + 0] = 0.219597039;
        CP2[8 * 4 + 1] = 0.105638419;
        CP2[8 * 4 + 2] = 2.0000000;
        CP2[8 * 4 + 3] = 1.0000000;

        CP2[9 * 4 + 0] = 1.4923261;
        CP2[9 * 4 + 1] = -0.288356341;
        CP2[9 * 4 + 2] = 2.0000000;
        CP2[9 * 4 + 3] = 1.0000000;

        CP2[10 * 4 + 0] = 1.15830572;
        CP2[10 * 4 + 1] = -0.810892925;
        CP2[10 * 4 + 2] = 2.0000000;
        CP2[10 * 4 + 3] = 1.0000000;

        CP2[11 * 4 + 0] = 1.0000000;
        CP2[11 * 4 + 1] = -1.0000000;
        CP2[11 * 4 + 2] = 2.0000000;
        CP2[11 * 4 + 3] = 1.0000000;

        int dofIndex2[noUCP2*noVCP2];
        for(int i = 0; i < noUCP2*noVCP2; i++)
            dofIndex2[i] = i;
	
	// Trimming Loops
	const int noTrimmingLoops2 = 1;
	const int noTrimmingCurves2 = 4;
	
	// Trimming Curves
	const int direction2_1 = 1;
	const int p2_1 = 2;
	const int noUKnots2_1 = 6;
	const int noCP2_1 = noUKnots2_1 - p2_1 - 1;
	double knotVectorU2_1[noUKnots2_1] = {0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
	double CP2_1[4*noCP2_1] = {1.0, 1.0, 0.0, 1.0,
				   0.5, 1.0, 0.0, 1.0,
				   2.68391769566e-17, 1.0, 0.0, 1.0};
	
	const int direction2_2 = 1;
	const int p2_2 = 2;
	const int noUKnots2_2 = 6;
	const int noCP2_2 = noUKnots2_2 - p2_2 - 1;
	double knotVectorU2_2[noUKnots2_2] = {0.0, 0.0, 0.0, 2.0, 2.0, 2.0};
	double CP2_2[4*noCP2_2] = {2.68391769566e-17, 1.0, 0.0, 1.0,
				   1.23207278139e-17, 0.5, 0.0, 1.0,
				   0.0, 2.66686007337e-18, 0.0, 1.0};

	const int direction2_3 = 1;
	const int p2_3 = 2;
	const int noUKnots2_3 = 6;
	const int noCP2_3 = noUKnots2_3 - p2_3 - 1;
	double knotVectorU2_3[noUKnots2_3] = {0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
	double CP2_3[4*noCP2_3] = {1.0,  1.0, 0.0, 1.0,
				   0.75, 1.0, 0.0, 1.0,
				   0.5,  1.0, 0.0, 1.0};
				   
	const int direction2_4 = 1;
	const int p2_4 = 8;
	const int noUKnots2_4 = 41;
	const int noCP2_4 = noUKnots2_4 - p2_4 - 1;
	double knotVectorU2_4[noUKnots2_4] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.4375, 0.4375, 0.4375, 0.4375, 0.4375, 0.4375, 0.4375, 0.4375, 0.71875, 0.71875, 0.71875, 0.71875, 0.71875, 0.71875, 0.71875, 0.71875, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
	double CP2_4[4*noCP2_4] = {-1.30104260698e-18,	0.0,	0.0,	1.0,
				  2.45106608077e-06,	0.0565081811119,	0.0,	1.0,
				  0.00705074428197,	0.108094322017,		0.0,	1.0,
				  0.0212852905443,	0.170831896976,		0.0,	1.0,
				  0.0407969358585,	0.207258034447,		0.0,	1.0,
				  0.0635009003927,	0.262397610619,		0.0,	1.0,
				  0.0893660563479,	0.297751092839,		0.0,	1.0,
				  0.117966117914,	0.336278660229,		0.0,	1.0,
				  0.17402819912,	0.394421101818,		0.0,	1.0,
				  0.200434767026,	0.41666940025,		0.0,	1.0,
				  0.228085712907,	0.437503136511,		0.0,	1.0,
				  0.258597695163,	0.456694607912,		0.0,	1.0,
				  0.294148232641,	0.473096228015,		0.0,	1.0,
				  0.335245119684,	0.485411857299,		0.0,	1.0,
				  0.378495841181,	0.492968133829,		0.0,	1.0,
				  0.414376987617,	0.496485803925,		0.0,	1.0,
				  0.469957729003,	0.502183260934,		0.0,	1.0,
				  0.513025838484,	0.5011073499,		0.0,	1.0,
				  0.564046962592,	0.499757983673,		0.0,	1.0,
				  0.601724586028,	0.492743734166,		0.0,	1.0,
				  0.642201591882,	0.487433067147,		0.0,	1.0,
				  0.704476379679,	0.479846158533,		0.0,	1.0,
				  0.765871154184,	0.445781979382,		0.0,	1.0,
				  0.811389997818,	0.411175336737,		0.0,	1.0,
				  0.857543515679,	0.375411583451,		0.0,	1.0,
				  0.895395538906,	0.335474224691,		0.0,	1.0,
				  0.92701825489,	0.290482422218,		0.0,	1.0,
				  0.953566439811,	0.240142803282,		0.0,	1.0,
				  0.975277458635,	0.184749460621,		0.0,	1.0,
				  0.991471265117,	0.125183952463,		0.0,	1.0,
				  1.0005504018,		0.0629153025212,	0.0,	1.0,
				  1.0,	0.0,	0.0,	1.0};

        // Total number of Gauss Points
        int noCPs = noUCP1*noVCP1 + noUCP2*noVCP2;

        // Construct the NURBS multipatch geometry
        theIGAMesh = new IGAMesh("myMesh", noCPs);

	// Add the underlying patches of the multipatch geometry
        theIGAMesh->addPatch(p1, noUKnots1, knotVectorU1, 
					 q1, noVKnots1, knotVectorV1,
					 noUCP1, noVCP1, CP1, dofIndex1);
	
	
	
	theIGAMesh->getSurfacePatch(0)->addTrimLoop(0, noTrimmingCurves1);
	theIGAMesh->getSurfacePatch(0)->addTrimCurve(direction1_1, p1_1, noUKnots1_1, knotVectorU1_1, noCP1_1, CP1_1);
	theIGAMesh->getSurfacePatch(0)->addTrimCurve(direction1_2, p1_2, noUKnots1_2, knotVectorU1_2, noCP1_2, CP1_2);
	theIGAMesh->getSurfacePatch(0)->addTrimCurve(direction1_3, p1_3, noUKnots1_3, knotVectorU1_3, noCP1_3, CP1_3);
	theIGAMesh->getSurfacePatch(0)->addTrimCurve(direction1_4, p1_4, noUKnots1_4, knotVectorU1_4, noCP1_4, CP1_4);
	theIGAMesh->getSurfacePatch(0)->addTrimCurve(direction1_5, p1_5, noUKnots1_5, knotVectorU1_5, noCP1_5, CP1_5);
	theIGAMesh->getSurfacePatch(0)->linearizeTrimming();

        theIGAMesh->addPatch(p2, noUKnots2, knotVectorU2, 
					 q2, noVKnots2, knotVectorV2,
					 noUCP2, noVCP2, CP2, dofIndex2);
	
	theIGAMesh->getSurfacePatch(1)->addTrimLoop(0, noTrimmingCurves2);
	theIGAMesh->getSurfacePatch(1)->addTrimCurve(direction2_1, p2_1, noUKnots2_1, knotVectorU2_1, noCP2_1, CP2_1);
	theIGAMesh->getSurfacePatch(1)->addTrimCurve(direction2_2, p2_2, noUKnots2_2, knotVectorU2_2, noCP2_2, CP2_2);
	theIGAMesh->getSurfacePatch(1)->addTrimCurve(direction2_3, p2_3, noUKnots2_3, knotVectorU2_3, noCP2_3, CP2_3);
	theIGAMesh->getSurfacePatch(1)->addTrimCurve(direction2_4, p2_4, noUKnots2_4, knotVectorU2_4, noCP2_4, CP2_4);
	theIGAMesh->getSurfacePatch(1)->linearizeTrimming();
	
	// Initialize auxiliary variables
        int connectionCtr = 0;
        int masterPatchIndex = 0;
        int masterPatchBLIndex = 0;
        int masterPatchBLTrCurveIndex = 3;
        int slavePatchIndex = 1;
        int slavePatchBLIndex = 0;
        int slavePatchBLTrCurveIndex = 3;
        int isGPProvided = 0;

        // Add the coupling data corresponding to the coupling between the patches
        theIGAMesh->addWeakContinuityCondition(connectionCtr,
					       masterPatchIndex, masterPatchBLIndex, masterPatchBLTrCurveIndex,
					       slavePatchIndex, slavePatchBLIndex, slavePatchBLTrCurveIndex,
					       isGPProvided);
	
	
	
    }

    void tearDown() {
      delete theIGAMesh;
    }
    /***********************************************************************************************
     * \brief Test case: Test the constructor
     ***********/

    void testCreateWeakContinuityConditionGPData() {
        theIGAMesh->createWeakContinuityConditionGPData(0);
    }

    void test4Leakage() {
        for (int i = 0; i < 100000000; i++) {
            testCreateWeakContinuityConditionGPData();
        }
    }

// Make the tests
    CPPUNIT_TEST_SUITE (TestIGAMortarMapperCreateInterfaceGPs);
    CPPUNIT_TEST (testCreateWeakContinuityConditionGPData);
//    CPPUNIT_TEST (test4Leakage);
    CPPUNIT_TEST_SUITE_END();
}
;

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION (EMPIRE::TestIGAMortarMapperCreateInterfaceGPs);
