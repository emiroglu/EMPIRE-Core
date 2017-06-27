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
class TestIGAMortarMapperWeakContinuityConditions: public CppUnit::TestFixture {

private:
    IGAMortarMapper* theMapper;
    IGAMesh* theIGAMesh;
    FEMesh* theFEMesh;
    bool isMappingIGA2FEM;
    double Tol;
    double TolRel100;

public:
    void setUp() {
        isMappingIGA2FEM = false;

        // Assign a tolerance value (corresponding to maximum accuracy provided by MATLAB)
        Tol = 1e-15;
        TolRel100 = Tol*1e2;

        // Patch 1

        // Polynomial orders
        int p1 = 2;
        int q1 = 2;

        // Number of knots at each parametric direction
        int noUKnots1 = 6;
        int noVKnots1 = 6;

        // Knot vectors
        double knotVectotU1[noUKnots1];
        knotVectotU1[0] = 0.0000000000e+00;
        knotVectotU1[1] = 0.0000000000e+00;
        knotVectotU1[2] = 0.0000000000e+00;
        knotVectotU1[3] = 5.0000000000e+00;
        knotVectotU1[4] = 5.0000000000e+00;
        knotVectotU1[5] = 5.0000000000e+00;

        double knotVectotV1[noVKnots1];
        knotVectotV1[0] = 0.00000000;
        knotVectotV1[1] = 0.00000000;
        knotVectotV1[2] = 0.00000000;
        knotVectotV1[3] = 1.00000000;
        knotVectotV1[4] = 1.00000000;
        knotVectotV1[5] = 1.00000000;

        // Number of Control Points at each parametric direction
        int noUCP1 = noUKnots1 - p1 - 1;
        int noVCP1 = noVKnots1 - q1 - 1;

        // Contol Point net
        double CP1[4*noUCP1*noVCP1];
        CP1[0 * 4 + 0] = -0.003908007399212;
        CP1[0 * 4 + 1] = -0.005234065460828;
        CP1[0 * 4 + 2] = -1.332504672648300e-10;
        CP1[0 * 4 + 3] = 1.0000000;

        CP1[1 * 4 + 0] = 0.114028181791499;
        CP1[1 * 4 + 1] = 0.203869643922430;
        CP1[1 * 4 + 2] = 4.657780423778900e-10;
        CP1[1 * 4 + 3] = 1.0000000;

        CP1[2 * 4 + 0] = 0.316023042789807;
        CP1[2 * 4 + 1] = 0.256681222103740;
        CP1[2 * 4 + 2] = 2.586590240316700e-10;
        CP1[2 * 4 + 3] = 1.0000000;

        CP1[3 * 4 + 0] = -0.003612528584176;
        CP1[3 * 4 + 1] = -0.002702058212826;
        CP1[3 * 4 + 2] = 0.049999999866791;
        CP1[3 * 4 + 3] = 0.707106781186547;

        CP1[4 * 4 + 0] = 0.138291045935836;
        CP1[4 * 4 + 1] = 0.197683752569160;
        CP1[4 * 4 + 2] = 0.050000000461007;
        CP1[4 * 4 + 3] = 0.707106781186547;

        CP1[5 * 4 + 0] = 0.366668497560549;
        CP1[5 * 4 + 1] = 0.259368994510050;
        CP1[5 * 4 + 2] = 0.050000000266021;
        CP1[5 * 4 + 3] = 0.707106781186547;

        CP1[6 * 4 + 0] = -0.003771817327604;
        CP1[6 * 4 + 1] = -0.005554345778163;
        CP1[6 * 4 + 2] = 0.099999999860943;
        CP1[6 * 4 + 3] = 1.0000000;

        CP1[7 * 4 + 0] = 0.113620343377848;
        CP1[7 * 4 + 1] = 0.204199374706320;
        CP1[7 * 4 + 2] = 0.100000000473276;
        CP1[7 * 4 + 3] = 1.0000000;

        CP1[8 * 4 + 0] = 0.316359281162531;
        CP1[8 * 4 + 1] = 0.257017472117850;
        CP1[8 * 4 + 2] = 0.100000000256655;
        CP1[8 * 4 + 3] = 1.0000000;

        int dofIndex1[noUCP1*noVCP1];
        for(int i = 0; i < noUCP1*noVCP1; i++)
            dofIndex1[i] = i;

        // Patch 2

        // Polynomial orders
        int p2 = 2;
        int q2 = 2;

        // Number of knots at each parametric direction
        int noUKnots2 = 6;
        int noVKnots2 = 6;

        // Knot vectors
        double knotVectotU2[noUKnots2];
        knotVectotU2[0] = 0.0000000000e+00;
        knotVectotU2[1] = 0.0000000000e+00;
        knotVectotU2[2] = 0.0000000000e+00;
        knotVectotU2[3] = 5.0000000000e+00;
        knotVectotU2[4] = 5.0000000000e+00;
        knotVectotU2[5] = 5.0000000000e+00;

        double knotVectotV2[noVKnots2];
        knotVectotV2[0] = 0.00000000;
        knotVectotV2[1] = 0.00000000;
        knotVectotV2[2] = 0.00000000;
        knotVectotV2[3] = 1.00000000;
        knotVectotV2[4] = 1.00000000;
        knotVectotV2[5] = 1.00000000;

        // Number of Control Points at each parametric direction
        int noUCP2 = noUKnots2 - p2 - 1;
        int noVCP2 = noVKnots2 - q2 - 1;

        // Contol Point net
        double CP2[4*noUCP2*noVCP2];
        CP2[0 * 4 + 0] = 0.316003856904912;
        CP2[0 * 4 + 1] = 0.256688793621200;
        CP2[0 * 4 + 2] = 2.582876756262200e-10;
        CP2[0 * 4 + 3] = 1.0000000;

        CP2[1 * 4 + 0] = 0.736485900451057;
        CP2[1 * 4 + 1] = 0.352668582279320;
        CP2[1 * 4 + 2] = 5.144890391225400e-10;
        CP2[1 * 4 + 3] = 1.0000000;

        CP2[2 * 4 + 0] = 1.009344658935689;
        CP2[2 * 4 + 1] = -0.007427201607165;
        CP2[2 * 4 + 2] = -6.052522674573600e-11;
        CP2[2 * 4 + 3] = 1.0000000;

        CP2[3 * 4 + 0] = 0.366663956039134;
        CP2[3 * 4 + 1] = 0.259375572812560;
        CP2[3 * 4 + 2] = 0.050000000265932;
        CP2[3 * 4 + 3] = 0.707106781186547;

        CP2[4 * 4 + 0] = 0.761015516537076;
        CP2[4 * 4 + 1] = 0.350232764137380;
        CP2[4 * 4 + 2] = 0.050000000507262;
        CP2[4 * 4 + 3] = 0.707106781186547;

        CP2[5 * 4 + 0] = 1.009405620078563;
        CP2[5 * 4 + 1] = -0.007974181436323;
        CP2[5 * 4 + 2] = 0.049999999940135;
        CP2[5 * 4 + 3] = 0.707106781186547;

        CP2[6 * 4 + 0] = 0.316341098950630;
        CP2[6 * 4 + 1] = 0.257025224191480;
        CP2[6 * 4 + 2] = 0.100000000256298;
        CP2[6 * 4 + 3] = 1.0000000;

        CP2[7 * 4 + 0] = 0.736033287858653;
        CP2[7 * 4 + 1] = 0.353010531551170;
        CP2[7 * 4 + 2] = 0.100000000520956;
        CP2[7 * 4 + 3] = 1.0000000;

        CP2[8 * 4 + 0] = 1.009826851544427;
        CP2[8 * 4 + 1] = -0.008501105140788;
        CP2[8 * 4 + 2] = 0.099999999934003;
        CP2[8 * 4 + 3] = 1.0000000;

        int dofIndex2[noUCP2*noVCP2];
        for(int i = 0; i < noUCP2*noVCP2; i++)
            dofIndex2[i] = i;

        // Total number of Gauss Points
        int noCPs = noUCP1*noVCP1 + noUCP2*noVCP2;

        // Construct the NURBS multipatch geometry
        theIGAMesh = new IGAMesh("modifiedCavity2Patches", noCPs);

        // Add the underlying patches of the multipatch geometry
        theIGAMesh->addPatch(p1, noUKnots1, knotVectotU1, q1, noVKnots1, knotVectotV1,
                             noUCP1, noVCP1, CP1, dofIndex1);

        theIGAMesh->addPatch(p2, noUKnots2, knotVectotU2, q2, noVKnots2, knotVectotV2,
                             noUCP2, noVCP2, CP2, dofIndex2);

        // Initialize auxiliary variables
        int noGPsConnection = 3;
        int connectionCtr = 0;
        int masterPatchIndex = 0;
        int masterPatchBLIndex = 0;
        int masterPatchBLTrCurveIndex = 1;
        int slavePatchIndex = 1;
        int slavePatchBLIndex = 0;
        int slavePatchBLTrCurveIndex = 3;
        int isGPProvided = 1;

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
        trCurveMasterGPs[2*0 + 1] = 0.112701665379258;
        trCurveMasterGPs[2*1 + 0] = 5;
        trCurveMasterGPs[2*1 + 1] = 0.500000000000000;
        trCurveMasterGPs[2*2 + 0] = 5;
        trCurveMasterGPs[2*2 + 1] = 0.887298334620742;

        // The parametric coordinates of the Gauss Points on the slave
        trCurveSlaveGPs[2*0 + 0] = 0;
        trCurveSlaveGPs[2*0 + 1] = 0.112701665379258;
        trCurveSlaveGPs[2*1 + 0] = 0;
        trCurveSlaveGPs[2*1 + 1] = 0.500000000000000;
        trCurveSlaveGPs[2*2 + 0] = 0;
        trCurveSlaveGPs[2*2 + 1] = 0.887298334620742;

        // The Gauss Point weights for each Gauss Point
        trCurveGPWeights[0] = 0.555555555555556;
        trCurveGPWeights[1] = 0.888888888888889;
        trCurveGPWeights[2] = 0.555555555555556;

        // The parametric components of the tangent to the trimming curve vectors on the master curve
        trCurveMasterGPTangents[3*0 + 0] = 0.586948988282539;
        trCurveMasterGPTangents[3*0 + 1] = 0.031882855063476;
        trCurveMasterGPTangents[3*0 + 2] = 0.808995901539127;
        trCurveMasterGPTangents[3*1 + 0] = 0.003362345712789;
        trCurveMasterGPTangents[3*1 + 1] = 0.003362462125332;
        trCurveMasterGPTangents[3*1 + 2] = 0.999988694175971;
        trCurveMasterGPTangents[3*2 + 0] = -0.583458967239414;
        trCurveMasterGPTangents[3*2 + 1] = -0.026531040888276;
        trCurveMasterGPTangents[3*2 + 2] = 0.811709145825831;

        // The parametric components of the tangent to the trimming curve vectors on the slave curve
        trCurveSlaveGPTangents[3*0 + 0] = -0.587061796109393;
        trCurveSlaveGPTangents[3*0 + 1] = -0.031868572585806;
        trCurveSlaveGPTangents[3*0 + 2] = -0.808914607131159;
        trCurveSlaveGPTangents[3*1 + 0] = -0.003372382194799;
        trCurveSlaveGPTangents[3*1 + 1] = -0.003364267532487;
        trCurveSlaveGPTangents[3*1 + 2] = -0.999988654306789;
        trCurveSlaveGPTangents[3*2 + 0] = 0.583562122978769;
        trCurveSlaveGPTangents[3*2 + 1] = 0.026514587745393;
        trCurveSlaveGPTangents[3*2 + 2] = -0.811635524888607;

        // The products of the Jacobian transformations with the Gauss weights at each Gauss Point
        trCurveGPJacobianProducts[0] = 0.029664256375105;
        trCurveGPJacobianProducts[1] = 0.052070494265638;
        trCurveGPJacobianProducts[2] = 0.029565099700221;

        // Add the coupling data corresponding to the coupling between the patches
        WeakIGAPatchContinuityCondition* theWeakContCond = theIGAMesh->addWeakContinuityCondition(connectionCtr,
                                                                                                  masterPatchIndex,
                                                                                                  masterPatchBLIndex,
                                                                                                  masterPatchBLTrCurveIndex,
                                                                                                  slavePatchIndex,
                                                                                                  slavePatchBLIndex,
                                                                                                  slavePatchBLTrCurveIndex,
                                                                                                  isGPProvided);

        // Add the Gauss Point data of the weak continuity condition
        theWeakContCond->addWeakContinuityConditionGPData(noGPsConnection,
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
        theFEMesh->elems[0 * 4 + 2] = 870;
        theFEMesh->elems[0 * 4 + 3] = 869;
        theFEMesh->elems[1 * 4 + 0] = 2;
        theFEMesh->elems[1 * 4 + 1] = 3;
        theFEMesh->elems[1 * 4 + 2] = 871;
        theFEMesh->elems[1 * 4 + 3] = 870;
        theFEMesh->elems[2 * 4 + 0] = 3;
        theFEMesh->elems[2 * 4 + 1] = 4;
        theFEMesh->elems[2 * 4 + 2] = 872;
        theFEMesh->elems[2 * 4 + 3] = 871;
        theFEMesh->elems[3 * 4 + 0] = 4;
        theFEMesh->elems[3 * 4 + 1] = 5;
        theFEMesh->elems[3 * 4 + 2] = 873;
        theFEMesh->elems[3 * 4 + 3] = 872;
        theFEMesh->elems[4 * 4 + 0] = 5;
        theFEMesh->elems[4 * 4 + 1] = 6;
        theFEMesh->elems[4 * 4 + 2] = 874;
        theFEMesh->elems[4 * 4 + 3] = 873;
        theFEMesh->elems[5 * 4 + 0] = 6;
        theFEMesh->elems[5 * 4 + 1] = 7;
        theFEMesh->elems[5 * 4 + 2] = 875;
        theFEMesh->elems[5 * 4 + 3] = 874;
        theFEMesh->elems[6 * 4 + 0] = 7;
        theFEMesh->elems[6 * 4 + 1] = 8;
        theFEMesh->elems[6 * 4 + 2] = 876;
        theFEMesh->elems[6 * 4 + 3] = 875;
        theFEMesh->elems[7 * 4 + 0] = 8;
        theFEMesh->elems[7 * 4 + 1] = 9;
        theFEMesh->elems[7 * 4 + 2] = 877;
        theFEMesh->elems[7 * 4 + 3] = 876;
        theFEMesh->elems[8 * 4 + 0] = 9;
        theFEMesh->elems[8 * 4 + 1] = 10;
        theFEMesh->elems[8 * 4 + 2] = 878;
        theFEMesh->elems[8 * 4 + 3] = 877;
        theFEMesh->elems[9 * 4 + 0] = 10;
        theFEMesh->elems[9 * 4 + 1] = 11;
        theFEMesh->elems[9 * 4 + 2] = 879;
        theFEMesh->elems[9 * 4 + 3] = 878;
        theFEMesh->elems[10 * 4 + 0] = 11;
        theFEMesh->elems[10 * 4 + 1] = 12;
        theFEMesh->elems[10 * 4 + 2] = 880;
        theFEMesh->elems[10 * 4 + 3] = 879;
        theFEMesh->elems[11 * 4 + 0] = 12;
        theFEMesh->elems[11 * 4 + 1] = 13;
        theFEMesh->elems[11 * 4 + 2] = 881;
        theFEMesh->elems[11 * 4 + 3] = 880;
        theFEMesh->elems[12 * 4 + 0] = 13;
        theFEMesh->elems[12 * 4 + 1] = 14;
        theFEMesh->elems[12 * 4 + 2] = 882;
        theFEMesh->elems[12 * 4 + 3] = 881;
        theFEMesh->elems[13 * 4 + 0] = 14;
        theFEMesh->elems[13 * 4 + 1] = 15;
        theFEMesh->elems[13 * 4 + 2] = 883;
        theFEMesh->elems[13 * 4 + 3] = 882;
        theFEMesh->elems[14 * 4 + 0] = 15;
        theFEMesh->elems[14 * 4 + 1] = 16;
        theFEMesh->elems[14 * 4 + 2] = 884;
        theFEMesh->elems[14 * 4 + 3] = 883;
        theFEMesh->elems[15 * 4 + 0] = 16;
        theFEMesh->elems[15 * 4 + 1] = 17;
        theFEMesh->elems[15 * 4 + 2] = 885;
        theFEMesh->elems[15 * 4 + 3] = 884;
        theFEMesh->elems[16 * 4 + 0] = 17;
        theFEMesh->elems[16 * 4 + 1] = 18;
        theFEMesh->elems[16 * 4 + 2] = 886;
        theFEMesh->elems[16 * 4 + 3] = 885;
        theFEMesh->elems[17 * 4 + 0] = 18;
        theFEMesh->elems[17 * 4 + 1] = 19;
        theFEMesh->elems[17 * 4 + 2] = 887;
        theFEMesh->elems[17 * 4 + 3] = 886;
        theFEMesh->elems[18 * 4 + 0] = 19;
        theFEMesh->elems[18 * 4 + 1] = 20;
        theFEMesh->elems[18 * 4 + 2] = 888;
        theFEMesh->elems[18 * 4 + 3] = 887;
        theFEMesh->elems[19 * 4 + 0] = 20;
        theFEMesh->elems[19 * 4 + 1] = 21;
        theFEMesh->elems[19 * 4 + 2] = 889;
        theFEMesh->elems[19 * 4 + 3] = 888;
        theFEMesh->elems[20 * 4 + 0] = 21;
        theFEMesh->elems[20 * 4 + 1] = 22;
        theFEMesh->elems[20 * 4 + 2] = 890;
        theFEMesh->elems[20 * 4 + 3] = 889;
        theFEMesh->elems[21 * 4 + 0] = 22;
        theFEMesh->elems[21 * 4 + 1] = 23;
        theFEMesh->elems[21 * 4 + 2] = 891;
        theFEMesh->elems[21 * 4 + 3] = 890;
        theFEMesh->elems[22 * 4 + 0] = 23;
        theFEMesh->elems[22 * 4 + 1] = 24;
        theFEMesh->elems[22 * 4 + 2] = 892;
        theFEMesh->elems[22 * 4 + 3] = 891;
        theFEMesh->elems[23 * 4 + 0] = 24;
        theFEMesh->elems[23 * 4 + 1] = 25;
        theFEMesh->elems[23 * 4 + 2] = 893;
        theFEMesh->elems[23 * 4 + 3] = 892;
        theFEMesh->elems[24 * 4 + 0] = 25;
        theFEMesh->elems[24 * 4 + 1] = 26;
        theFEMesh->elems[24 * 4 + 2] = 894;
        theFEMesh->elems[24 * 4 + 3] = 893;
        theFEMesh->elems[25 * 4 + 0] = 26;
        theFEMesh->elems[25 * 4 + 1] = 27;
        theFEMesh->elems[25 * 4 + 2] = 895;
        theFEMesh->elems[25 * 4 + 3] = 894;
        theFEMesh->elems[26 * 4 + 0] = 27;
        theFEMesh->elems[26 * 4 + 1] = 28;
        theFEMesh->elems[26 * 4 + 2] = 896;
        theFEMesh->elems[26 * 4 + 3] = 895;
        theFEMesh->elems[27 * 4 + 0] = 28;
        theFEMesh->elems[27 * 4 + 1] = 29;
        theFEMesh->elems[27 * 4 + 2] = 897;
        theFEMesh->elems[27 * 4 + 3] = 896;
        theFEMesh->elems[28 * 4 + 0] = 29;
        theFEMesh->elems[28 * 4 + 1] = 30;
        theFEMesh->elems[28 * 4 + 2] = 898;
        theFEMesh->elems[28 * 4 + 3] = 897;
        theFEMesh->elems[29 * 4 + 0] = 30;
        theFEMesh->elems[29 * 4 + 1] = 31;
        theFEMesh->elems[29 * 4 + 2] = 899;
        theFEMesh->elems[29 * 4 + 3] = 898;

        // Initialize the isogeometric mortar-based mapper
        theMapper = new IGAMortarMapper("Test IGA Mortar Mapper", theIGAMesh, theFEMesh, isMappingIGA2FEM);

        // Assign the mapper's properties
        theMapper->setParametersProjection(0.05, 2, 1e-3);
        theMapper->setParametersNewtonRaphson(10, 1e-6);
        theMapper->setParametersNewtonRaphsonBoundary(0, 1e-6);
        theMapper->setParametersBisection(100, 1e-6);
    }

    void tearDown() {

        delete theFEMesh;
        delete theIGAMesh;
//        delete theMapper;

    }
    /***********************************************************************************************
     * \brief Test case: Test the constructor
     ***********/

    void testIGAPatchContinuityConditions() {

        // Get the weak patch continuity conditions
        std::vector<WeakIGAPatchContinuityCondition*> weakIGAPatchContinuityConditions = theIGAMesh->getWeakIGAPatchContinuityConditions();

        // Define tolerances
        const double tolAngle = 1e-1;
        const double tolVct = 1e-4;

        // Initialize constant array sizes
        const int noCoord = 3;

        // Initialize varying array sizes
        int indexMaster;
        int indexSlave;
        int counter;
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
        int iGP;
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
        double* BDisplacementsGCMaster1 = new double[noCoord*noDOFsLocMaster];
        double* BDisplacementsGCSlave1 = new double[noCoord*noDOFsLocSlave];
        double* BOperatorOmegaTMaster1 = new double[noDOFsLocMaster];
        double* BOperatorOmegaTSlave1 = new double[noDOFsLocSlave];
        double* BOperatorOmegaNMaster1 = new double[noDOFsLocMaster];
        double* BOperatorOmegaNSlave1 = new double[noDOFsLocSlave];
        double* BDisplacementsGCMaster2 = new double[noCoord*noDOFsLocMaster];
        double* BDisplacementsGCSlave2 = new double[noCoord*noDOFsLocSlave];
        double* BOperatorOmegaTMaster2 = new double[noDOFsLocMaster];
        double* BOperatorOmegaTSlave2 = new double[noDOFsLocSlave];
        double* BOperatorOmegaNMaster2 = new double[noDOFsLocMaster];
        double* BOperatorOmegaNSlave2 = new double[noDOFsLocSlave];
        double* BDisplacementsGCMaster3 = new double[noCoord*noDOFsLocMaster];
        double* BDisplacementsGCSlave3 = new double[noCoord*noDOFsLocSlave];
        double* BOperatorOmegaTMaster3 = new double[noDOFsLocMaster];
        double* BOperatorOmegaTSlave3 = new double[noDOFsLocSlave];
        double* BOperatorOmegaNMaster3 = new double[noDOFsLocMaster];
        double* BOperatorOmegaNSlave3 = new double[noDOFsLocSlave];

        // Compute the B-operator matrices at the first Gauss Point
        iGP = 0;

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
        theMapper->computeIGAPatchContinuityConditionBOperatorMatrices(BDisplacementsGCMaster1, BOperatorOmegaTMaster1, BOperatorOmegaNMaster1, normalTrCurveVctMaster,
                                                                       patchMaster, tangentTrCurveVctMaster, uGPMaster, vGPMaster, uKnotSpanMaster, vKnotSpanMaster);

        // Expected values for the displacement B-operator matrix of the master patch by MATLAB
        double CorrectBOperatorDisplacementMatrixMaster1[] = { 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                               0000000000000000, 0000000000000000, 8.3628688621153702e-01, 0000000000000000,
                                                               0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                               0000000000000000, 0000000000000000, 0000000000000000, 1.5022110482233481e-01,
                                                               0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                               0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                               1.3492008966128155e-02, 0000000000000000, 0000000000000000, 0000000000000000,
                                                               0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                               0000000000000000, 0000000000000000, 8.3628688621153702e-01, 0000000000000000,
                                                               0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                               0000000000000000, 0000000000000000, 0000000000000000, 1.5022110482233481e-01,
                                                               0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                               0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                               1.3492008966128155e-02, 0000000000000000, 0000000000000000, 0000000000000000,
                                                               0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                               0000000000000000, 0000000000000000, 8.3628688621153702e-01, 0000000000000000,
                                                               0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                               0000000000000000, 0000000000000000, 0000000000000000, 1.5022110482233481e-01,
                                                               0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                               0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                               1.3492008966128155e-02};

        // Compare the computed values against the expected ones
        for(int i = 0; i < noCoord*noDOFsLocMaster; i++)
            CPPUNIT_ASSERT(fabs( BDisplacementsGCMaster1[i] - CorrectBOperatorDisplacementMatrixMaster1[i]) <= Tol);

        // Expected values for the bending rotation B-operator matrix of the master patch by MATLAB
        double CorrectBOperatorBendingRotationMatrixMaster1[] = { 0000000000000000, 0000000000000000, 0000000000000000, 1.2082071334501523e+00,
                                                                  -4.5960593897818907e+00, -6.9545526517053780e-01, -3.6665344714560595e+00,
                                                                  1.3947616860507358e+01, 2.1104913491297439e+00, 0000000000000000, 0000000000000000,
                                                                  0000000000000000, 2.1702864583147174e-01, -8.2558405583734096e-01, -1.2492370741540243e-01,
                                                                  1.8333340057510987e+00, -6.9740624256018142e+00, -1.0552841080112678e+00, 0000000000000000,
                                                                  0000000000000000, 0000000000000000, 1.9492283983185778e-02, -7.4149284794727072e-02,
                                                                  -1.1219939984624483e-02, 3.8847240244015213e-01, -1.4777617044915861e+00,
                                                                  -2.2360832854791138e-01};

        // Compare the computed values against the expected ones
        for(int i = 0; i < noDOFsLocMaster; i++)
            CPPUNIT_ASSERT(fabs( BOperatorOmegaTMaster1[i] + CorrectBOperatorBendingRotationMatrixMaster1[i]) <= TolRel100);

        // Expected values for the twisting rotation B-operator matrix of the master patch by MATLAB
        double CorrectBOperatorTwistingRotationMatrixMaster1[] = {0000000000000000, 0000000000000000, 0000000000000000, 2.1955203747055016e-16,
                                                                  -8.3518312003402813e-16, -1.2637619511630939e-16, -3.4909183754991653e+00,
                                                                  1.3279567496724303e+01, 2.0094050906558723e+00, 0000000000000000, 0000000000000000,
                                                                  0000000000000000, 3.9437841461593997e-17, -1.5002283676697640e-16, -2.2700788409604323e-17,
                                                                  2.9115929950373385e+00, -1.1075794831513168e+01, -1.6759400125789063e+00, 0000000000000000,
                                                                  0000000000000000, 0000000000000000, 3.5420835922746796e-18, -1.3474201652144170e-17,
                                                                  -2.0388562653882316e-18, 5.7932538046182791e-01, -2.2037726652111371e+00,
                                                                  -3.3346507807696679e-01};

        // Compare the computed values against the expected ones
        for(int i = 0; i < noDOFsLocMaster; i++)
            CPPUNIT_ASSERT(fabs( BOperatorOmegaNMaster1[i] - CorrectBOperatorTwistingRotationMatrixMaster1[i]) <= TolRel100);

        // Compute the B-operator matrices needed for the computation of the patch weak continuity contributions at the slave patch
        theMapper->computeIGAPatchContinuityConditionBOperatorMatrices(BDisplacementsGCSlave1, BOperatorOmegaTSlave1, BOperatorOmegaNSlave1, normalTrCurveVctSlave,
                                                                       patchSlave, tangentTrCurveVctSlave, uGPSlave, vGPSlave, uKnotSpanSlave, vKnotSpanSlave);

        // Expected values by MATLAB
        double CorrectBOperatorDisplacementMatrixSlave1[] = { 8.3628688621153702e-01, 0000000000000000, 0000000000000000, 0000000000000000,
                                                              0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                              0000000000000000, 1.5022110482233481e-01, 0000000000000000, 0000000000000000,
                                                              0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                              0000000000000000, 0000000000000000, 1.3492008966128155e-02, 0000000000000000,
                                                              0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                              0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                              8.3628688621153702e-01, 0000000000000000, 0000000000000000, 0000000000000000,
                                                              0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                              0000000000000000, 1.5022110482233481e-01, 0000000000000000, 0000000000000000,
                                                              0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                              0000000000000000, 0000000000000000, 1.3492008966128155e-02, 0000000000000000,
                                                              0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                              0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                              8.3628688621153702e-01, 0000000000000000, 0000000000000000, 0000000000000000,
                                                              0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                              0000000000000000, 1.5022110482233481e-01, 0000000000000000, 0000000000000000,
                                                              0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                              0000000000000000, 0000000000000000, 1.3492008966128155e-02, 0000000000000000,
                                                              0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                              0000000000000000};

        // Compare the computed values against the expected ones
        for(int i = 0; i < noCoord*noDOFsLocSlave; i++)
            CPPUNIT_ASSERT(fabs( BDisplacementsGCSlave1[i] - CorrectBOperatorDisplacementMatrixSlave1[i]) <= Tol);

        // Expected values for the bending rotation B-operator matrix of the master patch by MATLAB
        double CorrectBOperatorBendingRotationMatrixSlave1[] = { 1.6503648110537932e+00, -7.2203598335196046e+00, -9.1327757250282193e-01,
                                                                 5.3105615910724480e-01, -2.3233751318975284e+00, -2.9387543687534928e-01,
                                                                 0000000000000000, 0000000000000000, 0000000000000000, -1.9148024024100738e+00,
                                                                 8.3772765044963737e+00, 1.0596118374452201e+00, 9.5392913914013930e-02,
                                                                 -4.1734479517128364e-01, -5.2788455176606186e-02, 0000000000000000, 0000000000000000,
                                                                 0000000000000000, -3.7057913300368456e-01, 1.6212868022627205e+00, 2.0507078722413091e-01,
                                                                 8.5676513387059738e-03, -3.7483546170675550e-02, -4.7411601145735420e-03, 0000000000000000,
                                                                 0000000000000000, 0000000000000000};

        // Compare the computed values against the expected ones
        for(int i = 0; i < noDOFsLocSlave; i++)
            CPPUNIT_ASSERT(fabs( BOperatorOmegaTSlave1[i] + CorrectBOperatorBendingRotationMatrixSlave1[i]) <= TolRel100);

        // Expected values for the twisting rotation B-operator matrix of the master patch by MATLAB
        double CorrectBOperatorTwistingRotationMatrixSlave1[] = { 3.0685833881628630e+00, -1.3425077954461102e+01, -1.6980902458617186e+00, 8.9318362842006234e-17,
                                                                  -3.9076858349176899e-16, -4.9426924913763718e-17, 0000000000000000, 0000000000000000,
                                                                  0000000000000000, -2.5593454033096381e+00, 1.1197157517167645e+01, 1.4162885329810238e+00,
                                                                  1.6044139120524749e-17, -7.0193242665703890e-17, -8.8784930039403261e-18, 0000000000000000,
                                                                  0000000000000000, 0000000000000000, -5.0923798485322713e-01, 2.2279204372934602e+00,
                                                                  2.8180171288069539e-01, 1.4409937213811724e-18, -6.3043595673680579e-18, -7.9741596466449943e-19,
                                                                  0000000000000000, 0000000000000000, 0000000000000000};

        // Compare the computed values against the expected ones
        for(int i = 0; i < noDOFsLocMaster; i++)
            CPPUNIT_ASSERT(fabs( BOperatorOmegaNSlave1[i] - CorrectBOperatorTwistingRotationMatrixSlave1[i]) <= TolRel100);

        // Compute the B-operator matrices at the second Gauss Point
        iGP = 1;

        // Get the parametric coordinates of the Gauss Point on the master patch
        uGPMaster = trCurveMasterGPs[2*iGP];
        vGPMaster = trCurveMasterGPs[2*iGP + 1];

        std::cout << std::endl;
        std::cout << "uGPMaster = " << trCurveMasterGPs[2*iGP] << std::endl;
        std::cout << "vGPMaster = " << trCurveMasterGPs[2*iGP + 1] << std::endl;
        std::cout << std::endl;

        // Get the parametric coordinates of the Gauss Point on the slave patch
        uGPSlave = trCurveSlaveGPs[2*iGP];
        vGPSlave = trCurveSlaveGPs[2*iGP + 1];

        std::cout << std::endl;
        std::cout << "uGPSlave = " << trCurveSlaveGPs[2*iGP] << std::endl;
        std::cout << "vGPSlave = " << trCurveSlaveGPs[2*iGP + 1] << std::endl;
        std::cout << std::endl;

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
        theMapper->computeIGAPatchContinuityConditionBOperatorMatrices(BDisplacementsGCMaster2, BOperatorOmegaTMaster2, BOperatorOmegaNMaster2, normalTrCurveVctMaster,
                                                                       patchMaster, tangentTrCurveVctMaster, uGPMaster, vGPMaster, uKnotSpanMaster, vKnotSpanMaster);

        // Expected values for the displacement B-operator matrix of the master patch by MATLAB
        double CorrectBOperatorDisplacementMatrixMaster2[] = {  0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                                2.9289321881345248e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                                0000000000000000, 0000000000000000, 0000000000000000, 4.1421356237309503e-01, 0000000000000000, 0000000000000000,
                                                                0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                                2.9289321881345248e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                                0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 2.9289321881345248e-01, 0000000000000000,
                                                                0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                                0000000000000000, 4.1421356237309503e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                                0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 2.9289321881345248e-01, 0000000000000000,
                                                                0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                                0000000000000000, 0000000000000000, 2.9289321881345248e-01, 0000000000000000, 0000000000000000, 0000000000000000,
                                                                0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 4.1421356237309503e-01,
                                                                0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                                0000000000000000, 0000000000000000, 2.9289321881345248e-01};

        // Compare the computed values against the expected ones
        for(int i = 0; i < noCoord*noDOFsLocMaster; i++)
            CPPUNIT_ASSERT(fabs( BDisplacementsGCMaster2[i] - CorrectBOperatorDisplacementMatrixMaster2[i]) <= Tol);

        // Expected values for the bending rotation B-operator matrix of the master patch by MATLAB
        double CorrectBOperatorBendingRotationMatrixMaster2[] = { 0000000000000000, 0000000000000000, 0000000000000000, 3.4029929963324351e-01, -1.2839945950373786e+00,
                                                                  3.1732151795696960e-03, -3.5083249889435997e-01, 1.3237377591705441e+00, -3.2714349167858158e-03,
                                                                  0000000000000000, 0000000000000000, 0000000000000000, 4.8125588480739845e-01, -1.8158425703156105e+00,
                                                                  4.4876039432756405e-03, -4.8125588480739845e-01, 1.8158425703156105e+00, -4.4876039432756396e-03,
                                                                  0000000000000000, 0000000000000000, 0000000000000000, 3.4029929963324351e-01, -1.2839945950373786e+00,
                                                                  3.1732151795696960e-03, -3.2976610037212706e-01, 1.2442514309042130e+00, -3.0749954423535744e-03};

        // Compare the computed values against the expected ones
        for(int i = 0; i < noDOFsLocMaster; i++){
            CPPUNIT_ASSERT(fabs( BOperatorOmegaTMaster2[i] + CorrectBOperatorBendingRotationMatrixMaster2[i]) <= TolRel100);
        }

        // Expected values for the twisting rotation B-operator matrix of the master patch by MATLAB
        double CorrectBOperatorTwistingRotationMatrixMaster2[] = { 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                                   -2.5618326305314834e+00, 9.6661358237819819e+00, -2.3888518723020162e-02, 0000000000000000, 0000000000000000,
                                                                   0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 4.6270698887313225e-18,
                                                                   4.7381195660608742e-16, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                                   0000000000000000, 2.5618326305314834e+00, -9.6661358237819819e+00, 2.3888518723020637e-02};

        // Compare the computed values against the expected ones
        for(int i = 0; i < noDOFsLocMaster; i++)
            CPPUNIT_ASSERT(fabs( BOperatorOmegaNMaster2[i] - CorrectBOperatorTwistingRotationMatrixMaster2[i]) <= TolRel100);

        // Compute the B-operator matrices needed for the computation of the patch weak continuity contributions at the slave patch
        theMapper->computeIGAPatchContinuityConditionBOperatorMatrices(BDisplacementsGCSlave2, BOperatorOmegaTSlave2, BOperatorOmegaNSlave2, normalTrCurveVctSlave,
                                                                       patchSlave, tangentTrCurveVctSlave, uGPSlave, vGPSlave, uKnotSpanSlave, vKnotSpanSlave);

        // Expected values for the displacement B-operator matrix of the master patch by MATLAB
        double CorrectBOperatorDisplacementMatrixSlave2[] = { 2.9289321881345248e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                              0000000000000000, 0000000000000000, 0000000000000000, 4.1421356237309503e-01, 0000000000000000, 0000000000000000,
                                                              0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                              2.9289321881345248e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                              0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 2.9289321881345248e-01, 0000000000000000,
                                                              0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                              0000000000000000, 4.1421356237309503e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                              0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 2.9289321881345248e-01, 0000000000000000,
                                                              0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                              0000000000000000, 0000000000000000, 2.9289321881345248e-01, 0000000000000000, 0000000000000000, 0000000000000000,
                                                              0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 4.1421356237309503e-01,
                                                              0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                              0000000000000000, 0000000000000000, 2.9289321881345248e-01, 0000000000000000, 0000000000000000, 0000000000000000,
                                                              0000000000000000, 0000000000000000, 0000000000000000};

        // Compare the computed values against the expected ones
        for(int i = 0; i < noCoord*noDOFsLocSlave; i++)
            CPPUNIT_ASSERT(fabs( BDisplacementsGCSlave2[i] - CorrectBOperatorDisplacementMatrixSlave2[i]) <= Tol);

        // Expected values for the bending rotation B-operator matrix of the master patch by MATLAB
        double CorrectBOperatorBendingRotationMatrixSlave2[] = { -1.4678423016427308e-01, 6.4029113709351637e-01, -1.6591169834116850e-03, 1.5580898895687181e-01, -6.7965826163299092e-01,
                                                                 1.7611247438314361e-03, 0000000000000000, 0000000000000000, 0000000000000000, -2.2034718532244793e-01, 9.6118193138029717e-01,
                                                                 -2.4906064977572584e-03, 2.2034718532244793e-01, -9.6118193138029706e-01, 2.4906064977572593e-03, 0000000000000000,
                                                                 0000000000000000, 0000000000000000, -1.6483374774947077e-01, 7.1902538617246547e-01, -1.8631325042511861e-03, 1.5580898895687181e-01,
                                                                 -6.7965826163299092e-01, 1.7611247438314361e-03, 0000000000000000, 0000000000000000, 0000000000000000};

        // Compare the computed values against the expected ones
        for(int i = 0; i < noDOFsLocSlave; i++){
            CPPUNIT_ASSERT(fabs( BOperatorOmegaTSlave2[i] + CorrectBOperatorBendingRotationMatrixSlave2[i]) <= TolRel100);
        }

        // Expected values for the twisting rotation B-operator matrix of the master patch by MATLAB
        double CorrectBOperatorTwistingRotationMatrixSlave2[] = { 2.2344644649049124e+00, -9.7470129551917903e+00, 2.5256377598633858e-02, -4.5412764855402270e-20, 1.9809614980631200e-19,
                                                                  -5.1330506929087047e-22, 0000000000000000, 0000000000000000, 0000000000000000, 9.8963728881358436e-19, -1.2055642025601931e-18,
                                                                  2.3690669478064548e-16, -6.4223347963370136e-20, 2.8015026170997876e-19, -7.2592299062601025e-22, 0000000000000000, 0000000000000000,
                                                                  0000000000000000, -2.2344644649049124e+00, 9.7470129551917903e+00, -2.5256377598633618e-02, -4.5412764855402270e-20,
                                                                  1.9809614980631200e-19, -5.1330506929087047e-22, 0000000000000000, 0000000000000000, 0000000000000000};

        // Compare the computed values against the expected ones
        for(int i = 0; i < noDOFsLocSlave; i++)
            CPPUNIT_ASSERT(fabs( BOperatorOmegaNSlave2[i] - CorrectBOperatorTwistingRotationMatrixSlave2[i]) <= TolRel100);

        // Compute the B-operator matrices at the second Gauss Point
        iGP = 2;

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
        theMapper->computeIGAPatchContinuityConditionBOperatorMatrices(BDisplacementsGCMaster3, BOperatorOmegaTMaster3, BOperatorOmegaNMaster3, normalTrCurveVctMaster,
                                                                       patchMaster, tangentTrCurveVctMaster, uGPMaster, vGPMaster, uKnotSpanMaster, vKnotSpanMaster);

        // Expected values for the displacement B-operator matrix of the master patch by MATLAB
        double CorrectBOperatorDisplacementMatrixMaster3[] = {  0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                                1.3492008966128155e-02, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                                0000000000000000, 0000000000000000, 0000000000000000, 1.5022110482233481e-01, 0000000000000000, 0000000000000000,
                                                                0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                                8.3628688621153702e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                                0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 1.3492008966128155e-02, 0000000000000000,
                                                                0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                                1.5022110482233481e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                                0000000000000000, 0000000000000000, 0000000000000000, 8.3628688621153702e-01, 0000000000000000, 0000000000000000,
                                                                0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                                1.3492008966128155e-02, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                                0000000000000000, 0000000000000000, 0000000000000000, 1.5022110482233481e-01, 0000000000000000, 0000000000000000,
                                                                0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 8.3628688621153702e-01};

        // Compare the computed values against the expected ones
        for(int i = 0; i < noCoord*noDOFsLocMaster; i++)
            CPPUNIT_ASSERT(fabs( BDisplacementsGCMaster3[i] - CorrectBOperatorDisplacementMatrixMaster3[i]) <= Tol);

        // Expected values for the bending rotation B-operator matrix of the master patch by MATLAB
        double CorrectBOperatorBendingRotationMatrixMaster3[] = { 0000000000000000, 0000000000000000, 0000000000000000, 1.9297090468641944e-02, -7.3617664372221947e-02,
                                                                  1.1464583416185073e-02, 3.8378454811742241e-01, -1.4641234180077731e+00, 2.2801001906919766e-01, 0000000000000000,
                                                                  0000000000000000, 0000000000000000, 2.1485534565930742e-01, -8.1966495161681097e-01, 1.2764758691094841e-01,
                                                                  1.8109659197576788e+00, -6.9087659347868247e+00, 1.0759119300740561e+00, 0000000000000000, 0000000000000000,
                                                                  0000000000000000, 1.1961082846503659e+00, -4.5631075003413031e+00, 7.1061921103846160e-01, -3.6250111886534158e+00,
                                                                  1.3829279469124934e+01, -2.1536533305088490e+00};

        // Compare the computed values against the expected ones
        for(int i = 0; i < noDOFsLocMaster; i++){
            CPPUNIT_ASSERT(fabs( BOperatorOmegaTMaster3[i] + CorrectBOperatorBendingRotationMatrixMaster3[i]) <= TolRel100);
        }

        // Expected values for the twisting rotation B-operator matrix of the master patch by MATLAB
        double CorrectBOperatorTwistingRotationMatrixMaster3[] = { 0000000000000000, 0000000000000000, 0000000000000000, -2.5538947003129771e-18, 9.7430109059787522e-18,
                                                                   -1.5172929243126293e-18, -5.7934837759982427e-01, 2.2101919709627977e+00, -3.4419633430324637e-01, 0000000000000000,
                                                                   0000000000000000, 0000000000000000, -2.8435267456764674e-17, 1.0847946115857058e-16, -1.6893660536512676e-17,
                                                                   -2.9117085748274776e+00, 1.1108057194409538e+01, -1.7298735212946448e+00, 0000000000000000, 0000000000000000,
                                                                   0000000000000000, -1.5830026884793869e-16, 6.0390948993152508e-16, -9.4047682471140848e-17, 3.4910569524273001e+00,
                                                                   -1.3318249165372338e+01, 2.0740698555978900e+00};

        // Compare the computed values against the expected ones
        for(int i = 0; i < noDOFsLocMaster; i++)
            CPPUNIT_ASSERT(fabs( BOperatorOmegaNMaster3[i] - CorrectBOperatorTwistingRotationMatrixMaster3[i]) <= TolRel100);

        // Compute the B-operator matrices needed for the computation of the patch weak continuity contributions at the slave patch
        theMapper->computeIGAPatchContinuityConditionBOperatorMatrices(BDisplacementsGCSlave3, BOperatorOmegaTSlave3, BOperatorOmegaNSlave3, normalTrCurveVctSlave,
                                                                       patchSlave, tangentTrCurveVctSlave, uGPSlave, vGPSlave, uKnotSpanSlave, vKnotSpanSlave);

        // Expected values for the displacement B-operator matrix of the master patch by MATLAB
        double CorrectBOperatorDisplacementMatrixSlave3[] = { 1.3492008966128155e-02, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                              0000000000000000, 0000000000000000, 1.5022110482233481e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                              0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 8.3628688621153702e-01, 0000000000000000, 0000000000000000,
                                                              0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                              1.3492008966128155e-02, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                              0000000000000000, 0000000000000000, 1.5022110482233481e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                              0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 8.3628688621153702e-01, 0000000000000000, 0000000000000000,
                                                              0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                              1.3492008966128155e-02, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                              0000000000000000, 0000000000000000, 1.5022110482233481e-01, 0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000,
                                                              0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000, 8.3628688621153702e-01, 0000000000000000, 0000000000000000,
                                                              0000000000000000, 0000000000000000, 0000000000000000, 0000000000000000};

        // Compare the computed values against the expected ones
        for(int i = 0; i < noCoord*noDOFsLocSlave; i++)
            CPPUNIT_ASSERT(fabs( BDisplacementsGCSlave3[i] - CorrectBOperatorDisplacementMatrixSlave3[i]) <= Tol);

        // Expected values for the bending rotation B-operator matrix of the master patch by MATLAB
        double CorrectBOperatorBendingRotationMatrixSlave3[] = { -3.6774185675646665e-01, 1.6062869507087991e+00, -2.1193032717698823e-01, 8.5534159708889959e-03, -3.7361100471959334e-02,
                                                                 4.9293497922153658e-03, 0000000000000000, 0000000000000000, 0000000000000000, -1.9004557708238417e+00, 8.3011418172481051e+00,
                                                                 -1.0952362530839634e+00, 9.5234416192407967e-02, -4.1598147498760529e-01, 5.4883774069625468e-02, 0000000000000000,
                                                                 0000000000000000, 0000000000000000, 1.6342359994519522e+00, -7.1383007184755591e+00, 9.4181329561738636e-01, 5.3017379596506020e-01,
                                                                 -2.3157854740217849e+00, 3.0554016078172525e-01, 0000000000000000, 0000000000000000, 0000000000000000};

        // Compare the computed values against the expected ones
        for(int i = 0; i < noDOFsLocSlave; i++){
            CPPUNIT_ASSERT(fabs( BOperatorOmegaTSlave3[i] + CorrectBOperatorBendingRotationMatrixSlave3[i]) <= TolRel100);
        }

        // Expected values for the twisting rotation B-operator matrix of the master patch by MATLAB
        double CorrectBOperatorTwistingRotationMatrixSlave3[] = { 5.1139632877411267e-01, -2.2337659813749196e+00, 2.9471867094526083e-01, 1.1770085829963533e-18, -5.1411431497484843e-18,
                                                                  6.7831227123468878e-19, 0000000000000000, 0000000000000000, 0000000000000000, 2.5701928808289587e+00, -1.1226536249351257e+01,
                                                                  1.4812070155581185e+00, 1.3104907517254869e-17, -5.7241898218707404e-17, 7.5523829738947936e-18, 0000000000000000, 0000000000000000,
                                                                  0000000000000000, -3.0815892096030719e+00, 1.3460302230726183e+01, -1.7759256865033803e+00, 7.2955543195191520e-17,
                                                                  -3.1866793203775688e-16, 4.2044417448435991e-17, 0000000000000000, 0000000000000000, 0000000000000000};

        // Compare the computed values against the expected ones
        for(int i = 0; i < noDOFsLocSlave; i++)
            CPPUNIT_ASSERT(fabs( BOperatorOmegaNSlave3[i] - CorrectBOperatorTwistingRotationMatrixSlave3[i]) <= TolRel100);
    }

    void testComputePenaltyFactors() {

//        double alphaPrim, alphaSec;

//        std::vector<double> resultingPenaltyFactors;

//        for(int patchCounter = 0 ; patchCounter < theIGAMeshM->getIGAPatchCouplingData()->getNumPatches() ; patchCounter++) {

//            // get master patch and polynomial degrees
//            IGAPatchSurface* masterPatch = theIGAMeshM->getSurfacePatch(patchCounter);

//            // loop through all BReps where current patch is master
//            int* numOfBRepsPerPatch = theIGAMeshM->getIGAPatchCouplingData()->getNumBrepsPerPatch();
//            for(int BRepCounter = 0 ; BRepCounter < numOfBRepsPerPatch[patchCounter] ; BRepCounter++) {

//                int slaveID = theIGAMeshM->getIGAPatchCouplingData()->getSlaveID(patchCounter, BRepCounter);

//                // get slave patch and polynomial degrees
//                IGAPatchSurface* slavePatch = theIGAMeshM->getSurfacePatch(slaveID);

//                double* gausspoints_master = theIGAMeshM->getIGAPatchCouplingData()->getGPs_master(patchCounter, BRepCounter);
//                double* gausspoints_slave = theIGAMeshM->getIGAPatchCouplingData()->getGPs_slave(patchCounter, BRepCounter);
//                double* gausspoints_weight = theIGAMeshM->getIGAPatchCouplingData()->getGPs_weight(patchCounter, BRepCounter);
//                double* mappings = theIGAMeshM->getIGAPatchCouplingData()->getMappings(patchCounter, BRepCounter);

//                int numElemsPerBRep = theIGAMeshM->getIGAPatchCouplingData()->getNumElemsOfBRep(patchCounter, BRepCounter);
//                int numGPsPerElem = theIGAMeshM->getIGAPatchCouplingData()->getNumGPsOfElem(patchCounter, BRepCounter);

//                theMapperM->computePenaltyFactorsForPatchCoupling(alphaPrim, alphaSec, masterPatch,
//                        slavePatch, gausspoints_master, gausspoints_slave, gausspoints_weight,
//                        mappings, numElemsPerBRep, numGPsPerElem);

//                resultingPenaltyFactors.push_back(alphaPrim);
//                resultingPenaltyFactors.push_back(alphaSec);
//            }
//        }

//        double correctPenaltyFactors[] = {
//            1.273239485397344, 1.128379140802126,
//            1.273239485397344, 1.128379140802126,
//            1.273239485397344, 1.128379140802126,
//            1.273239485397344, 1.128379140802126,
//            0.954176621909422, 0.976819646561955,
//            0.954176621909422, 0.976819646561955,
//            0.954176621909422, 0.976819646561955,
//            0.954176621909422, 0.976819646561955};

//        for(int i = 0 ; i < resultingPenaltyFactors.size() ; i++) {
//            CPPUNIT_ASSERT(fabs(resultingPenaltyFactors[i] - correctPenaltyFactors[i])<=1e-6);
//        }
    }

// Make the tests
    CPPUNIT_TEST_SUITE (TestIGAMortarMapperWeakContinuityConditions);
    CPPUNIT_TEST (testIGAPatchContinuityConditions);
    CPPUNIT_TEST (testComputePenaltyFactors);
    CPPUNIT_TEST_SUITE_END();
}
;

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION (EMPIRE::TestIGAMortarMapperWeakContinuityConditions);
