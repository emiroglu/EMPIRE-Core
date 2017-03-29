/*
 * TestIGABarycentricMapper.cpp
 *
 *  Created on: May 8, 2013
 *      Author: Chenshen
 */

#include "cppunit/TestFixture.h"
#include "cppunit/TestAssert.h"
#include "cppunit/extensions/HelperMacros.h"
#include "IGABarycentricMapper.h"
#include "FEMesh.h"
#include "IGAPatchSurface.h"
#include "IGAMesh.h"
#include "DataField.h"
#include <iostream>
#include <math.h>
#include <stdlib.h>

using namespace std;

namespace EMPIRE {
/*
 * Test the projection of a point on a multipatch structure. The benchmark consists of a 3-patch structure perpendicular to each other, and two fluid elements matching the two structural patches.
 */
class TestIGABarycentricMapperCurve: public CppUnit::TestFixture {
private:
    IGABarycentricMapper* theMapper;
    IGAMesh* theIGAMesh;
    FEMesh* theFEMesh;
    double Tol;
    int a;
    int b;

public:
    void setUp() {

        Tol = 1e-13;

        // Provide an id for the basis
        int id_basis = 1;

        int numControlPoints = 10;

        int* controlPointID = new int[numControlPoints];
        for (int i = 0; i < numControlPoints; i++)
            controlPointID[i] = i + 1;

        theIGAMesh = new IGAMesh("IGAMesh", numControlPoints);

        // The polynomial degrees
        int p = 2;
        int q = 1;

        // Number of knots in both directions
        int uNoKnots = 6;
        int vNoKnots = 4;

        // The Control Point net
        int uNoControlPoints = uNoKnots - p - 1;
        int vNoControlPoints = vNoKnots - q - 1;

        // The knot vectors in each direction
        double* uKnotVector1 = new double[uNoKnots];
        double* vKnotVector1 = new double[vNoKnots];
        double* uKnotVector2 = new double[uNoKnots];
        double* vKnotVector2 = new double[vNoKnots];

        uKnotVector1[0] = 0;
        uKnotVector1[1] = 0;
        uKnotVector1[2] = 0;
        uKnotVector1[3] = 1;
        uKnotVector1[4] = 1;
        uKnotVector1[5] = 1;

        vKnotVector1[0] = 0;
        vKnotVector1[1] = 0;
        vKnotVector1[2] = 1;
        vKnotVector1[3] = 1;
        for (int i = 0; i < 6; i++) {
            uKnotVector2[i] = uKnotVector1[i];
            if (i < 5) {
                vKnotVector2[i] = vKnotVector1[i];
            }
        }

        double* controlPointNet1 = new double[6 * 4];

        controlPointNet1[0] = 0.0; // 0
        controlPointNet1[1] = 0.0;
        controlPointNet1[2] = 0.0;
        controlPointNet1[3] = 1.0;

        controlPointNet1[4] = 2.0; // 1
        controlPointNet1[5] = 1.0;
        controlPointNet1[6] = 0.0;
        controlPointNet1[7] = 1.0;

        controlPointNet1[8] = 5.0; // 2
        controlPointNet1[9] = 2.0;
        controlPointNet1[10] = 0.0;
        controlPointNet1[11] = 1.0;

        controlPointNet1[12] = 0.0; // 5
        controlPointNet1[13] = 0.0;
        controlPointNet1[14] = 1.0;
        controlPointNet1[15] = 1.0;

        controlPointNet1[16] = 2.0; // 6
        controlPointNet1[17] = 1.0;
        controlPointNet1[18] = 1.0;
        controlPointNet1[19] = 1.0;

        controlPointNet1[20] = 5.0; // 7
        controlPointNet1[21] = 2.0;
        controlPointNet1[22] = 1.0;
        controlPointNet1[23] = 1.0;

        double* controlPointNet2 = new double[6 * 4];

        controlPointNet2[0] = 5.0; // 2
        controlPointNet2[1] = 2.0;
        controlPointNet2[2] = 0.0;
        controlPointNet2[3] = 1.0;

        controlPointNet2[4] = 6.0; // 3
        controlPointNet2[5] = 2.0;
        controlPointNet2[6] = 0.0;
        controlPointNet2[7] = 1.0;

        controlPointNet2[8] = 6.0; // 4
        controlPointNet2[9] = 0.0;
        controlPointNet2[10] = 0.0;
        controlPointNet2[11] = 1.0;

        controlPointNet2[12] = 5.0; // 7
        controlPointNet2[13] = 2.0;
        controlPointNet2[14] = 1.0;
        controlPointNet2[15] = 1.0;

        controlPointNet2[16] = 6.0; // 8
        controlPointNet2[17] = 2.0;
        controlPointNet2[18] = 1.0;
        controlPointNet2[19] = 1.0;

        controlPointNet2[20] = 6.0; // 9
        controlPointNet2[21] = 0.0;
        controlPointNet2[22] = 1.0;
        controlPointNet2[23] = 1.0;

        int* dofIndexNet1 = new int[6];
        int* dofIndexNet2 = new int[6];
        dofIndexNet1[0] = 0;
        dofIndexNet1[1] = 1;
        dofIndexNet1[2] = 2;
        dofIndexNet1[3] = 5;
        dofIndexNet1[4] = 6;
        dofIndexNet1[5] = 7;

        dofIndexNet2[0] = 2;
        dofIndexNet2[1] = 3;
        dofIndexNet2[2] = 4;
        dofIndexNet2[3] = 7;
        dofIndexNet2[4] = 8;
        dofIndexNet2[5] = 9;

        theIGAMesh->addPatch(p, uNoKnots, uKnotVector1, q, vNoKnots, vKnotVector1, uNoControlPoints,
                vNoControlPoints, controlPointNet1, dofIndexNet1);
        theIGAMesh->addPatch(p, uNoKnots, uKnotVector2, q, vNoKnots, vKnotVector2, uNoControlPoints,
                vNoControlPoints, controlPointNet2, dofIndexNet2);

        // std::vector<IGAControlPoint*> igaCpNet;
        // igaCpNet.resize(theIGAMesh->getNumNodes());
        // std::vector<IGAPatchSurface*> patches = theIGAMesh->getSurfacePatches();
        // for (int i = 0; i < theIGAMesh->getNumPatches(); ++i) {
        //     IGAControlPoint** cpNet = patches[i]->getControlPointNet();
        //     int num = patches[i]->getNoControlPoints();
        //     for (int j = 0; j < num; ++j) {
        //         igaCpNet[ cpNet[j]->getDofIndex() ] = cpNet[j];
        //     }
        // }

        int numNodes = 16;
        int numElems = 4;
        theFEMesh = new FEMesh("Fluid", numNodes, numElems);
        theFEMesh->numNodesPerElem[0] = 4;
        theFEMesh->numNodesPerElem[1] = 4;
        theFEMesh->numNodesPerElem[2] = 4;
        theFEMesh->numNodesPerElem[3] = 4;
        theFEMesh->initElems();
        theFEMesh->elems[0] = 1;
        theFEMesh->elems[1] = 2;
        theFEMesh->elems[2] = 3;
        theFEMesh->elems[3] = 4;
        theFEMesh->elems[4] = 5;
        theFEMesh->elems[5] = 6;
        theFEMesh->elems[6] = 7;
        theFEMesh->elems[7] = 8;
        theFEMesh->elems[8] = 9;
        theFEMesh->elems[9] = 10;
        theFEMesh->elems[10] = 11;
        theFEMesh->elems[11] = 12;
        theFEMesh->elems[12] = 13;
        theFEMesh->elems[13] = 14;
        theFEMesh->elems[14] = 15;
        theFEMesh->elems[15] = 16;
        // theFEMesh->elems[4] = 4;
        // theFEMesh->elems[5] = 5;
        for (int i = 0; i < numNodes; i++)
            theFEMesh->nodeIDs[i] = i + 1;

        theFEMesh->nodes[0 * 3 + 0] = 0.0;
        theFEMesh->nodes[0 * 3 + 1] = 0.0;
        theFEMesh->nodes[0 * 3 + 2] = 0.0;

        theFEMesh->nodes[1 * 3 + 0] = 1.0;
        theFEMesh->nodes[1 * 3 + 1] = 1.0;
        theFEMesh->nodes[1 * 3 + 2] = 0.0;

        theFEMesh->nodes[2 * 3 + 0] = 1.0;
        theFEMesh->nodes[2 * 3 + 1] = 1.0;
        theFEMesh->nodes[2 * 3 + 2] = 1.0;

        theFEMesh->nodes[3 * 3 + 0] = 0.0;
        theFEMesh->nodes[3 * 3 + 1] = 0.0;
        theFEMesh->nodes[3 * 3 + 2] = 1.0;

        theFEMesh->nodes[4 * 3 + 0] = 1.0;
        theFEMesh->nodes[4 * 3 + 1] = 1.0;
        theFEMesh->nodes[4 * 3 + 2] = 0.0;

        theFEMesh->nodes[5 * 3 + 0] = 3.0;
        theFEMesh->nodes[5 * 3 + 1] = 2.5;
        theFEMesh->nodes[5 * 3 + 2] = 0.0;

        theFEMesh->nodes[6 * 3 + 0] = 3.0;
        theFEMesh->nodes[6 * 3 + 1] = 2.5;
        theFEMesh->nodes[6 * 3 + 2] = 1.0;

        theFEMesh->nodes[7 * 3 + 0] = 1.0;
        theFEMesh->nodes[7 * 3 + 1] = 1.0;
        theFEMesh->nodes[7 * 3 + 2] = 1.0;

        theFEMesh->nodes[8 * 3 + 0] = 3.0;
        theFEMesh->nodes[8 * 3 + 1] = 2.5;
        theFEMesh->nodes[8 * 3 + 2] = 0.0;

        theFEMesh->nodes[9 * 3 + 0] = 7.0;
        theFEMesh->nodes[9 * 3 + 1] = 2.5;
        theFEMesh->nodes[9 * 3 + 2] = 0.0;

        theFEMesh->nodes[10 * 3 + 0] = 7.0;
        theFEMesh->nodes[10 * 3 + 1] = 2.5;
        theFEMesh->nodes[10 * 3 + 2] = 1.0;

        theFEMesh->nodes[11 * 3 + 0] = 3.0;
        theFEMesh->nodes[11 * 3 + 1] = 2.5;
        theFEMesh->nodes[11 * 3 + 2] = 1.0;

        theFEMesh->nodes[12 * 3 + 0] = 7.0;
        theFEMesh->nodes[12 * 3 + 1] = 2.5;
        theFEMesh->nodes[12 * 3 + 2] = 0.0;

        theFEMesh->nodes[13 * 3 + 0] = 7.0;
        theFEMesh->nodes[13 * 3 + 1] = 0.0;
        theFEMesh->nodes[13 * 3 + 2] = 0.0;

        theFEMesh->nodes[14 * 3 + 0] = 7.0;
        theFEMesh->nodes[14 * 3 + 1] = 0.0;
        theFEMesh->nodes[14 * 3 + 2] = 1.0;

        theFEMesh->nodes[15 * 3 + 0] = 7.0;
        theFEMesh->nodes[15 * 3 + 1] = 2.5;
        theFEMesh->nodes[15 * 3 + 2] = 1.0;

        bool isMappingIGA2FEM = false;
        theMapper = new IGABarycentricMapper("Test IGA Barycentric Mapper for Curve", theIGAMesh, theFEMesh, isMappingIGA2FEM);
        theMapper->setParametersProjection(7.5, 10, 1e-3);
        theMapper->buildCouplingMatrices();
	
	delete controlPointID;
	delete uKnotVector1;
	delete uKnotVector2;
	delete vKnotVector1;
	delete vKnotVector2;
	delete controlPointNet1;
	delete controlPointNet2;
	delete dofIndexNet1;
	delete dofIndexNet2;
    }

    void tearDown() {

        delete theMapper;
 	// delete theIGAMesh;
 	delete theFEMesh;

    }
    /***********************************************************************************************
     * \brief Test case: Test the constructor
     ***********/

    void testMapping() {

        int nS = theIGAMesh->getNumNodes();
        int nF = theFEMesh->numNodes;
        double fieldS[nS];
        double fieldF[nF];

        for (int i = 0; i < nF; i++)
            fieldF[i] = 1;
        theMapper->consistentMapping(fieldF, fieldS);
        for (int i = 0; i < nS; i++) 
            CPPUNIT_ASSERT(fabs(fieldS[i] - 1.0) < Tol);

        /* cout << "Barycentric:\n";
        for (int i = 0; i < nS; ++i)
        {
            cout << "[" << i << "] = " << fieldS[i] << endl;
        } */
    }

// Make the tests
    CPPUNIT_TEST_SUITE (TestIGABarycentricMapperCurve);

    CPPUNIT_TEST (testMapping);
	// CPPUNIT_TEST(testMappingPrint);
    CPPUNIT_TEST_SUITE_END();
};

} /* namespace EMPIRE */
CPPUNIT_TEST_SUITE_REGISTRATION (EMPIRE::TestIGABarycentricMapperCurve);
