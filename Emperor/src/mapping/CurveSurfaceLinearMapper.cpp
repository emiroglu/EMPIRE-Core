#include "CurveSurfaceLinearMapper.h"
#include "KinematicMotion.h"

#include <assert.h>
#include <map>
#include <iostream>
#include <math.h>

using namespace std;

namespace EMPIRE {

CurveSurfaceLinearMapper::CurveSurfaceLinearMapper(int _curveNumNodes, int _curveNumElements,
        const double *_curveNodeCoors, const int *_curveNodeIDs, const int *_curveElems,
        int _surfaceNumNodes, const double *_surfaceNodeCoors, int _surfaceNumSections,
        int _surfaceNumRootSectionNodes, int _surfaceNumNormalSectionNodes,
        int _surfaceNumTipSectionNodes, const double *rotation_O_Q, const double *translation_O_Q) :
        curveNumNodes(_curveNumNodes), curveNumElements(_curveNumElements), curveElems(_curveElems), surfaceNumNodes(
                _surfaceNumNodes), surfaceNodeCoors(_surfaceNodeCoors), surfaceNumSections(
                _surfaceNumSections), surfaceNumRootSectionNodes(_surfaceNumRootSectionNodes), surfaceNumNormalSectionNodes(
                _surfaceNumNormalSectionNodes), surfaceNumTipSectionNodes(
                _surfaceNumTipSectionNodes) {
    /*
     * Coordinate systems:
     *   O --- global
     *   Q --- beam root (x is the direction orthogonal to root section, y-z is the root plane)
     *         origin is an arbitrary point in the root section
     *   P --- origin is the cross point of the section with the beam element, orientation is the same as the global system
     */
    KinematicMotion *KM_O_Q = new KinematicMotion();
    KM_O_Q->addRotation(rotation_O_Q);
    KM_O_Q->addTranslation(translation_O_Q);

    // surface node coordinates in Q
    double *surfaceNodeCoorsInQ = new double[surfaceNumNodes * 3];
    for (int i = 0; i < surfaceNumNodes * 3; i++)
        surfaceNodeCoorsInQ[i] = surfaceNodeCoors[i];
    KinematicMotion *KM_Q_O = KM_O_Q->newInverse();
    for (int i = 0; i < surfaceNumNodes; i++) {
        KM_Q_O->move(&surfaceNodeCoorsInQ[i * 3]);
    }

    assert(
            surfaceNumNodes
                    == surfaceNumRootSectionNodes + surfaceNumTipSectionNodes
                            + (surfaceNumSections - 2) * surfaceNumNormalSectionNodes);

    // sort surface nodes according to x in Q system
    sortedPosToUnsortedPos = new int[surfaceNumNodes];

    multimap<double, int> *mapCoorX2Pos = new multimap<double, int>; // sort by x coordinate
    for (int i = 0; i < surfaceNumNodes; i++) {
        mapCoorX2Pos->insert(pair<double, int>(surfaceNodeCoorsInQ[i * 3 + 0], i));
    }

    int count = 0;
    for (multimap<double, int>::iterator it = mapCoorX2Pos->begin(); it != mapCoorX2Pos->end();
            it++) {
        sortedPosToUnsortedPos[count] = it->second;
        count++;
    }

    delete mapCoorX2Pos;

    // curve node coordinates in Q
    double *curveNodeCoorsInQ = new double[_curveNumNodes * 3];
    for (int i = 0; i < _curveNumNodes * 3; i++)
        curveNodeCoorsInQ[i] = _curveNodeCoors[i];
    for (int i = 0; i < _curveNumNodes; i++) {
        KM_Q_O->move(&curveNodeCoorsInQ[i * 3]);
    }

    // map curve node ID to curve node position
    curveNodeIDToPos = new map<int, int>;
    for (int i = 0; i < _curveNumNodes; i++) {
        curveNodeIDToPos->insert(pair<int, int>(_curveNodeIDs[i], i));
    }

    double *curveElementLength = new double[curveNumElements];
    for (int i = 0; i < curveNumElements; i++) {
        const double *node1 = &_curveNodeCoors[curveNodeIDToPos->at(curveElems[i * 2 + 0]) * 3];
        const double *node2 = &_curveNodeCoors[curveNodeIDToPos->at(curveElems[i * 2 + 1]) * 3];
        double tmp[3];
        for (int j = 0; j < 3; j++)
            tmp[j] = node2[j] - node1[j];
        curveElementLength[i] = sqrt(tmp[0] * tmp[0] + tmp[1] * tmp[1] + tmp[2] * tmp[2]);
    }

    // map the x coordinate of the right node to an element
    map<double, int> *rightXToElemPos = new map<double, int>;

    for (int i = 0; i < curveNumElements; i++) {
        int leftID = curveElems[i * 2 + 0];
        int rightID = curveElems[i * 2 + 1];
        double leftX = curveNodeCoorsInQ[curveNodeIDToPos->at(leftID) * 3 + 0];
        double rightX = curveNodeCoorsInQ[curveNodeIDToPos->at(rightID) * 3 + 0];

        if (leftX > rightX) {
            rightX = leftX;
        }
        rightXToElemPos->insert(pair<double, int>(rightX, i));
    }

    // relate a section to a curve element: section center P, shape function + derivative, section to curve element
    sectionP = new double[surfaceNumSections * 3]; // cross point between section and beam
    sectionToCurveElem = new int[surfaceNumSections];
    shapeFuncOfSection = new double[10 * surfaceNumSections];
    for (int i = 0; i < surfaceNumSections; i++) {
        // compute the x of a section in the Q system
        double sectionX = 0.0;
        if (i == 0) { // root
            for (int j = 0; j < surfaceNumRootSectionNodes; j++) {
                sectionX += surfaceNodeCoorsInQ[sortedPosToUnsortedPos[j] * 3 + 0];
            }
            sectionX /= (double) surfaceNumRootSectionNodes;
        } else if (i == surfaceNumSections - 1) { // tip
            for (int j = 0; j < surfaceNumTipSectionNodes; j++) {
                int pos = sortedPosToUnsortedPos[surfaceNumNodes - j - 1];
                sectionX += surfaceNodeCoorsInQ[pos * 3 + 0];
            }
            sectionX /= (double) surfaceNumTipSectionNodes;
        } else { // normal
            for (int j = 0; j < surfaceNumNormalSectionNodes; j++) {
                int pos = sortedPosToUnsortedPos[surfaceNumRootSectionNodes
                        + (i - 1) * surfaceNumNormalSectionNodes + j];
                sectionX += surfaceNodeCoorsInQ[pos * 3 + 0];
            }
            sectionX /= (double) surfaceNumNormalSectionNodes;
        }

        // use section x and rightXToElemPos to relate a section to the corresponding curve element
        map<double, int>::iterator it = rightXToElemPos->lower_bound(sectionX); // find the first one bigger than section x
        if (it != rightXToElemPos->end()) {
            sectionToCurveElem[i] = it->second;
        } else {
            sectionToCurveElem[i] = rightXToElemPos->rbegin()->second;
            /*cout << "sectionX: " << sectionX << endl;
             cout << "rightX: " << rightXToElemPos->rbegin()->first << endl;
             cout << "rightX - sectionX:" << rightXToElemPos->rbegin()->first - sectionX << endl;*/
        }

        // compute the local coordinate xi of the section in the beam/curve
        int node1ID = curveElems[sectionToCurveElem[i] * 2 + 0];
        int node2ID = curveElems[sectionToCurveElem[i] * 2 + 1];
        double node1X = curveNodeCoorsInQ[curveNodeIDToPos->at(node1ID) * 3 + 0];
        double node2X = curveNodeCoorsInQ[curveNodeIDToPos->at(node2ID) * 3 + 0];

        double length = curveElementLength[sectionToCurveElem[i]];

        double diff = node2X - node1X; // can be negative
        double xi = 2.0 * (sectionX - node1X) / diff - 1.0;

        if (diff < 0) // node2 is on the left
            xi = -xi;
        //cout << "xi:" << xi <<endl;
        //cout << "sectionX:" << sectionX <<endl;
        //cout << "node1X:" << node1X <<endl;
        //cout << "diff:" << diff <<endl;
        // compute the shape functions and their derivatives
        double linearShapeFunc1 = 0.5 * (1.0 - xi);
        double linearShapeFunc2 = 0.5 * (1.0 + xi);

        double cubicShapeFuncDisp1 = 1.0 / 4.0 * (1.0 - xi) * (1.0 - xi) * (2.0 + xi);
        double cubicShapeFuncRot1 = 1.0 / 8.0 * length * (1.0 - xi) * (1.0 - xi) * (1.0 + xi);
        double cubicShapeFuncDisp2 = 1.0 / 4.0 * (1.0 + xi) * (1.0 + xi) * (2.0 - xi);
        double cubicShapeFuncRot2 = -1.0 / 8.0 * length * (1.0 + xi) * (1.0 + xi) * (1.0 - xi);
        double Dxi_Dx = 2.0 / length;
        double cubicShapeFuncDispDeriv1 = Dxi_Dx * (-3.0) / 4.0 * (1.0 - xi) * (1.0 + xi);
        double cubicShapeFuncRotDeriv1 = Dxi_Dx * (-1.0) / 8.0 * length * (1.0 - xi)
                * (1.0 + 3 * xi);
        double cubicShapeFuncDispDeriv2 = Dxi_Dx * 3.0 / 4.0 * (1.0 - xi) * (1.0 + xi);
        double cubicShapeFuncRotDeriv2 = Dxi_Dx * (-1.0) / 8.0 * length * (1.0 + xi)
                * (1.0 - 3 * xi);

        if (diff > 0) {
            shapeFuncOfSection[i * 10 + 0] = linearShapeFunc1;
            shapeFuncOfSection[i * 10 + 1] = cubicShapeFuncDisp1;
            shapeFuncOfSection[i * 10 + 2] = cubicShapeFuncRot1;
            shapeFuncOfSection[i * 10 + 3] = cubicShapeFuncDispDeriv1;
            shapeFuncOfSection[i * 10 + 4] = cubicShapeFuncRotDeriv1;
            shapeFuncOfSection[i * 10 + 5 + 0] = linearShapeFunc2;
            shapeFuncOfSection[i * 10 + 5 + 1] = cubicShapeFuncDisp2;
            shapeFuncOfSection[i * 10 + 5 + 2] = cubicShapeFuncRot2;
            shapeFuncOfSection[i * 10 + 5 + 3] = cubicShapeFuncDispDeriv2;
            shapeFuncOfSection[i * 10 + 5 + 4] = cubicShapeFuncRotDeriv2;
        } else {
            shapeFuncOfSection[i * 10 + 0] = linearShapeFunc2;
            shapeFuncOfSection[i * 10 + 1] = cubicShapeFuncDisp2;
            shapeFuncOfSection[i * 10 + 2] = cubicShapeFuncRot2;
            shapeFuncOfSection[i * 10 + 3] = cubicShapeFuncDispDeriv2;
            shapeFuncOfSection[i * 10 + 4] = cubicShapeFuncRotDeriv2;
            shapeFuncOfSection[i * 10 + 5 + 0] = linearShapeFunc1;
            shapeFuncOfSection[i * 10 + 5 + 1] = cubicShapeFuncDisp1;
            shapeFuncOfSection[i * 10 + 5 + 2] = cubicShapeFuncRot1;
            shapeFuncOfSection[i * 10 + 5 + 3] = cubicShapeFuncDispDeriv1;
            shapeFuncOfSection[i * 10 + 5 + 4] = cubicShapeFuncRotDeriv1;
        }

        //for (int j=0; j<10; j++)
        // cout << "shapeFunction " << j << ": " << shapeFuncOfSection[i*10 +j] << endl;

        // compute P, which is the cross point between section and beam/curve
        for (int j = 0; j < 3; j++) {
            sectionP[i * 3 + j] = shapeFuncOfSection[i * 10 + 0]
                    * _curveNodeCoors[curveNodeIDToPos->at(node1ID) * 3 + j]
                    + shapeFuncOfSection[i * 10 + 5 + 0]
                            * _curveNodeCoors[curveNodeIDToPos->at(node2ID) * 3 + j];
        }
    }

    // transformation from local beam to global system
    ROT_O_ELEM = new KinematicMotion*[curveNumElements];
    for (int i = 0; i < curveNumElements; i++) {
        // For the definition of local axes, see carat ElementBeam1::calc_transformation_matrix, or
        // carat.st.bv.tum.de/caratuserswiki/index.php/Users:General_FEM_Analysis/Elements_Reference/Beam1
        ROT_O_ELEM[i] = new KinematicMotion;
        int node1ID = curveElems[i * 2 + 0];
        int node2ID = curveElems[i * 2 + 1];
        const double *node1 = &_curveNodeCoors[curveNodeIDToPos->at(node1ID) * 3];
        const double *node2 = &_curveNodeCoors[curveNodeIDToPos->at(node2ID) * 3];
        double tmp[3];
        for (int j = 0; j < 3; j++)
            tmp[j] = node2[j] - node1[j];
        //double length = sqrt(tmp[0] * tmp[0] + tmp[1] * tmp[1] + tmp[2] * tmp[2]);
        double length = curveElementLength[i];
        // check if local x-axis is parallel to z-axis
        if (fabs(fabs(node2[2] - node1[2]) / length - 1.0) < 1E-6) {
            // if true, x_loc -> z_gl; y_loc -> y_gl; z_loc -> -x_gl
            double x_loc[] = { 0.0, 0.0, 1.0 };
            double y_loc[] = { 0.0, 1.0, 0.0 };
            double z_loc[] = { -1.0, 0.0, 0.0 };
            if (node1[2] > node2[2]) {
                x_loc[2] = -1.0;
                z_loc[0] = 1.0;
            }
            ROT_O_ELEM[i]->addRotation(x_loc, y_loc, z_loc, true);
        } else {
            double x_loc[3];
            double y_loc[3];
            double z_loc[3];
            // mapping to global coordinates through direction cosines
            x_loc[0] = tmp[0] / length; //direction cosine x
            x_loc[1] = tmp[1] / length; //direction cosine y
            x_loc[2] = tmp[2] / length; //direction cosine z

            // using the projection L_x of x_loc to x_gl-y_gl-plane as perpendicular to y_loc
            double lengthXY = sqrt(tmp[0] * tmp[0] + tmp[1] * tmp[1]);
            y_loc[0] = -tmp[1] / lengthXY;
            y_loc[1] = tmp[0] / lengthXY;
            y_loc[2] = 0;

            // z-direction via cross product of x and y vectors
            z_loc[0] = x_loc[1] * y_loc[2] - x_loc[2] * y_loc[1];
            z_loc[1] = x_loc[2] * y_loc[0] - x_loc[0] * y_loc[2];
            z_loc[2] = x_loc[0] * y_loc[1] - x_loc[1] * y_loc[0];

            ROT_O_ELEM[i]->addRotation(x_loc, y_loc, z_loc, true);
        }
    }

    delete KM_O_Q, KM_Q_O;
    delete[] surfaceNodeCoorsInQ;
    delete[] curveNodeCoorsInQ;
    delete rightXToElemPos;
    delete[] curveElementLength;
}

CurveSurfaceLinearMapper::~CurveSurfaceLinearMapper() {
    delete[] sortedPosToUnsortedPos;
    delete curveNodeIDToPos;
    delete[] sectionToCurveElem;
    delete[] shapeFuncOfSection;
    for (int i = 0; i < curveNumElements; i++) {
        delete ROT_O_ELEM[i];
    }
    delete[] ROT_O_ELEM;
    delete[] sectionP;
}

void CurveSurfaceLinearMapper::consistentMapping(const double *curveDispRot, double *surfaceDisp) {
    double *node1DispLocal = new double[curveNumElements * 3];
    double *node1RotLocal = new double[curveNumElements * 3];
    double *node2DispLocal = new double[curveNumElements * 3];
    double *node2RotLocal = new double[curveNumElements * 3];
    // get displacements and rotations in the curve element local system
    for (int i = 0; i < curveNumElements; i++) {
        int node1ID = curveElems[i * 2 + 0];
        int node2ID = curveElems[i * 2 + 1];
        for (int j = 0; j < 3; j++) {
            node1DispLocal[i * 3 + j] = curveDispRot[curveNodeIDToPos->at(node1ID) * 6 + 0 + j];
            node1RotLocal[i * 3 + j] = curveDispRot[curveNodeIDToPos->at(node1ID) * 6 + 3 + j];
            node2DispLocal[i * 3 + j] = curveDispRot[curveNodeIDToPos->at(node2ID) * 6 + 0 + j];
            node2RotLocal[i * 3 + j] = curveDispRot[curveNodeIDToPos->at(node2ID) * 6 + 3 + j];
        }
        KinematicMotion *ROT_ELEM_O = ROT_O_ELEM[i]->newInverse();
        ROT_ELEM_O->move(&node1DispLocal[i * 3]);
        ROT_ELEM_O->move(&node1RotLocal[i * 3]);
        ROT_ELEM_O->move(&node2DispLocal[i * 3]);
        ROT_ELEM_O->move(&node2RotLocal[i * 3]);

        delete ROT_ELEM_O;
        //cout << "node1RotLocal" << endl;
        //cout << node1RotLocal[i * 3 + 0] <<" "<< node1RotLocal[i * 3 + 1] <<" "<< node1RotLocal[i * 3 + 2] << endl;
        //cout << "node2RotLocal" << endl;
        //cout << node2RotLocal[i * 3 + 0] <<" "<< node2RotLocal[i * 3 + 1] <<" "<< node2RotLocal[i * 3 + 2] << endl;
    }

    // compute displacements of nodes through sections
    for (int i = 0; i < surfaceNumSections; i++) {
        // interpolate section rot and disp on the element local system
        int elem = sectionToCurveElem[i];
        double sectionDisp[3];
        double sectionRot[3];

        double linearShapeFunc1 = shapeFuncOfSection[i * 10 + 0];
        double cubicShapeFuncDisp1 = shapeFuncOfSection[i * 10 + 1];
        double cubicShapeFuncRot1 = shapeFuncOfSection[i * 10 + 2];
        double cubicShapeFuncDispDeriv1 = shapeFuncOfSection[i * 10 + 3];
        double cubicShapeFuncRotDeriv1 = shapeFuncOfSection[i * 10 + 4];
        double linearShapeFunc2 = shapeFuncOfSection[i * 10 + 5 + 0];
        double cubicShapeFuncDisp2 = shapeFuncOfSection[i * 10 + 5 + 1];
        double cubicShapeFuncRot2 = shapeFuncOfSection[i * 10 + 5 + 2];
        double cubicShapeFuncDispDeriv2 = shapeFuncOfSection[i * 10 + 5 + 3];
        double cubicShapeFuncRotDeriv2 = shapeFuncOfSection[i * 10 + 5 + 4];
        //disp_x = NL_1*disp_x1 + NL_2*disp_x2
        sectionDisp[0] = linearShapeFunc1 * node1DispLocal[elem * 3 + 0]
                + linearShapeFunc2 * node2DispLocal[elem * 3 + 0];
        //disp_y = NCd_1*disp_y1 + NCd_2*disp_y2 + NCr_1*rot_z1 + NCr_2*rot_z2
        sectionDisp[1] = cubicShapeFuncDisp1 * node1DispLocal[elem * 3 + 1]
                + cubicShapeFuncDisp2 * node2DispLocal[elem * 3 + 1]
                + cubicShapeFuncRot1 * node1RotLocal[elem * 3 + 2]
                + cubicShapeFuncRot2 * node2RotLocal[elem * 3 + 2];
        //disp_z = NCd_1*disp_z1 + NCd_2*disp_z2 + NCr_1*(-rot_y1) + NCr_2*(-rot_y2), since disp_z' = -rot_y
        sectionDisp[2] = cubicShapeFuncDisp1 * node1DispLocal[elem * 3 + 2]
                + cubicShapeFuncDisp2 * node2DispLocal[elem * 3 + 2]
                + cubicShapeFuncRot1 * (-node1RotLocal[elem * 3 + 1])
                + cubicShapeFuncRot2 * (-node2RotLocal[elem * 3 + 1]);
        //rot_x = NL_1*rot_x1 + NL_2*rot_x2
        sectionRot[0] = linearShapeFunc1 * node1RotLocal[elem * 3 + 0]
                + linearShapeFunc2 * node2RotLocal[elem * 3 + 0];
        //rot_y = - disp_z' = - [NCd'_1*disp_z1 + NCd'_2*disp_z2 + NCr'_1*(-rot_y1) + NCr'_2*(-rot_y2)]
        sectionRot[1] = -(cubicShapeFuncDispDeriv1 * node1DispLocal[elem * 3 + 2]
                + cubicShapeFuncDispDeriv2 * node2DispLocal[elem * 3 + 2]
                + cubicShapeFuncRotDeriv1 * (-node1RotLocal[elem * 3 + 1])
                + cubicShapeFuncRotDeriv2 * (-node2RotLocal[elem * 3 + 1]));
        //rot_z = disp_y' = NCd'_1*disp_y1 + NCd'_2*disp_y2 + NCr'_1*rot_z1 + NCr'_2*rot_z2
        sectionRot[2] = cubicShapeFuncDispDeriv1 * node1DispLocal[elem * 3 + 1]
                + cubicShapeFuncDispDeriv2 * node2DispLocal[elem * 3 + 1]
                + cubicShapeFuncRotDeriv1 * node1RotLocal[elem * 3 + 2]
                + cubicShapeFuncRotDeriv2 * node2RotLocal[elem * 3 + 2];

        // Compute KM_Oo_Od
        KinematicMotion KM_Eo_Ed;
        double localX[] = { 1.0, 0.0, 0.0 };
        double localY[] = { 0.0, 1.0, 0.0 };
        double localZ[] = { 0.0, 0.0, 1.0 };
        // bending first, then torsion
        KM_Eo_Ed.addRotation(localY, true, sectionRot[1]);
        KM_Eo_Ed.addRotation(localZ, true, sectionRot[2]);
        KM_Eo_Ed.addRotation(localX, true, sectionRot[0]);
        KM_Eo_Ed.addTranslation(sectionDisp);

        /*cout << "sectionDisp" << endl;
         cout << sectionDisp[0] <<" "<< sectionDisp[1] <<" "<< sectionDisp[2] << endl;
         cout << "sectionRot" << endl;
         cout << sectionRot[0] <<" "<< sectionRot[1] <<" "<< sectionRot[2] << endl;*/

        KinematicMotion TRL_O_P;
        TRL_O_P.addTranslation(&sectionP[i * 3]);
        KinematicMotion *TRL_P_O = TRL_O_P.newInverse();

        KinematicMotion KM_Oo_Od;
        KinematicMotion *ROT_ELEM_O = ROT_O_ELEM[elem]->newInverse();
        KM_Oo_Od.addKinematicMotion(TRL_P_O);
        KM_Oo_Od.addKinematicMotion(ROT_ELEM_O);
        KM_Oo_Od.addKinematicMotion(&KM_Eo_Ed);
        KM_Oo_Od.addKinematicMotion(ROT_O_ELEM[elem]);
        KM_Oo_Od.addKinematicMotion(&TRL_O_P);

        delete ROT_ELEM_O;
        delete TRL_P_O;

        // move the section nodes, get the displacements
        int numSectionNodes;
        if (i == 0) { // root
            numSectionNodes = surfaceNumRootSectionNodes;
        } else if (i == surfaceNumSections - 1) { // tip
            numSectionNodes = surfaceNumTipSectionNodes;
        } else { // normal
            numSectionNodes = surfaceNumNormalSectionNodes;
        }
        for (int j = 0; j < numSectionNodes; j++) {
            int pos;
            if (i == 0) { // root
                pos = sortedPosToUnsortedPos[j];
            } else if (i == surfaceNumSections - 1) { // tip
                pos = sortedPosToUnsortedPos[surfaceNumNodes - j - 1];
            } else { // normal
                pos = sortedPosToUnsortedPos[surfaceNumRootSectionNodes
                        + (i - 1) * surfaceNumNormalSectionNodes + j];
            }
            double nodalOldCoor[3];
            double nodalNewCoor[3];
            for (int k = 0; k < 3; k++) {
                nodalOldCoor[k] = surfaceNodeCoors[pos * 3 + k];
                nodalNewCoor[k] = surfaceNodeCoors[pos * 3 + k];
            }
            KM_Oo_Od.move(nodalNewCoor);
            for (int k = 0; k < 3; k++)
                surfaceDisp[pos * 3 + k] = nodalNewCoor[k] - nodalOldCoor[k];
        }
    }
}

void CurveSurfaceLinearMapper::conservativeMapping(const double *surfaceForce, double *curveForceMoment) {
    // set forces and moments on the curve/beam to 0
    for (int i = 0; i < curveNumNodes * 6; i++) {
        curveForceMoment[i] = 0.0;
    }

    for (int i = 0; i < surfaceNumSections; i++) {
        int elem = sectionToCurveElem[i];
        int node1ID = curveElems[elem * 2 + 0];
        int node2ID = curveElems[elem * 2 + 1];
        int node1Pos = curveNodeIDToPos->at(node1ID);
        int node2Pos = curveNodeIDToPos->at(node2ID);
        double linearShapeFunc1 = shapeFuncOfSection[i * 10 + 0];
        double linearShapeFunc2 = shapeFuncOfSection[i * 10 + 5 + 0];
        int numSectionNodes;
        if (i == 0) { // root
            numSectionNodes = surfaceNumRootSectionNodes;
        } else if (i == surfaceNumSections - 1) { // tip
            numSectionNodes = surfaceNumTipSectionNodes;
        } else { // normal
            numSectionNodes = surfaceNumNormalSectionNodes;
        }
        for (int j = 0; j < numSectionNodes; j++) {
            int pos;
            if (i == 0) { // root
                pos = sortedPosToUnsortedPos[j];
            } else if (i == surfaceNumSections - 1) { // tip
                pos = sortedPosToUnsortedPos[surfaceNumNodes - j - 1];
            } else { // normal
                pos = sortedPosToUnsortedPos[surfaceNumRootSectionNodes
                        + (i - 1) * surfaceNumNormalSectionNodes + j];
            }
            // distance from P to a section node
            double R[3];
            for (int k = 0; k < 3; k++) {
                R[k] = surfaceNodeCoors[pos * 3 + k] - sectionP[i * 3 + k];
            }
            // force vector on the nodal
            double F[3];
            for (int k = 0; k < 3; k++) {
                F[k] = surfaceForce[pos * 3 + k];
            }
            // compute moment around P with the nodal force: M = R X F
            double M[3];
            M[0] = R[1] * F[2] - R[2] * F[1];
            M[1] = R[2] * F[0] - R[0] * F[2];
            M[2] = R[0] * F[1] - R[1] * F[0];

            /*cout << "F" << endl;
            cout << F[0] << " " << F[1] << " " << F[2] << endl;
            cout << "M" << endl;
            cout << M[0] << " " << M[1] << " " << M[2] << endl;
            cout << "linearShapeFunc1: " << linearShapeFunc1 << endl;
            cout << "linearShapeFunc1: " << linearShapeFunc1 << endl;*/
            // split F and M to curve/beam end nodes
            for (int k = 0; k < 3; k++) {
                curveForceMoment[node1Pos * 6 + 0 + k] += linearShapeFunc1 * F[k];
                curveForceMoment[node1Pos * 6 + 3 + k] += linearShapeFunc1 * M[k];
                curveForceMoment[node2Pos * 6 + 0 + k] += linearShapeFunc2 * F[k];
                curveForceMoment[node2Pos * 6 + 3 + k] += linearShapeFunc2 * M[k];
            }
        }
    }
}

} /* namespace EMPIRE */
