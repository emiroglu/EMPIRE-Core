#include "CurveSurfaceMapper.h"
#include "KinematicMotion.h"
#include <assert.h>
#include <map>
#include <iostream>
#include <math.h>

using namespace std;

namespace EMPIRE {

CurveSurfaceMapper::CurveSurfaceMapper(EMPIRE_CurveSurfaceMapper_type _type, int _curveNumNodes,
        int _curveNumElements, const double *_curveNodeCoors, const int *_curveNodeIDs,
        const int *_curveElems, int _surfaceNumNodes, const double *_surfaceNodeCoors,
        int _surfaceNumSections, int _surfaceNumRootSectionNodes, int _surfaceNumNormalSectionNodes,
        int _surfaceNumTipSectionNodes, const double *rotation_O_Q, const double *translation_O_Q) :
        type(_type), curveNumNodes(_curveNumNodes), curveNumElements(_curveNumElements), curveNodeCoors(
                _curveNodeCoors), curveElems(_curveElems), surfaceNumNodes(_surfaceNumNodes), surfaceNodeCoors(
                _surfaceNodeCoors), surfaceNumSections(_surfaceNumSections), surfaceNumRootSectionNodes(
                _surfaceNumRootSectionNodes), surfaceNumNormalSectionNodes(
                _surfaceNumNormalSectionNodes), surfaceNumTipSectionNodes(
                _surfaceNumTipSectionNodes) {
    /*
     * Coordinate systems:
     *   O --- global
     *   Q --- beam root, origin is an arbitrary point in the root section
     *   P --- origin is the cross point of the section with the beam element, orientation is the same as the global system
     */
    // construct KM_O_Q
    KM_O_Q = new KinematicMotion();
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

    // check number of surface nodes
    assert(
            surfaceNumNodes
                    == surfaceNumRootSectionNodes + surfaceNumTipSectionNodes
                            + (surfaceNumSections - 2) * surfaceNumNormalSectionNodes);

    // sort surface nodes according to x in Q system
    sortedPosToUnsortedPos = new int[surfaceNumNodes];

    cout << "Andreas started debugging" << endl;
    
    multimap<double, int> *mapCoorX2Pos = new multimap<double, int>; // sort by x coordinate
    for (int i = 0; i < surfaceNumNodes; i++) {
        mapCoorX2Pos->insert(pair<double, int>(surfaceNodeCoorsInQ[i * 3 + 0], i));
    }
    
    cout << "Andreas ended debugging" << endl;

    int count = 0;
    for (multimap<double, int>::iterator it = mapCoorX2Pos->begin(); it != mapCoorX2Pos->end();
            it++) {
        sortedPosToUnsortedPos[count] = it->second;
        count++;
    }
    
    delete mapCoorX2Pos;
    
    // curve node coordinates in Q
    curveNodeCoorsInQ = new double[curveNumNodes * 3];
    for (int i = 0; i < curveNumNodes * 3; i++)
        curveNodeCoorsInQ[i] = curveNodeCoors[i];
    for (int i = 0; i < curveNumNodes; i++) {
        KM_Q_O->move(&curveNodeCoorsInQ[i * 3]);
    }

    // map curve node ID to curve node position
    curveNodeIDToPos = new map<int, int>;
    for (int i = 0; i < curveNumNodes; i++) {
        curveNodeIDToPos->insert(pair<int, int>(_curveNodeIDs[i], i));
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
        //double cubicShapeFuncRot1 = 1.0 / 8.0 * length * (1.0 - xi) * (1.0 - xi) * (1.0 + xi);
        double cubicShapeFuncRot1 = 1.0 / 8.0 * (1.0 - xi) * (1.0 - xi) * (1.0 + xi); // *=length
        double cubicShapeFuncDisp2 = 1.0 / 4.0 * (1.0 + xi) * (1.0 + xi) * (2.0 - xi);
        //double cubicShapeFuncRot2 = -1.0 / 8.0 * length * (1.0 + xi) * (1.0 + xi) * (1.0 - xi);
        double cubicShapeFuncRot2 = -1.0 / 8.0 * (1.0 + xi) * (1.0 + xi) * (1.0 - xi); // *=length
        //double Dxi_Dx = 2.0 / length;
        //double cubicShapeFuncDispDeriv1 = Dxi_Dx * (-3.0) / 4.0 * (1.0 - xi) * (1.0 + xi);
        double cubicShapeFuncDispDeriv1 = 2.0 * (-3.0) / 4.0 * (1.0 - xi) * (1.0 + xi); // /=length
        double cubicShapeFuncRotDeriv1 = 2.0 * (-1.0) / 8.0 * (1.0 - xi) * (1.0 + 3 * xi);
        //double cubicShapeFuncDispDeriv2 = Dxi_Dx * 3.0 / 4.0 * (1.0 - xi) * (1.0 + xi);
        double cubicShapeFuncDispDeriv2 = 2.0 * 3.0 / 4.0 * (1.0 - xi) * (1.0 + xi); // /=length
        double cubicShapeFuncRotDeriv2 = 2.0 * (-1.0) / 8.0 * (1.0 + xi) * (1.0 - 3 * xi);

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
                    * curveNodeCoors[curveNodeIDToPos->at(node1ID) * 3 + j]
                    + shapeFuncOfSection[i * 10 + 5 + 0]
                            * curveNodeCoors[curveNodeIDToPos->at(node2ID) * 3 + j];
        }
    }

    // transformation from local beam to the global system
    ROT_O_ELEM = new KinematicMotion*[curveNumElements];
    // compute ROT_O_ELEM
    for (int i = 0; i < curveNumElements; i++) {
        // For the definition of local axes, see carat ElementBeam1::calc_transformation_matrix, or
        // carat.st.bv.tum.de/caratuserswiki/index.php/Users:General_FEM_Analysis/Elements_Reference/Beam1
        ROT_O_ELEM[i] = new KinematicMotion;
        int node1ID = curveElems[i * 2 + 0];
        int node2ID = curveElems[i * 2 + 1];
        const double *node1 = &curveNodeCoors[curveNodeIDToPos->at(node1ID) * 3];
        const double *node2 = &curveNodeCoors[curveNodeIDToPos->at(node2ID) * 3];
        double tmp[3];
        for (int j = 0; j < 3; j++)
            tmp[j] = node2[j] - node1[j];
        double length = sqrt(tmp[0] * tmp[0] + tmp[1] * tmp[1] + tmp[2] * tmp[2]); // do not use normalizeVector here
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

    // initialize section rotation for conservative mapping
    sectionRot = new double[surfaceNumSections * 3];

    // delete
    delete[] surfaceNodeCoorsInQ;
    delete rightXToElemPos;
    if (type != EMPIRE_CurveSurfaceMapper_corotate2D) {
        delete KM_Q_O;
        delete KM_O_Q;
        delete[] curveNodeCoorsInQ;
    }

    // special initialization for different algorithms
    if (type == EMPIRE_CurveSurfaceMapper_linear) {
        curveElemLength = new double[curveNumElements];
        for (int i = 0; i < curveNumElements; i++) {
            int node1ID = curveElems[i * 2 + 0];
            int node2ID = curveElems[i * 2 + 1];
            const double *node1 = &curveNodeCoors[curveNodeIDToPos->at(node1ID) * 3];
            const double *node2 = &curveNodeCoors[curveNodeIDToPos->at(node2ID) * 3];

            // corotate axis in element system
            double elemVec[3];
            for (int j = 0; j < 3; j++)
                elemVec[j] = node2[j] - node1[j];

            curveElemLength[i] = normalizeVector(elemVec);
            assert(curveElemLength[i] != 0.0);
        }
    } else if (type == EMPIRE_CurveSurfaceMapper_corotate2D) {
        // construct KM_O_Q
        const double EPS = 1E-10;
        for (int i = 0; i < 9; i++) {
            double tmp = fabs(rotation_O_Q[i]);
            if (tmp > 0.5)
                assert(fabs(tmp - 1.0) < EPS); // should be 1
            else
                assert(fabs(tmp - 0.0) < EPS); // should be 0
        }

        // z coordinates should be 0 in Q
        for (int i = 0; i < curveNumNodes; i++) {
            assert(fabs(curveNodeCoorsInQ[i * 3 + 2] - 0.0) < EPS);
        }

        // transformation from local beam to Q
        ROT_Q_ELEM = new KinematicMotion*[curveNumElements];
        // rotation angle from local beam to Q
        angle_Q_ELEM = new double[curveNumElements];
        // compute ROT_Q_ELEM and angle_Q_ELEM
        for (int i = 0; i < curveNumElements; i++) {
            // For the definition of local axes, see carat ElementBeam1::calc_transformation_matrix, or
            // carat.st.bv.tum.de/caratuserswiki/index.php/Users:General_FEM_Analysis/Elements_Reference/Beam1
            // x and y are defined on Q instead of O, which is different than the classic carat definition
            ROT_Q_ELEM[i] = new KinematicMotion;
            int node1ID = curveElems[i * 2 + 0];
            int node2ID = curveElems[i * 2 + 1];
            const double *node1 = &curveNodeCoorsInQ[curveNodeIDToPos->at(node1ID) * 3];
            const double *node2 = &curveNodeCoorsInQ[curveNodeIDToPos->at(node2ID) * 3];
            double tmp[2];
            for (int j = 0; j < 2; j++) // use only Qx-Qy
                tmp[j] = node2[j] - node1[j];

            double lengthXY = sqrt(tmp[0] * tmp[0] + tmp[1] * tmp[1]);
            {
                double x_loc[3];
                double y_loc[3];
                double z_loc[3];
                // mapping to global coordinates through direction cosines
                x_loc[0] = tmp[0] / lengthXY; //direction cosine x
                x_loc[1] = tmp[1] / lengthXY; //direction cosine y
                x_loc[2] = 0.0; //direction cosine z

                // using the projection L_x of x_loc to x_gl-y_gl-plane as perpendicular to y_loc
                y_loc[0] = -x_loc[1];
                y_loc[1] = x_loc[0];
                y_loc[2] = 0.0;

                // z-direction via cross product of x and y vectors
                z_loc[0] = 0.0;
                z_loc[1] = 0.0;
                z_loc[2] = 1.0;

                ROT_Q_ELEM[i]->addRotation(x_loc, y_loc, z_loc, true);
                angle_Q_ELEM[i] = acos(x_loc[0]);
                if (x_loc[1] < 0.0) {
                    angle_Q_ELEM[i] = 2.0 * M_PI - angle_Q_ELEM[i];
                }
            }
        }
        delete KM_Q_O;
    } else if (type == EMPIRE_CurveSurfaceMapper_corotate3D) {
        // do nothing
    } else {
        assert(false);
    }
    
}

CurveSurfaceMapper::~CurveSurfaceMapper() {
    delete[] sortedPosToUnsortedPos;
    delete curveNodeIDToPos;
    delete[] sectionToCurveElem;
    delete[] shapeFuncOfSection;
    for (int i = 0; i < curveNumElements; i++) {
        delete ROT_O_ELEM[i];
    }
    delete[] ROT_O_ELEM;
    delete[] sectionP;
    delete[] sectionRot;

    // special destruction for different algorithms
    if (type == EMPIRE_CurveSurfaceMapper_linear) {
        delete[] curveElemLength;
    } else if (type == EMPIRE_CurveSurfaceMapper_corotate2D) {

        for (int i = 0; i < curveNumElements; i++) {
            delete ROT_Q_ELEM[i];
        }
        delete[] ROT_Q_ELEM;
        delete KM_O_Q;
        delete[] angle_Q_ELEM;
        delete[] curveNodeCoorsInQ;
    } else if (type == EMPIRE_CurveSurfaceMapper_corotate3D) {
        // do nothing
    } else {
        assert(false);
    }
}

void CurveSurfaceMapper::consistentMapping(const double *curveDispRot, double *surfaceDisp) {
    if (type == EMPIRE_CurveSurfaceMapper_linear) {
        consistentMappingLinear(curveDispRot, surfaceDisp);
    } else if (type == EMPIRE_CurveSurfaceMapper_corotate2D) {
        consistentMappingCorotate2D(curveDispRot, surfaceDisp);
    } else if (type == EMPIRE_CurveSurfaceMapper_corotate3D) {
        consistentMappingCorotate3D(curveDispRot, surfaceDisp);
    } else {
        assert(false);
    }
}
void CurveSurfaceMapper::consistentMappingLinear(const double *curveDispRot, double *surfaceDisp) {
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
        double tmpSectionDisp[3];
        double tmpSectionRot[3];

        double length = curveElemLength[elem];
        double linearShapeFunc1 = shapeFuncOfSection[i * 10 + 0];
        double cubicShapeFuncDisp1 = shapeFuncOfSection[i * 10 + 1];
        double cubicShapeFuncRot1 = shapeFuncOfSection[i * 10 + 2] * length;
        double cubicShapeFuncDispDeriv1 = shapeFuncOfSection[i * 10 + 3] / length;
        double cubicShapeFuncRotDeriv1 = shapeFuncOfSection[i * 10 + 4];
        double linearShapeFunc2 = shapeFuncOfSection[i * 10 + 5 + 0];
        double cubicShapeFuncDisp2 = shapeFuncOfSection[i * 10 + 5 + 1];
        double cubicShapeFuncRot2 = shapeFuncOfSection[i * 10 + 5 + 2] * length;
        double cubicShapeFuncDispDeriv2 = shapeFuncOfSection[i * 10 + 5 + 3] / length;
        double cubicShapeFuncRotDeriv2 = shapeFuncOfSection[i * 10 + 5 + 4];
        //disp_x = NL_1*disp_x1 + NL_2*disp_x2
        tmpSectionDisp[0] = linearShapeFunc1 * node1DispLocal[elem * 3 + 0]
                + linearShapeFunc2 * node2DispLocal[elem * 3 + 0];
        //disp_y = NCd_1*disp_y1 + NCd_2*disp_y2 + NCr_1*rot_z1 + NCr_2*rot_z2
        tmpSectionDisp[1] = cubicShapeFuncDisp1 * node1DispLocal[elem * 3 + 1]
                + cubicShapeFuncDisp2 * node2DispLocal[elem * 3 + 1]
                + cubicShapeFuncRot1 * node1RotLocal[elem * 3 + 2]
                + cubicShapeFuncRot2 * node2RotLocal[elem * 3 + 2];
        //disp_z = NCd_1*disp_z1 + NCd_2*disp_z2 + NCr_1*(-rot_y1) + NCr_2*(-rot_y2), since disp_z' = -rot_y
        tmpSectionDisp[2] = cubicShapeFuncDisp1 * node1DispLocal[elem * 3 + 2]
                + cubicShapeFuncDisp2 * node2DispLocal[elem * 3 + 2]
                + cubicShapeFuncRot1 * (-node1RotLocal[elem * 3 + 1])
                + cubicShapeFuncRot2 * (-node2RotLocal[elem * 3 + 1]);
        //rot_x = NL_1*rot_x1 + NL_2*rot_x2
        tmpSectionRot[0] = linearShapeFunc1 * node1RotLocal[elem * 3 + 0]
                + linearShapeFunc2 * node2RotLocal[elem * 3 + 0];
        //rot_y = - disp_z' = - [NCd'_1*disp_z1 + NCd'_2*disp_z2 + NCr'_1*(-rot_y1) + NCr'_2*(-rot_y2)]
        tmpSectionRot[1] = -(cubicShapeFuncDispDeriv1 * node1DispLocal[elem * 3 + 2]
                + cubicShapeFuncDispDeriv2 * node2DispLocal[elem * 3 + 2]
                + cubicShapeFuncRotDeriv1 * (-node1RotLocal[elem * 3 + 1])
                + cubicShapeFuncRotDeriv2 * (-node2RotLocal[elem * 3 + 1]));
        //rot_z = disp_y' = NCd'_1*disp_y1 + NCd'_2*disp_y2 + NCr'_1*rot_z1 + NCr'_2*rot_z2
        tmpSectionRot[2] = cubicShapeFuncDispDeriv1 * node1DispLocal[elem * 3 + 1]
                + cubicShapeFuncDispDeriv2 * node2DispLocal[elem * 3 + 1]
                + cubicShapeFuncRotDeriv1 * node1RotLocal[elem * 3 + 2]
                + cubicShapeFuncRotDeriv2 * node2RotLocal[elem * 3 + 2];

        // Compute KM_Oo_Od
        KinematicMotion KM_Eo_Ed;
        double localX[] = { 1.0, 0.0, 0.0 };
        double localY[] = { 0.0, 1.0, 0.0 };
        double localZ[] = { 0.0, 0.0, 1.0 };
        // bending first, then torsion
        KM_Eo_Ed.addRotation(localY, true, tmpSectionRot[1]);
        KM_Eo_Ed.addRotation(localZ, true, tmpSectionRot[2]);
        KM_Eo_Ed.addRotation(localX, true, tmpSectionRot[0]);
        KM_Eo_Ed.addTranslation(tmpSectionDisp);

        /*cout << "tmpSectionDisp" << endl;
         cout << tmpSectionDisp[0] <<" "<< tmpSectionDisp[1] <<" "<< tmpSectionDisp[2] << endl;
         cout << "tmpSectionRot" << endl;
         cout << tmpSectionRot[0] <<" "<< tmpSectionRot[1] <<" "<< tmpSectionRot[2] << endl;*/

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
        // compute the final section displacement and rotation from KM_Oo_Ord
        KM_Oo_Od.getRotationVector(&sectionRot[i * 3]);
    }
    delete[] node1DispLocal;
    delete[] node1RotLocal;
    delete[] node2DispLocal;
    delete[] node2RotLocal;
}
void CurveSurfaceMapper::consistentMappingCorotate2D(const double *curveDispRot,
        double *surfaceDisp) {
    // Rotation ROT_O_Q
    KinematicMotion ROT_O_Q;
    ROT_O_Q.addRotation(KM_O_Q->getRotationMatrix());
    KinematicMotion *ROT_Q_O = ROT_O_Q.newInverse();
    // DOFs on curve/beam element local system
    double *node1Disp = new double[curveNumElements * 3];
    double *node1Rot = new double[curveNumElements * 3];
    double *node2Disp = new double[curveNumElements * 3];
    double *node2Rot = new double[curveNumElements * 3];
    // Corotate angles
    double *angle_elem_corotate = new double[curveNumElements];
    // Corotate transformations
    KinematicMotion **ROT_elem_corotate = new KinematicMotion*[curveNumElements];
    // length of the element after deformation
    double *corotateElementLength = new double[curveNumElements];
    // compute DOFs in the curve element local system, compute corotate transformations, current length of the element
    for (int i = 0; i < curveNumElements; i++) {
        int node1ID = curveElems[i * 2 + 0];
        int node2ID = curveElems[i * 2 + 1];
        for (int j = 0; j < 3; j++) {
            node1Disp[i * 3 + j] = curveDispRot[curveNodeIDToPos->at(node1ID) * 6 + 0 + j];
            node1Rot[i * 3 + j] = curveDispRot[curveNodeIDToPos->at(node1ID) * 6 + 3 + j];
            node2Disp[i * 3 + j] = curveDispRot[curveNodeIDToPos->at(node2ID) * 6 + 0 + j];
            node2Rot[i * 3 + j] = curveDispRot[curveNodeIDToPos->at(node2ID) * 6 + 3 + j];
        }

        // transform the DOFs to Q system
        ROT_Q_O->move(&node1Disp[i * 3]);
        ROT_Q_O->move(&node1Rot[i * 3]);
        ROT_Q_O->move(&node2Disp[i * 3]);
        ROT_Q_O->move(&node2Rot[i * 3]);
        { // compute corotate angle in Q
            const double *node1 = &curveNodeCoorsInQ[curveNodeIDToPos->at(node1ID) * 3];
            const double *node2 = &curveNodeCoorsInQ[curveNodeIDToPos->at(node2ID) * 3];
            double tmp[2];
            for (int j = 0; j < 2; j++) // use only Qx-Qy
                tmp[j] = node2[j] + node2Disp[i * 3 + j] - node1[j] - node1Disp[i * 3 + j];

            double lengthXY = sqrt(tmp[0] * tmp[0] + tmp[1] * tmp[1]);
            tmp[0] /= lengthXY; //direction cosine x
            tmp[1] /= lengthXY; //direction cosine y

            double elem_Q_corotate = acos(tmp[0]);
            if (tmp[1] < 0.0) {
                elem_Q_corotate = 2.0 * M_PI - elem_Q_corotate;
            }
            angle_elem_corotate[i] = elem_Q_corotate - angle_Q_ELEM[i];

            double z_loc[] = { 0.0, 0.0, 1.0 };
            ROT_elem_corotate[i] = new KinematicMotion;
            ROT_elem_corotate[i]->addRotation(z_loc, true, angle_elem_corotate[i]);

            corotateElementLength[i] = lengthXY;
        }

        // transform the DOFs to element local system
        KinematicMotion *ROT_ELEM_Q = ROT_Q_ELEM[i]->newInverse();
        ROT_ELEM_Q->move(&node1Disp[i * 3]);
        ROT_ELEM_Q->move(&node1Rot[i * 3]);
        ROT_ELEM_Q->move(&node2Disp[i * 3]);
        ROT_ELEM_Q->move(&node2Rot[i * 3]);

        delete ROT_ELEM_Q;
        //cout << "node1Rot" << endl;
        //cout << node1Rot[i * 3 + 0] <<" "<< node1Rot[i * 3 + 1] <<" "<< node1Rot[i * 3 + 2] << endl;
        //cout << "node2Rot" << endl;
        //cout << node2Rot[i * 3 + 0] <<" "<< node2Rot[i * 3 + 1] <<" "<< node2Rot[i * 3 + 2] << endl;
    }

    // compute displacements of nodes according to rotations and displacements of sections
    for (int i = 0; i < surfaceNumSections; i++) {
        // interpolate section rot and disp on the element local system
        int elem = sectionToCurveElem[i];
        double tmpSectionDisp[3];
        double tmpSectionRot[3];

        double length = corotateElementLength[elem];
        double linearShapeFunc1 = shapeFuncOfSection[i * 10 + 0];
        double cubicShapeFuncDisp1 = shapeFuncOfSection[i * 10 + 1];
        double cubicShapeFuncRot1 = shapeFuncOfSection[i * 10 + 2] * length;
        double cubicShapeFuncDispDeriv1 = shapeFuncOfSection[i * 10 + 3] / length;
        double cubicShapeFuncRotDeriv1 = shapeFuncOfSection[i * 10 + 4];
        double linearShapeFunc2 = shapeFuncOfSection[i * 10 + 5 + 0];
        double cubicShapeFuncDisp2 = shapeFuncOfSection[i * 10 + 5 + 1];
        double cubicShapeFuncRot2 = shapeFuncOfSection[i * 10 + 5 + 2] * length;
        double cubicShapeFuncDispDeriv2 = shapeFuncOfSection[i * 10 + 5 + 3] / length;
        double cubicShapeFuncRotDeriv2 = shapeFuncOfSection[i * 10 + 5 + 4];

        // compute disp_y and rot_z in corotate system with origin P
        tmpSectionDisp[0] = 0.0;
        tmpSectionDisp[1] = cubicShapeFuncRot1
                * (node1Rot[elem * 3 + 2] - angle_elem_corotate[elem])
                + cubicShapeFuncRot2 * (node2Rot[elem * 3 + 2] - angle_elem_corotate[elem]);
        tmpSectionDisp[2] = 0.0;
        tmpSectionRot[0] = 0.0;
        tmpSectionRot[1] = 0.0;
        tmpSectionRot[2] = cubicShapeFuncRotDeriv1
                * (node1Rot[elem * 3 + 2] - angle_elem_corotate[elem])
                + cubicShapeFuncRotDeriv2 * (node2Rot[elem * 3 + 2] - angle_elem_corotate[elem]);

        // transform the section disp to the element local system (rot does not change)
        ROT_elem_corotate[elem]->move(tmpSectionDisp);

        // compute the DOFs on the element local system with origin P
        //disp_x += NL_1*disp_x1 + NL_2*disp_x2
        tmpSectionDisp[0] += linearShapeFunc1 * node1Disp[elem * 3 + 0]
                + linearShapeFunc2 * node2Disp[elem * 3 + 0];
        //disp_y += NL_1*disp_y1 + NL_2*disp_y2
        tmpSectionDisp[1] += linearShapeFunc1 * node1Disp[elem * 3 + 1]
                + linearShapeFunc2 * node2Disp[elem * 3 + 1];
        //disp_z += NCd_1*disp_z1 + NCd_2*disp_z2 + NCr_1*(-rot_y1) + NCr_2*(-rot_y2), since disp_z' = -rot_y
        tmpSectionDisp[2] += cubicShapeFuncDisp1 * node1Disp[elem * 3 + 2]
                + cubicShapeFuncDisp2 * node2Disp[elem * 3 + 2]
                + cubicShapeFuncRot1 * (-node1Rot[elem * 3 + 1])
                + cubicShapeFuncRot2 * (-node2Rot[elem * 3 + 1]);
        //rot_x += NL_1*rot_x1 + NL_2*rot_x2
        tmpSectionRot[0] += linearShapeFunc1 * node1Rot[elem * 3 + 0]
                + linearShapeFunc2 * node2Rot[elem * 3 + 0];
        //rot_y += - disp_z' = - [NCd'_1*disp_z1 + NCd'_2*disp_z2 + NCr'_1*(-rot_y1) + NCr'_2*(-rot_y2)]
        tmpSectionRot[1] += -(cubicShapeFuncDispDeriv1 * node1Disp[elem * 3 + 2]
                + cubicShapeFuncDispDeriv2 * node2Disp[elem * 3 + 2]
                + cubicShapeFuncRotDeriv1 * (-node1Rot[elem * 3 + 1])
                + cubicShapeFuncRotDeriv2 * (-node2Rot[elem * 3 + 1]));
        //rot_z += corotate angle
        tmpSectionRot[2] += angle_elem_corotate[elem];

        // Compute KM_Oo_Ord
        KinematicMotion KM_Eo_Erd;
        double localX[] = { 1.0, 0.0, 0.0 };
        double localY[] = { 0.0, 1.0, 0.0 };
        double localZ[] = { 0.0, 0.0, 1.0 };
        // bending first, then torsion, the corotate (the order of first two does not matter, the corotate must be the last)
        KM_Eo_Erd.addRotation(localY, true, tmpSectionRot[1]);
        KM_Eo_Erd.addRotation(localX, true, tmpSectionRot[0]);
        KM_Eo_Erd.addRotation(localZ, true, tmpSectionRot[2]);
        KM_Eo_Erd.addTranslation(tmpSectionDisp);

        /*cout << "tmpSectionDisp" << endl;
         cout << tmpSectionDisp[0] <<" "<< tmpSectionDisp[1] <<" "<< tmpSectionDisp[2] << endl;
         cout << "tmpSectionRot" << endl;
         cout << tmpSectionRot[0] <<" "<< tmpSectionRot[1] <<" "<< tmpSectionRot[2] << endl;*/

        KinematicMotion TRL_O_P;
        TRL_O_P.addTranslation(&sectionP[i * 3]);
        KinematicMotion *TRL_P_O = TRL_O_P.newInverse();

        KinematicMotion KM_Oo_Ord;

        KinematicMotion ROT_O_ELEM;
        ROT_O_ELEM.addKinematicMotion(ROT_Q_ELEM[elem]);
        ROT_O_ELEM.addKinematicMotion(&ROT_O_Q);
        KinematicMotion *ROT_ELEM_O = ROT_O_ELEM.newInverse();

        KM_Oo_Ord.addKinematicMotion(TRL_P_O);
        KM_Oo_Ord.addKinematicMotion(ROT_ELEM_O);
        KM_Oo_Ord.addKinematicMotion(&KM_Eo_Erd);
        KM_Oo_Ord.addKinematicMotion(&ROT_O_ELEM);
        KM_Oo_Ord.addKinematicMotion(&TRL_O_P);

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
            KM_Oo_Ord.move(nodalNewCoor);
            for (int k = 0; k < 3; k++)
                surfaceDisp[pos * 3 + k] = nodalNewCoor[k] - nodalOldCoor[k];
        }

        // compute the final section displacement and rotation from KM_Oo_Ord
        KM_Oo_Ord.getRotationVector(&sectionRot[i * 3]);
    }
    delete[] node1Disp;
    delete[] node1Rot;
    delete[] node2Disp;
    delete[] node2Rot;

    delete[] angle_elem_corotate;
    for (int i = 0; i < curveNumElements; i++) {
        delete ROT_elem_corotate[i];
    }
    delete[] ROT_elem_corotate;
    delete ROT_Q_O;
    delete[] corotateElementLength;
}
void CurveSurfaceMapper::consistentMappingCorotate3D(const double *curveDispRot,
        double *surfaceDisp) {
    // motion definition: o---original, r---corotating, d---deformed, rd---r+d
    // DOFs on curve/beam element local system
    double *node1Disp = new double[curveNumElements * 3];
    double *node1Rot = new double[curveNumElements * 3];
    double *node2Disp = new double[curveNumElements * 3];
    double *node2Rot = new double[curveNumElements * 3];
    // Corotate transformations
    KinematicMotion **ROT_Eo_Er = new KinematicMotion*[curveNumElements];
    KinematicMotion **ROT_torsion = new KinematicMotion*[curveNumElements];
    // length of the element after deformation
    double *corotateElementLength = new double[curveNumElements];
    // compute DOFs in the curve element local system, compute corotate transformations, current length of the element
    for (int i = 0; i < curveNumElements; i++) {
        int node1ID = curveElems[i * 2 + 0];
        int node2ID = curveElems[i * 2 + 1];
        for (int j = 0; j < 3; j++) {
            node1Disp[i * 3 + j] = curveDispRot[curveNodeIDToPos->at(node1ID) * 6 + 0 + j];
            node1Rot[i * 3 + j] = curveDispRot[curveNodeIDToPos->at(node1ID) * 6 + 3 + j];
            node2Disp[i * 3 + j] = curveDispRot[curveNodeIDToPos->at(node2ID) * 6 + 0 + j];
            node2Rot[i * 3 + j] = curveDispRot[curveNodeIDToPos->at(node2ID) * 6 + 3 + j];
        }
        KinematicMotion *ROT_ELEM_O = ROT_O_ELEM[i]->newInverse();
        { // compute ROT_elem_corotate, see Non-linear Modeling and Analysis of Solids and Structures (Krenk2009) P136
            const double *node1 = &curveNodeCoors[curveNodeIDToPos->at(node1ID) * 3];
            const double *node2 = &curveNodeCoors[curveNodeIDToPos->at(node2ID) * 3];

            // corotate axis in element system
            double corotateXAxis[3];
            for (int j = 0; j < 3; j++)
                corotateXAxis[j] = node2[j] + node2Disp[i * 3 + j] - node1[j]
                        - node1Disp[i * 3 + j]; // the disp. now is global instead of local

            corotateElementLength[i] = normalizeVector(corotateXAxis);
            assert(corotateElementLength[i] != 0.0);

            ROT_ELEM_O->move(corotateXAxis);

            // average axis
            double elemXAxis[] = { 1.0, 0.0, 0.0 };
            double averageAxis[3];
            for (int j = 0; j < 3; j++)
                averageAxis[j] = corotateXAxis[j] + elemXAxis[j];
            double averageAxisLength = normalizeVector(averageAxis);
            assert(averageAxisLength != 0.0);

            // new corotateAxis in the element local system, only the computation result is coded here
            double corotateYAxis[3];
            double corotateZAxis[3];
            corotateYAxis[0] = -2.0 * averageAxis[0] * averageAxis[1];
            corotateYAxis[1] = 1.0 - 2.0 * averageAxis[1] * averageAxis[1]; // diagonal
            corotateYAxis[2] = -2.0 * averageAxis[2] * averageAxis[1];
            corotateZAxis[0] = -2.0 * averageAxis[0] * averageAxis[2];
            corotateZAxis[1] = -2.0 * averageAxis[1] * averageAxis[2];
            corotateZAxis[2] = 1.0 - 2.0 * averageAxis[2] * averageAxis[2]; // diagonal

            // compute ROT_elem_corotate
            ROT_Eo_Er[i] = new KinematicMotion;
            ROT_Eo_Er[i]->addRotation(corotateXAxis, corotateYAxis, corotateZAxis, true);
        }

        { // transform the DOFs to element local system
            KinematicMotion *ROT_Er_Eo = ROT_Eo_Er[i]->newInverse();

            // step 1: transform rotation to element local system
            double angleNode1 = normalizeVector(&node1Rot[i * 3]);
            KinematicMotion ROT_Eo_Ed__node1;
            if (angleNode1 != 0.0) {
                KinematicMotion ROT_Oo_Ord__node1;
                KinematicMotion ROT_Eo_Erd__node1;
                ROT_Oo_Ord__node1.addRotation(&node1Rot[i * 3], true, angleNode1);
                ROT_Eo_Erd__node1.addRotation(ROT_O_ELEM[i]->getRotationMatrix());
                ROT_Eo_Erd__node1.addRotation(ROT_Oo_Ord__node1.getRotationMatrix());
                ROT_Eo_Erd__node1.addRotation(ROT_ELEM_O->getRotationMatrix()); // O->E
                ROT_Eo_Ed__node1.addRotation(ROT_Eo_Erd__node1.getRotationMatrix());
            }

            double angleNode2 = normalizeVector(&node2Rot[i * 3]);
            KinematicMotion ROT_Eo_Ed__node2;
            if (angleNode2 != 0.0) {
                KinematicMotion ROT_Oo_Ord__node2;
                KinematicMotion ROT_Eo_Erd__node2;
                ROT_Oo_Ord__node2.addRotation(&node2Rot[i * 3], true, angleNode2);
                ROT_Eo_Erd__node2.addRotation(ROT_O_ELEM[i]->getRotationMatrix());
                ROT_Eo_Erd__node2.addRotation(ROT_Oo_Ord__node2.getRotationMatrix());
                ROT_Eo_Erd__node2.addRotation(ROT_ELEM_O->getRotationMatrix()); // O->E
                ROT_Eo_Ed__node2.addRotation(ROT_Eo_Erd__node2.getRotationMatrix());
            }
            /*cout << "transform rotation to element local system" << endl;
             cout << "node1Rot" << endl;
             cout << node1Rot[i * 3 + 0] << " " << node1Rot[i * 3 + 1] << " "
             << node1Rot[i * 3 + 2] << endl;
             cout << "node2Rot" << endl;
             cout << node2Rot[i * 3 + 0] << " " << node2Rot[i * 3 + 1] << " "
             << node2Rot[i * 3 + 2] << endl;*/

            // step 2: compute rotations without corotate
            ROT_Eo_Ed__node1.addRotation(ROT_Er_Eo->getRotationMatrix()); // o_rd -> o_d
            ROT_Eo_Ed__node2.addRotation(ROT_Er_Eo->getRotationMatrix()); // o_rd -> o_d
            /*cout << "compute rotations without corotate" << endl;
             cout << "node1Rot" << endl;
             cout << node1Rot[i * 3 + 0] << " " << node1Rot[i * 3 + 1] << " "
             << node1Rot[i * 3 + 2] << endl;
             cout << "node2Rot" << endl;
             cout << node2Rot[i * 3 + 0] << " " << node2Rot[i * 3 + 1] << " "
             << node2Rot[i * 3 + 2] << endl;*/

            // step 3: compute rotations without nonlinear torsion
            ROT_Eo_Ed__node1.getRotationVector(&node1Rot[i * 3]);
            ROT_Eo_Ed__node2.getRotationVector(&node2Rot[i * 3]);
            double torsionAngleAverage = (node1Rot[i * 3] + node2Rot[i * 3]) / 2.0;
            double xAxis[] = { 1.0, 0.0, 0.0 };
            ROT_torsion[i] = new KinematicMotion;
            ROT_torsion[i]->addRotation(xAxis, true, torsionAngleAverage);
            KinematicMotion *ROT_torsion_inv = ROT_torsion[i]->newInverse();
            ROT_Eo_Ed__node1.addRotation(ROT_torsion_inv->getRotationMatrix());
            ROT_Eo_Ed__node2.addRotation(ROT_torsion_inv->getRotationMatrix());

            // rotations related to linear motion (bending and small torsion)
            ROT_Eo_Ed__node1.getRotationVector(&node1Rot[i * 3]);
            ROT_Eo_Ed__node2.getRotationVector(&node2Rot[i * 3]);

            /*cout << "rotations related to linear motion" << endl;
             cout << "node1Rot" << endl;
             cout << node1Rot[i * 3 + 0] << " " << node1Rot[i * 3 + 1] << " "
             << node1Rot[i * 3 + 2] << endl;
             cout << "node2Rot" << endl;
             cout << node2Rot[i * 3 + 0] << " " << node2Rot[i * 3 + 1] << " "
             << node2Rot[i * 3 + 2] << endl;*/

            // displacements relative to the element local system
            ROT_ELEM_O->move(&node1Disp[i * 3]);
            ROT_ELEM_O->move(&node2Disp[i * 3]);
            delete ROT_Er_Eo;
            delete ROT_torsion_inv;
        }
        delete ROT_ELEM_O;
        //cout << "node1Disp" << endl;
        //cout << node1Disp[i * 3 + 0] <<" "<< node1Disp[i * 3 + 1] <<" "<< node1Disp[i * 3 + 2] << endl;
        //cout << "node2Disp" << endl;
        //cout << node2Disp[i * 3 + 0] <<" "<< node2Disp[i * 3 + 1] <<" "<< node2Disp[i * 3 + 2] << endl;
        //cout << "node1Rot" << endl;
        //cout << node1Rot[i * 3 + 0] <<" "<< node1Rot[i * 3 + 1] <<" "<< node1Rot[i * 3 + 2] << endl;
        //cout << "node2Rot" << endl;
        //cout << node2Rot[i * 3 + 0] <<" "<< node2Rot[i * 3 + 1] <<" "<< node2Rot[i * 3 + 2] << endl;
    }

    // compute displacements of nodes according to rotations and displacements of sections
    for (int i = 0; i < surfaceNumSections; i++) {
        // interpolate section rot and disp on the element local system
        int elem = sectionToCurveElem[i];
        double tmpSectionDisp[3];
        double tmpSectionRot[3];

        double length = corotateElementLength[elem];
        double linearShapeFunc1 = shapeFuncOfSection[i * 10 + 0];
        double cubicShapeFuncDisp1 = shapeFuncOfSection[i * 10 + 1];
        double cubicShapeFuncRot1 = shapeFuncOfSection[i * 10 + 2] * length;
        double cubicShapeFuncDispDeriv1 = shapeFuncOfSection[i * 10 + 3] / length;
        double cubicShapeFuncRotDeriv1 = shapeFuncOfSection[i * 10 + 4];
        double linearShapeFunc2 = shapeFuncOfSection[i * 10 + 5 + 0];
        double cubicShapeFuncDisp2 = shapeFuncOfSection[i * 10 + 5 + 1];
        double cubicShapeFuncRot2 = shapeFuncOfSection[i * 10 + 5 + 2] * length;
        double cubicShapeFuncDispDeriv2 = shapeFuncOfSection[i * 10 + 5 + 3] / length;
        double cubicShapeFuncRotDeriv2 = shapeFuncOfSection[i * 10 + 5 + 4];

        // DOFs related to linear motion
        //disp_x = 0.0
        tmpSectionDisp[0] = 0.0;
        //disp_y = NCr_1*rot_z1 + NCr_2*rot_z2
        tmpSectionDisp[1] = cubicShapeFuncRot1 * node1Rot[elem * 3 + 2]
                + cubicShapeFuncRot2 * node2Rot[elem * 3 + 2];
        //disp_z = NCr_1*(-rot_y1) + NCr_2*(-rot_y2), since disp_z' = -rot_y
        tmpSectionDisp[2] = cubicShapeFuncRot1 * (-node1Rot[elem * 3 + 1])
                + cubicShapeFuncRot2 * (-node2Rot[elem * 3 + 1]);
        //rot_x = NL_1*rot_x1 + NL_2*rot_x2
        tmpSectionRot[0] = linearShapeFunc1 * node1Rot[elem * 3 + 0]
                + linearShapeFunc2 * node2Rot[elem * 3 + 0];
        //rot_y = - disp_z' = - [NCr'_1*(-rot_y1) + NCr'_2*(-rot_y2)]
        tmpSectionRot[1] = -(cubicShapeFuncRotDeriv1 * (-node1Rot[elem * 3 + 1])
                + cubicShapeFuncRotDeriv2 * (-node2Rot[elem * 3 + 1]));
        //rot_z = disp_y' = NCr'_1*rot_z1 + NCr'_2*rot_z2
        tmpSectionRot[2] = cubicShapeFuncRotDeriv1 * node1Rot[elem * 3 + 2]
                + cubicShapeFuncRotDeriv2 * node2Rot[elem * 3 + 2];

        KinematicMotion KM_Eo_Ed;
        double localX[] = { 1.0, 0.0, 0.0 };
        double localY[] = { 0.0, 1.0, 0.0 };
        double localZ[] = { 0.0, 0.0, 1.0 };
        KM_Eo_Ed.addRotation(localY, true, tmpSectionRot[1]);
        KM_Eo_Ed.addRotation(localX, true, tmpSectionRot[0]);
        KM_Eo_Ed.addRotation(localZ, true, tmpSectionRot[2]);
        KM_Eo_Ed.addTranslation(tmpSectionDisp);

        // Add big torsion
        KM_Eo_Ed.addRotation(ROT_torsion[elem]->getRotationMatrix());

        // Add corotate
        // disp = NL_1*disp + NL_2*disp
        for (int j = 0; j < 3; j++) {
            tmpSectionDisp[j] = linearShapeFunc1 * node1Disp[elem * 3 + j]
                    + linearShapeFunc2 * node2Disp[elem * 3 + j];
        }

        KinematicMotion KM_Eo_Erd;
        KM_Eo_Erd.addKinematicMotion(&KM_Eo_Ed);
        KM_Eo_Erd.addRotation(ROT_Eo_Er[elem]->getRotationMatrix());
        KM_Eo_Erd.addTranslation(tmpSectionDisp);

        /*cout << "tmpSectionDisp" << endl;
         cout << tmpSectionDisp[0] <<" "<< tmpSectionDisp[1] <<" "<< tmpSectionDisp[2] << endl;
         cout << "tmpSectionRot" << endl;
         cout << tmpSectionRot[0] <<" "<< tmpSectionRot[1] <<" "<< tmpSectionRot[2] << endl;*/

        KinematicMotion TRL_O_P;
        TRL_O_P.addTranslation(&sectionP[i * 3]);
        KinematicMotion *TRL_P_O = TRL_O_P.newInverse();

        KinematicMotion KM_Oo_Ord;
        KinematicMotion *ROT_ELEM_O = ROT_O_ELEM[elem]->newInverse();
        KM_Oo_Ord.addKinematicMotion(TRL_P_O);
        KM_Oo_Ord.addKinematicMotion(ROT_ELEM_O);
        KM_Oo_Ord.addKinematicMotion(&KM_Eo_Erd);
        KM_Oo_Ord.addKinematicMotion(ROT_O_ELEM[elem]);
        KM_Oo_Ord.addKinematicMotion(&TRL_O_P);

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
            KM_Oo_Ord.move(nodalNewCoor);
            for (int k = 0; k < 3; k++)
                surfaceDisp[pos * 3 + k] = nodalNewCoor[k] - nodalOldCoor[k];
        }

        // compute the final section displacement and rotation from KM_Oo_Ord
        KM_Oo_Ord.getRotationVector(&sectionRot[i * 3]);
    }

    delete[] node1Disp;
    delete[] node1Rot;
    delete[] node2Disp;
    delete[] node2Rot;

    for (int i = 0; i < curveNumElements; i++) {
        delete ROT_Eo_Er[i];
        delete ROT_torsion[i];
    }
    delete[] ROT_Eo_Er;
    delete[] ROT_torsion;
    delete[] corotateElementLength;
}

void CurveSurfaceMapper::conservativeMapping(const double *surfaceForce, double *curveForceMoment) {
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

        KinematicMotion ROT;
        double tmpSectionRot[3];
        for (int j = 0; j < 3; j++)
            tmpSectionRot[j] = sectionRot[i * 3 + j];
        double angle = normalizeVector(tmpSectionRot);

        ROT.addRotation(tmpSectionRot, true, angle);

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
            // distance from P to a section node in the original configuration
            double P2N[3];
            for (int k = 0; k < 3; k++) {
                P2N[k] = surfaceNodeCoors[pos * 3 + k] - sectionP[i * 3 + k];
            }

            // distance from P to a section node in the deformed configuration
            ROT.move(P2N);

            // force vector on the nodal
            double F[3];
            for (int k = 0; k < 3; k++) {
                F[k] = surfaceForce[pos * 3 + k];
            }
            // compute moment around P with the nodal force: M = R X F
            double M[3];
            M[0] = P2N[1] * F[2] - P2N[2] * F[1];
            M[1] = P2N[2] * F[0] - P2N[0] * F[2];
            M[2] = P2N[0] * F[1] - P2N[1] * F[0];

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

double CurveSurfaceMapper::normalizeVector(double *vector) {
    double length = sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]);
    if (length < 1E-20) // 1E-20 is a physical "small" length, or a physical small angle
        return 0.0;
    for (int j = 0; j < 3; j++)
        vector[j] /= length;
    return length;
}

} /* namespace EMPIRE */
