#include "AbstractCurveSurfaceMapper.h"
#include "KinematicMotion.h"
#include <assert.h>
#include <map>
#include <iostream>
#include <math.h>

using namespace std;

namespace EMPIRE {

AbstractCurveSurfaceMapper::AbstractCurveSurfaceMapper(int _curveNumNodes, int _curveNumElements,
        const double *_curveNodeCoors, const int *_curveNodeIDs, const int *_curveElems,
        int _surfaceNumNodes, const double *_surfaceNodeCoors, int _surfaceNumSections,
        int _surfaceNumRootSectionNodes, int _surfaceNumNormalSectionNodes,
        int _surfaceNumTipSectionNodes, const double *rotation_O_Q, const double *translation_O_Q) :
        curveNumNodes(_curveNumNodes), curveNumElements(_curveNumElements), curveNodeCoors(
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

    // check number of surface nodes
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
    double *curveNodeCoorsInQ = new double[curveNumNodes * 3];
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

    delete KM_Q_O;
    delete[] surfaceNodeCoorsInQ;
    delete rightXToElemPos;
    delete KM_O_Q;
    delete[] curveNodeCoorsInQ;
}

AbstractCurveSurfaceMapper::~AbstractCurveSurfaceMapper() {
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
}

void AbstractCurveSurfaceMapper::conservativeMapping(const double *surfaceForce,
        double *curveForceMoment) {
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
        for (int j=0; j<3; j++)
            tmpSectionRot[j] = sectionRot[i*3 +j];
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

double AbstractCurveSurfaceMapper::normalizeVector(double *vector) {
    double length = sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]);
    if (length < 1E-20) // 1E-20 is a physical "small" length, or a physical small angle
        return 0.0;
    for (int j = 0; j < 3; j++)
        vector[j] /= length;
    return length;
}

} /* namespace EMPIRE */
