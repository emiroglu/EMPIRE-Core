#include "CurveSurfaceCorotate2DMapper.h"
#include "KinematicMotion.h"

#include <assert.h>
#include <map>
#include <iostream>
#include <math.h>

using namespace std;

namespace EMPIRE {

CurveSurfaceCorotate2DMapper::CurveSurfaceCorotate2DMapper(int _curveNumNodes,
        int _curveNumElements, const double *_curveNodeCoors, const int *_curveNodeIDs,
        const int *_curveElems, int _surfaceNumNodes, const double *_surfaceNodeCoors,
        int _surfaceNumSections, int _surfaceNumRootSectionNodes, int _surfaceNumNormalSectionNodes,
        int _surfaceNumTipSectionNodes, const double *rotation_O_Q, const double *translation_O_Q) :
        AbstractCurveSurfaceMapper(_curveNumNodes, _curveNumElements, _curveNodeCoors,
                _curveNodeIDs, _curveElems, _surfaceNumNodes, _surfaceNodeCoors,
                _surfaceNumSections, _surfaceNumRootSectionNodes, _surfaceNumNormalSectionNodes,
                _surfaceNumTipSectionNodes, rotation_O_Q, translation_O_Q) {
    /*
     * Coordinate systems:
     *   O --- global
     *   Q --- beam root (x is the direction orthogonal to root section, y is flapwise direction, z is edgewise direction)
     *         origin is an arbitrary point in the root section
     *         Qx,Qy,Qz = +-Ox or Oy or Oz!!!
     *   P --- origin is the cross point of the section with the beam element, orientation is the same as the global system
     */
    // construct KM_O_Q
    const double EPS = 1E-10;
    for (int i = 0; i < 9; i++) {
        double tmp = fabs(rotation_O_Q[i]);
        if (tmp > 0.5)
            assert(fabs(tmp - 1.0) < EPS); // should be 1
        else
            assert(fabs(tmp - 0.0) < EPS); // should be 0
    }
    KM_O_Q = new KinematicMotion();
    KM_O_Q->addRotation(rotation_O_Q);
    KM_O_Q->addTranslation(translation_O_Q);

    // curve node coordinates in Q
    KinematicMotion *KM_Q_O = KM_O_Q->newInverse();
    curveNodeCoorsInQ = new double[curveNumNodes * 3];
    for (int i = 0; i < curveNumNodes * 3; i++)
        curveNodeCoorsInQ[i] = _curveNodeCoors[i];
    for (int i = 0; i < curveNumNodes; i++) {
        KM_Q_O->move(&curveNodeCoorsInQ[i * 3]);
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
}

CurveSurfaceCorotate2DMapper::~CurveSurfaceCorotate2DMapper() {
    for (int i = 0; i < curveNumElements; i++) {
        delete ROT_Q_ELEM[i];
    }
    delete[] ROT_Q_ELEM;
    delete KM_O_Q;
    delete[] angle_Q_ELEM;
    delete[] curveNodeCoorsInQ;
}

void CurveSurfaceCorotate2DMapper::consistentMapping(const double *curveDispRot,
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
        tmpSectionDisp[1] = cubicShapeFuncRot1 * (node1Rot[elem * 3 + 2] - angle_elem_corotate[elem])
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

} /* namespace EMPIRE */
