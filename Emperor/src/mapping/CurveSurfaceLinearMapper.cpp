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
        AbstractCurveSurfaceMapper(_curveNumNodes, _curveNumElements, _curveNodeCoors,
                _curveNodeIDs, _curveElems, _surfaceNumNodes, _surfaceNodeCoors,
                _surfaceNumSections, _surfaceNumRootSectionNodes, _surfaceNumNormalSectionNodes,
                _surfaceNumTipSectionNodes, rotation_O_Q, translation_O_Q) {
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
}

CurveSurfaceLinearMapper::~CurveSurfaceLinearMapper() {
    delete[] curveElemLength;
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
}

} /* namespace EMPIRE */
