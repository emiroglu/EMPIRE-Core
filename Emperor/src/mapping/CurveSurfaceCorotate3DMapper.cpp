#include "CurveSurfaceCorotate3DMapper.h"
#include "KinematicMotion.h"

#include <assert.h>
#include <map>
#include <iostream>
#include <math.h>

using namespace std;

namespace EMPIRE {

CurveSurfaceCorotate3DMapper::CurveSurfaceCorotate3DMapper(int _curveNumNodes,
        int _curveNumElements, const double *_curveNodeCoors, const int *_curveNodeIDs,
        const int *_curveElems, int _surfaceNumNodes, const double *_surfaceNodeCoors,
        int _surfaceNumSections, int _surfaceNumRootSectionNodes, int _surfaceNumNormalSectionNodes,
        int _surfaceNumTipSectionNodes, const double *rotation_O_Q, const double *translation_O_Q) :
        AbstractCurveSurfaceMapper(_curveNumNodes, _curveNumElements, _curveNodeCoors,
                _curveNodeIDs, _curveElems, _surfaceNumNodes, _surfaceNodeCoors,
                _surfaceNumSections, _surfaceNumRootSectionNodes, _surfaceNumNormalSectionNodes,
                _surfaceNumTipSectionNodes, rotation_O_Q, translation_O_Q) {

}

CurveSurfaceCorotate3DMapper::~CurveSurfaceCorotate3DMapper() {

}

void CurveSurfaceCorotate3DMapper::consistentMapping(const double *curveDispRot,
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

} /* namespace EMPIRE */
