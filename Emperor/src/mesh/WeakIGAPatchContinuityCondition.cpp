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

// Inclusion of standard libraries
#include <iostream>
#include <assert.h>
#include <math.h>
#include <algorithm>

// Inclusion of user defined libraries
#include "IGAPatchSurface.h"
#include "WeakIGAPatchContinuityCondition.h"
#include "MathLibrary.h"
#include "Message.h"

using namespace std;

namespace EMPIRE {

WeakIGAPatchContinuityCondition::WeakIGAPatchContinuityCondition(int _ID, 
                                                                 int _masterPatchIndex, int _masterPatchBLIndex, int _masterPatchBLTrCurveIndex,
                                                                 int _slavePatchIndex, int _slavePatchBLIndex, int _slavePatchBLTrCurveIndex) :
    AbstractCondition(_ID),
    masterPatchIndex(_masterPatchIndex), masterPatchBLIndex(_masterPatchBLIndex), masterPatchBLTrCurveIndex(_masterPatchBLTrCurveIndex),
    slavePatchIndex(_slavePatchIndex), slavePatchBLIndex(_slavePatchBLIndex), slavePatchBLTrCurveIndex(_slavePatchBLTrCurveIndex)
{

    type = EMPIRE_WeakIGAPatchContinuityCondition;

    isGPDataInitialized = false;
}

void WeakIGAPatchContinuityCondition::addWeakContinuityConditionGPData(int _trCurveNumGP,
                                                                       double* _trCurveMasterGPs, double* _trCurveSlaveGPs, double* _trCurveGPWeights,
                                                                       double* _trCurveMasterGPTangents, double* _trCurveSlaveGPTangents,
                                                                       double* _trCurveGPJacobianProducts){

    if (isGPDataInitialized) assert(false);

    const int noCoordParam = 2;
    const int noCoord = 3;

    trCurveNumGP = _trCurveNumGP;

    trCurveMasterGPs = new double[trCurveNumGP*noCoordParam];
    trCurveSlaveGPs = new double[trCurveNumGP*noCoordParam];
    trCurveGPWeights = new double[trCurveNumGP];
    trCurveMasterGPTangents = new double[trCurveNumGP*noCoord];
    trCurveSlaveGPTangents = new double[trCurveNumGP*noCoord];
    trCurveGPJacobianProducts = new double[trCurveNumGP];

    for (int i = 0; i < trCurveNumGP; i++){
        for (int j = 0; j < noCoordParam; j++){
            trCurveMasterGPs[i*noCoordParam+j]=_trCurveMasterGPs[i*noCoordParam+j];
            trCurveSlaveGPs[i*noCoordParam+j]=_trCurveSlaveGPs[i*noCoordParam+j];
        }
        trCurveGPWeights[i] = _trCurveGPWeights[i];
        for (int j = 0; j < noCoord; j++){
            trCurveMasterGPTangents[i*noCoord+j]=_trCurveMasterGPTangents[i*noCoord+j];
            trCurveSlaveGPTangents[i*noCoord+j]=_trCurveSlaveGPTangents[i*noCoord+j];
        }
        trCurveGPJacobianProducts[i] = _trCurveGPJacobianProducts[i];
    }

    isGPDataInitialized = true;
}

void WeakIGAPatchContinuityCondition::createGPData(IGAPatchSurface* _masterPatch, IGAPatchSurface* _slavePatch) {

    if (isGPDataInitialized) assert(false);

    const int noCoordParam = 2;
    const int noCoord = 3;

    // Initialize parameter and coordinate sets to store the intersections and transfer them
    std::vector<double> masterUTilde;
    std::vector<double> masterUV;
    std::vector<double> masterXYZ;
    std::vector<double> slaveUTilde;
    std::vector<double> slaveUV;
    std::vector<double> slaveXYZ;
    std::vector<double> masterUTildeFromSlave;
    std::vector<double> masterUVFromSlave;
    std::vector<double> masterXYZFromSlave;

    // Compute the knot intersections of the master trimming curve and their respective parameters/coordinates on the master patch
    _masterPatch->computeKnotIntersectionsWithTrimmingCurve(masterUTilde, masterUV, masterXYZ,
                                                            masterPatchBLIndex, masterPatchBLTrCurveIndex);

    // Compute the knot intersections of slave the trimming curve and their respective parameters/coordinates on the slave patch
    _slavePatch->computeKnotIntersectionsWithTrimmingCurve(slaveUTilde, slaveUV, slaveXYZ,
                                                           slavePatchBLIndex, slavePatchBLTrCurveIndex);

    // Compute the coordinates of the knot intersections of the master trimming curve on the slave trimming curve and slave patch
    _masterPatch->computePointProjectionOnTrimmingCurve(masterUTildeFromSlave, masterUVFromSlave, masterXYZFromSlave,
                                                        slaveXYZ, masterPatchBLIndex, masterPatchBLTrCurveIndex);

    // Sort the masterUTildeFromSlave, masterUVFromSlave and masterXYZFromSlave ascending order with respect to masterUTildeFromSlave
    if ( *masterUTildeFromSlave.begin() > *masterUTildeFromSlave.end()){
        std::sort(masterUTildeFromSlave.begin(),masterUTildeFromSlave.end());
        masterUVFromSlave.clear();
        masterXYZFromSlave.clear();
        for (int iUTilde = 0; iUTilde < masterUTildeFromSlave.size(); iUTilde++) {
            double uv[2];
            double xyz[3];
            _masterPatch->getTrimming().getLoop(masterPatchBLIndex).getIGACurve(masterPatchBLTrCurveIndex).computeCartesianCoordinates(uv, masterUTildeFromSlave[iUTilde]);
            masterUVFromSlave.push_back(uv[0]);
            masterUVFromSlave.push_back(uv[1]);
            _masterPatch->computeCartesianCoordinates(xyz, uv);
            masterXYZFromSlave.push_back(xyz[0]);
            masterXYZFromSlave.push_back(xyz[1]);
            masterXYZFromSlave.push_back(xyz[2]);
        }
    }

    // Sort the uTilde, uvParams and xyzCoords of both intersections in ascending order with respect to UTilde

    // Initiate iterators
    std::vector<double>::iterator iMasterUTilde = masterUTilde.begin();
    std::vector<double>::iterator iMasterUTildeFromSlave = masterUTildeFromSlave.begin();

    // Loop over each uTilde of master curve knot intersections on the slave curve
    while (iMasterUTildeFromSlave != masterUTildeFromSlave.end()){

        // Loop over each uTilde of slave curve knot intersections
        while (iMasterUTilde != masterUTilde.end()) {

            // If the slaveUTildeFromMaster is greater, advance the iterator on the slaveUTilde
            if (*iMasterUTildeFromSlave > *iMasterUTilde)  iMasterUTilde++;

            // If the uTildes are equal break loop
            if (fabs (*iMasterUTildeFromSlave - *iMasterUTilde) < 1e-2)   break;

            // If the slaveUTildeFromMaster is smaller than slaveUTilde, insert into place
            else if (*iMasterUTildeFromSlave < *iMasterUTilde && fabs(*iMasterUTildeFromSlave - *iMasterUTilde) > 1e-6){

                // Compute index to insert uTilde as an integer
                int indexUTilde = std::distance(masterUTilde.begin(),iMasterUTilde);

                // Compute index to insert uvParam as an integer
                int indexUV = indexUTilde*noCoordParam;

                // Compute index to insert xyzCoords as an integer
                int indexXYZCoords = indexUTilde*noCoord;

                // Compute index to retrieve from slaveUVFromMaster as an integer
                int indexUVFromSlave = std::distance(masterUTildeFromSlave.begin(),iMasterUTildeFromSlave)*noCoordParam;

                // Compute index to retrieve from slaveXYZFromMaster as an integer
                int indexXYZFromSlave = std::distance(masterUTildeFromSlave.begin(),iMasterUTildeFromSlave)*noCoord;

                // Insert UTilde
                masterUTilde.insert(iMasterUTilde,*iMasterUTildeFromSlave);

                // Insert UV
                masterUV.insert(masterUV.begin()+indexUV,*(masterUTildeFromSlave.begin()+indexUVFromSlave+1));
                masterUV.insert(masterUV.begin()+indexUV,*(masterUTildeFromSlave.begin()+indexUVFromSlave));

                // Insert XYZ
                masterXYZ.insert(masterXYZ.begin()+indexXYZCoords,*(masterXYZFromSlave.begin()+indexXYZFromSlave+2));
                masterXYZ.insert(masterXYZ.begin()+indexXYZCoords,*(masterXYZFromSlave.begin()+indexXYZFromSlave+1));
                masterXYZ.insert(masterXYZ.begin()+indexXYZCoords,*(masterXYZFromSlave.begin()+indexXYZFromSlave));

                // Reset the iterator after each insert
                iMasterUTilde = masterUTilde.begin();
            }
            iMasterUTilde++;
        }
        iMasterUTildeFromSlave++;
    }

    // Getting the polynomial orders
    int pMaster = _masterPatch->getIGABasis(0)->getPolynomialDegree();
    int qMaster = _masterPatch->getIGABasis(1)->getPolynomialDegree();
    int pSlave = _slavePatch->getIGABasis(0)->getPolynomialDegree();
    int qSlave = _slavePatch->getIGABasis(1)->getPolynomialDegree();
    int p = std::max(std::max(pMaster,qMaster),std::max(pSlave,qSlave));

    // Get the number of Gauss points
    int numGPPerSection = p+1;
    trCurveNumGP = numGPPerSection * (masterUTilde.size()-1);

    // Create a Gauss quadrature rule
    MathLibrary::IGAGaussQuadratureOnBiunitInterval* theGPQuadrature = new MathLibrary::IGAGaussQuadratureOnBiunitInterval(numGPPerSection);

    // Initialize the GP data
    trCurveGPWeights = new double[trCurveNumGP];
    trCurveGPJacobianProducts = new double[trCurveNumGP];
    trCurveMasterGPs = new double[trCurveNumGP*noCoordParam];
    trCurveMasterGPTangents = new double[trCurveNumGP*noCoord];
    trCurveSlaveGPs = new double[trCurveNumGP*noCoordParam];
    trCurveSlaveGPTangents = new double[trCurveNumGP*noCoord];

    // Initialize auxiliary variables
    double GP;
    double GW;
    double uTildeGP;
    double detJ1;
    double detJ2;
    int pTrCurveMaster = _masterPatch->getTrimming().getLoop(masterPatchBLIndex).getIGACurve(masterPatchBLTrCurveIndex).getIGABasis()->getPolynomialDegree();

    // Initialize pointers
    int noDeriv = 1;
    int noDerivBaseVct = 0;
    int derivDegree = 1;
    int noLocalBasisFunctionsTrCurveMaster = pTrCurveMaster + 1;
    int noLocalBasisFunctionsMaster = (pMaster + 1)*(qMaster + 1);
    int knotSpanIndexTrCurveMaster;
    int uKnotSpanMaster;
    int vKnotSpanMaster;
    int counterGP;
    double* uv = new double[noCoordParam];
    double* baseVectorTrCurveMaster = new double[(noDerivBaseVct + 1)*noCoord];
    double* localBasisFunctionsAndDerivativesMaster = new double[(derivDegree + 1) * (derivDegree + 2) * noLocalBasisFunctionsMaster / 2];
    double* localBasisFunctionsAndDerivativesTrCurveMaster = new double[noLocalBasisFunctionsTrCurveMaster*(noDeriv + 1)];
    double* baseVectorsMaster = new double[6];
    double* A1 = new double[3];
    double* A2 = new double[3];
    double* XYZ = new double[3];
    std::vector<double> masterGPUTilde;
    std::vector<double> masterGPXYZ;

    // Initialize GP counter on the curve
    counterGP = 0;

    // Loop over the sections of the trimming curve
    for (int iSection = 0; iSection < masterUTilde.size() - 1; iSection++) {

        // Determinant of the Jacobian of the transformation from the parent space to the parameter space of the patch
        detJ1 = (masterUTilde[iSection+1]-masterUTilde[iSection])/2.0;

        // Loop over the GPs of the section
        for (int iGP = 0; iGP < numGPPerSection; iGP++) {

            // Get GP coordinates and weights
            GP = theGPQuadrature->gaussPoints[iGP];
            GW = theGPQuadrature->weights[iGP];
            trCurveGPWeights[counterGP] = GW;

            // Compute the image of the GP in the curve parameter space
            uTildeGP = ((1.0 - GP)*masterUTilde[iSection] + (1.0 + GP)*masterUTilde[iSection + 1])/2.0;
            masterGPUTilde.push_back(uTildeGP);

            /// computeGPDataTangentVectorAndJacobianProducts

            // Find the knot span in the parameter space of the curve
            knotSpanIndexTrCurveMaster = _masterPatch->getTrimming().getLoop(masterPatchBLIndex).getIGACurve(masterPatchBLTrCurveIndex).getIGABasis()->findKnotSpan(uTildeGP);

            // Compute the basis functions at the parametric location of the trimming curve
            _masterPatch->getTrimming().getLoop(masterPatchBLIndex).getIGACurve(masterPatchBLTrCurveIndex).getIGABasis()->computeLocalBasisFunctionsAndDerivatives
                    (localBasisFunctionsAndDerivativesTrCurveMaster, noDeriv, uTildeGP, knotSpanIndexTrCurveMaster);

            // Get the corresponding u,v parameters of the uTildeGP parametric location
            _masterPatch->getTrimming().getLoop(masterPatchBLIndex).getIGACurve(masterPatchBLTrCurveIndex).computeCartesianCoordinates
                    (uv, localBasisFunctionsAndDerivativesTrCurveMaster, knotSpanIndexTrCurveMaster);
            trCurveMasterGPs[noCoordParam*counterGP] = uv[0];
            trCurveMasterGPs[noCoordParam*counterGP + 1] = uv[1];

            // Compute the base vector of the trimming curve in the master
            _masterPatch->getTrimming().getLoop(masterPatchBLIndex).getIGACurve(masterPatchBLTrCurveIndex).computeBaseVectorAndDerivatives
                    (baseVectorTrCurveMaster, knotSpanIndexTrCurveMaster, localBasisFunctionsAndDerivativesTrCurveMaster, noDerivBaseVct);

            // Find the knot span indices for on the master patch
            uKnotSpanMaster = _masterPatch->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(uv[0]);
            vKnotSpanMaster = _masterPatch->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(uv[1]);

            // Compute the basis functions of the master patch at the (u,v) parametric location
            _masterPatch->getIGABasis()->computeLocalBasisFunctionsAndDerivatives
                    (localBasisFunctionsAndDerivativesMaster, derivDegree, uv[0], uKnotSpanMaster, uv[1], vKnotSpanMaster);

            // Compute the base vectors on the master patch in the physical space
            _masterPatch->computeBaseVectors(baseVectorsMaster, localBasisFunctionsAndDerivativesMaster, uKnotSpanMaster, vKnotSpanMaster);

            // Compute the physical coordinates of the parametric location (u,v)
            _masterPatch->computeCartesianCoordinates(XYZ, localBasisFunctionsAndDerivativesMaster, derivDegree, uKnotSpanMaster, vKnotSpanMaster);
            for(int iCoord = 0; iCoord < noCoord; iCoord++)
                masterGPXYZ.push_back(XYZ[iCoord]);

            for(int iCoord = 0; iCoord < noCoord; iCoord++){
                A1[iCoord] = baseVectorsMaster[iCoord];
                A2[iCoord] = baseVectorsMaster[noCoord+iCoord];
            }

            // Compute the tangent vector on the physical space
            MathLibrary::computeDenseVectorMultiplicationScalar(A1, baseVectorTrCurveMaster[0], noCoord);
            MathLibrary::computeDenseVectorMultiplicationScalar(A2, baseVectorTrCurveMaster[1], noCoord);
            MathLibrary::computeDenseVectorAddition(A1, A2, 1.0, noCoord);

            // Compute the determinant of the Jacobian of the transformation from the parameter space to the physical space
            detJ2 = MathLibrary::vector2norm(A1, noCoord);
            MathLibrary::computeDenseVectorMultiplicationScalar(A1, 1.0/detJ2, noCoord);

            for(int iCoord = 0; iCoord < noCoord; iCoord++)
                trCurveMasterGPTangents[noCoord*counterGP + iCoord] = A1[iCoord];

            // Compute the element length on the Gauss point
            trCurveGPJacobianProducts[counterGP] = detJ1 * detJ2 * GW;

            // Update Gauss point counter
            counterGP++;
        }
    }

    // Initialize parameter and coordinate sets to store the GP coordinates
    std::vector<double> slaveGPUTilde;
    std::vector<double> slaveGPUV;
    std::vector<double> slaveGPXYZ;

    _slavePatch->computePointProjectionOnTrimmingCurve(slaveGPUTilde,slaveGPUV,slaveGPXYZ,
                                                       masterGPXYZ,slavePatchBLIndex,slavePatchBLTrCurveIndex);

    int knotSpanIndexTrCurveSlave;
    int uKnotSpanSlave;
    int vKnotSpanSlave;

    // Initialize auxiliary variables
    int pTrCurveSlave = _slavePatch->getTrimming().getLoop(slavePatchBLIndex).getIGACurve(slavePatchBLTrCurveIndex).getIGABasis()->getPolynomialDegree();
    int noLocBasisFunctionsTrCurveSlave = pTrCurveSlave + 1;
    int noLocalBasisFunctionsSlave = (pSlave + 1)*(qSlave + 1);

    // Initialize pointers
    double* localBasisFunctionsAndDerivativesTrCurveSlave = new double[noLocBasisFunctionsTrCurveSlave*(noDeriv + 1)];
    double* localBasisFunctionsAndDerivativesSlave = new double[(derivDegree + 1) * (derivDegree + 2) * noLocalBasisFunctionsSlave / 2];
    double* baseVectorTrCurveSlave = new double[noDeriv*noCoord];
    double* baseVectorsSlave = new double[6];

    // Loop over the GPs
    for (int iSlaveGP = 0; iSlaveGP < slaveGPUTilde.size(); iSlaveGP++) {

        // Uodate the curve parameter of the GP
        uTildeGP = slaveGPUTilde[iSlaveGP];

        // Find the knot span
        knotSpanIndexTrCurveSlave = _slavePatch->getTrimming().getLoop(slavePatchBLIndex).getIGACurve(slavePatchBLTrCurveIndex).getIGABasis()->findKnotSpan(uTildeGP);

        // Compute the basis functions at the parametric location
        _slavePatch->getTrimming().getLoop(slavePatchBLIndex).getIGACurve(slavePatchBLTrCurveIndex).getIGABasis()->computeLocalBasisFunctionsAndDerivatives
                (localBasisFunctionsAndDerivativesTrCurveSlave, noDeriv, uTildeGP, knotSpanIndexTrCurveSlave);

        // Get the corresponding u,v parameters of the uTildeGP parametric location
        _slavePatch->getTrimming().getLoop(slavePatchBLIndex).getIGACurve(slavePatchBLTrCurveIndex).computeCartesianCoordinates
                (uv, localBasisFunctionsAndDerivativesTrCurveSlave, knotSpanIndexTrCurveSlave);
        trCurveSlaveGPs[noCoordParam*iSlaveGP] = uv[0];
        trCurveSlaveGPs[noCoordParam*iSlaveGP + 1] = uv[1];

        // Compute the base vector
        _slavePatch->getTrimming().getLoop(slavePatchBLIndex).getIGACurve(slavePatchBLTrCurveIndex).computeBaseVectorAndDerivatives
                (baseVectorTrCurveSlave, knotSpanIndexTrCurveSlave, localBasisFunctionsAndDerivativesTrCurveSlave, noDerivBaseVct);

        // Find the knot span indices for on the slave patch
        uKnotSpanSlave = _slavePatch->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(uv[0]);
        vKnotSpanSlave = _slavePatch->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(uv[1]);

        // Compute the basis functions of the slave patch at the (u,v) parametric location
        _slavePatch->getIGABasis()->computeLocalBasisFunctionsAndDerivatives(
                    localBasisFunctionsAndDerivativesSlave, derivDegree, uv[0], uKnotSpanSlave, uv[1], vKnotSpanSlave);

        // Compute the base vectors on the slave patch
        _slavePatch->computeBaseVectors(baseVectorsSlave, localBasisFunctionsAndDerivativesSlave, uKnotSpanSlave, vKnotSpanSlave);
        for(int iCoord = 0; iCoord < noCoord; iCoord++){
            A1[iCoord] = baseVectorsSlave[iCoord];
            A2[iCoord] = baseVectorsSlave[noCoord+iCoord];
        }

        // Compute the tangent vector on the physical space
        MathLibrary::computeDenseVectorMultiplicationScalar(A1, baseVectorTrCurveSlave[0], noCoord);
        MathLibrary::computeDenseVectorMultiplicationScalar(A2, baseVectorTrCurveSlave[1], noCoord);
        MathLibrary::computeDenseVectorAddition(A1, A2, 1.0, noCoord);
        // Compute the determinant of the Jacobian of the transformation from the parameter space to the physical space
        detJ2 = MathLibrary::vector2norm(A1, noCoord);
        MathLibrary::computeDenseVectorMultiplicationScalar(A1, 1.0/detJ2, noCoord);

        for(int iCoord = 0; iCoord < noCoord; iCoord++)
            trCurveSlaveGPTangents[noCoord*iSlaveGP + iCoord] = A1[iCoord];

    }

    // Set the initialized flag to true
    isGPDataInitialized = true;

    // delete pointers
    delete theGPQuadrature;

    delete uv;
    delete XYZ;

    delete baseVectorsMaster;   delete baseVectorsSlave;
    delete A1;    delete A2;

    delete localBasisFunctionsAndDerivativesTrCurveMaster;
    delete localBasisFunctionsAndDerivativesMaster;
    delete baseVectorTrCurveMaster;

    delete localBasisFunctionsAndDerivativesTrCurveSlave;
    delete localBasisFunctionsAndDerivativesSlave;
    delete baseVectorTrCurveSlave;

}

WeakIGAPatchContinuityCondition::~WeakIGAPatchContinuityCondition() {
    if (isGPDataInitialized) {
        delete trCurveMasterGPs;
        delete trCurveSlaveGPs;
        delete trCurveGPWeights;
        delete trCurveMasterGPTangents;
        delete trCurveSlaveGPTangents;
        delete trCurveGPJacobianProducts;
    }
}


} /* namespace EMPIRE */
