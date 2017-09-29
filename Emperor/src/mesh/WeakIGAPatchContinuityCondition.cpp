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
#include <string>

// Inclusion of user defined libraries
#include "IGAPatchCurve.h"
#include "IGAPatchSurface.h"
#include "WeakIGAPatchContinuityCondition.h"
#include "MathLibrary.h"
#include "Message.h"

using namespace std;

namespace EMPIRE {

WeakIGAPatchContinuityCondition::WeakIGAPatchContinuityCondition(int _ID,
                                                                 int _masterPatchIndex, int _pMaster, int _uNoKnotsMaster, double* _uKnotVectorMaster, int _uNoControlPointsMaster, double* _controlPointNetMaster,
                                                                 int _slavePatchIndex,  int _pSlave, int _uNoKnotsSlave, double* _uKnotVectorSlave, int _uNoControlPointsSlave, double* _controlPointNetSlave) :
    AbstractCondition(_ID),
    masterPatchIndex(_masterPatchIndex), slavePatchIndex(_slavePatchIndex)
{

    type = EMPIRE_WeakIGAPatchContinuityCondition;

    isTrimmingCurve = false;

    // Initialize the masterCurve and the slaveCurve
    masterCurve = new IGAPatchCurve(0, _pMaster, _uNoKnotsMaster, _uKnotVectorMaster, _uNoControlPointsMaster, _controlPointNetMaster);
    slaveCurve = new IGAPatchCurve(0, _pSlave, _uNoKnotsSlave, _uKnotVectorSlave, _uNoControlPointsSlave, _controlPointNetSlave);

    // Linearize the curves with the default algorithm and direction (combined algorithm and forward direction)
    masterCurve->linearize();
    slaveCurve->linearize();

    isGPDataInitialized = false;
}

WeakIGAPatchContinuityCondition::WeakIGAPatchContinuityCondition(int _ID,
                                                                 int _masterPatchIndex, IGAPatchCurve* _masterCurve,
                                                                 int _slavePatchIndex, IGAPatchCurve* _slaveCurve) :
    AbstractCondition(_ID),
    masterPatchIndex(_masterPatchIndex), masterCurve(_masterCurve),
    slavePatchIndex(_slavePatchIndex), slaveCurve(_slaveCurve)
{

    type = EMPIRE_WeakIGAPatchContinuityCondition;

    isTrimmingCurve = false;

    // Linearize the curves with the default algorithm and direction (combined algorithm and forward direction)
    masterCurve->linearize();
    slaveCurve->linearize();

    isGPDataInitialized = false;
}

WeakIGAPatchContinuityCondition::WeakIGAPatchContinuityCondition(int _ID,
                                                                 int _masterPatchIndex, int _masterPatchBLIndex, int _masterPatchBLTrCurveIndex,
                                                                 int _slavePatchIndex, int _slavePatchBLIndex, int _slavePatchBLTrCurveIndex) :
    AbstractCondition(_ID),
    masterPatchIndex(_masterPatchIndex), masterPatchBLIndex(_masterPatchBLIndex), masterPatchBLTrCurveIndex(_masterPatchBLTrCurveIndex),
    slavePatchIndex(_slavePatchIndex), slavePatchBLIndex(_slavePatchBLIndex), slavePatchBLTrCurveIndex(_slavePatchBLTrCurveIndex)
{

    type = EMPIRE_WeakIGAPatchContinuityCondition;

    isTrimmingCurve = true;

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

    // Set the master and slave curves in case they are trimming curves of the given patches
    if (isTrimmingCurve) {
        masterCurve = &_masterPatch->getTrimming().getLoop(masterPatchBLIndex).getIGACurve(masterPatchBLTrCurveIndex);
        slaveCurve = &_slavePatch->getTrimming().getLoop(slavePatchBLIndex).getIGACurve(slavePatchBLTrCurveIndex);
    }

    // Initialize the coordinates
    const int noCoordParam = 2;
    const int noCoord = 3;

    // Initialize auxiliary variables
    bool isProjectedOnSlave = false;

    // Initialize parameter and coordinate sets to store the intersections and transfer them
    std::vector<double> masterUTildes;
    std::vector<double> slaveUTildes;
    std::vector<double> masterUTildesFromSlave;
    std::vector<double> masterUTildesMerged;

    // Initialize temporary variables to use recursively in the loops
    double tmpUTilde;
    double tmpUVMaster[noCoordParam];
    double tmpUVSlave[noCoordParam];
    double tmpXYZ[noCoord];

    // Compute the knot intersections of the master trimming curve
    _masterPatch->computeKnotIntersectionsWithTrimmingCurve(masterUTildes, masterCurve);

    // Compute the knot intersections of slave the trimming curve
    _slavePatch->computeKnotIntersectionsWithTrimmingCurve(slaveUTildes, slaveCurve);

    // Project the knot intersections onto the master curve
    for (std::vector<double>::iterator iUTilde = slaveUTildes.begin(); iUTilde != slaveUTildes.end(); iUTilde++) {

        // Compute the physical coordinates of the knot intersection
        _slavePatch->computeCartesianCoordinates(tmpXYZ, *iUTilde, slaveCurve);

        // Compute the point projection of the slave knot intersection on the master trimming curve
        isProjectedOnSlave = false;
        isProjectedOnSlave = _masterPatch->computePointProjectionOnTrimmingCurve(tmpUTilde, tmpXYZ, masterCurve);

        // Consider the intersection only if the projection is successful
        if (isProjectedOnSlave)     masterUTildesFromSlave.push_back(tmpUTilde);
        else    WARNING_OUT("In \"WeakIGAPatchContinuityCondition::createGPData\"; a knot intersection could not be projected and discarded!");
    }

    // Merge the knot intersections. Here an assert for the vector sizes is not necessary since at least the beginning and the end knots are added
    EMPIRE::MathLibrary::mergeSortRemoveDuplicates(masterUTildesMerged, masterUTildes, masterUTildesFromSlave);

    // Second check to kick out very close knots
    EMPIRE::MathLibrary::sortRemoveSimilar(masterUTildesMerged, 1e-4);

    /// Create the Gauss points on the master and the slave sides
    // Getting the polynomial orders
    int pMaster = _masterPatch->getIGABasis(0)->getPolynomialDegree();
    int qMaster = _masterPatch->getIGABasis(1)->getPolynomialDegree();
    int pSlave = _slavePatch->getIGABasis(0)->getPolynomialDegree();
    int qSlave = _slavePatch->getIGABasis(1)->getPolynomialDegree();
    int p = std::max(std::max(pMaster,qMaster),std::max(pSlave,qSlave));

    // Get the number of Gauss points
    int numGPPerSection = p+1;

    // Create a Gauss quadrature rule
    MathLibrary::IGAGaussQuadratureOnBiunitInterval* theGPQuadrature = new MathLibrary::IGAGaussQuadratureOnBiunitInterval(numGPPerSection);

    // Initialize auxiliary variables
    double GP;
    double GW;
    double masterGPUTilde;
    double slaveGPUTilde;
    int counterValidGP = 0;

    std::vector<double> masterValidGPUTildes;
    std::vector<double> slaveValidGPUTildes;
    std::vector<double> validGPWeights;
    std::vector<int> validGPsPerSection;

    // Initialize GP counter on the curve
    trCurveNumGP = 0;

    /// Define and count the valid GPs
    // Loop over the sections of the trimming curve
    for (int iSection = 0; iSection < masterUTildesMerged.size() - 1; iSection++) {

        // Reset the valid GP per section
        counterValidGP = 0;

        // Loop over the GPs of the section
        for (int iGP = 0; iGP < numGPPerSection; iGP++) {

            // Get GP coordinates and weights
            GP = *theGPQuadrature->getGaussPoint(iGP);
            GW = theGPQuadrature->getGaussWeight(iGP);

            // Compute the image of the GP in the curve parameter space
            masterGPUTilde = ((1.0 - GP)*masterUTildesMerged[iSection] + (1.0 + GP)*masterUTildesMerged[iSection + 1])/2.0;

            // Compute the physical coordinates of the GP
            _masterPatch->computeCartesianCoordinates(tmpXYZ, masterGPUTilde, masterCurve);

            // Project the GP onto the slave side
            isProjectedOnSlave = false;
            isProjectedOnSlave = _slavePatch->computePointProjectionOnTrimmingCurve(slaveGPUTilde, tmpXYZ, slaveCurve);

            // If the GP is projected onto the slave side successfully then store the GP pair
            if (isProjectedOnSlave) {

                // Store the valid GP curve parameters
                masterValidGPUTildes.push_back(masterGPUTilde);
                slaveValidGPUTildes.push_back(slaveGPUTilde);

                // Store the valid GP weights
                validGPWeights.push_back(GW);

                // Update Gauss point counters
                counterValidGP++;
                trCurveNumGP++;

            } else  WARNING_OUT("In \"WeakIGAPatchContinuityCondition::createGPData\"; an interface GP could not be projected and discarded!");
        }

        // Store number of valid GPs per section
        validGPsPerSection.push_back(counterValidGP);
    }

    // Initialize the GP data
    trCurveGPWeights = new double[trCurveNumGP];
    trCurveGPJacobianProducts = new double[trCurveNumGP];
    trCurveMasterGPs = new double[trCurveNumGP*noCoordParam];
    trCurveMasterGPTangents = new double[trCurveNumGP*noCoord];
    trCurveSlaveGPs = new double[trCurveNumGP*noCoordParam];
    trCurveSlaveGPTangents = new double[trCurveNumGP*noCoord];

    /// Compute the GP data for the valid GPs
    // Initialize variables
    int noDeriv = 1;
    int noDerivBaseVct = 0;
    int derivDegree = 1;
    int uKnotSpanMaster;
    int vKnotSpanMaster;
    int uKnotSpanSlave;
    int vKnotSpanSlave;
    int knotSpanIndexTrCurveMaster;
    int knotSpanIndexTrCurveSlave;
    int pTrCurveMaster = masterCurve->getIGABasis()->getPolynomialDegree();
    int pTrCurveSlave = slaveCurve->getIGABasis()->getPolynomialDegree();
    int noLocalBasisFunctionsTrCurveMaster = pTrCurveMaster + 1;
    int noLocalBasisFunctionsMaster = (pMaster + 1)*(qMaster + 1);
    int noLocalBasisFunctionsTrCurveSlave = pTrCurveSlave + 1;
    int noLocalBasisFunctionsSlave = (pSlave + 1)*(qSlave + 1);

    double detJ1Master;
    double detJ2Master;
    double detJ2Slave;

    // Initialize pointers
    double baseVectorTrCurveMaster[(noDerivBaseVct + 1)*noCoord];
    double baseVectorTrCurveSlave[(noDerivBaseVct + 1)*noCoord];
    double localBasisFunctionsAndDerivativesTrCurveMaster[(derivDegree + 1) * (derivDegree + 2) * noLocalBasisFunctionsTrCurveMaster / 2];
    double localBasisFunctionsAndDerivativesTrCurveSlave[(derivDegree + 1) * (derivDegree + 2) * noLocalBasisFunctionsTrCurveSlave / 2];
    double localBasisFunctionsAndDerivativesMaster[(derivDegree + 1) * (derivDegree + 2) * noLocalBasisFunctionsMaster / 2];
    double localBasisFunctionsAndDerivativesSlave[(derivDegree + 1) * (derivDegree + 2) * noLocalBasisFunctionsSlave / 2];
    double baseVectorsMaster[6];
    double baseVectorsSlave[6];
    double A1Master[3];
    double A2Master[3];
    double A1Slave[3];
    double A2Slave[3];

    int counterGP = 0;
    // Loop over the sections
    for (int iSection = 0; iSection < masterUTildesMerged.size() - 1; iSection++) {

        // Determinant of the Jacobian of the transformation from the parent space to the parameter space of the patch
        detJ1Master = (masterUTildesMerged[iSection+1]-masterUTildesMerged[iSection])/2.0;

        // Loop over the GPs in this section
        for (int iGP = 0; iGP < validGPsPerSection[iSection]; iGP++) {

            // Get the curve parameters of the GPs
            masterGPUTilde = masterValidGPUTildes[counterGP];
            slaveGPUTilde = slaveValidGPUTildes[counterGP];

            // Find the knot spans in the parameter space of the curves
            knotSpanIndexTrCurveMaster = masterCurve->getIGABasis()->findKnotSpan(masterGPUTilde);
            knotSpanIndexTrCurveSlave = slaveCurve->getIGABasis()->findKnotSpan(slaveGPUTilde);

            // Compute the basis functions at the parametric locations of the trimming curves
            masterCurve->getIGABasis()->computeLocalBasisFunctionsAndDerivatives
                    (localBasisFunctionsAndDerivativesTrCurveMaster, noDeriv, masterGPUTilde, knotSpanIndexTrCurveMaster);
            slaveCurve->getIGABasis()->computeLocalBasisFunctionsAndDerivatives
                    (localBasisFunctionsAndDerivativesTrCurveSlave, noDeriv, slaveGPUTilde, knotSpanIndexTrCurveSlave);

            // Get the corresponding u,v parameters of the parametric locations
            masterCurve->computeCartesianCoordinates
                    (tmpUVMaster, localBasisFunctionsAndDerivativesTrCurveMaster, knotSpanIndexTrCurveMaster);
            slaveCurve->computeCartesianCoordinates
                    (tmpUVSlave, localBasisFunctionsAndDerivativesTrCurveSlave, knotSpanIndexTrCurveSlave);
            for (int iCoord = 0; iCoord < noCoordParam; iCoord++) {
                _masterPatch->getIGABasis(iCoord)->clampKnot(tmpUVMaster[iCoord]);
                _slavePatch->getIGABasis(iCoord)->clampKnot(tmpUVSlave[iCoord]);
                trCurveMasterGPs[noCoordParam*counterGP + iCoord] = tmpUVMaster[iCoord];
                trCurveSlaveGPs[noCoordParam*counterGP + iCoord] = tmpUVSlave[iCoord];
            }

            // Compute the base vectors of the trimming curves
            masterCurve->computeBaseVectorAndDerivatives
                    (baseVectorTrCurveMaster, knotSpanIndexTrCurveMaster, localBasisFunctionsAndDerivativesTrCurveMaster, noDerivBaseVct);
            slaveCurve->computeBaseVectorAndDerivatives
                    (baseVectorTrCurveSlave, knotSpanIndexTrCurveSlave, localBasisFunctionsAndDerivativesTrCurveSlave, noDerivBaseVct);

            // Find the knot span indices on the patch
            uKnotSpanMaster = _masterPatch->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(tmpUVMaster[0]);
            vKnotSpanMaster = _masterPatch->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(tmpUVMaster[1]);
            uKnotSpanSlave = _slavePatch->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(tmpUVSlave[0]);
            vKnotSpanSlave = _slavePatch->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(tmpUVSlave[1]);

            // Compute the basis functions of the patches at the (u,v) parametric locations
            _masterPatch->getIGABasis()->computeLocalBasisFunctionsAndDerivatives
                    (localBasisFunctionsAndDerivativesMaster, derivDegree, tmpUVMaster[0], uKnotSpanMaster, tmpUVMaster[1], vKnotSpanMaster);
            _slavePatch->getIGABasis()->computeLocalBasisFunctionsAndDerivatives
                    (localBasisFunctionsAndDerivativesSlave, derivDegree, tmpUVSlave[0], uKnotSpanSlave, tmpUVSlave[1], vKnotSpanSlave);

            // Compute the base vectors on the patches in the physical space
            _masterPatch->computeBaseVectors(baseVectorsMaster, localBasisFunctionsAndDerivativesMaster, uKnotSpanMaster, vKnotSpanMaster);
            _slavePatch->computeBaseVectors(baseVectorsSlave, localBasisFunctionsAndDerivativesSlave, uKnotSpanSlave, vKnotSpanSlave);
            for(int iCoord = 0; iCoord < noCoord; iCoord++){
                A1Master[iCoord] = baseVectorsMaster[iCoord];
                A2Master[iCoord] = baseVectorsMaster[noCoord+iCoord];
                A1Slave[iCoord] = baseVectorsSlave[iCoord];
                A2Slave[iCoord] = baseVectorsSlave[noCoord+iCoord];
            }

            // Compute the tangent vector on the physical space
            MathLibrary::computeDenseVectorMultiplicationScalar(A1Master, baseVectorTrCurveMaster[0], noCoord);
            MathLibrary::computeDenseVectorMultiplicationScalar(A2Master, baseVectorTrCurveMaster[1], noCoord);
            MathLibrary::computeDenseVectorAddition(A1Master, A2Master, 1.0, noCoord);
            MathLibrary::computeDenseVectorMultiplicationScalar(A1Slave, baseVectorTrCurveSlave[0], noCoord);
            MathLibrary::computeDenseVectorMultiplicationScalar(A2Slave, baseVectorTrCurveSlave[1], noCoord);
            MathLibrary::computeDenseVectorAddition(A1Slave, A2Slave, 1.0, noCoord);

            // Compute the determinant of the Jacobian of the transformation from the parameter space to the physical space
            detJ2Master = MathLibrary::vector2norm(A1Master, noCoord);
            detJ2Slave = MathLibrary::vector2norm(A1Slave, noCoord);

            // Normalize the tangent vectors
            MathLibrary::computeDenseVectorMultiplicationScalar(A1Master, 1.0/detJ2Master, noCoord);
            MathLibrary::computeDenseVectorMultiplicationScalar(A1Slave, 1.0/detJ2Slave, noCoord);

            for(int iCoord = 0; iCoord < noCoord; iCoord++) {
                trCurveMasterGPTangents[noCoord*counterGP + iCoord] = A1Master[iCoord];
                trCurveSlaveGPTangents[noCoord*counterGP + iCoord] = A1Slave[iCoord];
            }

            // Store the GP Weights
            trCurveGPWeights[counterGP] = validGPWeights[counterGP];

            // Store the Jacobian and Gauss weight products
            trCurveGPJacobianProducts[counterGP] = detJ1Master * detJ2Master * trCurveGPWeights[counterGP];

            // Update the GP counter
            counterGP++;
        }
    }

    // Set the initialized flag to true
    isGPDataInitialized = true;

    // Delete pointers
    delete theGPQuadrature;

}

WeakIGAPatchContinuityCondition::~WeakIGAPatchContinuityCondition() {
    if (!isTrimmingCurve) {
        delete masterCurve;
        delete slaveCurve;
    }
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
