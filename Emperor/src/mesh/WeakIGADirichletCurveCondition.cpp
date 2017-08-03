
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
#include "WeakIGADirichletCurveCondition.h"
#include "IGAPatchSurface.h"
#include "IGAPatchCurve.h"
#include "MathLibrary.h"
#include "Message.h"

using namespace std;

namespace EMPIRE {

WeakIGADirichletCurveCondition::WeakIGADirichletCurveCondition(int _ID,
                                                               int _patchIndex, int _p, int _uNoKnots, double* _uKnotVector, int _uNoControlPoints, double* _controlPointNet) :
    AbstractCondition(_ID),
    patchIndex(_patchIndex)
{

    type = EMPIRE_WeakIGADirichletCurveCondition;

    isTrimmingCurve = false;

    // Initialize the dirichletCurve
    dirichletCurve = new IGAPatchCurve(0, _p, _uNoKnots, _uKnotVector, _uNoControlPoints, _controlPointNet);

    // Linearize the curve with the default algorithm and direction (combined algorithm and forward direction)
    dirichletCurve->linearize();

    isGPDataInitialized = false;
}

WeakIGADirichletCurveCondition::WeakIGADirichletCurveCondition(int _ID,
                                 int _patchIndex, IGAPatchCurve* _dirichletCurve) :
        AbstractCondition(_ID),
        patchIndex(_patchIndex), dirichletCurve(_dirichletCurve)
        {

    type = EMPIRE_WeakIGADirichletCurveCondition;

    isTrimmingCurve = false;

    // Linearize the curve with the default algorithm and direction (combined algorithm and forward direction)
    dirichletCurve->linearize();

    isGPDataInitialized = false;

}

WeakIGADirichletCurveCondition::WeakIGADirichletCurveCondition(int _ID,
                                 int _patchIndex, int _patchBLIndex, int _patchBLTrCurveIndex) :
        AbstractCondition(_ID),
        patchIndex(_patchIndex), patchBLIndex(_patchBLIndex), patchBLTrCurveIndex(_patchBLTrCurveIndex)
        {

    type = EMPIRE_WeakIGADirichletCurveCondition;

    isTrimmingCurve = true;

    isGPDataInitialized = false;

}

void WeakIGADirichletCurveCondition::addWeakDirichletCurveConditionGPData(int _curveNumGP,
                            double* _curveGPs, double* _curveGPWeights,
                            double* _curveGPTangents,
                            double* _curveGPJacobianProducts){

    if (isGPDataInitialized) assert(false);

    const int noCoordParam = 2;
    const int noCoord = 3;

    curveNumGP = _curveNumGP;

    curveGPs = new double[curveNumGP*noCoordParam];
    curveGPWeights = new double[curveNumGP];
    curveGPTangents = new double[curveNumGP*noCoord];
    curveGPJacobianProducts = new double[curveNumGP];

    for (int i = 0; i < curveNumGP; i++){
        for (int j = 0; j < noCoordParam; j++){
            curveGPs[i*noCoordParam+j]=_curveGPs[i*noCoordParam+j];
        }
        curveGPWeights[i] = _curveGPWeights[i];
        for (int j = 0; j < noCoord; j++){
            curveGPTangents[i*noCoord+j]= _curveGPTangents[i*noCoord+j];
        }
        curveGPJacobianProducts[i] = _curveGPJacobianProducts[i];
    }

    isGPDataInitialized = true;
}

void WeakIGADirichletCurveCondition::createGPData(IGAPatchSurface* _patch) {

    if (isGPDataInitialized) assert(false);

    // Initialize the coordinates
    const int noCoordParam = 2;
    const int noCoord = 3;

    // Initialize parameter set to store the intersections and transfer them
    std::vector<double> UTildes;

    // Check to verify if the curve is a trimming curve and assign it to the member variable
    if (isTrimmingCurve)
        dirichletCurve = &_patch->getTrimming().getLoop(patchBLIndex).getIGACurve(patchBLTrCurveIndex);

    // Compute the knot intersections of the trimming curve
    _patch->computeKnotIntersectionsWithTrimmingCurve(UTildes, dirichletCurve);

    // Sort and remove duplicates if they exist
    EMPIRE::MathLibrary::sortRemoveDuplicates(UTildes);

    /// Create the Gauss points on the master and the slave sides
    // Getting the polynomial orders
    int p = _patch->getIGABasis(0)->getPolynomialDegree();
    int q = _patch->getIGABasis(1)->getPolynomialDegree();
    int pMax = std::max(p, q);

    // Get the number of Gauss points
    int numGPPerSection = pMax + 1;
    curveNumGP = numGPPerSection * (UTildes.size()-1);

    // Initialize the GP data
    curveGPWeights = new double[curveNumGP];
    curveGPJacobianProducts = new double[curveNumGP];
    curveGPs = new double[curveNumGP*noCoordParam];
    curveGPTangents = new double[curveNumGP*noCoord];

    // Create a Gauss quadrature rule
    MathLibrary::IGAGaussQuadratureOnBiunitInterval* theGPQuadrature = new MathLibrary::IGAGaussQuadratureOnBiunitInterval(numGPPerSection);

    // Initialize variables
    int noDeriv = 1;
    int noDerivBaseVct = 0;
    int derivDegree = 1;
    int uKnotSpan;
    int vKnotSpan;
    int knotSpanIndexCurve;
    int pCurve = dirichletCurve->getIGABasis()->getPolynomialDegree();
    int noLocalBasisFunctionsCurve = pCurve + 1;
    int noLocalBasisFunctions = (p + 1)*(q + 1);
    double GP;
    double GW;
    double GPUTilde;
    double detJ1;
    double detJ2;

    // Initialize pointers
    double uv[2];
    double baseVectorCurve[(noDerivBaseVct + 1)*noCoord];
    double localBasisFunctionsAndDerivatives[(derivDegree + 1) * (derivDegree + 2) * noLocalBasisFunctions / 2];
    double localBasisFunctionsAndDerivativesCurve[noLocalBasisFunctionsCurve * (noDeriv + 1)];
    double baseVectors[6];
    double A1[3];
    double A2[3];
    // Initialize GP counter on the curve
    int counterGP = 0;

    // Loop over the sections of the curve
        for (int iSection = 0; iSection < UTildes.size() - 1; iSection++) {

            // Determinant of the Jacobian of the transformation from the parent space to the parameter space of the patch
            detJ1 = (UTildes[iSection+1]-UTildes[iSection])/2.0;

            // Loop over the GPs of the section
            for (int iGP = 0; iGP < numGPPerSection; iGP++) {

                // Get GP coordinates and weights
                GP = theGPQuadrature->gaussPoints[iGP];
                GW = theGPQuadrature->weights[iGP];
                curveGPWeights[counterGP] = GW;

                // Compute the image of the GP in the curve parameter space
                GPUTilde = ((1.0 - GP)*UTildes[iSection] + (1.0 + GP)*UTildes[iSection + 1])/2.0;

                // Find the knot span in the parameter space of the curve
                knotSpanIndexCurve = dirichletCurve->getIGABasis()->findKnotSpan(GPUTilde);

                // Compute the basis functions at the parametric location of the curve
                dirichletCurve->getIGABasis()->computeLocalBasisFunctionsAndDerivatives
                        (localBasisFunctionsAndDerivativesCurve, noDeriv, GPUTilde, knotSpanIndexCurve);

                // Get the corresponding u,v parameters of the GPUTilde parametric location
                dirichletCurve->computeCartesianCoordinates
                        (uv, localBasisFunctionsAndDerivativesCurve, knotSpanIndexCurve);
                for (int iCoord = 0; iCoord < noCoordParam; iCoord++) {
                    _patch->getIGABasis(iCoord)->clampKnot(uv[iCoord]);
                    curveGPs[noCoordParam*counterGP + iCoord] = uv[iCoord];
                }

                // Compute the base vector of the curve
                dirichletCurve->computeBaseVectorAndDerivatives
                        (baseVectorCurve, knotSpanIndexCurve, localBasisFunctionsAndDerivativesCurve, noDerivBaseVct);

                // Find the knot span indices for on the patch
                uKnotSpan = _patch->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(uv[0]);
                vKnotSpan = _patch->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(uv[1]);

                // Compute the basis functions of the patch at the (u,v) parametric location
                _patch->getIGABasis()->computeLocalBasisFunctionsAndDerivatives
                        (localBasisFunctionsAndDerivatives, derivDegree, uv[0], uKnotSpan, uv[1], vKnotSpan);

                // Compute the base vectors on the patch in the physical space
                _patch->computeBaseVectors(baseVectors, localBasisFunctionsAndDerivatives, uKnotSpan, vKnotSpan);
                for(int iCoord = 0; iCoord < noCoord; iCoord++){
                    A1[iCoord] = baseVectors[iCoord];
                    A2[iCoord] = baseVectors[noCoord+iCoord];
                }

                // Compute the tangent vector on the physical space
                MathLibrary::computeDenseVectorMultiplicationScalar(A1, baseVectorCurve[0], noCoord);
                MathLibrary::computeDenseVectorMultiplicationScalar(A2, baseVectorCurve[1], noCoord);
                MathLibrary::computeDenseVectorAddition(A1, A2, 1.0, noCoord);

                // Compute the determinant of the Jacobian of the transformation from the parameter space to the physical space
                detJ2 = MathLibrary::vector2norm(A1, noCoord);
                MathLibrary::computeDenseVectorMultiplicationScalar(A1, 1.0/detJ2, noCoord);

                for(int iCoord = 0; iCoord < noCoord; iCoord++)
                    curveGPTangents[noCoord*counterGP + iCoord] = A1[iCoord];

                // Compute the element length on the Gauss point
                curveGPJacobianProducts[counterGP] = detJ1 * detJ2 * GW;

                // Update Gauss point counter
                counterGP++;
            }
        }

    // Set the initialized flag to true
    isGPDataInitialized = true;

    // Delete pointers
    delete theGPQuadrature;

}

void WeakIGADirichletCurveCondition::getCurveGPData(double* _curveGP, double& _curveGPWeight, double* _curveGPTangent, double& _curveGPJacobianProduct, int _iGP) {
    // Check if the GP data is initialized
    if (!isGPDataInitialized) assert(false);
    // Check if the given arrays are initialized
    if (_curveGP == NULL || _curveGPWeight == NULL || _curveGPTangent == NULL || _curveGPJacobianProduct == NULL) assert(false);
    // Check if the counter exceeds the number of GPs
    if (_iGP == curveNumGP) assert(false);

    // Initialize coordinates
    int noCoordParam = 2;
    int noCoord = 3;
    int iCoord;

    _curveGPWeight = curveGPWeights[_iGP];
    _curveGPJacobianProduct = curveGPJacobianProducts[_iGP];
    for (iCoord = 0; iCoord < noCoordParam; iCoord++)
        _curveGP[iCoord] = curveGPs[_iGP*2 + iCoord];
    for (iCoord = 0; iCoord < noCoord; iCoord++)
        _curveGPTangent[iCoord] = curveGPTangents[_iGP*3 + iCoord];

}

void WeakIGADirichletCurveCondition::getCurveAllGPData(double* _curveGPs, double* _curveGPWeights, double* _curveGPTangents, double* _curveGPJacobianProducts) {

    // Check if the GP data is initialized
    if (!isGPDataInitialized) assert(false);
    // Check if the given arrays are initialized
    if (_curveGPs == NULL || _curveGPWeights == NULL || _curveGPTangents == NULL || _curveGPJacobianProducts == NULL) assert(false);

    // Initialize coordinates
    const int noCoordParam = 2;
    const int noCoord = 3;
    int iCoord;

    // Copy the variables
    for (int iGP = 0; iGP < curveNumGP; iGP++) {
        _curveGPWeights[iGP] = curveGPWeights[iGP];
        _curveGPJacobianProducts[iGP] = curveGPJacobianProducts[iGP];
        for (iCoord = 0; iCoord < noCoordParam; iCoord++)
            _curveGPs[iGP*2 + iCoord] = curveGPs[iGP*2 + iCoord];
        for (iCoord = 0; iCoord < noCoord; iCoord++)
            _curveGPTangents[iGP*3 + iCoord] = curveGPTangents[iGP*3 + iCoord];
    }
}

WeakIGADirichletCurveCondition::~WeakIGADirichletCurveCondition() {
    if (!isTrimmingCurve)
        delete dirichletCurve;
    if (isGPDataInitialized) {
        delete curveGPs;
        delete curveGPWeights;
        delete curveGPTangents;
        delete curveGPJacobianProducts;
    }
}

} /* namespace EMPIRE */
