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
#include "IGAPatchSurface.h"
#include "WeakIGADirichletBoundaryCondition.h"
#include "MathLibrary.h"
#include "Message.h"

using namespace std;

namespace EMPIRE {

WeakIGADirichletBoundaryCondition::WeakIGADirichletBoundaryCondition(int _ID,
                                 int _patchIndex, int _patchBLIndex, int _patchBLTrCurveIndex) :
        AbstractCondition(_ID),
        patchIndex(_patchIndex), patchBLIndex(_patchBLIndex), patchBLTrCurveIndex(_patchBLTrCurveIndex)
        {

    type = EMPIRE_WeakIGADirichletBoundaryCondition;

    isGPDataInitialized = false;

}

void WeakIGADirichletBoundaryCondition::addWeakDirichletBoundaryConditionGPData(int _trCurveNumGP,
                            double* _trCurveGPs, double* _trCurveGPWeights,
                            double* _trCurveGPTangents,
                            double* _trCurveGPJacobianProducts){

    if (isGPDataInitialized) assert(false);

    const int noCoordParam = 2;
    const int noCoord = 3;

    trCurveNumGP = _trCurveNumGP;

    trCurveGPs = new double[trCurveNumGP*noCoordParam];
    trCurveGPWeights = new double[trCurveNumGP];
    trCurveGPTangents = new double[trCurveNumGP*noCoord];
    trCurveGPJacobianProducts = new double[trCurveNumGP];

    for (int i = 0; i < trCurveNumGP; i++){
        for (int j = 0; j < noCoordParam; j++){
            trCurveGPs[i*noCoordParam+j]=_trCurveGPs[i*noCoordParam+j];
        }
        trCurveGPWeights[i] = _trCurveGPWeights[i];
        for (int j = 0; j < noCoord; j++){
            trCurveGPTangents[i*noCoord+j]=_trCurveGPTangents[i*noCoord+j];
        }
        trCurveGPJacobianProducts[i] = _trCurveGPJacobianProducts[i];
    }

    isGPDataInitialized = true;
}

void WeakIGADirichletBoundaryCondition::createGPData(IGAPatchSurface* _patch) {

    if (isGPDataInitialized) assert(false);

    // Initialize the coordinates
    const int noCoordParam = 2;
    const int noCoord = 3;

    // Initialize parameter and coordinate sets to store the intersections and transfer them
    std::vector<double> UTildes;

    // Compute the knot intersections of the trimming curve
    _patch->computeKnotIntersectionsWithTrimmingCurve(UTildes,
                                                      patchBLIndex, patchBLTrCurveIndex);

    // Sort and remove duplicates if they exist
    EMPIRE::MathLibrary::sortRemoveDuplicates(UTildes);

    /// Create the Gauss points on the master and the slave sides
    // Getting the polynomial orders
    int p = _patch->getIGABasis(0)->getPolynomialDegree();
    int q = _patch->getIGABasis(1)->getPolynomialDegree();
    int pMax = std::max(p, q);

    // Get the number of Gauss points
    int numGPPerSection = pMax+1;
    trCurveNumGP = numGPPerSection * (UTildes.size()-1);

    // Initialize the GP data
    trCurveGPWeights = new double[trCurveNumGP];
    trCurveGPJacobianProducts = new double[trCurveNumGP];
    trCurveGPs = new double[trCurveNumGP*noCoordParam];
    trCurveGPTangents = new double[trCurveNumGP*noCoord];

    // Create a Gauss quadrature rule
    MathLibrary::IGAGaussQuadratureOnBiunitInterval* theGPQuadrature = new MathLibrary::IGAGaussQuadratureOnBiunitInterval(numGPPerSection);

    // Initialize variables
    int noDeriv = 1;
    int noDerivBaseVct = 0;
    int derivDegree = 1;
    int uKnotSpan;
    int vKnotSpan;
    int knotSpanIndexTrCurve;
    int pTrCurve = _patch->getTrimming().getLoop(patchBLIndex).getIGACurve(patchBLTrCurveIndex).getIGABasis()->getPolynomialDegree();
    int noLocalBasisFunctionsTrCurve = pTrCurve + 1;
    int noLocalBasisFunctions = (p + 1)*(q + 1);
    double GP;
    double GW;
    double GPUTilde;
    double detJ1;
    double detJ2;

    // Initialize pointers
    double uv[2];
    double baseVectorTrCurve[(noDerivBaseVct + 1)*noCoord];
    double localBasisFunctionsAndDerivatives[(derivDegree + 1) * (derivDegree + 2) * noLocalBasisFunctions / 2];
    double localBasisFunctionsAndDerivativesTrCurve[noLocalBasisFunctionsTrCurve * (noDeriv + 1)];
    double baseVectors[6];
    double A1[3];
    double A2[3];

    // Initialize GP counter on the curve
    int counterGP = 0;

    // Loop over the sections of the trimming curve
        for (int iSection = 0; iSection < UTildes.size() - 1; iSection++) {

            // Determinant of the Jacobian of the transformation from the parent space to the parameter space of the patch
            detJ1 = (UTildes[iSection+1]-UTildes[iSection])/2.0;

            // Loop over the GPs of the section
            for (int iGP = 0; iGP < numGPPerSection; iGP++) {

                // Get GP coordinates and weights
                GP = theGPQuadrature->gaussPoints[iGP];
                GW = theGPQuadrature->weights[iGP];
                trCurveGPWeights[counterGP] = GW;

                // Compute the image of the GP in the curve parameter space
                GPUTilde = ((1.0 - GP)*UTildes[iSection] + (1.0 + GP)*UTildes[iSection + 1])/2.0;

                // Find the knot span in the parameter space of the curve
                knotSpanIndexTrCurve = _patch->getTrimming().getLoop(patchBLIndex).getIGACurve(patchBLTrCurveIndex).getIGABasis()->findKnotSpan(GPUTilde);

                // Compute the basis functions at the parametric location of the trimming curve
                _patch->getTrimming().getLoop(patchBLIndex).getIGACurve(patchBLTrCurveIndex).getIGABasis()->computeLocalBasisFunctionsAndDerivatives
                        (localBasisFunctionsAndDerivativesTrCurve, noDeriv, GPUTilde, knotSpanIndexTrCurve);

                // Get the corresponding u,v parameters of the GPUTilde parametric location
                _patch->getTrimming().getLoop(patchBLIndex).getIGACurve(patchBLTrCurveIndex).computeCartesianCoordinates
                        (uv, localBasisFunctionsAndDerivativesTrCurve, knotSpanIndexTrCurve);
                trCurveGPs[noCoordParam*counterGP] = uv[0];
                trCurveGPs[noCoordParam*counterGP + 1] = uv[1];

                // Compute the base vector of the trimming curve
                _patch->getTrimming().getLoop(patchBLIndex).getIGACurve(patchBLTrCurveIndex).computeBaseVectorAndDerivatives
                        (baseVectorTrCurve, knotSpanIndexTrCurve, localBasisFunctionsAndDerivativesTrCurve, noDerivBaseVct);

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
                MathLibrary::computeDenseVectorMultiplicationScalar(A1, baseVectorTrCurve[0], noCoord);
                MathLibrary::computeDenseVectorMultiplicationScalar(A2, baseVectorTrCurve[1], noCoord);
                MathLibrary::computeDenseVectorAddition(A1, A2, 1.0, noCoord);

                // Compute the determinant of the Jacobian of the transformation from the parameter space to the physical space
                detJ2 = MathLibrary::vector2norm(A1, noCoord);
                MathLibrary::computeDenseVectorMultiplicationScalar(A1, 1.0/detJ2, noCoord);

                for(int iCoord = 0; iCoord < noCoord; iCoord++)
                    trCurveGPTangents[noCoord*counterGP + iCoord] = A1[iCoord];

                // Compute the element length on the Gauss point
                trCurveGPJacobianProducts[counterGP] = detJ1 * detJ2 * GW;

                // Update Gauss point counter
                counterGP++;
            }
        }

    // Set the initialized flag to true
    isGPDataInitialized = true;

    // Delete pointers
    delete theGPQuadrature;

}

WeakIGADirichletBoundaryCondition::~WeakIGADirichletBoundaryCondition() {
    if (isGPDataInitialized) {
        delete trCurveGPs;
        delete trCurveGPWeights;
        delete trCurveGPTangents;
        delete trCurveGPJacobianProducts;
    }
}


} /* namespace EMPIRE */
