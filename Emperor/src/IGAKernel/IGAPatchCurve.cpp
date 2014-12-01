/*  Copyright &copy; 2013, TU Muenchen, Chair of Structural Analysis,
 *  Fabien Pean, Andreas Apostolatos, Munich
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

#include <assert.h>
#include "IGAPatchCurve.h"

using namespace std;

namespace EMPIRE {

IGAPatchCurve::IGAPatchCurve(int _IDBasis, int _pDegree, int _uNoKnots, double* _uKnotVector,
		int _uNoControlPoints, double* _controlPointNet):
		uNoControlPoints(_uNoControlPoints) {
	// Read input
	bool ucondition = _uNoControlPoints != _uNoKnots - _pDegree - 1;
	if (ucondition) {
		ERROR_BLOCK_OUT("IGAPatchCurve","IGAPatchCurve","Number of Control Points, number of knots and polynomial degree do not match!");
	}
	// On the Control Point net
	assert(_controlPointNet!=NULL);
	ControlPointNet.reserve(uNoControlPoints);
    for (int i = 0; i < uNoControlPoints; i++) {
            ControlPointNet.push_back(IGAControlPoint(i, &_controlPointNet[i * 4]));
    }
	// Figure out whether the patch has a B-Spline or a NURBS underlying basis
	int isNurbs = 0;
	for (int i = 0; i < uNoControlPoints; i++) {
		if (ControlPointNet[i].getW() != 1.0) {
			isNurbs = 1;
			break;
		}
	}
	// Create the NURBS or the B-Spline underlying basis
	if (!isNurbs) {
		IGABasis = new BSplineBasis1D(_IDBasis, _pDegree, _uNoKnots, _uKnotVector);
	} else {
		double* controlPointWeights = new double[uNoControlPoints];
		for (int i = 0; i < uNoControlPoints; i++)
			controlPointWeights[i] = ControlPointNet[i].getW();
		IGABasis = new NurbsBasis1D(_IDBasis, _pDegree, _uNoKnots, _uKnotVector, _uNoControlPoints, controlPointWeights);
	}
}

IGAPatchCurve::~IGAPatchCurve() {
	delete IGABasis;
}

void IGAPatchCurve::computeCartesianCoordinates(double* _cartesianCoordinates, double _uPrm,
        int _uKnotSpanIndex) const {
    // Read input
    assert(_cartesianCoordinates != NULL);

    // Initialize the coordinates of the point
    for (int i = 0; i < 2; i++)
        _cartesianCoordinates[i] = 0;
    // Compute the local basis functions in the vicinity of the point
    int pDegree = IGABasis->getPolynomialDegree();
    int noLocalBasisFunctions = IGABasis->computeNoBasisFunctions();
    double localBasisFunctions[noLocalBasisFunctions];
    IGABasis->computeLocalBasisFunctions(localBasisFunctions,_uPrm,_uKnotSpanIndex);
    // Initialize the Control Point index
    int CPindex = 0;
    // Initialize a basis functions counter
    int counter_basis = 0;
    // Loop over all the non-zero contributions
	for (int i = 0; i <= pDegree; i++) {

		// Update the correct index for the Control Points in 2D. Pattern A[i][j] = V[j*n+i]
		CPindex =(_uKnotSpanIndex - pDegree + i);

		// Compute iteratively the x-coordinate of the point
		_cartesianCoordinates[0] += localBasisFunctions[counter_basis]
				* ControlPointNet[CPindex].getX();
		// Compute iteratively the y-coordinate of the point
		_cartesianCoordinates[1] += localBasisFunctions[counter_basis]
				* ControlPointNet[CPindex].getY();

		// Update basis function's counter
		counter_basis++;
	}
}

void IGAPatchCurve::computeCartesianCoordinates(double* _cartesianCoordinates, double _localCoordinates) const {
    int _uKnotSpanIndex = IGABasis->findKnotSpan(_localCoordinates);
    computeCartesianCoordinates(_cartesianCoordinates, _localCoordinates, _uKnotSpanIndex);
}

} /* namespace EMPIRE */
