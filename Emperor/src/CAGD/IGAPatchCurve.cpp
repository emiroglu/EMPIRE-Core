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
#include "MathLibrary.h"
#include "MatrixVectorMath.h"
#include "Message.h"
#include <math.h>

using namespace std;

namespace EMPIRE {

int IGAPatchCurve::MAX_NUM_ITERATIONS = 50;
double IGAPatchCurve::TOL_CONVERGENCE = 1e-6;
int IGAPatchCurve::MAX_NUM_ITERATIONS_NEWTONRAPHSON = 20;
double IGAPatchCurve::TOL_CONVERGENCE_NEWTONRAPHSON = 1e-9;

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

void IGAPatchCurve::computeCartesianCoordinates(double* _cartesianCoordinates, int _knotSpanIndex, double* _locBasisFunctions) {
    /*
     * Returns the Cartesian coordinates of the parametric location where the basis functions are computed.
     */

    // Initialize auxiliary variables
    int indexCP;

    // Number of Cartesian coordinates
    int noCoord = 3;

    // Initialize output array
    for(int i = 0; i < noCoord; i++)
        _cartesianCoordinates[i] = 0.0;

    // Get the polynomial order of the basis
    int p = this->getIGABasis()->getPolynomialDegree();

    // Get the number of the basis functions at the parametric location
    int noLocBasisFunctions = p + 1;

    // Loop over all the local basis functions
    for (int iBF = 0; iBF < noLocBasisFunctions; iBF++) {
        // Get the index of the Control Point associated with the given basis function
        indexCP = _knotSpanIndex - p + iBF;

        // Compute the Cartesian coordinates of the point iteratively
        _cartesianCoordinates[0] += _locBasisFunctions[iBF]*ControlPointNet[indexCP].getX();
        _cartesianCoordinates[1] += _locBasisFunctions[iBF]*ControlPointNet[indexCP].getY();
        _cartesianCoordinates[2] += _locBasisFunctions[iBF]*ControlPointNet[indexCP].getZ();
    }
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

void IGAPatchCurve::computeBaseVectorAndDerivatives(double* _baseVectorAndDerivatives, int _knotSpanIndex, double* _basisFunctionsAndDerivatives, int _derivOrder) const{
    /*
     * Returns the base vector and its derivatives in an array of constant size [(polynomialOrder + 1)*(_derivOrder + 1)]. The coordinates of the n-th derivative
     * of the base vector can be retrieved as follows,
     *
     * G,n[0] = _baseVectorAndDerivatives[n*_derivOrder + 0]
     * G,n[1] = _baseVectorAndDerivatives[n*_derivOrder + 1]
     * G,n[2] = _baseVectorAndDerivatives[n*_derivOrder + 2]
     */

    // Initialize auxiliary variables
    int indexCP;

    // Number of Cartesian coordinates
    int noCoord = 3;

    // Initialize output array
    for(int i = 0; i < (_derivOrder + 1)*noCoord; i++)
        _baseVectorAndDerivatives[i] = 0.0;

    // Get the polynomial order of the basis
    int p = this->getIGABasis()->getPolynomialDegree();

    // Compute the number of the local basis functions
    int noLocBasisFunctions = p + 1;

    // Loop over all basis functions
    for(int iBF = 0; iBF < noLocBasisFunctions; iBF++){
        // Get the index of the Control Point associated with the current basis function
        indexCP = _knotSpanIndex - p + iBF;

        // Loop over all the derivative orders and add the corresponding contributions
        for(int iDeriv = 0; iDeriv < _derivOrder + 1; iDeriv++){
            _baseVectorAndDerivatives[noCoord*iDeriv + 0] += _basisFunctionsAndDerivatives[(iDeriv + 1)*noLocBasisFunctions + iBF]*ControlPointNet[indexCP].getX();
            _baseVectorAndDerivatives[noCoord*iDeriv + 1] += _basisFunctionsAndDerivatives[(iDeriv + 1)*noLocBasisFunctions + iBF]*ControlPointNet[indexCP].getY();
            _baseVectorAndDerivatives[noCoord*iDeriv + 2] += _basisFunctionsAndDerivatives[(iDeriv + 1)*noLocBasisFunctions + iBF]*ControlPointNet[indexCP].getZ();
        }
    }
}


bool IGAPatchCurve::computeIntersectionsWithKnotBisection(std::vector<double>& _uTilde, std::vector<double>& _uvSurface, unsigned int _dir, double _knot) const{

    /*
     * Found intersections are appended to the given _uvSurface and _uTilde vectors
     * Patch parameters of the intersections should be computed outside this function.
     * The bisection algorithm assumes one of the coordinates is constant and considers only one of the coordinates for convergence
     */

    // Initiate a single corresponding coordinate for both points
    double coordP1;
    double coordP2;

    // Check if the direction of the parameter curve is given correctly
    if (_dir<2)
        coordP1 = polyline.at(0+(_dir+1)%2);
    else
        ERROR_BLOCK_OUT("IGAPatchCurve","computeIntersectionWithKnotBisection","input variable \"dir\" can only be 0 for u or 1 for v direction!");

    double uvP[2];
    double uTilde = 0.0;
    bool isIntersecting;
    double uTildeP1 = 0.0;
    double uTildeP2 = 0.0;

    for (int iVertexCtr=1; iVertexCtr<polyline.size()/2; iVertexCtr++) {
        // Reset intersection flag at each iteration
        isIntersecting = false;

        // Get the corresponding coordinate of the P2
        coordP2 = polyline.at((iVertexCtr*2)+(_dir+1)%2);

        // Check if both points are aligned with a parameter line, if so skip this curve section
        if (fabs(coordP1-coordP2)<TOL_CONVERGENCE){
            coordP1 = coordP2;
            continue;
        }

        // Check if the knot is coinciding with P1 within given tolerance
        if (fabs(coordP1-_knot)<TOL_CONVERGENCE) {
            // Get the curve parameter and compute the coordinates
            uTilde = polylineKnots.at(iVertexCtr-1);
            computeCartesianCoordinates(uvP, uTilde);
            isIntersecting = true;
        }
        // Check if the knot is coinciding with P2 within given tolerance
        else if (fabs(coordP2-_knot)<TOL_CONVERGENCE) {
            // Get the corresponding curve parameter and compute the coordinates
            uTilde = polylineKnots.at(iVertexCtr);
            computeCartesianCoordinates(uvP, uTilde);
            isIntersecting = true;
        }
        // Check if the knot is between the two considered vertices on P1->P2 direction
        else if (coordP1>_knot && coordP2<_knot) {
            // Get the curve parameters of P1 and P2 in forward sense
            uTildeP1 = polylineKnots.at(iVertexCtr-1);
            uTildeP2 = polylineKnots.at(iVertexCtr);
            isIntersecting = solveIntersectionWithKnotBisection(uvP,uTilde,uTildeP1,uTildeP2,_dir,_knot);
        }
        // Check if the knot is between the two considered vertices on P2->P1 direction
        else if (coordP2>_knot && coordP1<_knot){
            // Get the curve parameters of P1 and P2 in backward sense
            uTildeP1 = polylineKnots.at(iVertexCtr);
            uTildeP2 = polylineKnots.at(iVertexCtr-1);
            isIntersecting = solveIntersectionWithKnotBisection(uvP,uTilde,uTildeP1,uTildeP2,_dir,_knot);
        }
        // If none of the cases holds then switch to the next vertex pair
        else {
            isIntersecting = false;
        }

        // Advance the first point
        coordP1 = coordP2;

        // If intersection found then add it to the list considering the ordering
        if (isIntersecting) {
            _uvSurface.push_back(uvP[0]);
            _uvSurface.push_back(uvP[1]);
            _uTilde.push_back(uTilde);
        }
    }

    if (_uvSurface.size()>0)    return true;
    else                        return false;
}

bool IGAPatchCurve::computeIntersectionsWithKnotBisection(std::vector<double>& _uTilde, unsigned int _dir, double _knot) const{

    /*
     * This function expects either an empty vector of curve parameters or a nonempty but sorted vector of curve parameters.
     * Found intersections are inserted into their corresponding positions in order to keep an ascending order.
     * Patch parameters of the intersections should be computed outside this function.
     * The bisection algorithm assumes one of the coordinates is constant and considers only the varying coordinate for convergence.
     */

    // Initiate a single corresponding coordinate for both points
    double coordP1;
    double coordP2;

    int noCoordParam = 2;

    // Check if the direction of the parameter curve is given correctly
    if (_dir<2)
        coordP1 = polyline.at(0+(_dir+1)%noCoordParam);
    else
        ERROR_BLOCK_OUT("IGAPatchCurve","computeIntersectionWithKnotBisection","input variable \"dir\" can only be 0 for u or 1 for v direction!");

    double uvP[2];
    double uTilde = 0.0;
    bool isIntersecting;
    double uTildeP1 = 0.0;
    double uTildeP2 = 0.0;

    for (int iVertexCtr=1; iVertexCtr<polyline.size()/noCoordParam; iVertexCtr++) {
        // Reset intersection flag at each iteration
        isIntersecting = false;

        // Get the corresponding coordinate of the P2
        coordP2 = polyline.at((iVertexCtr*noCoordParam)+(_dir+1)%noCoordParam);

        // Check if both points are aligned with a parameter line, if so skip this curve section
        if ( fabs(coordP1 - coordP2) < TOL_CONVERGENCE ) {
            coordP1 = coordP2;
            continue;
        }

        // Check if the knot is coinciding with P1 within given tolerance
        if ( fabs(coordP1 - _knot) < TOL_CONVERGENCE ) {
            // Get the curve parameter and compute the coordinates
            uTilde = polylineKnots.at(iVertexCtr-1);
            computeCartesianCoordinates(uvP, uTilde);
            isIntersecting = true;
        }
        // Check if the knot is coinciding with P2 within given tolerance
        else if ( fabs(coordP2 - _knot) < TOL_CONVERGENCE ) {
            // Get the corresponding curve parameter and compute the coordinates
            uTilde = polylineKnots.at(iVertexCtr);
            computeCartesianCoordinates(uvP, uTilde);
            isIntersecting = true;
        }
        // Check if the knot is between the two considered vertices on P1->P2 direction
        else if (coordP1 > _knot && coordP2 < _knot) {
            // Get the curve parameters of P1 and P2 in forward sense
            uTildeP1 = polylineKnots.at(iVertexCtr-1);
            uTildeP2 = polylineKnots.at(iVertexCtr);
            isIntersecting = solveIntersectionWithKnotBisection(uvP, uTilde, uTildeP1, uTildeP2, _dir, _knot);
        }
        // Check if the knot is between the two considered vertices on P2->P1 direction
        else if (coordP2 > _knot && coordP1 < _knot){
            // Get the curve parameters of P1 and P2 in backward sense
            uTildeP1 = polylineKnots.at(iVertexCtr);
            uTildeP2 = polylineKnots.at(iVertexCtr-1);
            isIntersecting = solveIntersectionWithKnotBisection(uvP, uTilde, uTildeP1, uTildeP2, _dir, _knot);
        }
        // If none of the cases holds then switch to the next vertex pair
        else {
            isIntersecting = false;
        }

        // Advance the first point
        coordP1 = coordP2;

        std::vector<double>::iterator iUTilde = _uTilde.begin();
        // If intersection found then add it to the list considering the ordering in the given _uTilde variable
        if (isIntersecting) {

            iUTilde = _uTilde.begin();

            if (_uTilde.empty())    _uTilde.push_back(uTilde);
            else {
                while( (uTilde > *iUTilde) && (iUTilde != _uTilde.end()) )  iUTilde++;
                if ( uTilde != *iUTilde )   _uTilde.insert(iUTilde, uTilde);
                else if ( iUTilde == _uTilde.end())   _uTilde.push_back(uTilde);
            }
        }
    }

    if (_uTilde.size()>0)    return true;
    else                     return false;
}

bool IGAPatchCurve::solveIntersectionWithKnotBisection(double* _uvP, double& _uTildeP, double _uTildeP1, double _uTildeP2, unsigned int _dir, double _knot) const{

    /*
     * The bisection algorithm assumes one of the coordinates is constant and considers only the varying coordinate for convergence.
     */

    int noCoordParam = 2;
    double distance = 1.0;
    unsigned int iteration = 0;

    while (true){
        // Bisection point coordinate in curve parameter space
        _uTildeP = 0.5*(_uTildeP1 + _uTildeP2);

        // Bisection point coordinate in patch parameter space
        computeCartesianCoordinates(_uvP, _uTildeP);

        // Distance of the bisection point to the knot (checks the distance in v when dir is u(0) and vice versa)
        distance = _uvP[(_dir+1)%noCoordParam]-_knot;

        // If the distance is within a tolerance the intersection is found at the P
        if (fabs(distance) < TOL_CONVERGENCE){
            return true;
        } // If the distance is greater than 0 then the intersection is between P and P1
        if (distance > 0.0){
            _uTildeP1 = _uTildeP;
        } // If the distance is less than 0 then the intersection is between P and P2
        else if (distance < 0.0) {
            _uTildeP2 = _uTildeP;
        }

        iteration++;
        if (iteration > MAX_NUM_ITERATIONS)     return false;
    }
}

bool IGAPatchCurve::computePointProjectionOn2DCurve(double& _uPrm, double* _P, int _noCoord, int _maxIt, double _tol) {
    /*
     * Returns the parametric and Cartesian coordinates of the projected point onto the 2D curve.
     */

    // Initialize auxiliary variables
    int noCoord = 3;
    int knotSpanIndex;
    double normEuclideanDistance;
    double jacobian;
    double residual;
    double normBaseVectorSquare;
    double distanceTimesCurvatureVector;
    double P[3];
    double cartesianCoordinates[3];
    double euclideanDistance[3];
    double baseVector[3];
    double curvatureVector[3];

    // Check input
    if(_noCoord == 2){
        P[0] = _P[0];
        P[1] = _P[1];
        P[3] = 0.0;
    }else if(_noCoord == 3){
        for(int iCoord = 0; iCoord < noCoord; iCoord++)
            P[iCoord] = _P[iCoord];
    }else{
        ERROR_OUT() << "The number of coordinates of the input point must be 2 or 3" << std::endl;
        exit(EXIT_FAILURE);
    }

    // Initialize Newton-Raphson iteration counter
    int counterNR = 1;

    // Get the polynomial order of the curve
    int p = this->getIGABasis()->getPolynomialDegree();

    // Compute the number of local basis functions
    int noLocBasisFunctions = p + 1;

    // Number of derivatives for the basis functions
    int noDeriv = 2;

    // Initialize the basis functions and their derivatives as well as the base vector and its derivative
    double* locBasisFunctionsAndDerivatives = new double[noLocBasisFunctions*(noDeriv + 1)];
    double* baseVectorAndDerivatives = new double[noDeriv*noCoord];

    // Initialize convergence flag
    bool isConvergent = false;

    // Loop over all Newton-Rapshon iterations
    while(!isConvergent && counterNR <= _maxIt){
        // Find the knot span of the parameter on the curve
        knotSpanIndex = this->findKnotSpan(_uPrm);

        // Compute the basis functions and up to their second derivatives
        this->getIGABasis()->computeLocalBasisFunctionsAndDerivatives(locBasisFunctionsAndDerivatives, noDeriv, _uPrm, knotSpanIndex);

        // Compute the base vector and its derivative
        this->computeBaseVectorAndDerivatives(baseVectorAndDerivatives, knotSpanIndex, locBasisFunctionsAndDerivatives, noDeriv - 1);

        // Compute the Cartesian coordinates of the parametric coordinate
        this->computeCartesianCoordinates(cartesianCoordinates, knotSpanIndex, locBasisFunctionsAndDerivatives);

        // Compute the Euclidean distance of the point to be projected and the projected point
        for(int iCoord = 0; iCoord < noCoord; iCoord++)
            euclideanDistance[iCoord] = P[iCoord] - cartesianCoordinates[iCoord];

        // Compute the norm of the Euclidean distance between the two points
        normEuclideanDistance = MathLibrary::vector2norm(euclideanDistance, noCoord);

        // Check for point coincidence
        if(normEuclideanDistance < _tol){
            isConvergent = true;
            break;
        }

        // Compute the square of the norm of the base vector and the product of the curvature vector with the distance vector
        for(int iCoord = 0; iCoord < noCoord; iCoord++){
            baseVector[iCoord] = baseVectorAndDerivatives[0*noCoord + iCoord];
            curvatureVector[iCoord] = baseVectorAndDerivatives[1*noCoord + iCoord];
        }
        normBaseVectorSquare = MathLibrary::computeDenseDotProduct(noCoord, baseVector, baseVector);
        distanceTimesCurvatureVector = MathLibrary::computeDenseDotProduct(noCoord, euclideanDistance, curvatureVector);

        // Compute the Jacobian of the system
        jacobian = - normBaseVectorSquare + distanceTimesCurvatureVector;

        // Compute the residual of the system
        residual = MathLibrary::computeDenseDotProduct(noCoord, euclideanDistance, baseVector);

        // Check the orthogonality condition
        if(fabs(residual) < _tol){
            isConvergent = true;
            break;
        }

        // Update the parametric coordinate
        _uPrm -= residual/jacobian;

        // Clamp parameter in case it exceeds the limits of the knot vector
        if(_uPrm < this->getIGABasis()->getFirstKnot())
            _uPrm = this->getIGABasis()->getFirstKnot();
        if(_uPrm > this->getIGABasis()->getLastKnot())
            _uPrm = this->getIGABasis()->getLastKnot();

        // Update Newton-Raphson iteration counter
        counterNR++;
    }

    // Copy the Cartesian and the parametric coordinates of the projected point
    for(int iCoord = 0; iCoord < _noCoord; iCoord++)
        _P[iCoord] = cartesianCoordinates[iCoord];

    // Delete pointers
    delete[] locBasisFunctionsAndDerivatives;
    delete[] baseVectorAndDerivatives;

    // Return the convergence flag
    return isConvergent;
}

} /* namespace EMPIRE */
