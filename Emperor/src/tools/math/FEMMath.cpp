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

#include "FEMMath.h"
#include "ConstantsAndVariables.h"
#include "GeometryMath.h"
#include "MatrixVectorMath.h"
#include "DebugMath.h"
#include "Message.h"

using namespace std;
namespace EMPIRE {
namespace MathLibrary {



// Methods

/***********************************************************************************************
 * \brief Compute mass matrix of a triangle element
 * \param[in] triangle the triangle
 * \param[in] numGaussPoints number of Gauss points used in the Gauss quadrature
 * \param[in] dual whether dual or not
 * \param[out] mass matrix (3x3)
 * \author Tianyang Wang
 ***********/
void computeMassMatrixOfTrianlge(const double *triangle, int numGaussPoints, bool dual,
        double *massMatrix) {
    const double *gaussPointsLocal;
    const double *weights;
    if (numGaussPoints == 3) {
        gaussPointsLocal = EMPIRE::MathLibrary::triGaussPoints3;
        weights = EMPIRE::MathLibrary::triWeights3;
    } else if (numGaussPoints == 6) {
        gaussPointsLocal = EMPIRE::MathLibrary::triGaussPoints6;
        weights = EMPIRE::MathLibrary::triWeights6;
    } else if (numGaussPoints == 7) {
        gaussPointsLocal = EMPIRE::MathLibrary::triGaussPoints7;
        weights = EMPIRE::MathLibrary::triWeights7;
    } else if (numGaussPoints == 12) {
        gaussPointsLocal = EMPIRE::MathLibrary::triGaussPoints12;
        weights = EMPIRE::MathLibrary::triWeights12;
    } else {
        assert(false);
    }
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            massMatrix[i * 3 + j] = 0.0;
        }
    }
    double area = EMPIRE::MathLibrary::computeAreaOfTriangle(triangle);
    if (!dual) {
        // since the mass matrix is symmetric, only calculate the upper part
        for (int i = 0; i < 3; i++) {
            for (int j = i; j < 3; j++) {
                for (int k = 0; k < numGaussPoints; k++) {
                    massMatrix[i * 3 + j] += weights[k] * gaussPointsLocal[k * 3 + i]
                            * gaussPointsLocal[k * 3 + j];
                }
                massMatrix[i * 3 + j] *= area;
            }
        }
        // set the lower part
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < i; j++) {
                massMatrix[i * 3 + j] = massMatrix[j * 3 + i];
            }
        }
    } else {
        for (int i = 0; i < 3; i++) {
            for (int k = 0; k < numGaussPoints; k++) {
                massMatrix[i * 3 + i] += weights[k] * gaussPointsLocal[k * 3 + i];
            }
            massMatrix[i * 3 + i] *= area;
        }
    }
}

/***********************************************************************************************
 * \brief Compute mass matrix of a quad element
 * \param[in] quad the quad
 * \param[in] numGaussPoints number of Gauss points used in the Gauss quadrature
 * \param[in] dual whether dual or not
 * \param[out] mass matrix (4x4)
 * \author Tianyang Wang
 ***********/
void computeMassMatrixOfQuad(const double *quad, int numGaussPoints, bool dual,
        double *massMatrix) {
    const double *gaussPointsLocal;
    const double *weights;
    if (numGaussPoints == 1) {
        gaussPointsLocal = EMPIRE::MathLibrary::quadGaussPoints1;
        weights = EMPIRE::MathLibrary::quadWeights1;
    } else if (numGaussPoints == 4) {
        gaussPointsLocal = EMPIRE::MathLibrary::quadGaussPoints4;
        weights = EMPIRE::MathLibrary::quadWeights4;
    } else if (numGaussPoints == 9) {
        gaussPointsLocal = EMPIRE::MathLibrary::quadGaussPoints9;
        weights = EMPIRE::MathLibrary::quadWeights9;
    } else {
        assert(false);
    }
    double GPShapeFunc[numGaussPoints * 4];
    for (int i = 0; i < numGaussPoints; i++) {
        for (int j = 0; j < 4; j++) {
            computeShapeFuncOfQuad(&gaussPointsLocal[i * 2], &GPShapeFunc[i * 4]);
        }
    }
    double detJ[numGaussPoints];
    for (int i = 0; i < numGaussPoints; i++) {
        detJ[i] = computeDetJOfQuad(quad, &gaussPointsLocal[i * 2]);
    }
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            massMatrix[i * 4 + j] = 0.0;
        }
    }
    if (!dual) {
        // since the mass matrix is symmetric, only calculate the upper part
        for (int i = 0; i < 4; i++) {
            for (int j = i; j < 4; j++) {
                for (int k = 0; k < numGaussPoints; k++) {
                    massMatrix[i * 4 + j] += weights[k] * detJ[k] * GPShapeFunc[k * 4 + i]
                            * GPShapeFunc[k * 4 + j];
                }
            }
        }
        // set the lower part
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < i; j++) {
                massMatrix[i * 4 + j] = massMatrix[j * 4 + i];
            }
        }
    } else {
        for (int i = 0; i < 4; i++) {
            for (int k = 0; k < numGaussPoints; k++) {
                massMatrix[i * 4 + i] += weights[k] * detJ[k] * GPShapeFunc[k * 4 + i];
            }
        }
    }
}

/***********************************************************************************************
 * \brief Compute the shape function value by local coordinates in a quadrilateral
 * \param[in] xi_eta local coordinates
 * \param[out] shapeFuncValues shape function values of xi_eta
 * \author Tianyang Wang
 ***********/
void computeShapeFuncOfQuad(const double *xi_eta, double *shapeFuncValues) {
    shapeFuncValues[0] = 0.25 * (1.0 - xi_eta[0]) * (1.0 - xi_eta[1]);
    shapeFuncValues[1] = 0.25 * (1.0 + xi_eta[0]) * (1.0 - xi_eta[1]);
    shapeFuncValues[2] = 0.25 * (1.0 + xi_eta[0]) * (1.0 + xi_eta[1]);
    shapeFuncValues[3] = 0.25 * (1.0 - xi_eta[0]) * (1.0 + xi_eta[1]);
}

/***********************************************************************************************
 * \brief Compute the determinant of Jocobian matrix by local coordinates in a quadrilateral
 * \param[in] quad the quadrilateral
 * \param[in] xi_eta local coordinates
 * \return determinant of Jocobian matrix
 * \author Tianyang Wang
 ***********/
double computeDetJOfQuad(const double *quad, const double *xi_eta) {
    // d_N_d_xi[4] contains the partial derivative w.r.t. xi of the shape functions
    double d_N_d_xi[4];
    // d_N_d_eta[4] contains the partial derivative w.r.t. eta of the shape functions
    double d_N_d_eta[4];

    d_N_d_xi[0] = -0.25 * (1 - xi_eta[1]);
    d_N_d_xi[1] = -d_N_d_xi[0];
    d_N_d_xi[2] = 0.25 * (1 + xi_eta[1]);
    d_N_d_xi[3] = -d_N_d_xi[2];

    d_N_d_eta[0] = -0.25 * (1 - xi_eta[0]);
    d_N_d_eta[1] = -0.25 * (1 + xi_eta[0]);
    d_N_d_eta[2] = -d_N_d_eta[1];
    d_N_d_eta[3] = -d_N_d_eta[0];

    // g1 and g2 are the local basis vectors, and det(J)=||g1 x g2||
    double g1[3] = { 0, 0, 0 };
    double g2[3] = { 0, 0, 0 };

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 4; j++) {
            g1[i] += quad[j * 3 + i] * d_N_d_xi[j];
            g2[i] += quad[j * 3 + i] * d_N_d_eta[j];
        }
    }

    double crossProduct[3];
    EMPIRE::MathLibrary::computeVectorCrossProduct(g1, g2, crossProduct);
    return EMPIRE::MathLibrary::computeVectorLength(crossProduct);
}

/***********************************************************************************************
 * \brief Compute global coordinates of a point in a triangle
 * \param[in] triangle the triangle
 * \param[in] localCoor local coordinates of the point
 * \param[out] globalCoor global coordinates of the point
 * \author Tianyang Wang
 ***********/
void computeGlobalCoorInTriangle(const double *triangle, const double *localCoor,
        double *globalCoor) {
    for (int i = 0; i < 3; i++) {
        globalCoor[i] = 0;
        for (int j = 0; j < 3; j++)
            globalCoor[i] += localCoor[j] * triangle[i + j * 3];
    }
}

/***********************************************************************************************
 * \brief Compute global coordinates of a point in a quad
 * \param[in] quad the quad
 * \param[in] localCoor local coordinates of the point
 * \param[out] globalCoor global coordinates of the point
 * \author Tianyang Wang
 ***********/
void computeGlobalCoorInQuad(const double *quad, const double *localCoor, double *globalCoor) {
    double shapeFuncValues[4];
    computeShapeFuncOfQuad(localCoor, shapeFuncValues);
    for (int i = 0; i < 3; i++)
        globalCoor[i] = 0.0;

    for (int i = 0; i < 3; i++) { // x, y, z of globalCoor
        for (int j = 0; j < 4; j++) { // 4 end nodes of quad
            globalCoor[i] += shapeFuncValues[j] * quad[j * 3 + i];
        }
    }
}

/***********************************************************************************************
 * \brief Compute local coordinates of a point in a triangle
 * \param[in] triangle the triangle
 * \param[in] planeToProject project to plane (case: {2:x-y ; 0:y-z ;1: z-x} )
 * \param[in] point the point
 * \param[out] localCoor local coordinates of the point
 * \return a boolean saying whether the point is inside the triangle or not (true of false)
 * \author Tianyang Wang
 ***********/
bool computeLocalCoorInTriangle(const double *triangle, int planeToProject, const double *point,
        double *localCoor) {
    /* Normally, the local coordinates can be solved by:
     *  |x_1 x_2 x_3|   |xi_1|    |x_4|
     *  |y_1 y_2 y_3| * |xi_2| =  |y_4|
     *  |z_1 z_2 z_3|   |xi_3|    |z_4|
     *
     * But if x_1 = x_2 = x_3 = 0 than the system cannot be solved.
     * So we remove one equation by xi_1 + xi_2 + xi_3 = 1.
     * This indicates projection.
     * Choose among planes (xy or yz or zx) the one which has the smallest angle
     * with the triangle normal.
     * For example, if the angle between of normals of yz-plane and the triangle
     * is 0, then the system is replaced by
     *  |1.0 1.0 1.0|   |xi_1|    |1.0|
     *  |y_1 y_2 y_3| * |xi_2| =  |y_4|
     *  |z_1 z_2 z_3|   |xi_3|    |z_4|
     */
    double A[9]; // Attention! A is the transpose of the above matrix
    for (int i = 0; i < 9; i++)
        A[i] = triangle[i];

    for (int i = 0; i < 3; i++)
        localCoor[i] = point[i];

    A[planeToProject] = A[planeToProject + 3] = A[planeToProject + 6] = localCoor[planeToProject] =
            1.0;

    EMPIRE::MathLibrary::solve3x3LinearSystem(A, planeToProject, localCoor);

    // make sure the sum is 1.0
    assert(fabs(localCoor[0] + localCoor[1] + localCoor[2] -1.0) < 1E-12);
    localCoor[0] = 1.0 - localCoor[1] - localCoor[2];

    for (int i = 0; i < 3; i++) {
        if (localCoor[i] > 1.0)
        	return false;
        if (localCoor[i] < 0.0)
            return false;
    }
    return true;
}


/***********************************************************************************************
 * \brief Compute local coordinates of a point in a quadrilateral
 * \param[in] quad the quadrilateral
 * \param[in] planeToProject project to plane (case: {2:x-y ; 0:y-z ;1: z-x} )
 * \param[in] point the point
 * \param[out] localCoor local coordinates of the point
 * \return a boolean saying whether the point is inside the quadrilateral or not (true of false)
 * \author Tianyang Wang
 ***********/
bool computeLocalCoorInQuad(const double *quad, int planeToProject, const double *point,
        double *localCoor) {
    /*
     * So we use two coordinates among x, y, z.
     * This indicates projection.
     * Choose among planes (xy or yz or zx) the one which has the smallest angle
     * with the quad normal.
     */
    double x[4];
    double y[4];
    int x_direc = (planeToProject + 1) % 3;
    int y_direc = (planeToProject + 2) % 3;
    for (int i = 0; i < 4; i++) {
        x[i] = quad[i * 3 + x_direc];
        y[i] = quad[i * 3 + y_direc];
    }
    double x0 = point[x_direc];
    double y0 = point[y_direc];

    double a1 = x[0] + x[1] + x[2] + x[3] - 4.0 * x0;
    double b1 = -x[0] + x[1] + x[2] - x[3];
    double c1 = -x[0] - x[1] + x[2] + x[3];
    double d1 = x[0] - x[1] + x[2] - x[3];

    double a2 = y[0] + y[1] + y[2] + y[3] - 4.0 * y0;
    double b2 = -y[0] + y[1] + y[2] - y[3];
    double c2 = -y[0] - y[1] + y[2] + y[3];
    double d2 = y[0] - y[1] + y[2] - y[3];

    double delta[2];
    double J_T[4]; // transpose of Jacobian --- to use column major in lapack
    double F[2]; // -F
    localCoor[0] = 0;
    localCoor[1] = 0;
    //int dummy[2];
    const double EPS = 1E-13; // should be smaller than 1E-15 from experience

    const int MAX_ITER_NUM = 100;
    for (int i = 0; i < MAX_ITER_NUM; i++) {
        J_T[0] = b1 + d1 * localCoor[1];
        J_T[2] = c1 + d1 * localCoor[0];
        J_T[1] = b2 + d2 * localCoor[1];
        J_T[3] = c2 + d2 * localCoor[0];
        F[0] = a1 + b1 * localCoor[0] + c1 * localCoor[1] + d1 * localCoor[0] * localCoor[1];
        F[1] = a2 + b2 * localCoor[0] + c2 * localCoor[1] + d2 * localCoor[0] * localCoor[1];
        delta[0] = -F[0];
        delta[1] = -F[1];
        /*int info = LAPACKE_dgesv(LAPACK_COL_MAJOR, 2, 1, J_T, 2, dummy,
         delta, 2);
         if (info != 0) {
         cerr << "ERROR in MortarMath::computeLocalCoorInQuad()!" << endl;
         exit(EXIT_FAILURE);
         }*/
        solve2x2LinearSystem(J_T, delta, EPS);
        if (i >= 10) { // normally should find a solution within 10 iterations
        	WARNING_BLOCK_OUT("FEMMath","computeLocalCoordInQuad", "More than 10 iterations are necessary for computing local coordinates in quad");
            WARNING_OUT() << "iteration #: " << i << endl;
            EMPIRE::MathLibrary::printElem(quad, 4);
            EMPIRE::MathLibrary::printPoint(point);
            DEBUG_OUT() << "plane to project:  " << planeToProject << endl;
            DEBUG_OUT() << "xi:  " << localCoor[0] << " ita: " << localCoor[1] << endl;
            DEBUG_OUT() << "delta-xi:  " << delta[0] << " delta-ita: " << delta[1] << endl;
            // if point is far out of quad, return false
            for (int i = 0; i < 2; i++) { // do not care accuracy if point is far outside the quad
                if (localCoor[i] > 2.0)
                    return false;
                if (localCoor[i] < -2.0)
                    return false;
            }
            //assert(false);
        }
        if (fabs(delta[0]) < EPS && fabs(delta[1]) < EPS) {
            break;
        }

        localCoor[0] += delta[0];
        localCoor[1] += delta[1];
    }
    for (int i = 0; i < 2; i++) {
        if (localCoor[i] > 1.0)
            return false;
        if (localCoor[i] < -1.0)
            return false;
    }
    return true;
}


// Class Methods
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

/***********************************************************************************************
 * \brief constructor of the class GaussQuadratureOnTriangle
 * \param[in] _triangle triangle to integrate on
 * \param[in] _nGaussPoints number of gauss points
 * \author Tianyang Wang
 ***********/
GaussQuadratureOnTriangle::GaussQuadratureOnTriangle(double *_triangle, int _numGaussPoints) :
        triangle(_triangle), numGaussPoints(_numGaussPoints) {
    if (numGaussPoints == 3) {
        gaussPointsLocal = EMPIRE::MathLibrary::triGaussPoints3;
        weights = EMPIRE::MathLibrary::triWeights3;
    } else if (numGaussPoints == 6) {
        gaussPointsLocal = EMPIRE::MathLibrary::triGaussPoints6;
        weights = EMPIRE::MathLibrary::triWeights6;
    } else if (numGaussPoints == 7) {
        gaussPointsLocal = EMPIRE::MathLibrary::triGaussPoints7;
        weights = EMPIRE::MathLibrary::triWeights7;
    } else if (numGaussPoints == 12) {
        gaussPointsLocal = EMPIRE::MathLibrary::triGaussPoints12;
        weights = EMPIRE::MathLibrary::triWeights12;
    } else {
        assert(false);
    }
    gaussPointsGlobal = new double[numGaussPoints * 3];
    for (int i = 0; i < numGaussPoints; i++) {
    	EMPIRE::MathLibrary::computeGlobalCoorInTriangle(triangle, &gaussPointsLocal[i * 3], &gaussPointsGlobal[i * 3]);
    }
    area = EMPIRE::MathLibrary::computeAreaOfTriangle(triangle);
}


/***********************************************************************************************
 * \brief Destructor of the class.
 * \author Tianyang Wang
 ***********/
GaussQuadratureOnTriangle::~GaussQuadratureOnTriangle() {
    delete[] gaussPointsGlobal;
}

/***********************************************************************************************
 * \brief To define the integrand function
 * \param[in] _integrandFunc Object of the class integrand functions
 * \author Tianyang Wang
 ***********/
void GaussQuadratureOnTriangle::setIntegrandFunc(IntegrandFunction *_integrandFunc) {
    integrandFunc = _integrandFunc;
}

/***********************************************************************************************
 * \brief To Compute integral on the triangle
 * \author Tianyang Wang
 ***********/
double GaussQuadratureOnTriangle::computeIntegral() {
    double toReturn = 0;
    for (int i = 0; i < numGaussPoints; i++)
        toReturn += weights[i] * (*integrandFunc)(&gaussPointsGlobal[i * 3]);
    toReturn *= area;
    return toReturn;
}

/***********************************************************************************************
 * \brief constructor of the class GaussQuadratureOnQuad
 * \param[in] _quad quad to integrate on
 * \param[in] _nGaussPoints number of gauss points
 * \author Tianyang Wang
 ***********/
GaussQuadratureOnQuad::GaussQuadratureOnQuad(double *_quad, int _numGaussPoints) :
        quad(_quad), numGaussPoints(_numGaussPoints) {
    if (numGaussPoints == 1) {
        gaussPointsLocal = EMPIRE::MathLibrary::quadGaussPoints1;
        weights = EMPIRE::MathLibrary::quadWeights1;
    } else if (numGaussPoints == 4) {
        gaussPointsLocal = EMPIRE::MathLibrary::quadGaussPoints4;
        weights = EMPIRE::MathLibrary::quadWeights4;
    } else if (numGaussPoints == 9) {
        gaussPointsLocal = EMPIRE::MathLibrary::quadGaussPoints9;
        weights = EMPIRE::MathLibrary::quadWeights9;
    } else {
        assert(false);
    }
    gaussPointsGlobal = new double[numGaussPoints * 3];
    detJ = new double[numGaussPoints];
    for (int i = 0; i < numGaussPoints; i++) {
    	EMPIRE::MathLibrary::computeGlobalCoorInQuad(quad, &gaussPointsLocal[i * 2], &gaussPointsGlobal[i * 3]);
    }
    for (int i = 0; i < numGaussPoints; i++) {
        detJ[i] = EMPIRE::MathLibrary::computeDetJOfQuad(quad, &gaussPointsLocal[i * 2]);
    }
}

/***********************************************************************************************
 * \brief Destructor of the class.
 * \author Tianyang Wang
 ***********/
GaussQuadratureOnQuad::~GaussQuadratureOnQuad() {
    delete[] gaussPointsGlobal;
    delete[] detJ;
}

/***********************************************************************************************
 * \brief To define the integrand function
 * \param[in] _integrandFunc Object of the class integrand functions
 * \author Tianyang Wang
 ***********/
void GaussQuadratureOnQuad::setIntegrandFunc(IntegrandFunction *_integrandFunc) {
    integrandFunc = _integrandFunc;
}

/***********************************************************************************************
 * \brief To Compute integral on the quad
 * \author Tianyang Wang
 ***********/
double GaussQuadratureOnQuad::computeIntegral() {
    double toReturn = 0;
    for (int i = 0; i < numGaussPoints; i++)
        toReturn += detJ[i] * weights[i] * (*integrandFunc)(&gaussPointsGlobal[i * 3]);
    return toReturn;
}

/***********************************************************************************************
 * \brief Returns the scalar/vector value from the linear combination of the values on the vertices of the element with the shape functions
 * \param[in] _nNodes The number of nodes of the element
 * \param[in] _nValue Takes value 1 for a scalar, or the values 2-3 for a vector
 * \param[in] _values The values on the nodes of the element
 * \param[in] _coords The coordinates in the element
 * \param[in/out] _returnValue The resulting linear combination of the nodal values
 * \author Chenshen Wu
 ***********/
void computeLinearCombinationValueFromVerticesValues(int _nNodes, int _nValue,
        const double *_values, const double* _coords, double *_returnValue) {
    double shapeFuncs[4];
    computeLowOrderShapeFunc(_nNodes, _coords, shapeFuncs);
    computeLinearCombination(_nNodes, _nValue, _values, shapeFuncs, _returnValue);
}

/***********************************************************************************************
 * \brief Computes the value of a data field in the interior of an element
 * \param[in] _nNodes The number of nodes of the element
 * \param[in] _nValue Takes value 1 for a scalar, or the values 2-3 for a vector
 * \param[in] _values The values on the nodes of the element
 * \param[in] _shapeFuncs The values of the shape functions at the interior of the element
 * \param[in/out] _returnValue The resulting linear combination of the nodal values
 * \author Chenshen Wu
 ***********/
void computeLinearCombination(int _nNodes, int _nValue, const double *_values,
        const double *_shapeFuncs, double *_returnValue) {

    for (int i = 0; i < _nValue; i++) {
        _returnValue[i] = 0;
        for (int j = 0; j < _nNodes; j++)
            _returnValue[i] += _values[j * _nValue + i] * _shapeFuncs[j];
    }

}

/***********************************************************************************************
 * \brief Computes the values of the low order shape functions (linear for triangle and bilinear
 *        for the quadrilateral)
 * \param[in] _nNodes The number of nodes in the element level
 * \param[in] _coords The coordinates of the point where to evaluate the shape functions
 * \param[in/out] _shapeFuncs The evaluated shape functions
 ***********/
void computeLowOrderShapeFunc(int _nNodes, const double *_coords, double *_shapeFuncs) {
    assert(_coords!=NULL);
    assert(_shapeFuncs!=NULL);
    if (_nNodes == 3) {
        _shapeFuncs[0] = 1 - _coords[0] - _coords[1];
        _shapeFuncs[1] = _coords[0];
        _shapeFuncs[2] = _coords[1];
    } else {
        _shapeFuncs[0] = (1 - _coords[0]) / 2 * (1 - _coords[1]) / 2;
        _shapeFuncs[1] = (1 + _coords[0]) / 2 * (1 - _coords[1]) / 2;
        _shapeFuncs[2] = (1 + _coords[0]) / 2 * (1 + _coords[1]) / 2;
        _shapeFuncs[3] = (1 - _coords[0]) / 2 * (1 + _coords[1]) / 2;
    }

}

/***********************************************************************************************
 * \brief Compute local coordinates of a point in a triangle in a 2D space
 * \param[in] _coordsTriangle, coordinates of the triangle. double[6].
 * \param[in] _coordsNode, coordinates of the point. double[2].
 * \param[out] _localCoords, local coordinates of the point. double[3]
 * \return a boolean saying whether the point is inside the triangle or not (true of false)
 * \author Chenshen Wu
 ***********/
bool computeLocalCoordsInTriangle(const double *_coordsTri, const double *_coordsNode,
        double* _localCoords) {
    assert(_coordsTri!=NULL);
    assert(_coordsNode!=NULL);

    double area = computeAreaTriangle(_coordsTri[2] - _coordsTri[0], _coordsTri[3] - _coordsTri[1],
            0, _coordsTri[4] - _coordsTri[0], _coordsTri[5] - _coordsTri[1], 0);
    double area1 = computeAreaTriangle(_coordsTri[2] - _coordsNode[0],
            _coordsTri[3] - _coordsNode[1], 0, _coordsTri[4] - _coordsNode[0],
            _coordsTri[5] - _coordsNode[1], 0);
    double area2 = computeAreaTriangle(_coordsTri[0] - _coordsNode[0],
            _coordsTri[1] - _coordsNode[1], 0, _coordsTri[4] - _coordsNode[0],
            _coordsTri[5] - _coordsNode[1], 0);
    _localCoords[0] = area1 / area;
    _localCoords[1] = area2 / area;
    if (_localCoords[0] < 0 || _localCoords[0] > 1 || _localCoords[1] < 0 || _localCoords[1] > 1
            || _localCoords[0] + _localCoords[1] > 1)
        return false;
    return true;
}

/***********************************************************************************************
 * \brief Compute local coordinates of a point in a quadriliteral in a 2D space by solving a nonlinear
 *        system using the Newton-Raphson scheme
 * \param[in] _coordsQuad Coordinates of the quadriliteral. double[8].
 * \param[in] _coordsNode Coordinates of the point. double[2].
 * \param[out] _localCoords local coordinates of the point. double[2]
 * \return a boolean saying whether the point is inside the quadriliteral or not (true of false)
 * \author Chenshen Wu
 ***********/
bool computeLocalCoordsInQuad(const double *_coordsQuad, const double *_coordsNode,
        double* _localCoords) {
    /*
     * So we use two coordinates among x, y, z.
     * This indicates projection.
     * Choose among planes (xy or yz or zx) the one which has the smallest angle
     * with the quad normal.
     */
    assert(_coordsQuad!=NULL);
    assert(_coordsNode!=NULL);
    double x[4];
    double y[4];
    for (int i = 0; i < 4; i++) {
        x[i] = _coordsQuad[i * 2];
        y[i] = _coordsQuad[i * 2 + 1];
        if (x[i] == _coordsNode[0] & y[i] == _coordsNode[1]) {
            _localCoords[0] = (i == 0 || i == 3) ? (-1) : (1);
            _localCoords[1] = (i == 0 || i == 1) ? (-1) : (1);
            return true;
        }
    }
    double x0 = _coordsNode[0];
    double y0 = _coordsNode[1];

    double a1 = x[0] + x[1] + x[2] + x[3] - 4.0 * x0;
    double b1 = -x[0] + x[1] + x[2] - x[3];
    double c1 = -x[0] - x[1] + x[2] + x[3];
    double d1 = x[0] - x[1] + x[2] - x[3];

    double a2 = y[0] + y[1] + y[2] + y[3] - 4.0 * y0;
    double b2 = -y[0] + y[1] + y[2] - y[3];
    double c2 = -y[0] - y[1] + y[2] + y[3];
    double d2 = y[0] - y[1] + y[2] - y[3];

    double delta[2];
    double J_T[4]; // transpose of Jacobian --- to use column major in lapack
    double F[2]; // -F
    _localCoords[0] = 0;
    _localCoords[1] = 0;
//int dummy[2];
    const double EPS = 1E-15;

    const int MAX_ITER_NUM = 100;
    for (int i = 0; i < MAX_ITER_NUM; i++) {
        J_T[0] = b1 + d1 * _localCoords[1];
        J_T[2] = c1 + d1 * _localCoords[0];
        J_T[1] = b2 + d2 * _localCoords[1];
        J_T[3] = c2 + d2 * _localCoords[0];
        F[0] = a1 + b1 * _localCoords[0] + c1 * _localCoords[1]
                + d1 * _localCoords[0] * _localCoords[1];
        F[1] = a2 + b2 * _localCoords[0] + c2 * _localCoords[1]
                + d2 * _localCoords[0] * _localCoords[1];
        delta[0] = -F[0];
        delta[1] = -F[1];

        solve2x2LinearSystem(J_T, delta, EPS);
        if (fabs(delta[0]) < EPS && fabs(delta[1]) < EPS) {
            assert(i < 100);
            break;
        }
        _localCoords[0] += delta[0];
        _localCoords[1] += delta[1];
    }
    for (int i = 0; i < 2; i++) {
        if (_localCoords[i] > 1.0)
            return false;
        if (_localCoords[i] < -1.0)
            return false;
    }
    return true;
}




// IGA Integration

/***********************************************************************************************
 * \brief Constructor
 * param[in] _numGaussPoints, number of Gauss points
 * \author Chenshen Wu
 ***********/
IGAGaussQuadratureOnTriangle::IGAGaussQuadratureOnTriangle(int _numGaussPoints) :
        IGAGaussQuadrature(_numGaussPoints) {
    switch (_numGaussPoints) {
    case 1:
        gaussPoints = IGAtriGaussPoints1;
        weights = IGAtriWeights1;
        break;
    case 3:
        gaussPoints = IGAtriGaussPoints3;
        weights = IGAtriWeights3;
        break;
    case 4:
        gaussPoints = IGAtriGaussPoints4;
        weights = IGAtriWeights4;
        break;
    case 6:
        gaussPoints = IGAtriGaussPoints6;
        weights = IGAtriWeights6;
        break;
    case 7:
        gaussPoints = IGAtriGaussPoints7;
        weights = IGAtriWeights7;
        break;
    case 12:
        gaussPoints = IGAtriGaussPoints12;
        weights = IGAtriWeights12;
        break;
    case 13:
        gaussPoints = IGAtriGaussPoints13;
        weights = IGAtriWeights13;
        break;
    case 16:
        gaussPoints = IGAtriGaussPoints16;
        weights = IGAtriWeights16;
        break;
    default:
        ERROR_OUT()<<"Number of Gauss Points for Triangle = " << numGaussPoints << "doesn't exist! Please choose from 1,3,4,6,7,12,13,16." << std::endl;
        exit(EXIT_FAILURE);

    }
}

/***********************************************************************************************
 * \brief Constructor
 * param[in] _numGaussPoints, number of Gauss points
 * \author Chenshen Wu
 ***********/
IGAGaussQuadratureOnQuad::IGAGaussQuadratureOnQuad(int _numGaussPoints) :
        IGAGaussQuadrature(_numGaussPoints) {
    switch (_numGaussPoints) {
    case 1:
        gaussPoints = IGAquadGaussPoints1;
        weights = IGAquadWeights1;
        break;
    case 4:
        gaussPoints = IGAquadGaussPoints4;
        weights = IGAquadWeights4;
        break;
    case 9:
        gaussPoints = IGAquadGaussPoints9;
        weights = IGAquadWeights9;
        break;
    case 16:
        gaussPoints = IGAquadGaussPoints16;
        weights = IGAquadWeights16;
        break;
    case 25:
        gaussPoints = IGAquadGaussPoints25;
        weights = IGAquadWeights25;
        break;
    default:
        ERROR_OUT()<<"Number of Gauss Points for Quadrilateral = " << numGaussPoints << "doesn't exist! Please choose from 1,4,9,16,25." << std::endl;
        exit(EXIT_FAILURE);
    }
}



}
}
