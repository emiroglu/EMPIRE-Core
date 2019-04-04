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

void computeShapeFuncOfQuad(const double *xi_eta, double *shapeFuncValues) {
    shapeFuncValues[0] = 0.25 * (1.0 - xi_eta[0]) * (1.0 - xi_eta[1]);
    shapeFuncValues[1] = 0.25 * (1.0 + xi_eta[0]) * (1.0 - xi_eta[1]);
    shapeFuncValues[2] = 0.25 * (1.0 + xi_eta[0]) * (1.0 + xi_eta[1]);
    shapeFuncValues[3] = 0.25 * (1.0 - xi_eta[0]) * (1.0 + xi_eta[1]);
}

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

void computeGlobalCoorInTriangle(const double *triangle, const double *localCoor,
        double *globalCoor) {
    for (int i = 0; i < 3; i++) {
        globalCoor[i] = 0;
        for (int j = 0; j < 3; j++)
            globalCoor[i] += localCoor[j] * triangle[i + j * 3];
    }
}

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
    if (fabs(localCoor[0] + localCoor[1] + localCoor[2] -1.0) > 1E-9) {
        cout << "Error in computing local coordinates in triangle!" << endl;
        cout << "Triangle: " << endl;
        for (int i=0; i<3; i++) {
            cout << "   " << triangle[i*3+0] << "   " << triangle[i*3+1] << "   " << triangle[i*3+2] << endl;
        }
        cout << "Point: " << endl;
        cout << "   " << point[0] << "   " << point[1] << "   " << point[2] << endl;
        cout << "Local coordinates: " << endl;
        cout << "   " << localCoor[0] << "   " << localCoor[1] << "   " << localCoor[2] << endl;
    }
    assert(fabs(localCoor[0] + localCoor[1] + localCoor[2] -1.0) < 1E-9);
    localCoor[0] = 1.0 - localCoor[1] - localCoor[2];

    for (int i = 0; i < 3; i++) {
        if (localCoor[i] > 1.0)
        	return false;
        if (localCoor[i] < 0.0)
            return false;
    }
    return true;
}

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

GaussQuadratureOnTriangle::~GaussQuadratureOnTriangle() {
    delete[] gaussPointsGlobal;
}

void GaussQuadratureOnTriangle::setIntegrandFunc(IntegrandFunction *_integrandFunc) {
    integrandFunc = _integrandFunc;
}

double GaussQuadratureOnTriangle::computeIntegral() {
    double toReturn = 0;
    for (int i = 0; i < numGaussPoints; i++)
        toReturn += weights[i] * (*integrandFunc)(&gaussPointsGlobal[i * 3]);
    toReturn *= area;
    return toReturn;
}

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

void GaussQuadratureOnQuad::setIntegrandFunc(IntegrandFunction *_integrandFunc) {
    integrandFunc = _integrandFunc;
}

double GaussQuadratureOnQuad::computeIntegral() {
    double toReturn = 0;
    for (int i = 0; i < numGaussPoints; i++)
        toReturn += detJ[i] * weights[i] * (*integrandFunc)(&gaussPointsGlobal[i * 3]);
    return toReturn;
}

void computeLinearCombinationValueFromVerticesValues(int _nNodes, int _nValue,
        const double *_values, const double* _coords, double *_returnValue) {
    double shapeFuncs[4];
    computeLowOrderShapeFunc(_nNodes, _coords, shapeFuncs);
    computeLinearCombination(_nNodes, _nValue, _values, shapeFuncs, _returnValue);
}

void computeLinearCombination(int _nNodes, int _nValue, const double *_values,
        const double *_shapeFuncs, double *_returnValue) {

    for (int i = 0; i < _nValue; i++) {
        _returnValue[i] = 0;
        for (int j = 0; j < _nNodes; j++)
            _returnValue[i] += _values[j * _nValue + i] * _shapeFuncs[j];
    }

}

void computeLowOrderShapeFunc(int _nNodes, const double *_coords, double *_shapeFuncs) {
    assert(_coords!=NULL);
    assert(_shapeFuncs!=NULL);
    if (_nNodes == 3) {
        _shapeFuncs[0] = 1 - _coords[0] - _coords[1];
        _shapeFuncs[1] = _coords[0];
        _shapeFuncs[2] = _coords[1];
    } else if (_nNodes == 4) {
        _shapeFuncs[0] = (1 - _coords[0]) / 2 * (1 - _coords[1]) / 2;
        _shapeFuncs[1] = (1 + _coords[0]) / 2 * (1 - _coords[1]) / 2;
        _shapeFuncs[2] = (1 + _coords[0]) / 2 * (1 + _coords[1]) / 2;
        _shapeFuncs[3] = (1 - _coords[0]) / 2 * (1 + _coords[1]) / 2;
    } else {
        ERROR_OUT() << "Low order basis functions are computed for only triangles or quadrilaterals" << endl;
        exit(-1);
    }
}

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

IGAGaussQuadratureOnTriangle::IGAGaussQuadratureOnTriangle(int _numGaussPoints) :
        IGAGaussQuadrature(_numGaussPoints, 2) {
    switch (_numGaussPoints) {
    case 1:
        setGaussPoints(IGACanonicalTriangleGaussPoints1);
        setGaussWeights(IGACanonicalTriangleWeights1);
        break;
    case 3:
        setGaussPoints(IGACanonicalTriangleGaussPoints3);
        setGaussWeights(IGACanonicalTriangleWeights3);
        break;
    case 4:
        setGaussPoints(IGACanonicalTriangleGaussPoints4);
        setGaussWeights(IGACanonicalTriangleWeights4);
        break;
    case 6:
        setGaussPoints(IGACanonicalTriangleGaussPoints6);
        setGaussWeights(IGACanonicalTriangleWeights6);
        break;
    case 7:
        setGaussPoints(IGACanonicalTriangleGaussPoints7);
        setGaussWeights(IGACanonicalTriangleWeights7);
        break;
    case 12:
        setGaussPoints(IGACanonicalTriangleGaussPoints12);
        setGaussWeights(IGACanonicalTriangleWeights12);
        break;
    case 13:
        setGaussPoints(IGACanonicalTriangleGaussPoints13);
        setGaussWeights(IGACanonicalTriangleWeights13);
        break;
    case 16:
        setGaussPoints(IGACanonicalTriangleGaussPoints16);
        setGaussWeights(IGACanonicalTriangleWeights16);
        break;
    default:
        ERROR_OUT() << "Selected number of Gauss points for triangle with symmetric locations = " << getNumGaussPoints() << " doesn't exist! Please choose from 1,3,4,6,7,12,13,16." << std::endl;
        exit(EXIT_FAILURE);
    }
}

IGAGaussQuadratureOnTriangleUsingDegeneratedQuadrilateral::IGAGaussQuadratureOnTriangleUsingDegeneratedQuadrilateral(int _numGaussPoints) :
        IGAGaussQuadrature(_numGaussPoints, 2) {
    switch (_numGaussPoints) {
    case 1:
        setGaussPoints(IGACanonicalTriangleUsingDegeneratedQuadrilateralGaussPoints1);
        setGaussWeights(IGACanonicalTriangleUsingDegeneratedQuadrilateralGaussWeights1);
        break;
    case 4:
        setGaussPoints(IGACanonicalTriangleUsingDegeneratedQuadrilateralGaussPoints4);
        setGaussWeights(IGACanonicalTriangleUsingDegeneratedQuadrilateralGaussWeights4);
        break;
    case 9:
        setGaussPoints(IGACanonicalTriangleUsingDegeneratedQuadrilateralGaussPoints9);
        setGaussWeights(IGACanonicalTriangleUsingDegeneratedQuadrilateralGaussWeights9);
        break;
    case 16:
        setGaussPoints(IGACanonicalTriangleUsingDegeneratedQuadrilateralGaussPoints16);
        setGaussWeights(IGACanonicalTriangleUsingDegeneratedQuadrilateralGaussWeights16);
        break;
    case 25:
        setGaussPoints(IGACanonicalTriangleUsingDegeneratedQuadrilateralGaussPoints25);
        setGaussWeights(IGACanonicalTriangleUsingDegeneratedQuadrilateralGaussWeights25);
        break;
    case 36:
        setGaussPoints(IGACanonicalTriangleUsingDegeneratedQuadrilateralGaussPoints36);
        setGaussWeights(IGACanonicalTriangleUsingDegeneratedQuadrilateralGaussWeights36);
        break;
    case 49:
        setGaussPoints(IGACanonicalTriangleUsingDegeneratedQuadrilateralGaussPoints49);
        setGaussWeights(IGACanonicalTriangleUsingDegeneratedQuadrilateralGaussWeights49);
        break;
    case 64:
        setGaussPoints(IGACanonicalTriangleUsingDegeneratedQuadrilateralGaussPoints64);
        setGaussWeights(IGACanonicalTriangleUsingDegeneratedQuadrilateralGaussWeights64);
        break;
    case 81:
        setGaussPoints(IGACanonicalTriangleUsingDegeneratedQuadrilateralGaussPoints81);
        setGaussWeights(IGACanonicalTriangleUsingDegeneratedQuadrilateralGaussWeights81);
        break;
    case 100:
        setGaussPoints(IGACanonicalTriangleUsingDegeneratedQuadrilateralGaussPoints100);
        setGaussWeights(IGACanonicalTriangleUsingDegeneratedQuadrilateralGaussWeights100);
        break;
    default:
        ERROR_OUT() << "Selected number of Gauss points for triangle using degenerated quadrilateral = " << getNumGaussPoints() << " doesn't exist! Please choose from 1, 4, 9, 16, 25, 36, 49, 64, 81, 100." << std::endl;
        exit(EXIT_FAILURE);
    }
}

IGAGaussQuadratureOnBiunitInterval::IGAGaussQuadratureOnBiunitInterval(int _numGaussPoints) :
    IGAGaussQuadrature(_numGaussPoints, 1) {
    switch (_numGaussPoints) {
    case 1:
        setGaussPoints(IGABiunitIntervalGaussPoints1);
        setGaussWeights(IGABiunitIntervalGaussWeights1);
        break;
    case 2:
        setGaussPoints(IGABiunitIntervalGaussPoints2);
        setGaussWeights(IGABiunitIntervalGaussWeights2);
        break;
    case 3:
        setGaussPoints(IGABiunitIntervalGaussPoints3);
        setGaussWeights(IGABiunitIntervalGaussWeights3);
        break;
    case 4:
        setGaussPoints(IGABiunitIntervalGaussPoints4);
        setGaussWeights(IGABiunitIntervalGaussWeights4);
        break;
    case 5:
        setGaussPoints(IGABiunitIntervalGaussPoints5);
        setGaussWeights(IGABiunitIntervalGaussWeights5);
        break;
    case 6:
        setGaussPoints(IGABiunitIntervalGaussPoints6);
        setGaussWeights(IGABiunitIntervalGaussWeights6);
        break;
    case 7:
        setGaussPoints(IGABiunitIntervalGaussPoints7);
        setGaussWeights(IGABiunitIntervalGaussWeights7);
        break;
    case 8:
        setGaussPoints(IGABiunitIntervalGaussPoints8);
        setGaussWeights(IGABiunitIntervalGaussWeights8);
        break;
    case 9:
        setGaussPoints(IGABiunitIntervalGaussPoints9);
        setGaussWeights(IGABiunitIntervalGaussWeights9);
        break;
    case 10:
        setGaussPoints(IGABiunitIntervalGaussPoints10);
        setGaussWeights(IGABiunitIntervalGaussWeights10);
        break;
    case 11:
        setGaussPoints(IGABiunitIntervalGaussPoints11);
        setGaussWeights(IGABiunitIntervalGaussWeights11);
        break;
    case 12:
        setGaussPoints(IGABiunitIntervalGaussPoints12);
        setGaussWeights(IGABiunitIntervalGaussWeights12);
        break;
    case 13:
        setGaussPoints(IGABiunitIntervalGaussPoints13);
        setGaussWeights(IGABiunitIntervalGaussWeights13);
        break;
    case 14:
        setGaussPoints(IGABiunitIntervalGaussPoints14);
        setGaussWeights(IGABiunitIntervalGaussWeights14);
        break;
    case 15:
        setGaussPoints(IGABiunitIntervalGaussPoints15);
        setGaussWeights(IGABiunitIntervalGaussWeights15);
        break;
    case 16:
        setGaussPoints(IGABiunitIntervalGaussPoints16);
        setGaussWeights(IGABiunitIntervalGaussWeights16);
        break;
    case 17:
        setGaussPoints(IGABiunitIntervalGaussPoints17);
        setGaussWeights(IGABiunitIntervalGaussWeights17);
        break;
    case 18:
        setGaussPoints(IGABiunitIntervalGaussPoints18);
        setGaussWeights(IGABiunitIntervalGaussWeights18);
        break;
    case 19:
        setGaussPoints(IGABiunitIntervalGaussPoints19);
        setGaussWeights(IGABiunitIntervalGaussWeights19);
        break;
    case 20:
        setGaussPoints(IGABiunitIntervalGaussPoints20);
        setGaussWeights(IGABiunitIntervalGaussWeights20);
        break;
    case 21:
        setGaussPoints(IGABiunitIntervalGaussPoints21);
        setGaussWeights(IGABiunitIntervalGaussWeights21);
        break;
    case 22:
        setGaussPoints(IGABiunitIntervalGaussPoints22);
        setGaussWeights(IGABiunitIntervalGaussWeights22);
        break;
    case 23:
        setGaussPoints(IGABiunitIntervalGaussPoints23);
        setGaussWeights(IGABiunitIntervalGaussWeights23);
        break;
    case 24:
        setGaussPoints(IGABiunitIntervalGaussPoints24);
        setGaussWeights(IGABiunitIntervalGaussWeights24);
        break;
    case 25:
        setGaussPoints(IGABiunitIntervalGaussPoints25);
        setGaussWeights(IGABiunitIntervalGaussWeights25);
        break;
    case 26:
        setGaussPoints(IGABiunitIntervalGaussPoints26);
        setGaussWeights(IGABiunitIntervalGaussWeights26);
        break;
    case 27:
        setGaussPoints(IGABiunitIntervalGaussPoints27);
        setGaussWeights(IGABiunitIntervalGaussWeights27);
        break;
    case 28:
        setGaussPoints(IGABiunitIntervalGaussPoints28);
        setGaussWeights(IGABiunitIntervalGaussWeights28);
        break;
    case 29:
        setGaussPoints(IGABiunitIntervalGaussPoints29);
        setGaussWeights(IGABiunitIntervalGaussWeights29);
        break;
    case 30:
        setGaussPoints(IGABiunitIntervalGaussPoints30);
        setGaussWeights(IGABiunitIntervalGaussWeights30);
        break;
    case 31:
        setGaussPoints(IGABiunitIntervalGaussPoints31);
        setGaussWeights(IGABiunitIntervalGaussWeights31);
        break;
    case 32:
        setGaussPoints(IGABiunitIntervalGaussPoints32);
        setGaussWeights(IGABiunitIntervalGaussWeights32);
        break;
    case 33:
        setGaussPoints(IGABiunitIntervalGaussPoints33);
        setGaussWeights(IGABiunitIntervalGaussWeights33);
        break;
    case 34:
        setGaussPoints(IGABiunitIntervalGaussPoints34);
        setGaussWeights(IGABiunitIntervalGaussWeights34);
        break;
    case 35:
        setGaussPoints(IGABiunitIntervalGaussPoints35);
        setGaussWeights(IGABiunitIntervalGaussWeights35);
        break;
    case 36:
        setGaussPoints(IGABiunitIntervalGaussPoints36);
        setGaussWeights(IGABiunitIntervalGaussWeights36);
        break;
    case 37:
        setGaussPoints(IGABiunitIntervalGaussPoints37);
        setGaussWeights(IGABiunitIntervalGaussWeights37);
        break;
    case 38:
        setGaussPoints(IGABiunitIntervalGaussPoints38);
        setGaussWeights(IGABiunitIntervalGaussWeights38);
        break;
    case 39:
        setGaussPoints(IGABiunitIntervalGaussPoints39);
        setGaussWeights(IGABiunitIntervalGaussWeights39);
        break;
    case 40:
        setGaussPoints(IGABiunitIntervalGaussPoints40);
        setGaussWeights(IGABiunitIntervalGaussWeights40);
        break;
    case 41:
        setGaussPoints(IGABiunitIntervalGaussPoints41);
        setGaussWeights(IGABiunitIntervalGaussWeights41);
        break;
    case 42:
        setGaussPoints(IGABiunitIntervalGaussPoints42);
        setGaussWeights(IGABiunitIntervalGaussWeights42);
        break;
    case 43:
        setGaussPoints(IGABiunitIntervalGaussPoints43);
        setGaussWeights(IGABiunitIntervalGaussWeights43);
        break;
    case 44:
        setGaussPoints(IGABiunitIntervalGaussPoints44);
        setGaussWeights(IGABiunitIntervalGaussWeights44);
        break;
    case 45:
        setGaussPoints(IGABiunitIntervalGaussPoints45);
        setGaussWeights(IGABiunitIntervalGaussWeights45);
        break;
    case 46:
        setGaussPoints(IGABiunitIntervalGaussPoints46);
        setGaussWeights(IGABiunitIntervalGaussWeights46);
        break;
    case 47:
        setGaussPoints(IGABiunitIntervalGaussPoints47);
        setGaussWeights(IGABiunitIntervalGaussWeights47);
        break;
    case 48:
        setGaussPoints(IGABiunitIntervalGaussPoints48);
        setGaussWeights(IGABiunitIntervalGaussWeights48);
        break;
    case 49:
        setGaussPoints(IGABiunitIntervalGaussPoints49);
        setGaussWeights(IGABiunitIntervalGaussWeights49);
        break;
    case 50:
        setGaussPoints(IGABiunitIntervalGaussPoints50);
        setGaussWeights(IGABiunitIntervalGaussWeights50);
        break;
    default:
        ERROR_OUT() << "Number of Gauss Points for biunit interval = " << getNumGaussPoints() << "doesn't exist! Please choose from 1,...,50." << std::endl;
        exit(EXIT_FAILURE);
    }
}

IGAGaussQuadratureOnBiunitQuadrilateral::IGAGaussQuadratureOnBiunitQuadrilateral(int _numGaussPoints) :
        IGAGaussQuadrature(_numGaussPoints, 2) {
    switch (_numGaussPoints) {
    case 1:
        setGaussPoints(IGABiunitQuadrilateralGaussPoints1);
        setGaussWeights(IGABiunitQuadrilateralWeights1);
        break;
    case 4:
        setGaussPoints(IGABiunitQuadrilateralGaussPoints4);
        setGaussWeights(IGABiunitQuadrilateralWeights4);
        break;
    case 9:
        setGaussPoints(IGABiunitQuadrilateralGaussPoints9);
        setGaussWeights(IGABiunitQuadrilateralWeights9);
        break;
    case 16:
        setGaussPoints(IGABiunitQuadrilateralGaussPoints16);
        setGaussWeights(IGABiunitQuadrilateralWeights16);
        break;
    case 25:
        setGaussPoints(IGABiunitQuadrilateralGaussPoints25);
        setGaussWeights(IGABiunitQuadrilateralWeights25);
        break;
    case 36:
        setGaussPoints(IGABiunitQuadrilateralGaussPoints36);
        setGaussWeights(IGABiunitQuadrilateralWeights36);
        break;
    case 49:
        setGaussPoints(IGABiunitQuadrilateralGaussPoints49);
        setGaussWeights(IGABiunitQuadrilateralWeights49);
        break;
    case 64:
        setGaussPoints(IGABiunitQuadrilateralGaussPoints64);
        setGaussWeights(IGABiunitQuadrilateralWeights64);
        break;
    case 81:
        setGaussPoints(IGABiunitQuadrilateralGaussPoints81);
        setGaussWeights(IGABiunitQuadrilateralWeights81);
        break;
    case 100:
        setGaussPoints(IGABiunitQuadrilateralGaussPoints100);
        setGaussWeights(IGABiunitQuadrilateralWeights100);
        break;
    default:
        ERROR_OUT() << "Selected number of Gauss points for quadrilateral = " << getNumGaussPoints() << "doesn't exist! Please choose from 1, 4, 9, 16, 25, 49, 64, 81, 100." << std::endl;
        exit(EXIT_FAILURE);
    }
}

}
}
