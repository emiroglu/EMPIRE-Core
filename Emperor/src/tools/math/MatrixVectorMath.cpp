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
#ifdef USE_INTEL_MKL
#include <mkl.h>
#endif

#ifndef USE_INTEL_MKL
#include "cblas.h"
#endif

#include "MatrixVectorMath.h"
#include "ConstantsAndVariables.h"
#include "DebugMath.h"
using namespace std;

namespace EMPIRE {
namespace MathLibrary {

// Variables
// %%%%%%%%%%%%%%%%%%%%%%%


// Methods
// %%%%%%%%%%%%%%%%%%%%%%%
/***********************************************************************************************
 * \brief Compute the dot product of two dense vectors
 * \param[in] vec1 the 1st vector
 * \param[in] vec2 the 2nd vector
 * \return dot product
 * \author Stefan Sicklinger
 ***********/
double computeDenseDotProduct(const double *vec1, const double *vec2, const int elements) {
    return cblas_ddot(elements, vec1, 1, vec2, 1);
}

/***********************************************************************************************
 * \brief Compute the dot product of two dense vectors
 * \param[in] vec1 the 1st vector
 * \param[in] vec2 the 2nd vector
 * \param[in] elements number of elements in vec1 (vec2)
 * \return dot product
 * \author Stefan Sicklinger
 ***********/
double computeDenseDotProduct(const std::vector<double> &vec1, const std::vector<double> &vec2) {
    return cblas_ddot(vec1.size(), &vec1[0], 1, &vec2[0], 1);
}

/***********************************************************************************************
 * \brief Copy dense vector vec1 <- vec2
 * \param[in] vec1 the 1st vector
 * \param[in] vec2 the 2nd vector
 * \author Stefan Sicklinger
 ***********/
void copyDenseVector(double *vec1, const double *vec2, const int elements){
	cblas_dcopy(elements, vec2, 1, vec1, 1);
}

/***********************************************************************************************
 * \brief Compute Euclidean norm of vector
 * \param[in] vec1 the 1st vector
 * \param[in] elements number of elements in vec1
 * \return Euclidean norm
 * \author Stefan Sicklinger
 ***********/
double computeDenseEuclideanNorm(const double *vec1, const int elements){
	return cblas_dnrm2 (elements, vec1, 1);
}

/***********************************************************************************************
 * \brief Computes a vector-scalar product and adds the result to a vector. vec1 <- a*vec1 + vec2
 * \param[in] vec1 the 1st vector
 * \param[in] vec2 the 2nd vector
 * \param[in] a    scalar
 * \param[in] elements number of elements in vec1
 * \author Stefan Sicklinger
 ***********/
void computeDenseVectorAddition(double *vec1, const double *vec2 ,const double a, const int elements){
	cblas_daxpy (elements, a, vec2, 1, vec1, 1);
}

/***********************************************************************************************
 * \brief Computes vector scales by scalar vec1 <- vec1*a
 * \param[in] vec1 the 1st vector
 * \param[in] a    scalar
 * \param[in] elements number of elements in vec1
 * \author Stefan Sicklinger
 ***********/
void computeDenseVectorMultiplicationScalar(double *vec1 ,const double a, const int elements){
	cblas_dscal (elements, a, vec1, 1);
}

/***********************************************************************************************
 * \brief Computes the index function for the binomial coefficients (_i;_j)
 * \param[out] The index function for the binomial coefficients (_i;_j)
 * \param[in] _i The integer on the nominator
 * \param[in] _j The integer on the denominator
 * \author Andreas Apostolatos
 ***********/
int indexBinomialCoefficients(int _i, int _j) {
    return _i * 49 + _j;
}

/***********************************************************************************************
 * \brief Compute the square of the Euclidean distance of two points in n-D space.
 * \param[out] The square of the Euclidean distance between two points in the n-D space
 * \param[in] _length The length of the n-dimensional space
 * \param[in] _Pi The first point
 * \param[in] _Pj The second point
 * \author Andreas Apostolatos
 ***********/
double squareEuclideanDistance(int _length, double* _Pi, double* _Pj) {
    /*
     * Returns the square of the Euclidean distance of two points in n-D space.
     * The input arguments are 1D arrays holding the coordinate information of the points:
     *  _Pi = _Pj = double[_length]
     */

    double squareDistance = 0.0;
    for (int i = 0; i < _length; i++)
        squareDistance += pow(_Pi[i] - _Pj[i], 2.0);

    return squareDistance;
}

/***********************************************************************************************
 * \brief Compute the square of the 2-norm of a vector in the n-D space
 * \param[out] The square of the 2-norm of a vector in the n-D space
 * \param[in] _length The dimensinality of the n-D space
 * \param[in] _vector A vector in the nD space
 * \author Andreas Apostolatos
 ***********/
double square2normVector(int _length, double* _vector) {
    /*
     * Returns the square of the 2-norm of a vector in _length-dimensional space
     * _vector = double[_length]
     */

    // Initialize the square of the 2-norm of the _length-dimensional vector
    double vector2norm = 0.0;

    // Loop over all spatial dimensions
    for (int i = 0; i < _length; i++)
        vector2norm += pow(_vector[i], 2.0);

    // Return the 2-norm of the _length-dimensional vector
    return vector2norm;
}

/***********************************************************************************************
 * \brief Compute the dot product between two vectors in the n-D space
 * \param[out] The dot product between two vectors in the n-D space
 * \param[in] _length The dimensinality of the n-D space
 * \param[in] _vecI The 1st vector
 * \param[in] _vecJ The 2nd vector
 * \author Andreas Apostolatos
 ***********/
double dotProduct(int _length, double* _vecI, double* _vecJ) {
    // Initialize the dot product
    double dotProduct = 0.0;

    // Loop over all the entries of the vectors
    for (int i = 0; i < _length; i++)
        dotProduct += _vecI[i] * _vecJ[i];

    // Return the dot product value
    return dotProduct;
}

/***********************************************************************************************
* \brief Compute the cross product between two vectors in the 3-D space
* \param[in/out] _product The product of vector1 and vector 2
* \param[in] _v1 The 1st vector
* \param[in] _v2 The 2nd vector
* \author Chenshen Wu
***********/
void crossProduct(double* _product, double* _v1, double* _v2) {
    // Check input
    assert(_product != NULL);
    assert(_v1 != NULL);
    assert(_v2 != NULL);

    // Compute the cross product using the permutation tensor
    _product[0] = _v1[1] * _v2[2] - _v1[2] * _v2[1];
    _product[1] = _v1[2] * _v2[0] - _v1[0] * _v2[2];
    _product[2] = _v1[0] * _v2[1] - _v1[1] * _v2[0];

    // No return value
    return;
}

/***********************************************************************************************

 * \brief Solve a 2x2 linear equation system
 * \param[out] Flag on whether the linear system is solvable up to tolerance EPS or not
 * \param[in/out] _b The right-hand side vector where the solution is also stored, _b = double[2]
 * \param[in] _A The 2x2 matrix stored in a vector format namely A[i][j] = V[2 * i + j]
 * \author Andreas Apostolatos
 ***********/
bool solve2x2linearSystem(double* _b, double* _A) {
    /*
     * Solves a 2x2 linear equation system and stores the solution into the given right-hand side vector
     * _b = double[2]
     * _A = double[4], access of entry (i,j) by the rule (i,j) --> i * 2 + j
     */

    // Compute the partial determinants of the system
    double detA = _A[0] * _A[3] - _A[2] * _A[1];
    double det0 = _A[3] * _b[0] - _A[2] * _b[1];
    double det1 = _A[0] * _b[1] - _A[1] * _b[0];

    // Check the tolerance criterion
    if (fabs(detA) < EMPIRE::MathLibrary::EPS * fabs(det0))
    return false;
    if (fabs(detA) < EMPIRE::MathLibrary::EPS * fabs(det1))
    return false;
    if (detA == 0)
    return false;

    // Return the solution to the right-hand side vector
    _b[0] = det0 / detA;
    _b[1] = det1 / detA;

    // Return success
    return true;
}

/***********************************************************************************************
 * \brief Solve a 2x2 linear system by close form formula
 * \param[in] A the left hand side matrix
 * \param[in] b the right hand side vector
 * \param[out] b the solution is written to b
 * \return whether the determinant is zero or not
 * \author Tianyang Wang
 ***********/
bool solve2x2LinearSystem(const double *A, double *b, double EPS) {
    /*
     * The function is called by PolygonClipper::intersect() and computeLocalCoorInQuad()
     * understand the system in this way:
     * A = [v1, v2], Ax = b means solving alpha and beta in:
     * alpha * v1 + beta * v2 = b
     * detA is the area of the quad formed by v1 and v2
     */
    // A is column major
    double detA = A[0] * A[3] - A[2] * A[1];
    double det0 = A[3] * b[0] - A[2] * b[1];
    double det1 = A[0] * b[1] - A[1] * b[0];

    double v1LengthSquare = A[0] * A[0] + A[1] * A[1];
    double v2LengthSquare = A[2] * A[2] + A[3] * A[3];
    double max = (v1LengthSquare > v2LengthSquare) ? v1LengthSquare : v2LengthSquare;
    assert(max > 1E-20); // assume the length of the vector cannot be too small

    if (fabs(detA) < EPS * fabs(det0))
        return false;
    if (fabs(detA) < EPS * fabs(det1))
        return false;
    if (fabs(detA) < 1E-15 * max) // Check detA relative to max. 1E-15 is chosen to pass TestMortarMath::testClipping()
        return false; // det0 and det1 could be 0!!! Without this, it could fail!!!!!

    b[0] = det0 / detA;
    b[1] = det1 / detA;
    return true;
}

/***********************************************************************************************
 * \brief Solve a 3x3 linear system by close form formula, the system has one row which has all entries 1.
 * \brief Therefore, it cannot solve general 3x3 linear system.
 * \param[in] A the left hand side matrix
 * \param[in] b the right hand side vector
 * \param[in] planeToProject project to plane (case: {2:x-y ; 0:y-z ;1: z-x} )
 * \param[out] b the solution is written to b
 * \author Tianyang Wang
 ***********/
void solve3x3LinearSystem(const double *A, int planeToProject, double *b) {
    int XX = (planeToProject + 1) % 3;
    int YY = (planeToProject + 2) % 3;
    double detA1A2 = A[XX + 3] * A[YY + 6] - A[YY + 3] * A[XX + 6];
    double detA0A2 = A[XX] * A[YY + 6] - A[YY] * A[XX + 6];
    double detA0A1 = A[XX] * A[YY + 3] - A[YY] * A[XX + 3];
    double detA0b = A[XX] * b[YY] - A[YY] * b[XX];
    double detbA2 = b[XX] * A[YY + 6] - b[YY] * A[XX + 6];
    double detA1b = A[XX + 3] * b[YY] - A[YY + 3] * b[XX];
    double detA = detA1A2 - detA0A2 + detA0A1;
    double det0 = detA1A2 - detbA2 - detA1b;
    double det1 = detbA2 - detA0A2 + detA0b;
    double det2 = detA1b - detA0b + detA0A1;
    b[0] = det0 / detA;
    b[1] = det1 / detA;
    b[2] = det2 / detA;
}

/***********************************************************************************************
 * \brief Compute the matrix product between two matrices (for mortar, not so general)
 * \param[in] n n
 * \param[in] m m
 * \param[in] A the matrix with size n*n
 * \param[in] B the matrix with size n*m
 * \param[out] B B is overwritten by A*B (n*m)
 * \author Tianyang Wang
 ***********/
void computeMatrixProduct(int n, int m, const double *A, double *B) {
    double *C = new double[n * m];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            double tmp = 0.0;
            for (int k = 0; k < n; k++) {
                tmp += A[i * n + k] * B[j + k * m];
            }
            C[i * m + j] = tmp;
        }
    }
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            B[i * m + j] = C[i * m + j];
    delete[] C;
}

/***********************************************************************************************
 * \brief My own sparse matrix (csr format, non-symmetric) vector multiplication routine (A * x = y).
 *        The interface is simplified and compatible with mkl_dcsrmv.
 *        This routine supports only one-based indexing of the input arrays.
 * \param[in] trans 'N' --- no transpose, 'T' --- transpose
 * \param[in] numRows number of rows of matrix A
 * \param[in] numCols number of columns of matrix A
 * \param[in] A sparse matrix A in CSR format
 * \param[in] JA JA of A
 * \param[in] IA IA of A
 * \param[in] x vector x
 * \param[in] y vector y
 * \author Tianyang Wang
 ***********/
void dcsrmv(char trans, int numRows, int numCols, const double *A, const int *JA, const int *IA,
        const double *x, double *y) {
    if (trans == 'N') {
        for (int i = 0; i < numRows; i++)
            y[i] = 0.0;
        for (int i = 0; i < numRows; i++) {
            int JA_row_begin = IA[i] - 1;
            int JA_row_end = IA[i + 1] - 1;
            for (int j = JA_row_begin; j < JA_row_end; j++) {
                int col = JA[j] - 1;
                y[i] += A[j] * x[col];
            }
        }
    } else if (trans == 'T') {
        for (int i = 0; i < numCols; i++)
            y[i] = 0.0;
        for (int i = 0; i < numRows; i++) {
            int JA_row_begin = IA[i] - 1;
            int JA_row_end = IA[i + 1] - 1;
            for (int j = JA_row_begin; j < JA_row_end; j++) {
                int col = JA[j] - 1;
                y[col] += A[j] * x[i];
            }
        }
    } else {
        assert(false);
    }
}

/***********************************************************************************************
 * \brief My own sparse matrix (csr format, symmetric) vector multiplication routine (A * x = y).
 *        The interface is simplified and compatible with mkl_dcsrsymv.
 *        This routine supports only one-based indexing of the input arrays.
 * \param[in] n size of matrix A
 * \param[in] A sparse matrix A in CSR format (symmetric)
 * \param[in] IA IA of A
 * \param[in] JA JA of A
 * \param[in] x vector x
 * \param[in] y vector y
 * \author Tianyang Wang
 ***********/
void dcsrsymv(int n, const double *A, const int *IA, const int *JA, const double *x, double *y) {
    for (int i = 0; i < n; i++)
        y[i] = 0.0;
    // use the upper part with diagonal
    for (int i = 0; i < n; i++) {
        int JA_row_begin = IA[i] - 1;
        int JA_row_end = IA[i + 1] - 1;
        for (int j = JA_row_begin; j < JA_row_end; j++) {
            int col = JA[j] - 1;
            y[i] += A[j] * x[col];
        }
    }
    // use the lower part without diagonal
    for (int i = 0; i < n; i++) {
        int JA_row_begin = IA[i];
        int JA_row_end = IA[i + 1] - 1;
        for (int j = JA_row_begin; j < JA_row_end; j++) {
            int col = JA[j] - 1;
            y[col] += A[j] * x[i];
        }
    }
}

/***********************************************************************************************
 * \brief Solve a 2x2 linear system by close form formula
 * \param[in] A the left hand side matrix
 * \param[in] b the right hand side vector
 * \param[out] b the solution is written to b
 * \return whether the determinant is zero or not
 * \author Tianyang Wang
 ***********/
bool IGAsolve2x2LinearSystem(const double* _A, double* _b, double _EPS) {
// A is column major
    double detA = _A[0] * _A[3] - _A[2] * _A[1];
    double det0 = _A[3] * _b[0] - _A[2] * _b[1];
    double det1 = _A[0] * _b[1] - _A[1] * _b[0];
    if (fabs(detA) < _EPS * fabs(det0))
        return false;
    if (fabs(detA) < _EPS * fabs(det1))
        return false;
    if (detA == 0) // det0 and det1 could be 0!!!
        return false; // without this, it could fail!!!!!
    _b[0] = det0 / detA;
    _b[1] = det1 / detA;
    return true;
}

/***********************************************************************************************
 * \brief Solve 3x3 Linear system, Ax = b
 * \param[in] _A Square 3x3 matirx
 * \param[in/out] _b Right hand side vector which stores also the solution
 * \param[in] _EPS tolerance up to which matrix _A is regular
 * \param[out] Flag on the solvability of the system
 * \author Chenshen Wu
 ***********/
bool solve3x3LinearSystem(const double* _A, double* _b, double _EPS) {
// A is column major
    double A[9];
    double b[3];
    double detA = det3x3(_A);
    if (fabs(detA) < _EPS)
        return false;
    for (int i = 0; i < 3; i++)
        b[i] = _b[i];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 9; j++)
            A[j] = _A[j];
        for (int j = 0; j < 3; j++)
            A[j * 3 + i] = b[j];
        _b[i] = det3x3(A) / detA;
    }
    return true;

}

/***********************************************************************************************
 * \brief Computes the Determinant of a given 3x3 matrix
 * \param[in] _A Matrix
 * \param[out] The Determinant of the given matrix
 * \author Chenshen Wu
 ***********/
double det3x3(const double* _A) {
    return _A[0] * _A[4] * _A[8] + _A[1] * _A[5] * _A[6] + _A[2] * _A[3] * _A[7]
            - _A[0] * _A[5] * _A[7] - _A[1] * _A[3] * _A[8] - _A[2] * _A[4] * _A[6];
}

/***********************************************************************************************
 * \brief Computes the cross product of between two vectors with zero component on z direction
 * \param[in] _x1, x component of the first vector
 * \param[in] _y1, y component of the first vector
 * \param[in] _x2, x component of the second vector
 * \param[in] _y2, y component of the second vector
 * \param[out] Z component of the cross product between (x1,y1,0) and (x2,y2,0)
 * \author Chenshen Wu
 ***********/
double computeCrossProduct2D(double _x1, double _y1, double _x2, double _y2) {
    return _x1 * _y2 - _x2 * _y1;
}





//Classes
// %%%%%%%%%%%%%%%%%%%%%%%



}
}




