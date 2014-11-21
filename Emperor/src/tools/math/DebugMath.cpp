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

#include "DebugMath.h"
#include "Message.h"
#include <math.h>
#include <assert.h>
#include <iostream>
#include <iomanip>
#include <list>
using namespace std;

namespace EMPIRE {
namespace MathLibrary {

/***********************************************************************************************
 * \brief Print the coordinates of an element, used in debugging
 * \param[in] elem the element
 * \param[in] size 3---triangle, 4---quadrilateral
 * \author Tianyang Wang
 ***********/
void printElem(const double *elem, int size) {
    stringstream s;
    DEBUG_OUT() << "Element: " << endl;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < 3; j++)
            s << elem[i * 3 + j] << "  ";
        s << endl;
    }
    DEBUG_OUT() << s;
}


/***********************************************************************************************
 * \brief Print the coordinates of a point, used in debugging
 * \param[in] p the point
 * \author Tianyang Wang
 ***********/
void printPoint(const double *p) {
    stringstream s;
    s << "Point: ";
    for (int i = 0; i < 3; i++) {
        s << p[i] << "  ";
    }
    s << endl;
    DEBUG_OUT()<<s;
}

/***********************************************************************************************
 * \brief print a diagonal matrix
 * \param[in] A the matrix
 * \param[in] numRows number of rows of A
 * \author Tianyang Wang
 ***********/
void printDiagonalMatrix(const double *A, int numRows) {
    stringstream s;
    DEBUG_OUT() << "======================================================" << endl;
    s << "diagonal matrix:" << endl;
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numRows; j++) {
            if (i==j)
                s << setw(15) << A[i];
            else
                s << setw(15) << 0.0;
        }
        s << endl;
    }
    DEBUG_OUT() << s;
    DEBUG_OUT() << "======================================================" << endl;
}

/***********************************************************************************************
 * \brief print a general matrix
 * \param[in] A the matrix
 * \param[in] numRows number of rows of A
 * \param[in] numCols number of columns of A
 * \author Tianyang Wang
 ***********/
void printGeneralMatrix(const double *A, int numRows, int numCols) {
    stringstream s;
    DEBUG_OUT() << "======================================================" << endl;
    s << "general matrix:" << endl;
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numCols; j++) {
            s << setw(15) << A[i * numCols + j];
        }
        s << endl;
    }
    DEBUG_OUT() << s;
    DEBUG_OUT() << "======================================================" << endl;
}

/***********************************************************************************************
 * \brief print an unsymmetric CSR formatted matrix
 * \param[in] A A
 * \param[in] IA IA
 * \param[in] JA JA
 * \param[in] numRows number of rows of A
 * \param[in] numCols number of columns of A
 * \author Tianyang Wang
 ***********/
void printCSRMatrixUnsymmetric(const double *A, const int *IA, const int *JA, int numRows,
        int numCols) {
    double *matrix = new double[numRows * numCols];
    for (int i = 0; i < numRows * numCols; i++) {
        matrix[i] = 0.0;
    }
    for (int i = 0; i < numRows; i++) {
        int JA_row_begin = IA[i] - 1;
        int JA_row_end = IA[i + 1] - 1;
        for (int j = JA_row_begin; j < JA_row_end; j++) {
            int col = JA[j] - 1;
            matrix[i * numCols + col] = A[j];
        }
    }
    stringstream s;
    DEBUG_OUT() << "======================================================" << endl;
    s << "unsymmetric CSR matrix:" << endl;
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numCols; j++) {
            s << setw(15) << matrix[i * numCols + j];
        }
        s << endl;
    }
    DEBUG_OUT() << s;
    DEBUG_OUT() << "======================================================" << endl;
    delete[] matrix;
}

/***********************************************************************************************
 * \brief print a symmetric CSR formatted matrix
 * \param[in] A A
 * \param[in] IA IA
 * \param[in] JA JA
 * \param[in] n nxn matrix
 * \author Tianyang Wang
 ***********/
void printCSRMatrixSymmetric(const double *A, const int *IA, const int *JA, int n) {
    double *matrix = new double[n * n];
    for (int i = 0; i < n * n; i++) {
        matrix[i] = 0.0;
    }
    // 1. set up upper part of the matrix (same as unsymmetric matrix)
    for (int i = 0; i < n; i++) {
        int JA_row_begin = IA[i] - 1;
        int JA_row_end = IA[i + 1] - 1;
        for (int j = JA_row_begin; j < JA_row_end; j++) {
            int col = JA[j] - 1;
            matrix[i * n + col] = A[j];
        }
    }
    // 2. copy the upper part to the lower part
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < i; j++) {
            matrix[i * n + j] = matrix[j * n + i];
        }
    }
    stringstream s;
    DEBUG_OUT() << "======================================================" << endl;
    s << "symmetric CSR matrix:" << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            s << setw(15) << matrix[i * n + j];
        }
        s << endl;
    }
    DEBUG_OUT() << s;
    DEBUG_OUT() << "======================================================" << endl;
    delete[] matrix;
}

} // end namespace MathLibrary
} // end namespace EMPIRE


