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

#ifndef DEBUGMATH_H_
#define DEBUGMATH_H_
#include <string>

namespace EMPIRE {
namespace MathLibrary {

/***********************************************************************************************
 * \brief Print the coordinates of an element, used in debugging
 * \param[in] elem the element
 * \param[in] size 3---triangle, 4---quadrilateral
 * \author Tianyang Wang
 ***********/
void printElem(const double *elem, int size);
/***********************************************************************************************
 * \brief Print the coordinates of a point, used in debugging
 * \param[in] p the point
 * \author Tianyang Wang
 ***********/
void printPoint(const double *p);
/***********************************************************************************************
 * \brief print a diagonal matrix
 * \param[in] A the matrix
 * \param[in] numRows number of rows of A
 * \author Tianyang Wang
 ***********/
void printDiagonalMatrix(const double *A, int numRows);
/***********************************************************************************************
 * \brief print a diagonal matrix
 * \param[in] filename the file name to which the matrix should be written to.
 * \param[in] A the matrix
 * \param[in] numRows number of rows of A
 * \author Aditya Ghantasala
 ***********/
void printFullToFile(std::string filename, const double *A, int numRows);
/***********************************************************************************************
 * \brief print a general matrix
 * \param[in] A the matrix
 * \param[in] numRows number of rows of A
 * \param[in] numCols number of columns of A
 * \author Tianyang Wang
 ***********/
void printGeneralMatrix(const double *A, int numRows, int numCols);
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
        int numCols);
/***********************************************************************************************
 * \brief print a symmetric CSR formatted matrix
 * \param[in] A A
 * \param[in] IA IA
 * \param[in] JA JA
 * \param[in] n nxn matrix
 * \author Tianyang Wang
 ***********/
void printCSRMatrixSymmetric(const double *A, const int *IA, const int *JA, int n);



}
}


#endif /* DEBUGMATH_H_ */
