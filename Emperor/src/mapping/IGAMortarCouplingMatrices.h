/*  Copyright &copy; 2013, TU Muenchen, Chair of Structural Analysis,
 *  Ragnar Björnsson, Munich
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
/***********************************************************************************************//**
 * \file IGAMortarCouplingMatrices.h
 * This file holds the class IGAMortarCouplingMatrices.h
 * \date 6/8/2013
 **************************************************************************************************/

#ifndef IGAMORTARCOUPLINGMATRICES_H
#define IGAMORTARCOUPLINGMATRICES_H

#include "MathLibrary.h"

#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <set>
#include <algorithm>

namespace EMPIRE {

namespace MathLibrary {
template<class T> class SparseMatrix;
}

class IGAMortarCouplingMatrices {

private:
    // CNN matrices
    MathLibrary::SparseMatrix<double> *C_NN;
    MathLibrary::SparseMatrix<double> *C_NN_expanded;
    MathLibrary::SparseMatrix<double> *C_NN_BCs;

    // CNR matrices
    MathLibrary::SparseMatrix<double> *C_NR;
    MathLibrary::SparseMatrix<double> *C_NR_expanded;
    MathLibrary::SparseMatrix<double> *C_NR_BCs;

    size_t size_N;
    size_t size_R;

    bool isIGAWeakDirichletCurveConditions;
    bool isIGAWeakDirichletBoundaryConditions;
    bool isIGAPatchContinuityConditions;
    bool isDirichletBCs;
    bool isClampedDofs;

    std::vector<int> indexEmptyRowCnn;

public:

    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _size_N number of master nodes
     * \param[in] _size_R number of slave nodes
     * \author Ragnar Björnsson
     ***********/
    IGAMortarCouplingMatrices(int _size_N , int _size_R);

    /***********************************************************************************************
     * \brief Destructor
     * \author Ragnar Björnsson
     ***********/
    ~IGAMortarCouplingMatrices();

    /***********************************************************************************************
     * \brief add value to a location in the CNN matrix
     * \param[in] _row row of added value
     * \param[in] _column column of added value
     * \param[in] value value to be added
     * \author Ragnar Björnsson
     ***********/
    void addCNNValue(int _row , int _column , double value) {
        (*C_NN)(_row , _column) += value;
        C_NN->setFactorization(false);
    }

    /***********************************************************************************************
     * \brief add value to a location in the CNR matrix
     * \param[in] _row row of added value
     * \param[in] _column column of added value
     * \param[in] value value to be added
     * \author Ragnar Björnsson
     ***********/
    void addCNRValue(int _row , int _column , double value) {
        (*C_NR)(_row , _column) += value;
    }

    /***********************************************************************************************
     * \brief add value to a location in the C_NN_expanded matrix
     * \param[in] _row row of added value
     * \param[in] _column column of added value
     * \param[in] value value to be added
     * \author Ragnar Björnsson
     ***********/
    void addCNN_expandedValue(int _row , int _column , double value) {
        (*C_NN_expanded)(_row , _column) += value;
        C_NN_expanded->setFactorization(false);
    }

    /***********************************************************************************************
     * \brief set isIGAPatchCoupling and isClampedDofs and expand if any of them is true
     * \param[in] _isIGAPatchContinuityConditions Flag on whether weak patch continuity conditions are applied
     * \param[in] _isIGAWeakDirichletBoundaryConditions Flag on whether weak Dirichlet conditions are applied
     * \param[in] _isClampedDofs are any clamped nodes where not all directions are clamped
     * \author Ragnar Björnsson, Andreas Apostolatos, Altug Emiroglu
     ***********/
    void setIsIGAConditions(bool _isIGAWeakDirichletCurveConditions, bool _isIGAWeakDirichletBoundaryConditions, bool _isIGAPatchContinuityConditions, bool _isClampedDofs);

    /***********************************************************************************************
     * \brief expand CNN and CNR matrices to account for all 3 directions
     * \author Ragnar Björnsson
     ***********/
    void expandMatrices();

    /***********************************************************************************************
     * \brief set value of the correct CNN matrix
     * \param[in] _row row of value
     * \param[in] _column column of value
     * \param[in] value value to be set
     * \author Ragnar Björnsson
     ***********/
    void setValue(int _row , int _column , double value) {
        if(isDirichletBCs) {
            (*C_NN_BCs)(_row,_column) = value;
            C_NN_BCs->setFactorization(false);
        }
        else if( isIGAWeakDirichletCurveConditions || isIGAWeakDirichletBoundaryConditions || isIGAPatchContinuityConditions ) {
            (*C_NN_expanded)(_row,_column) = value;
            C_NN_expanded->setFactorization(false);
        }
        else {
            (*C_NN)(_row,_column) = value;
            C_NN->setFactorization(false);
        }
    }

    /***********************************************************************************************
     * \brief set isDirichletBCs and initialize BCs matrices if that is true
     * \param[in] _isDirichletBCs will Dirichlet Boundary conditions be applied
     * \author Ragnar Björnsson
     ***********/
    void setIsDirichletBCs(bool _isDirichletBCs);

    /***********************************************************************************************
     * \brief apply Dirichlet boundary conditions
     * \param[in] clampedIds clamped Dofs
     * \author Ragnar Björnsson
     ***********/
    void applyDirichletBCs(std::vector<int> clampedIds);

    /***********************************************************************************************
     * \brief get correct master size, either number of nodes or dofs
     * \author Ragnar Björnsson
     ***********/
    int getCorrectSizeN() {
        if( isIGAWeakDirichletCurveConditions || isIGAWeakDirichletBoundaryConditions || isIGAPatchContinuityConditions || isClampedDofs)
            return 3*size_N;
        else
            return size_N;
    }

    /***********************************************************************************************
     * \brief get correct slave size, either number of nodes or dofs
     * \author Ragnar Björnsson
     ***********/
    int getCorrectSizeR() {
        if( isIGAWeakDirichletCurveConditions || isIGAWeakDirichletBoundaryConditions || isIGAPatchContinuityConditions || isClampedDofs)
            return 3*size_R;
        else
            return size_R;
    }

    /***********************************************************************************************
     * \brief get correct CNN matrix
     * \author Ragnar Björnsson
     ***********/
    MathLibrary::SparseMatrix<double>* getCorrectCNN() {
        if(isDirichletBCs)
            return C_NN_BCs;
        else if( isIGAWeakDirichletCurveConditions || isIGAWeakDirichletBoundaryConditions || isIGAPatchContinuityConditions )
                return C_NN_expanded;
        else
            return C_NN;
    }

    /***********************************************************************************************
     * \brief get correct CNR matrix
     * \author Ragnar Björnsson
     ***********/
    MathLibrary::SparseMatrix<double>* getCorrectCNR() {
        if(isDirichletBCs)
            return C_NR_BCs;
        else if( isIGAWeakDirichletCurveConditions || isIGAWeakDirichletBoundaryConditions || isIGAPatchContinuityConditions )
            return C_NR_expanded;
        else
            return C_NR;
    }

    /***********************************************************************************************
     * \brief get correct CNN matrix when using conservative mapping (cannot use Dirichlet boundary conditions here)
     * \author Ragnar Björnsson
     ***********/
    MathLibrary::SparseMatrix<double>* getCorrectCNN_conservative() {
        if( isIGAWeakDirichletCurveConditions || isIGAWeakDirichletBoundaryConditions || isIGAPatchContinuityConditions )
                return C_NN_expanded;
        else {
            return C_NN;
        }
    }

    /***********************************************************************************************
     * \brief get correct CNR matrix when using conservative mapping (cannot use Dirichlet boundary conditions here)
     * \author Ragnar Björnsson
     ***********/
    MathLibrary::SparseMatrix<double>* getCorrectCNR_conservative() {
        if( isIGAWeakDirichletCurveConditions || isIGAWeakDirichletBoundaryConditions || isIGAPatchContinuityConditions )
            return C_NR_expanded;
        else
            return C_NR;
    }

    /***********************************************************************************************
     * \brief factorize the correct CNN matrix
     * \author Ragnar Björnsson
     ***********/
    void factorizeCorrectCNN();

    /***********************************************************************************************
     * \brief enforce consistency on the correct CNN matrix
     * \author Ragnar Björnsson
     ***********/
    void enforceCnn();

    /***********************************************************************************************
     * \brief get indices of rows that are empty for the correct CNN matrix
     * \author Ragnar Björnsson
     ***********/
    std::vector<int> getIndexEmptyRowCnn() {
        return indexEmptyRowCnn;
    }

    /***********************************************************************************************
     * \brief delete a row of the correct CNN matrix
     * \author Ragnar Björnsson
     ***********/
    void deleterow(int row) {
        if(isDirichletBCs) {
            C_NN_BCs->deleteRow(row);
            C_NN_BCs->setFactorization(false);
        }
        else if( isIGAWeakDirichletCurveConditions || isIGAWeakDirichletBoundaryConditions || isIGAPatchContinuityConditions ) {
            C_NN_expanded->deleteRow(row);
            C_NN_expanded->setFactorization(false);
        }
        else {
            C_NN->deleteRow(row);
            C_NN->setFactorization(false);
        }
    }
};
}


#endif // IGAMORTARCOUPLINGMATRICES_H
