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
    // Cnn and Cnr matrices
    MathLibrary::SparseMatrix<double> *Cnn;

    // CNR matrices
    MathLibrary::SparseMatrix<double> *Cnr;

    // Master and slave sizes of Cnr
    size_t size_N;
    size_t size_R;

    // Flag on whether the expanded version of the coupling matrices is assumed
    bool isExpanded;

    // Vector containing all indices of empty rows in Cnn
    std::vector<int> indexEmptyRowCnn;

public:

    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _size_N Master size of Cnr
     * \param[in] _size_R Slave size of Cnr
     * \author Andreas Apostolatos
     ***********/
    IGAMortarCouplingMatrices(int _size_N , int _size_R, bool _isExpanded);

    /***********************************************************************************************
     * \brief Destructor
     * \author Andreas Apostolatos
     ***********/
    ~IGAMortarCouplingMatrices();

    /***********************************************************************************************
     * \brief Add value in the Cnn matrix
     * \param[in] _row row of added value
     * \param[in] _column column of added value
     * \param[in] value value to be added
     * \author Andreas Apostolatos
     ***********/
    void addCNNValue(int _row , int _column , double value) {
        (*Cnn)(_row, _column) += value;
        Cnn->setFactorization(false);
        Cnn->setDetermined(false);
    }

    /***********************************************************************************************
     * \brief Add value in the Cnr matrix
     * \param[in] _row row of added value
     * \param[in] _column column of added value
     * \param[in] value value to be added
     * \author Andreas Apostolatos
     ***********/
    void addCNRValue(int _row , int _column , double value) {
        (*Cnr)(_row, _column) += value;
    }

    /***********************************************************************************************
     * \brief Set value in the Cnn matrix
     * \param[in] _row row of value
     * \param[in] _column column of value
     * \param[in] value value to be set
     * \author Andreas Apostolatos
     ***********/
    void setValue(int _row , int _column , double value) {
        (*Cnn)(_row, _column) = value;
        Cnn->setFactorization(false);
        Cnn->setDetermined(false);
    }

    /***********************************************************************************************
     * \brief apply Dirichlet boundary conditions
     * \param[in] clampedIds clamped Dofs
     * \author Ragnar Björnsson
     ***********/
//    void applyDirichletBCs(std::vector<int> clampedIds);

    /***********************************************************************************************
     * \brief Get the size of C_nn
     * \author Andreas Apostolatos
     ***********/
    int getSizeN() {
        return size_N;
    }

    /***********************************************************************************************
     * \brief Get the slave size of C_nr
     * \author Andreas Apostolatos
     ***********/
    int getSizeR() {
        return size_R;
    }

    /***********************************************************************************************
     * \brief Get the flag on whether the expanded version of the coupling matrices is assumed
     * \author Andreas Apostolatos
     ***********/
    int getIsExpanded() {
        return isExpanded;
    }

    /***********************************************************************************************
     * \brief get CNN matrix
     * \author Andreas Apostolatos
     ***********/
    MathLibrary::SparseMatrix<double>* getCnn() {
        return Cnn;
    }

    /***********************************************************************************************
     * \brief get CNR matrix
     * \author Andreas Apostolatos
     ***********/
    MathLibrary::SparseMatrix<double>* getCnr() {
        return Cnr;
    }

    /***********************************************************************************************
     * \brief factorize the Cnn matrix
     * \author Andreas Apostolatos
     ***********/
    void factorizeCnn();

    /***********************************************************************************************
     * \brief Enforce consistency on the correct CNN matrix
     * \author Andreas Apostolatos
     ***********/
    void enforceCnn();

    /***********************************************************************************************
     * \brief Get indices of rows that are empty for the correct CNN matrix
     * \author Andreas Apostolatos
     ***********/
    std::vector<int> getIndexEmptyRowCnn() {
        return indexEmptyRowCnn;
    }

    /***********************************************************************************************
     * \brief Delete a row of the CNN matrix
     * \author Andreas Apostolatos
     ***********/
    void deleterow(int row) {
        Cnn->deleteRow(row);
        Cnn->setFactorization(false);
        Cnn->setDetermined(false);
    }
};
}


#endif // IGAMORTARCOUPLINGMATRICES_H
