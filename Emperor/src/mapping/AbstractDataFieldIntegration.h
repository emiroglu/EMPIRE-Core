/*  Copyright &copy; 2013, TU Muenchen, Chair of Structural Analysis,
 *  Fabien Pean, Munich
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
 * \file AbstractDataFieldIntegration.h
 * This file holds the class AbstractDataFieldIntegration
 * \date 1/4/2015
 **************************************************************************************************/

#ifndef ABSTRACTDATAFIELDINTEGRATION_H_
#define ABSTRACTDATAFIELDINTEGRATION_H_

#include "MathLibrary.h"

namespace EMPIRE {

class AbstractDataFieldIntegration {
public:
    /***********************************************************************************************
     * \brief Destructor
     * \author Fabien Pean
     ***********/
	virtual ~AbstractDataFieldIntegration() {
		delete massMatrix;
	};
    /***********************************************************************************************
     * \brief Do integration (massMatrix*tractions=forces).
     * \param[in] tractions tractions (in one direction)
     * \param[out] forces forces (in one direction)
     * \author Tianyang Wang
     ***********/
    void integrate(const double *tractions, double *forces) {
        // This routine supports only one-based indexing of the input arrays.
        // Edit Aditya
        massMatrix->determineCSR();
        massMatrix->mulitplyVec(false,const_cast<double *>(tractions),forces,numNodes);
    };
    /***********************************************************************************************
     * \brief Do deintegration (massMatrix^(-1)*forces=tractions).
     * \param[in] forces forces (in one direction)
     * \param[out] tractions tractions (in one direction)
     * \author Tianyang Wang
     ***********/
    void deIntegrate(const double *forces, double *tractions) {
        // Edit Aditya
        massMatrix->solve(tractions,const_cast<double*>(forces));
    };
protected:
    /***********************************************************************************************
     * \brief Constructor protected to prevent explicit instantiation
     * \author Fabien Pean
     ***********/
	AbstractDataFieldIntegration() {
	};
    /// massMatrix csr format
    EMPIRE::MathLibrary::SparseMatrix<double> *massMatrix;
    /// number of nodes
    int numNodes;
};

} /* namespace EMPIRE */
#endif /* ABSTRACTDATAFIELDINTEGRATION_H_ */
