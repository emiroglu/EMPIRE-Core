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
/***********************************************************************************************//**
 * \file DataFieldIntegration.h
 * This file holds the class DataFieldIntegration
 * \date 7/17/2013
 **************************************************************************************************/
#ifndef DATAFIELDINTEGRATION_H_
#define DATAFIELDINTEGRATION_H_

#include "AbstractDataFieldIntegration.h"

namespace EMPIRE {

namespace MathLibrary
{
	template<class T>
	class SparseMatrix;
}
class FEMesh;

/********//**
 * \brief Class DataFieldIntegration is an operator from traction to force or vice versa
 * \author Tianyang Wang
 ***********/
class DataFieldIntegration : public AbstractDataFieldIntegration {
public:
    /***********************************************************************************************
     * \brief Legacy constructor for backward compatibility
     * \param[in] _numNodes number of nodes
     * \param[in] _numElems number of elements
     * \param[in] _numNodesPerElem number of nodes per element
     * \param[in] _nodes coordinates of nodes
     * \param[in] _nodeIDs IDs of nodes
     * \param[in] _elems element table
     * \author Tianyang Wang
     ***********/
    DataFieldIntegration(int _numNodes, int _numElems, const int *_numNodesPerElem,
            const double *_nodes, const int *_nodeIDs, const int *_elems);
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _mesh FE mesh on which integration is proceeding
     * \author Fabien Pean
     ***********/
    DataFieldIntegration(FEMesh* mesh);
    /***********************************************************************************************
     * \brief Destructor
     * \author Tianyang Wang
     ***********/
    virtual ~DataFieldIntegration();

private:
    /// number of Gauss points used for computing triangle element mass matrix
    static const int numGPsMassMatrixTri;
    /// number of Gauss points used for computing quad element mass matrix
    static const int numGPsMassMatrixQuad;

};

} /* namespace EMPIRE */
#endif /* DATAFIELDINTEGRATION_H_ */
