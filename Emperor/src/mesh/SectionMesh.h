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
 * \file SectionMesh.h
 * This file holds the class DataField
 * \date 2/10/2015
 **************************************************************************************************/

#ifndef SECTIONMESH_H_
#define SECTIONMESH_H_

#include "FEMesh.h"
#include <string>

namespace EMPIRE {
class Message;

/********//**
 * \brief Class SectionMesh stores an FE mesh having parallel sections in the beam length direction
 ***********/
class SectionMesh : public FEMesh {
public:
    /***********************************************************************************************
     * \brief Constructor, allocating the storage of the mesh
     * \param[in] _name name of the mesh
     * \param[in] _numNodes number of nodes
     * \param[in] _numElems number of elements
     * \param[in] _triangulateAll triangulate all elements
     * \author Tianyang Wang
     ***********/
    SectionMesh(std::string _name, int _numNodes, int _numElems, bool _triangulateAll = false);
    /***********************************************************************************************
     * \brief Destructor
     * \author Tianyang Wang
     ***********/
    virtual ~SectionMesh();
    /***********************************************************************************************
     * \brief get number of sections
     * \return number of sections
     * \author Tianyang Wang
     ***********/
    int getNumSections();
    /***********************************************************************************************
     * \brief get number of nodes of the root section
     * \return number of nodes of the root section
     * \author Tianyang Wang
     ***********/
    int getNumRootSectionNodes();
    /***********************************************************************************************
     * \brief get number of nodes of the normal sections
     * \return number of nodes of the normal sections
     * \author Tianyang Wang
     ***********/
    int getNumNormalSectionNodes();
    /***********************************************************************************************
     * \brief get number of nodes of the tip section
     * \return number of nodes of the tip section
     * \author Tianyang Wang
     ***********/
    int getNumTipSectionNodes();
    /***********************************************************************************************
     * \brief get the translation from the global system to the root section system
     * \return a translation vector
     * \author Tianyang Wang
     ***********/
    const double *getTranslationGlobal2Root();
    /***********************************************************************************************
     * \brief get the rotation from the global system to the root section system
     * \return a rotation matrix
     * \author Tianyang Wang
     ***********/
    const double *getRotationGlobal2Root();
    /***********************************************************************************************
     * \brief set number of sections
     * \param[in] number of sections
     * \author Tianyang Wang
     ***********/
    void setNumSections(int _numSections);
    /***********************************************************************************************
     * \brief set number of nodes of the root section
     * \param[in] number of nodes of the root section
     * \author Tianyang Wang
     ***********/
    void setNumRootSectionNodes(int _numRootSectionNodes);
    /***********************************************************************************************
     * \brief set number of nodes of the normal sections
     * \param[in] number of nodes of the normal sections
     * \author Tianyang Wang
     ***********/
    void setNumNormalSectionNodes(int _numNormalSectionNodes);
    /***********************************************************************************************
     * \brief set number of nodes of the tip section
     * \param[in] number of nodes of the tip section
     * \author Tianyang Wang
     ***********/
    void setNumTipSectionNodes(int _numTipSectionNodes);
    /***********************************************************************************************
     * \brief set the translation from the global system to the root section system
     * \param[in] a translation vector
     * \author Tianyang Wang
     ***********/
    void setTranslationGlobal2Root(double *_translationGlobal2Root);
    /***********************************************************************************************
     * \brief set the rotation from the global system to the root section system
     * \param[in] a rotation matrix
     * \author Tianyang Wang
     ***********/
    void setRotationGlobal2Root(double *_rotationGlobal2Root);

private:
    /// number of sections
    int numSections;
    /// number of nodes of the root section
    int numRootSectionNodes;
    /// number of nodes of the normal sections
    int numNormalSectionNodes;
    /// number of nodes of the tip sections
    int numTipSectionNodes;
    /// the rotation from the global system to the root section system
    double *rotationGlobal2Root;
    /// the translation from the global system to the root section system
    double *translationGlobal2Root;
};
/***********************************************************************************************
 * \brief Allows for nice debug output later
 * \author Stefan Sicklinger
 ***********/
Message &operator<<(Message &message, SectionMesh &mesh);

} /* namespace EMPIRE */

#endif /* SECTIONMESH_H_ */
