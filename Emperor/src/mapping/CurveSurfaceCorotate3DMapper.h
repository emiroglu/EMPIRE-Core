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
 * \file CurveSurfaceCorotate3DMapper.h
 * This file holds the class CurveSurfaceCorotate3DMapper
 * \date 2/6/2015
 **************************************************************************************************/
#ifndef CURVESURFACECOROTATE3DMAPPER_H_
#define CURVESURFACECOROTATE3DMAPPER_H_

#include "AbstractCurveSurfaceMapper.h"
#include <map>

namespace EMPIRE {

class KinematicMotion;
/********//**
 * \brief Class CurveSurfaceCorotate3DMapper performs 3D corotate algorithm to reconstruct the deformed
 *        surface of a deformed curve. See Non-linear Modeling and Analysis of Solids and Structures (Krenk2009).
 * \author Tianyang Wang
 ***********/
class CurveSurfaceCorotate3DMapper: public AbstractCurveSurfaceMapper {
public:
    /***********************************************************************************************
     * \brief Constructor. Builds the relation between curve elements and surface sections
     * \param[in] _curveNumNodes number of nodes of curve
     * \param[in] _curveNumElements number of elements of curve
     * \param[in] _curveNodeCoors nodal coordinates of curve
     * \param[in] _curveNodeIDs nodal IDs of curve
     * \param[in] _curveElems connectivity table of curve elements
     * \param[in] _surfaceNumNodes number of nodes of surface
     * \param[in] _surfaceNodeCoors nodal coordinates of surface
     * \param[in] _surfaceNumSections number of sections of surface
     * \param[in] _surfaceNumRootSectionNodes number of nodes of surface root section
     * \param[in] _surfaceNumNormalSectionNodes number of nodes of a surface normal section
     * \param[in] _surfaceNumTipSectionNodes number of nodes of surface tip section
     * \param[in] rotation_O_Q rotation from global system O to beam root system Q
     * \param[in] translation_O_Q translation from global system O to beam root system Q
     * \author Tianyang Wang
     ***********/
    CurveSurfaceCorotate3DMapper(int _curveNumNodes, int _curveNumElements,
            const double *_curveNodeCoors, const int *_curveNodeIDs, const int *_curveElems,
            int _surfaceNumNodes, const double *_surfaceNodeCoors, int _surfaceNumSections,
            int _surfaceNumRootSectionNodes, int _surfaceNumNormalSectionNodes,
            int _surfaceNumTipSectionNodes, const double *rotation_O_Q,
            const double *translation_O_Q);
    /***********************************************************************************************
     * \brief Destructor
     * \author Tianyang Wang
     ***********/
    virtual ~CurveSurfaceCorotate3DMapper();
    /***********************************************************************************************
     * \brief Map deformation from curve to surface
     * \param[in] curveDispRot displacements and rotations on curve nodes
     * \param[out] surfaceDisp displacements on surface nodes
     * \author Tianyang Wang
     ***********/
    void consistentMapping(const double *curveDispRot, double *surfaceDisp);
};

} /* namespace EMPIRE */

#endif /* CURVESURFACECOROTATE3DMAPPER_H_ */
