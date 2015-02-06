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
 * \file CurveSurfaceMapper.h
 * This file holds the class CurveSurfaceMapper
 * \date 2/6/2015
 **************************************************************************************************/
#ifndef ABSTRACTCURVESURFACEMAPPER_H_
#define ABSTRACTCURVESURFACEMAPPER_H_

#include "AbstractMapper.h"
#include "EMPEROR_Enum.h"
#include <map>

namespace EMPIRE {

class KinematicMotion;

/********//**
 * \brief Class CurveSurfaceMapper maps DOFs between beam elements and surface mesh with sections
 * \author Tianyang Wang
 ***********/
class CurveSurfaceMapper: public AbstractMapper {
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
    CurveSurfaceMapper(EMPIRE_CurveSurfaceMapper_type _type, int _curveNumNodes,
            int _curveNumElements, const double *_curveNodeCoors, const int *_curveNodeIDs,
            const int *_curveElems, int _surfaceNumNodes, const double *_surfaceNodeCoors,
            int _surfaceNumSections, int _surfaceNumRootSectionNodes,
            int _surfaceNumNormalSectionNodes, int _surfaceNumTipSectionNodes,
            const double *rotation_O_Q, const double *translation_O_Q);
    /***********************************************************************************************
     * \brief Destructor
     * \author Tianyang Wang
     ***********/
    virtual ~CurveSurfaceMapper();
    /***********************************************************************************************
     * \brief Map deformation from curve to surface / reconstruct the deformed surface according to beam DOFs.
     * \param[in] curveDispRot displacements and rotations on curve nodes
     * \param[out] surfaceDisp displacements on surface nodes
     * \author Tianyang Wang
     ***********/
    void consistentMapping(const double *curveDispRot, double *surfaceDisp);
    /***********************************************************************************************
     * \brief Map deformation from curve to surface with the linear Bernoulli beam theory
     *        It allows only small displacements.
     * \param[in] curveDispRot displacements and rotations on curve nodes
     * \param[out] surfaceDisp displacements on surface nodes
     * \author Tianyang Wang
     ***********/
    void consistentMappingLinear(const double *curveDispRot, double *surfaceDisp);
    /***********************************************************************************************
     * \brief Map deformation from curve to surface with 2D corotating algorithm.
     *        Corotate2D allows large translation/rotation in 2D and small value for the rest DOFs.
     *        See Non-linear Modeling and Analysis of Solids and Structures (Krenk2009).
     * \param[in] curveDispRot displacements and rotations on curve nodes
     * \param[out] surfaceDisp displacements on surface nodes
     * \author Tianyang Wang
     ***********/
    void consistentMappingCorotate2D(const double *curveDispRot, double *surfaceDisp);
    /***********************************************************************************************
     * \brief Map deformation from curve to surface with 3D corotating algorithm
     *        See Non-linear Modeling and Analysis of Solids and Structures (Krenk2009).
     * \param[in] curveDispRot displacements and rotations on curve nodes
     * \param[out] surfaceDisp displacements on surface nodes
     * \author Tianyang Wang
     ***********/
    void consistentMappingCorotate3D(const double *curveDispRot, double *surfaceDisp);
    /***********************************************************************************************
     * \brief Map loads form surface to curve
     * \param[in] surfaceForce forces on surface nodes
     * \param[out] curveForceMoment forces and moments on curve nodes
     * \author Tianyang Wang
     * ***********/
    void conservativeMapping(const double *surfaceForce, double *curveForceMoment);

protected:
    /// type of the CurveSurfaceMapper
    EMPIRE_CurveSurfaceMapper_type type;
    /// number of nodes of a curve/beam
    int curveNumNodes;
    /// number of elements of a curve/beam
    int curveNumElements;
    /// coordinates of the curve points
    const double *curveNodeCoors;
    //const int *curveNodeIDs;
    /// connectivity table of curve/beam elements
    const int *curveElems;
    /// number of surface nodes
    int surfaceNumNodes;
    /// node coordinates of surface nodes
    const double *surfaceNodeCoors;
    /// number of sections of the surface
    int surfaceNumSections;
    /// number of root section nodes of the surface
    int surfaceNumRootSectionNodes;
    /// number of normal section nodes of the surface
    int surfaceNumNormalSectionNodes;
    /// number of the tip section nodes of the surface
    int surfaceNumTipSectionNodes;
    /// After sorting the surface nodes into sections, map from sorted position to its position in node array
    int *sortedPosToUnsortedPos;
    /// Map a curve/beam node ID to the its position in node array
    std::map<int, int> *curveNodeIDToPos;
    /// Map a section to the corresponding curve/beam element
    int *sectionToCurveElem;
    /// cross point between section and curve/beam
    double *sectionP;
    /// Shape functions and their derivatives of a section in the corresponding curve/beam element
    double *shapeFuncOfSection;
    /// Rotation from the global system the curve/beam element local system
    KinematicMotion **ROT_O_ELEM;
    /// rotation of a section needed by conservative mapping
    double *sectionRot;
    /***********************************************************************************************
     * \brief Normalize a rotation (or length) vector and return the rotation angle (or length)
     * \param[in] vector a rotation (or length) vector
     * \return length or rotation angle
     * \author Tianyang Wang
     * ***********/
    double normalizeVector(double *vector);

    /// unit test class
    friend class TestCurveSurfaceMapper;

    /// length of curve elements
    double *curveElemLength; // linear

    /// Kinematic motion from the global system O to the root section system Q
    KinematicMotion *KM_O_Q; // corotate2D
    /// Nodal coordinates of curve in root section system Q
    double *curveNodeCoorsInQ; // corotate2D
    /// Rotation from the root section system Q to the curve/beam element local system
    KinematicMotion **ROT_Q_ELEM; // corotate2D
    /// Rotation angle the root section system Q to the curve/beam element local system
    double *angle_Q_ELEM; // corotate2D
};

} /* namespace EMPIRE */

#endif /* ABSTRACTCURVESURFACEMAPPER_H_ */
