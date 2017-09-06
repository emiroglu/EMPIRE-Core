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
 * \file MapperAdapter.h
 * This file holds the class MapperAdapter
 * \date 3/5/2012
 **************************************************************************************************/

#ifndef MAPPERADAPTER_H_
#define MAPPERADAPTER_H_

#include <string>
#include "EMPEROR_Enum.h"

namespace EMPIRE {

class MortarMapper;
class AbstractMesh;
class DataField;
class AbstractMapper;

/********//**
 * \brief Class MapperAdapter is the adaptor of the mapper.
 ***********/
class MapperAdapter {
public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _name name of the mapper
     * \param[in] _meshA meshA
     * \param[in] _meshB meshB
     * \author Tianyang Wang
     ***********/
    MapperAdapter(std::string _name, AbstractMesh *_meshA, AbstractMesh *_meshB);
    /***********************************************************************************************
     * \brief Initialize MortarMapper
     * \param[in] oppositeSurfaceNormal whether the surface normal of A and B are opposite
     * \param[in] dual use dual mortar or not
     * \param[in] enforceConsistency enforce consistensy in consistent mapping
     * \author Tianyang Wang
     ***********/
    void initMortarMapper(bool oppositeSurfaceNormal, bool dual, bool enforceConsistency);
    /***********************************************************************************************
     * \brief Initialize IGA Mortar Mapper
     * \param[in] _enforceConsistency A flag on whether the consistency should be enforced
     * \param[in] _tolConsistency The tolerance up to which the consistency in the mapping is assumed
     * \param[in] _maxProjectionDistance The max distance allowed between FE mesh and NURBS surface
     * \param[in] _noInitialGuess The number of test point to find initial guess for Newton-Raphson scheme
     * \param[in] _maxProjectionDistanceOnDifferentPatches The max authorized distance between two projected points from a same physical node
     * \param[in] _noIterationsNewton The number of iteration for Newton-Raphson scheme of projecting a node on a NURBS patch
     * \param[in] _tolProjectionNewtonRaphson The tolerance for Newton-Raphson scheme of projecting a node on a NURBS patch
     * \param[in] _noIterationsNewtonRaphsonBoundary The number of iteration for Newton-Raphson scheme of projecting a node on a NURBS patch boundary
     * \param[in] _tolProjectionNewtonRaphsonBoundary The tolerance for Newton-Raphson scheme of projecting a node on a NURBS patch boundary
     * \param[in] _noIterationsBisection The number of iteration for bisection scheme of projecting a node on a NURBS patch boundary
     * \param[in] _tolProjectionBisection The tolerance for bisection scheme of projecting a node on a NURBS patch boundary
     * \param[in] _noGPTriangle The number of Gauss points when performs integration on triangle
     * \param[in] _noGPQuad The number of Gauss points when performs integration on quadrilateral
     * \param[in] _isWeakCurveDirichletConditions Flag on whether weak Dirichlet conditions along trimming curves are assumed
     * \param[in] _isAutomaticPenaltyParametersWeakCurveDirichletConditions Flag on whether the Penalty parameters for the application of weak Dirichlet conditions along trimming curves is automatic
     * \param[in] _isPrimPrescribedWeakCurveDirichletConditions Flag on whether the primary field is prescribed along the trimming curves
     * \param[in] _isSecBendingPrescribedWeakCurveDirichletConditions Flag on whether the bending rotation of the primary field is prescribed along the trimming curves
     * \param[in] _isSecTwistingPrescribedWeakCurveDirichletConditions Flag on whether the twisting rotation of the primary field is prescribed along the trimming curves
     * \param[in] _alphaPrimWeakCurveDirichletConditions Penalty parameter for the application of Dirichlet condition along trimming curves based on the primary field condition
     * \param[in] _alphaSecBendingWeakCurveDirichletConditions Penalty parameter for the application of Dirichlet condition along trimming curves based on the bending rotation of the primary field condition
     * \param[in] _alphaSecTwistingWeakCurveDirichletConditions Penalty parameter for the application of Dirichlet condition along trimming curves based on the twisting rotation of the primary field condition
     * \param[in] _isWeakSurfaceDirichletConditions Flag on whether weak Dirichlet conditions across a surface are assumed
     * \param[in] _isPrimPrescribedWeakSurfaceDirichletConditions Flag on whether the primary field is prescribed along a surface
     * \param[in] _isAutomaticPenaltyParametersWeakSurfaceDirichletConditions Flag on whether the Penalty parameters for the application of weak Dirichlet conditions across a surface is automatic
     * \param[in] _alphaPrimWeakCurveDirichletConditions Penalty parameter for the application of Dirichlet condition across a surface based on the primary field condition
     * \param[in] _isWeakPatchContinuityConditions Flag on whether across patch continuity conditions are assumed
     * \param[in] _isPrimCoupledWeakContinuityConditions Flag on whether the primary field is coupled between multipatches
     * \param[in] _isSecBendingCoupledWeakContinuityConditions Flag on whether the bending rotation of the primary field is coupled between multipatches
     * \param[in] _isSecTwistingCoupledWeakContinuityConditions Flag on whether the twisting rotation of the primary field is coupled between multipatches
     * \param[in] _isAutomaticPenaltyParametersWeakContinuityConditions Flag on whether the Penalty parameters for the weak enforcement of the continuity conditions across the patches are automatically computed
     * \param[in] _alphaPrimWeakContinuityConditions Penalty parameter for the continuity of the field across the patches
     * \param[in] _alphaSecBendingWeakContinuityConditions Penalty parameter for the continuity of the bending rotation of the field across the patches
     * \param[in] _alphaSecTwistingWeakContinuityConditions Penalty parameter for the continuity of the twisting rotation of the field across the patches
     * \param[in] _isStrongCurveDirichletConditions Flag on whether strong Dirichlet conditions are applied
     * \param[in] _isErrorComputation Flag on whether error computation is assumed
     * \param[in] _isDomainError Flag on the computation of the domain error from the mapping
     * \param[in] _isInterfaceError Flag on the computation of the interface error from the mapping
     * \param[in] _isCurveError Flag on the computation of the error along the trimming curves where Dirichlet conditions are applied
     * \author Andreas Apostolatos, Altug Emiroglu, Fabien Pean
     ***********/
    void initIGAMortarMapper(bool _enforceConsistency, double _tolConsistency,
                             double _maxProjectionDistance, int _noInitialGuess, double _maxProjectionDistanceOnDifferentPatches,
                             int _noIterationsNewton, double _tolProjectionNewtonRaphson,
                             int _noIterationsNewtonRaphsonBoundary, double _tolProjectionNewtonRaphsonBoundary,
                             int _noIterationsBisection, double _tolProjectionBisection,
                             int _noGPTriangle, int _noGPQuad,
                             bool _isWeakCurveDirichletConditions, bool _isAutomaticPenaltyParametersWeakCurveDirichletConditions, bool _isPrimPrescribedWeakCurveDirichletConditions, bool _isSecBendingPrescribedWeakCurveDirichletConditions, bool _isSecTwistingPrescribedWeakCurveDirichletConditions, double _alphaPrimWeakCurveDirichletConditions, double _alphaSecBendingWeakCurveDirichletConditions, double _alphaSecTwistingWeakCurveDirichletConditions,
                             bool _isWeakSurfaceDirichletConditions, bool _isAutomaticPenaltyParametersWeakSurfaceDirichletConditions, bool _isPrimPrescribedWeakSurfaceDirichletConditions, double _alphaPrimWeakSurfaceDirichletConditions,
                             bool _isWeakPatchContinuityConditions, bool _isAutomaticPenaltyParametersWeakContinuityConditions, bool _isPrimCoupledWeakContinuityConditions, bool _isSecBendingCoupledWeakContinuityConditions, bool _isSecTwistingCoupledWeakContinuityConditions, double _alphaPrimWeakContinuityConditions, double _alphaSecBendingWeakContinuityConditions, double _alphaSecTwistingWeakContinuityConditions,
                             bool _isStrongCurveDirichletConditions,
                             bool _isErrorComputation, bool _isDomainError, bool _isCurveError, bool _isInterfaceError);

    /***********************************************************************************************
     * \brief Initialize IGA Barycentric Mapper
     * \param[in] _maxProjectionDistance The max distance allowed between FE mesh and NURBS surface
     * \param[in] _numRefinementForIntialGuess The number of test point to find initial guess for Newton-Raphson scheme
     * \param[in] _maxDistanceForProjectedPointsOnDifferentPatches The max authorized distance between two projected points from a same physical node
     * \param[in] _newtonRaphsonMaxIt The number of iteration for Newton-Raphson scheme of projecting a node on a NURBS patch
     * \param[in] _newtonRaphsonTol The tolerance for Newton-Raphson scheme of projecting a node on a NURBS patch
     * \author Apostolos Petalas, Fabien Pean, Chenshen Wu
     ***********/
    void initIGABarycentricMapper(double _maxProjectionDistance, int _numRefinementForIntialGuess, double _maxDistanceForProjectedPointsOnDifferentPatches,
                             int _newtonRaphsonMaxIt, double _newtonRaphsonTol);
    /***********************************************************************************************
     * \brief Initialize NearestNeighborMapper
     * \author Tianyang Wang
     ***********/
    void initNearestNeighborMapper();
    /***********************************************************************************************
     * \brief Initialize BarycentricInterpolationMapper
     * \author Tianyang Wang
     ***********/
    void initBarycentricInterpolationMapper();
    /***********************************************************************************************
     * \brief Initialize NearestElementMapper
     * \author Tianyang Wang
     ***********/
    void initNearestElementMapper();
    /***********************************************************************************************
     * \brief Initialize CurveSurfaceMapper
     * \param[in] type type of the CurveSurfaceMapper
     * \author Tianyang Wang
     ***********/
    void initCurveSurfaceMapper(EMPIRE_CurveSurfaceMapper_type type);
    /***********************************************************************************************
     * \brief Destructor
     * \author Tianyang Wang
     ***********/
    virtual ~MapperAdapter();
    /***********************************************************************************************
     * \brief Do consistent mapping from A to B (map displacements)
     * \param[in] fieldA is the input data
     * \param[out] fieldB is the output data
     * \author Tianyang Wang
     ***********/
    void consistentMapping(const DataField *fieldA, DataField *fieldB);
    /***********************************************************************************************
     * \brief Do conservative mapping from B to A (map forces)
     * \param[in] fieldB is the input data
     * \param[out] fieldA is the output data
     * \author Tianyang Wang
     ***********/
    void conservativeMapping(const DataField *fieldB, DataField *fieldA);
    /***********************************************************************************************
     * \brief is it meshA or not
     * \param[in] mesh mesh
     * \return ture if it is meshA
     * \author Tianyang Wang
     ***********/
    bool isMeshA(AbstractMesh *mesh) {
        return meshA == mesh;
    }
    /***********************************************************************************************
     * \brief is it meshB or not
     * \param[in] mesh mesh
     * \return ture if it is meshB
     * \author Tianyang Wang
     ***********/
    bool isMeshB(AbstractMesh *mesh) {
        return meshB == mesh;
    }
private:
    /// the adapted mapper
    AbstractMapper *mapperImpl;
    /// name of the mapper
    std::string name;
    /// mesh A
    AbstractMesh *meshA;
    /// mesh B
    AbstractMesh *meshB;
};

} /* namespace EMPIRE */
#endif /* MAPPERADAPTER_H_ */
