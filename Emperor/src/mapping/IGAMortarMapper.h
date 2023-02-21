/*  Copyright &copy; 2013, TU Muenchen, Chair of Structural Analysis,
 *  Fabien Pean, Andreas Apostolatos
 *  Ragnar Björnsson, Stefan Sicklinger, Tianyang Wang, Munich
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
/******************************************************************************//**
 * \file IGAMortarMapper.h
 * The header file of class IGAMortarMapper.
 * \date 3/6/2013
 *********************************************************************************/

#ifndef IGAMORTARMAPPER_H_
#define IGAMORTARMAPPER_H_

#include <map>
#include <set>
#include <vector>
#include <iostream>
#include <assert.h>
#include "AbstractMapper.h"
#include "IGAMortarCouplingMatrices.h"

namespace flann {
template<typename Distance> class Index;
template<class T> struct L2;
template<typename T> class Matrix;
}

namespace EMPIRE {

namespace MathLibrary {
template<class T> class SparseMatrix;
}

class IGAPatchSurface;
class IGAPatchSurfaceTrimmingLoop;
class IGAMesh;
class Message;
class FEMesh;
class DataField;
class IGAMortarCouplingMatrices;

namespace MathLibrary {

class IGAGaussQuadratureOnTriangle;
class IGAGaussQuadratureOnQuad;

}

/***********************************************************************************************
 * \brief This class is computing coupling matrices to perform the Mortar method between a NURBS geometry and a polygon mesh
 * ***********/
class IGAMortarMapper: public AbstractMapper {
private:
    /// Type definitions
    typedef std::pair<double,double> Point2D;
    typedef std::vector<Point2D> Polygon2D;
    typedef std::vector<Polygon2D> ListPolygon2D;
private:
    /// Name of the mapper
    std::string name;

    /// IGA Mesh
    IGAMesh *meshIGA;

    /// Fluid Mesh
    FEMesh *meshFE;

    /// Flag on whether the meshFEDirectElemTable was created
    bool isMeshFEDirectElemTable;

    /// The element freedom table for the fluid mesh
    int **meshFEDirectElemTable;

    /// The reverse element freedom table for the finite element mesh
    std::map<int, std::vector<int> > meshFENodeToElementTable;

    /// Indices for rows which are identically zero in the mass matrix
    std::vector<int> indexEmptyRowCnn;

    /// Flag on whether the coupling matrices are initialized
    bool isCouplingMatrices;

    /// The isogeometric coupling matrices
    IGAMortarCouplingMatrices *couplingMatrices;

    /// Number of weak IGA Dirichlet curve conditions
    int noWeakIGADirichletCurveConditions;

    /// Penalty factors for the primary field to the application of weak Dirichlet curve conditions
    double* weakDirichletCCAlphaPrimary;

    /// Penalty factors for the bending secondary field to the application of weak Dirichlet curve conditions
    double* weakDirichletCCAlphaSecondaryBending;

    /// Penalty factors for the twisting secondary field to the application of weak Dirichlet curve conditions
    double* weakDirichletCCAlphaSecondaryTwisting;

    /// Number of weak IGA Dirichlet surface conditions
    int noWeakIGADirichletSurfaceConditions;

    /// Penalty factors for the primary field to the application of weak Dirichlet surface conditions
    double* weakDirichletSCAlphaPrimary;

    /// Penalty factors for the bending secondary field to the application of weak Dirichlet surface conditions
    double* weakDirichletSCAlphaSecondaryBending;

    /// Penalty factors for the twisting secondary field to the application of weak Dirichlet surface conditions
    double* weakDirichletSCAlphaSecondaryTwisting;

    /// Number of weak IGA patch continuity conditions
    int noWeakIGAPatchContinuityConditions;

    /// Penalty factors for the primary field to the application of weak patch continuity conditions
    double* weakPatchContinuityAlphaPrimaryIJ;

    /// Penalty factors for the bending secondary field to the application of weak patch continuity conditions
    double* weakPatchContinuityAlphaSecondaryBendingIJ;

    /// Penalty factors for the twisting secondary field to the application of weak patch continuity conditions
    double* weakPatchContinuityAlphaSecondaryTwistingIJ;

    /// Integration area
    double areaIntegration;

    /// Minimum element area in the multipatch geometry
    double minElArea;

    /// Minimum edge size of the Finite Element mesh
    double minEdgeSize;

    /// Minimum edge size at the Dirichlet boundary
    double minElEdgeSizeDirichlet;

    /// Minimum edge size at the interface between the patches in the multipatch geometry
    double minElEdgeSizeInterface;

    /// Flag on whether the expanded version of the coupling matrices is computed
    bool isExpanded;

    /// Flag on whether the gauss rules have been created
    bool isGaussQuadature;

    /// Quadrature rule over the triangulated subdomains
    EMPIRE::MathLibrary::IGAGaussQuadrature **gaussRuleOnTriangle;

    /// Quadrature rule over the non-triangulated subdomains
    EMPIRE::MathLibrary::IGAGaussQuadrature **gaussRuleOnQuadrilateral;

    /// The parametric coordinates of the projected nodes on the surface
    /// For each node i, for each possible patch j, store parametric coordinates of i in j
    std::vector<std::map<int, std::vector<double> > > projectedCoords;

    /// Polygon reconstructed in 2D parametric space stored for each patch
    std::map<int,ListPolygon2D> trimmedProjectedPolygons;

    /// Triangulated polygon in 2D parametric space stored for each patch
    std::map<int,ListPolygon2D> triangulatedProjectedPolygons2;

    /// List of all the projected polygons
    std::vector<std::map<int,Polygon2D> > projectedPolygons;

    /// List of all the triangulated polygons
    std::vector<std::map<int,ListPolygon2D> > triangulatedProjectedPolygons;

    /// Stream of gauss points stored in line with format
    /// Weight / Jacobian / NumOfFENode / Node1 / ShapeValue1 / Node2 / ShapeValue2 ... NumOfIGANode / Node1 / ShapeValue1/ ... cartesianCoordinatesGP
    std::vector<std::vector<double> > streamGPs;

    /// Stream of interface gauss points stored in line with format
    std::vector<std::vector<double> > streamInterfaceGPs;

    /// Stream of gauss points for the trimming curves where conditions are applied stored in line with format
    std::vector<std::vector<double> > streamCurveGPs;

    /// Flag on the mapping direction
    bool isMappingIGA2FEM;

    /// Number of nodes for the slave side
    size_t numNodesSlave;

    /// Number of nodes for the master side
    size_t numNodesMaster;

    /// Consistency properties
    struct propConsistency {
        bool enforceConsistency;
        bool tolConsistency;
    } propConsistency;

    /// Properties for the projection schemes
    struct propProjection {
        double maxProjectionDistance;
        int noInitialGuess;
        double maxProjectionDistanceOnDifferentPatches;
    } propProjection;

    /// Properties for the nonlinear solution schemes
    struct propNonlinearSchemes {
        int noIterations;
        double tolProjection;
    } propNewtonRaphson, propNewtonRaphsonBoundary, propBisection;

    /// Properties for the integration
    struct propIntegration {
        int isAutomaticNoGPTriangle;
        int noGPTriangle;
        int isAutomaticNoGPQuadrilateral;
        int noGPQuadrilateral;
    } propIntegration;

    /// Properties for the application of weak Dirichlet conditions along trimming curves
    struct propWeakCurveDirichletConditions {
        bool isWeakCurveDirichletConditions;
        bool isAutomaticPenaltyParameters;
        bool isPrimPrescribed;
        bool isSecBendingPrescribed;
        bool isSecTwistingPrescribed;
        double alphaPrim;
        double alphaSecBending;
        double alphaSecTwisting;
    } propWeakCurveDirichletConditions;

    /// Properties for the application of weak Dirichlet conditions across surfaces
    struct propWeakSurfaceDirichletConditions {
        bool isWeakSurfaceDirichletConditions;
        bool isAutomaticPenaltyParameters;
        bool isPrimPrescribed;
        bool alphaPrim;
    } propWeakSurfaceDirichletConditions;

    /// Properties for the application of weak continuity conditions across the multipatches
    struct propWeakPatchContinuityConditions {
        bool isWeakPatchContinuityConditions;
        bool isAutomaticPenaltyParameters;
        bool isPrimCoupled;
        bool isSecBendingCoupled;
        bool isSecTwistingCoupled;
        double alphaPrim;
        double alphaSecBending;
        double alphaSecTwisting;
    } propWeakPatchContinuityConditions;

    /// On the strong application of Dirichlet boundary conditions
    struct propStrongCurveDirichletConditions {
        int isStrongCurveDirichletConditions;
    } propStrongCurveDirichletConditions;

    /// On the error computation
    struct propErrorComputation {
        bool isErrorComputation;
        bool isDomainError;
        bool isCurveError;
        bool isInterfaceError;
    } propErrorComputation;

public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _name The name of the mapper
     * \param[in] _meshA Mesh A (either IGA or FEM mesh)
     * \param[in] _meshB Mesh B (either IGA or FEM mesh)
     * \param[in] _isExpanded Flag on whether or not the coupling matrices will be expanded
     * \param[in] _isMappingIGA2FEM
     * \author Andreas Apostolatos, Altug Emiroglu, Fabien Pean
     ***********/
    IGAMortarMapper(std::string _name, AbstractMesh *_meshA, AbstractMesh *_meshB);

    /***********************************************************************************************
     * \brief Destructor
     * \author Andreas Apostolatos, Fabien Pean
     ***********/
    virtual ~IGAMortarMapper();

    /***********************************************************************************************
     * \brief Initialization of the mapper
     * \author Andreas Apostolatos
     ***********/
    void initialize();

    /***********************************************************************************************
     * \brief Set the the parameters for the consistency enforcement
     * \param[in] _enforceConsistency The consistency flag
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    void setParametersConsistency(bool _enforceConsistency = false, double _tolConsistency = 1e-6);

    /***********************************************************************************************
     * \brief Set the parameters for the projection of mesh onto the NURBS surface
     * \param[in] _maxProjectionDistance The max distance allowed between FE mesh and NURBS surface
     * \param[in] _noInitialGuess The number of test point to find initial guess for Newton-Raphson scheme
     * \param[in] _maxProjectionDistanceOnDifferentPatches The max authorized distance between two projected points from a same physical node
     * \author Andreas Apostolatos, Fabien Pean
     ***********/
    void setParametersProjection(double _maxProjectionDistance = 1e-2, int _noInitialGuess = 10,
                                 double _maxProjectionDistanceOnDifferentPatches = 1e-3);

    /***********************************************************************************************
     * \brief Set the parameters for Newton-Raphson scheme of projection on NURBS patch
     * \param[in] _noIterations The number of iteration for Newton-Raphson scheme of projecting a node on a NURBS patch
     * \param[in] _tolProjection The tolerance for Newton-Raphson scheme of projecting a node on a NURBS patch
     * \author Andreas Apostolatos, Fabien Pean
     ***********/
    void setParametersNewtonRaphson(int _noIterations = 20, double _tolProjection = 1e-6);

    /***********************************************************************************************
     * \brief Set the parameters for Newton-Raphson scheme of projection on NURBS patch boundary
     * \param[in] _noIterations The number of iteration for Newton-Raphson scheme of projecting a node on a NURBS patch boundary
     * \param[in] _tolProjection The tolerance for Newton-Raphson scheme of projecting a node on a NURBS patch boundary
     * \author Andreas Apostolatos, Fabien Pean
     ***********/
    void setParametersNewtonRaphsonBoundary(int _noIterations = 20, double _tolProjection = 1e-6);

    /***********************************************************************************************
     * \brief Set the parameters for bisection scheme of projection on NURBS patch boundary
     * \param[in] _bisectionMaxIt The number of iteration for bisection scheme of projecting a node on a NURBS patch boundary
     * \param[in] _bisectionTol The tolerance for bisection scheme of projecting a node on a NURBS patch boundary
     * \author Andreas Apostolatos, Fabien Pean
     ***********/
    void setParametersBisection(int _noIterations = 40, double _tolProjection = 1e-6);

    /***********************************************************************************************
     * \brief Set the parameters for integration
     * \param[in] _isAutomaticNoGPTriangle Flag on whether the Gauss rule on a triangle is automatically created or not
     * \param[in] _noGPTriangle The number of Gauss points when performs integration on triangle
     * \param[in] _isAutomaticNoGPQuadrilateral Flag on whether the Gauss rule on a triangle is automatically created or not
     * \param[in] _noGPQuadrilateral The number of Gauss points when performs integration on quadrilateral
     * \author Andreas Apostolatos, Fabien Pean
     ***********/
    void setParametersIntegration(bool _isAutomaticNoGPTriangle = false, int _noGPTriangle = 16, bool _isAutomaticNoGPQuadrilateral = false, int _noGPQuadrilateral = 25);

    /***********************************************************************************************
     * \brief Set the parameters for the application of weak Dirichlet boundary conditions along trimming curves
     * \param[in] _isWeakCurveDirichletConditions Flag on the application of weak Dirichlet conditions along trimming curves
     * \param[in] _isAutomaticPenaltyParameters Flag on whether the Penalty parameters are automatically computed
     * \param[in] _isPrimPrescribed Flag on whether the primary field is prescribed along the trimming curves
     * \param[in] _isSecBendingPrescribed Flag on whether the bending rotation of the primary field is prescribed along the trimming curves
     * \param[in] _isSecTwistingPrescribed Flag on whether the twisting rotation of the primary field is prescribed along the trimming curves
     * \param[in] _alphaPrim Penalty parameter for the primary field
     * \param[in] _alphaSecBending Penalty parameter for the bending rotation of the primary field
     * \param[in] _alphaSecTwisting Penalty parameter for the twisting rotation of the primary field
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    void setParametersWeakCurveDirichletConditions(bool _isWeakCurveDirichletConditions = false, bool _isAutomaticPenaltyParameters = false,
                                                   bool _isPrimPrescribed = false, bool _isSecBendingPrescribed = false, bool _isSecTwistingPrescribed = false,
                                                   double _alphaPrim = 0.0, double _alphaSecBending = 0.0, double _alphaSecTwisting = 0.0);

    /***********************************************************************************************
     * \brief Set the parameters for the application of weak Dirichlet boundary conditions across surfaces
     * \param[in] _isWeakSurfaceDirichletConditions Flag on the application of weak Dirichlet conditions across surfaces
     * \param[in] _isAutomaticPenaltyFactors Flag on whether the Penalty parameters are automatically computed
     * \param[in] _isPrimPrescribed Flag on whether the primary field is prescribed across surfaces
     * \param[in] _alphaPrim Penalty parameter for the primary field
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    void setParametersWeakSurfaceDirichletConditions(bool _isWeakSurfaceDirichletConditions = false, bool _isAutomaticPenaltyParameters = false,
                                                     bool _isPrimPrescribed = false, double _alphaPrim = 0.0);

    /***********************************************************************************************
     * \brief Set parameter for the application of weak continuity conditions across the multipatches
     * \param[in] _isWeakCurveDirichletConditions Flag on the application of weak continuity conditions between the multipatches
     * \param[in] _isAutomaticPenaltyFactors Flag on whether the Penalty parameters are automatically computed
     * \param[in] _isPrimCoupled Flag on whether the primary fields are coupled between the multipatches
     * \param[in] _isSecBendingCoupled Flag on whether the bending rotation of the primary fields is coupled between the multipatches
     * \param[in] _isSecTwistingCoupled Flag on whether the twisting roation of the primary fields is coupled between the multipatches
     * \param[in] _alphaPrim Penalty parameter for the primary field
     * \param[in] _alphaSecBending Penalty parameter for the bending rotation of the primary field
     * \param[in] _alphaSecTwisting Penalty parameter for the twisting rotation of the primary field
     ***********/
    void setParametersWeakPatchContinuityConditions(bool _isWeakPatchContinuityConditions = false, bool _isAutomaticPenaltyParameters = false,
                                                    bool _isPrimCoupled = false, bool _isSecBendingCoupled = false, bool _isSecTwistingCoupled = false,
                                                    double _alphaPrim = 0.0, double _alphaSecBending = 0.0, double _alphaSecTwisting = 0.0);

    /***********************************************************************************************
     * \brief Set the parameters for the application of strong Dirichlet boundary conditions along curves
     * \param[in] _isStrongCurveDirichletConditions Flag on whether strong Dirichlet boundary conditions along curves are assumed
     * \author Andreas Apostolatos
     ***********/
    void setParametersStrongCurveDirichletConditions(bool _isStrongCurveDirichletConditions = false);

    /***********************************************************************************************
     * \brief Set parameters for the error computation
     * \param[in] _isErrorComputation Flag if error computation is assumed
     * \param[in] _isDomainError Flag on the computation of the error from the mapping in the domain
     * \param[in] _isCurveError Flag on the computation of the error of the constaint satisfaction along the trimming curves
     * \param[in] _isInterfaceError Flag on the computation of the interface error between the patches
     * \author Andreas Apostolatos
     ***********/
    void setParametersErrorComputation(bool _isErrorComputation = false, bool _isDomainError = false,
                                       bool _isCurveError = false, bool _isInterfaceError = false);

    /***********************************************************************************************
     * \brief Build the coupling matrcies Cnn and Cnr
     * \author Fabien Pean
     ***********/
    void buildCouplingMatrices();

    /***********************************************************************************************
     * \brief Perform consistent mapping
     * \param[in] _slaveField The reference field
     * \param[in/out] _masterField the target field
     * \author Andreas Apostolatos
     ***********/
    void consistentMapping(const double *_slaveField, double *_masterField);

    /***********************************************************************************************
     * \brief Perform conservative mapping (use it only when mapping forces)
     * \param[in] _masterField The target field
     * \param[in/out] _slaveField The reference field
     * \author Andreas Apostolatos
     ***********/
    void conservativeMapping(const double *_masterField, double *_slaveField);

    /***********************************************************************************************
     * \brief Compute the mapping errors
     * \param[in] _slaveField The reference field
     * \param[in] _masterField The target field
     * \author Andreas Apostolatos
     ***********/
    void computeErrorsConsistentMapping(const double* _slaveField, const double *_masterField);

    /***********************************************************************************************
     * \brief Compute the relative error in the L2 norm for the consistent mapping
     * \param[in] _slaveField The field to be mapped
     * \param[in] _masterField The mapped field
     * \param[out] The relative error in the L2 norm for the consistent mapping
     * \author Andreas Apostolatos
     ***********/
    double computeDomainErrorInL2Norm4ConsistentMapping(const double *_slaveField, const double *_masterField);

    /***********************************************************************************************
     * \brief Compute the relative error in the L2 norm of the error along the trimming curves where conditions are applied in terms of the primary and the secondary fields
     * \param[in/out] _errorL2Interface Double array of constant size 2 containing the L2 norm of the error along the curves where conditions are applied in terms of the primary and the secondary fields
     * \param[in] _fieldIGA The field on the isogeometric discretization
     * \author Andreas Apostolatos
     ***********/
    void computeIGADirichletCurveErrorInL2Norm(double* _errorL2Curve, const double *_fieldIGA);

    /***********************************************************************************************
     * \brief Compute the relative error in the L2 norm of the interface error in terms of the displacements and the rotations across the patch interfaces
     * \param[in/out] _errorL2Interface Double array of constant size 2 containing the L2 norm of the patch interface error in terms of the displacements and the rotations
     * \param[in] _fieldIGA The field on the isogeometric discretization
     * \author Andreas Apostolatos
     ***********/
    void computeIGAPatchInterfaceErrorInL2Norm(double* _errorL2Interface, const double *_fieldIGA);

    /// intern function used for mapping
private:
    /***********************************************************************************************
     * \brief Initialization of the element freedom tables
     * \author Andreas Apostolatos
     ***********/
    void initTables();

    /***********************************************************************************************
     * \brief Initialization of the coupling matrices
     * \author Andreas Apostolatos
     ***********/
    void initCouplingMatrices();

    /***********************************************************************************************
     * \brief Fills up the array projectedCoords by performing closest point projection
     * \author Andreas Apostolatos
     ***********/
    void projectPointsToSurface();

    /***********************************************************************************************
     * \brief Compute the initial guess for a node of an element
     * \param[in] _patchCount The index of the patch we are working on
     * \param[in] _elem The index of the element we are working with
     * \param[in] _localNode The index of the node in the element we are working with
     * \param[out] _u The output guess in u direction
     * \param[out] _v The output guess in v direction
     * \author Fabien Pean
     ***********/

    void computeInitialGuessForProjection(const int _patchCount, const int _elem, const int _localNode, double& _u, double& _v);
    /***********************************************************************************************
     * \brief Compute the projection of a point on a patch using Newton-Raphson
     * \param[in] _patchIndex The index of the patch we are working on
     * \param[in] _nodeIndex The global index of the node in the element we are working with
     * \param[in] _u The initial guess in u direction
     * \param[in] _v The initial guess in v direction
     * \param[out] _minProjectionDistance The previous distance computed
     * \param[out] _minProjectionPoint The previous point computed
     * \author Fabien Pean
     ***********/
    bool projectPointOnPatch(const int patchIndex, const int nodeIndex, const double u0, const double v0, double& minProjectionDistance, std::vector<double>& minProjectionPoint);

    /***********************************************************************************************
     * \brief Compute the projection of a point on a patch using Newton-Raphson
     * \param[in] _patchIndex The index of the patch we are working on
     * \param[in] _nodeIndex The global index of the node in the element we are working with
     * \param[in] _u The initial guess in u direction
     * \param[in] _v The initial guess in v direction
     * \param[out] _minProjectionDistance The previous distance computed
     * \param[out] _minProjectionPoint The previous point computed
     * \author Altug Emiroglu
     ***********/
    bool forceProjectPointOnPatchByRelaxation(const int patchIndex, const int nodeIndex, const double u0, const double v0, double& minProjectionDistance, std::vector<double>& minProjectionPoint);

    /***********************************************************************************************
     * \brief Compute the projection of a point on a patch using a brute force method
     * \param[in] _patchCount The index of the patch we are working on
     * \param[in] _nodeIndex The global index of the node in the element we are working with
     * \param[in] _u The initial guess in u direction
     * \param[in] _v The initial guess in v direction
     * \author Fabien Pean
     ***********/
    bool forceProjectPointOnPatchBySampling(const int patchIndex, const int nodeIndex, double& minProjectionDistance, std::vector<double>& minProjectionPoint);

    /***********************************************************************************************
     * \brief Compute matrices Cnn and Cnr by looping over the FE elements and processing them
     * \author Andreas Apostolatos, Fabien Pean
     ***********/
    void computeCouplingMatrices();

public:
    /***********************************************************************************************
     * \brief Create a Gauss quadrature rule for triangles and for quadrilaterals for each patch of the B-Spline multipatch geometry
     * \author Andreas Apostolatos
     ***********/
    void createGaussQuadratureRules();

    /***********************************************************************************************
     * \brief Compute and assemble the IGA weak Dirichlet curve condition matrices
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    void computeIGAWeakDirichletCurveConditionMatrices();

    /***********************************************************************************************
     * \brief Compute and assemble the IGA weak Dirichlet surface condition matrices
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    void computeIGAWeakDirichletSurfaceConditionMatrices();

    /***********************************************************************************************
     * \brief Compute and assemble the IGA Patch weak continuity condition matrices
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    void computeIGAPatchWeakContinuityConditionMatrices();

    /***********************************************************************************************
     * \brief Compute the B-operator matrices for continuity condition enforcement on a patch
     * \param[in/out] _BDisplacementsGC Pointer to the B-operator matrix for the displacement field in the global Cartesian space
     * \param[in/out] _BOperatorOmegaT Pointer to the B-operator matrix for the bending rotational field in the global Cartesian space
     * \param[in/out] _BOperatorOmegaN Pointer to the B-operator matrix for the twisting rotational field in the global Cartesian space
     * \param[in/out] _normalTrCurveVct Pointer to the normal to the trimming curve vector which is also tangent to the surface patch
     * \param[in/out] _surfaceNormalVct Pointer to the surface normal vector
     * \param[in] _patch Pointer to the patch for which the B-operator matrices are computed
     * \param[in] _tangentTrCurveVct Pointer to the tangent along the trimming curve vector
     * \param[in] _u The u parametric coordinate of the curve in the surface parameter space
     * \param[in] _v The v parametric coordinate of the curve in the surface parameter space
     * \param[in] _uKnotSpan The knot span index in the u-parametric direction
     * \param[in] _vKnotSpan The knot span index in the v-parametric direction
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    void computeDisplacementAndRotationBOperatorMatrices(double* _BDisplacementsGC, double* _BOperatorOmegaT,
                                                         double* _BOperatorOmegaN, double* _normalTrCurveVct,
                                                         double *_surfaceNormalVct, IGAPatchSurface* _patch,
                                                         double* _tangentTrCurveVct, double _u, double _v,
                                                         int _uKnotSpan, int _vKnotSpan);

    /***********************************************************************************************
     * \brief Compute and write the penalty factors for the primary and the secondary field for the weak Dirichlet curve conditions
     * [in] _filename The name of the file where the Penalty parameters are written out
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    void computePenaltyParametersForWeakDirichletCurveConditions(std::string _filename);

    /***********************************************************************************************
     * \brief Compute the penalty factors for the primary and the secondary field for the weak Dirichlet surface conditions
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    void computePenaltyParametersForWeakDirichletSurfaceConditions();

    /***********************************************************************************************
     * \brief Compute and write the penalty factors for the primary and the secondary field for each interface
     * [in] _filename The name of the file where the Penalty parameters are written out
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    void computePenaltyParametersForPatchContinuityConditions(std::string _filename);

private:
    /***********************************************************************************************
     * \brief Get the patches index on which the FE-side element is projected
     * \param[in] _elemIndex The index of the element one is getting the patches for
     * \param[out] _patchWithFullElt The set of patch indexes on which the element is fully projected (All nodes inside)
     * \param[out] _patchWithSplitElement The set of patch indexes on which the element is partially projected (At least 1 outside)
     * \author Fabien Pean
     ***********/
    void getPatchesIndexElementIsOn(int _elemIndex, std::set<int>& _patchWithFullElt, std::set<int>& _patchWithSplitElt);

    /***********************************************************************************************
     * \brief Build the polygon that further steps will work on in case where all nodes of element lies in patch.
     * \param[in] elemIndex The id of current element
     * \param[in] numNodesElementFE The nuöber of node in the current element
     * \param[in] patchIndex The index of the patch we are working on
     * \param[out] polygonUV The output polygon containing parametric coordinates on the patch from the element
     * \author Fabien Pean
     ***********/
    void buildFullParametricElement(int elemIndex, int numNodesElementFE, int patchIndex, Polygon2D& polygonUV);

    /***********************************************************************************************
     * \brief Build the polygon that further steps will work on in case where some nodes of element are out of patch.
     * \param[in] elemIndex The element index
     * \param[in] numNodesElementFE The nuöber of node in the current element
     * \param[in] patchIndex The index of the patch we are working on
     * \param[out] polygonUV The output polygon containing parametric coordinates on the patch from the element
     * \author Fabien Pean
     ***********/
    void buildBoundaryParametricElement(int elemIndex, int numNodesElementFE, int patchIndex, Polygon2D& polygonUV);

    /***********************************************************************************************
     * \brief Compute the projection of a line on patch boundary, display warnings and backup solution
     * \param[in] _thePatch The patch to compute the coupling matrices for
     * \param[in/out] _u The parameter value of the inside patch node
     * \param[in/out] _v The parameter value of the inside patch node
     * \param[in/out] _lambda The ratio on the line
     * \param[in/out] _distance The distance of the line to the patch
     * \param[in] _Pin The point of the line that could have been projected in the patch
     * ]param[in] _Pout The point of the line that could not have been projected in the patch
     * \return Flag if it has converged
     * \author Fabien Pean
     ***********/
    bool projectLineOnPatchBoundary(IGAPatchSurface* _thePatch, double& _u, double& _v, double& _lambda, double& _distance, double* _Pin, double* _Pout);

    /***********************************************************************************************
     * \brief Compute local matrices Cnn and Cnr for element _elemIndex  on patch _patchIndex
     * \param[in] _elemIndex The element index
     * \param[in] _patchIndex The patch index
     * \param[in/out] _projectedElement The polygon containing parametric coordinates of projected element on patch
     * \author Fabien Pean
     ***********/
    bool computeLocalCouplingMatrix(const int _elemIndex, const int _patchIndex, Polygon2D& _projectedElement);

    /***********************************************************************************************
     * \brief Clip the input polygon by the patch parametric quad
     * \param[in] _thePatch The patch for which boundaries are used for clipping
     * \param[in] _polygonUV An input polygon defined in parametric (i.e. 2D) space
     * \author Fabien Pean
     ***********/
    void clipByPatch(const IGAPatchSurface* _thePatch, Polygon2D& _polygonUV);

    /***********************************************************************************************
     * \brief Clip the input polygon by the trimming window of the patch
     * \param[in] _thePatch The patch for which trimming curves are used
     * \param[in] _polygonUV An input polygon defined in parametric (i.e. 2D) space
     * \param[out] _listPolygon	A set of polygons after application of trimming polygon
     * \author Andreas Apostolatos, Fabien Pean
     ***********/
    void clipByTrimming(const IGAPatchSurface* _thePatch, const Polygon2D& _polygonUV, ListPolygon2D& _listPolygonUV);

    /***********************************************************************************************
     * \brief Clip the input polygon by the trimming window of the trimming loop
     * This function does the same as above function but for a specific given trimming loop
     * \param[in] _theTrimmingLoop The trimming loop for which trimming curves are used
     * \param[in] _polygonUV An input polygon defined in parametric (i.e. 2D) space
     * \param[out] _listPolygon	A set of polygons after application of trimming polygon
     * \author Altug Emiroglu
     ***********/
    void clipByTrimming(const IGAPatchSurfaceTrimmingLoop* _theTrimmingLoop, const Polygon2D& _polygonUV, ListPolygon2D& _listPolygonUV);

    /***********************************************************************************************
     * \brief Clip the input polygon for every knot span of the patch it is crossing
     * \param[in] _thePatch The patch for which trimming curves are used
     * \param[in] _polygonUV An input polygon defined in parametric (i.e. 2D) space
     * \param[out] _listPolygon	A set of polygons after application of knot span clipping
     * \param[out] _listSpan The list of span index every polygon of the list above is linked to
     * \author Fabien Pean
     ***********/
    void clipByKnotSpan(const IGAPatchSurface* _thePatch, const Polygon2D& _polygonUV, ListPolygon2D& _listPolygon, Polygon2D& _listSpan);

    /***********************************************************************************************
     * \brief triangulate optimally a 2D input polygon
     * \param[in] _polygonUV 	An input polygon defined in parametric (i.e. 2D) space
     * \author Fabien Pean
     ***********/
    ListPolygon2D triangulatePolygon(const Polygon2D& _polygonUV);

    /***********************************************************************************************
     * \brief Subdivides the input polygon according to member numDivision and compute the canonical element
     * \param[in] _elementIndex	The element index for which the canonical space is related to
     * \param[in] _theElement The polygon of the projected element
     * \param[in] _polygonUV The polygon of a subelement
     * \return The polygon in canonical space of the polygonUV which is defined in nurbs parametric space
     * \author Fabien Pean
     ***********/
    Polygon2D computeCanonicalElement(const int _elementIndex, const Polygon2D& _theElement, const Polygon2D& _polygonUV);

    /***********************************************************************************************
     * \brief Computes the minimum element area size in the multipatch geometry neglecting trimmed elements
     * \author Andreas Apostolatos
     ***********/
    void computeMinimumElementAreaSize();

    /***********************************************************************************************
     * \brief Computes the minimum edge size in the Finite Element mesh
     * \author Andreas Apostolatos
     ***********/
    void computeMinimumEdgeSize();

    /***********************************************************************************************
     * \brief Integrate the element coupling matrices and assemble them to the global one
     * \param[in] _thePatch The patch to compute the coupling matrices for
     * \param[in] _patchIndex The index of the patch as stored in IGAMesh
     * \param[in] _polygonIGA The resulting from the clipping polygon at each knot span in the NURBS space
     * \param[in] _spanU The knot span index in the u-direction where basis will be evaluated
     * \param[in] _spanV The knot span index in the v-direction where basis will be evaluated
     * \param[in] _polygonFE The resulting from the clipping polygon at each knot span in the bilinear/linear space
     * \param[in] _elementIndex The global numbering of the element from the FE mesh the shape functions are evaluated for
     * \author Andreas Apostolatos, Fabien Pean
     ***********/
    void integrate(IGAPatchSurface* _thePatch, int _patchIndex, Polygon2D _polygonIGA,
                   int _spanU, int _spanV, Polygon2D _polygonFE, int _elementIndex);

    /// helper functions used for the computation of coupling matrix process
private:
    /***********************************************************************************************
     * \brief Compute the span of the projected element living in _thePatch
     * \param[in] _thePatch The patch to compute the coupling matrices for
     * \param[in] _polygonUV The resulting from the clipping polygon at each knot span in the NURBS space
     * \param[out] _span An array size 4 containing [minSpanU maxSpanU minSpanV maxSpanV]
     * \return True if inside a single knot span, false otherwise
     * \author Andreas Apostolatos, Fabien Pean
     ***********/
    bool computeKnotSpanOfProjElement(const IGAPatchSurface* _thePatch, const Polygon2D& _polygonUV, int* _span=NULL);

    /***********************************************************************************************
     * \brief Get the element id in the FE mesh of the neighbor of edge made up by node1 and node2
     * \param[in] _element The element for which we want the neighbor
     * \param[in] _node1 The first node index of the edge
     * \param[in] _node2 The second node of the edge
     * \return The index of the neighbor element, or -1 if does not exist
     * \author Fabien Pean
     ***********/
    int getNeighbourElementofEdge(int _element, int _node1, int _node2);

    /// Writing output functions
public:
    void writeGaussPointData();
    /***********************************************************************************************
     * \brief Writes the projected nodes of the FE mesh onto the IGA surface into a file
     * \author Andreas Apostolatos
     ***********/
    void writeProjectedNodesOntoIGAMesh();

    /***********************************************************************************************
     * \brief Writes the back projection of projected FE element in a Paraview (polydata vtk) format
     * 		Opens a file filename.csv, process it and write filename.vtk
     * \param[in] _filename		The substring to append to open csv file and write vtk file
     * \author Fabien Pean
     ***********/
    void writeCartesianProjectedPolygon(const std::string _filename, std::map<int,ListPolygon2D>& _data);

    /***********************************************************************************************
     * \brief Writes all FE mesh in parametric coordinates
     * \author Fabien Pean
     ***********/
    void writeParametricProjectedPolygons(std::string _filename);

    /***********************************************************************************************
     * \brief Writes all triangulated polygons to be integrated
     * \author Fabien Pean
     ***********/
    void writeTriangulatedParametricPolygon(std::string _filename);

    /***********************************************************************************************
     * \brief Print both coupling matrices Cnn and Cnr in file in csv format with space delimiter
     * \author Fabien Pean
     ***********/
    void writeCouplingMatricesToFile();

    /// Debugging functions
public:
    /***********************************************************************************************
     * \brief Print a polygon in debug stream
     * \author Fabien Pean
     ***********/
    void debugPolygon(const Polygon2D& _polygon, std::string _name="");

    /***********************************************************************************************
     * \brief Print a set of polygon
     * \author Fabien Pean
     ***********/
    void debugPolygon(const ListPolygon2D& _listPolygon, std::string _name="");

    /***********************************************************************************************
     * \brief Print both coupling matrices Cnn and Cnr
     * \author Andreas Apostolatos
     ***********/
    void printCouplingMatrices();

    /***********************************************************************************************
     * \brief Check and enforce consistency of the mapper based on the mapping of a unit field
     * \author Andreas Apostolatos, Fabien Pean
     ***********/
    void enforceConsistency();

    /// Get functions
public:
    /***********************************************************************************************
     * \brief Get the IGAMesh containing all B-Spline surface patches which form the geometry
     * \author Andreas Apostolatos
     ***********/
    IGAMesh* getIGAMesh() {
        return meshIGA;
    }

    /***********************************************************************************************
     * \brief Get couplingMatrices object
     * \author Andreas Apostolatos
     ***********/
    IGAMortarCouplingMatrices* getCouplingMatrices() {
        return couplingMatrices;
    }

    /***********************************************************************************************
     * \brief Get flag on whether a IGA weak Dirichlet curve condition was used or not
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    bool getIsIGAWeakDirichletCurveConditionPenalties() {
        return propWeakCurveDirichletConditions.isWeakCurveDirichletConditions;
    }

    /***********************************************************************************************
     * \brief Get flag on whether a IGA weak Dirichlet surface condition was used or not
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    bool getIsIGAWeakDirichletSurfaceConditionPenalties() {
        return propWeakSurfaceDirichletConditions.isWeakSurfaceDirichletConditions;
    }

    /***********************************************************************************************
     * \brief Get flag on whether the IGA patch coupling was used or not
     * \author Andreas Apostolatos
     ***********/
    bool getIsIGAPatchCouplingPenalties() {
        return propWeakPatchContinuityConditions.isWeakPatchContinuityConditions;
    }

    /***********************************************************************************************
     * \brief Get flag on whether the expanded version of the coupling matrices is assumed
     * \author Andreas Apostolatos
     ***********/
    bool getIsExpanded() {
        return isExpanded;
    }

    /***********************************************************************************************
     * \brief Get the number of the weak patch continuity conditions
     * \author Altug Emiroglu
     ***********/
    int getNoWeakIGADirichletCurveConditions(){
        return noWeakIGADirichletCurveConditions;
    }

    /***********************************************************************************************
     * \brief Get the number of the weak patch continuity conditions
     * \author Andreas Apostolatos
     ***********/
    int getNoWeakIGAPatchContinuityConditions(){
        return noWeakIGAPatchContinuityConditions;
    }

    /***********************************************************************************************
     * \brief Get the flag on whether error computation is enabled
     * \author Andreas Apostolatos
     ***********/
    bool getIsErrorComputation() {
        return propErrorComputation.isErrorComputation;
    }

    /***********************************************************************************************
     * \brief Get the minimum element area of the multipatch geometry
     * \return The minimum element area in the multipatch geometry
     * \author Andreas Apostolatos
     ***********/
    double getMinElArea() {
        return minElArea;
    }

    /***********************************************************************************************
     * \brief Get the minimum edge size of the Finite Element mesh
     * \return The minimum edge size of the Finite Element mesh
     * \author Andreas Apostolatos
     ***********/
    double getMinEdgeSize() {
        return minEdgeSize;
    }

    /***********************************************************************************************
     * \brief Get the minimum edge size at the Dirichlet boundary
     * \return The minimum edge size at the Dirichlet boundary
     * \author Andreas Apostolatos
     ***********/
    double getMinElEdgeSizeDirichlet() {
        return minElEdgeSizeDirichlet;
    }

    /***********************************************************************************************
     * \brief Get the minimum edge size at the interface between the patches in the multipatch geometry
     * \return The minimum edge size between the patch interfaces in the multipatch geometry
     * \author Andreas Apostolatos
     ***********/
    double getMinElEdgeSizeInterface() {
        return minElEdgeSizeInterface;
    }

    /***********************************************************************************************
     * \brief Get the penalty parameters of the primary field
     * \param[in/out] The vector of the penalty factors for the primary field
     * \author Altug Emiroglu
     ***********/
    void getPenaltyParameterForWeakDirichletCCPrimaryField(double* _alphaPrim);

    /***********************************************************************************************
     * \brief Get the penalty parameters for the bending rotation of the secondary field
     * \param[in/out] The vector of the penalty factors for the bending rotation of the primary field
     * \author Altug Emiroglu
     ***********/
    void getPenaltyParameterForWeakDirichletCCSecondaryFieldBending(double* _alphaSecBending);

    /***********************************************************************************************
     * \brief Get the penalty parameters of the twisting rotation of the secondary field
     * \param[in/out] The vector of the penalty factors for the twisting rotation of the primary field
     * \author Altug Emiroglu
     ***********/
    void getPenaltyParameterForWeakDirichletCCSecondaryFieldTwisting(double* _alphaSecTwisting);

    /***********************************************************************************************
     * \brief Get the penalty parameters of the primary field
     * \param[in/out] The vector of the penalty factors for the primary field
     * \author Andreas Apostolatos
     ***********/
    void getPenaltyParameterForPatchContinuityPrimaryField(double* _alphaPrim);

    /***********************************************************************************************
     * \brief Get the penalty parameters of the bending rotation of the secondary field
     * \param[in/out] The vector of the penalty factors for the bending rotation of the secondary field
     * \author Andreas Apostolatos
     ***********/
    void getPenaltyParameterForPatchContinuitySecondaryFieldBending(double* _alphaSecBending);

    /***********************************************************************************************
     * \brief Get the penalty parameters of the twisting rotation of the secondary field
     * \param[in/out] The vector of the penalty factors for the twisting rotation of the secondary field
     * \author Andreas Apostolatos
     ***********/
    void getPenaltyParameterForPatchContinuitySecondaryFieldTwisting(double* _alphaSecTwisting);

    /// Print functions
public:
    /***********************************************************************************************
     * \brief Print the domain and the interface error from the mapping
     * \param[in] message Reference to the message
     * \param[in] _domainError The computed error from the mapping in the domain
     * \param[in] _errorL2Curve The computed errors across the curves where conditions are applied from the mapped field
     * \param[in] _errorL2Interface The computed interface errors accross the patch interfaces from the mapped field
     * \author Andreas Apostolatos
     ***********/
    void printErrorMessage(Message &message, double _errorL2Domain, double* _errorL2Curve, double *_errorL2Interface);

    // Constant members of the class
private:
    /// Tolerance for cleaning a triangle before integrating
    static double EPS_CLEANTRIANGLE;
    static double EPS_CLIPPING;

    /// unit test class
    friend class TestIGAMortarMapperTube;
    friend class TestIGAMortarMapperMultiPatchPlanarSurface;
    friend class TestIGAMortarMapperCylinder;

};

extern Message infoOut;

}

#endif /* IGAMORTARMAPPER_H_ */
