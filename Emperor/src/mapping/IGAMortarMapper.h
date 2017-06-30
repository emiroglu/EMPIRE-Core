/*  Copyright &copy; 2013, TU Muenchen, Chair of Structural Analysis,
 *  Fabien Pean, Andreas Apostolatos, Chenshen Wu,
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

namespace EMPIRE {

namespace MathLibrary {
template<class T> class SparseMatrix;
}

class IGAPatchSurface;
class IGAMesh;
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

    /// The isogeometric coupling matrices
    IGAMortarCouplingMatrices *couplingMatrices;

    /// Flag on the application of patch continuity conditions using Penalty
    bool isIGAPatchContinuityConditions;

    /// Number of weak IGA patch continuity conditions
    int noWeakIGAPatchContinuityConditions;

    /// Penalty factors for the primary field to the application of weak patch continuity conditions
    double* alphaPrimaryIJ;

    /// Penalty factors for the secondary field to the application of weak patch continuity conditions
    double* alphaSecondaryIJ;

    /// Flag on the initialization of the quadrature
    bool isGaussQuadrature;

    /// Quadrature rule over the triangulated subdomains
    EMPIRE::MathLibrary::IGAGaussQuadratureOnTriangle *gaussTriangle;

    /// Quadrature rule over the non-triangulated subdomains
    EMPIRE::MathLibrary::IGAGaussQuadratureOnQuad *gaussQuad;

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
    /// Weight / Jacobian / NumOfFENode / Node1 / ShapeValue1 / Node2 / ShapeValue2 ... NumOfIGANode / Node1 / ShapeValue1/ ...
    std::vector<std::vector<double> > streamGP;

    /// Flag on the mapping direction
    bool isMappingIGA2FEM;

    /// Number of nodes for the slave side
    size_t numNodesSlave;

    /// Number of nodes for the master side
    size_t numNodesMaster;

    /// Number of Gauss Points for the integration over a triangle or a quadrilateral
    struct integration {
        int numGPTriangle;
        int numGPQuad;
    } integration;

    /// Properties for the nonlinear solution schemes
    struct nonlinearSchemeProperties {
        int maxNumOfIterations;
        double tolerance;
    } newtonRaphson, newtonRaphsonBoundary, bisection;

    /// Properties for the projection schemes
    struct projectionProperties {
        double maxProjectionDistance;
        int numRefinementForIntialGuess;
        double maxDistanceForProjectedPointsOnDifferentPatches;
    } projectionProperties;

    /// Properties for the application of weak patch continuity conditions with Penalty
    struct IgaPatchCoupling {
            double dispPenalty;
            double rotPenalty;
            int isAutomaticPenaltyFactors;
    } IgaPatchCoupling;

    /// On the strong application of Dirichlet boundary conditions
    struct dirichletBCs {
        int isDirichletBCs;
    }dirichletBCs;

public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _name The name of the mapper
     * \param[in] _meshIGA The IGAMesh
     * \param[in] _meshFE The FEMesh
     * \param[in] _disTol Tolerance up to which projection is trusted
     * \param[in] _numGPsTri The number of Gauss points used for computing triangle element
     * \param[in] _numGPsQuad The number of Gauss points used for computing quad element
     * \author Fabien Pean, Chenshen Wu
     ***********/
    IGAMortarMapper(std::string _name, IGAMesh *_meshIGA, FEMesh *_meshFE, bool _isMappingIGA2FEM);

    /***********************************************************************************************
     * \brief Destructor
     * \author Fabien Pean
     ***********/
    virtual ~IGAMortarMapper();

    /***********************************************************************************************
     * \brief Set parameter for the projection of mesh onto the NURBS surface
     * \param[in] _maxProjectionDistance The max distance allowed between FE mesh and NURBS surface
     * \param[in] _numRefinementForIntialGuess The number of test point to find initial guess for Newton-Raphson scheme
     * \param[in] _maxDistanceForProjectedPointsOnDifferentPatches The max authorized distance between two projected points from a same physical node
     * \author Fabien Pean
     ***********/
    void setParametersProjection(double _maxProjectionDistance = 1e-2, int _numRefinementForIntialGuess = 10,
                                 double _maxDistanceForProjectedPointsOnDifferentPatches = 1e-3);

    /***********************************************************************************************
     * \brief Set parameter for Newton-Raphson scheme of projection on NURBS patch
     * \param[in] _newtonRaphsonMaxIt The number of iteration for Newton-Raphson scheme of projecting a node on a NURBS patch
     * \param[in] _newtonRaphsonTol The tolerance for Newton-Raphson scheme of projecting a node on a NURBS patch
     * \author Fabien Pean
     ***********/
    void setParametersNewtonRaphson(int _maxNumOfIterations = 20, double _tolerance = 1e-6);

    /***********************************************************************************************
     * \brief Set parameter for Newton-Raphson scheme of projection on NURBS patch boundary
     * \param[in] _newtonRaphsonBoundaryMaxIt The number of iteration for Newton-Raphson scheme of projecting a node on a NURBS patch boundary
     * \param[in] _newtonRaphsonBoundaryTol The tolerance for Newton-Raphson scheme of projecting a node on a NURBS patch boundary
     * \author Fabien Pean
     ***********/
    void setParametersNewtonRaphsonBoundary(int _maxNumOfIterations = 20, double _tolerance = 1e-6);

    /***********************************************************************************************
     * \brief Set parameter for bisection scheme of projection on NURBS patch boundary
     * \param[in] _bisectionMaxIt The number of iteration for bisection scheme of projecting a node on a NURBS patch boundary
     * \param[in] _bisectionTol The tolerance for bisection scheme of projecting a node on a NURBS patch boundary
     * \author Fabien Pean
     ***********/
    void setParametersBisection(int _maxNumOfIterations = 40, double _tolerance = 1e-6);

    /***********************************************************************************************
     * \brief Set parameter for integration
     * \param[in] _numGPsTriangle The number of Gauss points when performs integration on triangle
     * \param[in] _numGPsQuad The number of Gauss points when performs integration on quadrilateral
     ***********/
    void setParametersIntegration(int _numGPTriangle = 16, int _numGPQuad = 25);

    /***********************************************************************************************
     * \brief Set parameter for penalty coupling
     * \param[in] _dispPenalty The displacement penalty coupling factor
     * \param[in] _rotPenalty The rotational penalty coupling factor
     * \param[in] isAutomaticPenaltyFactors flag whether to compute penalty factors automatically or not
     ***********/
    void setParametersIgaPatchCoupling(double _dispPenalty = 0, double _rotPenalty = 0, int isAutomaticPenaltyFactors = 0);

    /***********************************************************************************************
     * \brief Set parameter for penalty coupling
     * \param[in] _isDirichletBCs flag if dirichlet boundary conditions are used or not
     ***********/
    void setParametersDirichletBCs(int _isDirichletBCs = 0);

    /***********************************************************************************************
     * \brief Set the flag regarding the computation of weak continuity conditions
     * \author Ragnar Björnsson
     ***********/
    bool setUseIGAPatchCouplingPenalties(bool _isIGAPatchContinuityConditions) {
        isIGAPatchContinuityConditions = _isIGAPatchContinuityConditions;
    }

    /***********************************************************************************************
     * \brief Build the coupling matrcies C_NN and C_NR
     * \author Fabien Pean
     ***********/
    void buildCouplingMatrices();

    /***********************************************************************************************
     * \brief Perform consistent mapping from IGA to FE (map displacements)
     * \param[in] fieldIGA is the input data
     * \param[in/out] fieldFE is the output data
     * \author Chenshen Wu
     ***********/
    void consistentMapping(const double *fieldIGA, double *fieldFE);

    /***********************************************************************************************
     * \brief Perform conservative mapping from FE to IGA (map forces)
     * \param[in] fieldFE is the input data
     * \param[in/out] fieldIGA is the output data
     * \author Chenshen Wu
     ***********/
    void conservativeMapping(const double *fieldFE, double *fieldIGA);

    /***********************************************************************************************
     * \brief Compute the relative error in the L2 norm for the consistent mapping
     * \param[in] fieldFE The field on the Finite Element mesh
     * \param[in] fieldIGA The field on the isogeometric discretization
     * \param[out] The relative error in the L2 norm for the consistent mapping
     * \author Andreas Apostolatos
     ***********/
    double computeDomainErrorInL2Norm4ConsistentMapping(const double *fieldIGA, const double *fieldFE);

    /// intern function used for mapping
private:
    /***********************************************************************************************
     * \brief Initialization of the element freedom tables
     * \author Chenshen Wu
     ***********/
    void initTables();

    /***********************************************************************************************
     * \brief Fills up the array projectedCoords by performing closest point projection
     * \author Chenshen Wu
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
     * \brief Compute the projection of a point on a patch using a brute force method
     * \param[in] _patchCount The index of the patch we are working on
     * \param[in] _nodeIndex The global index of the node in the element we are working with
     * \param[in] _u The initial guess in u direction
     * \param[in] _v The initial guess in v direction
     * \author Fabien Pean
     ***********/
    bool forceProjectPointOnPatch(const int patchIndex, const int nodeIndex, double& minProjectionDistance, std::vector<double>& minProjectionPoint);

    /***********************************************************************************************
     * \brief Compute matrices C_NN and C_NR by looping over the FE elements and processing them
     * \author Fabien Pean, Chenshen Wu
     ***********/
    void computeCouplingMatrices();

public:

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
     * \param[in] _patch Pointer to the patch for which the B-operator matrices are computed
     * \param[in] _tangentTrCurveVct Pointer to the tangent along the trimming curve vector
     * \param[in] _u The u parametric coordinate of the curve in the surface parameter space
     * \param[in] _v The v parametric coordinate of the curve in the surface parameter space
     * \param[in] _uKnotSpan The knot span index in the u-parametric direction
     * \param[in] _vKnotSpan The knot span index in the v-parametric direction
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    void computeIGAPatchContinuityConditionBOperatorMatrices(double* _BDisplacementsGC, double* _BOperatorOmegaT,
                                                             double* _BOperatorOmegaN, double* _normalTrCurveVct,
                                                             IGAPatchSurface* _patch, double* _tangentTrCurveVct,
                                                             double _u, double _v, int _uKnotSpan, int _vKnotSpan);

    /***********************************************************************************************
     * \brief Compute the penalty factors for the primary and the secondary field for each interface
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    void computePenaltyFactorsForPatchContinuityConditions();

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
     * \param[in] _thePatch 	The patch to compute the coupling matrices for
     * \param[in/out] _u 		The parameter value of the inside patch node
     * \param[in/out] _v 		The parameter value of the inside patch node
     * \param[in/out] _div		The ratio on the line
     * \param[in/out] _dis		The distance of the line to the patch
     * \param[in] _Pin			The point of the line that could have been projected in the patch
     * ]param[in] _Pout			The point of the line that could not have been projected in the patch
     * \return 					Flag if it has converged
     * \author Fabien Pean
     ***********/
    bool projectLineOnPatchBoundary(IGAPatchSurface* _thePatch, double& _u, double& _v, double& _div, double& _dis, double* _Pin, double* _Pout);

    /***********************************************************************************************
     * \brief Compute local matrices C_NN and C_NR for element _elemIndex  on patch _patchIndex
     * \param[in] _elemIndex The element index
     * \param[in] _patchIndex The patch index
     * \param[in/out] _projectedElement The polygon containing parametric coordinates of projected element on patch
     * \author Fabien Pean
     ***********/
    bool computeLocalCouplingMatrix(const int _elemIndex, const int _patchIndex, Polygon2D& _projectedElement);

    /***********************************************************************************************
     * \brief Clip the input polygon by the patch parametric quad
     * \param[in] _thePatch 	The patch for which boundaries are used for clipping
     * \param[in] _polygonUV 	An input polygon defined in parametric (i.e. 2D) space
     * \author Fabien Pean
     ***********/
    void clipByPatch(const IGAPatchSurface* _thePatch, Polygon2D& _polygonUV);

    /***********************************************************************************************
     * \brief Clip the input polygon by the trimming window of the patch
     * \param[in] _thePatch 	The patch for which trimming curves are used
     * \param[in] _polygonUV 	An input polygon defined in parametric (i.e. 2D) space
     * \param[out] _listPolygon	A set of polygons after application of trimming polygon
     * \author Fabien Pean
     ***********/
    void clipByTrimming(const IGAPatchSurface* _thePatch, const Polygon2D& _polygonUV, ListPolygon2D& _listPolygonUV);

    /***********************************************************************************************
     * \brief Clip the input polygon for every knot span of the patch it is crossing
     * \param[in] _thePatch 	The patch for which trimming curves are used
     * \param[in] _polygonUV 	An input polygon defined in parametric (i.e. 2D) space
     * \param[out] _listPolygon	A set of polygons after application of knot span clipping
     * \param[out] _listSpan	The list of span index every polygon of the list above is linked to
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
     * \param[in] _theElement 	The polygon of the projected element
     * \param[in] _polygonUV 	The polygon of a subelement
     * \return					The polygon in canonical space of the polygonUV which is defined in nurbs parametric space
     * \author Fabien Pean
     ***********/
    Polygon2D computeCanonicalElement(const int _elementIndex, const Polygon2D& _theElement, const Polygon2D& _polygonUV);

    /***********************************************************************************************
     * \brief Integrate the element coupling matrices and assemble them to the global one
     * \param[in] _igaPatchSurface 	The patch to compute the coupling matrices for
     * \param[in] _polygonIGA 		The resulting from the clipping polygon at each knot span in the NURBS space
     * \param[in] _spanU 			The knot span index in the u-direction where basis will be evaluated
     * \param[in] _spanV 			The knot span index in the v-direction where basis will be evaluated
     * \param[in] _polygonFE 		The resulting from the clipping polygon at each knot span in the bilinear/linear space
     * \param[in] _elementIndex 	The global numbering of the element from the FE mesh the shape functions are evaluated for
     * \author Fabien Pean, Chenshen Wu
     ***********/
    void integrate(IGAPatchSurface* _igaPatchSurface, Polygon2D _polygonIGA,
            int _spanU, int _spanV, Polygon2D _polygonFE, int _elementIndex);

    /// helper functions used for the computation of coupling matrix process
private:
    /***********************************************************************************************
     * \brief Compute the span of the projected element living in _thePatch
     * \param[in] _thePatch 	The patch to compute the coupling matrices for
     * \param[in] _polygonUV 	The resulting from the clipping polygon at each knot span in the NURBS space
     * \param[out] _span 		An array size 4 containing [minSpanU maxSpanU minSpanV maxSpanV]
     * \return 					True if inside a single knot span, false otherwise
     * \author Fabien Pean
     ***********/
    bool computeKnotSpanOfProjElement(const IGAPatchSurface* _thePatch, const Polygon2D& _polygonUV, int* _span=NULL);

/***********************************************************************************************
     * \brief Get the element id in the FE mesh of the neighbor of edge made up by node1 and node2
     * \param[in] _element		The element for which we want the neighbor
     * \param[in] _node1	 	The first node index of the edge
     * \param[in] _node2 		The second node of the edge
     * \return 					The index of the neighbor element, or -1 if does not exist
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
     * \brief Print both coupling matrices C_NN and C_NR in file in csv format with space delimiter
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
     * \brief Print both coupling matrices C_NN and C_NR
     * \author Chenshen Wu
     ***********/
    void printCouplingMatrices();

    /***********************************************************************************************
     * \brief Check consistency of the coupling, constant field gives constant field
     * \author Fabien Pean
     ***********/
    void checkConsistency();

public:
    /***********************************************************************************************
     * \brief Get boolean whether the IGA patch coupling was used or not
     * \author Ragnar Björnsson
     ***********/
    bool getUseIGAPatchCouplingPenalties() {
        return isIGAPatchContinuityConditions;
    }

    /***********************************************************************************************
     * \brief Get couplingMatrices object
     * \author Ragnar Björnsson
     ***********/
    IGAMortarCouplingMatrices* getCouplingMatrices() {
        return couplingMatrices;
    }

    /***********************************************************************************************
     * \brief Get the number of the weak patch continuity conditions
     * \author Andreas Apostolatos
     ***********/
    int getNoWeakIGAPatchContinuityConditions(){
        return noWeakIGAPatchContinuityConditions;
    }

    /***********************************************************************************************
     * \brief Get the penalty parameters of the primary field
     * \param[in/out] The vector of the penalty factors for the primary field
     * \author Andreas Apostolatos
     ***********/
    void getPenaltyParameterForPrimaryField(double* _alphaPrim);

    /***********************************************************************************************
     * \brief Get the penalty parameters of the secondary field
     * \param[in/out] The vector of the penalty factors for the secondary field
     * \author Andreas Apostolatos
     ***********/
    void getPenaltyParameterForSecondaryField(double* _alphaSec);

    /// unit test class
    friend class TestIGAMortarMapperTube;
    friend class TestIGAMortarMapperMultiPatchPlanarSurface;
    friend class TestIGAMortarMapperCylinder;

};
}

#endif /* IGAMORTARMAPPER_H_ */
