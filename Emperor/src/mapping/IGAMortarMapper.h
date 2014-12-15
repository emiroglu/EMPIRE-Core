/*  Copyright &copy; 2013, TU Muenchen, Chair of Structural Analysis,
 *  Fabien Pean, Andreas Apostolatos, Chenshen Wu,
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

namespace EMPIRE {

namespace MathLibrary {
template<class T> class SparseMatrix;
}

class IGAPatchSurface;
class IGAMesh;
class FEMesh;
class DataField;

namespace MathLibrary {

class IGAGaussQuadratureOnTriangle;
class IGAGaussQuadratureOnQuad;

}

/***********************************************************************************************
 * \brief This class is related to IGA mortar mapping
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

    /// The element freedom table for the fluid mesh
    int **meshFEDirectElemTable;

    /// The mass-like matrix
    MathLibrary::SparseMatrix<double> *C_NN;

    /// The right-hand side matrix
    MathLibrary::SparseMatrix<double> *C_NR;

    /// Quadrature rule over the triangulated subdomains
    EMPIRE::MathLibrary::IGAGaussQuadratureOnTriangle *gaussTriangle;

    /// Quadrature rule over the non-triangulated subdomains
    EMPIRE::MathLibrary::IGAGaussQuadratureOnQuad *gaussQuad;

    /// The parametric coordinates of the projected nodes on the surface
    std::vector<std::map<int, std::vector<double> > > projectedCoords;

    /// Flag on the mapping direction
    bool isMappingIGA2FEM;

    size_t numNodesSlave;
    size_t numNodesMaster;

    struct integration {
        int numGPTriangle;
        int numGPQuad;
    } integration;
    struct nonlinearSchemeProperties {
        int maxNumOfIterations;
        double tolerance;
    } newtonRaphson, newtonRaphsonBoundary, bisection;
    struct projectionProperties {
        double maxProjectionDistance;
        int numRefinementForIntialGuess;
        int maxDistanceForProjectedPointsOnDifferentPatches;
    } projectionProperties;

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
     * \author Fabien Pean, Chenshen Wu
     ***********/
    virtual ~IGAMortarMapper();

    /***********************************************************************************************
     * \brief Set parameter for the projection of mesh onto the NURBS surface
     * \param[in] _maxProjectionDistance The max distance allowed between FE mesh and NURBS surface
     * \param[in] _numRefinementForIntialGuess The number of test point to find initial guess for Newton-Raphson scheme
     * \param[in] _maxDistanceForProjectedPointsOnDifferentPatches The max authorized distance between two projected points from a same physical node
     * \author Fabien Pean
     ***********/
    void setParametersProjection(double _maxProjectionDistance, int _numRefinementForIntialGuess,
                                 int _maxDistanceForProjectedPointsOnDifferentPatches);
    /***********************************************************************************************
     * \brief Set parameter for Newton-Raphson scheme of projection on NURBS patch
     * \param[in] _newtonRaphsonMaxIt The number of iteration for Newton-Raphson scheme of projecting a node on a NURBS patch
     * \param[in] _newtonRaphsonTol The tolerance for Newton-Raphson scheme of projecting a node on a NURBS patch
     * \author Fabien Pean
     ***********/
    void setParametersNewtonRaphson(int _maxNumOfIterations=20, double _tolerance=1e-6);
    /***********************************************************************************************
     * \brief Set parameter for Newton-Raphson scheme of projection on NURBS patch boundary
     * \param[in] _newtonRaphsonBoundaryMaxIt The number of iteration for Newton-Raphson scheme of projecting a node on a NURBS patch boundary
     * \param[in] _newtonRaphsonBoundaryTol The tolerance for Newton-Raphson scheme of projecting a node on a NURBS patch boundary
     * \author Fabien Pean
     ***********/
    void setParametersNewtonRaphsonBoundary(int _maxNumOfIterations=20, double _tolerance=1e-6);
    /***********************************************************************************************
     * \brief Set parameter for bisection scheme of projection on NURBS patch boundary
     * \param[in] _bisectionMaxIt The number of iteration for bisection scheme of projecting a node on a NURBS patch boundary
     * \param[in] _bisectionTol The tolerance for bisection scheme of projecting a node on a NURBS patch boundary
     * \author Fabien Pean
     ***********/
    void setParametersBisection(int _maxNumOfIterations=20, double _tolerance=1e-6);
    /***********************************************************************************************
     * \brief Set parameter for integration
     * \param[in] _numGPsTriangle The number of Gauss points when performs integration on triangle
     * \param[in] _numGPsQuad The number of Gauss points when performs integration on quadrilateral
     ***********/
    void setParametersIntegration(int _numGPTriangle=16, int _numGPQuad=25);

    /***********************************************************************************************
     * \brief Build the coupling matrcies C_NN and C_NR
     * \author Fabien Pean
     ***********/
    void buildCouplingMatrices();

    /***********************************************************************************************
     * \brief Perform consistent mapping from IGA to FE (map displacements)
     * \param[in] fieldIGA is the input data
     * \param[out] fieldFE is the output data
     * \author Chenshen Wu
     ***********/
    void consistentMapping(const double *fieldIGA, double *fieldFE);

    /***********************************************************************************************
     * \brief Perform conservative mapping from FE to IGA (map forces)
     * \param[in] fieldFE is the input data
     * \param[out] fieldIGA is the output data
     * \author Chenshen Wu
     ***********/
    void conservativeMapping(const double *fieldFE, double *fieldIGA);


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
     * \brief Compute matrices C_NN and C_NR
     * \author Fabien Pean, Chenshen Wu
     ***********/
    void computeCouplingMatrices();

    /***********************************************************************************************
     * \brief Compute local matrices C_NN and C_NR for element i  on patch j
     * \author Fabien Pean
     ***********/
    bool computeLocalCouplingMatrix(const int _elemIndex, const int _patchIndex, Polygon2D& _projectedElement);

    /***********************************************************************************************
     * \brief Get the patches index on which the FE-side element is projected
     * \param[in] _elemIndex The index of the element one is getting the patches for
     * \param[out] _patchWithFullElt The set of patch indexes on which the element is fully projected (All nodes inside)
     * \param[out] _patchWithSplitElement The set of patch indexes on which the element is partially projected (At least 1 outside)
     * \author Fabien Pean
     ***********/
    void getPatchesIndexElementIsOn(int _elemIndex, std::set<int>& _patchWithFullElt, std::set<int>& _patchWithSplitElt);

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
     * \brief Subdivides the input polygon according to member numDivision and compute the canonical element
     * \param[in] _thePatch 	The patch for which trimming curves are used
     * \param[in] _polygonUV 	An input polygon defined in parametric (i.e. 2D) space
     * \param[in] _elementIndex	The element for which the canonical space is related to
     * \return					The polygon in canonical space of the polygonUV which is defined in nurbs parametric space
     * \author Fabien Pean
     ***********/
    Polygon2D computeCanonicalElement(IGAPatchSurface* _thePatch, Polygon2D& _polygonUV, int _elementIndex);

    /***********************************************************************************************
     * \brief Integrate the element coupling matrices and assemble them to the global one
     * \param[in] _igaPatchSurface 	The patch to compute the coupling matrices for
     * \param[in] _polygonIGA 		The resulting from the clipping polygon at each knot span in the NURBS space
     * \param[in] _spanU 			The knot span index in the u-direction where basis will be evaluated
     * \param[in] _spanV 			The knot span index in the v-direction where basis will be evaluated
     * \param[in] _polygonFE 		The resulting from the clipping polygon at each knot span in the bilinear/linear space
     * \param[in] _elementIndex 	The global numbering of the element from the FE mesh the shape functions are evaluated for
     * \author Chenshen Wu, Fabien Pean
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

    /***********************************************************************************************
     * \brief Compute numDivision points projection between P1 and P2 and add them in parametric polygon
     * \param[in] _patchIndex		The index of the patch to project on
     * \param[in] _elemCount	 	The first of the element to project
     * \param[in] _nodeIndex1 		The index of the first node
     * \param[in] _nodeIndex2 		The index of the second node
     * \param[in] _P1				The first point of the line on FE element
     * \param[in] _isIn1			The first point was projected inside the patch
     * \param[in] _P2				The second point of the line on FE element
     * \param[in] _isIn2			The second point was projected inside the patch
     * \param[in/out]				The parametric polygon of the projection of the FE element
     * \return 1 if no problem have been found on the boundary
     * \author Fabien Pean
     ***********/
    bool computeIntermediatePoints(const int patchIndex, const int elemCount, const int nodeIndex1,const int nodeIndex2,
    		const double* P1, const bool isIn1, const double* P2, const bool isIn2,
    		Polygon2D& polygonUV, std::map<int,Polygon2D>* extraPolygonUV=NULL);
    /// Writing output functions
public:
    /***********************************************************************************************
     * \brief Writes the projected nodes of the FE mesh onto the IGA surface into a file
     * \author Andreas Apostolatos
     ***********/
    void writeProjectedNodesOntoIGAMesh();
    /***********************************************************************************************
     * \brief Writes a projected element of the FE mesh onto the IGA surface into a file
     * \param[in] _patchIndex	The id of the patch the parametric coordinates must be used with
     * \param[in] _polygonUV	The parametric coordinates of the polygon to write in the file
     * \author Fabien Pean
     ***********/
    void writeParametricProjectedPolygon(const int _patchIndex = -1, const Polygon2D* const _polygonUV = NULL);
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

    /// unit test class
    friend class TestIGAMortarMapperTube;
    friend class TestIGAMortarMapperMultiPatchPlanarSurface;
    friend class TestIGAMortarMapperCylinder;

};
}

#endif /* IGAMORTARMAPPER_H_ */
