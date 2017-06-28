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
 * \file IGABarycentricMapper.h
 * The header file of class IGABarycentricMapper.
 * \date 3/6/2013
 *********************************************************************************/

#ifndef IGABARYCENTRICMAPPER_H_
#define IGABARYCENTRICMAPPER_H_

#include <map>
#include <set>
#include <vector>
#include <iostream>
#include <assert.h>
#include "AbstractMapper.h"
#include "IGAControlPoint.h"

namespace EMPIRE {

namespace MathLibrary {
template<class T> class SparseMatrix;
}

class IGAPatchSurface;
class IGAMesh;
class FEMesh;
class DataField;


/***********************************************************************************************
 * \brief This class is computing coupling matrices to perform the Barycentric method between a NURBS geometry and a polygon mesh
 * ***********/
class IGABarycentricMapper: public AbstractMapper {
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

    /// The reverse element freedom table for the fluid mesh
    std::map<int, std::vector<int> > meshFENodeToElementTable;

    /// The IGA to FEM coupling matrix
    MathLibrary::SparseMatrix<double> *C_M;

    /// The FEM to IGA left hand side coupling matrix
    MathLibrary::SparseMatrix<double> *C_L;

    /// The FEM to IGA left hand side transpose coupling matrix
    MathLibrary::SparseMatrix<double> *C_LT;

    /// The FEM to IGA right hand side coupling matrix
    MathLibrary::SparseMatrix<double> *C_R;

    /// The parametric coordinates of the projection of the nodes on the surface
    /// Each element of vector (FE node) is a map from patch index to local patch coordinates
    std::vector<std::map<int, std::vector<double> > > projectedCoords;

    /// The parametric coordinates of the projection of the CPs on the surface
    /// Each element of vector (CP) is a pair of patch index and local patch coordinates
    std::vector<std::map<int, std::vector<double> > > projectedCPs;

    /// The parametric coordinates of the projection of the projectedCPs on the FE mesh
    /// Each element of vector (projectedCP) is a pair of element index and local element coordinates
    std::vector<std::pair<int, std::vector<double> > > projectedCPsOnFEMesh;

    /// The control point net for the whole IGA mesh
    std::vector<IGAControlPoint*> meshIGACpNet;

    std::vector<bool> isProjectionOfCpOnTrimmed;

    int numOfValidCPs;

    /// The number of nodes that the neighboring element to a control point has
    int *numNodesPerNeighborElem;

    /// number of FE neighbors to search
    static const int MAX_NUM_NEIGHBORS_TO_SEARCH;

    /// Flag on the mapping direction
    bool isMappingIGA2FEM;

    // Slave mesh has the known field, Master mesh has the unknown field
    size_t numNodesSlave;
    size_t numNodesMaster;

    struct nonlinearSchemeProperties {
        int maxNumOfIterations;
        double tolerance;
    } newtonRaphson;
    struct projectionProperties {
        double maxProjectionDistance;
        int numRefinementForIntialGuess;
        double maxDistanceForProjectedPointsOnDifferentPatches;
    } projectionProperties;

public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _name The name of the mapper
     * \param[in] _meshIGA The IGAMesh
     * \param[in] _meshFE The FEMesh
     * \param[in] _isMappingIGA2FEM Boolean; true if mapping is from IGA to FEM
     * \author Fabien Pean, Chenshen Wu, Apostolos Petalas
     ***********/
    IGABarycentricMapper(std::string _name, IGAMesh *_meshIGA, FEMesh *_meshFE, bool _isMappingIGA2FEM);

    /***********************************************************************************************
     * \brief Destructor
     * \author Fabien Pean, Chenshen Wu, Apostolos Petalas
     ***********/
    virtual ~IGABarycentricMapper();

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
     * \brief Build the coupling matrix C_M
     * \author Fabien Pean
     ***********/
    void buildCouplingMatrices();

    /***********************************************************************************************
     * \brief Perform consistent mapping from IGA to FE (map displacements)
     * \param[in] fieldIGA is the input data
     * \param[out] fieldFE is the output data
     * \author Chenshen Wu, Apostolos Petalas
     ***********/
    void consistentMapping(const double *fieldIGA, double *fieldFE);

    /***********************************************************************************************
     * \brief Perform conservative mapping from FE to IGA (map forces)
     * \param[in] fieldFE is the input data
     * \param[out] fieldIGA is the output data
     * \author Chenshen Wu, Apostolos Petalas
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
     * \brief Fills up the array projectedCPs by performing closest point projection
     * \author Chenshen Wu, Apostolos Petalas
     ***********/
    void projectCPsToSurface();

    /***********************************************************************************************
     * \brief Compute the initial guess for a CP
     * \param[in] _patchCount The index of the patch we are working on
     * \param[in] _cpIndex The index of the CP we are working with
     * \param[out] _u The output guess in u direction
     * \param[out] _v The output guess in v direction
     * \author Fabien Pean, Apostolos Petalas
     ***********/
    void computeInitialGuessForProjectionOfCPs(const int _patchIndex, const int _cpIndex, double& _u, double& _v);
    
    /***********************************************************************************************
     * \brief Compute the projection of a control point on a patch using Newton-Raphson
     * \param[in] _patchIndex The index of the patch we are working on
     * \param[in] _cpIndex The global index of the control point we are working with
     * \param[in] _u0 The initial guess in u direction
     * \param[in] _v0 The initial guess in v direction
     * \param[out] _minProjectionDistance The previous distance computed
     * \param[out] _minProjectionPoint The previous point computed
     * \author Fabien Pean, Apostolos Petalas
     ***********/
    bool projectCpOnPatch(const int _patchIndex, const int _cpIndex, const double _u0, const double _v0, double& _minProjectionDistance, std::vector<double>& _minProjectionPoint);

    /***********************************************************************************************
     * \brief Compute the knot span of a control point's shape function's support from the knot vector
     * \param[in] _thePatch The patch we are working on
     * \param[in] _cpIndex The index of the CP we are working with
     * \param[out] _u1 The start local coordinate in u direction
     * \param[out] _u2 The end local coordinate in u direction
     * \param[out] _v1 The start local coordinate in v direction
     * \param[out] _v2 The end local coordinate in v direction
     * \author Apostolos Petalas
     ***********/
    void knotSpanOfSupportCP(const int _patchIndex, const int _cpIndex, double& _u1, double& _u2, double& _v1, double& _v2);

    /***********************************************************************************************
     * \brief Fills up the array projectedCPsOnFEMesh by computing the projection on the FE mesh of the projections of the CPs on the IGA surface
     * \author Tianyang Wang, Apostolos Petalas
     ***********/
    void computeNeighborsAndWeights();

    /***********************************************************************************************
     * \brief Given the element index/id, return the element node cartesian coordinates
     * \param[in] elemIndex The element index/id
     * \param[out] elem The element node cartesian coordinates corresponding to the element index
     * \author Tianyang Wang
     ***********/
    void getElemCoorInFEMMesh(int elemIndex, double *elem);

    /***********************************************************************************************
     * \brief Determine whether a point given in local coordinates is inside the element or not
     * \param[in] numNodesThisElem Number of nodes of this element (3 or 4)
     * \param[out] localCoor Local coordinates of the point
     * \author Tianyang Wang
     ***********/
    bool insideElement(int numNodesThisElem, double *localCoor);

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
     * \param[in] _u0 The initial guess in u direction
     * \param[in] _v0 The initial guess in v direction
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
     * \brief Asserts if a point is not projected on a trimmed region
     * \param[in] _nodeIndex The global index of the patch on which the point is projected
     * \param[in] _u The local coordinate of the projection in the u direction
     * \param[in] _v The local coordinate of the projection in the v direction
     * \author Apostolos Petalas
     ***********/
    bool projectionInside(const int _patchIndex, const double _u, const double _v);
    
    /***********************************************************************************************
     * \brief Clip the trimmed regions off the patch
     * \param[in] _thePatch The patch that is to be clipped
     * \param[out] _listPolygon A set of polygons after application of clipping patch - the trimmed patch
     * \author Apostolos Petalas
     ***********/
    void clipPatch(const IGAPatchSurface* _thePatch, ListPolygon2D& _listPolygon);
    
    /***********************************************************************************************
     * \brief Compute matrices C_M, C_L and C_R by looping over the FE elements or over the control points and processing them
     * \author Apostolos Petalas
     ***********/
    void computeCouplingMatrices();

    /// Writing output functions
public:
    /***********************************************************************************************
     * \brief Writes the projected nodes of the FE mesh onto the IGA surface into a file
     * \author Andreas Apostolatos
     ***********/
    void writeProjectedNodesOntoIGAMesh();

    /***********************************************************************************************
     * \brief Writes the projected control points of the IGA surface onto the FE mesh into a file
     * \author Andreas Apostolatos, Apostolos Petalas
     ***********/
    void writeProjectedCPsOntoFEMMesh();
    
    /***********************************************************************************************
     * \brief Print coupling matrices C_M, C_L and C_R in file in csv format with space delimiter
     * \author Fabien Pean
     ***********/
    void writeCouplingMatricesToFile();

    /// Debugging functions
public:
    /***********************************************************************************************
     * \brief Print coupling matrices C_M, C_L and C_R
     * \author Chenshen Wu
     ***********/
    void printCouplingMatrices();
    
    /***********************************************************************************************
     * \brief Check consistency of the coupling, constant field gives constant field
     * \author Fabien Pean
     ***********/
    void checkConsistency();

    /// unit test class
    friend class TestIGABarycentricMapperCurve;
    friend class TestIGABarycentricMapperSemisphere;

};
}

#endif /* IGABARYCENTRICMAPPER_H_ */
