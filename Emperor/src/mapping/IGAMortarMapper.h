/*  Copyright &copy; 2013, TU Muenchen, Chair of Structural Analysis,
 *  Stefan Sicklinger, Tianyang Wang, Andreas Apostolatos, Munich
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
    std::vector<std::map<int, double*> > *projectedCoords;
//    std::vector<std::map<int, bool> > *isProjectionOrthogonal;

    /// Number of division made for linearization on boundary
    int numDivision;
/// Tolerance up to which projection is trusted
    double disTol;

    /// number of Gauss points used for computing triangle element
    int numGPsTri;

    /// number of Gauss points used for computing quad element
    int numGPsQuad;

    /// Flag on the mapping direction
    bool isMappingIGA2FEM;

    size_t numNodesSlave;
    size_t numNodesMaster;

public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _name The name of the mapper
     * \param[in] _meshIGA The IGAMesh
     * \param[in] _meshFE The FEMesh
     * \param[in] _disTol Tolerance up to which projection is trusted
     * \param[in] _numGPsTri The number of Gauss points used for computing triangle element
     * \param[in] _numGPsQuad The number of Gauss points used for computing quad element
     * \author Chenshen Wu
     ***********/
    IGAMortarMapper(std::string _name, IGAMesh *_meshIGA, FEMesh *_meshFE, double _disTol,
            int _numGPsTri, int _numGPsQuad, bool _isMappingIGA2FEM, int _numDivision=3);

    /***********************************************************************************************
     * \brief Destructor Chenshen Wu
     ***********/
    virtual ~IGAMortarMapper();

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
     * \author Chenshen Wu
     ***********/
    void computeCouplingMatrices();

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

    /// Writing output functions
public:
    /***********************************************************************************************
     * \brief Writes the projected nodes of the FE mesh onto the IGA surface into a file
     * \author Andreas Apostolatos
     ***********/
    void writeProjectedNodesOntoIGAMesh();
    /***********************************************************************************************
     * \brief Print both coupling matrices C_NN and C_NR in file in csv format with space delimiter
     * \author Fabien Pean
     ***********/
    void printCouplingMatricesToFile();

    /// Debugging functions
public:
    void debugPolygon(const Polygon2D& _polygon, std::string _name="");
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

    /// Number of refined parametric locations where to find the candidate closest points for the projection
    static const int REFINED_NUM_PARAMETRIC_LOCATIONS = 10;
};
}

#endif /* IGAMORTARMAPPER_H_ */
