/*  Copyright &copy; 2013, TU Muenchen, Chair of Structural Analysis,
 *  Fabien Pean, Munich
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
 * \file DataFieldIntegrationNURBS.h
 * This file holds the class DataFieldIntegrationNURBS
 * \date 1/4/2015
 **************************************************************************************************/

#ifndef DATAFIELDINTEGRATIONNURBS_H_
#define DATAFIELDINTEGRATIONNURBS_H_

#include "AbstractDataFieldIntegration.h"
#include <vector>
#include <utility>

namespace EMPIRE {

class IGAMesh;
class IGAPatchSurface;

/********//**
 * \brief Class DataFieldIntegrationNURBS is an operator from traction to force or vice versa for NURBS surfaces
 * \author Fabien Pean
 ***********/
class DataFieldIntegrationNURBS : public AbstractDataFieldIntegration {
private:
	/// Type definitions
    typedef std::pair<double,double> Point2D;
    typedef std::vector<Point2D> Polygon2D;
    typedef std::vector<Polygon2D> ListPolygon2D;

    /// The IGA multipatch geometry
    IGAMesh* meshIGA;

    /// Quadrature rule over the triangulated subdomains
    EMPIRE::MathLibrary::IGAGaussQuadrature **gaussRuleOnTriangle;

    /// Quadrature rule over the non-triangulated subdomains
    EMPIRE::MathLibrary::IGAGaussQuadrature **gaussRuleOnQuadrilateral;

public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _mesh Pointer to mesh class relies on.
     * \author Fabien Pean
     ***********/
    DataFieldIntegrationNURBS(IGAMesh* mesh);
    /***********************************************************************************************
     * \brief Destructor
     * \author Tianyang Wang
     ***********/
    virtual ~DataFieldIntegrationNURBS();

private:
    /***********************************************************************************************
     * \brief Compute the mass matrix
     * \author Andreas Apostolatos
     ***********/
    void computeMassMatrix();

    /***********************************************************************************************
     * \brief Create a Gauss quadrature rule for triangles and for quadrilaterals for each patch of the B-Spline multipatch geometry
     * \author Andreas Apostolatos
     ***********/
    void createGaussQuadratureRules();

    /***********************************************************************************************
     * \brief Clip the input polygon by the trimming window of the patch
     * \param[in] _thePatch 	The patch for which trimming curves are used
     * \param[in] _polygonUV 	An input polygon defined in parametric (i.e. 2D) space
     * \param[out] _listPolygon	A set of polygons after application of trimming polygon
     * \author Fabien Pean
     ***********/
    void clipByTrimming(const IGAPatchSurface* _thePatch, const Polygon2D& _polygonUV, ListPolygon2D& _listPolygonUV);

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
     * \brief Integrate the element coupling matrices and assemble them to the global one
     * \param[in] _thePatch	The patch to compute the coupling matrices for
     * \param[in] _indexPatch The index of the patch (ordering of storing the patches)
     * \param[in] _polygonIGA The resulting from the clipping polygon at each knot span in the NURBS space
     * \param[in] _spanU The knot span index in the u-direction where basis will be evaluated
     * \param[in] _spanV The knot span index in the v-direction where basis will be evaluated
     * \author Andreas Apostolatos, Fabien Pean
     ***********/
    void integrate(IGAPatchSurface* _thePatch, int _indexPatch, Polygon2D _polygonUV, int _spanU, int _spanV);

    /***********************************************************************************************
	 * \brief Fill empty rows with identity for ensuring numerical stability
     * \author Andreas Apostolatos, Fabien Pean
	 ***********/
    void enforceCnn();

    /// Get and set functions
public:
    /***********************************************************************************************
     * \brief Get the IGAMesh containing all B-Spline surface patches which form the geometry
     * \author Andreas Apostolatos
     ***********/
    IGAMesh* getIGAMesh() {
        return meshIGA;
    }

    // Constant members of the class
private:
    /// Tolerance for cleaning a triangle before integrating
    static double EPS_CLEANTRIANGLE;
    static double EPS_CLIPPING;
};

} /* namespace EMPIRE */
#endif /* DATAFIELDINTEGRATIONNURBS_H_ */
