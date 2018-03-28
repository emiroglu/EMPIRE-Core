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

#include "DataFieldIntegrationNURBS.h"
#include "IGAMesh.h"
#include "IGAPatchSurface.h"
#include "MathLibrary.h"
#include "ClipperAdapter.h"
#include "TriangulatorAdaptor.h"
#include "Message.h"

namespace EMPIRE {

double DataFieldIntegrationNURBS::EPS_CLEANTRIANGLE = 1e-6;
double DataFieldIntegrationNURBS::EPS_CLIPPING = 1e-9;

using std::make_pair;

DataFieldIntegrationNURBS::DataFieldIntegrationNURBS(IGAMesh* _mesh): meshIGA(_mesh) {

    // Get the number of patches
    int numPatches = meshIGA->getNumPatches();

    // Get the number of control points
    numNodes = meshIGA->getNumNodes();

    // Initialize the mass matrix
    massMatrix = new EMPIRE::MathLibrary::SparseMatrix<double>(numNodes,false);

    // Initialize the Gauss quadrature rules
    gaussRuleOnTriangle = new EMPIRE::MathLibrary::IGAGaussQuadrature*[numPatches];
    gaussRuleOnQuadrilateral = new EMPIRE::MathLibrary::IGAGaussQuadrature*[numPatches];

    // Initialize the integration area
    areaIntegration = 0.0;

    // Create the Gauss quadrature rules for all patches
    createGaussQuadratureRules();

    // Compute the mass matrix
    computeMassMatrix();

    // Print the integration area
    INFO_OUT() << "The integration area in the data integration filter is equal to: " << areaIntegration << std::endl;

    // Enforce flying nodes in Cnn
    enforceCnn();
}

DataFieldIntegrationNURBS::~DataFieldIntegrationNURBS() {
    // Initialize auxiliary arrays
    int numPatches = meshIGA->getNumPatches();

    // Delete the quadrature rules
    for (int iPatches = 0; iPatches < numPatches; iPatches++) {
        delete[] gaussRuleOnTriangle[iPatches];
        delete[] gaussRuleOnQuadrilateral[iPatches];
    }
    delete[] gaussRuleOnTriangle;
    delete[] gaussRuleOnQuadrilateral;
}

void DataFieldIntegrationNURBS::createGaussQuadratureRules() {
    /*
     * Creates a Gauss quadrature rule for quadrilaterals and for triangles at each patch
     */

    // Initialize variables
    int pDegree;
    int qDegree;
    int polOrder;
    int pTilde;
    int numGPs;
    int numGPsPerPolOrder[8] = {1, 3, 4, 6, 7, 12, 13, 16};

    // Number of patches
    int numPatches = getIGAMesh()->getNumPatches();

    // Loop over all patches
    for (int iPatches = 0; iPatches < numPatches; iPatches++) {
        // Get the polynomial order of the patches
        pDegree = meshIGA->getSurfacePatch(iPatches)->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree();
        qDegree = meshIGA->getSurfacePatch(iPatches)->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree();

        // Find the polynomial degree of the integrand \int_{\Omega} R_i*R_i d\Omega when a triangle is considered
        polOrder = 2*(pDegree + qDegree);
        if (polOrder > 8)
            pTilde = 2*std::max(pDegree, qDegree);

        // Find the number of Gauss points when a triangle is considered
        if (polOrder <= 8) { // Use the Gauss quadrature over the triangle with the symmetric rule
            if (polOrder == 0)
                numGPs = numGPsPerPolOrder[0];
            else
                numGPs = numGPsPerPolOrder[polOrder - 1];
        } else // Use the Gauss quadrature over the triangle with the degenerated quadrilateral
            numGPs = pow(std::ceil((pTilde + 1)/2.0), 2.0);

        // Instantiate the corresponding Gauss quadrature on triangle
        if (polOrder <= 8) // Use the Gauss quadrature over the triangle with the symmetric rule
            gaussRuleOnTriangle[iPatches] = new MathLibrary::IGAGaussQuadratureOnTriangle(numGPs);
        else // Use the Gauss quadrature over the triangle with the degenerated quadrilateral
            gaussRuleOnTriangle[iPatches] = new MathLibrary::IGAGaussQuadratureOnTriangleUsingDegeneratedQuadrilateral(numGPs);

        // Find the polynomial degree of the integrand \int_{\Omega} R_i*R_i d\Omega when a quadrilateral is considered
        pTilde = 2*std::max(pDegree, qDegree);

        // Find the number of Gauss points when a quadrilateral is considered
        numGPs = pow(std::ceil((pTilde + 1)/2.0), 2.0);

        // Instantiate the corresponding Gauss quadrature on quadrilateral
        gaussRuleOnQuadrilateral[iPatches] = new MathLibrary::IGAGaussQuadratureOnBiunitQuadrilateral(numGPs);
    }
}

void DataFieldIntegrationNURBS::computeMassMatrix() {
    /*
     * Compute the mass matrix of the multipatch B-Spline geometry
     *
     * Function layout :
     *
     * 1. Loop over all patches
     * ->
     *    1i. Get the surface patch
     *   1ii. Get the boundaries of the parameter space
     *  1iii. Clip the parameter space by the knot spans
     *   1iv. Loop over all the generated polygons
     *   ->
     *        1iv.1. Clean the polygon
     *        1iv.2. Check if the polygon has less than 3 vertices and if yes continue
     *        1iv.3. Initialize list of trimmed polygons in case patch is not trimmed
     *        1iv.4. Clip the polygon by trimming if the patch is trimmed
     *        1iv.5. Loop over all generated subpolygons
     *        ->
     *               1iv.5i. Clean the polygon
     *              1iv.5ii. Check if the polygon has less than 3 vertices and if yes continue
     *             1iv.5iii. Triangulate the generated polygon
     *              1iv.5iv. Loop over all triangles of the triangulation
     *              ->
     *                       1iv.5iv.1. Clean triangle
     *                       1iv.5iv.2. Check if the triangle after cleaning is no more a triangle
     *                       1iv.5iv.3. Pass the triangle into the integration function
     *              <-
     *        <-
     * <-
     */

    // Initialize auxiliary arrays
    int numPatches = meshIGA->getNumPatches();

    // 1. Loop over all patches
    for(unsigned int iPatches = 0; iPatches < numPatches; iPatches++) {
        // 1i. Get the surface patch
        IGAPatchSurface* patch = meshIGA->getSurfacePatch(iPatches);

        // 1ii. Get the boundaries of the parameter space
        Polygon2D polygonUV(4);
        {
            const double u0 = patch->getIGABasis()->getUBSplineBasis1D()->getFirstKnot();
            const double v0 = patch->getIGABasis()->getVBSplineBasis1D()->getFirstKnot();
            const double u1 = patch->getIGABasis()->getUBSplineBasis1D()->getLastKnot();
            const double v1 = patch->getIGABasis()->getVBSplineBasis1D()->getLastKnot();
            polygonUV[0] = make_pair(u0,v0);
            polygonUV[1] = make_pair(u1,v0);
            polygonUV[2] = make_pair(u1,v1);
            polygonUV[3] = make_pair(u0,v1);
        }

        // 1iii. Clip the parameter space by the knot spans
        Polygon2D listSpan;
        ListPolygon2D listKnotPolygonUV;
        clipByKnotSpan(patch,polygonUV,listKnotPolygonUV,listSpan);

        // 1iv. Loop over all the generated polygons
        for(int index = 0;index < listSpan.size(); index++) {
            // 1iv.1. Clean the polygon
            ClipperAdapter::cleanPolygon(listKnotPolygonUV[index]);

            // 1iv.2. Check if the polygon has less than 3 vertices and if yes continue
            if(listKnotPolygonUV[index].size() < 3)
                continue;

            // 1iv.3. Initialize list of trimmed polygons in case patch is not trimmed
            ListPolygon2D listTrimmedPolygonUV(1, listKnotPolygonUV[index]);

            // 1iv.4. Clip the polygon by trimming if the patch is trimmed
            if(patch->isTrimmed())
                clipByTrimming(patch,listKnotPolygonUV[index], listTrimmedPolygonUV);

            // 1iv.5. Loop over all generated subpolygons
            for(int trimmedPolygonIndex = 0; trimmedPolygonIndex < listTrimmedPolygonUV.size(); trimmedPolygonIndex++) {
                // 1iv.5i. Clean the polygon
                ClipperAdapter::cleanPolygon(listTrimmedPolygonUV[trimmedPolygonIndex]);

                // 1iv.5ii. Check if the polygon has less than 3 vertices and if yes continue
                if(listTrimmedPolygonUV[trimmedPolygonIndex].size()<3)
                    continue;

                // 1iv.5iii. Triangulate the generated polygon
                ListPolygon2D triangulatedPolygons = triangulatePolygon(listTrimmedPolygonUV[trimmedPolygonIndex]);

                // 1iv.5iv. Loop over all triangles of the triangulation
                for(ListPolygon2D::iterator triangulatedPolygon = triangulatedPolygons.begin(); triangulatedPolygon != triangulatedPolygons.end(); triangulatedPolygon++) {
                    // 1iv.5iv.1. Clean triangle
                    ClipperAdapter::cleanPolygon(*triangulatedPolygon,EPS_CLEANTRIANGLE);

                    // 1iv.5iv.2. Check if the triangle after cleaning is no more a triangle
                    if(triangulatedPolygon->size() < 3)
                        continue;

                    // 1iv.5iv.3. Pass the triangle into the integration function
                    integrate(patch, iPatches, *triangulatedPolygon, listSpan[index].first, listSpan[index].second);
                }
            }
        }
    }
}

void DataFieldIntegrationNURBS::clipByTrimming(const IGAPatchSurface* _thePatch, const Polygon2D& _polygonUV, ListPolygon2D& _listPolygonUV) {
    /*
     * Clips a given polygon by the tirmming curves within a patch
     *
     * Function layout:
     *
     * 1. Initialize auxiliary variables
     *
     * 2. Loop over all trimming loops within the patch
     * ->
     *    2i. Get the linearized trimming loop
     *   2ii. Add the trimming loop in the clipper
     * <-
     *
     * 3. Setup filling rule to have for sure clockwise loop as hole and counterclockwise as boundaries
     */

    // 1. Initialize auxiliary variables
	ClipperAdapter c;

    // 2. Loop over all trimming loops within the patch
    for(int iLoop = 0; iLoop < _thePatch->getTrimming().getNumOfLoops(); iLoop++) {
        // 2i. Get the linearized trimming loop
        const std::vector<double> clippingWindow = _thePatch->getTrimming().getLoop(iLoop).getPolylines();

        // 2ii. Add the trimming loop in the clipper
		c.addPathClipper(clippingWindow);
	}

    // 3. Setup filling rule to have for sure clockwise loop as hole and counterclockwise as boundaries
	c.setFilling(ClipperAdapter::POSITIVE, 0);
	c.addPathSubject(_polygonUV);
	c.clip();
	c.getSolution(_listPolygonUV);
}

bool DataFieldIntegrationNURBS::computeKnotSpanOfProjElement(const IGAPatchSurface* _thePatch, const Polygon2D& _polygonUV, int* _span) {
    /*
     * Clips a given polygon by the knot spans in the given patch
     *
     * Function layout:
     *
     * 1. Find minimum and maximum knot span indices in u- and v-directions
     *
     * 2. Loop over all nodes in the given polygon
     * ->
     *    2i. Find the knot span indices where the nodes belong into
     *   2ii. Clamp the knot span indices in case the image of the node is outside the knot limits
     * <-
     *
     * 3. Flag on whether all nodes of the polygon are in the same knot span
     *
     * 4. Return the flag
     */

    // 1. Find minimum and maximum knot span indices in u- and v-directions
    int minSpanU = _thePatch->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(_polygonUV[0].first);
    int minSpanV = _thePatch->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(_polygonUV[0].second);
    int maxSpanU = minSpanU;
    int maxSpanV = minSpanV;

    // 2. Loop over all nodes in the given polygon
    for (int iNodes = 1; iNodes < _polygonUV.size(); iNodes++) {
        // 2i. Find the knot span indices where the nodes belong into
        int spanU = _thePatch->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(_polygonUV[iNodes].first);
        int spanV = _thePatch->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(_polygonUV[iNodes].second);

        // 2ii. Clamp the knot span indices in case the image of the node is outside the knot limits
        if (spanU < minSpanU)
            minSpanU = spanU;
        if (spanU > maxSpanU)
            maxSpanU = spanU;
        if (spanV < minSpanV)
            minSpanV = spanV;
        if (spanV > maxSpanV)
            maxSpanV = spanV;
    }

    // 3. Flag on whether all nodes of the polygon are in the same knot span
    bool OnSameKnotSpan = (minSpanU == maxSpanU && minSpanV == maxSpanV);
    if(_span!=NULL) {
        _span[0] = minSpanU;
        _span[1] = maxSpanU;
        _span[2] = minSpanV;
        _span[3] = maxSpanV;
    }

    // 4. Return the flag
    return OnSameKnotSpan;
}

void DataFieldIntegrationNURBS::clipByKnotSpan(const IGAPatchSurface* _thePatch, const Polygon2D& _polygonUV, ListPolygon2D& _listPolygon, Polygon2D& _listSpan) {
    /*
     * Clips a given polygon by the knot spans of the given patch
     *
     * Function layout:
     *
     * 1. Initialize auxiliary variables
     *
     * 2. Get the knot vectors of the patch
     *
     * 3.find the knot span which the current element located in. (from minSpanu to maxSpanu in U-direction, and from minSpanV to max SpanV in V-direction)
     *
     * 4. Find the clipped polygon
     * ->
     *    ### If all vertices of the polygon are on same knot span then returned the same polygon as input ###
     *    4i. Add the polygon into the list
     *   4ii. Add the corresponding knot span indices of the knot span containing the polygon
     *    ### Else clip the polygon for every knot span window it is crossing ###
     *    4i. Loop over all spans in u-direction
     *    ->
     *        4i.1. Loop over all spans in v-direction
     *        ->
     *              4i.1i. Initialize a clipper adapter
     *             4i.1ii. Clip the polygon with the knot span window
     *             ### Check if the encountered knot span is collapsed ###
     *             4i.1ii.1. Create the knot span window
     *             4i.1ii.2. Clip the polygon with the knot span window (WARNING : here we assume to get only a single output polygon from the clipping!)
     *        <-
     *    <-
     * <-
     */

    // 1. Initialize auxiliary variables
    int span[4];

    // 2. Get the knot vectors of the patch
    const double *knotVectorU = _thePatch->getIGABasis()->getUBSplineBasis1D()->getKnotVector();
    const double *knotVectorV = _thePatch->getIGABasis()->getVBSplineBasis1D()->getKnotVector();

    // 3.find the knot span which the current element located in. (from minSpanu to maxSpanu in U-direction, and from minSpanV to max SpanV in V-direction)
	int isOnSameKnotSpan = computeKnotSpanOfProjElement(_thePatch, _polygonUV,span);
    int minSpanU = span[0];
    int maxSpanU = span[1];
    int minSpanV = span[2];
    int maxSpanV = span[3];

    // 4. Find the clipped polygon
    if (isOnSameKnotSpan) { // ### If all vertices of the polygon are on same knot span then returned the same polygon as input ###
        // 4i. Add the polygon into the list
		_listPolygon.push_back(_polygonUV);

        // 4ii. Add the corresponding knot span indices of the knot span containing the polygon
		_listSpan.push_back(make_pair(minSpanU,minSpanV));
    } else { // ### Else clip the polygon for every knot span window it is crossing ###
        // 4i. Loop over all spans in u-direction
		for (int spanU = minSpanU; spanU <= maxSpanU; spanU++) {
            // 4i.1. Loop over all spans in v-direction
			for (int spanV = minSpanV; spanV <= maxSpanV; spanV++) {
                // 4i.1i. Initialize a clipper adapter
                ClipperAdapter c(EPS_CLIPPING);

                // 4i.1ii. Clip the polygon with the knot span window
                if (knotVectorU[spanU] != knotVectorU[spanU + 1] && knotVectorV[spanV] != knotVectorV[spanV + 1]) { // ### Check if the encountered knot span is collapsed ###
                    // 4i.1ii.1. Create the knot span window
					Polygon2D knotSpanWindow(4);
                    knotSpanWindow[0] = make_pair(knotVectorU[spanU],knotVectorV[spanV]);
                    knotSpanWindow[1] = make_pair(knotVectorU[spanU + 1],knotVectorV[spanV]);
                    knotSpanWindow[2] = make_pair(knotVectorU[spanU + 1],knotVectorV[spanV + 1]);
                    knotSpanWindow[3] = make_pair(knotVectorU[spanU],knotVectorV[spanV + 1]);

                    // 4i.1ii.2. Clip the polygon with the knot span window (WARNING : here we assume to get only a single output polygon from the clipping!)
					Polygon2D solution = c.clip(_polygonUV,knotSpanWindow);
					/// Store polygon and its knot span for integration
					_listPolygon.push_back(solution);
					_listSpan.push_back(make_pair(spanU,spanV));
				}
			}
		}
	}
}

DataFieldIntegrationNURBS::ListPolygon2D DataFieldIntegrationNURBS::triangulatePolygon(const Polygon2D& _polygonUV) {
    /*
     * Triangulates given polygon
     *
     * Function layout :
     *
     * 1. If the polygon is already a triangle do nothing
     *
     * 2. Initialize triangulator
     *
     * 3. Loop over all segments of the given polygon and add its ends into the triangulator
     *
     * 4. Triangulate the polygon
     *
     * 5. Fill the output list of the newly generated polygons
     *
     * 6. Return the updated list of the triangulated polygons
     */

    // 1. If the polygon is already a triangle do nothing
    if(_polygonUV.size() < 4)
		return ListPolygon2D(1,_polygonUV);

    // 2. Initialize triangulator
	TriangulatorAdaptor triangulator;

    // 3. Loop over all segments of the given polygon and add its ends into the triangulator
    for(Polygon2D::const_iterator it = _polygonUV.begin(); it!=_polygonUV.end(); it++)
		triangulator.addPoint(it->first,it->second,0);

    // 4. Triangulate the polygon
	int numTriangles = _polygonUV.size() - 2;
	int triangleIndexes[3 * numTriangles];
    bool triangulated = triangulator.triangulate(triangleIndexes);
	if(!triangulated)
		return ListPolygon2D();

    // 5. Fill the output list of the newly generated polygons
	ListPolygon2D out(numTriangles, Polygon2D(3));
	for(int i = 0; i < numTriangles; i++)
        for(int j = 0; j < 3; j++)
            out[i][j] = _polygonUV[triangleIndexes[3*i + j]];

    // 6. Return the updated list of the triangulated polygons
	return out;
}

void DataFieldIntegrationNURBS::integrate(IGAPatchSurface* _thePatch, int _indexPatch, Polygon2D _polygonUV, int _spanU, int _spanV) {
    /*
     * Integrates the triangulated polygons and assembles the mass matrix
     *
     * Function layout :
     *
     * 1. Check input
     *
     * 2. Get the corresponding quadrature rule depending on the integration domain
     *
     * 3. Initialize auxiliary variables
     *
     * 4. Create an element freedom table
     *
     * 5. Copy input polygon into contiguous C format
     *
     * 6. Loop over all Gauss points
     * ->
     *    6i. Get the Gauss point coordinates in the integration space
     *   6ii. Get the Gauss point weight
     *  6iii. Compute the basis functions describing parametrically the integration domain
     *   6iv. Compute the image of the gauss point in the NURBS parameter space
     *    6v. Compute the local basis functions and their first derivatives at the Gauss point
     *   6vi. Compute the base vectors at the Gauss point
     *  6vii. Compute the surface normal vector at the Gauss point
     * 6viii. Compute the determinant of the Jacobian of the transformation from the physical space to the NURBS parameter space
     *   6ix. Compute the Jacobian of the transformation from the NURBS parameter space to the integration space
     *    6x. Compute the Jacobian products
     *   6xi. Update the integration area at the Gauss point
     *  6xii. Loop over all local basis functions in a nested loop to compute and assemble the local mass matrix
     *   ->
     *        6xii.1. Compute the local basis functions of the dual product RI*RJ
     *        6xii.2. Compute and add the local mass matrix contributions
     *   <-
     * <-
     */

    // 1. Check input
    assert(!_polygonUV.empty());
    int numNodesUV = _polygonUV.size();
    assert(numNodesUV > 2);
    assert(numNodesUV < 5);

    // 2. Get the corresponding quadrature rule depending on the integration domain
    EMPIRE::MathLibrary::IGAGaussQuadrature* theGaussQuadrature;
    int numNodesQuadrature = numNodesUV;
    if (numNodesQuadrature == 3)
        theGaussQuadrature = gaussRuleOnTriangle[_indexPatch];
    else if (numNodesQuadrature == 4)
        theGaussQuadrature = gaussRuleOnQuadrilateral[_indexPatch];
    else {
        ERROR_OUT() << "Only triangles and quadrilaterals are expected as integration domains";
        exit(-1);
    }

    // 3. Initialize auxiliary variables
    int indexBaseVctU, indexBaseVctV;
    int derivDegree = 1;
    int derivDegreeBaseVec = 0;
    int noBaseVec = 2;
    int noCoord = 3;
    int pDegree = _thePatch->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree();
    int qDegree = _thePatch->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree();
    int numBasisFunctions = (pDegree + 1) * (qDegree + 1);
    double shapeFuncs[numNodesQuadrature];
    double uv[2];
    double baseVectorU[3];
    double baseVectorV[3];
    double surfaceNormalTilde[3];
    double gaussWeight;
    double IGABasisFctsI;
    double IGABasisFctsJ;
    double JacobianUVToPhysical;
    double JacobianCanonicalToUV;
    double Jacobian;
    double dudx;
    double dudy;
    double dvdx;
    double dvdy;
    double localBasisFunctionsAndDerivatives[(derivDegree + 1) * (derivDegree + 2) * numBasisFunctions / 2];
    double baseVctsAndDerivs[(derivDegreeBaseVec + 1) * (derivDegreeBaseVec + 2) * noCoord * noBaseVec / 2];

    // 4. Create an element freedom table
    int dofIGA[numBasisFunctions];
    _thePatch->getIGABasis()->getBasisFunctionsIndex(_spanU, _spanV, dofIGA);
    for (int i = 0; i < numBasisFunctions; i++)
        dofIGA[i] = _thePatch->getControlPointNet()[dofIGA[i]]->getDofIndex();

    // 5. Copy input polygon into contiguous C format
	double nodesUV[8];
	for (int i = 0; i < numNodesUV; i++) {
        nodesUV[i*2] = _polygonUV[i].first;
        nodesUV[i*2+1] = _polygonUV[i].second;
	}

    // 6. Loop over all Gauss points
    for (int iGP = 0; iGP < theGaussQuadrature->getNumGaussPoints(); iGP++) {
        // 6i. Get the Gauss point coordinates in the integration space
        const double *gaussPoint = theGaussQuadrature->getGaussPoint(iGP);

        // 6ii. Get the Gauss point weight
        gaussWeight = theGaussQuadrature->getGaussWeight(iGP);

        // 6iii. Compute the basis functions describing parametrically the integration domain
        MathLibrary::computeLowOrderShapeFunc(numNodesQuadrature, gaussPoint, shapeFuncs);

        // 6iv. Compute the image of the gauss point in the NURBS parameter space
        MathLibrary::computeLinearCombination(numNodesQuadrature, 2, nodesUV, shapeFuncs, uv);

        // 6v. Compute the local basis functions and their first derivatives at the Gauss point
        _thePatch->getIGABasis()->computeLocalBasisFunctionsAndDerivatives
                (localBasisFunctionsAndDerivatives, derivDegree, uv[0], _spanU, uv[1], _spanV);

        // 6vi. Compute the base vectors at the Gauss point
        _thePatch->computeBaseVectorsAndDerivatives(baseVctsAndDerivs, localBasisFunctionsAndDerivatives, derivDegreeBaseVec,_spanU, _spanV);
        for (int iCoord = 0; iCoord < noCoord; iCoord++) {
            indexBaseVctU = _thePatch->indexDerivativeBaseVector(0 , 0, 0, iCoord, 0);
            baseVectorU[iCoord] = baseVctsAndDerivs[indexBaseVctU];
            indexBaseVctV = _thePatch->indexDerivativeBaseVector(0 , 0, 0, iCoord, 1);
            baseVectorV[iCoord] = baseVctsAndDerivs[indexBaseVctV];
        }

        // 6vii. Compute the surface normal vector at the Gauss point
        MathLibrary::computeVectorCrossProduct(baseVectorU, baseVectorV, surfaceNormalTilde);

        // 6viii. Compute the determinant of the Jacobian of the transformation from the physical space to the NURBS parameter space
        JacobianUVToPhysical = MathLibrary::vector2norm(surfaceNormalTilde, noCoord);

        // 6ix. Compute the Jacobian of the transformation from the NURBS parameter space to the integration space
        if (numNodesQuadrature == 3) {
            JacobianCanonicalToUV = MathLibrary::computeAreaTriangle(nodesUV[2] - nodesUV[0], nodesUV[3] - nodesUV[1], 0,
                                                                     nodesUV[4] - nodesUV[0], nodesUV[5] - nodesUV[1], 0)*2.0;
        } else {
            dudx = .25 * (-(1 - gaussPoint[2]) * nodesUV[0] + (1 - gaussPoint[2]) * nodesUV[2]
                       + (1 + gaussPoint[2]) * nodesUV[4] - (1 + gaussPoint[2]) * nodesUV[6]);
            dudy = .25 * (-(1 - gaussPoint[1]) * nodesUV[0] - (1 + gaussPoint[1]) * nodesUV[2]
                       + (1 + gaussPoint[1]) * nodesUV[4] + (1 - gaussPoint[1]) * nodesUV[6]);
            dvdx = .25 * (-(1 - gaussPoint[2]) * nodesUV[1] + (1 - gaussPoint[2]) * nodesUV[3]
                       + (1 + gaussPoint[2]) * nodesUV[5] - (1 + gaussPoint[2]) * nodesUV[7]);
            dvdy = .25 * (-(1 - gaussPoint[1]) * nodesUV[1] - (1 + gaussPoint[1]) * nodesUV[3]
                       + (1 + gaussPoint[1]) * nodesUV[5] + (1 - gaussPoint[1]) * nodesUV[7]);
            JacobianCanonicalToUV = fabs(dudx * dvdy - dudy * dvdx);
		}

        // 6x. Compute the Jacobian products
        Jacobian = JacobianUVToPhysical * JacobianCanonicalToUV;

        // 6xi. Update the integration area at the Gauss point
        areaIntegration += Jacobian * gaussWeight;

        // 6xii. Loop over all local basis functions in a nested loop to compute and assemble the local mass matrix
        for (int i = 0; i < numBasisFunctions; i++) {
            for (int j = i; j < numBasisFunctions; j++) { // Starts from i because of computing only the upper triangular entries of the matrix
                    // 6xii.1. Compute the local basis functions of the dual product RI*RJ
                    IGABasisFctsI = localBasisFunctionsAndDerivatives[_thePatch->getIGABasis()->indexDerivativeBasisFunction(1, 0, 0, i)];
                    IGABasisFctsJ = localBasisFunctionsAndDerivatives[_thePatch->getIGABasis()->indexDerivativeBasisFunction(1, 0, 0, j)];

                    // 6xii.2. Compute and add the local mass matrix contributions
                    (*massMatrix)(dofIGA[i], dofIGA[j]) += IGABasisFctsI * IGABasisFctsJ * Jacobian * gaussWeight;
                    if(dofIGA[i] != dofIGA[j]) // Because matrix not instantiated as symmetric
                        (*massMatrix)(dofIGA[j], dofIGA[i]) += IGABasisFctsI * IGABasisFctsJ * Jacobian * gaussWeight;
			}
		}
    }
}

void DataFieldIntegrationNURBS::enforceCnn() {
    /*
     * Fixes the flying nodes to zero value
     */
    for(int i = 0; i < numNodes; i++) {
		if(massMatrix->isRowEmpty(i))
            (*massMatrix)(i,i) = 1;
	}
}

} /* namespace EMPIRE */
