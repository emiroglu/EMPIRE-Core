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

// Inclusion of standard libraries
#include <iostream>
#include <assert.h>
#include <math.h>
#include <algorithm>
#include <string>

// Inclusion of user defined libraries
#include "WeakIGADirichletSurfaceCondition.h"
#include "IGAPatchCurve.h"
#include "IGAPatchSurface.h"
#include "ClipperAdapter.h"
#include "TriangulatorAdaptor.h"
#include "MathLibrary.h"
#include "Message.h"

using namespace std;

namespace EMPIRE {

WeakIGADirichletSurfaceCondition::WeakIGADirichletSurfaceCondition(int _ID,
                                 int _patchIndex, IGAPatchSurfaceTrimmingLoop* _conditionBoundaryLoop) :
        AbstractCondition(_ID),
        patchIndex(_patchIndex), conditionBoundaryLoop(_conditionBoundaryLoop) {

    type = EMPIRE_WeakIGADirichletSurfaceCondition;

    isBoundaryLoop = false;

    // Linearize the boundary loop with the default algorithm (combined algorithm)
    conditionBoundaryLoop->linearize();

    isGPDataInitialized = false;

}

WeakIGADirichletSurfaceCondition::WeakIGADirichletSurfaceCondition(int _ID,
                                 int _patchIndex, int _patchBLIndex) :
        AbstractCondition(_ID),
        patchIndex(_patchIndex), patchBLIndex(_patchBLIndex) {

    type = EMPIRE_WeakIGADirichletSurfaceCondition;

    isBoundaryLoop = true;

    isGPDataInitialized = false;

}

void WeakIGADirichletSurfaceCondition::createGPData(IGAPatchSurface* _patch) {

    if (isGPDataInitialized) assert(false);

    // Initialize coordinates
    int noCoordParam = 2;

    /// Create a quadrature rule for the integration triangles
    // Get the patch polynomial degrees
    int p = _patch->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree();
    int q = _patch->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree();
    int numGPonTria = p + q + 1;

    // Gauss quadrature for the triangles
    WARNING_OUT("In \"WeakIGADirichletSurfaceCondition::createGPData\", using a default number of 6 GPs per triangle.");
    MathLibrary::IGAGaussQuadratureOnTriangle* gaussTriangle = new MathLibrary::IGAGaussQuadratureOnTriangle(6);

    /// Make polygons out of knot spans
    // Get the number of knots in each direction
    int noUKnots = _patch->getIGABasis()->getUBSplineBasis1D()->getNoKnots();
    int noVKnots = _patch->getIGABasis()->getVBSplineBasis1D()->getNoKnots();

    // Get the knot vectors
    double* uKnotVector = _patch->getIGABasis()->getUBSplineBasis1D()->getKnotVector();
    double* vKnotVector = _patch->getIGABasis()->getVBSplineBasis1D()->getKnotVector();

    // Count the number of knot spans
    int noUSpans = 1;
    for (int iKnot = 0; iKnot < noUKnots; iKnot++)
        if(uKnotVector[iKnot] != uKnotVector[iKnot+1])
            noUSpans++;

    int noVSpans = 1;
    for (int iKnot = 0; iKnot < noVKnots; iKnot++)
        if(vKnotVector[iKnot] != vKnotVector[iKnot+1])
            noVSpans++;

    // Initiate a polygon list of knot span windows
    ListPolygon2D knotSpanPolygonList;
    // Loop over U knots
    for (int iU = 0; iU <= noUSpans; iU++) {
        // Loop over V knots
        for (int iV = 0; iV <= noVSpans; iV++) {
            if (uKnotVector[iU] != uKnotVector[iU + 1] && vKnotVector[iV] != vKnotVector[iV + 1]) {
                // Create the knot span polygon
                Polygon2D knotSpanPolygon(4);
                knotSpanPolygon[0] = make_pair(uKnotVector[iU],vKnotVector[iV]);
                knotSpanPolygon[1] = make_pair(uKnotVector[iU+1],vKnotVector[iV]);
                knotSpanPolygon[2] = make_pair(uKnotVector[iU+1],vKnotVector[iV+1]);
                knotSpanPolygon[3] = make_pair(uKnotVector[iU],vKnotVector[iV+1]);

                // Store the knot span polygon
                knotSpanPolygonList.push_back(knotSpanPolygon);
            }
        }
    }

    /// Create the GP data
    // In case the trimming loop is a boundary loop of a patch set the conditionBoundaryLoop object
    if (isBoundaryLoop)
        conditionBoundaryLoop = &_patch->getTrimming().getLoop(patchBLIndex);

    // Initialize variables
    int derivDegree = 1;
    int uKnotSpan;
    int vKnotSpan;
    int noLocalBasisFunctions = (p + 1) * (q + 1);
    double localBasisFunctionsAndDerivatives[(derivDegree + 1) * (derivDegree + 2)
            * noLocalBasisFunctions / 2];
    double baseVectors[6];
    std::vector<double> tmpGPData; // temporary vector to store the created GP data

    // Loop over knot span windows
    for (int iPolygon = 0; iPolygon < knotSpanPolygonList.size(); iPolygon++) {

        // Clip the knot span windows with the boundary loops of the patch
        ListPolygon2D trimClippedPolygonList;
        clipByTrimming(_patch, knotSpanPolygonList[iPolygon], trimClippedPolygonList);

        for (int iTCW = 0; iTCW < trimClippedPolygonList.size(); iTCW++) {

            // Clip the trimming clipped knot span windows with the bounday loop of the condition
            ListPolygon2D trimCondClippedPolygonList;
            clipByCondition(conditionBoundaryLoop, trimClippedPolygonList[iTCW], trimCondClippedPolygonList);

            // Loop over the each polygon after clipping
            for (int iTCCW = 0; iTCCW < trimCondClippedPolygonList.size(); iTCCW++) {

                // Triangulate each sub polygon
                ListPolygon2D trimCondClippedTrias = triangulatePolygon(trimCondClippedPolygonList[iTCCW]);

                // Loop over each triangle
                for (int iTCP = 0; iTCP < trimCondClippedTrias.size(); iTCP++) {

                    // Copy triangle vertices to C++ type
                    int numNodes = trimCondClippedTrias[iTCP].size();
                    double triaUV[numNodes*noCoordParam];
                    for (int iNode = 0; iNode < numNodes; iNode++) {
                        triaUV[iNode * noCoordParam] = trimCondClippedTrias[iTCP][iNode].first;
                        triaUV[iNode * noCoordParam + 1]=trimCondClippedTrias[iTCP][iNode].second;
                    }

                    // Loop over GPs on the triangle
                    for (int iGP = 0; iGP < gaussTriangle->getNumGaussPoints(); iGP++) {
                        const double* GP = gaussTriangle->getGaussPoint(iGP);
                        double GW = gaussTriangle->getGaussWeight(iGP);

                        // Get shape function values
                        double shapeFuncs[numNodes];
                        MathLibrary::computeLowOrderShapeFunc(numNodes, GP, shapeFuncs);

                        // Compute parametric coordinates (u,v) of the GP
                        double uvGP[2];
                        MathLibrary::computeLinearCombination(numNodes, 2, triaUV, shapeFuncs, uvGP);

                        // Compute base vectors on GP
                        uKnotSpan = _patch->findSpanU(uvGP[0]);
                        vKnotSpan = _patch->findSpanV(uvGP[1]);
                        _patch->getIGABasis()->computeLocalBasisFunctionsAndDerivatives(
                                    localBasisFunctionsAndDerivatives, derivDegree, uvGP[0], uKnotSpan, uvGP[1], vKnotSpan);
                        _patch->computeBaseVectors(baseVectors, localBasisFunctionsAndDerivatives, uKnotSpan, vKnotSpan);

                        // Compute Jacobian on GP
                        double JacobianGP = MathLibrary::computeAreaTriangle(baseVectors[0], baseVectors[1], baseVectors[2],
                                baseVectors[3], baseVectors[4], baseVectors[5]) * 2 * GW;

                        tmpGPData.push_back(uvGP[0]);
                        tmpGPData.push_back(uvGP[1]);
                        tmpGPData.push_back(GW);
                        tmpGPData.push_back(JacobianGP);
                    }
                }
            }
        }
    }

    /// Set GP data
    // Initialize and fill up the object members
    surfaceNumGP = tmpGPData.size()/4;
    surfaceGPs = new double[surfaceNumGP*noCoordParam];
    surfaceGPWeights = new double[surfaceNumGP];
    surfaceGPJacobians = new double[surfaceNumGP];
    for (int iGP = 0; iGP < surfaceNumGP; iGP++) {
        surfaceGPs[iGP*noCoordParam] = tmpGPData[iGP*4];
        surfaceGPs[iGP*noCoordParam+1] = tmpGPData[iGP*4+1];
        surfaceGPWeights[iGP] = tmpGPData[iGP*4 + 2];
        surfaceGPJacobians[iGP] = tmpGPData[iGP*4 + 3];
    }

    // Set the initialized flag to true
    isGPDataInitialized = true;
}

void WeakIGADirichletSurfaceCondition::clipByTrimming(const IGAPatchSurface* _thePatch, const Polygon2D& _polygonUV, ListPolygon2D& _listPolygonUV) {
    ClipperAdapter c;
    // Fill clipper with trimming clipping curves
    for(int loop=0;loop<_thePatch->getTrimming().getNumOfLoops();loop++) {
        const std::vector<double> patchClippingWindow=_thePatch->getTrimming().getLoop(loop).getPolylines();
        c.addPathClipper(patchClippingWindow);
    }
    // Setup filling rule to have for sure clockwise loop as hole and counterclockwise as boundaries
    c.setFilling(ClipperAdapter::POSITIVE, 0);
    c.addPathSubject(_polygonUV);
    c.clip();
    c.getSolution(_listPolygonUV);
}

void WeakIGADirichletSurfaceCondition::clipByCondition(const IGAPatchSurfaceTrimmingLoop* _theTrimmingLoop, const Polygon2D& _polygonUV, ListPolygon2D& _listPolygonUV) {
    ClipperAdapter c;
    // Fill clipper with trimming loop to clip with
    const std::vector<double> conditionClippingWindow = _theTrimmingLoop->getPolylines();
    c.addPathClipper(conditionClippingWindow);
    // Setup filling rule to have for sure clockwise loop as hole and counterclockwise as boundaries
    c.setFilling(ClipperAdapter::NEGATIVE, 0);
    c.addPathSubject(_polygonUV);
    c.clip();
    c.getSolution(_listPolygonUV);
}

WeakIGADirichletSurfaceCondition::ListPolygon2D WeakIGADirichletSurfaceCondition::triangulatePolygon(const Polygon2D& _polygonUV) {
    // If already easily integrable by quadrature rule, do nothing
    if(_polygonUV.size()<4)
        return ListPolygon2D(1,_polygonUV);
    // Otherwise triangulate polygon
    TriangulatorAdaptor triangulator;
    // Fill adapter
    for(Polygon2D::const_iterator it=_polygonUV.begin();it!=_polygonUV.end();it++)
        triangulator.addPoint(it->first,it->second,0);
    // Triangulate
    int numTriangles = _polygonUV.size() - 2;
    int triangleIndexes[3 * numTriangles];
    bool triangulated=triangulator.triangulate(triangleIndexes);
    if(!triangulated)
        return ListPolygon2D();
    // Fill output structure
    ListPolygon2D out(numTriangles, Polygon2D(3));
    for(int i = 0; i < numTriangles; i++)
        for(int j=0;j<3;j++)
            out[i][j]=_polygonUV[triangleIndexes[3*i + j]];
    return out;
}

WeakIGADirichletSurfaceCondition::~WeakIGADirichletSurfaceCondition() {
    if (!isBoundaryLoop)
        delete conditionBoundaryLoop;
    if (isGPDataInitialized) {
        delete surfaceGPs;
        delete surfaceGPWeights;
        delete surfaceGPJacobians;
    }
}

} /* namespace EMPIRE */
