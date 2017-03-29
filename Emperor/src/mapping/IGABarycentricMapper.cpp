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

#ifdef FLANN
#include "flann/flann.hpp"
#endif

#include "IGABarycentricMapper.h"
#include "IGAPatchSurface.h"
#include "IGAMesh.h"
#include "FEMesh.h"
#include "ClipperAdapter.h"
#include "TriangulatorAdaptor.h"
#include "MathLibrary.h"
#include "DataField.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <set>
#include <algorithm>

using namespace std;

namespace EMPIRE {

/// Declaration statement
static const string HEADER_DECLARATION = "Author: Andreas Apostolatos";

const int IGABarycentricMapper::MAX_NUM_NEIGHBORS_TO_SEARCH = 10;

IGABarycentricMapper::IGABarycentricMapper(std::string _name, IGAMesh *_meshIGA, FEMesh *_meshFE,
    bool _isMappingIGA2FEM) :
    name(_name), meshIGA(_meshIGA), isMappingIGA2FEM(_isMappingIGA2FEM) {

    assert(_meshIGA != NULL);
    assert(_meshFE != NULL);
    assert(_meshIGA->type == EMPIRE_Mesh_IGAMesh);
    assert(_meshFE->type == EMPIRE_Mesh_FEMesh);

    if (_meshFE->triangulate() == NULL)
        meshFE = _meshFE;
    else
        meshFE = _meshFE->triangulate();

    projectedCoords.resize(meshFE->numNodes);
    projectedCPs.resize(meshIGA->getNumNodes());
    projectedCPsOnFEMesh.resize(meshIGA->getNumNodes());
    isProjectionOfCpOnTrimmed.resize(meshIGA->getNumNodes(), false);
    meshIGACpNet.resize(meshIGA->getNumNodes());

    if (isMappingIGA2FEM) {
        numNodesSlave = meshIGA->getNumNodes();
        numNodesMaster = meshFE->numNodes;
        C_M = new MathLibrary::SparseMatrix<double>(numNodesMaster, numNodesSlave);
    } else {
        numNodesSlave = meshFE->numNodes;
        numNodesMaster = meshIGA->getNumNodes(); // number of CPs
        C_L = new MathLibrary::SparseMatrix<double>(numNodesMaster, numNodesMaster);
        C_LT = new MathLibrary::SparseMatrix<double>(numNodesMaster, numNodesMaster);
        C_R = new MathLibrary::SparseMatrix<double>(numNodesMaster, numNodesSlave);
    }

    setParametersProjection();
    setParametersNewtonRaphson();
}

IGABarycentricMapper::~IGABarycentricMapper() {

    for (int i = 0; i < meshFE->numElems; i++)
        delete[] meshFEDirectElemTable[i];
    delete[] meshFEDirectElemTable;

    if (isMappingIGA2FEM) {
        delete C_M;
    } else {
        delete C_L;
        delete C_LT;
        delete C_R;
        delete numNodesPerNeighborElem;
    }
}

void IGABarycentricMapper::setParametersNewtonRaphson(int _maxNumOfIterations, double _tolerance) {
    newtonRaphson.maxNumOfIterations=_maxNumOfIterations;
    newtonRaphson.tolerance=_tolerance;
}

void IGABarycentricMapper::setParametersProjection(double _maxProjectionDistance, int _numRefinementForIntialGuess, double _maxDistanceForProjectedPointsOnDifferentPatches) {
    projectionProperties.maxProjectionDistance=_maxProjectionDistance;
    projectionProperties.numRefinementForIntialGuess=_numRefinementForIntialGuess;
    projectionProperties.maxDistanceForProjectedPointsOnDifferentPatches=_maxDistanceForProjectedPointsOnDifferentPatches;
}

void IGABarycentricMapper::buildCouplingMatrices() {
	int nIG = meshIGA->getNumNodes();
	int nFE = meshFE->numNodes;
    
    HEADING_OUT(3, "IGABarycentricMapper", "Building coupling matrices for ("+ name +")...", infoOut);
    {
        INFO_OUT() << "Number of nodes in NURBS mesh is " << nIG << endl;
        INFO_OUT() << "Number of nodes in FE mesh is    " << nFE << endl;
        if (isMappingIGA2FEM)
            INFO_OUT() << "Size of coupling matrix will be " << nFE << "x" << nIG << endl;
        else
            INFO_OUT() << "Size of coupling matrices will be " << nIG << "x" << nIG << " and " << nIG << "x" << nFE << endl;
    }

    //Set default scheme values
    IGAPatchSurface::MAX_NUM_ITERATIONS = newtonRaphson.maxNumOfIterations;
    IGAPatchSurface::TOL_ORTHOGONALITY = newtonRaphson.tolerance;

    initTables();

    if (isMappingIGA2FEM) {
        projectPointsToSurface();
        // Write the projected points on to a file
        writeProjectedNodesOntoIGAMesh();

    } else {
        projectCPsToSurface();
        computeNeighborsAndWeights();
        // Write the projected control points on to a file
        writeProjectedCPsOntoFEMMesh();
    }

    computeCouplingMatrices();

    writeCouplingMatricesToFile();

    checkConsistency();
}

void IGABarycentricMapper::initTables() {
    /* using the map to store the nodeIDs
     * but here the "key" is the node ID, and the value is the position in nodeIDs
     * the map is sorted automatically, so it is efficient for searching
     */

    // compute direct element table for fluid mesh
    meshFEDirectElemTable = new int*[meshFE->numElems]; // deleted
    for (int i = 0; i < meshFE->numElems; i++)
        meshFEDirectElemTable[i] = new int[meshFE->numNodesPerElem[i]];

    map<int, int> meshFENodesMap;
    for (int i = 0; i < meshFE->numNodes; i++)
        meshFENodesMap.insert(meshFENodesMap.end(), pair<int, int>(meshFE->nodeIDs[i], i));
    int count = 0;

    for (int i = 0; i < meshFE->numElems; i++) {
        const int numNodesPerElem = meshFE->numNodesPerElem[i];

        for (int j = 0; j < numNodesPerElem; j++) {
            if (meshFENodesMap.find(meshFE->elems[count + j]) == meshFENodesMap.end()) {
                ERROR_OUT() << "Cannot find node ID " << meshFE->elems[count + j] << endl;
                exit(-1);
            }
            meshFEDirectElemTable[i][j] = meshFENodesMap.at(meshFE->elems[count + j]);
        }
        count += numNodesPerElem;
    }

    for (int node = 0; node < meshFE->numNodes; node++) {
        for (int elem = 0; elem < meshFE->numElems; elem++) {
            const int numNodesPerElem = meshFE->numNodesPerElem[elem];
            int* out = find(meshFEDirectElemTable[elem],meshFEDirectElemTable[elem]+numNodesPerElem,node);
            if(out != meshFEDirectElemTable[elem]+numNodesPerElem) {
            	meshFENodeToElementTable[node].push_back(elem);
            }
        }
    }
}

void IGABarycentricMapper::projectCPsToSurface() {
    // Time stamps
    time_t timeStart, timeEnd;

    // Initialization of variables

    // Array of booleans containing flags on the projection of the FE nodes onto the NURBS patch
    // A CP needs to be projected at least once
    vector<bool> isProjected(meshIGA->getNumNodes(), false);
    // Keep track of the minimum distance found between a CP and a patch
    vector<double> minProjectionDistance(meshIGA->getNumNodes(), 1e9);
    // Keep track of the point on patch related to minimum distance
    vector<vector<double> > minProjectionPoint(meshIGA->getNumNodes());
    // List of patch to try a projection for every control point
    vector<set<int> > patchToProcessPerCP(meshIGA->getNumNodes());

    // Initial guess for projection onto the NURBS patch
    double initialU, initialV;

    // Get the number of patches in the IGA mesh
    int numPatches = meshIGA->getNumPatches();

    // Create vector of the control point net
    std::vector<IGAPatchSurface*> patches = meshIGA->getSurfacePatches();
    for (int i = 0; i < numPatches; ++i) {
        IGAControlPoint** cpNet = patches[i]->getControlPointNet();
        int num = patches[i]->getNoControlPoints();
        for (int j = 0; j < num; ++j) {
            meshIGACpNet[ cpNet[j]->getDofIndex() ] = cpNet[j];
        }
    }

    // Bounding box preprocessing, assign to each node the patches to be visited
    INFO_OUT()<<"Bounding box preprocessing..."<<endl;
    time(&timeStart);

    for (int i = 0; i < meshIGA->getNumNodes(); i++) {
        double P[3];
        P[0] = meshIGACpNet[i]->getX();
        P[1] = meshIGACpNet[i]->getY();
        P[2] = meshIGACpNet[i]->getZ();
        for(int patchCount = 0; patchCount < numPatches; patchCount++) {
            IGAPatchSurface* thePatch = meshIGA->getSurfacePatch(patchCount);
            bool isInside = thePatch->getBoundingBox().isPointInside(P, projectionProperties.maxProjectionDistance);
            if(isInside)
                patchToProcessPerCP[i].insert(patchCount);
        }

        if(patchToProcessPerCP[i].empty()) {
            stringstream msg;
            msg << "Control point [" << i << "] is not in any bounding box of NURBS patches ! Increase maxProjectionDistance !";
            ERROR_BLOCK_OUT("IGABarycentricMapper", "projectCPsToSurface", msg.str());
        }
    }

    time(&timeEnd);
    INFO_OUT()<<"Bounding box preprocessing done in "<< difftime(timeEnd, timeStart) << " seconds."<<endl;

    // Project the control point of the patches
    INFO_OUT()<<"Control point first pass projection..."<<endl;

    time(&timeStart);
    for(int patchIndex = 0; patchIndex < numPatches; patchIndex++) {
        IGAPatchSurface* thePatch = meshIGA->getSurfacePatch(patchIndex);
        IGAControlPoint **cpNet = thePatch->getControlPointNet();
        int num = thePatch->getNoControlPoints();
        for(int i = 0; i < num; i++) {
            int counterCP = cpNet[i]->getDofIndex();
            // If support of CP is completely in trimmed region, ignore CP
            double u1, u2, v1, v2;
            knotSpanOfSupportCP(patchIndex, counterCP, u1, u2, v1, v2);
            if (!projectionInside(patchIndex, u1, v1) && !projectionInside(patchIndex, u2, v1) && !projectionInside(patchIndex, u1, v2) && !projectionInside(patchIndex, u2, v2) ) {
                isProjectionOfCpOnTrimmed[counterCP] = true;
                continue;
            }
            // If already projected, go to next CP
            if(projectedCPs[counterCP].find(patchIndex) != projectedCPs[counterCP].end())
                continue;
            // If node in BBox of patch
            if(patchToProcessPerCP[counterCP].find(patchIndex) != patchToProcessPerCP[counterCP].end()) {
                computeInitialGuessForProjectionOfCPs(patchIndex, counterCP, initialU, initialV);
                bool flagProjected = projectCpOnPatch(patchIndex, counterCP, initialU, initialV, minProjectionDistance[counterCP], minProjectionPoint[counterCP]);
                isProjected[counterCP] = isProjected[counterCP] || flagProjected;
            }
        }
    }
    time(&timeEnd);

    INFO_OUT() << "Control point first pass projection done in " << difftime(timeEnd, timeStart) << " seconds."<<endl;

    int missing = 0;
    numOfValidCPs = 0;
    for (int i = 0; i < meshIGA->getNumNodes(); i++) {
        if(!isProjectionOfCpOnTrimmed[i]) {
            numOfValidCPs++;
            if(!isProjected[i]) {
                missing++;
                WARNING_OUT() << "Control point not projected: [" << i << "] of coordinates " << meshIGACpNet[i]->getX() << "," << meshIGACpNet[i]->getY() << "," << meshIGACpNet[i]->getZ() << endl;
            }
        }
    }

    INFO_OUT()<<meshIGA->getNumNodes() - missing << " nodes over " << meshIGA->getNumNodes() <<" could be projected during first pass." << endl;

    double initialTolerance = newtonRaphson.tolerance;
    int initialNumRefs = projectionProperties.numRefinementForIntialGuess;
    if(missing) {
        INFO_OUT()<<"Control point second pass projection..."<<endl;
        time(&timeStart);
        missing = 0;
        numOfValidCPs = 0;
        for(int patchIndex = 0; patchIndex < numPatches; patchIndex++) {
            IGAPatchSurface* thePatch = meshIGA->getSurfacePatch(patchIndex);
            IGAControlPoint **cpNet = thePatch->getControlPointNet();
            int num = thePatch->getNoControlPoints();
            for(int i = 0; i < num; i++) {
                int counterCP = cpNet[i]->getDofIndex();
                if(!isProjectionOfCpOnTrimmed[counterCP]) {
                    numOfValidCPs++;
                    if(!isProjected[counterCP]) {
                        newtonRaphson.tolerance = 10*newtonRaphson.tolerance;
                        projectionProperties.numRefinementForIntialGuess = 2*projectionProperties.numRefinementForIntialGuess;
                        computeInitialGuessForProjectionOfCPs(patchIndex, counterCP, initialU, initialV);
                        bool flagProjected = projectCpOnPatch(patchIndex, counterCP, initialU, initialV, minProjectionDistance[counterCP], minProjectionPoint[counterCP]);
                        isProjected[counterCP] = isProjected[counterCP] || flagProjected;
                    }
                    if(!isProjected[counterCP]) {
                        WARNING_OUT()<<"Control point not projected at second pass: ["<<counterCP<<"] of coordinates "<<meshIGACpNet[counterCP]->getX() << "," << meshIGACpNet[counterCP]->getY() << "," << meshIGACpNet[counterCP]->getZ() << endl;
                        
                        // // Work-around for flying nodes: check if knot values at boundaries of area of effect of CP are in trimmed region.
                        double u1, u2, v1, v2;

                        int numCPsInPatch = thePatch->getNoControlPoints();
                        int pDegree = thePatch->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree();
                        int qDegree = thePatch->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree();
                        double * pKnotVector = thePatch->getIGABasis()->getUBSplineBasis1D()->getKnotVector();
                        double * qKnotVector = thePatch->getIGABasis()->getVBSplineBasis1D()->getKnotVector();
                        int pNumKnots = thePatch->getIGABasis()->getUBSplineBasis1D()->getNoKnots();
                        int qNumKnots = thePatch->getIGABasis()->getVBSplineBasis1D()->getNoKnots();
                        int pNumCP = pNumKnots - pDegree - 1;
                        int qNumCP = qNumKnots - qDegree - 1;

                        int ip, iq;
                        ip = i % pNumCP;
                        iq = i / qNumCP;

                        u1 = pKnotVector[ip + pDegree/2];
                        u2 = pKnotVector[ip + pDegree/2 + 1];
                        v1 = qKnotVector[iq + qDegree/2];
                        v2 = qKnotVector[iq + qDegree/2 + 1];

                        if ( !projectionInside(patchIndex, u1, v1) && !projectionInside(patchIndex, u2, v1) && !projectionInside(patchIndex, u1, v2) && !projectionInside(patchIndex, u2, v2) ) {
                            isProjectionOfCpOnTrimmed[counterCP] = true;
                            WARNING_OUT()<<"...but it's a flying node and it will be ignored.\n";
                        }
                        else
                            missing++;
                    }
                }
                newtonRaphson.tolerance = initialTolerance;
                projectionProperties.numRefinementForIntialGuess = initialNumRefs;
            }
        }
        newtonRaphson.tolerance = initialTolerance;
        projectionProperties.numRefinementForIntialGuess = initialNumRefs;
        time(&timeEnd);
        INFO_OUT()<<"Control point second pass projection done! It took "<< difftime(timeEnd, timeStart) << " seconds."<<endl;

        if(missing) {
            stringstream msg;
            msg << missing << " control points over " << meshIGA->getNumNodes() << " could NOT be projected during second pass !";
            ERROR_BLOCK_OUT("IGABarycentricMapper", "ProjectCPsToSurface", msg.str());
        }
    }
}

void IGABarycentricMapper::computeInitialGuessForProjectionOfCPs(const int _patchIndex, const int _cpIndex, double& _u, double& _v) {
    IGAPatchSurface* thePatch = meshIGA->getSurfacePatch(_patchIndex);
    // Get the Cartesian coordinates of the input CP
    double P[3];
    P[0] = meshIGACpNet[_cpIndex]->getX();
    P[1] = meshIGACpNet[_cpIndex]->getY();
    P[2] = meshIGACpNet[_cpIndex]->getZ();
    // Get an initial guess for the projection onto the NURBS patch
    thePatch->findInitialGuess4PointProjection(_u, _v, P, projectionProperties.numRefinementForIntialGuess, projectionProperties.numRefinementForIntialGuess);
}

bool IGABarycentricMapper::projectCpOnPatch(const int _patchIndex, const int _cpIndex, const double _u0, const double _v0, double& _minProjectionDistance, std::vector<double>& _minProjectionPoint) {
    IGAPatchSurface* thePatch = meshIGA->getSurfacePatch(_patchIndex);
    /// Get the Cartesian coordinates of the CP
    double P[3], projectedP[3];
    projectedP[0] = P[0] = meshIGACpNet[_cpIndex]->getX();
    projectedP[1] = P[1] = meshIGACpNet[_cpIndex]->getY();
    projectedP[2] = P[2] = meshIGACpNet[_cpIndex]->getZ();
    /// Get an initial guess for the parametric location of the projected CP on the NURBS patch
    double u = _u0;
    double v = _v0;
    /// Compute point projection on the NURBS patch using the Newton-Rapshon iteration method
    bool hasResidualConverged;
    bool hasConverged = thePatch->computePointProjectionOnPatch(u, v, projectedP,
            hasResidualConverged, newtonRaphson.maxNumOfIterations, newtonRaphson.tolerance);
    double distance = MathLibrary::computePointDistance(P, projectedP);

    /// Return false if the projection falls on a trimmed patch region
    if(!projectionInside(_patchIndex, u, v)) {
        isProjectionOfCpOnTrimmed[_cpIndex] = true;
        return false;
    }

    if(hasConverged && distance < projectionProperties.maxProjectionDistance) {
        /// Perform some validity checks to validate the projected point
        if(distance > _minProjectionDistance + projectionProperties.maxDistanceForProjectedPointsOnDifferentPatches) {
            return false;
        }
        if(!_minProjectionPoint.empty() &&
                MathLibrary::computePointDistance(projectedP, &_minProjectionPoint[0]) > projectionProperties.maxDistanceForProjectedPointsOnDifferentPatches &&
                distance > _minProjectionDistance) {
            return false;
        }
        if(distance < _minProjectionDistance - projectionProperties.maxDistanceForProjectedPointsOnDifferentPatches
                || MathLibrary::computePointDistance(projectedP, &_minProjectionPoint[0]) > projectionProperties.maxDistanceForProjectedPointsOnDifferentPatches) {
            projectedCPs[_cpIndex].clear();
        }
        /// Store result
        vector<double> uv(2);
        uv[0] = u;
        uv[1] = v;
        projectedCPs[_cpIndex].insert(make_pair(_patchIndex, uv));
        _minProjectionDistance = distance;
        _minProjectionPoint = vector<double>(projectedP, projectedP + 3);
        return true;
    }
    return false;
}

void IGABarycentricMapper::computeNeighborsAndWeights() {
    // Array that holds the number of nodes of each element
    numNodesPerNeighborElem = new int[meshIGA->getNumNodes()];
    // Compute elementCentroids
    double *elementCentroids = new double[meshFE->numElems * 3];
    for (int i = 0; i < meshFE->numElems; i++) { // compute the centroid of all elements in FE mesh
        int numNodesThisElem = meshFE->numNodesPerElem[i];
        double thisElem[numNodesThisElem * 3];
        getElemCoorInFEMMesh(i, thisElem);
        EMPIRE::MathLibrary::computePolygonCenter(thisElem, numNodesThisElem, &(elementCentroids[i * 3]));
    }

    int NUM_NEIGHBORS_TO_SEARCH = MAX_NUM_NEIGHBORS_TO_SEARCH;
    if (NUM_NEIGHBORS_TO_SEARCH > meshFE->numElems) {
        NUM_NEIGHBORS_TO_SEARCH = meshFE->numElems;
    }

    {
    #ifdef FLANN
        double *castedCPs = new double[meshIGA->getNumNodes() * 3];
        for (int i = 0; i < meshIGA->getNumNodes(); ++i) {
            if (!isProjectionOfCpOnTrimmed[i]) {
                IGAPatchSurface* thePatch = meshIGA->getSurfacePatch(projectedCPs[i].begin()->first);
                double cartesianCoordinatesOfCP[3]; // cartesian coordinates of projection of the control point on the patch
                double uv[2] = {projectedCPs[i].begin()->second[0], projectedCPs[i].begin()->second[1]};
                thePatch->computeCartesianCoordinates(cartesianCoordinatesOfCP, uv);
                castedCPs[3 * i + 0] = cartesianCoordinatesOfCP[0];
                castedCPs[3 * i + 1] = cartesianCoordinatesOfCP[1];
                castedCPs[3 * i + 2] = cartesianCoordinatesOfCP[2];
            }
        }

        /// Build binary tree for searching
        const double *constCPs = const_cast<double*>(castedCPs);
        flann::Matrix<double> *elementCentroids_FLANN = new flann::Matrix<double>(
                const_cast<double*>(elementCentroids), meshFE->numElems, 3);
        flann::Index<flann::L2<double> > *ANodesTree = new flann::Index<flann::L2<double> >(
                *elementCentroids_FLANN, flann::KDTreeSingleIndexParams(1));
        ANodesTree->buildIndex();

        /// Traverse through all valid CPs and project them to the FE mesh; store element index and local element coordinates
        for (int i = 0; i < meshIGA->getNumNodes(); i++) {
            if (!isProjectionOfCpOnTrimmed[i]) {
                flann::Matrix<double> nodeI(&(castedCPs[i * 3]), 1, 3);
                vector<vector<int> > indexes_tmp;
                vector<vector<double> > dists_tmp;
                ANodesTree->knnSearch(nodeI, indexes_tmp, dists_tmp, NUM_NEIGHBORS_TO_SEARCH,
                        flann::SearchParams(1));

                for (int j = 0; j < NUM_NEIGHBORS_TO_SEARCH; j++) {
                    int numNodesThisElem = meshFE->numNodesPerElem[indexes_tmp[0][j]];
                    if (numNodesThisElem == 3) {
                        double triangle[3 * 3];
                        getElemCoorInFEMMesh(indexes_tmp[0][j], triangle);
                        double normal[3];
                        EMPIRE::MathLibrary::computeNormalOfTriangle(triangle, true, normal);
                        int planeToProject = EMPIRE::MathLibrary::computePlaneToProject(normal);

                        double projection[3];
                        EMPIRE::MathLibrary::projectToPlane(&triangle[0], normal, &(constCPs[i * 3]), 1,
                                projection);

                        double localCoors[3];
                        EMPIRE::MathLibrary::computeLocalCoorInTriangle(triangle, planeToProject, projection,
                                localCoors);
                        bool inside = insideElement(3, localCoors);
                        if (inside) {
                            numNodesPerNeighborElem[i] = 3;
                            int elemIndex = indexes_tmp[0][j];
                            vector<double> localCoords(3);
                            for (int k = 0; k < 3; k++) {
                                localCoords[k] = localCoors[k];
                            }
                            projectedCPsOnFEMesh[i] = make_pair(elemIndex, localCoords);
                            break;
                        }
                    } else if (numNodesThisElem == 4) {
                        double quad[4 * 3];
                        getElemCoorInFEMMesh(indexes_tmp[0][j], quad);
                        double normal[3];
                        EMPIRE::MathLibrary::computeNormalOfQuad(quad, true, normal);
                        int planeToProject = EMPIRE::MathLibrary::computePlaneToProject(normal);

                        { // replace the element by the projection of it on its "element plane"
                            double quadCenter[3];
                            EMPIRE::MathLibrary::computePolygonCenter(quad, 4, quadCenter);
                            double quadPrj[12];
                            EMPIRE::MathLibrary::projectToPlane(quadCenter, normal, quad, 4, quadPrj);
                            for (int i = 0; i < 12; i++)
                            quad[i] = quadPrj[i];
                        }
                        double projection[3];
                        EMPIRE::MathLibrary::projectToPlane(&quad[0], normal, &(constCPs[i * 3]), 1,
                                projection);

                        double localCoors[2];
                        EMPIRE::MathLibrary::computeLocalCoorInQuad(quad, planeToProject, projection,
                                localCoors);
                        bool inside = insideElement(4, localCoors);
                        if (inside) {
                            numNodesPerNeighborElem[i] = 4;
                            int elemIndex = indexes_tmp[0][j];
                            vector<double> localCoords(2);
                            for (int k = 0; k < 2; k++) {
                                localCoords[k] = localCoors[k];
                            }
                            projectedCPsOnFEMesh[i] = make_pair(elemIndex, localCoords);
                            break;
                        }
                    } else {
                        assert(false);
                    }
                    if (j == NUM_NEIGHBORS_TO_SEARCH - 1) { // projections do not locate inside the neighboring elements, use the nearest one
                        int numNodesThisElem = meshFE->numNodesPerElem[indexes_tmp[0][0]];
                        if (numNodesThisElem == 3) {
                            double triangle[3 * 3];
                            getElemCoorInFEMMesh(indexes_tmp[0][0], triangle);
                            double normal[3];
                            EMPIRE::MathLibrary::computeNormalOfTriangle(triangle, true, normal);
                            int planeToProject = EMPIRE::MathLibrary::computePlaneToProject(normal);

                            double projection[3];
                            EMPIRE::MathLibrary::projectToPlane(&triangle[0], normal, &(constCPs[i * 3]), 1,
                                    projection);

                            double localCoors[3];
                            EMPIRE::MathLibrary::computeLocalCoorInTriangle(triangle, planeToProject, projection,
                                    localCoors);
                            {
                                numNodesPerNeighborElem[i] = 3;
                                int elemIndex = indexes_tmp[0][0];
                                vector<double> localCoords(3);
                                for (int k = 0; k < 3; k++) {
                                    localCoords[k] = localCoors[k];
                                }
                                projectedCPsOnFEMesh[i] = make_pair(elemIndex, localCoords);
                            }
                        } else if (numNodesThisElem == 4) {
                            double quad[4 * 3];
                            getElemCoorInFEMMesh(indexes_tmp[0][0], quad);
                            double normal[3];
                            EMPIRE::MathLibrary::computeNormalOfQuad(quad, true, normal);
                            int planeToProject = EMPIRE::MathLibrary::computePlaneToProject(normal);

                            { // replace the element by the projection of it on its "element plane"
                                double quadCenter[3];
                                EMPIRE::MathLibrary::computePolygonCenter(quad, 4, quadCenter);
                                double quadPrj[12];
                                EMPIRE::MathLibrary::projectToPlane(quadCenter, normal, quad, 4, quadPrj);
                                for (int i = 0; i < 12; i++)
                                quad[i] = quadPrj[i];
                            }
                            double projection[3];
                            EMPIRE::MathLibrary::projectToPlane(&quad[0], normal, &(constCPs[i * 3]), 1,
                                    projection);

                            double localCoors[2];
                            EMPIRE::MathLibrary::computeLocalCoorInQuad(quad, planeToProject, projection,
                                    localCoors);
                            bool inside = insideElement(4, localCoors);
                            {
                                numNodesPerNeighborElem[i] = 4;
                                int elemIndex = indexes_tmp[0][0];
                                vector<double> localCoords(2);
                                for (int k = 0; k < 2; k++) {
                                    localCoords[k] = localCoors[k];
                                }
                                projectedCPsOnFEMesh[i] = make_pair(elemIndex, localCoords);
                            }
                        } else {
                            assert(false);
                        }
                    }
                }
            }
        }
        delete elementCentroids_FLANN;
        delete ANodesTree;
        delete castedCPs;
    #endif
    }
    delete[] elementCentroids;
}

void IGABarycentricMapper::knotSpanOfSupportCP(const int _patchIndex, const int _cpIndex, double& _u1, double& _u2, double& _v1, double& _v2) {
    IGAPatchSurface* thePatch = meshIGA->getSurfacePatch(_patchIndex);
    IGAControlPoint **cpNet = thePatch->getControlPointNet();
    int numCPsInPatch = thePatch->getNoControlPoints();
    int pDegree = thePatch->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree();
    int qDegree = thePatch->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree();
    double * pKnotVector = thePatch->getIGABasis()->getUBSplineBasis1D()->getKnotVector();
    double * qKnotVector = thePatch->getIGABasis()->getVBSplineBasis1D()->getKnotVector();
    int pNumKnots = thePatch->getIGABasis()->getUBSplineBasis1D()->getNoKnots();
    int qNumKnots = thePatch->getIGABasis()->getVBSplineBasis1D()->getNoKnots();
    int pNumCP = pNumKnots - pDegree - 1;
    int qNumCP = qNumKnots - qDegree - 1;

    int i;
    for(i = 0; i < numCPsInPatch; i++) {
        if (cpNet[i]->getDofIndex() == _cpIndex)
            break;
    }

    if (i == numCPsInPatch) {
        ERROR_OUT() << "Trying to find support for foreign control point [" << _cpIndex << "] to patch [" << _patchIndex << "] !\n";
        return;
    }

    int ip, iq;
    ip = i % pNumCP;
    iq = i / qNumCP;

    _u1 = pKnotVector[ip];
    _u2 = pKnotVector[ip + pDegree + 1];
    _v1 = qKnotVector[iq];
    _v2 = qKnotVector[iq + qDegree + 1];
}

void IGABarycentricMapper::getElemCoorInFEMMesh(int elemIndex, double *elem) {
    // compute the coordinates of an element by its id
    int numNodesThisElem = meshFE->numNodesPerElem[elemIndex];
    for (int i = 0; i < numNodesThisElem; i++) {
        int nodePos = meshFEDirectElemTable[elemIndex][i]; // position of the node
        for (int j = 0; j < 3; j++) {
            elem[i * 3 + j] = meshFE->nodes[nodePos * 3 + j];
        }
    }
}

bool IGABarycentricMapper::insideElement(int numNodesThisElem, double *localCoor) {
    const double EPS = 1e-10;
    if (numNodesThisElem == 3) {
        for (int i = 0; i < 3; i++) {
            if (localCoor[i] > 1.0 + EPS)
                return false;
            if (localCoor[i] < 0.0 - EPS)
                return false;
        }
        return true;
    } else if (numNodesThisElem == 4) {
        for (int i = 0; i < 2; i++) {
            if (localCoor[i] > 1.0 + EPS)
                return false;
            if (localCoor[i] < -1.0 - EPS)
                return false;
        }
        return true;
    } else {
        assert(false);
        return false;
    }
}

void IGABarycentricMapper::projectPointsToSurface() {
	// Time stamps
	time_t timeStart, timeEnd;

    // Initialization of variables

    // Array of booleans containing flags on the projection of the FE nodes onto the NURBS patch
	// A node needs to be projected at least once
	vector<bool> isProjected(meshFE->numNodes);
	// Keep track of the minimum distance found between a node and a patch
    vector<double> minProjectionDistance(meshFE->numNodes, 1e9);
    // Keep track of the point on patch related to minimum distance
    vector<vector<double> > minProjectionPoint(meshFE->numNodes);
    // List of patch to try a projection for every node
    vector<set<int> > patchToProcessPerNode(meshFE->numNodes);

    // Initial guess for projection onto the NURBS patch
    double initialU, initialV;

    // Initialize the array of the Cartesian coordinates of a node in the FE side
    double P[3];

    // Get the number of patches in the IGA mesh
    int numPatches = meshIGA->getNumPatches();
    // Bounding box preprocessing, assign to each node the patches to be visited
    INFO_OUT()<<"Bounding box preprocessing..."<<endl;
    time(&timeStart);

    for (int i = 0; i < meshFE->numNodes; i++) {
    	double P[3];
    	P[0] = meshFE->nodes[3 * i + 0];
    	P[1] = meshFE->nodes[3 * i + 1];
    	P[2] = meshFE->nodes[3 * i + 2];
        for(int patchCount = 0; patchCount < numPatches; patchCount++) {
            IGAPatchSurface* thePatch = meshIGA->getSurfacePatch(patchCount);
            bool isInside = thePatch->getBoundingBox().isPointInside(P, projectionProperties.maxProjectionDistance);
            if(isInside)
            	patchToProcessPerNode[i].insert(patchCount);
        }

        if(patchToProcessPerNode[i].empty()) {
            stringstream msg;
            msg << "Node [" << i << "] is not in any bounding box of NURBS patches ! Increase maxProjectionDistance !";
            ERROR_BLOCK_OUT("IGABarycentricMapper", "projectPointsToSurface", msg.str());
        }
    }

    time(&timeEnd);
    INFO_OUT()<<"Bounding box preprocessing done in "<< difftime(timeEnd, timeStart) << " seconds."<<endl;
    // Project the node for every patch's bounding box the node lies into
    // or on every patch if not found in a single bounding box
    INFO_OUT()<<"First pass projection..."<<endl;
    time(&timeStart);
    for(int i = 0; i < meshFE->numElems; i++) {
        int numNodesInElem = meshFE->numNodesPerElem[i];
        for(int patchIndex = 0; patchIndex < numPatches; patchIndex++) {


            IGAPatchSurface* thePatch = meshIGA->getSurfacePatch(patchIndex);
            bool initialGuessComputed = false;
            for(int j = 0; j < numNodesInElem; j++) {
                int nodeIndex = meshFEDirectElemTable[i][j];
                // If already projected, go to next node
                if(projectedCoords[nodeIndex].find(patchIndex) != projectedCoords[nodeIndex].end())
                    continue;
                // If node in BBox of patch
                if(patchToProcessPerNode[nodeIndex].find(patchIndex) != patchToProcessPerNode[nodeIndex].end()) {
                    if(!initialGuessComputed) {
                        computeInitialGuessForProjection(patchIndex, i, nodeIndex, initialU, initialV);
                        initialGuessComputed = true;
                    }
                    bool flagProjected = projectPointOnPatch(patchIndex, nodeIndex, initialU, initialV, minProjectionDistance[nodeIndex], minProjectionPoint[nodeIndex]);
                    isProjected[nodeIndex] = isProjected[nodeIndex] || flagProjected;
                }
            }
        }
    }
    time(&timeEnd);
    INFO_OUT()<<"First pass projection done in "<< difftime(timeEnd, timeStart) << " seconds."<<endl;
    int missing = 0;
    for (int i = 0; i < meshFE->numNodes; i++) {
        if(!isProjected[i]) {
            missing++;
            WARNING_OUT()<<"Node not projected at first pass: ["<<i<<"] of coordinates "<<meshFE->nodes[3*i]<<","<<meshFE->nodes[3*i+1]<<","<<meshFE->nodes[3*i+2]<<endl;
        }
    }
    INFO_OUT()<<meshFE->numNodes - missing << " nodes over " << meshFE->numNodes <<" could be projected during first pass." << endl;

    double initialTolerance = newtonRaphson.tolerance;
    if(missing) {
        INFO_OUT()<<"Second pass projection..."<<endl;
        time(&timeStart);
        missing = 0;
        for (int i = 0; i < meshFE->numNodes; i++) {
            if(!isProjected[i]) {
                newtonRaphson.tolerance = 10*newtonRaphson.tolerance;
                for(set<int>::iterator patchIndex=patchToProcessPerNode[i].begin();patchIndex!=patchToProcessPerNode[i].end();patchIndex++) {
                    IGAPatchSurface* thePatch = meshIGA->getSurfacePatch(*patchIndex);
                    computeInitialGuessForProjection(*patchIndex, meshFENodeToElementTable[i][0], i, initialU, initialV);
                    bool flagProjected = projectPointOnPatch(*patchIndex, i, initialU, initialV, minProjectionDistance[i], minProjectionPoint[i]);
                    isProjected[i] = isProjected[i] || flagProjected;
                }
                if(!isProjected[i]) {
                    for(set<int>::iterator patchIndex=patchToProcessPerNode[i].begin();patchIndex!=patchToProcessPerNode[i].end();patchIndex++) {
                        bool flagProjected = forceProjectPointOnPatch(*patchIndex, i, minProjectionDistance[i], minProjectionPoint[i]);
                        isProjected[i] = isProjected[i] || flagProjected;
                    }
                }
            }
            if(!isProjected[i]) {
                ERROR_OUT()<<"Node not projected at second pass: ["<<i<<"] of coordinates "<<meshFE->nodes[3*i]<<","<<meshFE->nodes[3*i+1]<<","<<meshFE->nodes[3*i+2]<<endl;
                missing++;
            }
            newtonRaphson.tolerance = initialTolerance;
        }
        newtonRaphson.tolerance = initialTolerance;
        time(&timeEnd);
        INFO_OUT()<<"Second pass projection done! It took "<< difftime(timeEnd, timeStart) << " seconds."<<endl;

        if(missing) {
            stringstream msg;
            msg << missing << " nodes over " << meshFE->numNodes << " could NOT be projected during second pass !";
            ERROR_BLOCK_OUT("IGABarycentricMapper", "ProjectPointsToSurface", msg.str());
        }
    }
}

void IGABarycentricMapper::computeInitialGuessForProjection(const int _patchIndex, const int _elemIndex, const int _nodeIndex, double& _u, double& _v) {
    IGAPatchSurface* thePatch = meshIGA->getSurfacePatch(_patchIndex);
	/// 1. Initialize the flag to false and the node id to zero
    bool isNodeInsideElementProjected = false;
    int projectedNode = -1;
	/// 2. Loop over all nodes of the current element to check if there exist one node has been projected already
    for (int j = 0; j < meshFE->numNodesPerElem[_elemIndex]; j++) {
    	int nodeIndex = meshFEDirectElemTable[_elemIndex][j];
        if (projectedCoords[nodeIndex].find(_patchIndex) != projectedCoords[nodeIndex].end()) {
            /// 1iii.2i. If the node has already been projected set projection flag to true
            isNodeInsideElementProjected = true;
            /// 1iii.2ii. Get the global ID of the projected node
            projectedNode = nodeIndex;
            /// 1iii.2iii. Break the loop
            break;
        }
    }
    /// 3. Check if there exist one node in the current element has been successfully projected
    if (isNodeInsideElementProjected) {
        /// 3i. If so, use result of the projected node as the initial guess for the projection step
        _u = projectedCoords[projectedNode][_patchIndex][0];
        _v = projectedCoords[projectedNode][_patchIndex][1];
    } else {
        /// 3ii. Otherwise, find the nearest knot intersection as initial guess for the projection step
        // Get the Cartesian coordinates of that input node
        double P[3];
        P[0] = meshFE->nodes[_nodeIndex * 3 + 0];
        P[1] = meshFE->nodes[_nodeIndex * 3 + 1];
        P[2] = meshFE->nodes[_nodeIndex * 3 + 2];
        // Get accordingly an initial guess for the projection onto the NURBS patch
        thePatch->findInitialGuess4PointProjection(_u, _v, P, projectionProperties.numRefinementForIntialGuess, projectionProperties.numRefinementForIntialGuess);
    }
}

bool IGABarycentricMapper::projectPointOnPatch(const int patchIndex, const int nodeIndex, const double u0, const double v0, double& minProjectionDistance, vector<double>& minProjectionPoint) {
    IGAPatchSurface* thePatch = meshIGA->getSurfacePatch(patchIndex);
    /// Get the Cartesian coordinates of the node in the FE side
	double P[3], projectedP[3];
    projectedP[0] = P[0] = meshFE->nodes[nodeIndex * 3 + 0];
    projectedP[1] = P[1] = meshFE->nodes[nodeIndex * 3 + 1];
    projectedP[2] = P[2] = meshFE->nodes[nodeIndex * 3 + 2];
    /// Get an initial guess for the parametric location of the projected node of the FE side on the NURBS patch
    double u = u0;
    double v = v0;
    /// Compute point projection on the NURBS patch using the Newton-Rapshon iteration method
    bool hasResidualConverged;
    bool hasConverged = thePatch->computePointProjectionOnPatch(u, v, projectedP,
            hasResidualConverged, newtonRaphson.maxNumOfIterations, newtonRaphson.tolerance);
    double distance = MathLibrary::computePointDistance(P, projectedP);
    if(!projectionInside(patchIndex, u, v)) {
        return false;
    }
    if(hasConverged &&  distance < projectionProperties.maxProjectionDistance) {
        /// Perform some validity checks to validate the projected point
        if(distance > minProjectionDistance + projectionProperties.maxDistanceForProjectedPointsOnDifferentPatches) {
            return false;
        }
        if(!minProjectionPoint.empty() &&
                MathLibrary::computePointDistance(projectedP, &minProjectionPoint[0]) > projectionProperties.maxDistanceForProjectedPointsOnDifferentPatches &&
                distance > minProjectionDistance) {
            return false;
        }
        if(distance < minProjectionDistance - projectionProperties.maxDistanceForProjectedPointsOnDifferentPatches
                || MathLibrary::computePointDistance(projectedP, &minProjectionPoint[0]) > projectionProperties.maxDistanceForProjectedPointsOnDifferentPatches) {
            projectedCoords[nodeIndex].clear();
        }
        /// Store result
        vector<double> uv(2);
        uv[0] = u;
        uv[1] = v;
        projectedCoords[nodeIndex].insert(make_pair(patchIndex, uv));
        minProjectionDistance = distance;
		minProjectionPoint = vector<double>(projectedP, projectedP + 3);
        return true;
    }
    return false;
}

bool IGABarycentricMapper::forceProjectPointOnPatch(const int patchIndex, const int nodeIndex, double& minProjectionDistance, vector<double>& minProjectionPoint) {
    IGAPatchSurface* thePatch = meshIGA->getSurfacePatch(patchIndex);
    /// Get the Cartesian coordinates of the node in the FE side
	double P[3], projectedP[3];
    projectedP[0] = P[0] = meshFE->nodes[nodeIndex * 3 + 0];
    projectedP[1] = P[1] = meshFE->nodes[nodeIndex * 3 + 1];
    projectedP[2] = P[2] = meshFE->nodes[nodeIndex * 3 + 2];
    /// Get an initial guess for the parametric location of the projected node of the FE side on the NURBS patch
    double u = 0;
    double v = 0;
    for(vector<int>::const_iterator it=meshFENodeToElementTable[nodeIndex].begin();it!=meshFENodeToElementTable[nodeIndex].end();it++) {
    	const int numNodesPerElem = meshFE->numNodesPerElem[*it];
    	for(int i = 0; i< numNodesPerElem; i++) {
    		if(projectedCoords[meshFEDirectElemTable[*it][i]].find(patchIndex) != projectedCoords[meshFEDirectElemTable[*it][i]].end()) {
    			u=projectedCoords[meshFEDirectElemTable[*it][i]][patchIndex][0];
    			v=projectedCoords[meshFEDirectElemTable[*it][i]][patchIndex][1];
    		}
    	}
    }
    /// Compute approximate of parametric position based on brute sampling
    thePatch->findInitialGuess4PointProjection(u,v,P,200,200);
    double uv[2] = {u, v};
    thePatch->computeCartesianCoordinates(projectedP, uv);
    double distance = MathLibrary::computePointDistance(P, projectedP);
    /// Perform some validity checks
    if(!projectionInside(patchIndex, u, v)) {
        return false;
    }
	if(distance > minProjectionDistance + projectionProperties.maxDistanceForProjectedPointsOnDifferentPatches) {
		return false;
	}
	if(distance < minProjectionDistance - projectionProperties.maxDistanceForProjectedPointsOnDifferentPatches) {
		projectedCoords[nodeIndex].clear();
	}
	/// Store result
    vector<double> coordTmp(2);
    coordTmp[0] = u;
    coordTmp[1] = v;
    projectedCoords[nodeIndex].insert(make_pair(patchIndex, coordTmp));
	minProjectionDistance = distance;
	minProjectionPoint = vector<double>(projectedP, projectedP + 3);
    return true;
}

bool IGABarycentricMapper::projectionInside(const int _patchIndex, const double _u, const double _v) {
    /// Examine each projection, whether it is in trimmed region or not
    Point2D projectedPoint;
    projectedPoint = make_pair(_u, _v);
    IGAPatchSurface* thePatch = meshIGA->getSurfacePatch(_patchIndex);

    /// Save trimmed patch to list polygon
    ListPolygon2D listTrimmedPolygon;
    clipPatch(thePatch, listTrimmedPolygon);
    for (int i = 0; i < listTrimmedPolygon.size(); ++i)
        ClipperAdapter::cleanPolygon(listTrimmedPolygon[i]);
    
    bool isPointIn = false; // returns true if point is inside the trimmed patch
    bool isPointOn = false; // returns true if point is on the boundary of the trimmed patch

    /// Check if point is valid
    for (int i = 0; i < listTrimmedPolygon.size(); ++i) {
        int pointIn = ClipperAdapter::isPointIn(projectedPoint, listTrimmedPolygon[i]);
        bool inside = pointIn == 1;
        bool onBoundary = pointIn == -1;
        isPointIn = isPointIn || inside;
        isPointOn = isPointOn || onBoundary;
    }

    // Return false only if: point is not inside and point is not on boundary of the trimmed patch
    if (isPointIn == false && isPointOn == false)
        return false;
    return true;
}

void IGABarycentricMapper::clipPatch(const IGAPatchSurface* _thePatch, ListPolygon2D& _listPolygon) {
    /// Clip patch by its trimmed regions
    const double u0 = _thePatch->getIGABasis()->getUBSplineBasis1D()->getFirstKnot();
    const double v0 = _thePatch->getIGABasis()->getVBSplineBasis1D()->getFirstKnot();
    const double u1 = _thePatch->getIGABasis()->getUBSplineBasis1D()->getLastKnot();
    const double v1 = _thePatch->getIGABasis()->getVBSplineBasis1D()->getLastKnot();
    
    Polygon2D knotSpanWindow(4);
    knotSpanWindow[0]=make_pair(u0,v0);
    knotSpanWindow[1]=make_pair(u1,v0);
    knotSpanWindow[2]=make_pair(u1,v1);
    knotSpanWindow[3]=make_pair(u0,v1);

    ClipperAdapter c;

    if(_thePatch->isTrimmed()) {
        // Fill clipper with trimming clipping curves
        for(int loop=0;loop<_thePatch->getTrimming().getNumOfLoops();loop++) {
            const std::vector<double> clippingWindow=_thePatch->getTrimming().getLoop(loop).getPolylines();
            c.addPathClipper(clippingWindow);
        }
        // Setup filling rule to have for sure clockwise loop as hole and counterclockwise as boundaries
        c.setFilling(ClipperAdapter::POSITIVE, 0);
        c.addPathSubject(knotSpanWindow);
        c.clip();
        c.getSolution(_listPolygon);
    } else {
        _listPolygon.resize(1);
        for (int i = 0; i < 4; ++i) {
            _listPolygon[0].push_back(knotSpanWindow[i]);
        }
    }
}

void IGABarycentricMapper::computeCouplingMatrices() {
    /*
     * Computes the coupling matrices C_M or C_L and C_R.
     * 1. Mapping is from IGA surface to FE mesh
     * ->
     * 1i. Loop over all the elements in the FE side
     * 1ii. Compute the coupling matrix C_M
     * <-
     *
     * 2. Mapping is from FE mesh to IGA surface
     * ->
     * 2i. Loop over all IGA surface control points
     * 2ii. Compute the coupling matrix C_L
     * 2iii. Compute the coupling matrix C_R
     * <-
     */

	// Time stamps
	time_t timeStart, timeEnd;

    time(&timeStart);
    /// 1. Mapping is from IGA surface to FE mesh
    if (isMappingIGA2FEM) {
        /// 1i. Loop through all nodes of master surface (FE)
        for (int counterNodeFE = 0; counterNodeFE < meshFE->numNodes; ++counterNodeFE) {
            // get patch index
            int patchIndex = projectedCoords[counterNodeFE].begin()->first;
            
            // assign patch to local variable
            IGAPatchSurface *thePatch = meshIGA->getSurfacePatch(patchIndex);

            // get number of CPs in this patch, pol. degrees for u and v
            int numCPsInPatch = thePatch->getNoControlPoints();
            int pDegree = thePatch->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree();
            int qDegree = thePatch->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree();
            int nShapeFuncsIGA = (pDegree + 1) * (qDegree + 1);

            // projectedCoords[counterNodeFE].begin()->second[0] gives parametrical coordinates
            double u = projectedCoords[counterNodeFE].begin()->second[0];
            double v = projectedCoords[counterNodeFE].begin()->second[1];
            int knotSpanU = thePatch->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(u);
            int knotSpanV = thePatch->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(v);

            // get value of basis fcts at point (u, v)
            double basisFctsMaster[nShapeFuncsIGA];
            thePatch->getIGABasis()->computeLocalBasisFunctions(basisFctsMaster, u, knotSpanU, v, knotSpanV);

            // get the cp net (point 1, point 2, ...)
            IGAControlPoint **cpNet = thePatch->getControlPointNet();
            int dofIGA[nShapeFuncsIGA];
            thePatch->getIGABasis()->getBasisFunctionsIndex(knotSpanU, knotSpanV, dofIGA);

            /// 1ii. Compute the coupling matrix C_M
            for (int i = 0; i < nShapeFuncsIGA; ++i) {
                int position = cpNet[dofIGA[i]]->getDofIndex();
                (*C_M)(counterNodeFE, position) += basisFctsMaster[i];
            }
        }
    } else {
        /// 2. Mapping is from FE mesh to IGA surface
        int numNodesIGA = meshIGA->getNumNodes();
        int numNodesFEM = meshFE->numNodes;

        /// 2i. Loop over all IGA surface control points
        for (int counterCP = 0; counterCP < numNodesIGA; ++counterCP) {
            if(!isProjectionOfCpOnTrimmed[counterCP]) {
                /// Build left hand side matrix
                // get patch index
                int patchIndex = projectedCPs[counterCP].begin()->first;
                
                // assign patch to local variable
                IGAPatchSurface *thePatch = meshIGA->getSurfacePatch(patchIndex);

                // get number of CPs in this patch, pol. degrees for u and v
                int numCPsInPatch = thePatch->getNoControlPoints();
                int pDegree = thePatch->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree();
                int qDegree = thePatch->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree();
                int nShapeFuncsIGA = (pDegree + 1) * (qDegree + 1);

                // projectedCPs[counterCP].begin()->second[0] gives parametrical coordinates
                double u = projectedCPs[counterCP].begin()->second[0];
                double v = projectedCPs[counterCP].begin()->second[1];
                int knotSpanU = thePatch->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(u);
                int knotSpanV = thePatch->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(v);

                // get value of basis fcts at point (u, v)
                double basisFctsMaster[nShapeFuncsIGA];
                thePatch->getIGABasis()->computeLocalBasisFunctions(basisFctsMaster, u, knotSpanU, v, knotSpanV);

                // get the cp net (point 1, point 2, ...)
                IGAControlPoint **cpNet = thePatch->getControlPointNet();
                int dofIGA[nShapeFuncsIGA];
                thePatch->getIGABasis()->getBasisFunctionsIndex(knotSpanU, knotSpanV, dofIGA);

                /// 2ii. Compute the coupling matrix C_L
                for (int i = 0; i < nShapeFuncsIGA; ++i) {
                    int position1 = counterCP;
                    int position2 = cpNet[dofIGA[i]]->getDofIndex();
                    double value = basisFctsMaster[i];
                    (*C_L)(position1, position2) = value;
                    (*C_LT)(position2, position1) = value;
                }

                /// 2iii. Compute the coupling matrix C_R
                int elemIndex = projectedCPsOnFEMesh[counterCP].first;
                double weights[ numNodesPerNeighborElem[counterCP] ];
                if (numNodesPerNeighborElem[counterCP] == 4) {
                    double localCoords[2] = {projectedCPsOnFEMesh[counterCP].second[0], projectedCPsOnFEMesh[counterCP].second[1]};
                    EMPIRE::MathLibrary::computeShapeFuncOfQuad(localCoords, weights);
                } else if (numNodesPerNeighborElem[counterCP] == 3) {
                    for (int k = 0; k < 3; ++k) {
                        weights[k] = projectedCPsOnFEMesh[counterCP].second[k];
                    }
                } else {
                    assert(false);
                }
                for (int k = 0; k < numNodesPerNeighborElem[counterCP]; ++k) {
                    int position1 = counterCP;
                    int position2 = meshFEDirectElemTable[elemIndex][k];
                    (*C_R)(position1, position2) = weights[k];
                }
            } else {
                (*C_L)(counterCP, counterCP) = 1.0; // if CP is projected on a trimmed region, ignore it
                (*C_LT)(counterCP, counterCP) = 1.0;
                (*C_R)(counterCP, 0) = 1.0; // this basically assigns the stress of the first FEM node to every CP that's ignored; required for
                                            // consistency of mapping
            }
        }
    }

    time(&timeEnd);
    INFO_OUT() << "Computing coupling matrices done! It took " << difftime(timeEnd, timeStart) << " seconds." << endl;
}

void IGABarycentricMapper::consistentMapping(const double* _slaveField, double *_masterField) {
    /*
     * Mapping of the
     * x_master = C_M * x_slave (mapping IGA to FEM)
     * or
     * x_master = (C_L)^(-1) * C_R * x_slave (mapping FEM to IGA)
     */

    if(isMappingIGA2FEM) {
        C_M->mulitplyVec(false, const_cast<double *>(_slaveField), _masterField, numNodesMaster);
    } else {
        double* tmpVec = new double[numNodesMaster];
        C_R->mulitplyVec(false, const_cast<double *>(_slaveField), tmpVec, numNodesMaster);
        C_L->solve(_masterField, tmpVec);

        delete[] tmpVec;
    }
}

void IGABarycentricMapper::conservativeMapping(const double* _masterField, double *_slaveField) {
    /*
     * Mapping of the
     * f_slave = (C_M)^T * f_master (mapping IGA to FEM)
     * or
     * f_slave = (C_L^(-1) * C_R)^T * f_master (mapping FEM to IGA)
     */

    if(isMappingIGA2FEM) {
        C_M->transposeMulitplyVec(const_cast<double *>(_masterField), _slaveField, numNodesMaster);
    } else {
        double* tmpVec = new double[numNodesMaster];
        C_LT->solve(tmpVec, const_cast<double *>(_masterField));
        C_R->transposeMulitplyVec(tmpVec, _slaveField, numNodesMaster);

        delete[] tmpVec;
    }
}

void IGABarycentricMapper::writeProjectedNodesOntoIGAMesh() {
    // Initializations
    IGAPatchSurface* IGAPatch;
    int numXiKnots, numEtaKnots, numXiControlPoints, numEtaControlPoints;
    double xiKnot, etaKnot;

    // Open file for writing the projected nodes
    ofstream projectedNodesFile;
    const string UNDERSCORE = "_";
    string projectedNodesFileName = name + "_projectedNodesOntoNURBSSurface.m";
    projectedNodesFile.open(projectedNodesFileName.c_str());
    projectedNodesFile.precision(14);
    projectedNodesFile << std::dec;

    projectedNodesFile << HEADER_DECLARATION << endl << endl;

    // Loop over all the patches
    for (int patchCounter = 0; patchCounter < meshIGA->getNumPatches(); patchCounter++) {
        // Get the IGA patch
        IGAPatch = meshIGA->getSurfacePatches()[patchCounter];

        // Get number of knots in each parametric direction
        numXiKnots = IGAPatch->getIGABasis()->getUBSplineBasis1D()->getNoKnots();
        numEtaKnots = IGAPatch->getIGABasis()->getVBSplineBasis1D()->getNoKnots();

        // Write out the knot vector in u-direction
        projectedNodesFile << "Patch" << patchCounter << endl << endl;
        projectedNodesFile << "xiKnotVector" << endl;

        for (int xiCounter = 0; xiCounter < numXiKnots; xiCounter++) {
            xiKnot = IGAPatch->getIGABasis()->getUBSplineBasis1D()->getKnotVector()[xiCounter];
            projectedNodesFile << xiKnot << " ";
        }
        projectedNodesFile << endl << endl;

        // Write out the knot vector in v-direction
        projectedNodesFile << "etaKnotVector" << endl;
        for (int etaCounter = 0; etaCounter < numEtaKnots; etaCounter++) {
            etaKnot = IGAPatch->getIGABasis()->getVBSplineBasis1D()->getKnotVector()[etaCounter];
            projectedNodesFile << etaKnot << " ";
        }
        projectedNodesFile << endl << endl;

        // Loop over all the nodes
        for (int nodeIndex = 0; nodeIndex < meshFE->numNodes; nodeIndex++) {
            // Loop over all the patches on which this node has been successfully projected
            for (map<int, vector<double> >::iterator it = projectedCoords[nodeIndex].begin();
                    it != projectedCoords[nodeIndex].end(); it++)
                if (it->first == patchCounter) {
                    projectedNodesFile << nodeIndex << "\t" << it->first << "\t" << it->second[0]
                            << "\t" << it->second[1] << endl;
                }
        }
        projectedNodesFile << endl;
    }

    // Close file
    projectedNodesFile.close();
}

void IGABarycentricMapper::writeProjectedCPsOntoFEMMesh() {
    // Initializations
    int numXiKnots, numEtaKnots, numXiControlPoints, numEtaControlPoints;

    // Open file for writing the projected nodes
    ofstream projectedCPsFile;
    const string UNDERSCORE = "_";
    string projectedCPsFileName = name + "_projectedCPOntoFEMElements.m";
    projectedCPsFile.open(projectedCPsFileName.c_str());
    projectedCPsFile.precision(14);
    projectedCPsFile << std::dec;

    projectedCPsFile << HEADER_DECLARATION << endl << endl;

    projectedCPsFile << "CP\tElement\tu\tv\n";

    // Loop over all the control points
    for (int counterCP = 0; counterCP < meshIGACpNet.size(); counterCP++) {
        if(!isProjectionOfCpOnTrimmed[counterCP]) {
            // Loop over all the patches on which this node has been successfully projected
            projectedCPsFile << counterCP << "\t" << projectedCPsOnFEMesh[counterCP].first << "\t" << projectedCPsOnFEMesh[counterCP].second[0]
                << "\t" << projectedCPsOnFEMesh[counterCP].second[1] << endl;
        }
    }

    projectedCPsFile << endl;

    // Close file
    projectedCPsFile.close();
}

void IGABarycentricMapper::printCouplingMatrices() {
    if (isMappingIGA2FEM) {
        ERROR_OUT() << "C_M" << endl;
        C_M->printCSR();
    } else {
        ERROR_OUT() << "C_L" << endl;
        C_L->printCSR();
        ERROR_OUT() << "C_R" << endl;
        C_R->printCSR();
    }
}

void IGABarycentricMapper::writeCouplingMatricesToFile() {
	DEBUG_OUT()<<"### Printing matrices into file ###"<<endl;
    if (isMappingIGA2FEM) {
    	DEBUG_OUT()<<"Size of C_M is "<<numNodesMaster<<" by "<<numNodesSlave<<endl;
        // if(Message::isDebugMode())
    		C_M->printCSRToFile(name + "_Cm.dat",1);
    } else {
        DEBUG_OUT()<<"Size of C_L is "<<numNodesMaster<<" by "<<numNodesMaster<<endl;
        DEBUG_OUT()<<"Size of C_R is "<<numNodesMaster<<" by "<<numNodesSlave<<endl;
        // if(Message::isDebugMode()) {
            C_L->printCSRToFile(name + "_Cl.dat",1);
            C_R->printCSRToFile(name + "_Cr.dat",1);
        // }
    }
}

void IGABarycentricMapper::checkConsistency() {
    double unitField[numNodesSlave];
    for(int i=0;i<numNodesSlave;i++) {
    	unitField[i]=1.0;
    }
    double output[numNodesMaster];
    this->consistentMapping(unitField,output);
    double norm=0;
    for(int i=0;i<numNodesMaster;i++) {
    	norm+=output[i]*output[i];
    }
    int denom = numNodesMaster;
    norm=sqrt(norm/denom);

    DEBUG_OUT()<<"### Check consistency ###"<<endl;
    DEBUG_OUT()<<"Norm of output field = "<<norm<<endl;
    /// WARNING hard coded tolerance. Used to decide if mapping is valid or not
    if(fabs(norm-1.0)>1e-1) {
    	ERROR_OUT()<<"Coupling not consistent !"<<endl;
    	ERROR_OUT()<<"Coupling of unit field deviating from 1 by "<<fabs(norm-1.0)<<endl;
    	exit(-1);
    }
}

} // namespace EMPIRE