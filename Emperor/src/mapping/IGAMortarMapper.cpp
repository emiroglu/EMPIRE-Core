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
/*
 * IGAMortarMapper.cpp
 *
 *  Created on: May 8, 2013
 *      Author: chenshen
 */

#include "IGAMortarMapper.h"
#include "IGAPatchSurface.h"
#include "IGAMesh.h"
#include "FEMesh.h"
#include "ClipperInterface.h"
#include "MathLibrary.h"
#include "DataField.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <set>
#include <algorithm>

using namespace std;

namespace EMPIRE {

/// Declaration statement
static const string HEADER_DECLARATION = "Author: Andreas Apostolatos";

IGAMortarMapper::IGAMortarMapper(std::string _name, IGAMesh *_meshIGA, FEMesh *_meshFE,
        double _disTol, int _numGPsTri, int _numGPsQuad, bool _isMappingIGA2FEM, int _numDivision) :
        name(_name), meshIGA(_meshIGA), disTol(_disTol), numGPsTri(_numGPsTri), numGPsQuad(
                _numGPsQuad), isMappingIGA2FEM(_isMappingIGA2FEM),numDivision(_numDivision) {

    assert(_meshIGA != NULL);
    assert(_meshFE != NULL);
    assert(_meshIGA->type == EMPIRE_Mesh_IGAMesh);
    assert(_meshFE->type == EMPIRE_Mesh_FEMesh);

    DEBUG_OUT()<<"------------------------------"<<endl;
    DEBUG_OUT()<<"DEBUG for mapper "<<_name<<endl;
    DEBUG_OUT()<<"------------------------------"<<endl;
    DEBUG_OUT()<<"numDivision is "<<numDivision<<endl;

    if (_meshFE->triangulate() == NULL)
        meshFE = _meshFE;
    else
        meshFE = _meshFE->triangulate();

    projectedCoords = new vector<map<int, double*> >(meshFE->numNodes);

    if (isMappingIGA2FEM) {
        numNodesSlave = meshIGA->getNumNodes();
        numNodesMaster = meshFE->numNodes;
    } else {
        numNodesSlave = meshFE->numNodes;
        numNodesMaster = meshIGA->getNumNodes();
    }
    DEBUG_OUT()<<"NumNodesIGA is "<<meshIGA->getNumNodes()<<endl;
    DEBUG_OUT()<<"NumNodesFE is "<<meshFE->numNodes<<endl;
    DEBUG_OUT()<<"NumNodesMaster is "<<numNodesMaster<<endl;
    DEBUG_OUT()<<"NumNodesSlave is "<<numNodesSlave<<endl;


    C_NR = new MathLibrary::SparseMatrix<double>(numNodesMaster, numNodesSlave);
    C_NN = new MathLibrary::SparseMatrix<double>(numNodesMaster, true);

    gaussTriangle = new MathLibrary::IGAGaussQuadratureOnTriangle(numGPsTri);
    gaussQuad = new MathLibrary::IGAGaussQuadratureOnQuad(numGPsQuad);

    initTables();

    projectPointsToSurface();

    // Write the projected points on to a file
    writeProjectedNodesOntoIGAMesh();

    computeCouplingMatrices();

    printCouplingMatricesToFile();

    C_NN->factorize();

    checkConsistency();
}

IGAMortarMapper::~IGAMortarMapper() {

    for (int i = 0; i < meshFE->numElems; i++)
        delete[] meshFEDirectElemTable[i];
    delete[] meshFEDirectElemTable;

    for (int nodeIndex = 0; nodeIndex < meshFE->numNodes; nodeIndex++) {
        for (int patchCount = 0; patchCount < meshIGA->getSurfacePatches().size(); patchCount++) {
            if ((*projectedCoords)[nodeIndex].find(patchCount)
                    != (*projectedCoords)[nodeIndex].end()) {
                delete[] (*projectedCoords)[nodeIndex][patchCount];
            }
        }
    }
    delete projectedCoords;

    delete gaussTriangle;
    delete gaussQuad;

    C_NN->cleanPardiso();
    delete C_NR;
    delete C_NN;
}

void IGAMortarMapper::initTables() {
    /* using the map to store the nodeIDs
     * but here the "key" is the node ID, and the value is the position in nodeIDs
     * the map is sorted automatically, so it is efficient for searching
     */

// compute direct element table for fluid mesh
    meshFEDirectElemTable = new int*[meshFE->numElems]; // deleted
    for (int i = 0; i < meshFE->numElems; i++)
        meshFEDirectElemTable[i] = new int[meshFE->numNodesPerElem[i]];

    map<int, int> *meshFENodesMap = new map<int, int>(); // deleted
    for (int i = 0; i < meshFE->numNodes; i++)
        meshFENodesMap->insert(meshFENodesMap->end(), pair<int, int>(meshFE->nodeIDs[i], i));
    int count = 0;

    for (int i = 0; i < meshFE->numElems; i++) {
        const int numNodesPerElem = meshFE->numNodesPerElem[i];

        for (int j = 0; j < numNodesPerElem; j++) {
            if (meshFENodesMap->find(meshFE->elems[count + j]) == meshFENodesMap->end()) {
                ERROR_OUT() << "Cannot find node ID " << meshFE->elems[count + j] << endl;
                exit(-1);
            }
            meshFEDirectElemTable[i][j] = meshFENodesMap->at(meshFE->elems[count + j]);
        }
        count += numNodesPerElem;
    }

    delete meshFENodesMap;

}

void IGAMortarMapper::projectPointsToSurface() {
    /*  Projects all nodes of the FE side onto the IGA mesh.
     *
     *  0. Read input
     *
     *  1. Loop over all the patches in the IGA mesh
     *     1i. Get the IGA patch
     *    1ii. Initialize all projection flags to false
     *   1iii. Loop over all the elements on the fluid side
     *         1iii.1. Initialize the flag to false and the node id to zero
     *         1iii.2. Loop over all nodes of the current element to check if there exist one node has been projected already
     *                 1iii.2i. If the node has already been projected set projection flag to true
     *                1iii.2ii. Get the global ID of the projected node
     *               1iii.2iii. Break the loop
     *         1iii.3. Check if there exist one node in the current element has been successfully projected
     *                 1iii.3i. If so, use result of the projected node as the initial guess for the projection step
     *                1iii.3ii. Otherwise, find the nearest knot intersection as initial guess for the projection step
     *         1iii.4. Loop over each node at the current element in the FE side
     *                 1iii.4i. Get the node ID from the EFT
     *                1iii.4ii. Check if the node has been already projected. If not, compute the point projection on the IGA patch using the initial guess get from the last step the results are stored in the class member "projectedCoords"
     *                          1iii.4ii.1. get the Cartesian coordinates of the node in the FE side
     *                          1iii.4ii.2. Get an initial guess for the parametric location of the projected node of the FE side on the NURBS patch
     *                          1iii.4ii.3. Compute point projection on the NURBS patch using the Newton-Rapshon iteration method
     *                          1iii.4ii.4. Check if the Newton-Rapshon iterations have converged and if the points are coinciding
     *                                      1iii.4ii.4i. Set projection flag to true
     *                                     1iii.4ii.4ii. Insert the parametric locations of the projected FE nodes into the map projectedCoords
     *
     *   2. Loop over all the nodes in the FE side
     *      2i. Initialize projection flag to false
     *     2ii. Loop over all the patches in the IGA mesh and check if the node is has been projected to any patch
     *    2iii. If the node has not been projected to any patch of the IGA mesh, loop over all the patches in the mesh again and try to project it with better initial guess
     *          2iii.1. Get the NURBS patch
     *          2iii.2. Get the Cartesian coordinates of the node in the FE side
     *          2iii.3. Get a better initial guess for the Newton-Rapshon iterations using the variable REFINED_NUM_PARAMETRIC_LOCATIONS
     *          2iii.4. Compute point projection on the NURBS patch using the Newton-Rapshon iteration method
     *          2iii.5. Check if the Newton-Rapshon iterations have converged and if the points are coinciding
     *                  2iii.5i. Set projection flag to true
     *                 2iii.5ii. Insert the parametric locations of the projected FE nodes into the map projectedCoords
     *           2iv. If the node can still not be projected assert an error in the projection phase
     *
     *  3. Deallocate dynamically allocated memory
     */

    /// 0. Read input
    // Initialization of variables
    // Array of booleans containing flags on the projection of the FE nodes onto the NURBS patch
    bool *isProjected = new bool[meshFE->numNodes]; // deleted

    // Flag on whether a node is projected on the IGA mesh
    bool isNodeProjected = false;

    // Flag (??)
    bool isNodeInsideElementProjecteded;

    // id of the projected node
    int projectedNode;

    // Coordinates of the projected nodes on the NURBS patch
    double initialU, initialV;

    // Initialize the parametric coordinates of the FE node on the NURBS patch
    double projectedU;
    double projectedV;

    // Initialize flag on the convergence of the Newton-Rapshon iterations for the projection of a node on the NURBS patch
    bool isConverged, isConvergedInside;

    // Initialize the array of the Cartesian coordinates of a node in the FE side
    double cartesianCoords[3];

    // Initialize the node ID
    int nodeIndex;

    /// Initialize parametric coordinates of the projected node onto the NURBS patch
    double U, V;
    double P[3];

    // Get the number of patches in the IGA mesh
    int numPatches = meshIGA->getSurfacePatches().size();

    /// 1. Loop over all the patches in the IGA mesh
    for (int patchCount = 0; patchCount < numPatches; patchCount++) {
        /// 1i. Get the IGA patch
        IGAPatchSurface* thePatch = meshIGA->getSurfacePatches()[patchCount];

        /// 1ii. Initialize all projection flags to false
        for (int i = 0; i < meshFE->numNodes; i++)
            isProjected[i] = false;

        /// 1iii. Loop over all the elements on the fluid side
        for (int i = 0; i < meshFE->numElems; i++) {

            /// 1iii.1. Initialize the flag to false and the node id to zero
            isNodeInsideElementProjecteded = false;
            projectedNode = 0;

            /// 1iii.2. Loop over all nodes of the current element to check if there exist one node has been projected already
            for (int j = 0; j < meshFE->numNodesPerElem[i]; j++)
                if (isProjected[meshFEDirectElemTable[i][j]]) {
                    /// 1iii.2i. If the node has already been projected set projection flag to true
                    isNodeInsideElementProjecteded = true;

                    /// 1iii.2ii. Get the global ID of the projected node
                    projectedNode = meshFEDirectElemTable[i][j];

                    /// 1iii.2iii. Break the loop
                    break;
                }

            /// 1iii.3. Check if there exist one node in the current element has been successfully projected
            if (isNodeInsideElementProjecteded) {
                /// 1iii.3i. If so, use result of the projected node as the initial guess for the projection step
                initialU = (*projectedCoords)[projectedNode][patchCount][0];
                initialV = (*projectedCoords)[projectedNode][patchCount][1];

            } else {
                /// 1iii.3ii. Otherwise, find the nearest knot intersection as initial guess for the projection step

                // Get the node ID from the EFT
                nodeIndex = meshFEDirectElemTable[i][0];

                // Get the Cartesian coordinates of that node
                cartesianCoords[0] = meshFE->nodes[nodeIndex * 3];
                cartesianCoords[1] = meshFE->nodes[nodeIndex * 3 + 1];
                cartesianCoords[2] = meshFE->nodes[nodeIndex * 3 + 2];

                // Get accordingly an initial guess for the projection onto the NURBS patch
                thePatch->findInitialGuess4PointProjection(initialU, initialV, cartesianCoords);
            }

            /// 1iii.4. Loop over each node at the current element in the FE side
            for (int j = 0; j < meshFE->numNodesPerElem[i]; j++) {

                // 1iii.4i. Get the node ID from the EFT
                nodeIndex = meshFEDirectElemTable[i][j];

                /// 1iii.4ii. Check if the node has been already projected. If not, compute the point projection on the IGA patch using the initial guess get from the last step the results are stored in the class member "projectedCoords"
                if (!isProjected[nodeIndex]) {

                    /// 1iii.4ii.1. get the Cartesian coordinates of the node in the FE side
                    cartesianCoords[0] = meshFE->nodes[nodeIndex * 3];
                    cartesianCoords[1] = meshFE->nodes[nodeIndex * 3 + 1];
                    cartesianCoords[2] = meshFE->nodes[nodeIndex * 3 + 2];

                    /// 1iii.4ii.2. Get an initial guess for the parametric location of the projected node of the FE side on the NURBS patch
                    projectedU = initialU;
                    projectedV = initialV;

                    /// 1iii.4ii.3. Compute point projection on the NURBS patch using the Newton-Rapshon iteration method

                    // isConverged --> projected onto patch
                    isConvergedInside = thePatch->computePointProjectionOnPatch(projectedU,
                            projectedV, cartesianCoords, isConverged);

                    /// 1iii.4ii.4. Check if the Newton-Rapshon iterations have converged
                    if (isConvergedInside
                            && MathLibrary::computePointDistance(&meshFE->nodes[nodeIndex * 3],
                                    cartesianCoords) < disTol) {

                        /// 1iii.4ii.4i. Set projection flag to true
                        isProjected[nodeIndex] = true;

                        /// 1iii.4ii.4ii. Insert the parametric locations of the projected FE nodes into the map projectedCoords
                        double* coordTmp = new double(2);
                        coordTmp[0] = projectedU;
                        coordTmp[1] = projectedV;
                        (*projectedCoords)[nodeIndex].insert(
                                std::pair<int, double*>(patchCount, coordTmp));

                    }
                }
            }
        }
    }

    /// 2. Loop over all the nodes in the FE side
    for (int nodeIndex = 0; nodeIndex < meshFE->numNodes; nodeIndex++) {
        /// 2i.1 Initialize projection flag to false
        isNodeProjected = false;
        /// 2ii. Loop over all the patches in the IGA mesh and check if the node is has been projected to any patch
        for (int patchCount = 0; patchCount < numPatches; patchCount++) {
            if ((*projectedCoords)[nodeIndex].find(patchCount)
                    != (*projectedCoords)[nodeIndex].end()) {
                /// 2ii.1. If the node has been projected to a patch in the IGA mesh, set the flag to true
                isNodeProjected = true;

                /// 2ii.2. Break the loop
                break;
            }
        }

        /// 2iii. If the node has not been projected to any patch of the IGA mesh, loop over all the patches in the mesh again and try to project it with better initial guess
        if (!isNodeProjected) {
            /// 2iii.0 Get the Cartesian coordinates of the node in the FE side
            P[0] = meshFE->nodes[nodeIndex * 3];
            P[1] = meshFE->nodes[nodeIndex * 3 + 1];
            P[2] = meshFE->nodes[nodeIndex * 3 + 2];
            for (int patchCount = 0; patchCount < numPatches; patchCount++) {
                /// 2iii.1. Get the NURBS patch
                IGAPatchSurface* thePatch = meshIGA->getSurfacePatches()[patchCount];

                /// 2iii.2. Get the Cartesian coordinates of the node in the FE side
                double P_out[3];
                P_out[0] = P[0];
                P_out[1] = P[1];
				P_out[2] = P[2];

                /// 2iii.3. Get a better initial guess for the Newton-Rapshon iterations using the variable REFINED_NUM_PARAMETRIC_LOCATIONS
                thePatch->findInitialGuess4PointProjection(U, V, P,
                        REFINED_NUM_PARAMETRIC_LOCATIONS, REFINED_NUM_PARAMETRIC_LOCATIONS);

                /// 2iii.4. Compute point projection on the NURBS patch using the Newton-Rapshon iteration method
                isConvergedInside = thePatch->computePointProjectionOnPatch(U, V, P_out, isConverged);

                /// 2iii.5. Check if the Newton-Rapshon iterations have converged and if the points are coinciding
                bool hasConverged=isConvergedInside;
                if (hasConverged && MathLibrary::computePointDistance(P, P_out) < disTol) {
                    /// 2iii.5i. Set projection flag to true
                    isNodeProjected = true;

                    /// 2iii.5ii. Insert the parametric locations of the projected FE nodes into the map projectedCoords
                    double* coordTmp = new double(2);
                    coordTmp[0] = U;
                    coordTmp[1] = V;
                    (*projectedCoords)[nodeIndex].insert(make_pair(patchCount, coordTmp));
                }
            }

            /// 2iv. If the node can still not be projected assert an error in the projection phase
            if (!isNodeProjected) {
                ERROR_OUT() << " in IGAMortarMapper::projectPointsToSurface" << endl;
                ERROR_OUT() << "Cannot project node: " << nodeIndex << "  ("
                        << P[0] << ", " << P[1] << ", " << P[2] << ")" << endl;
                ERROR_OUT() << "Projection failed in IGA mapper " << name << endl;
                exit (EXIT_FAILURE);
            }
        }
    }

    /// 3. Deallocate dynamically allocated memory
    delete[] isProjected;

}

void IGAMortarMapper::computeCouplingMatrices() {
    /*
     * Computes the coupling matrices CNR and CNN.
     * Loop over all the elements in the FE side
     * ->
     * 1. Find whether the projected FE element is located on one patch or split
     *
     * 2. Compute the coupling matrices
     * ->
     * 2i. Loop over patches where element can be projected entirely on one patch
     * 2ii. Loop over patches where element is split
     * <-
     */

	//List of integrated element
	set<int> elementIntegrated;
/// The vertices of the canonical polygons
    double parentTriangle[6] = { 0, 0, 1, 0, 0, 1 };
    double parentQuadriliteral[8] = { -1, -1, 1, -1, 1, 1, -1, 1 };

/// Loop over all the elements in the FE side
    for (int elemCount = 0; elemCount < meshFE->numElems; elemCount++) {
    	DEBUG_OUT()<<"######################"<<endl;
    	DEBUG_OUT()<<"# ELEMENT ["<<elemCount<<"] #"<<endl;
    	DEBUG_OUT()<<"######################"<<endl;

        // Compute the number of shape functions. Depending on number of nodes in the current element
        int numNodesElementFE = meshFE->numNodesPerElem[elemCount];

        // The vertices of the FE element in the parent domain
        double* projectedElementFEWZ;

        // Decide whether the parent element is a triangle or a quadrilateral
        if (numNodesElementFE == 3)
            projectedElementFEWZ = parentTriangle;
        else
            projectedElementFEWZ = parentQuadriliteral;

        /// 1. Find whether the projected FE element is located on one patch

        // Initialize the patch that possibly contains that element
        IGAPatchSurface* thePatch = NULL;

        // Initialize the index of this patch
        int patchIndex = 0;

        set<int> patchWithFullElt;
        set<int> patchWithSplitElt;
        getPatchesIndexElementIsOn(elemCount, patchWithFullElt, patchWithSplitElt);
        DEBUG_OUT()<<"Number of patch in set FULL "<<patchWithFullElt.size()<<endl;
        DEBUG_OUT()<<"Number of patch in set SPLIT "<<patchWithSplitElt.size()<<endl;
        /////////////////////////////////////
        /// Compute the coupling matrices ///
        /////////////////////////////////////
        /// 1. If the current element can be projected on one patch
		for (set<int>::iterator it = patchWithFullElt.begin();
				it != patchWithFullElt.end(); ++it) {
			int patchIndex=*it;
			IGAPatchSurface* thePatch = meshIGA->getSurfacePatches()[patchIndex];
			/// Get the projected coordinates for the current element
			Polygon2D polygonUV;
			map<int,Polygon2D> extraPolygonUV;
			bool isProjectedOnPatchBoundary=true;
			bool isProjectedOnPatch=false;
			double u, v;
			double div, dis;
			int MAX_DIV=numDivision;
			/// 1.1 Get initial polygon from projection
			// For every point of polygon
			for (int nodeCounter = 0; nodeCounter < numNodesElementFE; nodeCounter++) {
				int nodeCount = (nodeCounter + 0) % numNodesElementFE;
				int nodeCountNext = (nodeCounter + 1) % numNodesElementFE;
				int nodeIndex = meshFEDirectElemTable[elemCount][nodeCount];
				int nodeIndexNext = meshFEDirectElemTable[elemCount][nodeCountNext];
				double* P1 = &(meshFE->nodes[nodeIndex * 3]);
				double* P2 = &(meshFE->nodes[nodeIndexNext * 3]);
				/// 1.1.1 First point, just add it to polygon
				u = (*projectedCoords)[nodeIndex][patchIndex][0];
				v = (*projectedCoords)[nodeIndex][patchIndex][1];
				polygonUV.push_back(make_pair(u,v));
				/// 1.1.2 Compute intermediate points on the line
				computeIntermediatePoints(patchIndex,elemCount,nodeIndex,nodeIndexNext,P1,P2,polygonUV,extraPolygonUV);
			}
			ClipperInterface::cleanPolygon(polygonUV);
			if(polygonUV.size()<3)
				continue;
			bool isIntegrated=false;
			ListPolygon2D listTrimmedPolygonUV(1, polygonUV);
			/// 1.2 Apply trimming
			if(thePatch->isTrimmed())
				clipByTrimming(thePatch,polygonUV,listTrimmedPolygonUV);
			/// 1.3 For each subelement output of the trimmed polygon, clip by knot span
			for(int trimmedPolygonIndex=0;trimmedPolygonIndex<listTrimmedPolygonUV.size();trimmedPolygonIndex++) {
				Polygon2D listSpan;
				ListPolygon2D listPolygonUV;
				/// 1.3.1 Clip by knot span
				clipByKnotSpan(thePatch,listTrimmedPolygonUV[trimmedPolygonIndex],listPolygonUV,listSpan);
				/// 1.3.2 For each subelement clipped by knot span, compute canonical element and integrate
				for(int index=0;index<listSpan.size();index++) {
					ClipperInterface::cleanPolygon(listPolygonUV[index]);
					if(listPolygonUV[index].size()<3)
						continue;
					isIntegrated=true;
					// Get canonical element
					Polygon2D polygonWZ = computeCanonicalElement(thePatch,listPolygonUV[index],elemCount);
					// Integrate
					integrate(thePatch,listPolygonUV[index],listSpan[index].first,listSpan[index].second,polygonWZ,elemCount);
				}
			}
			if(isIntegrated)
				elementIntegrated.insert(elemCount);

		}
        /// 2. If the current element is split in some patches
        // Loop over all the patches in the IGA Mesh having a part of the FE element projected inside
		for (set<int>::iterator it = patchWithSplitElt.begin(); it != patchWithSplitElt.end(); it++) {
			int patchIndex=*it;
			IGAPatchSurface* thePatch = meshIGA->getSurfacePatches()[patchIndex];

			// Cartesian coordinates of the low order element
			double elementFEXYZ[12];
			for (int i = 0; i < numNodesElementFE; i++) {
				int nodeIndex = meshFEDirectElemTable[elemCount][i];
				for (int j = 0; j < 3; j++)
					elementFEXYZ[i * 3 + j] = meshFE->nodes[nodeIndex * 3 + j];
			}
			// Stores points of the polygon clipped by the nurbs patch
			Polygon2D polygonUV;
			map<int,Polygon2D> extraPolygonUV;
			int nodeCounter=0;
			int numEdgesProcessed=0;
			int numEdges=numNodesElementFE;
			bool isProjectedOnPatchBoundary=true;
			bool isProjectedOnPatch=false;
			double u, v;
			double w, z;
			double div, dis;
			char edge1,edge2;
			/// 2.1 Compute initial polygon by computing point on ever edge case
			while(numEdgesProcessed!=numEdges) {
				// Get the node indexes
				int nodeCount = (nodeCounter + 0) % numNodesElementFE;
				int nodeCountNext = (nodeCounter + 1) % numNodesElementFE;
				int nodeIndex = meshFEDirectElemTable[elemCount][nodeCount];
				int nodeIndexNext = meshFEDirectElemTable[elemCount][nodeCountNext];

				// Flag on whether the node is inside the current NURBS patch
				bool isNodeInsidePatch = (*projectedCoords)[nodeIndex].find(patchIndex)
							!= (*projectedCoords)[nodeIndex].end();
				bool isNextNodeInsidePatch = (*projectedCoords)[nodeIndexNext].find(patchIndex)
							!= (*projectedCoords)[nodeIndexNext].end();

				// Get the 2 points of the line
				double* P1 = &(meshFE->nodes[nodeIndex * 3]);
				double* P2 = &(meshFE->nodes[nodeIndexNext * 3]);
				double P[3];

				// Loop until starting predefined starting point 1st point outside and 2nd point outside
				if(numEdgesProcessed==0 && (isNodeInsidePatch || !isNextNodeInsidePatch)) {
					nodeCounter++;
					continue;
				}
				/// 2.1.1 NODE.1 OUTSIDE AND NODE.2 INSIDE
				if(!isNodeInsidePatch && isNextNodeInsidePatch) {
					// 1. First point
					u = (*projectedCoords)[nodeIndexNext][patchIndex][0];
					v = (*projectedCoords)[nodeIndexNext][patchIndex][1];
					edge1=thePatch->computePointProjectionOnPatchBoundary_NewtonRhapson(u, v, div, dis, P2, P1);
					isProjectedOnPatchBoundary = edge1;
					if(isProjectedOnPatchBoundary && dis <= disTol) {
						polygonUV.push_back(make_pair(u,v));
					}
					// 2. Compute intermediate points
					for(int i=0;i<3;i++) {
						P[i]=P2[i]*(1-div)+P1[i]*div;
					}
					computeIntermediatePoints(patchIndex,elemCount,nodeIndexNext,nodeIndex,P,P2,polygonUV,extraPolygonUV);
					// 3. Second point
					u = (*projectedCoords)[nodeIndexNext][patchIndex][0];
					v = (*projectedCoords)[nodeIndexNext][patchIndex][1];
					polygonUV.push_back(make_pair(u,v));
					// 4. Update processed edge
					numEdgesProcessed++;

				}
				/// 2.1.2 NODE.1 INSIDE AND NODE.2 INSIDE
				if(isNodeInsidePatch && isNextNodeInsidePatch) {
					// 1. First point
					u = (*projectedCoords)[nodeIndex][patchIndex][0];
					v = (*projectedCoords)[nodeIndex][patchIndex][1];
					polygonUV.push_back(make_pair(u,v));
					// 2. Compute intermediate points
					computeIntermediatePoints(patchIndex,elemCount,nodeIndex,nodeIndexNext,P1,P2,polygonUV,extraPolygonUV);
					// 3. Second point
					u = (*projectedCoords)[nodeIndexNext][patchIndex][0];
					v = (*projectedCoords)[nodeIndexNext][patchIndex][1];
					polygonUV.push_back(make_pair(u,v));
					// 4. Update processed edge
					numEdgesProcessed++;
				}
				/// 2.1.3 NODE.1 INSIDE AND NODE.2 OUTSIDE
				if(isNodeInsidePatch && !isNextNodeInsidePatch) {
					// 1. First point
					u = (*projectedCoords)[nodeIndex][patchIndex][0];
					v = (*projectedCoords)[nodeIndex][patchIndex][1];
					polygonUV.push_back(make_pair(u,v));
					// 2. Second point, intersection
					u = (*projectedCoords)[nodeIndex][patchIndex][0];
					v = (*projectedCoords)[nodeIndex][patchIndex][1];
					edge2=thePatch->computePointProjectionOnPatchBoundary_NewtonRhapson(u, v, div, dis, P1, P2);
					isProjectedOnPatchBoundary = edge2;
					pair<double, double> uv_tmp,wz_tmp;
					if (isProjectedOnPatchBoundary && dis <= disTol) {
						uv_tmp=make_pair(u,v);
					}
					// 3. Compute intermediate points
					for(int i=0;i<3;i++) {
						P[i]=P1[i]*(1-div)+P2[i]*div;
					}
					computeIntermediatePoints(patchIndex,elemCount,nodeIndex,nodeIndexNext,P1,P,polygonUV,extraPolygonUV);
					// 3. Second point, store values
					polygonUV.push_back(uv_tmp);
					// Check if corners to be included
//					if(numDivision<=1) {
//						Polygon2D corners = thePatch->getCorner(edge1,edge2,ClipperInterface::isCounterclockwise(polygonUV));
//						polygonUV.insert(polygonUV.end(),corners.begin(),corners.end());
//					}
					// 4. Update processed edge
					numEdgesProcessed++;
				}
				/// 2.1.4 NODE.1 OUTSIDE AND NODE.2 OUTSIDE
				if(!isNodeInsidePatch && !isNextNodeInsidePatch) {
					// 1. Compute intermediate points
						computeIntermediatePoints(patchIndex,elemCount,nodeIndex,nodeIndexNext,P1,P2,polygonUV,extraPolygonUV);
					// 4. Update processed edge
					numEdgesProcessed++;
				}
				nodeCounter++;
				if(!isProjectedOnPatchBoundary) {
					if(thePatch->isTrimmed()) {
						WARNING_OUT() << "Warning in IGAMortarMapper::computeCouplingMatrices"
								<< endl;
						WARNING_OUT() << "Cannot find point projection on patch boundary" << endl;
						WARNING_OUT() << "But proceed anyway !" << endl;
						break;//break loop over node
					} else {
						ERROR_OUT() << "Error in IGAMortarMapper::computeCouplingMatrices"
								<< endl;
						ERROR_OUT() << "Cannot find point projection on patch boundary" << endl;
						ERROR_OUT()
								<< "Cannot find point projection on patch boundary between node ["
								<< nodeIndex << "]:(" << meshFE->nodes[nodeIndex * 3] << ","
								<< meshFE->nodes[nodeIndex * 3 + 1] << ","
								<< meshFE->nodes[nodeIndex * 3 + 2] << ") and node ["
								<< nodeIndexNext << "]:(" << meshFE->nodes[nodeIndexNext * 3]
								<< "," << meshFE->nodes[nodeIndexNext * 3 + 1] << ","
								<< meshFE->nodes[nodeIndexNext * 3 + 2] << ") on patch ["
								<< patchIndex << "] boundary" << endl;
						ERROR_OUT() << "Projection failed in IGA mapper " << name << endl;
						exit (EXIT_FAILURE);
					}
				}
			}
			ClipperInterface::cleanPolygon(polygonUV);
			// Proceed toward integration if the polygon is valid, i.e. at least a triangle
			if (polygonUV.size() >= 3) {
				bool isIntegrated=false;
				ListPolygon2D listTrimmedPolygonUV(1, polygonUV);
				/// Apply trimming window
				if(thePatch->isTrimmed())
					clipByTrimming(thePatch,polygonUV,listTrimmedPolygonUV);
				for(int trimmedPolygonIndex=0;trimmedPolygonIndex<listTrimmedPolygonUV.size();trimmedPolygonIndex++) {
					Polygon2D listSpan;
					ListPolygon2D listPolygonUV;
					/// Apply knot span window
					clipByKnotSpan(thePatch,listTrimmedPolygonUV[trimmedPolygonIndex],listPolygonUV,listSpan);
					for(int index=0;index<listSpan.size();index++) {
						ClipperInterface::cleanPolygon(listPolygonUV[index]);
						if(listPolygonUV[index].size()<3) continue;
						isIntegrated=true;
						/// Get FE canonical element
						Polygon2D polygonWZ = computeCanonicalElement(thePatch,listPolygonUV[index],elemCount);
						/// Integrate
						integrate(thePatch,listPolygonUV[index],listSpan[index].first,listSpan[index].second,polygonWZ,elemCount);
					}
				}
				if(isIntegrated)
					elementIntegrated.insert(elemCount);
			} //end of if polygon has more than 3 edges

//			DEBUG_OUT()<< "Number of polygon EXTRA "<<extraPolygonUV.size()<<endl;
//			/// Integrate the polygon appearing only on current patch from an external element.
//			/// Only the case when both points of an edge are outside and still the edge is partially projected in the patch
//			for (map<int,Polygon2D>::iterator it = extraPolygonUV.begin(); it != extraPolygonUV.end(); it++) {
//				debugPolygon(it->second,"extra polygon");
//				ClipperInterface::cleanPolygon(it->second);
//				// Proceed toward integration if the polygon is valid, i.e. at least a triangle
//				if (it->second.size() < 3)
//					continue;
//				bool isIntegrated=false;
//				ListPolygon2D listTrimmedPolygonUV(1, it->second);
//				/// Apply trimming window
//				if(thePatch->isTrimmed())
//					clipByTrimming(thePatch,it->second,listTrimmedPolygonUV);
//				for(int trimmedPolygonIndex=0;trimmedPolygonIndex<listTrimmedPolygonUV.size();trimmedPolygonIndex++) {
//					Polygon2D listSpan;
//					ListPolygon2D listPolygonUV;
//					/// Apply knot span window
//					clipByKnotSpan(thePatch,listTrimmedPolygonUV[trimmedPolygonIndex],listPolygonUV,listSpan);
//					for(int index=0;index<listSpan.size();index++) {
//						ClipperInterface::cleanPolygon(listPolygonUV[index]);
//						if(listPolygonUV[index].size()<3) continue;
//						isIntegrated=true;
//						/// Get FE canonical element for element in the map
//						Polygon2D polygonWZ = computeCanonicalElement(thePatch,listPolygonUV[index],it->first);
//						/// Integrate for element in the map
//						integrate(thePatch,listPolygonUV[index],listSpan[index].first,listSpan[index].second,polygonWZ,it->first);
//					}
//				}
//			}
		} // end of loop over set of split patch
    } // end of loop over all the element
    if(elementIntegrated.size()!=meshFE->numElems) {
    	ERROR_OUT()<<"Number of FE mesh integrated is "<<elementIntegrated.size()<<" over "<<meshFE->numElems<<endl;
    	for(int i=0;i<meshFE->numElems;i++) {
    		if(!elementIntegrated.count(i))
    			ERROR_OUT()<<"Missing element number "<<i<<endl;
    	}
    	ERROR_BLOCK_OUT("IGAMortarMapper","ComputeCouplingMatrices","Not all element in FE mesh integrated ! Coupling matrices invalid");
    	exit(-1);
    }
}


void IGAMortarMapper::getPatchesIndexElementIsOn(int elemIndex, set<int>& patchWithFullElt, set<int>& patchWithSplitElt) {
    // Initialize the flag whether the projected FE element is located on one patch
    bool isAllNodesOnPatch = true;
    // Initialize the flag whether the projected FE element is not located at all on the patch
    bool isAllNodesOut= true;
	// Loop over all the patches
	for (int patchCount = 0; patchCount < meshIGA->getSurfacePatches().size(); patchCount++) {
		isAllNodesOnPatch = true;
		isAllNodesOut= true;
		// Loop over all the nodes of the unclipped element
		for (int nodeCount = 0; nodeCount < meshFE->numNodesPerElem[elemIndex]; nodeCount++) {
			// Find the index of the node in the FE mesh
			int nodeIndex = meshFEDirectElemTable[elemIndex][nodeCount];
			// Find whether this index is in the projected nodes array
			bool isNodeOnPatch = (*projectedCoords)[nodeIndex].find(patchCount)
					!= (*projectedCoords)[nodeIndex].end();
			// Update flag
			if (!isNodeOnPatch) {
				isAllNodesOnPatch = false;
			} else {
				isAllNodesOut=false;
			}
		}
		// If all nodes are on the patch "patchCount", save this patch and go for next patch
		if (isAllNodesOnPatch) {
			patchWithFullElt.insert(patchCount);
			continue;
		}
		// If element is splitted for patch "patchCount", save this patch and go for next patch
		if(!isAllNodesOut) {
			patchWithSplitElt.insert(patchCount);
			continue;
		}
	}
}

void IGAMortarMapper::clipByTrimming(const IGAPatchSurface* _thePatch, const Polygon2D& _polygonUV, ListPolygon2D& _listPolygonUV) {
	ClipperInterface c;
	for(int loop=0;loop<_thePatch->getTrimming().getNumOfLoops();loop++) {
		const std::vector<double> clippingWindow=_thePatch->getTrimming().getLoop(loop).getPolylines();
		c.addPathClipper(clippingWindow);
	}
	c.setFilling(ClipperInterface::POSITIVE, 0);
	c.addPathSubject(_polygonUV);
	c.clip();
	c.getSolution(_listPolygonUV);
}

void IGAMortarMapper::clipByKnotSpan(const IGAPatchSurface* _thePatch, const Polygon2D& _polygonUV, ListPolygon2D& _listPolygon, Polygon2D& _listSpan) {
	/// 1.find the knot span which the current element located in.
	//      from minSpanu to maxSpanu in U-direction, and from minSpanV to max SpanV in V-direction
    double *knotVectorU = _thePatch->getIGABasis()->getUBSplineBasis1D()->getKnotVector();
    double *knotVectorV = _thePatch->getIGABasis()->getVBSplineBasis1D()->getKnotVector();

	int span[4];
	int isOnSameKnotSpan = computeKnotSpanOfProjElement(_thePatch, _polygonUV,span);
	int minSpanU = span[0];
	int maxSpanU = span[1];
	int minSpanV = span[2];
	int maxSpanV = span[3];
	/// If on same knot span then returned the same polygon as input
	if (isOnSameKnotSpan) {
		_listPolygon.push_back(_polygonUV);
		_listSpan.push_back(make_pair(minSpanU,minSpanV));
	/// Else clip the polygon for every knot span window it is crossing
	} else {
		for (int spanU = minSpanU; spanU <= maxSpanU; spanU++) {
			for (int spanV = minSpanV; spanV <= maxSpanV; spanV++) {
				ClipperInterface c(1e-12);
				if (knotVectorU[spanU] != knotVectorU[spanU + 1]
						&& knotVectorV[spanV] != knotVectorV[spanV + 1]) {
					Polygon2D knotSpanWindow(4);
					knotSpanWindow[0]=make_pair(knotVectorU[spanU],knotVectorV[spanV]);
					knotSpanWindow[1]=make_pair(knotVectorU[spanU+1],knotVectorV[spanV]);
					knotSpanWindow[2]=make_pair(knotVectorU[spanU+1],knotVectorV[spanV+1]);
					knotSpanWindow[3]=make_pair(knotVectorU[spanU],knotVectorV[spanV+1]);
					Polygon2D solution = c.clip(_polygonUV,knotSpanWindow);
					_listPolygon.push_back(solution);
					_listSpan.push_back(make_pair(spanU,spanV));
				}
			}
		}
	}
}

IGAMortarMapper::Polygon2D IGAMortarMapper::computeCanonicalElement(IGAPatchSurface* _thePatch, Polygon2D& _polygonUV, int _elementIndex) {
	// Subdivide the edges of the polygon in parametric space of Nurbs patch
	Polygon2D tmp;
	tmp.reserve(_polygonUV.size()*(numDivision));
	for(int node=0;node<_polygonUV.size();node++) {
		int nodeNext=(node+1)%_polygonUV.size();
		tmp.push_back(_polygonUV[node]);
		for(int sub=1;sub<numDivision;sub++){
			double t=(double)sub/(double)numDivision;
			double u=_polygonUV[node].first*(1.0-t)+_polygonUV[nodeNext].first*(t);
			double v=_polygonUV[node].second*(1.0-t)+_polygonUV[nodeNext].second*(t);
			tmp.push_back(make_pair(u,v));
		}
	}
	_polygonUV=tmp;
	int numNodesElementFE = meshFE->numNodesPerElem[_elementIndex];
	// Cartesian coordinates of the low order element
	double elementFEXYZ[12];
	for (int i = 0; i < numNodesElementFE; i++) {
		int nodeIndex = meshFEDirectElemTable[_elementIndex][i];
		for (int j = 0; j < 3; j++)
			elementFEXYZ[i * 3 + j] = meshFE->nodes[nodeIndex * 3 + j];
	}
	// Compute canonical coordinates polygon using parametric space
	Polygon2D polygonWZ;
	for(int i=0; i<_polygonUV.size(); i++) {
		double u = _polygonUV[i].first;
		double v = _polygonUV[i].second;
		double nodeXYZ[3];
		double normalVec[3];
		_thePatch->computeCartesianCoordinatesAndNormalVector(nodeXYZ, normalVec, u, v);
		// Compute intersection between normal of the nurbs surface and the low order element.
		// Must exist !
		double coordFE[2];
		if (numNodesElementFE == 3) {
			MathLibrary::computeIntersectionBetweenLineAndTriangle(
					elementFEXYZ, nodeXYZ, normalVec, coordFE);
			if(coordFE[0]<0) coordFE[0]=0;
			if(coordFE[1]<0) coordFE[1]=0;
			if(coordFE[0]>1) coordFE[0]=1;
			if(coordFE[1]>1) coordFE[1]=1;
		} else {
			MathLibrary::computeIntersectionBetweenLineAndQuad(
					elementFEXYZ, nodeXYZ, normalVec, coordFE);
			if(coordFE[0]<-1) coordFE[0]=-1;
			if(coordFE[1]<-1) coordFE[1]=-1;
			if(coordFE[0]>1) coordFE[0]=1;
			if(coordFE[1]>1) coordFE[1]=1;
		}
		double w = coordFE[0];
		double z = coordFE[1];
		polygonWZ.push_back(make_pair(w,z));
	}
	return polygonWZ;
}

void IGAMortarMapper::integrate(IGAPatchSurface* _thePatch, Polygon2D _polygonUV,
        int _spanU, int _spanV, Polygon2D _polygonWZ, int _elementIndex) {
    /*
     * 1. Divide the polygon into several quadratures(triangle or quadriliteral) for integration
     * 2. Loop through each quadrature
     *   2.1 Choose a Gauss quadrature (triangle or quadriliteral)
     *   2.2 Loop throught each Gauss point
     *       2.2.1 compute shape functions from Gauss points
     *       2.2.2 evaluate the coordinates in IGA patch from shape functions
     *       2.2.3 evaluate the coordinates in the linear element from shape functions
     *       2.2.4 compute the shape functions in the linear element for the current integration point
     *       2.2.5 Compute the local basis functions(shape functions of IGA) and their derivatives(for Jacobian)
     *       2.2.6 Compute the Jacobian from parameter space on IGA patch to physical
     *       2.2.7 Compute the Jacobian from the canonical space to the parameter space of IGA patch
     *       2.2.8 integrate the shape function product for C_NN(Linear shape function multiply linear shape function)
     *       2.2.9 integrate the shape function product for C_NR(Linear shape function multiply IGA shape function)
     * 3. Assemble the element coupling matrix to the global coupling matrix.
     */

    assert(!_polygonUV.empty());
    assert(!_polygonWZ.empty());

    int numNodesUV=_polygonUV.size();
    int numNodesWZ=_polygonWZ.size();
    assert(numNodesUV>2);
    assert(numNodesWZ>2);

    int _numNodesElementFE=meshFE->numNodesPerElem[_elementIndex];
	double elementFEXYZ[12];
	for (int i = 0; i < _numNodesElementFE; i++) {
		int nodeIndex = meshFEDirectElemTable[_elementIndex][i];
		for (int j = 0; j < 3; j++)
			elementFEXYZ[i * 3 + j] = meshFE->nodes[nodeIndex * 3 + j];
	}

    vector<double*> quadratureVecUV;
    vector<double*> quadratureVecWZ;
    vector<int> numNodesQuadrature;

    // Definitions
    double IGABasisFctsI = 0;
    double IGABasisFctsJ = 0;
    double basisFctsMaster = 0;
    double basisFctsSlave = 0;
    int numNodesElMaster = 0;
    int numNodesElSlave = 0;

    int pDegree = _thePatch->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree();
    int qDegree = _thePatch->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree();
    int nShapeFuncsIGA = (pDegree + 1) * (qDegree + 1);

    if (isMappingIGA2FEM) {
        numNodesElMaster = _numNodesElementFE;
        numNodesElSlave = nShapeFuncsIGA;
    } else {
        numNodesElMaster = nShapeFuncsIGA;
        numNodesElSlave = _numNodesElementFE;
    }

    double elementCouplingMatrixNN[numNodesElMaster * (numNodesElMaster + 1) / 2];
    double elementCouplingMatrixNR[numNodesElSlave * numNodesElMaster];

    for (int arrayIndex = 0; arrayIndex < numNodesElMaster * (numNodesElMaster + 1) / 2;
            arrayIndex++)
        elementCouplingMatrixNN[arrayIndex] = 0;

    for (int arrayIndex = 0; arrayIndex < numNodesElSlave * numNodesElMaster; arrayIndex++)
        elementCouplingMatrixNR[arrayIndex] = 0;

/// 1. Divide the polygon into several quadratures(triangle or quadriliteral) for integration
    if (numNodesUV <= 4) {
        double* nodesUV= new double[8];
        double* nodesWZ= new double[8];

        for (int i = 0; i < numNodesUV; i++) {
        	nodesUV[i*2]=_polygonUV[i].first;
        	nodesUV[i*2+1]=_polygonUV[i].second;
        	nodesWZ[i*2]=_polygonWZ[i].first;
        	nodesWZ[i*2+1]=_polygonWZ[i].second;
        }
        quadratureVecUV.push_back(nodesUV);
        quadratureVecWZ.push_back(nodesWZ);
        numNodesQuadrature.push_back(numNodesUV);
    } else {
        double centerIGA[2] = { 0, 0 };
        double centerFE[2] = { 0, 0 };
        for (int i = 0; i < numNodesUV; i++) {
            centerIGA[0] += _polygonUV[i].first  / numNodesUV;
            centerIGA[1] += _polygonUV[i].second / numNodesUV;
            centerFE[0]  += _polygonWZ[i].first  / numNodesUV;
            centerFE[1]  += _polygonWZ[i].second / numNodesUV;
        }

        double *triangleIGA;
        double *triangleFE;
        for (int i = 0; i < numNodesUV; i++) {
            int iNext = (i + 1) % numNodesUV;
            triangleIGA = new double[6]; // deleted
            triangleFE = new double[6]; // deleted
            triangleIGA[0] = _polygonUV[i].first;
            triangleIGA[1] = _polygonUV[i].second;
            triangleFE[0]  = _polygonWZ[i].first;
            triangleFE[1]  = _polygonWZ[i].second;

            triangleIGA[2] = _polygonUV[iNext].first;
            triangleIGA[3] = _polygonUV[iNext].second;
            triangleFE[2]  = _polygonWZ[iNext].first;
            triangleFE[3]  = _polygonWZ[iNext].second;


            triangleIGA[4] = centerIGA[0];
            triangleIGA[5] = centerIGA[1];
            triangleFE[4]  = centerFE[0];
            triangleFE[5]  = centerFE[1];

            quadratureVecUV.push_back(triangleIGA);
            quadratureVecWZ.push_back(triangleFE);
            numNodesQuadrature.push_back(3);
        }
    }

/// 2. Loop through each quadrature
    for (int quadratureCount = 0; quadratureCount < quadratureVecUV.size(); quadratureCount++) {

/// 2.1 Choose gauss triangle or gauss quadriliteral

        MathLibrary::IGAGaussQuadrature *theGaussQuadrature;

        int nNodesQuadrature;
        if (numNodesQuadrature[quadratureCount] == 3) {
            theGaussQuadrature = gaussTriangle;
            nNodesQuadrature = 3;
        } else {
            theGaussQuadrature = gaussQuad;
            nNodesQuadrature = 4;
        }

        double *quadratureUV = quadratureVecUV[quadratureCount];
        double *quadratureWZ = quadratureVecWZ[quadratureCount];

/// 2.2 Loop throught each Gauss point
        for (int GPCount = 0; GPCount < theGaussQuadrature->numGaussPoints; GPCount++) {

            /// 2.2.1 compute shape functions from Gauss points(in the quadrature).
            const double *GP = theGaussQuadrature->getGaussPoint(GPCount);

            double shapeFuncs[nNodesQuadrature];
            MathLibrary::computeLowOrderShapeFunc(nNodesQuadrature, GP, shapeFuncs);

            /// 2.2.2 evaluate the coordinates in IGA patch from shape functions
            double GPIGA[2];
            MathLibrary::computeLinearCombination(nNodesQuadrature, 2, quadratureUV, shapeFuncs,
                    GPIGA);


            /// 2.2.3 evaluate the coordinates in the linear element from shape functions
            double GPFE[2];
//    		double nodeXYZ[3];
//    		double normalVec[3];
//    		_thePatch->computeCartesianCoordinatesAndNormalVector(nodeXYZ, normalVec, GPIGA[0], GPIGA[1]);
//    		if (_numNodesElementFE == 3) {
//    			EMPIRE::MathLibrary::computeIntersectionBetweenLineAndTriangle(
//    					elementFEXYZ, nodeXYZ, normalVec, GPFE);
//    			if(coordFE[0]<0 || coordFE[0]>1) continue;
//        		if(coordFE[1]<0 || coordFE[1]>1) continue;
//    		} else {
//    			EMPIRE::MathLibrary::computeIntersectionBetweenLineAndQuad(
//    					elementFEXYZ, nodeXYZ, normalVec, GPFE);
//    			if(coordFE[0]<-1 || coordFE[0]>1) continue;
//        		if(coordFE[1]<-1 || coordFE[1]>1) continue;
//    		}
            MathLibrary::computeLinearCombination(nNodesQuadrature, 2, quadratureWZ, shapeFuncs,
                    GPFE);

            /// 2.2.4 compute the shape function(in the linear element) of the current integration point
            double shapeFuncsFE[_numNodesElementFE];
            MathLibrary::computeLowOrderShapeFunc(_numNodesElementFE, GPFE, shapeFuncsFE);
            int derivDegree = 1;

            /// 2.2.5 Compute the local basis functions(shape functions of IGA) and their derivatives(for Jacobian)
            double localBasisFunctionsAndDerivatives[(derivDegree + 1) * (derivDegree + 2)
                    * nShapeFuncsIGA / 2];

            _thePatch->getIGABasis()->computeLocalBasisFunctionsAndDerivatives(
                    localBasisFunctionsAndDerivatives, derivDegree, GPIGA[0], _spanU, GPIGA[1],
                    _spanV);

            /// 2.2.6 Compute the Jacobian from parameter space on IGA patch to physical
            double baseVectors[6];
            _thePatch->computeBaseVectors(baseVectors, localBasisFunctionsAndDerivatives, _spanU,
                    _spanV);

            double JacobianUVToPhysical = MathLibrary::computeAreaTriangle(baseVectors[0],
                    baseVectors[1], baseVectors[2], baseVectors[3], baseVectors[4], baseVectors[5])
                    * 2;

            /// 2.2.7 Compute the Jacobian from the canonical space to the parameter space of IGA patch
            double JacobianCanonicalToUV;
            if (nNodesQuadrature == 3) {
                JacobianCanonicalToUV = MathLibrary::computeAreaTriangle(
                        quadratureUV[2] - quadratureUV[0], quadratureUV[3] - quadratureUV[1], 0,
                        quadratureUV[4] - quadratureUV[0], quadratureUV[5] - quadratureUV[1], 0);
            } else {
                double dudx = .25
                        * (-(1 - GP[2]) * quadratureUV[0] + (1 - GP[2]) * quadratureUV[2]
                                + (1 + GP[2]) * quadratureUV[4] - (1 + GP[2]) * quadratureUV[6]);
                double dudy = .25
                        * (-(1 - GP[1]) * quadratureUV[0] - (1 + GP[1]) * quadratureUV[2]
                                + (1 + GP[1]) * quadratureUV[4] + (1 - GP[1]) * quadratureUV[6]);
                double dvdx = .25
                        * (-(1 - GP[2]) * quadratureUV[1] + (1 - GP[2]) * quadratureUV[3]
                                + (1 + GP[2]) * quadratureUV[5] - (1 + GP[2]) * quadratureUV[7]);
                double dvdy = .25
                        * (-(1 - GP[1]) * quadratureUV[1] - (1 + GP[1]) * quadratureUV[3]
                                + (1 + GP[1]) * quadratureUV[5] + (1 - GP[1]) * quadratureUV[7]);
                JacobianCanonicalToUV = fabs(dudx * dvdy - dudy * dvdx);
            }
            double Jacobian = JacobianUVToPhysical * JacobianCanonicalToUV;
            /// 2.2.8 integrate the shape function product for C_NN(Linear shape function multiply linear shape function)
            int count = 0;

            for (int i = 0; i < numNodesElMaster; i++)
                for (int j = i; j < numNodesElMaster; j++) {
                    if (isMappingIGA2FEM)
                        elementCouplingMatrixNN[count++] += shapeFuncsFE[i] * shapeFuncsFE[j]
                                * Jacobian * theGaussQuadrature->weights[GPCount];
                    else {
                        IGABasisFctsI =
                                localBasisFunctionsAndDerivatives[_thePatch->getIGABasis()->indexDerivativeBasisFunction(
                                        1, 0, 0, i)];
                        IGABasisFctsJ =
                                localBasisFunctionsAndDerivatives[_thePatch->getIGABasis()->indexDerivativeBasisFunction(
                                        1, 0, 0, j)];
                        elementCouplingMatrixNN[count++] += IGABasisFctsI * IGABasisFctsJ * Jacobian
                                * theGaussQuadrature->weights[GPCount];
                    }

                }

            /// 2.2.9 integrate the shape function product for C_NR(Linear shape function multiply IGA shape function)
            count = 0;

            for (int i = 0; i < numNodesElMaster; i++)
                for (int j = 0; j < numNodesElSlave; j++) {
                    if (isMappingIGA2FEM) {
                        basisFctsMaster = shapeFuncsFE[i];
                        basisFctsSlave =
                                localBasisFunctionsAndDerivatives[_thePatch->getIGABasis()->indexDerivativeBasisFunction(
                                        1, 0, 0, j)];
                    } else {
                        basisFctsMaster =
                                localBasisFunctionsAndDerivatives[_thePatch->getIGABasis()->indexDerivativeBasisFunction(
                                        1, 0, 0, i)];
                        basisFctsSlave = shapeFuncsFE[j];
                    }
                    elementCouplingMatrixNR[count++] += basisFctsMaster * basisFctsSlave * Jacobian
                            * theGaussQuadrature->weights[GPCount];
                }

        }
    }

/// 3.Assemble the element coupling matrix to the global coupling matrix.
    int dofIGA[nShapeFuncsIGA];
    _thePatch->getIGABasis()->getBasisFunctionsIndex(_spanU, _spanV, dofIGA);

    for (int i = 0; i < nShapeFuncsIGA; i++)
        dofIGA[i] = _thePatch->getControlPointNet()[dofIGA[i]]->getDofIndex();

    int count = 0;
    int dof1, dof2;
    for (int i = 0; i < numNodesElMaster; i++)
        for (int j = i; j < numNodesElMaster; j++) {
            if (isMappingIGA2FEM) {
                dof1 = meshFEDirectElemTable[_elementIndex][i];
                dof2 = meshFEDirectElemTable[_elementIndex][j];
            } else {
                dof1 = dofIGA[i];
                dof2 = dofIGA[j];
            }
            if (dof1 < dof2)
                (*C_NN)(dof1, dof2) += elementCouplingMatrixNN[count];
            else
                (*C_NN)(dof2, dof1) += elementCouplingMatrixNN[count];
            count++;
        }

    count = 0;
    for (int i = 0; i < numNodesElMaster; i++)
        for (int j = 0; j < numNodesElSlave; j++) {
            if (isMappingIGA2FEM) {
                dof1 = meshFEDirectElemTable[_elementIndex][i];
                dof2 = dofIGA[j];
            } else {
                dof1 = dofIGA[i];
                dof2 = meshFEDirectElemTable[_elementIndex][j];
            }
            (*C_NR)(dof1, dof2) += elementCouplingMatrixNR[count];
            count ++;
        }

    if (numNodesUV > 4)
        for (int i = 0; i < quadratureVecUV.size(); i++) {
            delete quadratureVecUV[i];
            delete quadratureVecWZ[i];
        }
}

bool IGAMortarMapper::computeKnotSpanOfProjElement(const IGAPatchSurface* _thePatch, const Polygon2D& _polygonUV, int* _span) {
    int minSpanU = _thePatch->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(
            _polygonUV[0].first);
    int minSpanV = _thePatch->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(
            _polygonUV[0].second);
    int maxSpanU = minSpanU;
    int maxSpanV = minSpanV;

    for (int nodeCount = 0; nodeCount < _polygonUV.size(); nodeCount++) {
        int spanU = _thePatch->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(_polygonUV[nodeCount].first);
        int spanV = _thePatch->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(_polygonUV[nodeCount].second);

        if (spanU < minSpanU)
            minSpanU = spanU;
        if (spanU > maxSpanU)
            maxSpanU = spanU;
        if (spanV < minSpanV)
            minSpanV = spanV;
        if (spanV > maxSpanV)
            maxSpanV = spanV;
    }

    bool OnSameKnotSpan=(minSpanU == maxSpanU && minSpanV == maxSpanV);
    if(_span!=NULL) {
		_span[0]=minSpanU;
		_span[1]=maxSpanU;
		_span[2]=minSpanV;
		_span[3]=maxSpanV;
    }
    return OnSameKnotSpan;
}

int IGAMortarMapper::getNeighbourElementofEdge(int _element, int _node1, int _node2) {
	for(int i=0; i<meshFE->numElems;i++) {
		bool isNode1=false;
		bool isNode2=false;
		for(int j=0;j<meshFE->numNodesPerElem[i];j++) {
			if(isNode1==false)
				isNode1=(meshFEDirectElemTable[i][j]==_node1)?true:false;
			if(isNode2==false)
				isNode2=(meshFEDirectElemTable[i][j]==_node2)?true:false;
			if(_element!=i && isNode1 && isNode2) {
				return i;
			}
		}
	}
	// If polygon is on boundary of mesh, can occur
	return -1;
}

bool IGAMortarMapper::computeIntermediatePoints(const int patchIndex, const int elemCount, const int nodeIndex1, const int nodeIndex2,
		const double* P1,const double* P2, Polygon2D& polygonUV, map<int,Polygon2D>& extraPolygonUV) {
	if(numDivision<=1)
		return 1;
	IGAPatchSurface* thePatch = meshIGA->getSurfacePatches()[patchIndex];
	// Get last point inside
	double u,v;
	if(polygonUV.empty()) {
		u = (*projectedCoords)[nodeIndex1][patchIndex][0];
		v = (*projectedCoords)[nodeIndex1][patchIndex][1];
	} else {
		u=polygonUV.back().first;
		v=polygonUV.back().second;
	}
	//Gets the neighbour element to the line
	int neighbourElement=getNeighbourElementofEdge(elemCount,nodeIndex1,nodeIndex2);
	//Init some flags
	bool isProjectedOnPatch=1, isProjectedOnPatchBoundary=1, outFlag=1;
	bool previousPointInside=0, currentPointInside=1;
	int MAX_DIV=numDivision;
	//Init data
	double div,dis;
	double P1P2[3];
	for(int i=0;i<3;i++)
		P1P2[i]=P2[i]-P1[i];
	//Loop in between the two point
	for(int sub=1;sub<MAX_DIV;sub++) {
		previousPointInside=currentPointInside;
		double t = (double)sub/(double)MAX_DIV;
		double P_projected[3];
		for(int i=0;i<3;i++)
			P_projected[i]=P1[i]+P1P2[i]*t;
		isProjectedOnPatch=thePatch->computePointProjectionOnPatch(u,v,P_projected);
		// If inside then compute boundary patch intersection between previous point and current point
		if(isProjectedOnPatch && neighbourElement!=-1) {
			currentPointInside=1;
			pair<double,double> uv_tmp1(u,v);
			pair<double,double> uv_tmp2;
			// Opposite flag means crossing the boundary, so we project on point on boundary
			if(currentPointInside != previousPointInside) {
				double P_inside[3], P_outside[3];
				for(int i=0;i<3;i++) {
					if(currentPointInside) {
						P_inside[i]= P1[i]+P1P2[i]*t;
						P_outside[i]=P1[i]+P1P2[i]*(t-1./MAX_DIV);
					} else {
						P_inside[i]= P1[i]+P1P2[i]*(t-1./MAX_DIV);
						P_outside[i]=P1[i]+P1P2[i]*t;
					}
				}
				isProjectedOnPatchBoundary=thePatch->computePointProjectionOnPatchBoundary_Bisection(u, v, div, dis, P_inside, P_outside);
				if (isProjectedOnPatchBoundary && dis <= disTol) {
					uv_tmp2=make_pair(u,v);
				}
			}
			if(currentPointInside && !previousPointInside) {
				polygonUV.push_back(uv_tmp2);
//				extraPolygonUV[neighbourElement].push_back(uv_tmp2);
			}
			polygonUV.push_back(uv_tmp1);
//			extraPolygonUV[neighbourElement].push_back(uv_tmp1);
			if(!currentPointInside && previousPointInside) {
				polygonUV.push_back(uv_tmp2);
//				extraPolygonUV[neighbourElement].push_back(uv_tmp2);
			}
		// If outside then compute patch intersection on boundary
		} else {
			currentPointInside=0;
		}
		outFlag = isProjectedOnPatchBoundary && outFlag;
	}
	return outFlag;
}


void IGAMortarMapper::consistentMapping(const double* _slaveField, double *_masterField) {
    /*
     * Mapping of the
     * C_NN * x_master = C_NR * x_slave
     */
    double* tmpVec = new double[numNodesMaster];

// 1. matrix vector product (x_tmp = C_NR * x_slave)
    C_NR->mulitplyVec(_slaveField, tmpVec, numNodesMaster);

// 2. solve C_NN * x_master = x_tmp
    C_NN->solve(_masterField, tmpVec);

    delete[] tmpVec;

}

void IGAMortarMapper::conservativeMapping(const double* _masterField, double *_slaveField) {
    /*
     * Mapping of the
     * f_slave = (C_NN^(-1) * C_NR)^T * f_master
     */

    double* tmpVec = new double[numNodesMaster];

// 1. solve C_NN * f_tmp = f_master;
    C_NN->solve(tmpVec, const_cast<double *>(_masterField));

// 2. matrix vector product (f_slave = C_NR^T * f_tmp)
    C_NR->transposeMulitplyVec(tmpVec, _slaveField, numNodesMaster);

    delete[] tmpVec;

}

void IGAMortarMapper::writeProjectedNodesOntoIGAMesh() {
    // Initializations
    IGAPatchSurface* IGAPatch;
    int numXiKnots, numEtaKnots, numXiControlPoints, numEtaControlPoints;
    double xiKnot, etaKnot;

    // Open file for writing the projected nodes
    ofstream projectedNodesFile;
    const string UNDERSCORE = "_";
    string projectedNodesFileName = "projectedNodesOntoNURBSSurface" + UNDERSCORE + name + ".m";
    projectedNodesFile.open(projectedNodesFileName.c_str());
    projectedNodesFile.precision(14);
    projectedNodesFile << std::dec;

    projectedNodesFile << HEADER_DECLARATION << endl << endl;

    // Loop over all the patches
    for (int patchCounter = 0; patchCounter < meshIGA->getSurfacePatches().size(); patchCounter++) {
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
            for (map<int, double*>::iterator it = (*projectedCoords)[nodeIndex].begin();
                    it != (*projectedCoords)[nodeIndex].end(); it++)
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


void IGAMortarMapper::debugPolygon(const Polygon2D& _polygon, string _name) {
	DEBUG_OUT()<<"----------------------------------"<<endl;
	if(_name!="")
		DEBUG_OUT()<<"Polygon name : "<<_name<<endl;
	for(int i=0; i<_polygon.size();i++)
		DEBUG_OUT()<<"\t"<<"u="<<_polygon[i].first<<" / v="<<_polygon[i].second<<endl;
	DEBUG_OUT()<<"----------------------------------"<<endl;
}
void IGAMortarMapper::debugPolygon(const ListPolygon2D& _listPolygon, string _name) {
	DEBUG_OUT()<<"++++++++++++++++++++++++++"<<endl;
	if(_name!="")
		DEBUG_OUT()<<"Polygon list name : "<<_name<<endl;
	for(int i=0; i<_listPolygon.size();i++) {
		DEBUG_OUT()<<"Polygon index : "<<i<<endl;
		debugPolygon(_listPolygon[i]);
	}
	DEBUG_OUT()<<"++++++++++++++++++++++++++"<<endl;
}

void IGAMortarMapper::printCouplingMatrices() {

    ERROR_OUT() << "C_NN" << endl;
    C_NN->printCSR();
    ERROR_OUT() << "C_NR" << endl;
    C_NR->printCSR();
}

void IGAMortarMapper::printCouplingMatricesToFile() {
	DEBUG_OUT()<<"Size of C_NR is "<<numNodesMaster<<" by "<<numNodesSlave<<endl;
    if(Message::userSetOutputLevel==Message::DEBUG) {
		C_NR->printToFile("C_NR.dat");
		C_NN->printToFile("C_NN.dat");
	}
}

void IGAMortarMapper::checkConsistency() {
    double ones[numNodesSlave];
    for(int i=0;i<numNodesSlave;i++) {
    	ones[i]=1.0;
    }
    double output[numNodesMaster];
    this->consistentMapping(ones,output);
    double norm=0;
    for(int i=0;i<numNodesMaster;i++) {
    	norm+=output[i]*output[i];
    }
    norm=sqrt(norm/numNodesMaster);
    DEBUG_OUT()<<"Norm of output field = "<<norm<<endl;
    if(fabs(norm-1.0)>1e-6) {
    	ERROR_OUT()<<"Coupling not consistent !"<<endl;
    	ERROR_OUT()<<"Coupling of unit field deviating from 1 of "<<fabs(norm-1.0)<<endl;
    	exit(-1);
    }
}

}


