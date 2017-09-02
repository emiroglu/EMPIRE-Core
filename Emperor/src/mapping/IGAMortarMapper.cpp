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

#include "IGAMortarMapper.h"
#include "IGAPatchSurface.h"
#include "WeakIGADirichletCurveCondition.h"
#include "WeakIGADirichletSurfaceCondition.h"
#include "WeakIGAPatchContinuityCondition.h"
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
#include "IGAMortarCouplingMatrices.h"
#include <iomanip>

using namespace std;

namespace EMPIRE {

/// Declaration statement
static const string HEADER_DECLARATION = "Author: Andreas Apostolatos";

IGAMortarMapper::IGAMortarMapper(std::string _name, IGAMesh *_meshIGA, FEMesh *_meshFE,
                                 bool _isWeakConditions, bool _isMappingIGA2FEM) :
    name(_name), meshIGA(_meshIGA), isMappingIGA2FEM(_isMappingIGA2FEM) {

    // Check input
    assert(_meshIGA != NULL);
    assert(_meshFE != NULL);
    assert(_meshIGA->type == EMPIRE_Mesh_IGAMesh);
    assert(_meshFE->type == EMPIRE_Mesh_FEMesh);

    // Number of Cartesian directions
    int noCoord = 3;

    // Check if the FE mesh is triangulated and if not triangulate it
    if (_meshFE->triangulate() == NULL)
        meshFE = _meshFE;
    else
        meshFE = _meshFE->triangulate();

    // Assign the mapper type
    mapperType = EMPIRE_IGAMortarMapper;

    projectedCoords.resize(meshFE->numNodes);
    projectedPolygons.resize(meshFE->numElems);
    triangulatedProjectedPolygons.resize(meshFE->numElems);

    // Find the mapping direction (meshIGA -> meshFEM or meshFEM -> meshIGA)
    if (isMappingIGA2FEM) {
        numNodesSlave = meshIGA->getNumNodes();
        numNodesMaster = meshFE->numNodes;
    } else {
        numNodesSlave = meshFE->numNodes;
        numNodesMaster = meshIGA->getNumNodes();
    }

    // Flag on whether the expanded version of the coupling matrices is assumed
    isExpanded = _isWeakConditions && ~isMappingIGA2FEM;

    // Initialize flag on whether the meshFEDirectElemTable was created
    isMeshFEDirectElemTable = false;

    // Initialize flag on whether the Gauss quadrature has been defined
    isGaussQuadrature = false;

    // Get the number of weak Dirichlet curve conditions
    noWeakIGADirichletCurveConditions = meshIGA->getWeakIGADirichletCurveConditions().size();

    // Initialize the penalty factors for the application of weak Dirichlet curve conditions
    weakDirichletCCAlphaPrimary = new double[noWeakIGADirichletCurveConditions];
    weakDirichletCCAlphaSecondaryBending = new double[noWeakIGADirichletCurveConditions];
    weakDirichletCCAlphaSecondaryTwisting = new double[noWeakIGADirichletCurveConditions];

    // Get the number of weak Dirichlet surface conditions
    noWeakIGADirichletSurfaceConditions = meshIGA->getWeakIGADirichletSurfaceConditions().size();

    // Initialize the penalty factors for the application of weak Dirichlet surface conditions
    weakDirichletSCAlphaPrimary = new double[noWeakIGADirichletSurfaceConditions];
    weakDirichletSCAlphaSecondaryBending = new double[noWeakIGADirichletSurfaceConditions];
    weakDirichletSCAlphaSecondaryTwisting = new double[noWeakIGADirichletSurfaceConditions];

    // Get the number of weak continuity conditions
    noWeakIGAPatchContinuityConditions = meshIGA->getWeakIGAPatchContinuityConditions().size();

    // Initialize the penalty factors for the application of weak patch continuity conditions
    weakPatchContinuityAlphaPrimaryIJ = new double[noWeakIGAPatchContinuityConditions];
    weakPatchContinuityAlphaSecondaryBendingIJ = new double[noWeakIGAPatchContinuityConditions];
    weakPatchContinuityAlphaSecondaryTwistingIJ = new double[noWeakIGAPatchContinuityConditions];

    // Initialize the problem parameters with default values
    setParametersConsistency();
    setParametersProjection();
    setParametersNewtonRaphson();
    setParametersNewtonRaphsonBoundary();
    setParametersBisection();
    setParametersIntegration();
    setParametersWeakCurveDirichletConditions();
    setParametersWeakSurfaceDirichletConditions();
    setParametersWeakPatchContinuityConditions();
    setParametersErrorComputation();

    // Initialize coupling matrices taking into consideration whether the expanded coupling matrices are assumed or not
    int size_N;
    int size_R;
    if (isExpanded){
        size_N = noCoord*numNodesMaster;
        size_R = noCoord*numNodesSlave;
    }else {
        size_N = numNodesMaster;
        size_R = numNodesSlave;
    }
    couplingMatrices = new IGAMortarCouplingMatrices(size_N, size_R, isExpanded);
}

void IGAMortarMapper::setParametersConsistency(bool _enforceConsistency, double _tolConsistency) {
    propConsistency.enforceConsistency = _enforceConsistency;
    propConsistency.tolConsistency = _tolConsistency;
}

void IGAMortarMapper::setParametersProjection(double _maxProjectionDistance, int _noInitialGuess,
                                              double _maxProjectionDistanceOnDifferentPatches) {
    propProjection.maxProjectionDistance = _maxProjectionDistance;
    propProjection.noInitialGuess = _noInitialGuess;
    propProjection.maxProjectionDistanceOnDifferentPatches = _maxProjectionDistanceOnDifferentPatches;
}

void IGAMortarMapper::setParametersNewtonRaphson(int _noIterations, double _tolProjection) {
    propNewtonRaphson.noIterations = _noIterations;
    propNewtonRaphson.tolProjection = _tolProjection;
}

void IGAMortarMapper::setParametersNewtonRaphsonBoundary(int _noIterations, double _tolProjection) {
    propNewtonRaphsonBoundary.noIterations = _noIterations;
    propNewtonRaphsonBoundary.tolProjection = _tolProjection;
}

void IGAMortarMapper::setParametersBisection(int _noIterations, double _tolProjection) {
    propBisection.noIterations = _noIterations;
    propBisection.tolProjection = _tolProjection;
}

void IGAMortarMapper::setParametersIntegration(int _noGPTriangle, int _noGPQuad) {
    propIntegration.noGPTriangle = _noGPTriangle;
    propIntegration.noGPQuad = _noGPQuad;
}

void IGAMortarMapper::setParametersWeakCurveDirichletConditions(bool _isWeakCurveDirichletConditions, bool _isAutomaticPenaltyParameters,
                                                                double _alphaPrim, double _alphaSecBending, double _alphaSecTwisting) {
    propWeakCurveDirichletConditions.isWeakCurveDirichletConditions = _isWeakCurveDirichletConditions;
    propWeakCurveDirichletConditions.isAutomaticPenaltyParameters = _isAutomaticPenaltyParameters;
    propWeakCurveDirichletConditions.alphaPrim = _alphaPrim;
    propWeakCurveDirichletConditions.alphaSecBending = _alphaSecBending;
    propWeakCurveDirichletConditions.alphaSecTwisting = _alphaSecTwisting;
}

void IGAMortarMapper::setParametersWeakSurfaceDirichletConditions(bool _isWeakSurfaceDirichletConditions, bool _isAutomaticPenaltyParameters,
                                                                  double _alphaPrim) {
    propWeakSurfaceDirichletConditions.isWeakSurfaceDirichletConditions = _isWeakSurfaceDirichletConditions;
    propWeakSurfaceDirichletConditions.isAutomaticPenaltyParameters = _isAutomaticPenaltyParameters;
    propWeakSurfaceDirichletConditions.alphaPrim = _alphaPrim;
}

void IGAMortarMapper::setParametersWeakPatchContinuityConditions(bool _isWeakPatchContinuityConditions, bool _isAutomaticPenaltyParameters,
                                                                 double _alphaPrim, double _alphaSecBending, double _alphaSecTwisting) {
    propWeakPatchContinuityConditions.isWeakPatchContinuityConditions = _isWeakPatchContinuityConditions;
    propWeakPatchContinuityConditions.isAutomaticPenaltyParameters = _isAutomaticPenaltyParameters;
    propWeakPatchContinuityConditions.alphaPrim = _alphaPrim;
    propWeakPatchContinuityConditions.alphaSecBending = _alphaSecBending;
    propWeakPatchContinuityConditions.alphaSecTwisting = _alphaSecTwisting;
}

void IGAMortarMapper::setParametersStrongCurveDirichletConditions(bool _isStrongCurveDirichletConditions) {
    propStrongCurveDirichletConditions.isStrongCurveDirichletConditions = _isStrongCurveDirichletConditions;
}

void IGAMortarMapper::setParametersErrorComputation(bool _isErrorComputation, bool _isDomainError, bool _isCurveError, bool _isInterfaceError){
    propErrorComputation.isErrorComputation = _isErrorComputation;
    propErrorComputation.isDomainError = _isDomainError;
    propErrorComputation.isCurveError = _isCurveError;
    propErrorComputation.isInterfaceError = _isInterfaceError;
}

void IGAMortarMapper::buildCouplingMatrices() {
    HEADING_OUT(3, "IGAMortarMapper", "Building coupling matrices for ("+ name +")...", infoOut);
    {
        int nIG = meshIGA->getNumNodes();
        int nFE = meshFE->numNodes;
        INFO_OUT() << "Number of nodes in NURBS mesh is " << nIG << endl;
        INFO_OUT() << "Number of nodes in FE mesh is    " << nFE << endl;
        INFO_OUT() << "Size of matrices will be " << (isMappingIGA2FEM?nFE:nIG) << "x" << (isMappingIGA2FEM?nFE:nIG) << " and "
                   << (isMappingIGA2FEM?nFE:nIG) << "x" << (isMappingIGA2FEM?nIG:nFE) << endl;
    }

    //Instantiate quadrature rules
    isGaussQuadrature = true;
    gaussTriangle = new MathLibrary::IGAGaussQuadratureOnTriangle(propIntegration.noGPTriangle);
    gaussQuad = new MathLibrary::IGAGaussQuadratureOnQuad(propIntegration.noGPQuad);

    // Initialize a string holding the filenames
    string filename;

    // Set default scheme values
    IGAPatchSurface::MAX_NUM_ITERATIONS = propNewtonRaphson.noIterations;
    IGAPatchSurface::TOL_ORTHOGONALITY = propNewtonRaphson.tolProjection;

    // Compute the EFT for the FE mesh
    initTables();
    isMeshFEDirectElemTable = true;

    // Project the FE nodes onto the multipatch trimmed geometry
    projectPointsToSurface();

    // Write the projected points on to a file only in DEBUG mode to be used in MATLAB
    if (Message::isDebugMode())
        writeProjectedNodesOntoIGAMesh();

    // Reserve some space for gauss point values in the domain
    if (propErrorComputation.isDomainError)
        streamGPs.reserve(8*meshFE->numElems*gaussQuad->numGaussPoints);

    // Reserve some space for the gauss point values along each trimming curve where conditions are applied
    if(propErrorComputation.isCurveError){
        int noCurveGPs = 0;
        std::vector<WeakIGADirichletCurveCondition*> weakIGADirichletCurveConditions = meshIGA->getWeakIGADirichletCurveConditions();
        for (int iWDC = 0; iWDC < weakIGADirichletCurveConditions.size(); iWDC++){
            noCurveGPs += weakIGADirichletCurveConditions[iWDC]->getCurveNumGP();
        }
        streamInterfaceGPs.reserve(noCurveGPs);
    }
    // Reserve some space for the interface gauss point values
    if(propErrorComputation.isInterfaceError){
        int noInterfaceGPs = 0;
        std::vector<WeakIGAPatchContinuityCondition*> weakIGAPatchContinuityConditions = meshIGA->getWeakIGAPatchContinuityConditions();
        for (int iWCC = 0; iWCC < weakIGAPatchContinuityConditions.size(); iWCC++){
            noInterfaceGPs += weakIGAPatchContinuityConditions[iWCC]->getTrCurveNumGP();
        }
        streamInterfaceGPs.reserve(noInterfaceGPs);
    }

    // Compute CNN and CNR
    computeCouplingMatrices();

    // Write the gauss point data into an csl file
    if(Message::isDebugMode())
        writeGaussPointData();

    // Write polygon net of projected elements to a vtk file
    writeCartesianProjectedPolygon("trimmedPolygonsOntoNURBSSurface", trimmedProjectedPolygons);
    writeCartesianProjectedPolygon("integratedPolygonsOntoNURBSSurface", triangulatedProjectedPolygons2);
    trimmedProjectedPolygons.clear();
    triangulatedProjectedPolygons2.clear();

    // On the application of strong Dirichlet boundary conditions
//    bool _isdirichletBCs = false;
//    bool isClampedDofs = false;
//    std::vector<int> clampedIDs;
//    if(dirichletBCs.isDirichletBCs==1) {
//        _isdirichletBCs = true;
//        clampedIDs = meshIGA->getClampedDofs();

//        int clampedDirections = meshIGA->getClampedDirections();
//        if(clampedDirections == 1 || clampedDirections == 2)
//            isClampedDofs = true;
//    }

    // Compute the Penalty parameters for the application of weak Dirichlet curve conditions
    if(propWeakCurveDirichletConditions.isWeakCurveDirichletConditions && !isMappingIGA2FEM) {
        filename = name + "_penaltyParametersWeakDirichletConditions.txt";
        computePenaltyParametersForWeakDirichletCurveConditions(filename);
    }

    // Compute the Penalty parameters for the application of weak Dirichlet surface conditions
    if(propWeakSurfaceDirichletConditions.isWeakSurfaceDirichletConditions && !isMappingIGA2FEM)
        computePenaltyParametersForWeakDirichletSurfaceConditions();

    // Compute the Penalty parameters for the application of weak patch continuity conditions
    if(propWeakPatchContinuityConditions.isWeakPatchContinuityConditions && !isMappingIGA2FEM) {
        filename = name + "_penaltyParametersWeakContinuityConditions.txt";
        computePenaltyParametersForPatchContinuityConditions(filename);
    }

    // Compute the Penalty matrices for the application of weak continuity conditions between the multipatches
    if (propWeakPatchContinuityConditions.isWeakPatchContinuityConditions && !isMappingIGA2FEM) {
        INFO_OUT() << "Application of weak patch continuity conditions started" << endl;
        if(!propWeakPatchContinuityConditions.isAutomaticPenaltyParameters) {
            INFO_OUT() << "Manual patch coupling penalties: alphaPrim = "<< propWeakPatchContinuityConditions.alphaPrim
                       << " alphaSecBending = " <<  propWeakPatchContinuityConditions.alphaSecBending <<" alphaSecTwisting = "
                       <<  propWeakPatchContinuityConditions.alphaSecTwisting << endl;
        } else {
            INFO_OUT() << "Automatic patch coupling penalties" << endl;
        }
        computeIGAPatchWeakContinuityConditionMatrices();
        INFO_OUT() << "Application of weak patch continuity conditions finished" << std::endl;
    } else
        INFO_OUT() << "No application of weak patch continuity conditions is assumed" << std::endl;

    // Remove empty rows and columns from system in case consistent mapping for the traction from FE Mesh to IGA multipatch geometry is required
    if(!isMappingIGA2FEM){
        INFO_OUT() << "Enforcing flying nodes in Cnn" << std::endl;
        couplingMatrices->enforceCnn();
    }

    // Write the coupling matrices in files
    if(Message::isDebugMode())
        writeCouplingMatricesToFile();

    // Check and enforce consistency in Cnn. This has to be done before the application of the Dirichlet conditions in general.
    if (propConsistency.enforceConsistency)
        enforceConsistency();

    // Compute the Penalty matrices for the application of weak Dirichlet conditions along trimming curves
    if (propWeakCurveDirichletConditions.isWeakCurveDirichletConditions) {
        INFO_OUT() << "Application of weak Dirichlet curve conditions started" << endl;
        if(!propWeakCurveDirichletConditions.isAutomaticPenaltyParameters) {
            INFO_OUT() << "Manual weak Dirichlet curve condition penalties: alphaPrim = " << propWeakCurveDirichletConditions.alphaPrim << " alphaSecBending = " <<  propWeakCurveDirichletConditions.alphaSecBending << " alphaSecTwisting = " <<  propWeakCurveDirichletConditions.alphaSecTwisting << endl;
        } else {
            INFO_OUT() << "Automatic weak Dirichlet curve condition penalties, use DEBUG mode to see the computed values" << endl;
        }
        computeIGAWeakDirichletCurveConditionMatrices();
        INFO_OUT() << "Application of weak Dirichlet curve conditions finished" << std::endl;
    } else
        INFO_OUT() << "No application of weak Dirichlet curve conditions are assumed" << std::endl;

    // Compute the Penalty matrices for the application of weak Dirichlet conditions across surfaces
    if (propWeakSurfaceDirichletConditions.isWeakSurfaceDirichletConditions) {
        ERROR_OUT() << "Function under construction" << endl;
        exit(-1);
        INFO_OUT() << "Application of weak Dirichlet surface conditions started" << endl;
        if(!propWeakSurfaceDirichletConditions.isAutomaticPenaltyParameters) {
            INFO_OUT() << "Manual weak Dirichlet surface condition penalties: alphaPrim = "<< propWeakSurfaceDirichletConditions.alphaPrim << " alphaSecBending = " << endl;
        } else {
            INFO_OUT() << "Automatic weak Dirichlet surface condition penalties, use DEBUG mode to see the computed values" << endl;
        }
        computeIGAWeakDirichletSurfaceConditionMatrices();
        INFO_OUT() << "Application of weak Dirichlet surface conditions finished" << std::endl;
    } else
        INFO_OUT() << "No application of weak Dirichlet surface conditions are assumed" << std::endl;

    // this is for MapperAdapter, this makes sure that the fields in all directions are sent at the same time
    // if it is not always clamped in all three dircetions
//    if(isClampedDofs)
//        isIGAPatchContinuityConditions = true;

//    couplingMatrices->setIsDirichletBCs(_isdirichletBCs);

//    if(dirichletBCs.isDirichletBCs == 1){
//        INFO_OUT() << "Applying strong Dirichlet boundary conditions" << std::endl;
//        couplingMatrices->applyDirichletBCs(clampedIDs);
//    }

//    if (isIGAWeakDirichletCurveConditions || dirichletBCs.isDirichletBCs == 1) {
//        couplingMatrices->factorizeCnn();
//        INFO_OUT() << "Factorize was successful" << std::endl;
//    }
}

IGAMortarMapper::~IGAMortarMapper() {

    // Delete the Finite Element direct element table
    if(isMeshFEDirectElemTable){
        for (int i = 0; i < meshFE->numElems; i++)
            delete[] meshFEDirectElemTable[i];
        delete[] meshFEDirectElemTable;
    }

    // Delete the quadrature rules
    if(isGaussQuadrature){
        delete gaussTriangle;
        delete gaussQuad;
    }

    // Delete the coupling matrices
    delete couplingMatrices;

    // Delete the penalty parameters for the application of weak Dirichlet conditions along trimming curves
    delete[] weakDirichletCCAlphaPrimary;
    delete[] weakDirichletCCAlphaSecondaryBending;
    delete[] weakDirichletCCAlphaSecondaryTwisting;

    // Delete the penalty parameters for the continuity enforcement across the patch interfaces
    delete[] weakPatchContinuityAlphaPrimaryIJ;
    delete[] weakPatchContinuityAlphaSecondaryBendingIJ;
    delete[] weakPatchContinuityAlphaSecondaryTwistingIJ;

    // Delete the stored GP values
    if (propErrorComputation.isDomainError)
        streamGPs.clear();
    if (propErrorComputation.isCurveError)
        streamCurveGPs.clear();
    if (propErrorComputation.isInterfaceError)
        streamInterfaceGPs.clear();
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

void IGAMortarMapper::projectPointsToSurface() {
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

    // Get the number of patches in the IGA mesh
    int numPatches = meshIGA->getNumPatches();

    // Bounding box preprocessing, assign to each node the patches to be visited
    INFO_OUT() << "Bounding box preprocessing started" << endl;
    time(&timeStart);
    for (int i = 0; i < meshFE->numNodes; i++) {
        double P[3];
        P[0] = meshFE->nodes[3 * i + 0];
        P[1] = meshFE->nodes[3 * i + 1];
        P[2] = meshFE->nodes[3 * i + 2];
        for(int patchCount = 0; patchCount < numPatches; patchCount++) {
            IGAPatchSurface* thePatch = meshIGA->getSurfacePatch(patchCount);
            bool isInside = thePatch->getBoundingBox().isPointInside(P, propProjection.maxProjectionDistance);
            if(isInside)
                patchToProcessPerNode[i].insert(patchCount);
        }
        if(patchToProcessPerNode[i].empty()) {
            stringstream msg;
            msg << "Node [" << i << "] is not in any bounding box of NURBS patches ! Increase maxProjectionDistance !";
            ERROR_BLOCK_OUT("IGAMortarMapper", "projectPointsToSurface", msg.str());
        }
    }
    time(&timeEnd);
    INFO_OUT()<<"Bounding box preprocessing done in "<< difftime(timeEnd, timeStart) << " seconds"<<endl;
    // Project the node for every patch's bounding box the node lies into
    // or on every patch if not found in a single bounding box
    INFO_OUT()<<"First pass projection started"<<endl;
    time(&timeStart);
    for(int i = 0; i < meshFE->numElems; i++) {
        int numNodesInElem = meshFE->numNodesPerElem[i];
        for(int patchIndex = 0; patchIndex < numPatches; patchIndex++) {
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
    INFO_OUT()<<"First pass projection done in "<< difftime(timeEnd, timeStart) << " seconds"<<endl;
    int missing = 0;
    for (int i = 0; i < meshFE->numNodes; i++) {
        if(!isProjected[i]) {
            missing++;
            WARNING_OUT()<<"Node not projected at first pass ["<<i<<"] of coordinates "<<meshFE->nodes[3*i]<<","<<meshFE->nodes[3*i+1]<<","<<meshFE->nodes[3*i+2]<<endl;
        }
    }
    INFO_OUT()<<meshFE->numNodes - missing << " nodes over " << meshFE->numNodes <<" could be projected during first pass" << endl;
    double initialTolerance = propNewtonRaphson.tolProjection;

    // Second pass projection --> relax Newton-Rapshon tolerance and if still fails refine the sampling points for the Newton-Raphson initial guesses
    if(missing) {
        INFO_OUT()<<"Second pass projection started"<<endl;
        time(&timeStart);
        missing = 0;
        for (int i = 0; i < meshFE->numNodes; i++) {
            if(!isProjected[i]) {
                propNewtonRaphson.tolProjection = 10*propNewtonRaphson.tolProjection;
                for(set<int>::iterator patchIndex = patchToProcessPerNode[i].begin();patchIndex != patchToProcessPerNode[i].end(); patchIndex++) {
                    computeInitialGuessForProjection(*patchIndex, meshFENodeToElementTable[i][0], i, initialU, initialV);
                    bool flagProjected = projectPointOnPatch(*patchIndex, i, initialU, initialV, minProjectionDistance[i], minProjectionPoint[i]);
                    isProjected[i] = isProjected[i] || flagProjected;
                }
                if(!isProjected[i]) {
                    for(set<int>::iterator patchIndex = patchToProcessPerNode[i].begin(); patchIndex != patchToProcessPerNode[i].end(); patchIndex++) {
                        bool flagProjected = forceProjectPointOnPatch(*patchIndex, i, minProjectionDistance[i], minProjectionPoint[i]);
                        isProjected[i] = isProjected[i] || flagProjected;
                    }
                }
            }
            if(!isProjected[i]) {
                ERROR_OUT()<<"Node not projected at second pass ["<<i<<"] of coordinates "<<meshFE->nodes[3*i]<<","<<meshFE->nodes[3*i+1]<<","<<meshFE->nodes[3*i+2]<<endl;
                missing++;
            }
            propNewtonRaphson.tolProjection = initialTolerance;
        }
        propNewtonRaphson.tolProjection = initialTolerance;
        time(&timeEnd);
        INFO_OUT()<<"Second pass projection done! It took "<< difftime(timeEnd, timeStart) << " seconds"<<endl;
        if(missing) {
            stringstream msg;
            msg << missing << " nodes over " << meshFE->numNodes << " could NOT be projected during second pass!" << endl;
            msg << "Treatment possibility 1." << endl;
            msg << "Possibly relax parameters in projectionProperties or newtonRaphson" << endl;
            msg << "Treatment possibility 2." << endl;
            msg << "Remesh with higher accuracy on coordinates of the FE nodes, i.e. more digits" << endl;
            ERROR_BLOCK_OUT("IGAMortarMapper", "ProjectPointsToSurface", msg.str());
        }
    }
}

void IGAMortarMapper::computeInitialGuessForProjection(const int _patchIndex, const int _elemIndex, const int _nodeIndex, double& _u, double& _v) {

    IGAPatchSurface* thePatch = meshIGA->getSurfacePatch(_patchIndex);
    /// 1iii.1. Initialize the flag to false and the node id to zero
    bool isNodeInsideElementProjected = false;
    int projectedNode = -1;
    /// 1iii.2. Loop over all nodes of the current element to check if there exist one node has been projected already
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
    /// 1iii.3. Check if there exist one node in the current element has been successfully projected
    if (isNodeInsideElementProjected) {
        /// 1iii.3i. If so, use result of the projected node as the initial guess for the projection step
        _u = projectedCoords[projectedNode][_patchIndex][0];
        _v = projectedCoords[projectedNode][_patchIndex][1];
    } else {
        /// 1iii.3ii. Otherwise, find the nearest knot intersection as initial guess for the projection step
        // Get the Cartesian coordinates of that input node
        double P[3];
        P[0] = meshFE->nodes[_nodeIndex * 3 + 0];
        P[1] = meshFE->nodes[_nodeIndex * 3 + 1];
        P[2] = meshFE->nodes[_nodeIndex * 3 + 2];
        // Get accordingly an initial guess for the projection onto the NURBS patch
        thePatch->findInitialGuess4PointProjection(_u, _v, P, propProjection.noInitialGuess, propProjection.noInitialGuess);
    }
}

bool IGAMortarMapper::projectPointOnPatch(const int patchIndex, const int nodeIndex, const double u0, const double v0, double& minProjectionDistance, vector<double>& minProjectionPoint) {

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
                                                                hasResidualConverged, propNewtonRaphson.noIterations, propNewtonRaphson.tolProjection);
    double distance = MathLibrary::computePointDistance(P, projectedP);
    if(hasConverged &&  distance < propProjection.maxProjectionDistance) {
        /// Perform some validity checks to validate the projected point
        if(distance > minProjectionDistance + propProjection.maxProjectionDistanceOnDifferentPatches) {
            return false;
        }
        if(!minProjectionPoint.empty() &&
                MathLibrary::computePointDistance(projectedP, &minProjectionPoint[0]) > propProjection.maxProjectionDistanceOnDifferentPatches &&
                distance > minProjectionDistance) {
            return false;
        }
        if(distance < minProjectionDistance - propProjection.maxProjectionDistanceOnDifferentPatches
                || MathLibrary::computePointDistance(projectedP, &minProjectionPoint[0]) > propProjection.maxProjectionDistanceOnDifferentPatches) {
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

bool IGAMortarMapper::forceProjectPointOnPatch(const int patchIndex, const int nodeIndex, double& minProjectionDistance, vector<double>& minProjectionPoint) {
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
    thePatch->findInitialGuess4PointProjection(u, v, P, 200, 200);
    double uv[2] = {u, v};
    thePatch->computeCartesianCoordinates(projectedP, uv);
    double distance = MathLibrary::computePointDistance(P, projectedP);
    /// Perform some validity checks
    if(distance > minProjectionDistance + propProjection.maxProjectionDistanceOnDifferentPatches) {
        return false;
    }
    if(distance < minProjectionDistance - propProjection.maxProjectionDistanceOnDifferentPatches) {
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
    // Time stamps
    time_t timeStart, timeEnd;

    /// List of integrated element
    set<int> elementIntegrated;

    /// The vertices of the canonical polygons
    double parentTriangle[6] = { 0, 0, 1, 0, 0, 1 };
    double parentQuadriliteral[8] = { -1, -1, 1, -1, 1, 1, -1, 1 };

    int elementStringLength;
    {
        stringstream ss;
        ss << meshFE->numElems;
        elementStringLength = ss.str().length();
    }

    INFO_OUT() << "Computing coupling matrices started" << endl;
    time(&timeStart);
    /// Loop over all the elements in the FE side
    for (int elemIndex = 0; elemIndex < meshFE->numElems; elemIndex++) {
        DEBUG_OUT()<< setfill ('#') << setw(18+elementStringLength) << "#" << endl;
        DEBUG_OUT()<< setfill (' ') << "### ELEMENT ["<< setw(elementStringLength) << elemIndex << "] ###"<<endl;
        DEBUG_OUT()<< setfill ('#') << setw(18+elementStringLength) << "#" << setfill (' ')<< endl;
        // Get the number of shape functions. Depending on number of nodes in the current element
        int numNodesElementFE = meshFE->numNodesPerElem[elemIndex];
        /// Find whether the projected FE element is located on one patch
        set<int> patchWithFullElt;
        set<int> patchWithSplitElt;
        getPatchesIndexElementIsOn(elemIndex, patchWithFullElt, patchWithSplitElt);
        DEBUG_OUT()<<"Element FULLY projected on \t" << patchWithFullElt.size() << " patch" << endl;
        DEBUG_OUT()<<"Element PARTLY projected on \t" << patchWithSplitElt.size() << " patch" << endl;
        /////////////////////////////////////
        /// Compute the coupling matrices ///
        /////////////////////////////////////
        /// 1. If the current element can be projected on one patch
        for (set<int>::iterator it = patchWithFullElt.begin();
             it != patchWithFullElt.end(); ++it) {
            int patchIndex=*it;
            /// Get the projected coordinates for the current element
            Polygon2D polygonUV;
            bool isProjectedOnPatchBoundary = true;
            /// 1.1 Get initial polygon from projection
            // For every point of polygon
            buildFullParametricElement(elemIndex, numNodesElementFE, patchIndex, polygonUV);
            ClipperAdapter::cleanPolygon(polygonUV);
            bool isIntegrated = computeLocalCouplingMatrix(elemIndex, patchIndex, polygonUV);
            if(isIntegrated) {
                elementIntegrated.insert(elemIndex);
                projectedPolygons[elemIndex][patchIndex]=polygonUV;
            }
        }
        /// 2. If the current element is split in more than one patches
        // Loop over all the patches in the IGA Mesh having a part of the FE element projected inside
        for (set<int>::iterator it = patchWithSplitElt.begin(); it != patchWithSplitElt.end(); it++) {
            int patchIndex=*it;
            IGAPatchSurface* thePatch = meshIGA->getSurfacePatch(patchIndex);
            // Stores points of the polygon clipped by the nurbs patch
            Polygon2D polygonUV;
            buildBoundaryParametricElement(elemIndex, numNodesElementFE, patchIndex, polygonUV);
            ClipperAdapter::cleanPolygon(polygonUV);
            bool isIntegrated = computeLocalCouplingMatrix(elemIndex, patchIndex, polygonUV);
            if(isIntegrated) {
                elementIntegrated.insert(elemIndex);
                projectedPolygons[elemIndex][patchIndex]=polygonUV;
            }
        } // end of loop over set of split patch
    } // end of loop over all the element
    time(&timeEnd);
    INFO_OUT() << "Computing coupling matrices done! It took " << difftime(timeEnd, timeStart) << " seconds" << endl;
    if(elementIntegrated.size() != meshFE->numElems) {
        WARNING_OUT()<<"Number of FE mesh integrated is "<<elementIntegrated.size()<<" over "<<meshFE->numElems<<endl;
        for(int i = 0; i < meshFE->numElems; i++) {
            if(!elementIntegrated.count(i))
                WARNING_OUT()<<"Missing element number "<< i <<endl;
        }
        WARNING_BLOCK_OUT("IGAMortarMapper","ComputeCouplingMatrices","Not all element in FE mesh integrated ! Coupling matrices invalid");
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
            bool isNodeOnPatch = projectedCoords[nodeIndex].find(patchCount)
                    != projectedCoords[nodeIndex].end();
            // Update flag
            if (!isNodeOnPatch) {
                isAllNodesOnPatch = false;
            } else {
                isAllNodesOut = false;
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

void IGAMortarMapper::buildFullParametricElement(int elemCount, int numNodesElementFE, int patchIndex, Polygon2D& polygonUV) {
    // Just look into projectedCoords structure and build the polygon
    for (int nodeCount = 0; nodeCount < numNodesElementFE; nodeCount++) {
        int nodeIndex = meshFEDirectElemTable[elemCount][nodeCount];
        double u = projectedCoords[nodeIndex][patchIndex][0];
        double v = projectedCoords[nodeIndex][patchIndex][1];
        polygonUV.push_back(make_pair(u,v));
    }
}

void IGAMortarMapper::buildBoundaryParametricElement(int elemIndex, int numNodesElementFE, int patchIndex, Polygon2D& polygonUV) {
    IGAPatchSurface* thePatch = meshIGA->getSurfacePatch(patchIndex);
    /// Split nodes from element into 2 subsets, node projected inside the NURBS patch, and node projected outside
    vector<int> insideNode, outsideNode;
    for(int nodeCount = 0; nodeCount < numNodesElementFE; nodeCount++) {
        int nodeIndex = meshFEDirectElemTable[elemIndex][nodeCount];
        /// Flag on whether the node is inside the current NURBS patch
        bool isNodeInsidePatch = projectedCoords[nodeIndex].find(patchIndex)
                != projectedCoords[nodeIndex].end();
        if(isNodeInsidePatch)
            insideNode.push_back(nodeIndex);
        else
            outsideNode.push_back(nodeIndex);
    }
    /// Tolerance validity projected point
    /// WARNING hard coded tolerance for valid
    double toleranceRatio = 1e-6;
    /// Process every node of the element
    for(int i = 0; i < numNodesElementFE; i++) {
        bool isProjectedOnPatchBoundary = true;
        double uIn, vIn;
        double u, v;
        double div = 0, dis = propProjection.maxProjectionDistance;

        /// Get the node indexes
        int nodeCount = (i + 0) % numNodesElementFE;
        int nodeCountPrev = (i + numNodesElementFE - 1) % numNodesElementFE;
        int nodeCountNext = (i + 1) % numNodesElementFE;
        int nodeIndex = meshFEDirectElemTable[elemIndex][nodeCount];
        int nodeIndexPrev = meshFEDirectElemTable[elemIndex][nodeCountPrev];
        int nodeIndexNext = meshFEDirectElemTable[elemIndex][nodeCountNext];

        /// Flag on whether the node and its neighbour is inside the current NURBS patch
        bool isNodeInsidePatch = projectedCoords[nodeIndex].find(patchIndex)
                != projectedCoords[nodeIndex].end();
        bool isPrevNodeInsidePatch = projectedCoords[nodeIndexPrev].find(patchIndex)
                != projectedCoords[nodeIndexPrev].end();
        bool isNextNodeInsidePatch = projectedCoords[nodeIndexNext].find(patchIndex)
                != projectedCoords[nodeIndexNext].end();

        /// Get the location of the node and its neighbour
        double* P0 = &(meshFE->nodes[nodeIndexPrev * 3]);
        double* P1 = &(meshFE->nodes[nodeIndex * 3]);
        double* P2 = &(meshFE->nodes[nodeIndexNext * 3]);

        /// Case selection
        /// Node inside
        if(isNodeInsidePatch) {
            u = projectedCoords[nodeIndex][patchIndex][0];
            v = projectedCoords[nodeIndex][patchIndex][1];
            polygonUV.push_back(make_pair(u,v));
            continue;
        }
        /// Node outside and both neighbors inside
        if(!isNodeInsidePatch && isPrevNodeInsidePatch && isNextNodeInsidePatch) {
            double u0In = projectedCoords[nodeIndexPrev][patchIndex][0];
            double v0In = projectedCoords[nodeIndexPrev][patchIndex][1];
            u = u0In;
            v = v0In;
            dis = propProjection.maxProjectionDistance;
            isProjectedOnPatchBoundary = projectLineOnPatchBoundary(thePatch, u, v, div, dis, P0, P1);
            double u0=u,v0=v,div0=div;
            double u2In = projectedCoords[nodeIndexNext][patchIndex][0];
            double v2In = projectedCoords[nodeIndexNext][patchIndex][1];
            u = u2In;
            v = v2In;
            dis = propProjection.maxProjectionDistance;
            isProjectedOnPatchBoundary = projectLineOnPatchBoundary(thePatch, u, v, div, dis, P2, P1);
            double u2 = u, v2 = v, div2 = div;
            double denominator = (u0In - u0)*(v2In - v2) - (v0In - v0)*(u2In - u2);
            /// If two valid line parameter found and the denominator is valid
            if(div0 >= toleranceRatio && div2 >= toleranceRatio && fabs(denominator) > toleranceRatio) {
                // Compute intersection of the two lines
                // See http://en.wikipedia.org/wiki/Line%E2%80%93line_intersection
                u = ((u0In*v0 - v0In*u0)*(u2In - u2) - (u0In - u0)*(u2In*v2 - v2In*u2))/denominator;
                v = ((u0In*v0 - v0In*u0)*(v2In - v2) - (v0In - v0)*(u2In*v2 - v2In*u2))/denominator;
                // Store the point
                polygonUV.push_back(make_pair(u,v));
                continue;
                // If only first line parameter valid
            } else if(div0 >= toleranceRatio) {
                /// Save data from first projection
                u = u0;
                v = v0;
                uIn = u0In;
                vIn = v0In;
                div = div0;
                // If only second line parameter valid
            } else if(div2 >= toleranceRatio) {
                // Save data from second projection
                u = u2;
                v = v2;
                uIn = u2In;
                vIn = v2In;
                div = div2;
            }
            /// Else no valid point found
        }
        /// Node outside and previous neighbor outside and next neighbor inside
        if(!isNodeInsidePatch && !isPrevNodeInsidePatch && isNextNodeInsidePatch) {
            // Set up initial guess from next node
            uIn = projectedCoords[nodeIndexNext][patchIndex][0];
            vIn = projectedCoords[nodeIndexNext][patchIndex][1];
            u = uIn;
            v = vIn;
            dis = propProjection.maxProjectionDistance;
            // Project on boundary
            isProjectedOnPatchBoundary = projectLineOnPatchBoundary(thePatch, u, v, div, dis, P2, P1);
        }
        /// Node outside and previous neighbor outside and next neighbor inside
        if(!isNodeInsidePatch && isPrevNodeInsidePatch && !isNextNodeInsidePatch) {
            // Set up initial guess from previous node
            uIn = projectedCoords[nodeIndexPrev][patchIndex][0];
            vIn = projectedCoords[nodeIndexPrev][patchIndex][1];
            u = uIn;
            v = vIn;
            dis = propProjection.maxProjectionDistance;
            // Project on boundary
            isProjectedOnPatchBoundary = projectLineOnPatchBoundary(thePatch, u, v, div, dis, P0, P1);
        }
        /// Node outside and both neighbor outside  or no valid line parameter div computed before
        if(div < toleranceRatio) {
            // Project on boundary with all possible inside node until line parameter is valid
            for(vector<int>::iterator it = insideNode.begin(); it != insideNode.end() && div < toleranceRatio; it++) {
                if(*it == nodeIndexPrev || *it == nodeIndexNext)
                    continue;
                double* P0 = &(meshFE->nodes[*it * 3]);
                uIn = projectedCoords[*it][patchIndex][0];
                vIn = projectedCoords[*it][patchIndex][1];
                u = uIn;
                v = vIn;
                dis = propProjection.maxProjectionDistance;
                isProjectedOnPatchBoundary = projectLineOnPatchBoundary(thePatch, u, v, div, dis, P0, P1);
            }
        }
        /// Add point in polygon if line parameter is valid
        if(div >= toleranceRatio) {
            u = uIn + (u - uIn)/div;
            v = vIn + (v - vIn)/div;
            polygonUV.push_back(make_pair(u,v));
        }
        /// Warning/Error output
        if(!isProjectedOnPatchBoundary) {
            if(thePatch->isTrimmed()) {
                DEBUG_OUT() << "Warning in IGAMortarMapper::buildBoundaryParametricElement"
                              << endl;
                DEBUG_OUT() << "Cannot find point projection on patch boundary. "
                              << "Element "<< elemIndex <<" on Patch "<< patchIndex <<" not integrated and skipped !" << endl;
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
                exit(EXIT_FAILURE);
            }
        }
    }
}

bool IGAMortarMapper::projectLineOnPatchBoundary(IGAPatchSurface* thePatch, double& u, double& v, double& div, double& dis, double* Pin, double* Pout) {
    double uIn = u;
    double vIn = v;
    bool isProjectedOnPatchBoundary = thePatch->computePointProjectionOnPatchBoundaryNewtonRhapson(u, v, div, dis, Pin, Pout,
                                                                                                   propNewtonRaphsonBoundary.noIterations, propNewtonRaphsonBoundary.tolProjection);
    if(!isProjectedOnPatchBoundary || dis > propProjection.maxProjectionDistance) {
        DEBUG_OUT() << "In IGAMortarMapper::projectLineOnPatchBoundary. Point projection on boundary using Newton-Rhapson did not converge. Trying bisection algorithm." << endl;
        // Reset initial guess
        u = uIn;
        v = vIn;
        isProjectedOnPatchBoundary = thePatch->computePointProjectionOnPatchBoundaryBisection(u, v, div, dis, Pin, Pout,
                                                                                              propBisection.noIterations, propBisection.tolProjection);
    }
    // Perform some validity check
    if(!isProjectedOnPatchBoundary)
        DEBUG_OUT() << "In IGAMortarMapper::projectLineOnPatchBoundary. Point projection on boundary did not converge. Relax newtonRaphsonBoundary and/or bisection parameters in XML input!"<<endl;
    if(isProjectedOnPatchBoundary && dis > propProjection.maxProjectionDistance) {
        DEBUG_OUT() << "IGAMortarMapper::projectLineOnPatchBoundary. Point projection on boundary found too far. Distance to edge is "<< dis << " for prescribed max of " <<
                       propProjection.maxProjectionDistance << ". Relax maxProjectionDistance in XML input!" << endl;
    }
    return isProjectedOnPatchBoundary;
}

bool IGAMortarMapper::computeLocalCouplingMatrix(const int _elemIndex, const int _patchIndex, Polygon2D& _projectedElement) {
    bool isIntegrated=false;
    // Proceed further if the polygon is valid, i.e. at least a triangle
    if (_projectedElement.size() < 3)
        return isIntegrated;
    IGAPatchSurface* thePatch = meshIGA->getSurfacePatch(_patchIndex);
    Polygon2D projectedElementOnPatch = _projectedElement;
    /// 1.0 Apply patch boundary on polygon
    clipByPatch(thePatch, projectedElementOnPatch);
    ClipperAdapter::cleanPolygon(projectedElementOnPatch);
    // Proceed further if the polygon is valid, i.e. at least a triangle
    if (projectedElementOnPatch.size() < 3)
        return isIntegrated;
    /// 1.1 Init list of trimmed polyggons in case patch is not trimmed
    ListPolygon2D listTrimmedPolygonUV(1, projectedElementOnPatch);
    /// 1.2 Apply trimming
    if(thePatch->isTrimmed())
        clipByTrimming(thePatch,projectedElementOnPatch,listTrimmedPolygonUV);
    /// Debug data
    trimmedProjectedPolygons[_patchIndex].insert(trimmedProjectedPolygons[_patchIndex].end(),listTrimmedPolygonUV.begin(), listTrimmedPolygonUV.end());
    /// 1.3 For each subelement output of the trimmed polygon, clip by knot span
    for(int trimmedPolygonIndex=0;trimmedPolygonIndex<listTrimmedPolygonUV.size();trimmedPolygonIndex++) {
        Polygon2D listSpan;
        ListPolygon2D listPolygonUV;
        /// 1.3.1 Clip by knot span
        clipByKnotSpan(thePatch,listTrimmedPolygonUV[trimmedPolygonIndex],listPolygonUV,listSpan);
        /// 1.3.2 For each subelement clipped by knot span, compute canonical element and integrate
        for(int index=0;index<listSpan.size();index++) {
            // ClipperAdapter::cleanPolygon(listPolygonUV[index],1e-9);
            if(listPolygonUV[index].size()<3)
                continue;
            isIntegrated=true;
            ListPolygon2D triangulatedPolygons = triangulatePolygon(listPolygonUV[index]);
            /// 1.3.3 For each triangle, compute canonical element and integrate
            for(ListPolygon2D::iterator triangulatedPolygon=triangulatedPolygons.begin();
                triangulatedPolygon != triangulatedPolygons.end(); triangulatedPolygon++) {
                /// WARNING hard coded tolerance. Cleaning of triangle. Avoid heavily distorted triangle to go further.
                ClipperAdapter::cleanPolygon(*triangulatedPolygon,1e-8);
                if(triangulatedPolygon->size()<3)
                    continue;
                triangulatedProjectedPolygons[_elemIndex][_patchIndex].push_back(*triangulatedPolygon);
                triangulatedProjectedPolygons2[_patchIndex].push_back(*triangulatedPolygon);
                // Get canonical element
                Polygon2D polygonWZ = computeCanonicalElement(_elemIndex, _projectedElement, *triangulatedPolygon);
                // Integrate
                integrate(thePatch,*triangulatedPolygon,listSpan[index].first,listSpan[index].second,polygonWZ,_elemIndex);
            }
        }
    }
    return isIntegrated;
}

void IGAMortarMapper::clipByPatch(const IGAPatchSurface* _thePatch, Polygon2D& _polygonUV) {
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
    _polygonUV = c.clip(_polygonUV,knotSpanWindow);
}

void IGAMortarMapper::clipByTrimming(const IGAPatchSurface* _thePatch, const Polygon2D& _polygonUV, ListPolygon2D& _listPolygonUV) {
    ClipperAdapter c;
    // Fill clipper with trimming clipping curves
    for(int loop=0;loop<_thePatch->getTrimming().getNumOfLoops();loop++) {
        const std::vector<double> clippingWindow=_thePatch->getTrimming().getLoop(loop).getPolylines();
        c.addPathClipper(clippingWindow);
    }
    // Setup filling rule to have for sure clockwise loop as hole and counterclockwise as boundaries
    c.setFilling(ClipperAdapter::POSITIVE, 0);
    c.addPathSubject(_polygonUV);
    c.clip();
    c.getSolution(_listPolygonUV);
}

void IGAMortarMapper::clipByTrimming(const IGAPatchSurfaceTrimmingLoop* _theTrimmingLoop, const Polygon2D& _polygonUV, ListPolygon2D& _listPolygonUV) {
    ClipperAdapter c;
    // Define the clipping window for the clipper
    const std::vector<double> clippingWindow = _theTrimmingLoop->getPolylines();
    c.addPathClipper(clippingWindow);
    // Setup filling rule to have for sure clockwise loop as hole and counterclockwise as boundaries
    c.setFilling(ClipperAdapter::NEGATIVE, 0);
    c.addPathSubject(_polygonUV);
    c.clip();
    ListPolygon2D tmpListPolygonUV;
    c.getSolution(tmpListPolygonUV);
    _listPolygonUV.insert(_listPolygonUV.end(),tmpListPolygonUV.begin(),tmpListPolygonUV.end());
}

void IGAMortarMapper::clipByKnotSpan(const IGAPatchSurface* _thePatch, const Polygon2D& _polygonUV, ListPolygon2D& _listPolygon, Polygon2D& _listSpan) {
    /// 1.find the knot span which the current element located in.
    //      from minSpanu to maxSpanu in U-direction, and from minSpanV to max SpanV in V-direction
    const double *knotVectorU = _thePatch->getIGABasis()->getUBSplineBasis1D()->getKnotVector();
    const double *knotVectorV = _thePatch->getIGABasis()->getVBSplineBasis1D()->getKnotVector();

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
                /// WARNING hard coded tolerance. Uses reduced clipping tolerance (initially 1e-12). Avoid numerical instability on knot span inside/outside.
                ClipperAdapter c(1e-9);
                if (knotVectorU[spanU] != knotVectorU[spanU + 1]
                        && knotVectorV[spanV] != knotVectorV[spanV + 1]) {
                    Polygon2D knotSpanWindow(4);
                    knotSpanWindow[0] = make_pair(knotVectorU[spanU],knotVectorV[spanV]);
                    knotSpanWindow[1] = make_pair(knotVectorU[spanU+1],knotVectorV[spanV]);
                    knotSpanWindow[2] = make_pair(knotVectorU[spanU+1],knotVectorV[spanV+1]);
                    knotSpanWindow[3] = make_pair(knotVectorU[spanU],knotVectorV[spanV+1]);
                    /// WARNING design. Here we assume to get only a single output polygon from the clipping !
                    Polygon2D solution = c.clip(_polygonUV,knotSpanWindow);

                    /// Store polygon and its knot span for integration
                    _listPolygon.push_back(solution);
                    _listSpan.push_back(make_pair(spanU,spanV));
                }
            }
        }
    }
}

IGAMortarMapper::ListPolygon2D IGAMortarMapper::triangulatePolygon(const Polygon2D& _polygonUV) {
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

IGAMortarMapper::Polygon2D IGAMortarMapper::computeCanonicalElement(const int _elementIndex, const Polygon2D& _theElement, const Polygon2D& _polygonUV) {
    int numNodesElementFE = meshFE->numNodesPerElem[_elementIndex];
    double elementFEUV[8], coordsNodeFEUV[2], coordsNodeFEWZ[2];
    for(int i = 0; i < numNodesElementFE; i++) {
        elementFEUV[2*i]   = _theElement[i].first;
        elementFEUV[2*i+1] = _theElement[i].second;
    }
    // Compute canonical coordinates polygon using parametric space
    Polygon2D polygonWZ;
    for(int i=0; i<_polygonUV.size(); i++) {
        coordsNodeFEUV[0] = _polygonUV[i].first;
        coordsNodeFEUV[1] = _polygonUV[i].second;
        if(numNodesElementFE == 3)
            MathLibrary::computeLocalCoordsInTriangle(elementFEUV, coordsNodeFEUV, coordsNodeFEWZ);
        else
            MathLibrary::computeLocalCoordsInQuad(elementFEUV, coordsNodeFEUV, coordsNodeFEWZ);
        polygonWZ.push_back(make_pair(coordsNodeFEWZ[0],coordsNodeFEWZ[1]));
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
     *       2.2.8 integrate the shape function product for Cnn(Linear shape function multiply linear shape function)
     *       2.2.9 integrate the shape function product for Cnr(Linear shape function multiply IGA shape function)
     * 3. Assemble the element coupling matrix to the global coupling matrix.
     */

    // Read input
    assert(!_polygonUV.empty());
    assert(!_polygonWZ.empty());
    int numNodesUV=_polygonUV.size();
    int numNodesWZ=_polygonWZ.size();
    assert(numNodesUV > 2);
    assert(numNodesUV < 5);
    assert(numNodesWZ > 2);
    assert(numNodesWZ < 5);

    // Number of Cartesian directions
    int noCoord = 3;

    // Definitions
    int numNodesElementFE = meshFE->numNodesPerElem[_elementIndex];
    int numNodesElMaster = 0;
    int numNodesElSlave = 0;

    int pDegree = _thePatch->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree();
    int qDegree = _thePatch->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree();
    int nShapeFuncsIGA = (pDegree + 1) * (qDegree + 1);

    if (isMappingIGA2FEM) {
        numNodesElMaster = numNodesElementFE;
        numNodesElSlave = nShapeFuncsIGA;
    } else {
        numNodesElMaster = nShapeFuncsIGA;
        numNodesElSlave = numNodesElementFE;
    }
    // Initialize element matrix
    double elementCouplingMatrixNN[numNodesElMaster * (numNodesElMaster + 1) / 2];
    double elementCouplingMatrixNR[numNodesElSlave * numNodesElMaster];

    for (int arrayIndex = 0; arrayIndex < numNodesElMaster * (numNodesElMaster + 1) / 2;
         arrayIndex++)
        elementCouplingMatrixNN[arrayIndex] = 0;

    for (int arrayIndex = 0; arrayIndex < numNodesElSlave * numNodesElMaster; arrayIndex++)
        elementCouplingMatrixNR[arrayIndex] = 0;

    int dofIGA[nShapeFuncsIGA];
    _thePatch->getIGABasis()->getBasisFunctionsIndex(_spanU, _spanV, dofIGA);

    for (int i = 0; i < nShapeFuncsIGA; i++)
        dofIGA[i] = _thePatch->getControlPointNet()[dofIGA[i]]->getDofIndex();

    /// 1. Copy input polygon into contiguous C format
    double nodesUV[8];
    double nodesWZ[8];
    for (int i = 0; i < numNodesUV; i++) {
        nodesUV[i*2]=_polygonUV[i].first;
        nodesUV[i*2+1]=_polygonUV[i].second;
        nodesWZ[i*2]=_polygonWZ[i].first;
        nodesWZ[i*2+1]=_polygonWZ[i].second;
    }

    /// 2. Loop through each quadrature
    /// 2.1 Choose gauss triangle or gauss quadriliteral
    MathLibrary::IGAGaussQuadrature *theGaussQuadrature;
    int nNodesQuadrature=numNodesUV;
    if (nNodesQuadrature == 3)
        theGaussQuadrature = gaussTriangle;
    else
        theGaussQuadrature = gaussQuad;

    double *quadratureUV = nodesUV;
    double *quadratureWZ = nodesWZ;

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
        MathLibrary::computeLinearCombination(nNodesQuadrature, 2, quadratureWZ, shapeFuncs,
                                              GPFE);

        /// 2.2.4 compute the shape function(in the linear element) of the current integration point
        double shapeFuncsFE[numNodesElementFE];
        MathLibrary::computeLowOrderShapeFunc(numNodesElementFE, GPFE, shapeFuncsFE);

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

        /// 2.2.8 integrate the shape function product for Cnn(Linear shape function multiply linear shape function)
        int count = 0;
        for (int i = 0; i < numNodesElMaster; i++) {
            for (int j = i; j < numNodesElMaster; j++) {
                if (isMappingIGA2FEM)
                    elementCouplingMatrixNN[count++] += shapeFuncsFE[i] * shapeFuncsFE[j]
                            * Jacobian * theGaussQuadrature->weights[GPCount];
                else {
                    double IGABasisFctsI =
                            localBasisFunctionsAndDerivatives[_thePatch->getIGABasis()->indexDerivativeBasisFunction(
                                1, 0, 0, i)];
                    double IGABasisFctsJ =
                            localBasisFunctionsAndDerivatives[_thePatch->getIGABasis()->indexDerivativeBasisFunction(
                                1, 0, 0, j)];
                    elementCouplingMatrixNN[count++] += IGABasisFctsI * IGABasisFctsJ * Jacobian
                            * theGaussQuadrature->weights[GPCount];
                }
            }
        }

        /// Save GP data
        if(propErrorComputation.isDomainError){
            std::vector<double> streamGP;
            // weight + jacobian + nShapeFuncsFE + (#dof, shapefuncvalue,...) + nShapeFuncsIGA + (#dof, shapefuncvalue,...)
            streamGP.reserve(1 + 1 + 1 + 2*numNodesElementFE + 1 + 2*nShapeFuncsIGA);
            streamGP.push_back(theGaussQuadrature->weights[GPCount]);
            streamGP.push_back(Jacobian);
            streamGP.push_back(numNodesElementFE);
            for (int i = 0; i < numNodesElementFE; i++) {
                streamGP.push_back(meshFEDirectElemTable[_elementIndex][i]);
                streamGP.push_back(shapeFuncsFE[i]);
            }
            streamGP.push_back(nShapeFuncsIGA);
            for (int i = 0; i < nShapeFuncsIGA; i++) {
                double IGABasisFctsI =
                        localBasisFunctionsAndDerivatives[_thePatch->getIGABasis()->indexDerivativeBasisFunction(
                            1, 0, 0, i)];
                streamGP.push_back(dofIGA[i]);
                streamGP.push_back(IGABasisFctsI);
            }
            streamGPs.push_back(streamGP);
        }

        /// 2.2.9 integrate the shape function product for C_NR(Linear shape function multiply IGA shape function)
        count = 0;
        for (int i = 0; i < numNodesElMaster; i++) {
            for (int j = 0; j < numNodesElSlave; j++) {
                double basisFctsMaster;
                double basisFctsSlave;
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
    ///Create new scope
    {
        int dofIGA[nShapeFuncsIGA];
        _thePatch->getIGABasis()->getBasisFunctionsIndex(_spanU, _spanV, dofIGA);

        for (int i = 0; i < nShapeFuncsIGA; i++)
            dofIGA[i] = _thePatch->getControlPointNet()[dofIGA[i]]->getDofIndex();

        int count = 0;
        int dof1, dof2;
        /// 3.1 Assemble Cnn
        for (int i = 0; i < numNodesElMaster; i++)
            for (int j = i; j < numNodesElMaster; j++) {
                if (isMappingIGA2FEM) {
                    dof1 = meshFEDirectElemTable[_elementIndex][i];
                    dof2 = meshFEDirectElemTable[_elementIndex][j];
                } else {
                    dof1 = dofIGA[i];
                    dof2 = dofIGA[j];
                }
                if (!isExpanded){
                    couplingMatrices->addCNNValue(dof1, dof2, elementCouplingMatrixNN[count]);
                    if (dof1 != dof2)
                        couplingMatrices->addCNNValue(dof2, dof1, elementCouplingMatrixNN[count]);
                } else {
                    for(int iCoord = 0; iCoord < noCoord; iCoord++){
                        couplingMatrices->addCNNValue(noCoord*dof1 + iCoord, noCoord*dof2 + iCoord, elementCouplingMatrixNN[count]);
                        if (dof1 != dof2)
                            couplingMatrices->addCNNValue(noCoord*dof2 + iCoord, noCoord*dof1 + iCoord, elementCouplingMatrixNN[count]);
                    }
                }
                count++;
            }

        count = 0;
        /// 3.1 Assemble Cnr
        for (int i = 0; i < numNodesElMaster; i++)
            for (int j = 0; j < numNodesElSlave; j++) {
                if (isMappingIGA2FEM) {
                    dof1 = meshFEDirectElemTable[_elementIndex][i];
                    dof2 = dofIGA[j];
                } else {
                    dof1 = dofIGA[i];
                    dof2 = meshFEDirectElemTable[_elementIndex][j];
                }
                if (!isExpanded){
                    couplingMatrices->addCNRValue(dof1, dof2, elementCouplingMatrixNR[count]);
                } else {
                    for(int iCoord = 0; iCoord < noCoord; iCoord++)
                        couplingMatrices->addCNRValue(noCoord*dof1 + iCoord, noCoord*dof2 + iCoord, elementCouplingMatrixNR[count]);
                }
                count ++;
            }
    }
}

bool IGAMortarMapper::computeKnotSpanOfProjElement(const IGAPatchSurface* _thePatch, const Polygon2D& _polygonUV, int* _span) {
    int minSpanU = _thePatch->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(
                _polygonUV[0].first);
    int minSpanV = _thePatch->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(
                _polygonUV[0].second);
    int maxSpanU = minSpanU;
    int maxSpanV = minSpanV;

    for (int nodeCount = 1; nodeCount < _polygonUV.size(); nodeCount++) {
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

void IGAMortarMapper::computeIGAWeakDirichletCurveConditionMatrices() {
    /*
     * Computes and assembles the patch weak Dirichlet curve conditions.
     */

    // Get the weak Dirichlet curve conditions
    std::vector<WeakIGADirichletCurveCondition*> weakIGADirichletCurveConditions = meshIGA->getWeakIGADirichletCurveConditions();

    // Initialize constant array sizes
    const int noCoordParam = 2;
    const int noCoord = 3;

    // Initialize varying array sizes
    int patchIndex;
    int counter;
    int p;
    int q;
    int noLocalBasisFcts;
    int noDOFsLoc;
    int noGPsOnCond;
    int uKnotSpan;
    int vKnotSpan;
    int indexCP;
    double uGP;
    double vGP;
    double tangentCurveVct[noCoord];
    double normalCurveVct[noCoord];
    double alphaPrimary;
    double alphaSecondaryBending;
    double alphaSecondaryTwisting;
    double elementLengthOnGP;

    // Initialize pointers
    double* curveGPs;
    double* curveGPWeights;
    double* curveGPTangents;
    double* curveGPJacobianProducts;
    IGAPatchSurface* thePatch;

    // Loop over all the conditions for the application of weak Dirichlet conditions
    for (int iDCC = 0; iDCC < weakIGADirichletCurveConditions.size(); iDCC++){
        // Get the penalty factors for the primary and the secondary field
        alphaPrimary = weakDirichletCCAlphaPrimary[iDCC];
        alphaSecondaryBending = weakDirichletCCAlphaSecondaryBending[iDCC];
        alphaSecondaryTwisting = weakDirichletCCAlphaSecondaryTwisting[iDCC];

        // Get the index of the patch
        patchIndex = weakIGADirichletCurveConditions[iDCC]->getPatchIndex();

        // Get the number of Gauss Points for the given condition
        noGPsOnCond = weakIGADirichletCurveConditions[iDCC]->getCurveNumGP();

        // Get the parametric coordinates of the Gauss Points
        curveGPs = weakIGADirichletCurveConditions[iDCC]->getCurveGPs();

        // Get the corresponding Gauss weights
        curveGPWeights = weakIGADirichletCurveConditions[iDCC]->getCurveGPWeights();

        // Get the tangent vectors at the trimming curve of the given condition in the Cartesian space
        curveGPTangents = weakIGADirichletCurveConditions[iDCC]->getCurveGPTangents();

        // Get the product of the Jacobian transformations
        curveGPJacobianProducts = weakIGADirichletCurveConditions[iDCC]->getCurveGPJacobianProducts();

        // Get the patch
        thePatch = meshIGA->getSurfacePatch(patchIndex);

        // Get the polynomial orders of the master and the slave patch
        p = thePatch->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree();
        q = thePatch->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree();

        // get the number of local basis functions
        noLocalBasisFcts = (p + 1)*(q + 1);

        // get the number of the local DOFs for the patch
        noDOFsLoc = noCoord*noLocalBasisFcts;

        // Initialize pointers
        double* BDisplacementsGC = new double[noCoord*noDOFsLoc];
        double* BOperatorOmegaT = new double[noDOFsLoc];
        double* BOperatorOmegaN = new double[noDOFsLoc];
        double* KPenaltyDisplacement = new double[noDOFsLoc*noDOFsLoc];
        double* KPenaltyBendingRotation = new double[noDOFsLoc*noDOFsLoc];
        double* KPenaltyTwistingRotation = new double[noDOFsLoc*noDOFsLoc];

        // Loop over all the Gauss Points of the given condition
        for(int iGP = 0; iGP < noGPsOnCond; iGP++){

            // Get the parametric coordinates of the Gauss Point on the patch
            uGP = curveGPs[iGP*noCoordParam];
            vGP = curveGPs[iGP*noCoordParam + 1];

            // Find the knot span indices of the Gauss point locations in the parameter space of the patch
            uKnotSpan = thePatch->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(uGP);
            vKnotSpan = thePatch->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(vGP);

            // Get the tangent to the boundary vector on the patch
            for(int iCoord = 0; iCoord < noCoord; iCoord++)
                tangentCurveVct[iCoord] = curveGPTangents[iGP*noCoord + iCoord];

            // compute elementLength on GP. The weight is already included in variable trCurveGPJacobianProducts
            elementLengthOnGP = curveGPJacobianProducts[iGP];

            // Compute the B-operator matrices needed for the computation of the patch weak Dirichlet conditions at the patch
            computeDisplacementAndRotationBOperatorMatrices(BDisplacementsGC, BOperatorOmegaT, BOperatorOmegaN, normalCurveVct,
                                                                thePatch, tangentCurveVct, uGP, vGP, uKnotSpan, vKnotSpan);

            // Compute the dual product matrices for the displacements
            EMPIRE::MathLibrary::computeTransposeMatrixProduct(noCoord,noDOFsLoc,noDOFsLoc,BDisplacementsGC,BDisplacementsGC,KPenaltyDisplacement);

            // Compute the dual product matrices for the bending rotations
            EMPIRE::MathLibrary::computeTransposeMatrixProduct(1, noDOFsLoc, noDOFsLoc, BOperatorOmegaT, BOperatorOmegaT, KPenaltyBendingRotation);

            // Compute the dual product matrices for the twisting rotations
            EMPIRE::MathLibrary::computeTransposeMatrixProduct(1, noDOFsLoc, noDOFsLoc, BOperatorOmegaN, BOperatorOmegaN, KPenaltyTwistingRotation);

            // Compute the element index tables for the patch
            int CPIndex[noLocalBasisFcts];
            thePatch->getIGABasis()->getBasisFunctionsIndex(uKnotSpan, vKnotSpan, CPIndex);

            // Compute the element freedom tables for the patch
            int EFT[noDOFsLoc];
            counter = 0;
            for (int i = 0; i < noLocalBasisFcts; i++){
                indexCP = thePatch->getControlPointNet()[CPIndex[i]]->getDofIndex();
                for (int j = 0; j < noCoord; j++){
                    EFT[counter] = noCoord*indexCP + j;
                    counter++;
                }
            }

            // Assemble KPenaltyDisplacement to the global coupling matrix CNN
            for(int i = 0; i < noDOFsLoc; i++){
                for(int j = 0; j < noDOFsLoc; j++){
                    // Assemble the displacement coupling entries
                    couplingMatrices->addCNNValue(EFT[i], EFT[i], alphaPrimary*KPenaltyDisplacement[i*noDOFsLoc + i]*elementLengthOnGP);

                    // Assemble the bending rotation coupling entries
                    couplingMatrices->addCNNValue(EFT[i], EFT[j], alphaSecondaryBending*KPenaltyBendingRotation[i*noDOFsLoc + j]*elementLengthOnGP);

                    // Assemble the twisting rotation coupling entries
                    couplingMatrices->addCNNValue(EFT[i], EFT[j], alphaSecondaryTwisting*KPenaltyTwistingRotation[i*noDOFsLoc + j]*elementLengthOnGP);
                }
            }

            // Compute the error along the trimming curves where constraints are applied
            if(propErrorComputation.isCurveError){
                // Initialize variable storing the Gauss Point data
                std::vector<double> streamCurveGP;

                // elementLengthOnGP + noBasisFuncs + (#indexCP, basisFuncValue,...) + (#indexDOF, BtValue, BnValue,...)
                streamCurveGP.reserve(1 + 1 + 2*noLocalBasisFcts + 3*noDOFsLoc);

                // Save the element length on the Gauss Point
                streamCurveGP.push_back(elementLengthOnGP);

                // Save the number of basis functions
                streamCurveGP.push_back(noLocalBasisFcts);

                // Save the Control Point index and the basis function's values
                for(int iBFs = 0; iBFs < noLocalBasisFcts; iBFs++){
                    indexCP = thePatch->getControlPointNet()[CPIndex[iBFs]]->getDofIndex();
                    streamCurveGP.push_back(indexCP);
                    streamCurveGP.push_back(BDisplacementsGC[0*noLocalBasisFcts + 3*iBFs]);
                }

                // Save the DOF index and the bending and twisting B-operator values
                for(int iDOFs = 0; iDOFs < noDOFsLoc; iDOFs++){
                    streamCurveGP.push_back(EFT[iDOFs]);
                    streamCurveGP.push_back(BOperatorOmegaT[iDOFs]);
                    streamCurveGP.push_back(BOperatorOmegaN[iDOFs]);
                }

                // Push back the Gauss Point values into the member variable
                streamCurveGPs.push_back(streamCurveGP);
            }

        } // End of Gauss Point loop

        // Delete pointers
        delete[] BDisplacementsGC;
        delete[] BOperatorOmegaT;
        delete[] BOperatorOmegaN;
        delete[] KPenaltyDisplacement;
        delete[] KPenaltyBendingRotation;
        delete[] KPenaltyTwistingRotation;

    } // End of weak Dirichlet curve condition loop
}

void IGAMortarMapper::computeIGAWeakDirichletSurfaceConditionMatrices() {
    /*
     * Computes and assembles the patch weak Dirichlet surface conditions.
     */
    ERROR_OUT() << "Function under construction" << endl;
    exit(-1);

    // Get the weak Dirichlet curve conditions
    std::vector<WeakIGADirichletSurfaceCondition*> weakIGADirichletSurfaceConditions = meshIGA->getWeakIGADirichletSurfaceConditions();

    // Initialize constant array sizes
    const int noCoordParam = 2;
    const int noCoord = 3;

    // Initialize varying array sizes
    int patchIndex;
    int counter;
    int p;
    int q;
    int noLocalBasisFcts;
    int noDOFsLoc;
    int noGPsOnCond;
    int uKnotSpan;
    int vKnotSpan;
    int indexCP;
    double uGP;
    double vGP;
    double tangentCurveVct[noCoord] = {0,0,0};
    double normalCurveVct[noCoord];
    double alphaPrimary;
    double alphaSecondaryBending;
    double alphaSecondaryTwisting;
    double jacobianOnGP;

    // Initialize pointers
    double* surfaceGPs;
    double* surfaceGPWeights;
    double* surfaceGPJacobians;
    IGAPatchSurface* thePatch;

    // Loop over all the conditions for the application of weak Dirichlet conditions
    for (int iDSC = 0; iDSC < weakIGADirichletSurfaceConditions.size(); iDSC++){
        // Get the penalty factors for the primary and the secondary field
        alphaPrimary = weakDirichletSCAlphaPrimary[iDSC];
        alphaSecondaryBending = weakDirichletSCAlphaSecondaryBending[iDSC];
        alphaSecondaryTwisting = weakDirichletSCAlphaSecondaryTwisting[iDSC];

        // Get the index of the patch
        patchIndex = weakIGADirichletSurfaceConditions[iDSC]->getPatchIndex();

        // Get the number of Gauss Points for the given condition
        noGPsOnCond = weakIGADirichletSurfaceConditions[iDSC]->getSurfaceNumGP();

        // Get the parametric coordinates of the Gauss Points
        surfaceGPs = weakIGADirichletSurfaceConditions[iDSC]->getSurfaceGPs();

        // Get the corresponding Gauss weights
        surfaceGPWeights = weakIGADirichletSurfaceConditions[iDSC]->getSurfaceGPWeights();

        // Get the Jacobians
        surfaceGPJacobians = weakIGADirichletSurfaceConditions[iDSC]->getSurfaceGPJacobians();

        // Get the patch
        thePatch = meshIGA->getSurfacePatch(patchIndex);

        // Get the polynomial orders of the master and the slave patch
        p = thePatch->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree();
        q = thePatch->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree();

        // get the number of local basis functions
        noLocalBasisFcts = (p + 1)*(q + 1);

        // get the number of the local DOFs for the patch
        noDOFsLoc = noCoord*noLocalBasisFcts;

        // Initialize pointers
        double* BDisplacementsGC = new double[noCoord*noDOFsLoc];
        double* BOperatorOmegaT = new double[noDOFsLoc];
        double* BOperatorOmegaN = new double[noDOFsLoc];
        double* KPenaltyDisplacement = new double[noDOFsLoc*noDOFsLoc];

        // Loop over all the Gauss Points of the given condition
        for(int iGP = 0; iGP < noGPsOnCond; iGP++){

            // Get the parametric coordinates of the Gauss Point on the patch
            uGP = surfaceGPs[iGP*noCoordParam];
            vGP = surfaceGPs[iGP*noCoordParam + 1];

            // Find the knot span indices of the Gauss point locations in the parameter space of the patch
            uKnotSpan = thePatch->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(uGP);
            vKnotSpan = thePatch->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(vGP);

            // compute elementLength on GP. The weight is already included in variable trCurveGPJacobianProducts
            jacobianOnGP = surfaceGPJacobians[iGP];

            // Compute the B-operator matrices needed for the computation of the patch weak Dirichlet conditions at the patch
            computeDisplacementAndRotationBOperatorMatrices(BDisplacementsGC, BOperatorOmegaT, BOperatorOmegaN, normalCurveVct,
                                                                thePatch, tangentCurveVct, uGP, vGP, uKnotSpan, vKnotSpan);

            // Compute the dual product matrices for the displacements
            EMPIRE::MathLibrary::computeTransposeMatrixProduct(noCoord,noDOFsLoc,noDOFsLoc,BDisplacementsGC,BDisplacementsGC,KPenaltyDisplacement);

            // Compute the element index tables for the patch
            int CPIndex[noLocalBasisFcts];
            thePatch->getIGABasis()->getBasisFunctionsIndex(uKnotSpan, vKnotSpan, CPIndex);

            // Compute the element freedom tables for the patch
            int EFT[noDOFsLoc];
            counter = 0;
            for (int i = 0; i < noLocalBasisFcts; i++){
                indexCP = thePatch->getControlPointNet()[CPIndex[i]]->getDofIndex();
                for (int j = 0; j < noCoord; j++){
                    EFT[counter] = noCoord*indexCP + j;
                    counter++;
                }
            }

            // Assemble KPenaltyDisplacement to the global coupling matrix CNN
            for(int i = 0; i < noDOFsLoc; i++){
                for(int j = 0; j < noDOFsLoc; j++){
                    // Assemble the displacement coupling entries
                    couplingMatrices->addCNNValue(EFT[i], EFT[i], alphaPrimary*KPenaltyDisplacement[i*noDOFsLoc + i]*jacobianOnGP);
                }
            }

            //// TODO
            // Store the GP data into array for later usage in the error computation
        } // End of Gauss Point loop

        // Delete pointers
        delete[] BDisplacementsGC;
        delete[] BOperatorOmegaT;
        delete[] BOperatorOmegaN;
        delete[] KPenaltyDisplacement;

    } // End of weak Dirichlet surface condition loop
}

void IGAMortarMapper::computeIGAPatchWeakContinuityConditionMatrices() {
    /*
     * Computes and assembles the patch weak continuity conditions.
     */

    // Get the weak patch continuity conditions
    std::vector<WeakIGAPatchContinuityCondition*> weakIGAPatchContinuityConditions = meshIGA->getWeakIGAPatchContinuityConditions();

    // Define tolerances
    const double tolAngle = 1e-1;
    const double tolVct = 1e-4;

    // Initialize constant array sizes
    const int noCoordParam = 2;
    const int noCoord = 3;

    // Initialize varying array sizes
    int indexMaster;
    int indexSlave;
    int counter;
    int pMaster;
    int qMaster;
    int pSlave;
    int qSlave;
    int noLocalBasisFctsMaster;
    int noLocalBasisFctsSlave;
    int noDOFsLocMaster;
    int noDOFsLocSlave;
    int noGPsOnContCond;
    int uKnotSpanMaster;
    int uKnotSpanSlave;
    int vKnotSpanMaster;
    int vKnotSpanSlave;
    int indexCP;
    double uGPMaster;
    double vGPMaster;
    double uGPSlave;
    double vGPSlave;
    double phiTangents;
    double phiNormals;
    double condAligned;
    double factorTangent;
    double factorNormal;
    double factorTwisting;
    double tangentTrCurveVctMaster[noCoord];
    double tangentTrCurveVctSlave[noCoord];
    double normalTrCurveVctMaster[noCoord];
    double normalTrCurveVctSlave[noCoord];
    double normTangentTrCurveVctMaster;
    double normTangentTrCurveVctSlave;
    double normNormalTrCurveVctMaster;
    double normNormalTrCurveVctSlave;
    double alphaPrimary;
    double alphaSecondaryBending;
    double alphaSecondaryTwisting;
    double elementLengthOnGP;

    // Initialize pointers
    double* trCurveMasterGPs;
    double* trCurveSlaveGPs;
    double* trCurveGPWeights;
    double* trCurveMasterGPTangents;
    double* trCurveSlaveGPTangents;
    double* trCurveGPJacobianProducts;
    IGAPatchSurface* patchMaster;
    IGAPatchSurface* patchSlave;

    // Loop over all the conditions for the application of weak continuity across patch interfaces
    for (int iWCC = 0; iWCC < weakIGAPatchContinuityConditions.size(); iWCC++){
        // Get the penalty factors for the primary and the secondary field
        alphaPrimary = weakPatchContinuityAlphaPrimaryIJ[iWCC];
        alphaSecondaryBending = weakPatchContinuityAlphaSecondaryBendingIJ[iWCC];
        alphaSecondaryTwisting = weakPatchContinuityAlphaSecondaryTwistingIJ[iWCC];

        // Get the index of the master and slave patches
        indexMaster = weakIGAPatchContinuityConditions[iWCC]->getMasterPatchIndex();
        indexSlave = weakIGAPatchContinuityConditions[iWCC]->getSlavePatchIndex();

        // Get the number of Gauss Points for the given condition
        noGPsOnContCond = weakIGAPatchContinuityConditions[iWCC]->getTrCurveNumGP();

        // Get the parametric coordinates of the Gauss Points
        trCurveMasterGPs = weakIGAPatchContinuityConditions[iWCC]->getTrCurveMasterGPs();
        trCurveSlaveGPs = weakIGAPatchContinuityConditions[iWCC]->getTrCurveSlaveGPs();

        // Get the corresponding Gauss weights
        trCurveGPWeights = weakIGAPatchContinuityConditions[iWCC]->getTrCurveGPWeights();

        // Get the tangent vectors at the trimming curve of the given condition in the Cartesian space
        trCurveMasterGPTangents = weakIGAPatchContinuityConditions[iWCC]->getTrCurveMasterGPTangents();
        trCurveSlaveGPTangents = weakIGAPatchContinuityConditions[iWCC]->getTrCurveSlaveGPTangents();

        // Get the product of the Jacobian transformations
        trCurveGPJacobianProducts = weakIGAPatchContinuityConditions[iWCC]->getTrCurveGPJacobianProducts();

        // Get the master and the slave patch
        patchMaster = meshIGA->getSurfacePatch(indexMaster);
        patchSlave = meshIGA->getSurfacePatch(indexSlave);

        // Get the polynomial orders of the master and the slave patch
        pMaster = patchMaster->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree();
        qMaster = patchMaster->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree();
        pSlave = patchSlave->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree();
        qSlave = patchSlave->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree();

        // get the number of local basis functions for master and slave patch
        noLocalBasisFctsMaster = (pMaster + 1)*(qMaster + 1);
        noLocalBasisFctsSlave = (pSlave + 1)*(qSlave + 1);

        // get the number of the local DOFs for the master and slave patch
        noDOFsLocMaster = noCoord*noLocalBasisFctsMaster;
        noDOFsLocSlave = noCoord*noLocalBasisFctsSlave;

        // Initialize pointers
        double* BDisplacementsGCMaster = new double[noCoord*noDOFsLocMaster];
        double* BDisplacementsGCSlave = new double[noCoord*noDOFsLocSlave];
        double* BOperatorOmegaTMaster = new double[noDOFsLocMaster];
        double* BOperatorOmegaTSlave = new double[noDOFsLocSlave];
        double* BOperatorOmegaNMaster = new double[noDOFsLocMaster];
        double* BOperatorOmegaNSlave = new double[noDOFsLocSlave];
        double* KPenaltyDisplacementMaster = new double[noDOFsLocMaster*noDOFsLocMaster];
        double* KPenaltyDisplacementSlave = new double[noDOFsLocSlave*noDOFsLocSlave];
        double* CPenaltyDisplacement = new double[noDOFsLocMaster*noDOFsLocSlave];
        double* KPenaltyBendingRotationMaster = new double[noDOFsLocMaster*noDOFsLocMaster];
        double* KPenaltyBendingRotationSlave = new double[noDOFsLocSlave*noDOFsLocSlave];
        double* CPenaltyBendingRotation = new double[noDOFsLocMaster*noDOFsLocSlave];
        double* KPenaltyTwistingRotationMaster = new double[noDOFsLocMaster*noDOFsLocMaster];
        double* KPenaltyTwistingRotationSlave = new double[noDOFsLocSlave*noDOFsLocSlave];
        double* CPenaltyTwistingRotation = new double[noDOFsLocMaster*noDOFsLocSlave];

        // Loop over all the Gauss Points of the given condition
        for(int iGP = 0; iGP < noGPsOnContCond; iGP++){
            // Get the parametric coordinates of the Gauss Point on the master patch
            uGPMaster = trCurveMasterGPs[iGP*noCoordParam];
            vGPMaster = trCurveMasterGPs[iGP*noCoordParam + 1];

            // Get the parametric coordinates of the Gauss Point on the slave patch
            uGPSlave = trCurveSlaveGPs[iGP*noCoordParam];
            vGPSlave = trCurveSlaveGPs[iGP*noCoordParam + 1];

            // Find the knot span indices of the Gauss point locations in the parameter space of the master patch
            uKnotSpanMaster = patchMaster->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(uGPMaster);
            vKnotSpanMaster = patchMaster->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(vGPMaster);

            // Find the knot span indices of the Gauss point locations in the parameter space of the slave patch
            uKnotSpanSlave = patchSlave->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(uGPSlave);
            vKnotSpanSlave = patchSlave->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(vGPSlave);

            // Get the tangent to the boundary vector on the master and the slave patch
            for(int iCoord = 0; iCoord < noCoord; iCoord++){
                tangentTrCurveVctMaster[iCoord] = trCurveMasterGPTangents[iGP*noCoord + iCoord];
                tangentTrCurveVctSlave[iCoord] = trCurveSlaveGPTangents[iGP*noCoord + iCoord];
            }

            // compute elementLength on GP. The weight is already included in variable trCurveGPJacobianProducts
            elementLengthOnGP = trCurveGPJacobianProducts[iGP];

            // Compute the B-operator matrices needed for the computation of the patch weak continuity contributions at the master patch
            computeDisplacementAndRotationBOperatorMatrices(BDisplacementsGCMaster, BOperatorOmegaTMaster, BOperatorOmegaNMaster, normalTrCurveVctMaster,
                                                            patchMaster, tangentTrCurveVctMaster, uGPMaster, vGPMaster, uKnotSpanMaster, vKnotSpanMaster);

            // Compute the B-operator matrices needed for the computation of the patch weak continuity contributions at the slave patch
            computeDisplacementAndRotationBOperatorMatrices(BDisplacementsGCSlave, BOperatorOmegaTSlave, BOperatorOmegaNSlave, normalTrCurveVctSlave,
                                                            patchSlave, tangentTrCurveVctSlave, uGPSlave, vGPSlave, uKnotSpanSlave, vKnotSpanSlave);

            // Determine the alignment of the tangent vectors from both patches at their common interface
            normTangentTrCurveVctMaster = EMPIRE::MathLibrary::computeDenseDotProduct(noCoord,tangentTrCurveVctMaster,tangentTrCurveVctMaster);
            normTangentTrCurveVctMaster = sqrt(normTangentTrCurveVctMaster);
            normTangentTrCurveVctSlave = EMPIRE::MathLibrary::computeDenseDotProduct(noCoord,tangentTrCurveVctSlave,tangentTrCurveVctSlave);
            normTangentTrCurveVctSlave = sqrt(normTangentTrCurveVctSlave);
            phiTangents = EMPIRE::MathLibrary::computeDenseDotProduct(noCoord,tangentTrCurveVctMaster,tangentTrCurveVctSlave);
            phiTangents = phiTangents/(normTangentTrCurveVctMaster*normTangentTrCurveVctSlave);

            // Check if the tangent vectors at the coupled trimming curves are zero and if yes go to the next Gauss point
            if(normTangentTrCurveVctMaster < tolVct && normTangentTrCurveVctSlave < tolVct)
                continue;
            else if((normTangentTrCurveVctMaster < tolVct && normTangentTrCurveVctSlave > tolVct) || (normTangentTrCurveVctMaster > tolVct && normTangentTrCurveVctSlave < tolVct))
                assert(false);

            // Check if the tangent vectors are aligned and if not assert error
            condAligned = phiTangents*phiTangents - 1;
            if (abs(condAligned) > tolAngle) {
                INFO_OUT() << "Found boundaries for which the tangent vectors are not aligned with angle = " << phiTangents << endl;
            }
            assert(abs(condAligned) < tolAngle);

            // Check if the angle between the tangent vectors is 0° (have the same direction cos(phi) = 1) or 180° (have opposite directions cos(phi) = - 1) to formulate the interface constraint
            if (abs(phiTangents - 1) < sqrt(tolAngle)) {
                factorTangent = -1.0;
            } else if (abs(phiTangents + 1) < sqrt(tolAngle)) {
                factorTangent = 1.0;
            } else {
                assert(false);
            }

            // Determine the alignment of the normal vectors from both patches at their common interface
            normNormalTrCurveVctMaster = EMPIRE::MathLibrary::computeDenseDotProduct(noCoord,normalTrCurveVctMaster,normalTrCurveVctMaster);
            normNormalTrCurveVctMaster = sqrt(normNormalTrCurveVctMaster);
            normNormalTrCurveVctSlave = EMPIRE::MathLibrary::computeDenseDotProduct(noCoord,normalTrCurveVctSlave,normalTrCurveVctSlave);
            normNormalTrCurveVctSlave = sqrt(normNormalTrCurveVctSlave);
            phiNormals = EMPIRE::MathLibrary::computeDenseDotProduct(noCoord,normalTrCurveVctMaster,normalTrCurveVctSlave);
            phiNormals = phiTangents/(normNormalTrCurveVctMaster*normNormalTrCurveVctSlave);

            // Check if the normal vectors at the coupled trimming curves are zero and if yes go to the next Gauss point
            if(normNormalTrCurveVctMaster < tolVct && normNormalTrCurveVctSlave < tolVct)
                continue;
            else if((normNormalTrCurveVctMaster < tolVct && normNormalTrCurveVctSlave > tolVct) || (normNormalTrCurveVctMaster > tolVct && normNormalTrCurveVctSlave < tolVct))
                assert(false);

            // Check if the normal vectors are not aligned and in this case neglect the interface twisting rotation continuity condition
            condAligned = phiNormals*phiNormals - 1;
            if (abs(condAligned) < tolAngle)
                factorTwisting = 1.0;
            else {
                INFO_OUT() << "Found boundaries for which the normal vectors are not aligned with angle = " << phiNormals << endl;
                factorTwisting = 0.0;
            }

            // Check if the angle between the normal vectors is 0° (have the same direction) or 180° (have opposite directions) to formulate the interface constraint
            if (abs(phiNormals - 1) < sqrt(tolAngle))
                factorNormal = -1.0;
            else if (abs(phiNormals + 1) < sqrt(tolAngle))
                factorNormal = 1.0;
            else
                assert(false);

            // Compute the dual product matrices for the displacements
            EMPIRE::MathLibrary::computeTransposeMatrixProduct(noCoord,noDOFsLocMaster,noDOFsLocMaster,BDisplacementsGCMaster,BDisplacementsGCMaster,KPenaltyDisplacementMaster);
            EMPIRE::MathLibrary::computeTransposeMatrixProduct(noCoord,noDOFsLocSlave,noDOFsLocSlave,BDisplacementsGCSlave,BDisplacementsGCSlave,KPenaltyDisplacementSlave);
            EMPIRE::MathLibrary::computeTransposeMatrixProduct(noCoord,noDOFsLocMaster,noDOFsLocSlave,BDisplacementsGCMaster,BDisplacementsGCSlave,CPenaltyDisplacement);

            // Compute the dual product matrices for the bending rotations
            EMPIRE::MathLibrary::computeTransposeMatrixProduct(1, noDOFsLocMaster, noDOFsLocMaster ,BOperatorOmegaTMaster, BOperatorOmegaTMaster, KPenaltyBendingRotationMaster);
            EMPIRE::MathLibrary::computeTransposeMatrixProduct(1, noDOFsLocSlave, noDOFsLocSlave, BOperatorOmegaTSlave, BOperatorOmegaTSlave, KPenaltyBendingRotationSlave);
            EMPIRE::MathLibrary::computeTransposeMatrixProduct(1, noDOFsLocMaster, noDOFsLocSlave, BOperatorOmegaTMaster, BOperatorOmegaTSlave,CPenaltyBendingRotation);

            // Compute the dual product matrices for the twisting rotations
            EMPIRE::MathLibrary::computeTransposeMatrixProduct(1, noDOFsLocMaster, noDOFsLocMaster ,BOperatorOmegaNMaster, BOperatorOmegaNMaster, KPenaltyTwistingRotationMaster);
            EMPIRE::MathLibrary::computeTransposeMatrixProduct(1, noDOFsLocSlave, noDOFsLocSlave, BOperatorOmegaNSlave, BOperatorOmegaNSlave, KPenaltyTwistingRotationSlave);
            EMPIRE::MathLibrary::computeTransposeMatrixProduct(1, noDOFsLocMaster, noDOFsLocSlave, BOperatorOmegaNMaster, BOperatorOmegaNSlave,CPenaltyTwistingRotation);

            // Compute the element index tables for the master and slave patch
            int CPIndexMaster[noLocalBasisFctsMaster];
            int CPIndexSlave[noLocalBasisFctsSlave];
            patchMaster->getIGABasis()->getBasisFunctionsIndex(uKnotSpanMaster, vKnotSpanMaster, CPIndexMaster);
            patchSlave->getIGABasis()->getBasisFunctionsIndex(uKnotSpanSlave, vKnotSpanSlave, CPIndexSlave);

            // Compute the element freedom tables for the master and the slave patch
            int EFTMaster[noDOFsLocMaster];
            counter = 0;
            for (int i = 0; i < noLocalBasisFctsMaster ; i++){
                indexCP = patchMaster->getControlPointNet()[CPIndexMaster[i]]->getDofIndex();
                for (int j = 0; j < noCoord; j++){
                    EFTMaster[counter] = noCoord*indexCP + j;
                    counter++;
                }
            }

            counter = 0;
            int EFTSlave[noDOFsLocSlave];
            for (int i = 0; i < noLocalBasisFctsSlave ; i++){
                indexCP = patchSlave->getControlPointNet()[CPIndexSlave[i]]->getDofIndex();
                for (int j = 0; j < noCoord; j++){
                    EFTSlave[counter] = noCoord*indexCP + j;
                    counter++;
                }
            }

            // Assemble KPenaltyDisplacementMaster to the global coupling matrix CNN
            for(int i = 0; i < noDOFsLocMaster; i++){
                for(int j = 0; j < noDOFsLocMaster; j++){
                    // Assemble the displacement coupling entries
                    couplingMatrices->addCNNValue(EFTMaster[i], EFTMaster[j], alphaPrimary*KPenaltyDisplacementMaster[i*noDOFsLocMaster + j]*elementLengthOnGP);

                    // Assemble the bending rotation coupling entries
                    couplingMatrices->addCNNValue(EFTMaster[i], EFTMaster[j], alphaSecondaryBending*KPenaltyBendingRotationMaster[i*noDOFsLocMaster + j]*elementLengthOnGP);

                    // Assemble the twisting rotation coupling entries
                    couplingMatrices->addCNNValue(EFTMaster[i], EFTMaster[j], factorTwisting*alphaSecondaryTwisting*KPenaltyTwistingRotationMaster[i*noDOFsLocMaster + j]*elementLengthOnGP);
                }
            }

            // Assemble KPenaltyDisplacementSlave to the global coupling matrix CNN
            for(int i = 0; i < noDOFsLocSlave; i++){
                for(int j = 0; j < noDOFsLocSlave; j++) {
                    // Assemble the displacement coupling entries
                    couplingMatrices->addCNNValue(EFTSlave[i], EFTSlave[j], alphaPrimary*KPenaltyDisplacementSlave[i*noDOFsLocSlave + j]*elementLengthOnGP);

                    // Assemble the bending rotation coupling entries
                    couplingMatrices->addCNNValue(EFTSlave[i], EFTSlave[j], alphaSecondaryBending*KPenaltyBendingRotationSlave[i*noDOFsLocSlave + j]*elementLengthOnGP);

                    // Assemble the twisting rotation coupling entries
                    couplingMatrices->addCNNValue(EFTSlave[i], EFTSlave[j], factorTwisting*alphaSecondaryTwisting*KPenaltyTwistingRotationSlave[i*noDOFsLocSlave + j]*elementLengthOnGP);
                }
            }

            // Assemble CPenaltyDisplacement to the global coupling matrix CNN
            for(int i = 0; i < noDOFsLocMaster; i++){
                for(int j = 0; j < noDOFsLocSlave; j++){
                    // Assemble the displacement coupling entries
                    couplingMatrices->addCNNValue(EFTMaster[i], EFTSlave[j], alphaPrimary*(-1.0)*CPenaltyDisplacement[i*noDOFsLocSlave + j]*elementLengthOnGP);
                    couplingMatrices->addCNNValue(EFTSlave[j], EFTMaster[i], alphaPrimary*(-1.0)*CPenaltyDisplacement[i*noDOFsLocSlave + j]*elementLengthOnGP);

                    // Assemble the bending rotation coupling entries
                    couplingMatrices->addCNNValue(EFTMaster[i], EFTSlave[j], alphaSecondaryBending*factorTangent*CPenaltyBendingRotation[i*noDOFsLocSlave + j]*elementLengthOnGP);
                    couplingMatrices->addCNNValue(EFTSlave[j], EFTMaster[i], alphaSecondaryBending*factorTangent*CPenaltyBendingRotation[i*noDOFsLocSlave + j]*elementLengthOnGP);

                    // Assemble the twisting rotation coupling entries
                    couplingMatrices->addCNNValue(EFTMaster[i], EFTSlave[j], factorTwisting*alphaSecondaryTwisting*factorNormal*CPenaltyTwistingRotation[i*noDOFsLocSlave + j]*elementLengthOnGP);
                    couplingMatrices->addCNNValue(EFTSlave[j], EFTMaster[i], factorTwisting*alphaSecondaryTwisting*factorNormal*CPenaltyTwistingRotation[i*noDOFsLocSlave + j]*elementLengthOnGP);
                }
            }


            if(propErrorComputation.isInterfaceError){
                // Initialize variable storing the Gauss Point data
                std::vector<double> streamInterfaceGP;

                // elementLengthOnGP + noBasisFuncsI + (#indexCP, basisFuncValueI,...) + (#indexDOF, BtValueI, BnValueI,...) + noBasisFuncsJ + (#indexCP, basisFuncValueJ,...) + (#indexDOF, BtValueJ, BnValueJ,...) + factorTangent + factorNormal + factorTwisting
                streamInterfaceGP.reserve(1 + 1 + 2*noLocalBasisFctsMaster + 3*noDOFsLocMaster + 1 + 2*noLocalBasisFctsSlave + 3*noDOFsLocSlave + 1 + 1);

                // Save the element length on the Gauss Point
                streamInterfaceGP.push_back(elementLengthOnGP);

                // Save the number of basis functions of the master patch
                streamInterfaceGP.push_back(noLocalBasisFctsMaster);

                // Save the Control Point index and the basis function's values of the master patch
                for(int iBFs = 0; iBFs < noLocalBasisFctsMaster; iBFs++){
                    indexCP = patchMaster->getControlPointNet()[CPIndexMaster[iBFs]]->getDofIndex();
                    streamInterfaceGP.push_back(indexCP);
                    streamInterfaceGP.push_back(BDisplacementsGCMaster[0*noLocalBasisFctsMaster + 3*iBFs]);
                }

                // Save the DOF index and the bending and twisting B-operator values of the master patch
                for(int iDOFs = 0; iDOFs < noDOFsLocMaster; iDOFs++){
                    streamInterfaceGP.push_back(EFTMaster[iDOFs]);
                    streamInterfaceGP.push_back(BOperatorOmegaTMaster[iDOFs]);
                    streamInterfaceGP.push_back(BOperatorOmegaNMaster[iDOFs]);
                }

                // Save the number of basis functions of the slave patch
                streamInterfaceGP.push_back(noLocalBasisFctsSlave);

                // Save the Control Point index and the basis function's values of the slave patch
                for(int iBFs = 0; iBFs < noLocalBasisFctsSlave; iBFs++){
                    indexCP = patchSlave->getControlPointNet()[CPIndexSlave[iBFs]]->getDofIndex();
                    streamInterfaceGP.push_back(indexCP);
                    streamInterfaceGP.push_back(BDisplacementsGCSlave[0*noLocalBasisFctsSlave + 3*iBFs]);
                }

                // Save the DOF index and the bending and twisting B-operator values of the slave patch
                for(int iDOFs = 0; iDOFs < noDOFsLocSlave; iDOFs++){
                    streamInterfaceGP.push_back(EFTSlave[iDOFs]);
                    streamInterfaceGP.push_back(BOperatorOmegaTSlave[iDOFs]);
                    streamInterfaceGP.push_back(BOperatorOmegaNSlave[iDOFs]);
                }

                // Save the factors
                streamInterfaceGP.push_back(factorTangent);
                streamInterfaceGP.push_back(factorNormal);
                streamInterfaceGP.push_back(factorTwisting);

                // Push back the Gauss Point values into the member variable
                streamInterfaceGPs.push_back(streamInterfaceGP);
            }
        } // End of Gauss Point loop

        // Delete pointers
        delete[] BDisplacementsGCMaster;
        delete[] BDisplacementsGCSlave;
        delete[] BOperatorOmegaTMaster;
        delete[] BOperatorOmegaTSlave;
        delete[] BOperatorOmegaNMaster;
        delete[] BOperatorOmegaNSlave;
        delete[] KPenaltyDisplacementMaster;
        delete[] KPenaltyDisplacementSlave;
        delete[] CPenaltyDisplacement;
        delete[] KPenaltyBendingRotationMaster;
        delete[] KPenaltyBendingRotationSlave;
        delete[] CPenaltyBendingRotation;
        delete[] KPenaltyTwistingRotationMaster;
        delete[] KPenaltyTwistingRotationSlave;
        delete[] CPenaltyTwistingRotation;
    } // End of weak continuity condition loop
}

void IGAMortarMapper::computeDisplacementAndRotationBOperatorMatrices(double* _BDisplacementsGC, double* _BOperatorOmegaT,
                                                                      double* _BOperatorOmegaN, double* _normalTrCurveVct,
                                                                      IGAPatchSurface* _patch, double* _tangentTrCurveVct,
                                                                      double _u, double _v, int _uKnotSpan, int _vKnotSpan) {
    /*
     * Returns the B-operator matrix for the displacement field, the bending and the twisting rotations in the global Cartesian space.
     * Additionally the normal to the boundary vector which is tangent to the surface patch is returned.
     *
     * Hints on the size of the [in/out] arrays:
     *
     * double* _BDisplacementsGC = new double[noCoord*noDOFsLoc]
     * double* _BOperatorOmegaT = new double[noDOFsLoc]
     * double* _BOperatorOmegaN = new double[noDOFsLoc]
     * double _normalTrCurveVct[noCoord]
     */

    // Initialize constant array sizes
    const int noCoord = 3;
    const int noParametricCoord = 2;

    // Initialize varying array sizes
    double surfaceNormalVct[noCoord];
    double tangentTrCurveVctCov[noCoord];
    double normalTrCurveVctCov[noCoord];
    double surfNormalVctAndDervs[3*noCoord];
    double covariantMetricTensor[4];
    double contravariantCurvatureTensor[4];
    double contravariantBaseVcts[6];
    double contravariantBaseVct[3];
    double dT3Cov2GC[noParametricCoord*noCoord];
    double BabTimesContravariantBaseVct[noParametricCoord*noCoord];
    int derivDegreeBasis = 2;
    int derivDegreeBaseVec = derivDegreeBasis - 1;
    int noBaseVec = 2;
    int indexBasis;
    int indexBasisDerivU;
    int indexBasisDerivV;

    // Get the polynomial orders of the basis
    int p = _patch->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree();
    int q = _patch->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree();

    // Get the number of the basis functions affecting the current knot span
    int noLocalBasisFcts = (p + 1)*(q + 1);

    // get the number of the local DOFs
    int noDOFsLoc = noCoord*noLocalBasisFcts;

    // Initialize pointers
    double* basisFctsAndDerivs = new double[(derivDegreeBasis + 1) * (derivDegreeBasis + 2)
            * noLocalBasisFcts / 2];
    double* baseVctsAndDerivs = new double[(derivDegreeBaseVec + 1)
            * (derivDegreeBaseVec + 2) * noCoord * noBaseVec / 2];
    double* BdDisplacementsdUGC = new double[noCoord*noDOFsLoc];
    double* BdDisplacementsdVGC = new double[noCoord*noDOFsLoc];
    double* commonBOperator1 = new double[noParametricCoord*noDOFsLoc];
    double* commonBOperator2Part1 = new double[noDOFsLoc];
    double* commonBOperator2Part2 = new double[noDOFsLoc];
    double* commonBOperator2 = new double[noParametricCoord*noDOFsLoc];
    double* commonBOperator3 = new double[noParametricCoord*noDOFsLoc];
    double* commonBOperator = new double[noParametricCoord*noDOFsLoc];

    // Compute the basis functions and their derivatives
    _patch->getIGABasis()->computeLocalBasisFunctionsAndDerivatives(basisFctsAndDerivs, derivDegreeBasis, _u, _uKnotSpan,
                                                                    _v, _vKnotSpan);

    // Compute the base vectors and their derivatives
    _patch->computeBaseVectorsAndDerivatives(baseVctsAndDerivs, basisFctsAndDerivs, derivDegreeBaseVec,
                                             _uKnotSpan, _vKnotSpan);

    // Compute the derivatives of the surface normal vector
    _patch->computeSurfaceNormalVectorAndDerivatives(surfNormalVctAndDervs, baseVctsAndDerivs, derivDegreeBaseVec);

    // Compute the normal to the boundary vector which is tangent to the surface
    for(int i = 0; i < noCoord; i++){
        surfaceNormalVct[i] = surfNormalVctAndDervs[i];
    }
    EMPIRE::MathLibrary::computeVectorCrossProduct(surfaceNormalVct, _tangentTrCurveVct, _normalTrCurveVct);

    // Compute the covariant metric tensor
    _patch->computeCovariantMetricTensor(covariantMetricTensor, baseVctsAndDerivs, derivDegreeBaseVec);

    // Compute the contravariant base vectors of the master patch
    _patch->computeContravariantBaseVectors(contravariantBaseVcts, covariantMetricTensor, baseVctsAndDerivs, derivDegreeBaseVec);

    // Initialize the B-operator matrix for the displacement field
    for(int iBF = 0; iBF < noCoord*noDOFsLoc; iBF++){
        _BDisplacementsGC[iBF] = 0.0;
        BdDisplacementsdUGC[iBF] = 0.0;
        BdDisplacementsdVGC[iBF] = 0.0;
    }

    // Compute the B-operator matrix for the displacement field
    for(int iBF = 0; iBF < noLocalBasisFcts; iBF++){
        indexBasis = _patch->getIGABasis()->indexDerivativeBasisFunction(derivDegreeBasis,0,0,iBF);
        _BDisplacementsGC[0*noDOFsLoc + 3*iBF + 0] = basisFctsAndDerivs[indexBasis];
        _BDisplacementsGC[1*noDOFsLoc + 3*iBF + 1] = basisFctsAndDerivs[indexBasis];
        _BDisplacementsGC[2*noDOFsLoc + 3*iBF + 2] = basisFctsAndDerivs[indexBasis];
        indexBasisDerivU = _patch->getIGABasis()->indexDerivativeBasisFunction(derivDegreeBasis,1,0,iBF);
        BdDisplacementsdUGC[0*noDOFsLoc + 3*iBF + 0] = basisFctsAndDerivs[indexBasisDerivU];
        BdDisplacementsdUGC[1*noDOFsLoc + 3*iBF + 1] = basisFctsAndDerivs[indexBasisDerivU];
        BdDisplacementsdUGC[2*noDOFsLoc + 3*iBF + 2] = basisFctsAndDerivs[indexBasisDerivU];
        indexBasisDerivV = _patch->getIGABasis()->indexDerivativeBasisFunction(derivDegreeBasis,0,1,iBF);
        BdDisplacementsdVGC[0*noDOFsLoc + 3*iBF + 0] = basisFctsAndDerivs[indexBasisDerivV];
        BdDisplacementsdVGC[1*noDOFsLoc + 3*iBF + 1] = basisFctsAndDerivs[indexBasisDerivV];
        BdDisplacementsdVGC[2*noDOFsLoc + 3*iBF + 2] = basisFctsAndDerivs[indexBasisDerivV];
    }

    // Transform the normal and the tangent vectors to the covariant base for the master patch
    for(int iCov = 0; iCov < noParametricCoord; iCov++){
        for (int iCoord = 0; iCoord < noCoord; iCoord++)
            contravariantBaseVct[iCoord] = contravariantBaseVcts[noCoord*iCov + iCoord];
        tangentTrCurveVctCov[iCov] = EMPIRE::MathLibrary::computeDenseDotProduct(noCoord,contravariantBaseVct,_tangentTrCurveVct);
        normalTrCurveVctCov[iCov] = EMPIRE::MathLibrary::computeDenseDotProduct(noCoord,contravariantBaseVct,_normalTrCurveVct);
    }

    // Compute the curvature tensor in the contravariant basis
    _patch->computeContravariantCurvatureTensor(contravariantCurvatureTensor, surfaceNormalVct, baseVctsAndDerivs, derivDegreeBaseVec);

    // Compute the B-operator matrices for the rotation vector
    for(int iCovCoord = 0; iCovCoord < noParametricCoord; iCovCoord++)
        for(int iCoord = 0; iCoord < noCoord; iCoord++)
            dT3Cov2GC[iCovCoord*noCoord + iCoord] = surfNormalVctAndDervs[noCoord*(iCovCoord + 1) + iCoord]; // Parametric derivatives of the transformation matrix from the covariant to the global Cartesian basis
    EMPIRE::MathLibrary::computeMatrixProduct(noParametricCoord, noCoord, noDOFsLoc, dT3Cov2GC, _BDisplacementsGC, commonBOperator1);
    EMPIRE::MathLibrary::computeMatrixProduct(1, noCoord, noDOFsLoc, surfaceNormalVct, BdDisplacementsdUGC, commonBOperator2Part1);
    EMPIRE::MathLibrary::computeMatrixProduct(1, noCoord, noDOFsLoc, surfaceNormalVct, BdDisplacementsdVGC, commonBOperator2Part2);
    for (int iDOFs = 0; iDOFs < noDOFsLoc; iDOFs++){
        commonBOperator2[0*noDOFsLoc + iDOFs] = commonBOperator2Part1[iDOFs];
        commonBOperator2[1*noDOFsLoc + iDOFs] = commonBOperator2Part2[iDOFs];
    }
    EMPIRE::MathLibrary::computeTransposeMatrixProduct(noParametricCoord, noParametricCoord, noCoord, contravariantCurvatureTensor, contravariantBaseVcts, BabTimesContravariantBaseVct);
    EMPIRE::MathLibrary::computeMatrixProduct(noParametricCoord, noCoord, noDOFsLoc, BabTimesContravariantBaseVct, _BDisplacementsGC, commonBOperator3);
    for(int i = 0; i < noParametricCoord*noDOFsLoc; i++)
        commonBOperator[i] = commonBOperator1[i] + commonBOperator2[i] + commonBOperator3[i];
    EMPIRE::MathLibrary::computeDenseVectorMultiplicationScalar(normalTrCurveVctCov, -1.0, noCoord);
    EMPIRE::MathLibrary::computeTransposeMatrixProduct(noParametricCoord, 1, noDOFsLoc, normalTrCurveVctCov, commonBOperator, _BOperatorOmegaT);
    EMPIRE::MathLibrary::computeTransposeMatrixProduct(noParametricCoord, 1, noDOFsLoc, tangentTrCurveVctCov, commonBOperator, _BOperatorOmegaN);

    // Delete pointers
    delete[] basisFctsAndDerivs;
    delete[] baseVctsAndDerivs;
    delete[] BdDisplacementsdUGC;
    delete[] BdDisplacementsdVGC;
    delete[] commonBOperator1;
    delete[] commonBOperator2Part1;
    delete[] commonBOperator2Part2;
    delete[] commonBOperator2;
    delete[] commonBOperator3;
    delete[] commonBOperator;
}

void IGAMortarMapper::computePenaltyParametersForWeakDirichletCurveConditions(std::string _filename){
    /*
     * Compute the penalty factors related to the application of weak Dirichlet curve conditions
     * as a function of the minimum element edge size across interface.
     */

    // Initialize a file stream to write out the Penalty parameters
    std::ofstream ofs;
    if (!_filename.empty()) {
        ofs.open(_filename.c_str(), std::ofstream::out);
        ofs << std::scientific;
        ofs << std::endl;
    }

    // Get the weak patch continuity conditions
    std::vector<WeakIGADirichletCurveCondition*> weakIGADirichletCurveConditions = meshIGA->getWeakIGADirichletCurveConditions();

    // Initialize constant array sizes
    const int noCoord = 3;

    // Initialize varying array sizes
    bool isElementChanged;
    int patchIndex;
    int p;
    int q;
    int pMax;
    int noLocalBasisFcts;
    int noDOFsLoc;
    int noGPsOnCond;
    int uKnotSpan;
    int vKnotSpan;
    int uKnotSpanSaved;
    int vKnotSpanSaved;
    double uGP;
    double vGP;
    double elEdgeSize;
    double minElEdgeSize;
    double alphaPrim;
    double alphaSec;

    // Initialize pointers
    double* curveGPs;
    double* curveGPWeights;
    double* curveGPTangents;
    double* curveGPJacobianProducts;
    IGAPatchSurface* thePatch;

    // Loop over all the conditions for the application of weak continuity across patch interfaces
    for (int iWDCC = 0; iWDCC < weakIGADirichletCurveConditions.size(); iWDCC++){
        // Check if penalty factors are to be assigned manually in the xml file
        if(!propWeakCurveDirichletConditions.isAutomaticPenaltyParameters){
            weakDirichletCCAlphaPrimary[iWDCC] = propWeakCurveDirichletConditions.alphaPrim;
            weakDirichletCCAlphaSecondaryBending[iWDCC] = propWeakCurveDirichletConditions.alphaSecBending;
            weakDirichletCCAlphaSecondaryTwisting[iWDCC] = propWeakCurveDirichletConditions.alphaSecTwisting;
            continue;
        }

        // Get the index of the patch
        patchIndex = weakIGADirichletCurveConditions[iWDCC]->getPatchIndex();

        // Get the number of Gauss Points for the given condition
        noGPsOnCond = weakIGADirichletCurveConditions[iWDCC]->getCurveNumGP();

        // Get the parametric coordinates of the Gauss Points
        curveGPs = weakIGADirichletCurveConditions[iWDCC]->getCurveGPs();

        // Get the corresponding Gauss weights
        curveGPWeights = weakIGADirichletCurveConditions[iWDCC]->getCurveGPWeights();

        // Get the tangent vectors at the curve of the given condition in the Cartesian space
        curveGPTangents = weakIGADirichletCurveConditions[iWDCC]->getCurveGPTangents();

        // Get the product of the Jacobian transformations
        curveGPJacobianProducts = weakIGADirichletCurveConditions[iWDCC]->getCurveGPJacobianProducts();

        // Get the master and the slave patch
        thePatch = meshIGA->getSurfacePatch(patchIndex);

        // Get the polynomial orders of the patch
        p = thePatch->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree();
        q = thePatch->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree();
        pMax = std::max(p, q);

        // get the number of local basis functions for the patch
        noLocalBasisFcts = (p + 1)*(q + 1);

        // get the number of the local DOFs for the patch
        noDOFsLoc = noCoord*noLocalBasisFcts;

        // Initialize flags on whether an element has been changed while looping over the Gauss Points
        isElementChanged = false;

        // Initialize the element edge sizes
        elEdgeSize = 0.0;
        minElEdgeSize = std::numeric_limits<double>::max();

        // Loop over all the Gauss Points of the given condition
        for(int iGP = 0; iGP < noGPsOnCond; iGP++){

            // Get the parametric coordinates of the Gauss Point on the patch
            uGP = curveGPs[2*iGP];
            vGP = curveGPs[2*iGP + 1];

            // Find the knot span indices of the Gauss point locations in the parameter space of the patch
            uKnotSpan = thePatch->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(uGP);
            vKnotSpan = thePatch->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(vGP);

            // Initialize the saved knot span indices
            if(iGP == 0){
                uKnotSpanSaved = uKnotSpan;
                vKnotSpanSaved = vKnotSpan;
            }

            // Initialize element edge sizes if an element has been crossed
            if(uKnotSpan != uKnotSpanSaved || vKnotSpan != vKnotSpanSaved){
                if(elEdgeSize < minElEdgeSize)
                    minElEdgeSize = elEdgeSize;
                elEdgeSize = 0.0;
            }

            // Add the contribution from the Gauss Point to the element edge sizes
            elEdgeSize += curveGPJacobianProducts[iGP];

            // Save the knot span indices
            uKnotSpanSaved = uKnotSpan;
            vKnotSpanSaved = vKnotSpan;

        } // End of Gauss Point loop

        // Check the element sizes for the last elements
        if(elEdgeSize < minElEdgeSize)
            minElEdgeSize = elEdgeSize;

        // Compute correspondingly the penalty factors
        alphaPrim = pMax/minElEdgeSize;
        alphaSec = pMax/sqrt(minElEdgeSize);
        weakDirichletCCAlphaPrimary[iWDCC] = alphaPrim;
        weakDirichletCCAlphaSecondaryBending[iWDCC] = alphaSec;
        weakDirichletCCAlphaSecondaryTwisting[iWDCC] = alphaSec;

        // Write the Penalty parameters into a file
        if (!_filename.empty()) {
            ofs << "Weak Dirichlet curve condition on patch[" << patchIndex << "] :" << std::endl;
            ofs << "weakDirichletCCAlphaPrimary[" << iWDCC << "] = " << scientific << setprecision(15) << alphaPrim << std::endl;
            ofs << "weakDirichletCCAlphaSecondary[" << iWDCC << "] = " << scientific << setprecision(15) << alphaSec << std::endl;
            ofs << std::endl;
        }
    } // End of weak Dirichlet curve condition loop

    // Close file
    if (!_filename.empty())
        ofs.close();
}

void IGAMortarMapper::computePenaltyParametersForWeakDirichletSurfaceConditions(){

    /*
     * Compute the penalty factors related to the application of weak Dirichlet surface conditions
     */
    ERROR_OUT() << "Function under construction" << endl;
    exit(-1);

    std::vector<WeakIGADirichletSurfaceCondition*> weakIGADirichletSurfaceConditions = meshIGA->getWeakIGADirichletSurfaceConditions();

    int patchIndex;

    for (int iWDSC = 0; iWDSC < weakIGADirichletSurfaceConditions.size(); iWDSC++){
        weakDirichletSCAlphaPrimary[iWDSC] = propWeakSurfaceDirichletConditions.alphaPrim;

        patchIndex = weakIGADirichletSurfaceConditions[iWDSC]->getPatchIndex();

        if (Message::isDebugMode()){
            DEBUG_OUT() << std::endl;
            DEBUG_OUT() << "Weak Dirichlet surface condition on patch[" << patchIndex << "] :" << std::endl;
            DEBUG_OUT() << "weakDirichletSCAlphaPrimary[" << iWDSC << "] = " << scientific << setprecision(15) << weakDirichletSCAlphaPrimary[iWDSC] << std::endl;
            DEBUG_OUT() << "weakDirichletSCAlphaSecondaryBending[" << iWDSC << "] = " << scientific << setprecision(15) << weakDirichletSCAlphaSecondaryBending[iWDSC] << std::endl;
            DEBUG_OUT() << "weakDirichletSCAlphaSecondaryTwisting[" << iWDSC << "] = " << scientific << setprecision(15) << weakDirichletSCAlphaSecondaryTwisting[iWDSC] << std::endl;
            DEBUG_OUT() << std::endl;
        }
    }

}

void IGAMortarMapper::computePenaltyParametersForPatchContinuityConditions(std::string _filename){
    /*
     * Compute the penalty factors related to the application of weak patch continuity conditions
     * as a function of the minimum element edge size across each interface.
     */

    // Initialize a file stream to write out the Penalty parameters
    std::ofstream ofs;
    if (!_filename.empty()) {
        ofs.open(_filename.c_str(), std::ofstream::out);
        ofs << std::scientific;
        ofs << std::endl;
    }

    // Get the weak patch continuity conditions
    std::vector<WeakIGAPatchContinuityCondition*> weakIGAPatchContinuityConditions = meshIGA->getWeakIGAPatchContinuityConditions();

    // Initialize constant array sizes
    const int noCoord = 3;

    // Initialize varying array sizes
    bool isElementMasterChanged;
    bool isElementSlaveChanged;
    int indexMaster;
    int indexSlave;
    int pMaster;
    int qMaster;
    int pSlave;
    int qSlave;
    int pMaxMaster;
    int pMaxSlave;
    int pMax;
    int noLocalBasisFctsMaster;
    int noLocalBasisFctsSlave;
    int noDOFsLocMaster;
    int noDOFsLocSlave;
    int noGPsOnContCond;
    int uKnotSpanMaster;
    int uKnotSpanSlave;
    int vKnotSpanMaster;
    int vKnotSpanSlave;
    int uKnotSpanMasterSaved;
    int uKnotSpanSlaveSaved;
    int vKnotSpanMasterSaved;
    int vKnotSpanSlaveSaved;
    double uGPMaster;
    double vGPMaster;
    double uGPSlave;
    double vGPSlave;
    double elEdgeSizeMaster;
    double elEdgeSizeSlave;
    double minElEdgeSizeMaster;
    double minElEdgeSizeSlave;
    double minElEdgeSize;
    double alphaSec;
    double alphaPrim;

    // Initialize pointers
    double* trCurveMasterGPs;
    double* trCurveSlaveGPs;
    double* trCurveGPWeights;
    double* trCurveMasterGPTangents;
    double* trCurveSlaveGPTangents;
    double* trCurveGPJacobianProducts;
    IGAPatchSurface* patchMaster;
    IGAPatchSurface* patchSlave;

    // Loop over all the conditions for the application of weak continuity across patch interfaces
    for (int iWCC = 0; iWCC < weakIGAPatchContinuityConditions.size(); iWCC++){
        // Check if penalty factors are to be assigned manually in the xml file
        if(!propWeakPatchContinuityConditions.isAutomaticPenaltyParameters){
            weakPatchContinuityAlphaPrimaryIJ[iWCC] = propWeakPatchContinuityConditions.alphaPrim;
            weakPatchContinuityAlphaSecondaryBendingIJ[iWCC] = propWeakPatchContinuityConditions.alphaSecBending;
            weakPatchContinuityAlphaSecondaryTwistingIJ[iWCC] = propWeakPatchContinuityConditions.alphaSecTwisting;
            continue;
        }

        // Get the index of the master and slave patches
        indexMaster = weakIGAPatchContinuityConditions[iWCC]->getMasterPatchIndex();
        indexSlave = weakIGAPatchContinuityConditions[iWCC]->getSlavePatchIndex();

        // Get the number of Gauss Points for the given condition
        noGPsOnContCond = weakIGAPatchContinuityConditions[iWCC]->getTrCurveNumGP();

        // Get the parametric coordinates of the Gauss Points
        trCurveMasterGPs = weakIGAPatchContinuityConditions[iWCC]->getTrCurveMasterGPs();
        trCurveSlaveGPs = weakIGAPatchContinuityConditions[iWCC]->getTrCurveSlaveGPs();

        // Get the corresponding Gauss weights
        trCurveGPWeights = weakIGAPatchContinuityConditions[iWCC]->getTrCurveGPWeights();

        // Get the tangent vectors at the trimming curve of the given condition in the Cartesian space
        trCurveMasterGPTangents = weakIGAPatchContinuityConditions[iWCC]->getTrCurveMasterGPTangents();
        trCurveSlaveGPTangents = weakIGAPatchContinuityConditions[iWCC]->getTrCurveSlaveGPTangents();

        // Get the product of the Jacobian transformations
        trCurveGPJacobianProducts = weakIGAPatchContinuityConditions[iWCC]->getTrCurveGPJacobianProducts();

        // Get the master and the slave patch
        patchMaster = meshIGA->getSurfacePatch(indexMaster);
        patchSlave = meshIGA->getSurfacePatch(indexSlave);

        // Get the polynomial orders of the master and the slave patch
        pMaster = patchMaster->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree();
        qMaster = patchMaster->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree();
        pSlave = patchSlave->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree();
        qSlave = patchSlave->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree();
        pMaxMaster = pMaster;
        if(pMaxMaster < qMaster)
            pMaxMaster = qMaster;
        pMaxSlave = pSlave;
        if(pMaxSlave < qSlave)
            pMaxSlave = qSlave;
        pMax = pMaxMaster;
        if(pMaxMaster < pMaxSlave)
            pMax = pMaxSlave;

        // get the number of local basis functions for master and slave patch
        noLocalBasisFctsMaster = (pMaster + 1)*(qMaster + 1);
        noLocalBasisFctsSlave = (pSlave + 1)*(qSlave + 1);

        // get the number of the local DOFs for the master and slave patch
        noDOFsLocMaster = noCoord*noLocalBasisFctsMaster;
        noDOFsLocSlave = noCoord*noLocalBasisFctsSlave;

        // Initialize flags on whether an element has been changed while looping over the Gauss Points
        isElementMasterChanged = false;
        isElementSlaveChanged = false;

        // Initialize the element edge sizes
        elEdgeSizeMaster = 0.0;
        elEdgeSizeSlave = 0.0;
        minElEdgeSizeMaster = std::numeric_limits<double>::max();
        minElEdgeSizeSlave = std::numeric_limits<double>::max();

        // Loop over all the Gauss Points of the given condition
        for(int iGP = 0; iGP < noGPsOnContCond; iGP++){

            // Get the parametric coordinates of the Gauss Point on the master patch
            uGPMaster = trCurveMasterGPs[2*iGP];
            vGPMaster = trCurveMasterGPs[2*iGP + 1];

            // Get the parametric coordinates of the Gauss Point on the slave patch
            uGPSlave = trCurveSlaveGPs[2*iGP];
            vGPSlave = trCurveSlaveGPs[2*iGP + 1];

            // Find the knot span indices of the Gauss point locations in the parameter space of the master patch
            uKnotSpanMaster = patchMaster->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(uGPMaster);
            vKnotSpanMaster = patchMaster->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(vGPMaster);

            // Find the knot span indices of the Gauss point locations in the parameter space of the slave patch
            uKnotSpanSlave = patchSlave->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(uGPSlave);
            vKnotSpanSlave = patchSlave->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(vGPSlave);

            // Initialize the saved knot span indices
            if(iGP == 0){
                uKnotSpanMasterSaved = uKnotSpanMaster;
                vKnotSpanMasterSaved = vKnotSpanMaster;
                uKnotSpanSlaveSaved = uKnotSpanSlave;
                vKnotSpanSlaveSaved = vKnotSpanSlave;
            }

            // Initialize element edge sizes if an element has been crossed
            if(uKnotSpanMaster != uKnotSpanMasterSaved || vKnotSpanMaster != vKnotSpanMasterSaved){
                if(elEdgeSizeMaster < minElEdgeSizeMaster)
                    minElEdgeSizeMaster = elEdgeSizeMaster;
                elEdgeSizeMaster = 0.0;
            }
            if(uKnotSpanSlave != uKnotSpanSlaveSaved || vKnotSpanSlave != vKnotSpanSlaveSaved){
                if(elEdgeSizeSlave < minElEdgeSizeSlave)
                    minElEdgeSizeSlave = elEdgeSizeSlave;
                elEdgeSizeSlave = 0.0;
            }

            // Add the contribution from the Gauss Point to the element edge sizes
            elEdgeSizeMaster += trCurveGPJacobianProducts[iGP];
            elEdgeSizeSlave += trCurveGPJacobianProducts[iGP];

            // Save the knot span indices
            uKnotSpanMasterSaved = uKnotSpanMaster;
            vKnotSpanMasterSaved = vKnotSpanMaster;
            uKnotSpanSlaveSaved = uKnotSpanSlave;
            vKnotSpanSlaveSaved = vKnotSpanSlave;
        } // End of Gauss Point loop

        // Check the element sizes for the last elements
        if(elEdgeSizeMaster < minElEdgeSizeMaster)
            minElEdgeSizeMaster = elEdgeSizeMaster;
        if(elEdgeSizeSlave < minElEdgeSizeSlave)
            minElEdgeSizeSlave = elEdgeSizeSlave;

        // Get the minimum of the minimum element edge sizes between both patches
        minElEdgeSize = minElEdgeSizeMaster;
        if(minElEdgeSizeSlave < minElEdgeSize){
            minElEdgeSize = minElEdgeSizeSlave;
        }

        // Compute correspondingly the penalty factors
        alphaPrim = pMax/minElEdgeSize;
        alphaSec = pMax/sqrt(minElEdgeSize);
        weakPatchContinuityAlphaPrimaryIJ[iWCC] = alphaPrim;
        weakPatchContinuityAlphaSecondaryBendingIJ[iWCC] = alphaSec;
        weakPatchContinuityAlphaSecondaryTwistingIJ[iWCC] = alphaSec;

        // Write the Penalty parameters into a file
        if (!_filename.empty()) {
            ofs << "Coupling between patch[" << indexMaster << "] and patch[" << indexSlave << "]:" << std::endl;
            ofs << "weakPatchContinuityAlphaPrimaryIJ[" << iWCC << "] = " << alphaPrim << std::endl;
            ofs << "weakPatchContinuityAlphaSecondaryIJ[" << iWCC << "] = " << alphaSec << std::endl;
            ofs << std::endl;
        }
    } // End of weak continuity condition loop

    // Close file
    if (!_filename.empty())
        ofs.close();
}

void IGAMortarMapper::consistentMapping(const double* _slaveField, double *_masterField) {
    /*
     * Mapping of the
     * Cnn * x_master = Cnr * x_slave
     */

    // Get the appropriate sizes of the matrices
    int size_N = couplingMatrices->getSizeN();

    // Initialize right hand side vector
    double* tmpVec = new double[size_N]();

    // C_NR * x_slave = tmpVec
    couplingMatrices->getCnr()->mulitplyVec(false,const_cast<double *>(_slaveField), tmpVec, size_N);

    // Solve for the master field Cnn * x_master = tmpVec
    couplingMatrices->getCnn()->solve(_masterField, tmpVec);

    // Compute the error in the relative L2-norm
    if (propErrorComputation.isErrorComputation){
        double errorL2Domain;
        double errorL2Interface[2];
        double errorL2Curve[2];
        if(propErrorComputation.isDomainError)
            errorL2Domain = computeDomainErrorInL2Norm4ConsistentMapping(_slaveField, _masterField);
        if(propErrorComputation.isInterfaceError)
            if(isMappingIGA2FEM){
                computeIGAPatchInterfaceErrorInL2Norm(errorL2Interface, _slaveField);
            }else{
                computeIGAPatchInterfaceErrorInL2Norm(errorL2Interface, _masterField);
            }
        if(propErrorComputation.isCurveError)
            if(isMappingIGA2FEM){
                computeIGADirichletCurveErrorInL2Norm(errorL2Curve, _slaveField);
            }else{
                computeIGADirichletCurveErrorInL2Norm(errorL2Curve, _masterField);
            }
        printErrorMessage(infoOut, errorL2Domain, errorL2Curve, errorL2Interface);
    }

    // Delete pointers
    delete[] tmpVec;
}

void IGAMortarMapper::conservativeMapping(const double* _masterField, double *_slaveField) {
    /*
     * Mapping of the
     * f_slave = (Cnn^(-1) * Cnr)^T * f_master
     */
    int size_N = couplingMatrices->getSizeN();

    double* tmpVec = new double[size_N];

    couplingMatrices->getCnn()->solve(tmpVec, const_cast<double *>(_masterField));

    couplingMatrices->getCnr()->transposeMulitplyVec(tmpVec, _slaveField, numNodesMaster);

    delete[] tmpVec;
}

double IGAMortarMapper::computeDomainErrorInL2Norm4ConsistentMapping(const double *_slaveField, const double *_masterField){
    /*
     * Returns the relative error in the L2 norm in the domain using the isogeometric mortar-based mapping.
     *
     * The values of the basis functions and other consituents necessary for the integration are provided in the array streamGPs
     * in the following sequence,
     *
     * weight + jacobian + nShapeFuncsFE + (#dof, shapefuncValue,...) + nShapeFuncsIGA + (#dof, shapefuncValue,...)
     *
     */

    // Initialize output array
    double errorL2Domain = 0.0;

    // Initialize auxiliary arrays
    double slaveFieldL2Domain = 0.0;
    double slaveFieldNorm;
    double JacobianProducts;
    double GW;
    double basisFctFEM;
    double basisFctIGA;
    double errorGPSquare;
    double fieldFEM[3];
    double fieldIGA[3];
    double errorVct[3];
    int noNodesFE;
    int noCPsIGA;
    int indexNode;
    int indexDOF;
    int indexCP;
    int noCoord = 3;

    // Define tolerance
    double tolNormSlaveField = 1e-6;

    // Loop over all the Gauss Points
    for(int iGP = 0; iGP < streamGPs.size(); iGP++){
        // Get the Gauss Point Weight
        GW = streamGPs[iGP][0];

        // Get the product of the Jacobian transformations at the Gauss point
        JacobianProducts = streamGPs[iGP][1];

        // Get the number of the basis functions of the finite element
        noNodesFE = streamGPs[iGP][2];

        // Initialize the field on the finite element mesh at the Gauss point
        for(int iCoord = 0; iCoord < noCoord; iCoord++)
            fieldFEM[iCoord] = 0.0;

        // Loop over the nodes of the finite element
        for(int iNodesFE = 0; iNodesFE < noNodesFE; iNodesFE++){
            // Get the value of the basis function
            basisFctFEM = streamGPs[iGP][3 + 2*iNodesFE + 1];

            // Get the index of the node
            indexNode = streamGPs[iGP][3 + 2*iNodesFE];
            for(int iCoord = 0; iCoord < noCoord; iCoord++)
                if(!isMappingIGA2FEM)
                    fieldFEM[iCoord] += basisFctFEM*_slaveField[noCoord*indexNode + iCoord];
                else
                    fieldFEM[iCoord] += basisFctFEM*_masterField[noCoord*indexNode + iCoord];
        }

        // Get the number of basis functions of the isogeometric discretization
        noCPsIGA = streamGPs[iGP][3 + 2*noNodesFE];

        // Initialize the field on the isogeometric discretization at the Gauss point
        for(int iCoord = 0; iCoord < noCoord; iCoord++)
            fieldIGA[iCoord] = 0.0;

        // Loop over the Control Points of the isogeometric discretization
        for(int iCPsIGA = 0; iCPsIGA < noCPsIGA; iCPsIGA++){
            // Get the value of the basis function
            basisFctIGA = streamGPs[iGP][3 + 2*noNodesFE + 2*iCPsIGA + 2];

            // Get the index of the CP
            indexCP = streamGPs[iGP][3 + 2*noNodesFE + 2*iCPsIGA + 1];
            for(int iCoord = 0; iCoord < noCoord; iCoord++)
                if(!isMappingIGA2FEM)
                    fieldIGA[iCoord] += basisFctIGA*_masterField[noCoord*indexCP + iCoord];
                else
                    fieldIGA[iCoord] += basisFctIGA*_slaveField[noCoord*indexCP + iCoord];
        }

        // Compute the difference of the vectors on the Gauss Point
        for(int iCoord = 0; iCoord < noCoord; iCoord++)
            errorVct[iCoord] = fieldFEM[iCoord] - fieldIGA[iCoord];

        // Compute the norm of the difference of the fields at the Gauss Point
        errorGPSquare = EMPIRE::MathLibrary::computeDenseDotProduct(noCoord, errorVct, errorVct);

        // Compute the norm of the reference field at the Gauss Point
        if(!isMappingIGA2FEM)
            slaveFieldNorm = EMPIRE::MathLibrary::computeDenseDotProduct(noCoord, fieldFEM, fieldFEM);
        else
            slaveFieldNorm = EMPIRE::MathLibrary::computeDenseDotProduct(noCoord, fieldIGA, fieldIGA);

        // Add the contributions from the Gauss Point
        errorL2Domain += errorGPSquare*JacobianProducts*GW;
        slaveFieldL2Domain += slaveFieldNorm*JacobianProducts*GW;
    }

    // Compute the relative L2 norm of the mapping error
    errorL2Domain = sqrt(errorL2Domain);
    slaveFieldL2Domain = sqrt(slaveFieldL2Domain);
    if (slaveFieldL2Domain > tolNormSlaveField)
        errorL2Domain /= slaveFieldL2Domain;
    else
        WARNING_OUT() << "The norm of the slave field is smaller than the tolerance, no division of the mapping error is made" << std::endl;

    // Return the relative L2 norm of the mapping error
    return errorL2Domain;
}

void IGAMortarMapper::computeIGADirichletCurveErrorInL2Norm(double* _errorL2Curve, const double *_fieldIGA){
    /*
     * Returns the error in the L2 norm along the trimming curves where conditions are applied in terms of the primary and the
     * secondary fields in a double array of constant size 2
     *
     * The values of the basis functions and other consituents necessary for the integration are provided in the array streamGPs
     * in the following sequence,
     *
     * elementLengthOnGP + noBasisFuncs + (#indexCP, basisFuncValue,...) + (#indexDOF, BtValue, BnValue,...)
     */

    // Initialize output array
    for(int i = 0; i < 2; i++){
        _errorL2Curve[i] = 0.0;
    }

    // Initialize auxiliary arrays
    int noCPs;
    int noDOFs;
    int indexCP;
    int indexDOF;
    int noCoord = 3;
    double elementLengthOnGP;
    double basisFct;
    double BoperatorT;
    double BoperatorN;
    double field[3];
    double omegaT;
    double omegaN;
    double errorBendingRotation;
    double errorTwistingRotation;
    double normRotationSquare;
    double normErrorFieldSquare;
    double errorField[3];

    // Loop over all the Gauss Points
    for(int iGP = 0; iGP < streamCurveGPs.size(); iGP++){
        // Get the element length on the Gauss Point
        elementLengthOnGP = streamCurveGPs[iGP][0];

        // Get the number of basis functions of patch
        noCPs = streamCurveGPs[iGP][1];
        noDOFs = 3*noCPs;

        // Initialize the field and its rotations on the isogeometric discretization at the Gauss point
        for(int iCoord = 0; iCoord < noCoord; iCoord++)
            field[iCoord] = 0.0;
        omegaT = 0.0;
        omegaN = 0.0;

        // Loop over all basis functions of the patch and compute the field at the Gauss Point
        for(int iBFs = 0; iBFs < noCPs; iBFs++){
            // Get the index of the Control Point
            indexCP = streamCurveGPs[iGP][1 + 1 + 2*iBFs];

            // Get the value of the basis function
            basisFct = streamCurveGPs[iGP][1 + 1 + 2*iBFs + 1];

            // Loop over all the Cartesian coordinates
            for(int iCoord = 0; iCoord < noCoord; iCoord++)
                field[iCoord] += basisFct*_fieldIGA[noCoord*indexCP + iCoord];
        }

        // Loop over all the DOFs of the patch and compute the rotations of the field at the Gauss point
        for(int iDOFs = 0; iDOFs < noDOFs; iDOFs++){
            // Get the index of the DOF
            indexDOF = streamCurveGPs[iGP][1 + 1 + 2*noCPs + 3*iDOFs];

            // Get the value of the B-operator for the bending rotation
            BoperatorT = streamCurveGPs[iGP][1 + 1 + 2*noCPs + 3*iDOFs + 1];

            // Get the value of the B-operator for the twisting rotation
            BoperatorN = streamCurveGPs[iGP][1 + 1 + 2*noCPs + 3*iDOFs + 2];

            // Compute the tangent and the bending rotations
            omegaT += BoperatorT*_fieldIGA[indexDOF];
            omegaN += BoperatorN*_fieldIGA[indexDOF];
        }

        // Compute the error vector for the displacements
        for(int iCoord = 0; iCoord < noCoord; iCoord++){
            errorField[iCoord] = field[iCoord] - 0.0;
        }
        normErrorFieldSquare = EMPIRE::MathLibrary::computeDenseDotProduct(noCoord, errorField, errorField);
        _errorL2Curve[0] += normErrorFieldSquare*elementLengthOnGP;

        // Compute the error in terms of the rotations
        errorBendingRotation = omegaT + 0.0;
        errorTwistingRotation = omegaN + 0.0;
        normRotationSquare = errorBendingRotation*errorBendingRotation + errorTwistingRotation*errorTwistingRotation;
        _errorL2Curve[1] += normRotationSquare*elementLengthOnGP;
    }

    // Take the necessary for the norm square roots
    for(int i = 0; i < 2; i++)
        _errorL2Curve[i] = sqrt(_errorL2Curve[i]);
}

void IGAMortarMapper::computeIGAPatchInterfaceErrorInL2Norm(double* _errorL2Interface, const double *_fieldIGA){
    /*
     * Returns the error in the L2 norm across the patch interfaces in terms of the displacements and the rotations in a double array of constant size 2
     *
     * The values of the basis functions and other consituents necessary for the integration are provided in the array streamGPs
     * in the following sequence,
     *
     * elementLengthOnGP + noBasisFuncsI + (#indexCP, basisFuncValueI,...) + (#indexDOF, BtValueI, BnValueI,...) + noBasisFuncsJ + (#indexCP, basisFuncValueJ,...) + (#indexDOF, BtValueJ, BnValueJ,...)
     */

    // Initialize output array
    for(int i = 0; i < 2; i++){
        _errorL2Interface[i] = 0.0;
    }

    // Initialize auxiliary arrays
    int noCPsI;
    int noCPsJ;
    int noDOFsI;
    int noDOFsJ;
    int indexCP;
    int indexDOF;
    int noCoord = 3;
    double factorTangent;
    double factorNormal;
    double factorTwisting;
    double elementLengthOnGP;
    double basisFct;
    double BoperatorT;
    double BoperatorN;
    double fieldI[3];
    double fieldJ[3];
    double omegaTI;
    double omegaTJ;
    double omegaNI;
    double omegaNJ;
    double errorBendingRotation;
    double errorTwistingRotation;
    double normRotationSquare;
    double normErrorFieldSquare;
    double errorField[3];

    // Loop over all the Gauss Points
    for(int iGP = 0; iGP < streamInterfaceGPs.size(); iGP++){
        // Get the element length on the Gauss Point
        elementLengthOnGP = streamInterfaceGPs[iGP][0];

        // Get the number of basis functions of patch I
        noCPsI = streamInterfaceGPs[iGP][1];
        noDOFsI = 3*noCPsI;

        // Initialize the field and its rotations on the isogeometric discretization at the Gauss point
        for(int iCoord = 0; iCoord < noCoord; iCoord++){
            fieldI[iCoord] = 0.0;
            fieldJ[iCoord] = 0.0;
        }
        omegaTI = 0.0;
        omegaTJ = 0.0;
        omegaNI = 0.0;
        omegaNJ = 0.0;

        // Loop over all basis functions of patch I and compute the field at the Gauss Point
        for(int iBFs = 0; iBFs < noCPsI; iBFs++){
            // Get the index of the Control Point
            indexCP = streamInterfaceGPs[iGP][1 + 1 + 2*iBFs];

            // Get the value of the basis function
            basisFct = streamInterfaceGPs[iGP][1 + 1 + 2*iBFs + 1];

            // Loop over all the Cartesian coordinates
            for(int iCoord = 0; iCoord < noCoord; iCoord++)
                fieldI[iCoord] += basisFct*_fieldIGA[noCoord*indexCP + iCoord];
        }

        // Loop over all the DOFs of patch I and compute the rotations of the field at the Gauss point
        for(int iDOFs = 0; iDOFs < noDOFsI; iDOFs++){
            // Get the index of the DOF
            indexDOF = streamInterfaceGPs[iGP][1 + 1 + 2*noCPsI + 3*iDOFs];

            // Get the value of the B-operator for the bending rotation
            BoperatorT = streamInterfaceGPs[iGP][1 + 1 + 2*noCPsI + 3*iDOFs + 1];

            // Get the value of the B-operator for the twisting rotation
            BoperatorN = streamInterfaceGPs[iGP][1 + 1 + 2*noCPsI + 3*iDOFs + 2];

            // Compute the tangent and the bending rotations
            omegaTI += BoperatorT*_fieldIGA[indexDOF];
            omegaNI += BoperatorN*_fieldIGA[indexDOF];
        }

        // Get the number of basis functions of patch J
        noCPsJ = streamInterfaceGPs[iGP][1 + 1 + 2*noCPsI + 3*noDOFsI];
        noDOFsJ = 3*noCPsJ;

        // Loop over all basis functions of patch J and compute the field at the Gauss Point
        for(int iBFs = 0; iBFs < noCPsJ; iBFs++){
            // Get the index of the Control Point
            indexCP = streamInterfaceGPs[iGP][1 + 1 + 2*noCPsI + 3*noDOFsI + 1 + 2*iBFs];

            // Get the value of the basis function
            basisFct = streamInterfaceGPs[iGP][1 + 1 + 2*noCPsI + 3*noDOFsI + 1 + 2*iBFs + 1];

            // Loop over all the Cartesian coordinates
            for(int iCoord = 0; iCoord < noCoord; iCoord++)
                fieldJ[iCoord] += basisFct*_fieldIGA[noCoord*indexCP + iCoord];
        }

        // Loop over all the DOFs of patch J and compute the rotations of the field at the Gauss point
        for(int iDOFs = 0; iDOFs < noDOFsJ; iDOFs++){
            // Get the index of the DOF
            indexDOF = streamInterfaceGPs[iGP][1 + 1 + 2*noCPsI + 3*noDOFsI + 1 + 2*noCPsJ + 3*iDOFs];

            // Get the value of the B-operator for the bending rotation
            BoperatorT = streamInterfaceGPs[iGP][1 + 1 + 2*noCPsI + 3*noDOFsI + 1 + 2*noCPsJ + 3*iDOFs + 1];

            // Get the value of the B-operator for the twisting rotation
            BoperatorN = streamInterfaceGPs[iGP][1 + 1 + 2*noCPsI + 3*noDOFsI + 1 + 2*noCPsJ + 3*iDOFs + 2];

            // Compute the tangent and the bending rotations
            omegaTJ += BoperatorT*_fieldIGA[indexDOF];
            omegaNJ += BoperatorN*_fieldIGA[indexDOF];
        }

        // Get the factors for the bending and the twisting rotation
        factorTangent = streamInterfaceGPs[iGP][1 + 1 + 2*noCPsI + 3*noDOFsI + 1 + 2*noCPsJ + 3*noDOFsJ];
        factorNormal = streamInterfaceGPs[iGP][1 + 1 + 2*noCPsI + 3*noDOFsI + 1 + 2*noCPsJ + 3*noDOFsJ + 1];

        // Get the factor on whether the twisting rotations are to be used
        factorTwisting = streamInterfaceGPs[iGP][1 + 1 + 2*noCPsI + 3*noDOFsI + 1 + 2*noCPsJ + 3*noDOFsJ + 2];

        // Compute the error vector for the displacements
        for(int iCoord = 0; iCoord < noCoord; iCoord++){
            errorField[iCoord] = fieldI[iCoord] - fieldJ[iCoord];
        }
        normErrorFieldSquare = EMPIRE::MathLibrary::computeDenseDotProduct(noCoord, errorField, errorField);
        _errorL2Interface[0] += normErrorFieldSquare*elementLengthOnGP;

        // Compute the error in terms of the rotations
        errorBendingRotation = omegaTI + factorTangent*omegaTJ;
        errorTwistingRotation = factorTwisting*(omegaNI + factorNormal*omegaNJ);
        normRotationSquare = errorBendingRotation*errorBendingRotation + errorTwistingRotation*errorTwistingRotation;
        _errorL2Interface[1] += normRotationSquare*elementLengthOnGP;
    }

    // Take the necessary for the norm square roots
    for(int i = 0; i < 2; i++)
        _errorL2Interface[i] = sqrt(_errorL2Interface[i]);
}

void IGAMortarMapper::writeGaussPointData() {
    string filename = name + "_GaussPointData.csv";
    ofstream filestream;
    filestream.open(filename.c_str());
    filestream.precision(12);
    filestream << std::dec;
    for(std::vector<std::vector<double> >::iterator it1 = streamGPs.begin(); it1!=streamGPs.end(); it1++) {
        for(std::vector<double>::iterator it2=it1->begin(); it2!=it1->end(); it2++) {
            filestream << *it2 << " ";
        }
        filestream << endl;
    }
    filestream.close();
}

void IGAMortarMapper::writeProjectedNodesOntoIGAMesh() {
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

void IGAMortarMapper::writeParametricProjectedPolygons(string _filename) {
    string filename = name + "_" + _filename + ".csv";
    ofstream out;
    out.open(filename.c_str(), ofstream::out);
    for(int i=0;i<projectedPolygons.size();i++) {
        for(map<int,Polygon2D>::const_iterator it=projectedPolygons[i].begin(); it!=projectedPolygons[i].end(); it++) {
            out << i << "\t" << it->first;
            for(int j=0; j<it->second.size(); j++) {
                out<<"\t"<<(it->second)[j].first<<"\t"<<(it->second)[j].second;
            }
            out << endl;
        }
    }
    out.close();
}

void IGAMortarMapper::writeTriangulatedParametricPolygon(string _filename) {
    string filename = name + "_" + _filename + ".csv";
    ofstream out;
    out.open(filename.c_str(), ofstream::out);
    for(int i=0;i<triangulatedProjectedPolygons.size();i++) {
        for(map<int,ListPolygon2D>::const_iterator it1=triangulatedProjectedPolygons[i].begin(); it1!=triangulatedProjectedPolygons[i].end(); it1++) {
            out << i << "\t" << it1->first;
            for(ListPolygon2D::const_iterator it2=it1->second.begin(); it2!=it1->second.end(); it2++) {
                for(int j=0; j<it2->size(); j++) {
                    out<<"\t"<<(*it2)[j].first<<"\t"<<(*it2)[j].second;
                }
            }
            out << endl;
        }
    }
    out.close();
}

void IGAMortarMapper::writeCartesianProjectedPolygon(const string _filename, std::map<int, ListPolygon2D>& _data) {
    ofstream out;
    string outName = name + "_" +_filename + ".vtk";
    out.open(outName.c_str(), ofstream::out);
    out << "# vtk DataFile Version 2.0\n";
    out << "Back projection of projected FE elements on NURBS mesh\n";
    out << "ASCII\nDATASET POLYDATA\n";
    string points, pointsHeader, polygons, polygonsHeader, patchColor, patchColorHeader;
    int pointsNumber=0, polygonsNumber=0, polygonsEntriesNumber=0;
    /// Loop over data
    for(map<int,ListPolygon2D>::iterator itPatch=_data.begin(); itPatch!=_data.end(); itPatch++) {
        int idPatch = itPatch->first;
        IGAPatchSurface* thePatch = meshIGA->getSurfacePatch(idPatch);
        for(ListPolygon2D::iterator itListPolygon=itPatch->second.begin(); itListPolygon!=itPatch->second.end(); itListPolygon++) {
            int nEdge=0;
            for(Polygon2D::iterator itPolygon=itListPolygon->begin(); itPolygon!=itListPolygon->end(); itPolygon++) {
                double local[2], global[3];
                local[0] = itPolygon->first;
                local[1] = itPolygon->second;
                thePatch->computeCartesianCoordinates(global,local);
                stringstream pointStream;
                pointStream << global[0];
                pointStream << " " << global[1];
                pointStream << " " << global[2];
                pointStream << "\n";
                points += pointStream.str();
                pointsNumber++;
                nEdge++;
            }
            stringstream colorStream;
            /// Concatenate new polygon color
            colorStream << idPatch << "\n";
            patchColor += colorStream.str();
            polygonsNumber++;
            polygonsEntriesNumber += nEdge + 1;
            /// Concatenate new polygon connectivity
            stringstream polygonStream;
            polygonStream << nEdge;
            for(int i=nEdge;i>0;i--) {
                polygonStream << " " << pointsNumber - i;
            }
            polygonStream << "\n";
            polygons += polygonStream.str();
        }
    }
    /// Write actually the file
    stringstream header;
    header << "POINTS " << pointsNumber << " float\n";
    pointsHeader = header.str();
    header.str("");
    header.clear();
    header << "POLYGONS " << polygonsNumber << " " << polygonsEntriesNumber <<"\n";
    polygonsHeader = header.str();
    header.str("");
    header.clear();
    header << "CELL_DATA " << polygonsNumber << "\nSCALARS patch_belonging int 1\nLOOKUP_TABLE default\n";
    patchColorHeader = header.str();
    out << pointsHeader << points;
    out << polygonsHeader << polygons;
    out << patchColorHeader << patchColor;
    out.close();
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

    ERROR_OUT() << "Cnn" << endl;
    couplingMatrices->getCnn()->printCSR();
    ERROR_OUT() << "Cnr" << endl;
    couplingMatrices->getCnr()->printCSR();
}

void IGAMortarMapper::writeCouplingMatricesToFile() {
    DEBUG_OUT()<<"### Printing matrices into file ###"<<endl;
    DEBUG_OUT()<<"Size of Cnr is "<<numNodesMaster<<" by "<<numNodesSlave<<endl;
    if(Message::isDebugMode()) {
        couplingMatrices->getCnr()->printCSRToFile(name + "_Cnr.dat",1);
        couplingMatrices->getCnn()->printCSRToFile(name + "_Cnn.dat",1);
    }
}

void IGAMortarMapper::enforceConsistency() {
    /*
     * The function checks the consistency up to the specified tolerance when mapping a unit field and if the consistency up
     * to the specified tolerance is violated then it is being enforced. In case after enforcing the consistency the mapping
     * of a unit field still violates the specified consistency tolerance an error is asserted.
     */

    INFO_OUT() << "Checking consistency" << std::endl;

    int size_N = couplingMatrices->getSizeN();
    int size_R = couplingMatrices->getSizeR();

    double ones[size_R];
    for(int i = 0; i < size_R; i++) {
        ones[i] = 1.0;
    }

    double output[size_N];
    this->consistentMapping(ones, output);

    double norm = 0;
    vector<int> inconsistentDoF;
    for(int i = 0; i < size_N; i++) {
        if(fabs(output[i] - 1) > propConsistency.tolConsistency && output[i] != 0)
            inconsistentDoF.push_back(i);
        norm += output[i]*output[i];
    }

    // Replace badly conditioned row of Cnn by sum value of Cnr
    if(!inconsistentDoF.empty()) {
        INFO_OUT() << "Mapping found inconsistent up to specified tolerance" << std::endl;
        INFO_OUT() << "inconsistent DOF size = " << inconsistentDoF.size() << std::endl;
        INFO_OUT() << "Enforcing consistency" << std::endl;
        for(vector<int>::iterator it = inconsistentDoF.begin(); it != inconsistentDoF.end(); it++) {
            couplingMatrices->deleterow(*it);       // deleterow might not be working for the new coupling matrix datastructure
            couplingMatrices->addCNNValue(*it ,*it , couplingMatrices->getCnr()->getRowSum(*it));
        }

        couplingMatrices->factorizeCnn();
        this->consistentMapping(ones, output);
        norm = 0;

        for(int i = 0; i < size_N; i++) {
            norm += output[i]*output[i];
        }
    } else
        INFO_OUT() << "Mapping found consistent up to specified tolerance" << std::endl;

    int denom = size_N - couplingMatrices->getIndexEmptyRowCnn().size();

    norm = sqrt(norm/denom);

    DEBUG_OUT() << "Checking consistency after enforcement"<<endl;
    DEBUG_OUT() << "Deviation from unit field : "<< fabs(norm - 1.0) << endl;
    if(fabs(norm - 1.0) > propConsistency.tolConsistency) {
        ERROR_OUT() << "Coupling matrices not consistent up to " << propConsistency.tolConsistency << "!"<<endl;
        exit(-1);
    }
}

void IGAMortarMapper::getPenaltyParameterForWeakDirichletCCPrimaryField(double* _alphaPrim){
    if( propWeakCurveDirichletConditions.isWeakCurveDirichletConditions ){
        for(int i = 0; i < noWeakIGADirichletCurveConditions; i++)
            _alphaPrim[i] = weakDirichletCCAlphaPrimary[i];
    }else{
        ERROR_OUT() << "Penalty parameters were not computed" << std::endl;
        exit(-1);
    }
}

void IGAMortarMapper::getPenaltyParameterForWeakDirichletCCSecondaryFieldBending(double* _alphaSecBending){
    if( propWeakCurveDirichletConditions.isWeakCurveDirichletConditions ){
        for(int i = 0; i < noWeakIGADirichletCurveConditions; i++)
            _alphaSecBending[i] = weakDirichletCCAlphaSecondaryBending[i];
    }else{
        ERROR_OUT() << "Penalty parameters were not computed" << std::endl;
        exit(-1);
    }
}

void IGAMortarMapper::getPenaltyParameterForWeakDirichletCCSecondaryFieldTwisting(double* _alphaSecTwisting){
    if( propWeakCurveDirichletConditions.isWeakCurveDirichletConditions ){
        for(int i = 0; i < noWeakIGADirichletCurveConditions; i++)
            _alphaSecTwisting[i] = weakDirichletCCAlphaSecondaryTwisting[i];
    }else{
        ERROR_OUT() << "Penalty parameters were not computed" << std::endl;
        exit(-1);
    }
}


void IGAMortarMapper::getPenaltyParameterForPatchContinuityPrimaryField(double* _alphaPrim){
    if(propWeakPatchContinuityConditions.isWeakPatchContinuityConditions){
        for(int i = 0; i < noWeakIGAPatchContinuityConditions; i++)
            _alphaPrim[i] = weakPatchContinuityAlphaPrimaryIJ[i];
    }else{
        ERROR_OUT() << "Penalty parameters were not computed" << std::endl;
        exit(-1);
    }
}

void IGAMortarMapper::getPenaltyParameterForPatchContinuitySecondaryFieldBending(double* _alphaSecBending){
    if(propWeakPatchContinuityConditions.isWeakPatchContinuityConditions){
        for(int i = 0; i < noWeakIGAPatchContinuityConditions; i++)
            _alphaSecBending[i] = weakPatchContinuityAlphaSecondaryBendingIJ[i];
    }else{
        ERROR_OUT() << "Penalty parameters were not computed" << std::endl;
        exit(-1);
    }
}

void IGAMortarMapper::getPenaltyParameterForPatchContinuitySecondaryFieldTwisting(double* _alphaSecTwisting){
    if(propWeakPatchContinuityConditions.isWeakPatchContinuityConditions){
        for(int i = 0; i < noWeakIGAPatchContinuityConditions; i++)
            _alphaSecTwisting[i] = weakPatchContinuityAlphaSecondaryTwistingIJ[i];
    }else{
        ERROR_OUT() << "Penalty parameters were not computed" << std::endl;
        exit(-1);
    }
}

void IGAMortarMapper::printErrorMessage(Message &message, double _errorL2Domain, double* _errorL2Curve, double* _errorL2Interface){
    message << std::endl;
    message() << "\t+" << "Mapping error: " << std::endl;
    if(propErrorComputation.isDomainError)
        message << "\t\t+" << '\t' << "L2 norm of the error in the domain: " << _errorL2Domain << std::endl;
    if(propErrorComputation.isCurveError){
        message << "\t\t+" << '\t' << "L2 norm of the field along the Dirichlet boundary: " << _errorL2Curve[0] << std::endl;
        message << "\t\t+" << '\t' << "L2 norm of the field rotation along the Dirichlet boundary: " << _errorL2Curve[1] << std::endl;
    }
    if(propErrorComputation.isInterfaceError){
        message << "\t\t+" << '\t' << "L2 norm of the field interface jump: " << _errorL2Interface[0] << std::endl;
        message << "\t\t+" << '\t' << "L2 norm of the field rotation interface jump: " << _errorL2Interface[1] << std::endl;
    }
    message() << "\t+" << "---------------------------------" << endl;
    message << std::endl;
}

}
