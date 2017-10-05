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


#include "flann/flann.hpp"

#include "IGAMortarMapper.h"
#include "IGAPatchSurface.h"
#include "WeakIGADirichletCurveCondition.h"
#include "WeakIGADirichletSurfaceCondition.h"
#include "WeakIGAPatchContinuityCondition.h"
#include "IGAMortarCouplingMatrices.h"
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
#include <iomanip>

using namespace std;

namespace EMPIRE {

double IGAMortarMapper::EPS_CLEANTRIANGLE = 1e-6;
double IGAMortarMapper::EPS_CLIPPING = 1e-9;

/// Declaration statement
static const string HEADER_DECLARATION = "Author: Andreas Apostolatos";

IGAMortarMapper::IGAMortarMapper(std::string _name, AbstractMesh *_meshA, AbstractMesh *_meshB): name(_name) {

    // Check input
    assert(_meshA != NULL);
    assert(_meshB != NULL);

    bool isMeshAIGA = (_meshA->type == EMPIRE_Mesh_IGAMesh);
    bool isMeshBIGA = (_meshB->type == EMPIRE_Mesh_IGAMesh);

    if (isMeshAIGA && !isMeshBIGA) {
        assert(_meshB->type == EMPIRE_Mesh_FEMesh || _meshB->type == EMPIRE_Mesh_SectionMesh);
        meshIGA = dynamic_cast<IGAMesh *>(_meshA);
        if (dynamic_cast<FEMesh *>(_meshB)->triangulate() == NULL)
            meshFE = dynamic_cast<FEMesh *>(_meshB);
        else
            meshFE = dynamic_cast<FEMesh *>(_meshB)->triangulate();
        isMappingIGA2FEM = true;
    } else if (!isMeshAIGA && isMeshBIGA) {
        assert(_meshA->type == EMPIRE_Mesh_FEMesh || _meshA->type == EMPIRE_Mesh_SectionMesh);
        if (dynamic_cast<FEMesh *>(_meshA)->triangulate() == NULL)
            meshFE = dynamic_cast<FEMesh *>(_meshA);
        else
            meshFE = dynamic_cast<FEMesh *>(_meshA)->triangulate();
        meshIGA = dynamic_cast<IGAMesh *>(_meshB);
        isMappingIGA2FEM = false;
    } else {
        ERROR_BLOCK_OUT("IGAMortarMapper","IGAMortarMapper","Wrong type of mesh! Put a NURBS mesh and a FE mesh!");
    }

    // Number of Patches
    int numPatches = getIGAMesh()->getNumPatches();

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

    // Initialize flag on whether the meshFEDirectElemTable was created
    isMeshFEDirectElemTable = false;

    // Initialize flag on the expansion of the coupling matrices
    isExpanded = false;

    // Initialize flag on the initialization of the coupling matrices
    isCouplingMatrices = false;

    // Initialize Gauss quadratures
    gaussRuleOnTriangle = new EMPIRE::MathLibrary::IGAGaussQuadrature*[numPatches];
    gaussRuleOnQuadrilateral = new EMPIRE::MathLibrary::IGAGaussQuadrature*[numPatches];
    isGaussQuadature = false;

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

void IGAMortarMapper::setParametersIntegration(bool _isAutomaticNoGPTriangle, int _noGPTriangle, bool _isAutomaticNoGPQuadrilateral, int _noGPQuadrilateral) {
    propIntegration.isAutomaticNoGPTriangle = _isAutomaticNoGPTriangle;
    propIntegration.noGPTriangle = _noGPTriangle;
    propIntegration.isAutomaticNoGPQuadrilateral = _isAutomaticNoGPQuadrilateral;
    propIntegration.noGPQuadrilateral = _noGPQuadrilateral;
}

void IGAMortarMapper::setParametersWeakCurveDirichletConditions(bool _isWeakCurveDirichletConditions, bool _isAutomaticPenaltyParameters,
                                                                bool _isPrimPrescribed, bool _isSecBendingPrescribed, bool _isSecTwistingPrescribed,
                                                                double _alphaPrim, double _alphaSecBending, double _alphaSecTwisting) {
    propWeakCurveDirichletConditions.isWeakCurveDirichletConditions = _isWeakCurveDirichletConditions;
    propWeakCurveDirichletConditions.isAutomaticPenaltyParameters = _isAutomaticPenaltyParameters;
    propWeakCurveDirichletConditions.isPrimPrescribed = _isPrimPrescribed;
    propWeakCurveDirichletConditions.isSecBendingPrescribed = _isSecBendingPrescribed;
    propWeakCurveDirichletConditions.isSecTwistingPrescribed = _isSecTwistingPrescribed;
    propWeakCurveDirichletConditions.alphaPrim = _alphaPrim;
    propWeakCurveDirichletConditions.alphaSecBending = _alphaSecBending;
    propWeakCurveDirichletConditions.alphaSecTwisting = _alphaSecTwisting;
}

void IGAMortarMapper::setParametersWeakSurfaceDirichletConditions(bool _isWeakSurfaceDirichletConditions, bool _isAutomaticPenaltyParameters,
                                                                  bool _isPrimPrescribed, double _alphaPrim) {
    propWeakSurfaceDirichletConditions.isWeakSurfaceDirichletConditions = _isWeakSurfaceDirichletConditions;
    propWeakSurfaceDirichletConditions.isAutomaticPenaltyParameters = _isAutomaticPenaltyParameters;
    propWeakSurfaceDirichletConditions.isPrimPrescribed = _isPrimPrescribed;
    propWeakSurfaceDirichletConditions.alphaPrim = _alphaPrim;
}

void IGAMortarMapper::setParametersWeakPatchContinuityConditions(bool _isWeakPatchContinuityConditions, bool _isAutomaticPenaltyParameters,
                                                                 bool _isPrimCoupled, bool _isSecBendingCoupled, bool _isSecTwistingCoupled,
                                                                 double _alphaPrim, double _alphaSecBending, double _alphaSecTwisting) {
    propWeakPatchContinuityConditions.isWeakPatchContinuityConditions = _isWeakPatchContinuityConditions;
    propWeakPatchContinuityConditions.isAutomaticPenaltyParameters = _isAutomaticPenaltyParameters;
    propWeakPatchContinuityConditions.isPrimCoupled = _isPrimCoupled;
    propWeakPatchContinuityConditions.isSecBendingCoupled = _isSecBendingCoupled;
    propWeakPatchContinuityConditions.isSecTwistingCoupled = _isSecTwistingCoupled;
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
    /*
     * Builds the coupling matrices Cnn and Cnr for the isogeometric mortar-based mapper.
     *
     * Function layout :
     *
     * 0. Print message
     *
     * 1. Check input
     *
     * 2. Initialize coupling matrices
     *
     * 3. Initialize auxiliary variables
     *
     * 4. Create the Gauss quadrature rules for each patch
     *
     * 5. Find the maximum number of Gauss points within an element
     *
     * 6. Project the FE nodes onto the multipatch trimmed geometry
     *
     * 7. Write the projected points on to a file only in DEBUG mode to be used in MATLAB
     *
     * 8. Reserve some space for gauss point values in the domain (the size is not known exactly in advance)
     *
     * 9. Reserve space for the gauss point values along each trimming curve where conditions are applied
     *
     * 10. Reserve space for the interface gauss point values
     *
     * 11. Compute mortar coupling matrices
     *
     * 12. Write the gauss point and the coupling matrices data in files
     *
     * 13. Write polygon net of projected elements to a vtk file
     *
     * 14. Compute the Penalty parameters for the application of weak Dirichlet curve conditions
     *
     * 15. Compute the Penalty parameters for the application of weak Dirichlet surface conditions
     *
     * 16. Compute the Penalty parameters for the application of weak patch continuity conditions
     *
     * 17. Compute the Penalty matrices for the application of weak continuity conditions between the multipatches
     *
     * 18. Remove empty rows and columns from system (flying nodes)
     *
     * 19. Check and enforce consistency in Cnn. This has to be done before the weak application of the Dirichlet conditions.
     *
     * 20. Compute the Penalty matrices for the application of weak Dirichlet conditions along trimming curves
     *
     * 21. Compute the Penalty matrices for the application of weak Dirichlet conditions across surfaces
     *
     * 22. Factorize Cnn matrix
     */

    // 0. Print message
    HEADING_OUT(3, "IGAMortarMapper", "Building coupling matrices for ("+ name +")...", infoOut);
    {
        int sizeN = couplingMatrices->getSizeN();
        int sizeR = couplingMatrices->getSizeR();
        if (isExpanded) {
            INFO_OUT() << "The number of DOFs in the finite element mesh is " << (isMappingIGA2FEM?sizeN:sizeR) << endl;
            INFO_OUT() << "The number of DOFs in NURBS multipatch geometry is " << (isMappingIGA2FEM?sizeR:sizeN) << endl;
            INFO_OUT() << "Isogeometric mortar-based mapping for all Cartesian components simultaneously" << endl;
        } else {
            INFO_OUT() << "The number of nodes in the finite element mesh is " << (isMappingIGA2FEM?sizeN:sizeR) << endl;
            INFO_OUT() << "The number of Control Points in the NURBS multipatch geometry is " << (isMappingIGA2FEM?sizeR:sizeN) << endl;
            INFO_OUT() << "Isogeometric mortar-based mapping per each Cartesian component seperately" << endl;
        }
        INFO_OUT() << "The size of matrices is " << "Cnn: " << (isMappingIGA2FEM?sizeN:sizeR) << " x " << (isMappingIGA2FEM?sizeN:sizeR) << " and "
                   << "Cnr: " << (isMappingIGA2FEM?sizeN:sizeR) << " x " << (isMappingIGA2FEM?sizeR:sizeN) << endl;
    }

    // 1. Check input
    if (propErrorComputation.isErrorComputation)
        if (propErrorComputation.isCurveError && !propWeakCurveDirichletConditions.isWeakCurveDirichletConditions){
            ERROR_OUT() << "Error in MapperAdapter::initIGAMortarMapper" << endl;
            ERROR_OUT() << "Error computation along trimming curves is requested but no conditions along trimming curves are prescribed" << endl;
            exit(-1);
        }
        if (propErrorComputation.isInterfaceError && !propWeakPatchContinuityConditions.isWeakPatchContinuityConditions) {
            ERROR_OUT() << "Error in MapperAdapter::initIGAMortarMapper" << endl;
            ERROR_OUT() << "Error computation along patch interfaces is requested but no conditions along the patch interfaces are prescribed" << endl;
            exit(-1);
        }
    if(propWeakCurveDirichletConditions.isWeakCurveDirichletConditions && !isMappingIGA2FEM) {
        assert(isExpanded);
    }
    if(propWeakPatchContinuityConditions.isWeakPatchContinuityConditions && !isMappingIGA2FEM) {
        assert(isExpanded);
    }

    // 2. Initialize coupling matrices
    initialize();

    // 3. Initialize auxiliary variables
    int numPatches = getIGAMesh()->getNumPatches();
    string filename; // String holding the mapper names
    IGAPatchSurface::MAX_NUM_ITERATIONS = propNewtonRaphson.noIterations; // Set default scheme values
    IGAPatchSurface::TOL_ORTHOGONALITY = propNewtonRaphson.tolProjection; // Set default scheme values

    // 4. Create the Gauss quadrature rules for each patch
    createGaussQuadratureRules();

    // 5. Find the maximum number of Gauss points within an element
    int maxNumGP = 0;
    int numGP;
    for (int iPatches = 0; iPatches < numPatches; iPatches++) {
        numGP = gaussRuleOnTriangle[iPatches]->getNumGaussPoints();
        if (numGP > maxNumGP)
            maxNumGP = numGP;
        numGP = gaussRuleOnQuadrilateral[iPatches]->getNumGaussPoints();
        if (numGP > maxNumGP)
            maxNumGP = numGP;
    }

    // 6. Project the FE nodes onto the multipatch trimmed geometry
    projectPointsToSurface();

    // 7. Write the projected points on to a file only in DEBUG mode to be used in MATLAB
    if (Message::isDebugMode())
        writeProjectedNodesOntoIGAMesh();

    // 8. Reserve some space for gauss point values in the domain (the size is not known exactly in advance)
    if (propErrorComputation.isDomainError)
        streamGPs.reserve(8*meshFE->numElems*maxNumGP);

    // 9. Reserve space for the gauss point values along each trimming curve where conditions are applied
    if(propErrorComputation.isCurveError){
        int noCurveGPs = 0;
        std::vector<WeakIGADirichletCurveCondition*> weakIGADirichletCurveConditions = meshIGA->getWeakIGADirichletCurveConditions();
        for (int iWDC = 0; iWDC < weakIGADirichletCurveConditions.size(); iWDC++){
            noCurveGPs += weakIGADirichletCurveConditions[iWDC]->getCurveNumGP();
        }
        streamInterfaceGPs.reserve(noCurveGPs);
    }
    // 10. Reserve space for the interface gauss point values
    if(propErrorComputation.isInterfaceError){
        int noInterfaceGPs = 0;
        std::vector<WeakIGAPatchContinuityCondition*> weakIGAPatchContinuityConditions = meshIGA->getWeakIGAPatchContinuityConditions();
        for (int iWCC = 0; iWCC < weakIGAPatchContinuityConditions.size(); iWCC++){
            noInterfaceGPs += weakIGAPatchContinuityConditions[iWCC]->getTrCurveNumGP();
        }
        streamInterfaceGPs.reserve(noInterfaceGPs);
    }

    // 11. Compute mortar coupling matrices
    computeCouplingMatrices();

    // 12. Write the gauss point and the coupling matrices data in files
    if(Message::isDebugMode()) {
        writeGaussPointData();
        writeCouplingMatricesToFile();
    }

    // 13. Write polygon net of projected elements to a vtk file
    writeCartesianProjectedPolygon("trimmedPolygonsOntoNURBSSurface", trimmedProjectedPolygons);
    writeCartesianProjectedPolygon("integratedPolygonsOntoNURBSSurface", triangulatedProjectedPolygons2);
    trimmedProjectedPolygons.clear();
    triangulatedProjectedPolygons2.clear();

    // 14. Compute the Penalty parameters for the application of weak Dirichlet curve conditions
    if(propWeakCurveDirichletConditions.isWeakCurveDirichletConditions && !isMappingIGA2FEM) {
        filename = name + "_penaltyParametersWeakDirichletConditions.txt";
        computePenaltyParametersForWeakDirichletCurveConditions(filename);
    }

    // 15. Compute the Penalty parameters for the application of weak Dirichlet surface conditions
    if(propWeakSurfaceDirichletConditions.isWeakSurfaceDirichletConditions && !isMappingIGA2FEM)
        computePenaltyParametersForWeakDirichletSurfaceConditions();

    // 16. Compute the Penalty parameters for the application of weak patch continuity conditions
    if(propWeakPatchContinuityConditions.isWeakPatchContinuityConditions && !isMappingIGA2FEM) {
        filename = name + "_penaltyParametersWeakContinuityConditions.txt";
        computePenaltyParametersForPatchContinuityConditions(filename);
    }

    // 17. Compute the Penalty matrices for the application of weak continuity conditions between the multipatches
    if (propWeakPatchContinuityConditions.isWeakPatchContinuityConditions && !isMappingIGA2FEM) {
        INFO_OUT() << "Application of weak patch continuity conditions started" << endl;
        if(!propWeakPatchContinuityConditions.isAutomaticPenaltyParameters) {
            INFO_OUT() << "Manual assignment of the penalty parameters: alphaPrim = "<< propWeakPatchContinuityConditions.alphaPrim
                       << " alphaSecBending = " <<  propWeakPatchContinuityConditions.alphaSecBending <<" alphaSecTwisting = "
                       <<  propWeakPatchContinuityConditions.alphaSecTwisting << endl;
        } else {
            INFO_OUT() << "Automatic determination of the penalty parameters are assumed" << endl;
        }
        computeIGAPatchWeakContinuityConditionMatrices();
        INFO_OUT() << "Application of weak patch continuity conditions finished" << std::endl;
    } else
        INFO_OUT() << "No application of weak patch continuity conditions is assumed" << std::endl;

    // 18. Remove empty rows and columns from system (flying nodes)
    if(!isMappingIGA2FEM){
        INFO_OUT() << "Enforcing flying nodes in Cnn" << std::endl;
        couplingMatrices->enforceCnn();
    }

    // 19. Check and enforce consistency in Cnn. This has to be done before the weak application of the Dirichlet conditions.
    if (propConsistency.enforceConsistency)
        enforceConsistency();

    // 20. Compute the Penalty matrices for the application of weak Dirichlet conditions along trimming curves
    if (propWeakCurveDirichletConditions.isWeakCurveDirichletConditions) {
        INFO_OUT() << "Application of weak Dirichlet curve conditions started" << endl;
        if(!propWeakCurveDirichletConditions.isAutomaticPenaltyParameters) {
            INFO_OUT() << "Manual assignment of the penalty parameters: alphaPrim = " << propWeakCurveDirichletConditions.alphaPrim << " alphaSecBending = " <<  propWeakCurveDirichletConditions.alphaSecBending << " alphaSecTwisting = " <<  propWeakCurveDirichletConditions.alphaSecTwisting << endl;
        } else {
            INFO_OUT() << "Automatic determination of the penalty parameters are assumed" << endl;
        }
        computeIGAWeakDirichletCurveConditionMatrices();
        INFO_OUT() << "Application of weak Dirichlet curve conditions finished" << std::endl;
    } else
        INFO_OUT() << "No application of weak Dirichlet curve conditions are assumed" << std::endl;

    // 21. Compute the Penalty matrices for the application of weak Dirichlet conditions across surfaces
    if (propWeakSurfaceDirichletConditions.isWeakSurfaceDirichletConditions) {
        INFO_OUT() << "Application of weak Dirichlet surface conditions started" << endl;
        if(!propWeakSurfaceDirichletConditions.isAutomaticPenaltyParameters) {
            INFO_OUT() << "Manual assignment of the penalty parameters: alphaPrim = "<< propWeakSurfaceDirichletConditions.alphaPrim << " alphaSecBending = " << endl;
        } else {
            INFO_OUT() << "Automatic determination of the penalty parameters are assumed" << endl;
        }
        computeIGAWeakDirichletSurfaceConditionMatrices();
        INFO_OUT() << "Application of weak Dirichlet surface conditions finished" << std::endl;
    } else
        INFO_OUT() << "No application of weak Dirichlet surface conditions are assumed" << std::endl;

    // 22. Factorize Cnn matrix
    couplingMatrices->factorizeCnn();
    INFO_OUT() << "Factorize was successful" << std::endl;
}

IGAMortarMapper::~IGAMortarMapper() {
    // Initialize auxiliary variables
    int numPatches = getIGAMesh()->getNumPatches();

    // Delete the Finite Element direct element table
    if(isMeshFEDirectElemTable){
        for (int i = 0; i < meshFE->numElems; i++)
            delete[] meshFEDirectElemTable[i];
        delete[] meshFEDirectElemTable;
    }

    // Delete the quadrature rules
    if (isGaussQuadature)
        for (int iPatches = 0; iPatches < numPatches; iPatches++) {
                delete gaussRuleOnTriangle[iPatches];
                delete gaussRuleOnQuadrilateral[iPatches];
        }
    delete[] gaussRuleOnTriangle;
    delete[] gaussRuleOnQuadrilateral;

    // Delete the coupling matrices
    if (isCouplingMatrices)
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
    /*
     * Create the map to store the nodeIDs but here the "key" is the node ID, and the value is the position in nodeIDs
     * the map is sorted automatically, so it is efficient for searching
     *
     * Function layout:
     *
     * 1. Initialize the direct element freedom tables for the Finite Element mesh
     *
     * 2. Make a map between the node ids and the local indexing
     *
     * 3. Loop over all elements in the Finite Element mesh
     * ->
     *    3i. Get the number of nodes per element for the given element
     *   3ii. Loop over all nodes in the Finite Element
     *   ->
     *        3ii.1. Assert an error if an ID is not found in the map
     *        3ii.2. Insert the id of the node in the map
     *   <-
     * <-
     *
     * 4. Loop over all nodes in the Finite Element mesh
     * ->
     *    4i. Loop over all elements in the Finite Element mesh
     *    ->
     *        4i.1. Get the number of nodes for the given element
     *        4i.2. Find the given node in the current element
     *        4i.3. Insert the element in the node map
     *    <-
     * <-
     *
     * 5. Set flag of the initialization of the flag
     */

    // 1. Initialize the direct element freedom tables for the Finite Element mesh
    meshFEDirectElemTable = new int*[meshFE->numElems]; // deleted
    for (int i = 0; i < meshFE->numElems; i++)
        meshFEDirectElemTable[i] = new int[meshFE->numNodesPerElem[i]];

    // 2. Make a map between the node ids and the local indexing
    map<int, int> meshFENodesMap;
    for (int i = 0; i < meshFE->numNodes; i++)
        meshFENodesMap.insert(meshFENodesMap.end(), pair<int, int>(meshFE->nodeIDs[i], i));
    int count = 0;

    // 3. Loop over all elements in the Finite Element mesh
    for (int i = 0; i < meshFE->numElems; i++) {
        // 3i. Get the number of nodes per element for the given element
        const int numNodesPerElem = meshFE->numNodesPerElem[i];

        // 3ii. Loop over all nodes in the Finite Element
        for (int j = 0; j < numNodesPerElem; j++) {
            // 3ii.1. Assert an error if an ID is not found in the map
            if (meshFENodesMap.find(meshFE->elems[count + j]) == meshFENodesMap.end()) {
                ERROR_OUT() << "Cannot find node ID " << meshFE->elems[count + j] << endl;
                exit(-1);
            }

            // 3ii.2. Insert the id of the node in the map
            meshFEDirectElemTable[i][j] = meshFENodesMap.at(meshFE->elems[count + j]);
        }

        // 3iii. Updated counter
        count += numNodesPerElem;
    }

    // 4. Loop over all nodes in the Finite Element mesh
    for (int iNodes = 0; iNodes < meshFE->numNodes; iNodes++) {
        // 4i. Loop over all elements in the Finite Element mesh
        for (int elem = 0; elem < meshFE->numElems; elem++) {
            // 4i.1. Get the number of nodes for the given element
            const int numNodesPerElem = meshFE->numNodesPerElem[elem];

            // 4i.2. Find the given node in the current element
            int* out = find(meshFEDirectElemTable[elem],meshFEDirectElemTable[elem]+numNodesPerElem,iNodes);

            // 4i.3. Insert the element in the node map
            if(out != meshFEDirectElemTable[elem] + numNodesPerElem) {
                meshFENodeToElementTable[iNodes].push_back(elem);
            }
        }
    }

    // 5. Set flag of the initialization of the flag
    isMeshFEDirectElemTable = true;
}

void IGAMortarMapper::initCouplingMatrices() {
    /*
     * Initialize the coupling matrices.
     *
     * Function layout:
     *
     * 1. Initialize auxiliary arrays
     *
     * 2. Get the number of rows and columns of the coupling matrices
     *
     * 3. Initialize coupling matrices
     *
     * 4. Set flag on the initialization of the coupling matrices accordingly
     */

    // 1. Initialize auxiliary arrays
    int noCoord = 3;

    // 2. Get the number of rows and columns of the coupling matrices
    int size_N;
    int size_R;
    if (isExpanded){
        size_N = noCoord*numNodesMaster;
        size_R = noCoord*numNodesSlave;
    } else {
        size_N = numNodesMaster;
        size_R = numNodesSlave;
    }

    // 3. Initialize coupling matrices
    couplingMatrices = new IGAMortarCouplingMatrices(size_N, size_R, isExpanded);

    // 4. Set flag on the initialization of the coupling matrices accordingly
    isCouplingMatrices = true;
}

void IGAMortarMapper::initialize() {
    /*
     * Initialize the freedom tables and the coupling matrices.
     *
     * Function layout :
     *
     * 1. Initialize freedom tables
     *
     * 2. Assign flag on whether the expanded version of the coupling matrices is assumed
     *
     * 3. Initialize coupling matrices
     */

    // 1. Initialize freedom tables
    initTables();

    // 2. Assign flag on whether the expanded version of the coupling matrices is assumed
    isExpanded = (propWeakCurveDirichletConditions.isWeakCurveDirichletConditions ||
                  propWeakSurfaceDirichletConditions.isWeakSurfaceDirichletConditions ||
                  propWeakPatchContinuityConditions.isWeakPatchContinuityConditions) && ~isMappingIGA2FEM;

    // 3. Initialize coupling matrices
    initCouplingMatrices();
}


void IGAMortarMapper::projectPointsToSurface() {

    int numCoord = 3;

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

    // Get the number of patches in the IGA mesh
    int numPatches = getIGAMesh()->getNumPatches();

    // Lists to keep track of projection nodes and patches
    vector<set<int> > patchIndicesToProcessPerNode(meshFE->numNodes);
    vector<set<int> > nodeIndicesToProcessPerPatch(numPatches);
    vector<vector< double > > nodeCoordsToProcessPerPatch(numPatches);
    set<int> notProjectedNodeIndicesFirstPass;
    set<int> notProjectedNodeIndicesSecondPass;

    // Initial guess for projection onto the NURBS patch
    double initialU, initialV;

    // Variables for generating initial guesses both in the patch parametric spaces and in the physical space
    int numUCP, numVCP;
    double* candidatesU;
    double* candidatesV;
    double* candidatesXYZ;
    double tmpUV[2];
    double* P;

    // Pointer to the patch to process
    IGAPatchSurface* thePatch;

    INFO_OUT() << "Bounding box preprocessing started" << endl;
    time(&timeStart);
    for (int iPatch = 0; iPatch < numPatches; iPatch++) {
        thePatch = meshIGA->getSurfacePatch(iPatch);
        for (int iNode = 0; iNode < meshFE->numNodes; iNode++) {
            P = &meshFE->nodes[numCoord * iNode];
            bool isInside = thePatch->getBoundingBox().isPointInside(P, propProjection.maxProjectionDistance);
            if(isInside) {
                nodeIndicesToProcessPerPatch[iPatch].insert(iNode);
                std::copy(P, P + numCoord, std::back_inserter(nodeCoordsToProcessPerPatch[iPatch]));
                patchIndicesToProcessPerNode[iNode].insert(iPatch);
            }
        }
        if(nodeIndicesToProcessPerPatch[iPatch].empty()) {
            stringstream msg;
            msg << "Patch [" << iPatch << "] does not have any nodes in its bounding box! Increase maxProjectionDistance !";
            ERROR_BLOCK_OUT("IGAMortarMapper", "projectPointsToSurface", msg.str());
        }
    }
    time(&timeEnd);
    INFO_OUT()<<"Bounding box preprocessing done in "<< difftime(timeEnd, timeStart) << " seconds"<<endl;

    INFO_OUT()<<"First pass projection started"<<endl;
    time(&timeStart);
    for (int iPatch = 0; iPatch < numPatches; iPatch++) {
        // Get the patch to project points onto
        thePatch = meshIGA->getSurfacePatch(iPatch);

        // Generate candidate points for initial guesses in the patch parameter space using the Greville abcissae
        numUCP = thePatch->getUNoControlPoints();
        numVCP = thePatch->getVNoControlPoints();

        candidatesU = new double[numUCP];
        for (int iUCP = 0; iUCP < numUCP; iUCP++)
            candidatesU[iUCP] = thePatch->getIGABasis(0)->computeGrevilleAbscissae(iUCP);

        candidatesV = new double[numVCP];
        for (int iVCP = 0; iVCP < numVCP; iVCP++)
            candidatesV[iVCP] = thePatch->getIGABasis(1)->computeGrevilleAbscissae(iVCP);

        // Compute the physical coordinates of the initial guess points
        candidatesXYZ = new double[numUCP*numVCP*numCoord];
        for (int iVCP = 0; iVCP < numVCP; iVCP++) {
            tmpUV[1] = candidatesV[iVCP];
            for (int iUCP = 0; iUCP < numUCP; iUCP++) {
                tmpUV[0] = candidatesU[iUCP];
                thePatch->computeCartesianCoordinates(&candidatesXYZ[(iVCP*numUCP+iUCP)*numCoord],tmpUV);
            }
        }

        #ifdef ANN
            double dummy;
            int* patchCandidateIndices[nodeIndicesToProcessPerPatch[iPatch].size()];
            ANNkd_tree ANNKdTree(candidatesXYZ, numUCP*numVCP, numCoord);
            ANNKdTree->annkSearch(&nodeCoordsToProcessPerPatch[iPatch][0], 1, patchCandidateIndices, &dummy);
        #endif

        #ifdef FLANN
            // Store the nodes inside the bounding box of the current patch and the initial guess candidates in FLANN types
            flann::Matrix<double> FLANNpatchNodesXYZ(&nodeCoordsToProcessPerPatch[iPatch][0], nodeIndicesToProcessPerPatch[iPatch].size(), numCoord);
            flann::Matrix<double> FLANNcandidatesXYZ(candidatesXYZ, numUCP*numVCP, numCoord);

            // Construct the FLANN KDTree search tree
            flann::Index<flann::L2<double> > FLANNKdTree(FLANNcandidatesXYZ, flann::KDTreeSingleIndexParams(1));
            FLANNKdTree.buildIndex();

            // Find the initial guess indices for each node inside the current patch bounding box
            vector<vector<int> > patchCandidateIndices;
            vector<vector<double> > dummyDistances;
            FLANNKdTree.knnSearch(FLANNpatchNodesXYZ, patchCandidateIndices, dummyDistances, 1, flann::SearchParams(1));
            dummyDistances.clear();
        #endif

        // Retrieve the initial guess coordinates from the patch parametric space(UV)
        std::div_t dv;
        int patchNodeIndex = 0;
        // Loop over the nodes in the current patch bounding box
        for(set<int>::iterator iNode = nodeIndicesToProcessPerPatch[iPatch].begin(); iNode != nodeIndicesToProcessPerPatch[iPatch].end(); iNode++) {
            // The quotient = iVCP, remainder = iUCP, since the ordering of the initial guesses were done in such order. See the ordering of "candidatesXYZ" above.
            #ifdef ANN
                dv = std::div(patchCandidateIndices[patchNodeIndex], numUCP);
            #endif
            #ifdef FLANN
                dv = std::div(patchCandidateIndices[patchNodeIndex][0], numUCP);
            #endif

            initialU = candidatesU[dv.rem];
            initialV = candidatesV[dv.quot];
            bool flagProjected = projectPointOnPatch(iPatch, *iNode, initialU, initialV, minProjectionDistance[*iNode], minProjectionPoint[*iNode]);
            isProjected[*iNode] = isProjected[*iNode] || flagProjected;

            patchNodeIndex++;
        }

        delete[] candidatesU;
        delete[] candidatesV;
        delete[] candidatesXYZ;
    }
    time(&timeEnd);
    INFO_OUT()<<"First pass projection done in "<< difftime(timeEnd, timeStart) << " seconds"<<endl;

    int missing = 0;
    for (int iNode = 0; iNode < meshFE->numNodes; iNode++) {
        if(!isProjected[iNode]) {
            missing++;
            notProjectedNodeIndicesFirstPass.insert(iNode);
            WARNING_OUT()<<"Node ["<<iNode<<"] not projected at first pass with coordinates "<<meshFE->nodes[iNode*numCoord]<<","<<meshFE->nodes[iNode*numCoord+1]<<","<<meshFE->nodes[iNode*numCoord+2]<<endl;

        }
    }

    // Second pass projection --> relax Newton-Rapshon tolerance and if still fails refine the sampling points for the Newton-Raphson initial guesses
    if(missing) {
        INFO_OUT()<< missing << " out of " << meshFE->numNodes <<" nodes could NOT be projected during first pass" << endl;
        INFO_OUT()<<"Second pass projection started"<<endl;
        time(&timeStart);
        missing = 0;
        for(set<int>::iterator iNode = notProjectedNodeIndicesFirstPass.begin(); iNode != notProjectedNodeIndicesFirstPass.end(); iNode++) {
            for(set<int>::iterator iPatch = patchIndicesToProcessPerNode[*iNode].begin();iPatch != patchIndicesToProcessPerNode[*iNode].end(); iPatch++) {
                computeInitialGuessForProjection(*iPatch, meshFENodeToElementTable[*iNode][0], *iNode, initialU, initialV);
                bool flagProjected = forceProjectPointOnPatchByRelaxation(*iPatch, *iNode, initialU, initialV, minProjectionDistance[*iNode], minProjectionPoint[*iNode]);
                isProjected[*iNode] = isProjected[*iNode] || flagProjected;
            }
            if(!isProjected[*iNode]) {
                notProjectedNodeIndicesSecondPass.insert(*iNode);
                WARNING_OUT()<<"Node ["<<*iNode<<"] not projected at second pass with coordinates "<<meshFE->nodes[(*iNode)*numCoord]<<","<<meshFE->nodes[(*iNode)*numCoord+1]<<","<<meshFE->nodes[(*iNode)*numCoord+2]<<endl;
                missing++;
            }
        }
        notProjectedNodeIndicesFirstPass.clear();
        time(&timeEnd);
        INFO_OUT()<<"Second pass projection done in "<< difftime(timeEnd, timeStart) << " seconds"<<endl;
    }

    // Third pass projection --> Closest point projection based on brute sampling
    if(missing) {
        INFO_OUT()<< missing << " out of " << meshFE->numNodes <<" nodes could NOT be projected during second pass" << endl;
        INFO_OUT()<<"Third pass projection started"<<endl;
        time(&timeStart);
        missing = 0;
        for(set<int>::iterator iNode = notProjectedNodeIndicesSecondPass.begin(); iNode != notProjectedNodeIndicesSecondPass.end(); iNode++) {
            for(set<int>::iterator iPatch = patchIndicesToProcessPerNode[*iNode].begin(); iPatch != patchIndicesToProcessPerNode[*iNode].end(); iPatch++) {
                bool flagProjected = forceProjectPointOnPatchBySampling(*iPatch, *iNode, minProjectionDistance[*iNode], minProjectionPoint[*iNode]);
                isProjected[*iNode] = isProjected[*iNode] || flagProjected;
            }
            if(!isProjected[*iNode]) {
                WARNING_OUT()<<"Node ["<<*iNode<<"] not projected at third pass with coordinates "<<meshFE->nodes[(*iNode)*numCoord]<<","<<meshFE->nodes[(*iNode)*numCoord+1]<<","<<meshFE->nodes[(*iNode)*numCoord+2]<<endl;
                missing++;
            }
        }
        notProjectedNodeIndicesSecondPass.clear();
        time(&timeEnd);
        INFO_OUT()<<"Third pass projection done in "<< difftime(timeEnd, timeStart) << " seconds"<<endl;
    }

    if(missing) {
        stringstream msg;
        msg << missing << " nodes over " << meshFE->numNodes << " could NOT be projected!" << endl;
        msg << "Treatment possibility 1." << endl;
        msg << "Possibly relax parameters in projectionProperties or newtonRaphson" << endl;
        msg << "Treatment possibility 2." << endl;
        msg << "Remesh with higher accuracy on coordinates of the FE nodes, i.e. more digits" << endl;
        ERROR_BLOCK_OUT("IGAMortarMapper", "projectPointsToSurface", msg.str());
    }

}

void IGAMortarMapper::computeInitialGuessForProjection(const int _patchIndex, const int _elemIndex, const int _nodeIndex, double& _u, double& _v) {
    /*
     * Finds an initial guess for the projection of a node on the given patch
     *
     * Function layout:
     *
     * 1. Initialize auxiliary arrays
     *
     * 2. Get the patch from the multipatch geometry
     *
     * 3. Loop over all nodes of the current element
     * ->
     *    3i. Get the index of the node from the direct element freedom table
     *   3ii. Check if the node has already been projected and get its index
     * <-
     *
     * 4. Check if there exist one node in the current element which was successfully projected
     * ### If so, use result of the projected node as the initial guess for the projection step ###
     * ### Otherwise, find the nearest knot intersection as initial guess for the projection step ###
     */

    // 1. Initialize auxiliary arrays
    bool isNodeInsideElementProjected = false;
    int projectedNode = -1;

    // 2. Get the patch from the multipatch geometry
    IGAPatchSurface* thePatch = meshIGA->getSurfacePatch(_patchIndex);

    // 3. Loop over all nodes of the current element
    for (int iNodesElmnt = 0; iNodesElmnt < meshFE->numNodesPerElem[_elemIndex]; iNodesElmnt++) {
        // 3i. Get the index of the node from the direct element freedom table
        int nodeIndex = meshFEDirectElemTable[_elemIndex][iNodesElmnt];

        // 3ii. Check if the node has already been projected and get its index
        if (projectedCoords[nodeIndex].find(_patchIndex) != projectedCoords[nodeIndex].end()) {
            isNodeInsideElementProjected = true;
            projectedNode = nodeIndex;
            break;
        }
    }

    // 4. Check if there exist one node in the current element which was successfully projected
    if (isNodeInsideElementProjected) { // ### If so, use result of the projected node as the initial guess for the projection step ###
        _u = projectedCoords[projectedNode][_patchIndex][0];
        _v = projectedCoords[projectedNode][_patchIndex][1];
    } else { // ### Otherwise, find the nearest knot intersection as initial guess for the projection step ###
        double P[3];
        P[0] = meshFE->nodes[_nodeIndex * 3 + 0];
        P[1] = meshFE->nodes[_nodeIndex * 3 + 1];
        P[2] = meshFE->nodes[_nodeIndex * 3 + 2];
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
                                                                hasResidualConverged, propNewtonRaphson.noIterations, propNewtonRaphson.tolProjection, propProjection.maxProjectionDistance);
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

bool IGAMortarMapper::forceProjectPointOnPatchByRelaxation(const int patchIndex, const int nodeIndex, const double u0, const double v0, double& minProjectionDistance, vector<double>& minProjectionPoint) {

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
    bool hasConverged = thePatch->computeForcedPointProjectionOnPatch(u, v, projectedP);
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

bool IGAMortarMapper::forceProjectPointOnPatchBySampling(const int patchIndex, const int nodeIndex, double& minProjectionDistance, vector<double>& minProjectionPoint) {
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
            int patchIndex = *it;

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
            if(listPolygonUV[index].size() < 3)
                continue;
            isIntegrated = true;
            ListPolygon2D triangulatedPolygons = triangulatePolygon(listPolygonUV[index]);
            /// 1.3.3 For each triangle, compute canonical element and integrate
            for(ListPolygon2D::iterator triangulatedPolygon = triangulatedPolygons.begin(); triangulatedPolygon != triangulatedPolygons.end(); triangulatedPolygon++) {
                /// WARNING hard coded tolerance. Cleaning of triangle. Avoid heavily distorted triangle to go further.
                ClipperAdapter::cleanPolygon(*triangulatedPolygon,1e-8);
                if(triangulatedPolygon->size() < 3)
                    continue;
                triangulatedProjectedPolygons[_elemIndex][_patchIndex].push_back(*triangulatedPolygon);
                triangulatedProjectedPolygons2[_patchIndex].push_back(*triangulatedPolygon);
                // Get canonical element
                Polygon2D polygonWZ = computeCanonicalElement(_elemIndex, _projectedElement, *triangulatedPolygon);
                // Integrate
                integrate(thePatch, _patchIndex, *triangulatedPolygon, listSpan[index].first, listSpan[index].second, polygonWZ, _elemIndex);
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
    knotSpanWindow[0] = make_pair(u0,v0);
    knotSpanWindow[1] = make_pair(u1,v0);
    knotSpanWindow[2] = make_pair(u1,v1);
    knotSpanWindow[3] = make_pair(u0,v1);
    ClipperAdapter c;
    _polygonUV = c.clip(_polygonUV,knotSpanWindow);
}

void IGAMortarMapper::clipByTrimming(const IGAPatchSurface* _thePatch, const Polygon2D& _polygonUV, ListPolygon2D& _listPolygonUV) {
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
    for(int loop = 0; loop < _thePatch->getTrimming().getNumOfLoops(); loop++) {
        // 2i. Get the linearized trimming loop
        const std::vector<double> clippingWindow = _thePatch->getTrimming().getLoop(loop).getPolylines();

        // 2ii. Add the trimming loop in the clipper
        c.addPathClipper(clippingWindow);
    }

    // 3. Setup filling rule to have for sure clockwise loop as hole and counterclockwise as boundaries
    c.setFilling(ClipperAdapter::POSITIVE, 0);
    c.addPathSubject(_polygonUV);
    c.clip();
    c.getSolution(_listPolygonUV);
}

void IGAMortarMapper::clipByTrimming(const IGAPatchSurfaceTrimmingLoop* _theTrimmingLoop, const Polygon2D& _polygonUV, ListPolygon2D& _listPolygonUV) {
    /*
     * Clips a given polygon by the tirmming curves within a trimming loop
     *
     * Function layout:
     *
     * 1. Initialize auxiliary variables
     *
     * 2. Define the clipping window for the clipper
     *
     * 3. Setup filling rule to have for sure clockwise loop as hole and counterclockwise as boundaries
     */

    // 1. Initialize auxiliary variables
    ClipperAdapter c;
    ListPolygon2D tmpListPolygonUV;

    // 2. Define the clipping window for the clipper
    const std::vector<double> clippingWindow = _theTrimmingLoop->getPolylines();
    c.addPathClipper(clippingWindow);

    // 3. Setup filling rule to have for sure clockwise loop as hole and counterclockwise as boundaries
    c.setFilling(ClipperAdapter::NEGATIVE, 0);
    c.addPathSubject(_polygonUV);
    c.clip();
    c.getSolution(tmpListPolygonUV);
    _listPolygonUV.insert(_listPolygonUV.end(),tmpListPolygonUV.begin(),tmpListPolygonUV.end());
}

void IGAMortarMapper::clipByKnotSpan(const IGAPatchSurface* _thePatch, const Polygon2D& _polygonUV, ListPolygon2D& _listPolygon, Polygon2D& _listSpan) {
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

IGAMortarMapper::ListPolygon2D IGAMortarMapper::triangulatePolygon(const Polygon2D& _polygonUV) {
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

IGAMortarMapper::Polygon2D IGAMortarMapper::computeCanonicalElement(const int _elementIndex, const Polygon2D& _theElement, const Polygon2D& _polygonUV) {
    /*
     *  Finds the vertices of the clipped in the parameter space element in the canonical space and stores them the returned polygon
     *
     * Function layout:
     *
     * 1. Define auxiliary variables
     *
     * 2. Loop over all vertices of the clipped element in the parameter space of the patch and create a polygon
     *
     * 3. Loop over all vertices of the polygon in the parameter space of the patch
     * ->
     *    3i. Get the coordinates of a segment in the polygon
     *   3ii. Compute the local coordinates in the patch given the type of the canonical element
     *  3iii. Push back the coordinates of the element in the canonical space
     * <-
     *
     * 4. Return the newly created polygon containing the image of the clipped Finite Element in the canonical space
     */

    // 1. Define auxiliary variables
    int numNodesElementFE = meshFE->numNodesPerElem[_elementIndex];
    double elementFEUV[8];
    double coordsNodeFEUV[2];
    double coordsNodeFEWZ[2];
    Polygon2D polygonWZ;

    // 2. Loop over all vertices of the clipped element in the parameter space of the patch and create a polygon
    for(int iNodes = 0; iNodes < numNodesElementFE; iNodes++) {
        elementFEUV[2*iNodes] = _theElement[iNodes].first;
        elementFEUV[2*iNodes + 1] = _theElement[iNodes].second;
    }

    // 3. Loop over all vertices of the polygon in the parameter space of the patch
    for(int iVertices = 0; iVertices < _polygonUV.size(); iVertices++) {
        // 3i. Get the coordinates of a segment in the polygon
        coordsNodeFEUV[0] = _polygonUV[iVertices].first;
        coordsNodeFEUV[1] = _polygonUV[iVertices].second;

        // 3ii. Compute the local coordinates in the patch given the type of the canonical element
        if(numNodesElementFE == 3)
            MathLibrary::computeLocalCoordsInTriangle(elementFEUV, coordsNodeFEUV, coordsNodeFEWZ);
        else if (numNodesElementFE == 4)
            MathLibrary::computeLocalCoordsInQuad(elementFEUV, coordsNodeFEUV, coordsNodeFEWZ);
        else {
            ERROR_OUT() << "The canonical element can only be a triangle or a quadrilateral" << endl;
            exit(-1);
        }

        // 3iii. Push back the coordinates of the element in the canonical space
        polygonWZ.push_back(make_pair(coordsNodeFEWZ[0],coordsNodeFEWZ[1]));
    }

    // 4. Return the newly created polygon containing the image of the clipped Finite Element in the canonical space
    return polygonWZ;
}

void IGAMortarMapper::integrate(IGAPatchSurface* _thePatch, int _patchIndex, Polygon2D _polygonUV,
                                int _spanU, int _spanV, Polygon2D _polygonWZ, int _elementIndex) {
    /*
     * Compute and assembles the element contributions to the global Cnn and Cnr matrices
     *
     * Function layout:
     *
     * 1. Read input
     *
     * 2. Get the corresponding quadrature rule depending on the integration domain
     *
     * 3. Find the number of the master and the slave DOFs
     *
     * 4. Create an element freedom table for the isogeometric field
     *
     * 5. Copy input polygon into contiguous C format
     *
     * 6. Loop over all Gauss points
     * ->
     *    6i. Get the Gauss point coordinates in the integration space
     *   6ii. Get the Gauss point weight
     *  6iii. Compute the basis functions describing parametrically the integration domain
     *   6iv. Compute the image of the gauss point in the NURBS parameter space
     *    6v. Compute the image of the Gauss point in the low order Finite Element
     *   6vi. Compute the basis functions of the low order Finite Element
     *  6vii. Compute the local NURBS basis functions and their first order derivatives
     * 6viii. Compute the base vectors and their first derivatives
     *   6ix. Compute the surface normal vector at the Gauss point
     *    6x. Compute the determinant of the Jacobian of the transformation from the physical space to the NURBS parameter space
     *   6xi. Compute the determinant of the Jacobian of the transformation from the NURBS parameter space to the canonical space
     *  6xii. Compute the product of the determinants of the Jacobian matrices
     * 6xiii. Loop over all local basis functions in the master side
     * ->
     *        6xiii.1. Loop over all local basis functions in the master side to compute and assemble the local Cnn matrix to the global one
     *        ->
     *                 6xiii.1i. Compute the product of the basis functions
     *                6xiii.1ii. Compute the integrand on the Gauss point times the Gauss weight
     *               6xiii.1iii. Find the DOF numbering of the dual basis functions product
     *                6xiii.1iv. Assemble the element contributions to the global Cnn matrix
     *        <-
     *
     *        6xiii.2. Loop over all local basis functions in the slave side to compute and assemble the local Cnr matrix to the global one
     *        ->
     *                 6xiii.2i. Compute the integrand on the Gauss point times the Gauss weight
     *                6xiii.2ii. Find the DOF numbering of the dual basis functions product
     *               6xiii.2iii. Assemble the element contributions to the global Cnr matrix
     *        <-
     * <-
     *
     *  6xiv. Save the gauss point data for the computation of the L2 norm of the error
     * <-
     */

    // 1. Read input
    assert(!_polygonUV.empty());
    assert(!_polygonWZ.empty());
    int numNodesUV=_polygonUV.size();
    int numNodesWZ=_polygonWZ.size();
    assert(numNodesUV > 2);
    assert(numNodesUV < 5);
    assert(numNodesWZ > 2);
    assert(numNodesWZ < 5);

    // 2. Get the corresponding quadrature rule depending on the integration domain
    EMPIRE::MathLibrary::IGAGaussQuadrature* theGaussQuadrature;
    int nNodesQuadrature = numNodesUV;
    if (nNodesQuadrature == 3)
        theGaussQuadrature = gaussRuleOnTriangle[_patchIndex];
    else if (nNodesQuadrature == 4)
        theGaussQuadrature = gaussRuleOnQuadrilateral[_patchIndex];
    else {
        ERROR_OUT() << "Only triangles and quadrilaterals are expected as integration domains";
        exit(-1);
    }

    // 2. Initialize auxiliary variables
    int noCoord = 3;
    int derivDegree = 1;
    int numNodesElementFE = meshFE->numNodesPerElem[_elementIndex];
    int numNodesElMaster = 0;
    int numNodesElSlave = 0;
    int dof1;
    int dof2;
    int pDegree = _thePatch->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree();
    int qDegree = _thePatch->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree();
    int numBasisFunctionsIGA = (pDegree + 1) * (qDegree + 1);
    double basisFunctions[nNodesQuadrature];
    double basisFunctionsFE[numNodesElementFE];
    double uv[2];
    double wz[2];
    double baseVectors[6];
    double baseVectorU[3];
    double baseVectorV[3];
    double surfaceNormalTilde[3];
    double gaussWeight;
    double JacobianUVToPhysical;
    double JacobianCanonicalToUV;
    double JacobianProduct;
    double dudx;
    double dudy;
    double dvdx;
    double dvdy;
    double IGABasisFctsI;
    double IGABasisFctsJ;
    double basisFctsMaster;
    double basisFctsSlave;
    double basisFunctionsProduct;
    double integrand;
    double localBasisFunctionsAndDerivatives[(derivDegree + 1) * (derivDegree + 2) * numBasisFunctionsIGA / 2];

    // 3. Find the number of the master and the slave DOFs
    if (isMappingIGA2FEM) {
        numNodesElMaster = numNodesElementFE;
        numNodesElSlave = numBasisFunctionsIGA;
    } else {
        numNodesElMaster = numBasisFunctionsIGA;
        numNodesElSlave = numNodesElementFE;
    }

    // 4. Create an element freedom table for the isogeometric field
    int dofIGA[numBasisFunctionsIGA];
    _thePatch->getIGABasis()->getBasisFunctionsIndex(_spanU, _spanV, dofIGA);
    for (int i = 0; i < numBasisFunctionsIGA; i++)
        dofIGA[i] = _thePatch->getControlPointNet()[dofIGA[i]]->getDofIndex();

    // 5. Copy input polygon into contiguous C format
    double nodesUV[8];
    double nodesWZ[8];
    for (int i = 0; i < numNodesUV; i++) {
        nodesUV[i*2] = _polygonUV[i].first;
        nodesUV[i*2 + 1] = _polygonUV[i].second;
        nodesWZ[i*2] = _polygonWZ[i].first;
        nodesWZ[i*2 + 1] = _polygonWZ[i].second;
    }

    // 6. Loop over all Gauss points
    for (int iGP = 0; iGP < theGaussQuadrature->getNumGaussPoints(); iGP++) {
        // 6i. Get the Gauss point coordinates in the integration space
        const double *gaussPoint = theGaussQuadrature->getGaussPoint(iGP);

        // 6ii. Get the Gauss point weight
        gaussWeight = theGaussQuadrature->getGaussWeight(iGP);

        // 6iii. Compute the basis functions describing parametrically the integration domain
        MathLibrary::computeLowOrderShapeFunc(nNodesQuadrature, gaussPoint, basisFunctions);

        // 6iv. Compute the image of the gauss point in the NURBS parameter space
        MathLibrary::computeLinearCombination(nNodesQuadrature, 2, nodesUV, basisFunctions, uv);

        // 6v. Compute the image of the Gauss point in the low order Finite Element
        MathLibrary::computeLinearCombination(nNodesQuadrature, 2, nodesWZ, basisFunctions, wz);

        // 6vi. Compute the basis functions of the low order Finite Element
        MathLibrary::computeLowOrderShapeFunc(numNodesElementFE, wz, basisFunctionsFE);

        // 6vii. Compute the local NURBS basis functions and their first order derivatives
        _thePatch->getIGABasis()->computeLocalBasisFunctionsAndDerivatives(localBasisFunctionsAndDerivatives,
                                                                           derivDegree, uv[0], _spanU, uv[1], _spanV);

        // 6viii. Compute the base vectors and their first derivatives
        _thePatch->computeBaseVectors(baseVectors, localBasisFunctionsAndDerivatives, _spanU, _spanV);
        for (int iCoord = 0; iCoord < noCoord; iCoord++) {
            baseVectorU[iCoord] = baseVectors[iCoord];
            baseVectorV[iCoord] = baseVectors[noCoord + 1 + iCoord];
        }

        // 6ix. Compute the surface normal vector at the Gauss point
        MathLibrary::computeVectorCrossProduct(baseVectorU, baseVectorV, surfaceNormalTilde);

        // 6x. Compute the determinant of the Jacobian of the transformation from the physical space to the NURBS parameter space
        JacobianUVToPhysical = MathLibrary::vector2norm(surfaceNormalTilde, noCoord);

        // 6xi. Compute the determinant of the Jacobian of the transformation from the NURBS parameter space to the canonical space
        if (nNodesQuadrature == 3) {
            JacobianCanonicalToUV = MathLibrary::computeAreaTriangle( nodesUV[2] - nodesUV[0], nodesUV[3] - nodesUV[1], 0,
                                                                      nodesUV[4] - nodesUV[0], nodesUV[5] - nodesUV[1], 0)*2.0;;
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

        // 6xii. Compute the product of the determinants of the Jacobian matrices
        JacobianProduct = JacobianUVToPhysical * JacobianCanonicalToUV;

        // 6xiii. Loop over all local basis functions in the master side
        for (int i = 0; i < numNodesElMaster; i++) {
            // 6xiii.1. Loop over all local basis functions in the master side to compute and assemble the local Cnn matrix to the global one
            for (int j = i; j < numNodesElMaster; j++) { // Starts from i because of computing only the upper triangular entries of the matrix
                // 6xiii.1i. Compute the product of the basis functions
                if (isMappingIGA2FEM)
                    basisFunctionsProduct = basisFunctionsFE[i] * basisFunctionsFE[j];
                else {
                    IGABasisFctsI = localBasisFunctionsAndDerivatives[_thePatch->getIGABasis()->indexDerivativeBasisFunction(1, 0, 0, i)];
                    IGABasisFctsJ = localBasisFunctionsAndDerivatives[_thePatch->getIGABasis()->indexDerivativeBasisFunction(1, 0, 0, j)];
                    basisFunctionsProduct = IGABasisFctsI*IGABasisFctsJ;
                }

                // 6xiii.1ii. Compute the integrand on the Gauss point times the Gauss weight
                integrand = basisFunctionsProduct*JacobianProduct*theGaussQuadrature->getGaussWeight(iGP);

                // 6xiii.1iii. Find the DOF numbering of the dual basis functions product
                if (isMappingIGA2FEM) {
                    dof1 = meshFEDirectElemTable[_elementIndex][i];
                    dof2 = meshFEDirectElemTable[_elementIndex][j];
                } else {
                    dof1 = dofIGA[i];
                    dof2 = dofIGA[j];
                }

                // 6xiii.1iv. Assemble the element contributions to the global Cnn matrix
                if (!isExpanded){
                    couplingMatrices->addCNNValue(dof1, dof2, integrand);
                    if (dof1 != dof2)
                        couplingMatrices->addCNNValue(dof2, dof1, integrand);
                } else {
                    for(int iCoord = 0; iCoord < noCoord; iCoord++){
                        couplingMatrices->addCNNValue(noCoord*dof1 + iCoord, noCoord*dof2 + iCoord, integrand);
                        if (dof1 != dof2)
                            couplingMatrices->addCNNValue(noCoord*dof2 + iCoord, noCoord*dof1 + iCoord, integrand);
                    }
                }
            }

            // 6xiii.2. Loop over all local basis functions in the slave side to compute and assemble the local Cnr matrix to the global one
            for (int j = 0; j < numNodesElSlave; j++) {
                // 6xiii.2i. Compute the integrand on the Gauss point times the Gauss weight
                if (isMappingIGA2FEM) {
                    basisFctsMaster = basisFunctionsFE[i];
                    basisFctsSlave = localBasisFunctionsAndDerivatives[_thePatch->getIGABasis()->indexDerivativeBasisFunction(1, 0, 0, j)];
                } else {
                    basisFctsMaster = localBasisFunctionsAndDerivatives[_thePatch->getIGABasis()->indexDerivativeBasisFunction(1, 0, 0, i)];
                    basisFctsSlave = basisFunctionsFE[j];
                }
                integrand = basisFctsMaster*basisFctsSlave*JacobianProduct*theGaussQuadrature->getGaussWeight(iGP);

                // 6xiii.2ii. Find the DOF numbering of the dual basis functions product
                if (isMappingIGA2FEM) {
                    dof1 = meshFEDirectElemTable[_elementIndex][i];
                    dof2 = dofIGA[j];
                } else {
                    dof1 = dofIGA[i];
                    dof2 = meshFEDirectElemTable[_elementIndex][j];
                }

                // 6xiii.2iii. Assemble the element contributions to the global Cnr matrix
                if (!isExpanded){
                    couplingMatrices->addCNRValue(dof1, dof2, integrand);
                } else {
                    for(int iCoord = 0; iCoord < noCoord; iCoord++)
                        couplingMatrices->addCNRValue(noCoord*dof1 + iCoord, noCoord*dof2 + iCoord, integrand);
                }
            }
        }

        // 6xiv. Save the gauss point data for the computation of the L2 norm of the error
        if(propErrorComputation.isDomainError){
            std::vector<double> streamGP;
            // weight + JacobianProduct + numBasisFuncsFE + (#dof, shapefuncvalue,...) + nShapeFuncsIGA + (#dof, shapefuncvalue,...)
            streamGP.reserve(1 + 1 + 1 + 2*numNodesElementFE + 1 + 2*numBasisFunctionsIGA);
            streamGP.push_back(theGaussQuadrature->getGaussWeight(iGP));
            streamGP.push_back(JacobianProduct);
            streamGP.push_back(numNodesElementFE);
            for (int i = 0; i < numNodesElementFE; i++) {
                streamGP.push_back(meshFEDirectElemTable[_elementIndex][i]);
                streamGP.push_back(basisFunctionsFE[i]);
            }
            streamGP.push_back(numBasisFunctionsIGA);
            for (int i = 0; i < numBasisFunctionsIGA; i++) {
                double IGABasisFctsI = localBasisFunctionsAndDerivatives[_thePatch->getIGABasis()->indexDerivativeBasisFunction(1, 0, 0, i)];
                streamGP.push_back(dofIGA[i]);
                streamGP.push_back(IGABasisFctsI);
            }
            streamGPs.push_back(streamGP);
        }
    }
}

bool IGAMortarMapper::computeKnotSpanOfProjElement(const IGAPatchSurface* _thePatch, const Polygon2D& _polygonUV, int* _span) {
    /*
     * Returns the knot span bounds in terms of knot span indices where the given polygon belongs to
     *
     * Function layout:
     *
     * 1. Initialize auxiliary variables
     *
     * 2. Loop over all the segments of the given polygon
     * ->
     *    2i. Get the (u,v) parameters of the vertex of the polygon
     *   2ii. Update the bounding values accordingly
     * <-
     *
     * 3. flag on whether all vertices of the polygon are found in the same knot span
     *
     * 4. Fill up the span pointer with the span bounding values
     *
     * 5. Return the flag
     */

    // 1. Initialize auxiliary variables
    int minSpanU = _thePatch->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(
                _polygonUV[0].first);
    int minSpanV = _thePatch->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(
                _polygonUV[0].second);
    int maxSpanU = minSpanU;
    int maxSpanV = minSpanV;

    // 2. Loop over all the segments of the given polygon
    for (int nodeCount = 1; nodeCount < _polygonUV.size(); nodeCount++) {
        // 2i. Get the (u,v) parameters of the vertex of the polygon
        int spanU = _thePatch->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(_polygonUV[nodeCount].first);
        int spanV = _thePatch->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(_polygonUV[nodeCount].second);

        // 2ii. Update the bounding values accordingly
        if (spanU < minSpanU)
            minSpanU = spanU;
        if (spanU > maxSpanU)
            maxSpanU = spanU;
        if (spanV < minSpanV)
            minSpanV = spanV;
        if (spanV > maxSpanV)
            maxSpanV = spanV;
    }

    // 3. flag on whether all vertices of the polygon are found in the same knot span
    bool OnSameKnotSpan = (minSpanU == maxSpanU && minSpanV == maxSpanV);

    // 4. Fill up the span pointer with the span bounding values
    if(_span != NULL) {
        _span[0] = minSpanU;
        _span[1] = maxSpanU;
        _span[2] = minSpanV;
        _span[3] = maxSpanV;
    }

    // 5. Return the flag
    return OnSameKnotSpan;
}

int IGAMortarMapper::getNeighbourElementofEdge(int _element, int _node1, int _node2) {
    /*
     * Return the index of the element neighbouring the segment defined by thw two given nodes
     *
     * Function layout:
     *
     * 1. Initialize auxilary variables
     *
     * 2. Loop over all the elements in the Finite Element mesh
     * ->
     *    2i. Initialize the flags to false
     *   2ii. Loop over all nodes of the element
     *   ->
     *        2ii.1. Check if first vertex is found in the given element
     *        2ii.2. Check if first vertex is found in the given element
     *        2ii.3. If both vertices are found return the index of the element
     *   <-
     * <-
     *
     * 3. If polygon is on boundary of mesh, can occur
     */

    // 1. Initialize auxilary variables
    bool isNode1;
    bool isNode2;

    // 2. Loop over all the elements in the Finite Element mesh
    for(int i = 0; i < meshFE->numElems; i++) {
        // 2i. Initialize the flags to false
        isNode1 = false;
        isNode2 = false;

        // 2ii. Loop over all nodes of the element
        for(int j = 0; j < meshFE->numNodesPerElem[i]; j++) {
            // 2ii.1. Check if first vertex is found in the given element
            if (isNode1 == false)
                isNode1 = (meshFEDirectElemTable[i][j] == _node1)?true:false;

            // 2ii.2. Check if first vertex is found in the given element
            if (isNode2 == false)
                isNode2 = (meshFEDirectElemTable[i][j] == _node2)?true:false;

            // 2ii.3. If both vertices are found return the index of the element
            if(_element != i && isNode1 && isNode2) {
                return i;
            }
        }
    }

    // 3. If polygon is on boundary of mesh, can occur
    return -1;
}

void IGAMortarMapper::createGaussQuadratureRules() {
    /*
     * Creates a Gauss quadrature rule for quadrilaterals and for triangles at each patch
     *
     * Function layout:
     *
     * 1. Initialize auxiliary variables
     *
     * 2. Set flag of Gauss quadratures to true (needed for the destructor of the class)
     *
     * 3. Get the number of the patches
     *
     * 4. Loop over all patches
     * ->
     *    4i. Get the polynomial order of the patch
     *   4ii. Find the polynomial degree of the integrands when a triangle is considered
     *  4iii. Find the number of Gauss points when a triangle is considered
     *   4iv. Instantiate the corresponding Gauss quadrature on triangle
     *    4v. Find the polynomial degree of the integrands when a quadrilateral is considered
     *   4vi. Compute the number of Gauss points when a quadrilateral is considered
     *  4vii. Instantiate the corresponding Gauss quadrature on quadrilateral
     * <-
     */

    // 1. Initialize auxiliary variables
    int pDegree;
    int qDegree;
    int polOrder;
    int pTilde;
    int numGPs;
    int numGPsPerPolOrder[8] = {1, 3, 4, 6, 7, 12, 13, 16};

    // 2. Set flag of Gauss quadratures to true (needed for the destructor of the class)
    isGaussQuadature = true;

    // 3. Get the number of the patches
    int numPatches = getIGAMesh()->getNumPatches();

    // 4. Loop over all patches
    for (int iPatches = 0; iPatches < numPatches; iPatches++) {
        // 4i. Get the polynomial order of the patch
        pDegree = meshIGA->getSurfacePatch(iPatches)->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree();
        qDegree = meshIGA->getSurfacePatch(iPatches)->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree();

        // 4ii. Find the polynomial degree of the integrands when a triangle is considered
        if (propIntegration.isAutomaticNoGPTriangle)  { // Automatic definition of the quadrature rule
            if (isMappingIGA2FEM) // When mapping from IGA to FEM integrals \int_{\Omega} N_i*N_i d\Omega and \int_{\Omega} N_i*R_i d\Omega are computed within the same loop
                if ((pDegree == 0 && qDegree == 0) || (pDegree == 1 && qDegree == 0) || (pDegree == 0 && qDegree == 1)) { // Integrand N_i*N_i has higher polynomial order than N_i*R_i
                    polOrder = 2;
                    if (polOrder > 8)
                        pTilde = 1;
                } else { // Integrand N_i*R_i has higher polynomial order than N_i*N_i
                    polOrder = 1 + pDegree + qDegree;
                    if (polOrder > 8)
                        pTilde = std::max(pDegree, qDegree) + 1;
                }
            else // When mapping from FEM to IGA integrals \int_{\Omega} R_i*R_i d\Omega and \int_{\Omega} R_i*N_i d\Omega are computed within the same loop
                if ((pDegree == 0 && qDegree == 0) || (pDegree == 1 && qDegree == 0) || (pDegree == 0 && qDegree == 1)) { // Integrand N_i*R_i has higher polynomial order than R_i*R_i
                    polOrder = pDegree + qDegree + 1;
                    if (polOrder > 8)
                        pTilde = std::max(pDegree, qDegree) + 1;
                } else { // Integrand R_i*R_i has higher polynomial order than N_i*R_i
                    polOrder = 2*(pDegree + qDegree);
                    if (polOrder > 8)
                        pTilde = 2*std::max(pDegree, qDegree);
                }

        // 4iii. Find the number of Gauss points when a triangle is considered
        if (propIntegration.isAutomaticNoGPTriangle)
            if (polOrder <= 8) { // Use the Gauss quadrature over the triangle with the symmetric rule
                if (polOrder == 0)
                    numGPs = numGPsPerPolOrder[0];
                else
                    numGPs = numGPsPerPolOrder[polOrder - 1];
            } else // Use the Gauss quadrature over the triangle with the degenerated quadrilateral
                numGPs = pow(std::ceil((pTilde + 1)/2.0), 2.0);
            }
        else if (!propIntegration.isAutomaticNoGPTriangle)
            numGPs = propIntegration.noGPTriangle;
        else
            ERROR_BLOCK_OUT("createGaussQuadratureRules","IGAMortarMapper","Found corner case not encountered before related to the integration over a triangle!");

        // 4iv. Instantiate the corresponding Gauss quadrature on triangle
        if ((propIntegration.isAutomaticNoGPTriangle && polOrder <= 8) || !propIntegration.isAutomaticNoGPTriangle && numGPs <= 16) // Use the Gauss quadrature over the triangle with the symmetric rule
            gaussRuleOnTriangle[iPatches] = new MathLibrary::IGAGaussQuadratureOnTriangle(numGPs);
        else if ((propIntegration.isAutomaticNoGPTriangle && polOrder > 8) || !propIntegration.isAutomaticNoGPTriangle && numGPs > 16) // Use the Gauss quadrature over the triangle with the degenerated quadrilateral
            gaussRuleOnTriangle[iPatches] = new MathLibrary::IGAGaussQuadratureOnTriangleUsingDegeneratedQuadrilateral(numGPs);
        else
            ERROR_BLOCK_OUT("createGaussQuadratureRules","IGAMortarMapper","Found corner case not encountered before related to the integration over a triangle!");

        // 4v. Find the polynomial degree of the integrands when a quadrilateral is considered
        if (propIntegration.isAutomaticNoGPQuadrilateral) // Automatic definition of the quadrature rule
            if (isMappingIGA2FEM) // When mapping from IGA to FEM integrals \int_{\Omega} N_i*N_i d\Omega and \int_{\Omega} N_i*R_i d\Omega are computed within the same loop
                if ((pDegree == 0 && qDegree == 0) || (pDegree == 1 && qDegree == 0) || (pDegree == 0 && qDegree == 1)) { // Integrand N_i*N_i has higher polynomial order than N_i*R_i
                    pTilde = 2;
                } else { // Integrand N_i*R_i has higher polynomial order than N_i*N_i
                    pTilde = std::max(pDegree, qDegree) + 2;
                }
            else // When mapping from FEM to IGA integrals \int_{\Omega} R_i*R_i d\Omega and \int_{\Omega} R_i*N_i d\Omega are computed within the same loop
                if ((pDegree == 0 && qDegree == 0) || (pDegree == 1 && qDegree == 0) || (pDegree == 0 && qDegree == 1)) { // Integrand R_i*N_i has higher polynomial order than R_i*R_i
                    pTilde = std::max(pDegree, qDegree) + 2;
                } else {
                    pTilde = 2*std::max(pDegree, qDegree);
                }

        // 4vi. Compute the number of Gauss points when a quadrilateral is considered
        if (propIntegration.isAutomaticNoGPQuadrilateral)
            numGPs = pow(std::ceil((pTilde + 1)/2.0), 2.0);
        else if (!propIntegration.isAutomaticNoGPQuadrilateral)
            numGPs = propIntegration.noGPQuadrilateral;
        else
            ERROR_BLOCK_OUT("createGaussQuadratureRules","IGAMortarMapper","Found corner case not encountered before related to the integration over a quadrilateral!");

        // 4vii. Instantiate the corresponding Gauss quadrature on quadrilateral
        gaussRuleOnQuadrilateral[iPatches] = new MathLibrary::IGAGaussQuadratureOnBiunitQuadrilateral(numGPs);
    }
}

void IGAMortarMapper::computeIGAWeakDirichletCurveConditionMatrices() {
    /*
     * Computes and assembles the patch weak Dirichlet curve conditions' contributions.
     *
     * Function layout:
     *
     * 1. Initialize auxiliary variables
     *
     * 2. Get the weak Dirichlet curve conditions
     *
     * 3. Loop over all the conditions for the application of weak Dirichlet conditions
     * ->
     *    3i. Get the penalty factors for the primary and the secondary field
     *   3ii. Get the index of the patch
     *  3iii. Get the number of Gauss Points for the given condition
     *   3iv. Get the parametric coordinates of the Gauss Points
     *    3v. Get the corresponding Gauss weights
     *   3vi. Get the tangent vectors at the trimming curve of the given condition in the Cartesian space
     *  3vii. Get the product of the Jacobian transformations
     * 3viii. Get the patch
     *   3ix. Get the polynomial orders of the master and the slave patch
     *    3x. get the number of local basis functions
     *   3xi. get the number of the local DOFs for the patch
     *  3xii. Initialize pointers
     * 3xiii. Loop over all the Gauss Points of the given condition
     * ->
     *        3xiii.1. Get the parametric coordinates of the Gauss Point on the patch
     *        3xiii.2. Find the knot span indices of the Gauss point locations in the parameter space of the patch
     *        3xiii.3. Get the tangent to the boundary vector on the patch
     *        3xiii.4. compute elementLength on GP. The weight is already included in variable trCurveGPJacobianProducts
     *        3xiii.5. Compute the B-operator matrices needed for the computation of the patch weak Dirichlet conditions at the patch
     *        3xiii.6. Compute the local penalty factor for scaling the rotational contributions
     *        3xiii.7. Scale the B-operator matrices
     *        3xiii.8. Compute the dual product matrices for the displacements
     *        3xiii.9. Compute the dual product matrices for the bending rotations
     *       3xiii.10. Compute the dual product matrices for the twisting rotations
     *       3xiii.11. Compute the element index tables for the patch
     *       3xiii.12. Compute the element freedom tables for the patch
     *       3xiii.13. Loop over all DOFs to assemble the contributions of the weak Dirichlet conditions into the Cnn
     *       ->
     *                 3xiii.13i. Assemble the displacement coupling entries
     *                3xiii.13ii. Assemble the bending rotation coupling entries
     *               3xiii.13iii. Assemble the twisting rotation coupling entries
     *       <-
     *       3xiii.14. Store the Gauss point values necessary for the error computation
     * <-
     *  3xiv. Delete pointers
     * <-
     */

    // 1. Initialize auxiliary variables
    const int noCoordParam = 2;
    const int noCoord = 3;
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
    double normBOperatorOmegaT;
    double normBOperatorOmegaN;
    double surfaceNormalVct[noCoord];
    double alphaLocal;
    double alphaPrimary;
    double alphaSecondaryBending;
    double alphaSecondaryTwisting;
    double elementLengthOnGP;
    double* curveGPs;
    double* curveGPWeights;
    double* curveGPTangents;
    double* curveGPJacobianProducts;
    IGAPatchSurface* thePatch;

    // 2. Get the weak Dirichlet curve conditions
    std::vector<WeakIGADirichletCurveCondition*> weakIGADirichletCurveConditions = meshIGA->getWeakIGADirichletCurveConditions();

    // 3. Loop over all the conditions for the application of weak Dirichlet conditions
    for (int iDCC = 0; iDCC < weakIGADirichletCurveConditions.size(); iDCC++){
        // 3i. Get the penalty factors for the primary and the secondary field
        alphaPrimary = weakDirichletCCAlphaPrimary[iDCC];
        alphaSecondaryBending = weakDirichletCCAlphaSecondaryBending[iDCC];
        alphaSecondaryTwisting = weakDirichletCCAlphaSecondaryTwisting[iDCC];

        // 3ii. Get the index of the patch
        patchIndex = weakIGADirichletCurveConditions[iDCC]->getPatchIndex();

        // 3iii. Get the number of Gauss Points for the given condition
        noGPsOnCond = weakIGADirichletCurveConditions[iDCC]->getCurveNumGP();

        // 3iv. Get the parametric coordinates of the Gauss Points
        curveGPs = weakIGADirichletCurveConditions[iDCC]->getCurveGPs();

        // 3v. Get the corresponding Gauss weights
        curveGPWeights = weakIGADirichletCurveConditions[iDCC]->getCurveGPWeights();

        // 3vi. Get the tangent vectors at the trimming curve of the given condition in the Cartesian space
        curveGPTangents = weakIGADirichletCurveConditions[iDCC]->getCurveGPTangents();

        // 3vii. Get the product of the Jacobian transformations
        curveGPJacobianProducts = weakIGADirichletCurveConditions[iDCC]->getCurveGPJacobianProducts();

        // 3viii. Get the patch
        thePatch = meshIGA->getSurfacePatch(patchIndex);

        // 3ix. Get the polynomial orders of the master and the slave patch
        p = thePatch->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree();
        q = thePatch->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree();

        // 3x. get the number of local basis functions
        noLocalBasisFcts = (p + 1)*(q + 1);

        // 3xi. get the number of the local DOFs for the patch
        noDOFsLoc = noCoord*noLocalBasisFcts;

        // 3xii. Initialize pointers
        double* BDisplacementsGC = new double[noCoord*noDOFsLoc];
        double* BOperatorOmegaT = new double[noDOFsLoc];
        double* BOperatorOmegaN = new double[noDOFsLoc];
        double* KPenaltyDisplacement = new double[noDOFsLoc*noDOFsLoc];
        double* KPenaltyBendingRotation = new double[noDOFsLoc*noDOFsLoc];
        double* KPenaltyTwistingRotation = new double[noDOFsLoc*noDOFsLoc];

        // 3xiii. Loop over all the Gauss Points of the given condition
        for(int iGP = 0; iGP < noGPsOnCond; iGP++){
            // 3xiii.1. Get the parametric coordinates of the Gauss Point on the patch
            uGP = curveGPs[iGP*noCoordParam];
            vGP = curveGPs[iGP*noCoordParam + 1];

            // 3xiii.2. Find the knot span indices of the Gauss point locations in the parameter space of the patch
            uKnotSpan = thePatch->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(uGP);
            vKnotSpan = thePatch->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(vGP);

            // 3xiii.3. Get the tangent to the boundary vector on the patch
            for(int iCoord = 0; iCoord < noCoord; iCoord++)
                tangentCurveVct[iCoord] = curveGPTangents[iGP*noCoord + iCoord];

            // 3xiii.4. compute elementLength on GP. The weight is already included in variable trCurveGPJacobianProducts
            elementLengthOnGP = curveGPJacobianProducts[iGP];

            // 3xiii.5. Compute the B-operator matrices needed for the computation of the patch weak Dirichlet conditions at the patch
            computeDisplacementAndRotationBOperatorMatrices(BDisplacementsGC, BOperatorOmegaT, BOperatorOmegaN, normalCurveVct,
                                                            surfaceNormalVct, thePatch, tangentCurveVct, uGP, vGP, uKnotSpan, vKnotSpan);

            // 3xiii.6. Compute the local penalty factor for scaling the rotational contributions
            if (propWeakCurveDirichletConditions.isSecBendingPrescribed) {
                normBOperatorOmegaT = EMPIRE::MathLibrary::vector2norm(BOperatorOmegaT, noDOFsLoc);
                alphaLocal = normBOperatorOmegaT;
            } else
                alphaLocal = - 1.0;

            if (propWeakCurveDirichletConditions.isSecTwistingPrescribed) {
                normBOperatorOmegaN = EMPIRE::MathLibrary::vector2norm(BOperatorOmegaN, noDOFsLoc);
                if (normBOperatorOmegaN > alphaLocal)
                    alphaLocal = normBOperatorOmegaN;
            }
            alphaLocal = 1.0/abs(alphaLocal);

            // 3xiii.7. Scale the B-operator matrices
            EMPIRE::MathLibrary::computeDenseVectorMultiplicationScalar(BOperatorOmegaT, alphaLocal, noDOFsLoc);
            EMPIRE::MathLibrary::computeDenseVectorMultiplicationScalar(BOperatorOmegaN, alphaLocal, noDOFsLoc);

            // 3xiii.6. Compute the dual product matrices for the displacements
            EMPIRE::MathLibrary::computeTransposeMatrixProduct(noCoord,noDOFsLoc,noDOFsLoc,BDisplacementsGC,BDisplacementsGC,KPenaltyDisplacement);

            // 3xiii.7. Compute the dual product matrices for the bending rotations
            EMPIRE::MathLibrary::computeTransposeMatrixProduct(1, noDOFsLoc, noDOFsLoc, BOperatorOmegaT, BOperatorOmegaT, KPenaltyBendingRotation);

            // 3xiii.8. Compute the dual product matrices for the twisting rotations
            EMPIRE::MathLibrary::computeTransposeMatrixProduct(1, noDOFsLoc, noDOFsLoc, BOperatorOmegaN, BOperatorOmegaN, KPenaltyTwistingRotation);

            // 3xiii.9. Compute the element index tables for the patch
            int CPIndex[noLocalBasisFcts];
            thePatch->getIGABasis()->getBasisFunctionsIndex(uKnotSpan, vKnotSpan, CPIndex);

            // 3xiii.10. Compute the element freedom tables for the patch
            int EFT[noDOFsLoc];
            counter = 0;
            for (int i = 0; i < noLocalBasisFcts; i++){
                indexCP = thePatch->getControlPointNet()[CPIndex[i]]->getDofIndex();
                for (int j = 0; j < noCoord; j++){
                    EFT[counter] = noCoord*indexCP + j;
                    counter++;
                }
            }

            // 3xiii.11. Loop over all DOFs to assemble the contributions of the weak Dirichlet conditions into the Cnn
            for(int i = 0; i < noDOFsLoc; i++){
                for(int j = 0; j < noDOFsLoc; j++){
                    // 3xiii.11i. Assemble the displacement coupling entries
                    if (propWeakCurveDirichletConditions.isPrimPrescribed)
                        couplingMatrices->addCNNValue(EFT[i], EFT[i], alphaPrimary*KPenaltyDisplacement[i*noDOFsLoc + i]*elementLengthOnGP);

                    // 3xiii.11ii. Assemble the bending rotation coupling entries
                    if (propWeakCurveDirichletConditions.isSecBendingPrescribed)
                        couplingMatrices->addCNNValue(EFT[i], EFT[j], alphaSecondaryBending*KPenaltyBendingRotation[i*noDOFsLoc + j]*elementLengthOnGP);

                    // 3xiii.11iii. Assemble the twisting rotation coupling entries
                    if (propWeakCurveDirichletConditions.isSecTwistingPrescribed)
                        couplingMatrices->addCNNValue(EFT[i], EFT[j], alphaSecondaryTwisting*KPenaltyTwistingRotation[i*noDOFsLoc + j]*elementLengthOnGP);
                }
            }

            // 3xiii.12. Store the Gauss point values necessary for the error computation
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

        // 3xiv. Delete pointers
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
     *
     * Function layout:
     *
     * 1. Initialize auxiliary variables
     *
     * 2. Get the weak Dirichlet surface conditions
     *
     * 3. Loop over all the conditions for the application of weak Dirichlet conditions over surfaces
     * ->
     *    3i. Get the penalty factors for the primary field
     *   3ii. Get the index of the patch
     *  3iii. Get the number of Gauss Points for the given condition
     *   3iv. Get the parametric coordinates of the Gauss Points
     *    3v. Get the corresponding Gauss weights
     *   3vi. Get the Jacobian products
     *  3vii. Get the patch
     * 3viii. Get the polynomial orders of the master and the slave patch
     *   3ix. get the number of local basis functions
     *    3x. get the number of the local DOFs for the patch
     *   3xi. Initialize pointers
     *  3xii. Loop over all the Gauss Points of the given condition
     *  ->
     *        3xii.1. Get the parametric coordinates of the Gauss Point on the patch
     *        3xii.2. Find the knot span indices of the Gauss point locations in the parameter space of the patch
     *        3xii.3. compute elementLength on GP. The weight is already included in variable trCurveGPJacobianProducts
     *        3xii.4. Compute the B-operator matrices needed for the computation of the patch weak Dirichlet conditions at the patch
     *        3xii.5. Compute the dual product matrices for the displacements
     *        3xii.6. Compute the element index tables for the patch
     *        3xii.7. Compute the element freedom tables for the patch
     *        3xii.8. Assemble KPenaltyDisplacement to the global coupling matrix CNN
     *  <-
     * 3xiii. Delete pointers
     * <-
     */

    // 1. Initialize auxiliary variables
    const int noCoordParam = 2;
    const int noCoord = 3;
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
    double surfaceNormalVct[noCoord];
    double alphaPrimary;
    double jacobianOnGP;
    double* surfaceGPs;
    double* surfaceGPWeights;
    double* surfaceGPJacobians;
    IGAPatchSurface* thePatch;

    // 2. Get the weak Dirichlet surface conditions
    std::vector<WeakIGADirichletSurfaceCondition*> weakIGADirichletSurfaceConditions = meshIGA->getWeakIGADirichletSurfaceConditions();

    // 3. Loop over all the conditions for the application of weak Dirichlet conditions over surfaces
    for (int iDSC = 0; iDSC < weakIGADirichletSurfaceConditions.size(); iDSC++){
        // 3i. Get the penalty factors for the primary field
        alphaPrimary = weakDirichletSCAlphaPrimary[iDSC];

        // 3ii. Get the index of the patch
        patchIndex = weakIGADirichletSurfaceConditions[iDSC]->getPatchIndex();

        // 3iii. Get the number of Gauss Points for the given condition
        noGPsOnCond = weakIGADirichletSurfaceConditions[iDSC]->getSurfaceNumGP();

        // 3iv. Get the parametric coordinates of the Gauss Points
        surfaceGPs = weakIGADirichletSurfaceConditions[iDSC]->getSurfaceGPs();

        // 3v. Get the corresponding Gauss weights
        surfaceGPWeights = weakIGADirichletSurfaceConditions[iDSC]->getSurfaceGPWeights();

        // 3vi. Get the Jacobian products
        surfaceGPJacobians = weakIGADirichletSurfaceConditions[iDSC]->getSurfaceGPJacobians();

        // 3vii. Get the patch
        thePatch = meshIGA->getSurfacePatch(patchIndex);

        // 3viii. Get the polynomial orders of the master and the slave patch
        p = thePatch->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree();
        q = thePatch->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree();

        // 3ix. get the number of local basis functions
        noLocalBasisFcts = (p + 1)*(q + 1);

        // 3x. get the number of the local DOFs for the patch
        noDOFsLoc = noCoord*noLocalBasisFcts;

        // 3xi. Initialize pointers
        double* BDisplacementsGC = new double[noCoord*noDOFsLoc];
        double* BOperatorOmegaT = new double[noDOFsLoc];
        double* BOperatorOmegaN = new double[noDOFsLoc];
        double* KPenaltyDisplacement = new double[noDOFsLoc*noDOFsLoc];

        // 3xii. Loop over all the Gauss Points of the given condition
        for(int iGP = 0; iGP < noGPsOnCond; iGP++){
            // 3xii.1. Get the parametric coordinates of the Gauss Point on the patch
            uGP = surfaceGPs[iGP*noCoordParam];
            vGP = surfaceGPs[iGP*noCoordParam + 1];

            // 3xii.2. Find the knot span indices of the Gauss point locations in the parameter space of the patch
            uKnotSpan = thePatch->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(uGP);
            vKnotSpan = thePatch->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(vGP);

            // 3xii.3. compute elementLength on GP. The weight is already included in variable trCurveGPJacobianProducts
            jacobianOnGP = surfaceGPJacobians[iGP];

            // 3xii.4. Compute the B-operator matrices needed for the computation of the patch weak Dirichlet conditions at the patch
            computeDisplacementAndRotationBOperatorMatrices(BDisplacementsGC, BOperatorOmegaT, BOperatorOmegaN, normalCurveVct,
                                                            surfaceNormalVct, thePatch, tangentCurveVct, uGP, vGP, uKnotSpan, vKnotSpan);

            // 3xii.5. Compute the dual product matrices for the displacements
            if (propWeakSurfaceDirichletConditions.isPrimPrescribed)
                EMPIRE::MathLibrary::computeTransposeMatrixProduct(noCoord,noDOFsLoc,noDOFsLoc,BDisplacementsGC,BDisplacementsGC,KPenaltyDisplacement);

            // 3xii.6. Compute the element index tables for the patch
            int CPIndex[noLocalBasisFcts];
            thePatch->getIGABasis()->getBasisFunctionsIndex(uKnotSpan, vKnotSpan, CPIndex);

            // 3xii.7. Compute the element freedom tables for the patch
            int EFT[noDOFsLoc];
            counter = 0;
            for (int i = 0; i < noLocalBasisFcts; i++){
                indexCP = thePatch->getControlPointNet()[CPIndex[i]]->getDofIndex();
                for (int j = 0; j < noCoord; j++){
                    EFT[counter] = noCoord*indexCP + j;
                    counter++;
                }
            }

            // 3xii.8. Assemble KPenaltyDisplacement to the global coupling matrix CNN
            if (propWeakSurfaceDirichletConditions.isPrimPrescribed)
                for(int i = 0; i < noDOFsLoc; i++){
                    for(int j = 0; j < noDOFsLoc; j++){
                        // Assemble the displacement coupling entries
                        couplingMatrices->addCNNValue(EFT[i], EFT[i], alphaPrimary*KPenaltyDisplacement[i*noDOFsLoc + i]*jacobianOnGP);
                    }
                }

            //// 3xii.9. TODO Store the GP data into array for later usage in the error computation
        } // End of Gauss Point loop

        // 3xiii. Delete pointers
        delete[] BDisplacementsGC;
        delete[] BOperatorOmegaT;
        delete[] BOperatorOmegaN;
        delete[] KPenaltyDisplacement;

    } // End of weak Dirichlet surface condition loop
}

void IGAMortarMapper::computeIGAPatchWeakContinuityConditionMatrices() {
    /*
     * Computes and assembles the patch weak continuity conditions.
     *
     * Function layout:
     *
     * 1. Initialize auxiliary variables
     *
     * 2. Get the weak patch continuity conditions
     *
     * 3. Loop over all the conditions for the application of weak continuity across patch interfaces
     * ->
     *    3i. Get the penalty factors for the primary and the secondary field
     *   3ii. Get the index of the master and slave patches
     *  3iii. Get the number of Gauss Points for the given condition
     *   3iv. Get the parametric coordinates of the Gauss Points
     *    3v. Get the corresponding Gauss weights
     *   3vi. Get the tangent vectors at the trimming curve of the given condition in the Cartesian space
     *  3vii. Get the product of the Jacobian transformations
     * 3viii. Get the master and the slave patch
     *   3ix. Get the polynomial orders of the master and the slave patch
     *    3x. get the number of local basis functions for master and slave patch
     *   3xi. get the number of the local DOFs for the master and slave patch
     *  3xii. Initialize pointers
     * 3xiii. Loop over all the Gauss Points of the given condition
     * ->
     *        3xiii.1. Get the parametric coordinates of the Gauss Point on the master patch
     *        3xiii.2. Get the parametric coordinates of the Gauss Point on the slave patch
     *        3xiii.3. Find the knot span indices of the Gauss point locations in the parameter space of the master patch
     *        3xiii.4. Find the knot span indices of the Gauss point locations in the parameter space of the slave patch
     *        3xiii.5. Get the tangent to the boundary vector on the master and the slave patch
     *        3xiii.6. Compute elementLength on GP. The weight is already included in variable trCurveGPJacobianProducts
     *        3xiii.7. Compute the B-operator matrices needed for the computation of the patch weak continuity contributions at the master patch
     *        3xiii.8. Compute the B-operator matrices needed for the computation of the patch weak continuity contributions at the slave patch
     *        3xiii.9. Compute the local penalty factor for scaling the rotational contributions
     *       3xiii.10. Scale the B-operator matrices
     *       3xiii.11. Compute the angle of the surface normal vectors
     *       3xiii.12. Determine the alignment of the tangent vectors from both patches at their common interface
     *       3xiii.13. Check if the tangent vectors at the coupled trimming curves are zero and if yes go to the next Gauss point
     *       3xiii.14. Check if the tangent vectors are aligned and if not assert error
     *       3xiii.15. Check if the angle between the tangent vectors is 0° (have the same direction cos(phi) = 1) or 180° (have opposite directions cos(phi) = - 1) to formulate the interface constraint
     *       3xiii.16. Determine the alignment of the normal vectors from both patches at their common interface
     *       3xiii.17. Check if the normal vectors at the coupled trimming curves are zero and if yes go to the next Gauss point
     *       3xiii.18. Check if the boundary normal vectors have an angle of [0,90]U[270,360] meaning the the vectors are in the positive quadrant or [90,270] meaning that the vectors are in the negative quadrant
     *       3xiii.19. Compute the dual product matrices for the displacements
     *       3xiii.20. Compute the dual product matrices for the bending rotations
     *       3xiii.21. Compute the dual product matrices for the twisting rotations
     *       3xiii.22. Compute the element index tables for the master and slave patch
     *       3xiii.23. Compute the element freedom tables for the master and the slave patch
     *       3xiii.24. Loop over all DOFs of the master patch to assemble KPenaltyDisplacementMaster to the global coupling matrix CNN
     *       ->
     *                 3xiii.24i. Assemble the displacement coupling entries
     *                3xiii.24ii. Assemble the bending rotation coupling entries
     *               3xiii.24iii. Assemble the twisting rotation coupling entries
     *       <-
     *       3xiii.25. Loop over all DOFs of the slave patch to assemble KPenaltyDisplacementSlave to the global coupling matrix CNN
     *       ->
     *                 3xiii.25i. Assemble the displacement coupling entries
     *                3xiii.25ii. Assemble the bending rotation coupling entries
     *               3xiii.25iii. Assemble the twisting rotation coupling entries
     *       <-
     *       3xiii.26. Loop over all DOFs of the master/slave patches to assemble CPenaltyDisplacement to the global coupling matrix CNN
     *       ->
     *                 3xiii.26i. Assemble the displacement coupling entries
     *                3xiii.26ii. Assemble the bending rotation coupling entries
     *               3xiii.26iii. Assemble the twisting rotation coupling entries
     *       <-
     *       3xiii.27. Save the Gauss point data for further computation of the error in the L2 norm
     * <-
     *  3xiv. Delete pointers
     * <-
     */

    // 1. Initialize auxiliary variables
    const double tolAngle = 1e-1;
    const double tolVct = 1e-4;
    const int noCoordParam = 2;
    const int noCoord = 3;
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
    double cosPhiTangents;
    double cosPhiNormals;
    double cosPhiSurfaceNormals;
    double condAligned;
    double factorTangent;
    double factorNormal;
    double tangentTrCurveVctMaster[noCoord];
    double tangentTrCurveVctSlave[noCoord];
    double normalTrCurveVctMaster[noCoord];
    double normalTrCurveVctSlave[noCoord];
    double surfaceNormalVctMaster[noCoord];
    double surfaceNormalVctSlave[noCoord];
    double normTangentTrCurveVctMaster;
    double normTangentTrCurveVctSlave;
    double normNormalTrCurveVctMaster;
    double normNormalTrCurveVctSlave;
    double normSurfaceNormalVctMaster;
    double normSurfaceNormalVctSlave;
    double normBOperatorOmegaTMaster;
    double normBOperatorOmegaNMaster;
    double normBOperatorOmegaTSlave;
    double normBOperatorOmegaNSlave;
    double alphaLocal;
    double alphaPrimary;
    double alphaSecondaryBending;
    double alphaSecondaryTwisting;
    double elementLengthOnGP;
    double* trCurveMasterGPs;
    double* trCurveSlaveGPs;
    double* trCurveGPWeights;
    double* trCurveMasterGPTangents;
    double* trCurveSlaveGPTangents;
    double* trCurveGPJacobianProducts;
    IGAPatchSurface* patchMaster;
    IGAPatchSurface* patchSlave;

    // 2. Get the weak patch continuity conditions
    std::vector<WeakIGAPatchContinuityCondition*> weakIGAPatchContinuityConditions = meshIGA->getWeakIGAPatchContinuityConditions();

    // 3. Loop over all the conditions for the application of weak continuity across patch interfaces
    for (int iWCC = 0; iWCC < weakIGAPatchContinuityConditions.size(); iWCC++){
        // 3i. Get the penalty factors for the primary and the secondary field
        alphaPrimary = weakPatchContinuityAlphaPrimaryIJ[iWCC];
        alphaSecondaryBending = weakPatchContinuityAlphaSecondaryBendingIJ[iWCC];
        alphaSecondaryTwisting = weakPatchContinuityAlphaSecondaryTwistingIJ[iWCC];

        // 3ii. Get the index of the master and slave patches
        indexMaster = weakIGAPatchContinuityConditions[iWCC]->getMasterPatchIndex();
        indexSlave = weakIGAPatchContinuityConditions[iWCC]->getSlavePatchIndex();

        // 3iii. Get the number of Gauss Points for the given condition
        noGPsOnContCond = weakIGAPatchContinuityConditions[iWCC]->getTrCurveNumGP();

        // 3iv. Get the parametric coordinates of the Gauss Points
        trCurveMasterGPs = weakIGAPatchContinuityConditions[iWCC]->getTrCurveMasterGPs();
        trCurveSlaveGPs = weakIGAPatchContinuityConditions[iWCC]->getTrCurveSlaveGPs();

        // 3v. Get the corresponding Gauss weights
        trCurveGPWeights = weakIGAPatchContinuityConditions[iWCC]->getTrCurveGPWeights();

        // 3vi. Get the tangent vectors at the trimming curve of the given condition in the Cartesian space
        trCurveMasterGPTangents = weakIGAPatchContinuityConditions[iWCC]->getTrCurveMasterGPTangents();
        trCurveSlaveGPTangents = weakIGAPatchContinuityConditions[iWCC]->getTrCurveSlaveGPTangents();

        // 3vii. Get the product of the Jacobian transformations
        trCurveGPJacobianProducts = weakIGAPatchContinuityConditions[iWCC]->getTrCurveGPJacobianProducts();

        // 3viii. Get the master and the slave patch
        patchMaster = meshIGA->getSurfacePatch(indexMaster);
        patchSlave = meshIGA->getSurfacePatch(indexSlave);

        // 3ix. Get the polynomial orders of the master and the slave patch
        pMaster = patchMaster->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree();
        qMaster = patchMaster->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree();
        pSlave = patchSlave->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree();
        qSlave = patchSlave->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree();

        // 3x. get the number of local basis functions for master and slave patch
        noLocalBasisFctsMaster = (pMaster + 1)*(qMaster + 1);
        noLocalBasisFctsSlave = (pSlave + 1)*(qSlave + 1);

        // 3xi. get the number of the local DOFs for the master and slave patch
        noDOFsLocMaster = noCoord*noLocalBasisFctsMaster;
        noDOFsLocSlave = noCoord*noLocalBasisFctsSlave;

        // 3xii. Initialize pointers
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

        // 3xiii. Loop over all the Gauss Points of the given condition
        for(int iGP = 0; iGP < noGPsOnContCond; iGP++){
            // 3xiii.1. Get the parametric coordinates of the Gauss Point on the master patch
            uGPMaster = trCurveMasterGPs[iGP*noCoordParam];
            vGPMaster = trCurveMasterGPs[iGP*noCoordParam + 1];

            // 3xiii.2. Get the parametric coordinates of the Gauss Point on the slave patch
            uGPSlave = trCurveSlaveGPs[iGP*noCoordParam];
            vGPSlave = trCurveSlaveGPs[iGP*noCoordParam + 1];

            // 3xiii.3. Find the knot span indices of the Gauss point locations in the parameter space of the master patch
            uKnotSpanMaster = patchMaster->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(uGPMaster);
            vKnotSpanMaster = patchMaster->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(vGPMaster);

            // 3xiii.4. Find the knot span indices of the Gauss point locations in the parameter space of the slave patch
            uKnotSpanSlave = patchSlave->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(uGPSlave);
            vKnotSpanSlave = patchSlave->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(vGPSlave);

            // 3xiii.5. Get the tangent to the boundary vector on the master and the slave patch
            for(int iCoord = 0; iCoord < noCoord; iCoord++){
                tangentTrCurveVctMaster[iCoord] = trCurveMasterGPTangents[iGP*noCoord + iCoord];
                tangentTrCurveVctSlave[iCoord] = trCurveSlaveGPTangents[iGP*noCoord + iCoord];
            }

            // 3xiii.6. Compute elementLength on GP. The weight is already included in variable trCurveGPJacobianProducts
            elementLengthOnGP = trCurveGPJacobianProducts[iGP];

            // 3xiii.7. Compute the B-operator matrices needed for the computation of the patch weak continuity contributions at the master patch
            computeDisplacementAndRotationBOperatorMatrices(BDisplacementsGCMaster, BOperatorOmegaTMaster, BOperatorOmegaNMaster, normalTrCurveVctMaster,
                                                            surfaceNormalVctMaster, patchMaster, tangentTrCurveVctMaster, uGPMaster, vGPMaster, uKnotSpanMaster,
                                                            vKnotSpanMaster);

            // 3xiii.8. Compute the B-operator matrices needed for the computation of the patch weak continuity contributions at the slave patch
            computeDisplacementAndRotationBOperatorMatrices(BDisplacementsGCSlave, BOperatorOmegaTSlave, BOperatorOmegaNSlave, normalTrCurveVctSlave,
                                                            surfaceNormalVctSlave, patchSlave, tangentTrCurveVctSlave, uGPSlave, vGPSlave, uKnotSpanSlave,
                                                            vKnotSpanSlave);

            // 3xiii.9. Compute the local penalty factor for scaling the rotational contributions
            if (propWeakPatchContinuityConditions.isSecBendingCoupled) {
                normBOperatorOmegaTMaster = EMPIRE::MathLibrary::vector2norm(BOperatorOmegaTMaster, noDOFsLocMaster);
                alphaLocal = normBOperatorOmegaTMaster;
                normBOperatorOmegaTSlave = EMPIRE::MathLibrary::vector2norm(BOperatorOmegaTSlave, noDOFsLocSlave);
                if (normBOperatorOmegaTSlave > alphaLocal)
                    alphaLocal = normBOperatorOmegaTSlave;
            } else
                alphaLocal = - 1.0;

            if (propWeakPatchContinuityConditions.isSecTwistingCoupled) {
                normBOperatorOmegaNMaster = EMPIRE::MathLibrary::vector2norm(BOperatorOmegaNMaster, noDOFsLocMaster);
                if (normBOperatorOmegaNMaster > alphaLocal)
                    alphaLocal = normBOperatorOmegaNMaster;
                normBOperatorOmegaNSlave = EMPIRE::MathLibrary::vector2norm(BOperatorOmegaNSlave, noDOFsLocSlave);
                if (normBOperatorOmegaNSlave > alphaLocal)
                    alphaLocal = normBOperatorOmegaNSlave;
            }
            alphaLocal = 1.0/abs(alphaLocal);

            // 3xiii.10. Scale the B-operator matrices
            EMPIRE::MathLibrary::computeDenseVectorMultiplicationScalar(BOperatorOmegaTMaster, alphaLocal, noDOFsLocMaster);
            EMPIRE::MathLibrary::computeDenseVectorMultiplicationScalar(BOperatorOmegaNMaster, alphaLocal, noDOFsLocMaster);
            EMPIRE::MathLibrary::computeDenseVectorMultiplicationScalar(BOperatorOmegaTSlave, alphaLocal, noDOFsLocSlave);
            EMPIRE::MathLibrary::computeDenseVectorMultiplicationScalar(BOperatorOmegaNSlave, alphaLocal, noDOFsLocSlave);

            // 3xiii.11. Compute the angle of the surface normal vectors
            normSurfaceNormalVctMaster = EMPIRE::MathLibrary::computeDenseDotProduct(noCoord, surfaceNormalVctMaster, surfaceNormalVctMaster);
            normSurfaceNormalVctMaster = sqrt(normSurfaceNormalVctMaster);
            normSurfaceNormalVctSlave = EMPIRE::MathLibrary::computeDenseDotProduct(noCoord, surfaceNormalVctSlave, surfaceNormalVctSlave);
            normSurfaceNormalVctSlave = sqrt(normSurfaceNormalVctSlave);
            cosPhiSurfaceNormals = EMPIRE::MathLibrary::computeDenseDotProduct(noCoord, surfaceNormalVctMaster, surfaceNormalVctSlave);
            cosPhiSurfaceNormals = cosPhiSurfaceNormals/(normSurfaceNormalVctMaster*normSurfaceNormalVctSlave);

            // 3xiii.12. Determine the alignment of the tangent vectors from both patches at their common interface
            normTangentTrCurveVctMaster = EMPIRE::MathLibrary::computeDenseDotProduct(noCoord,tangentTrCurveVctMaster,tangentTrCurveVctMaster);
            normTangentTrCurveVctMaster = sqrt(normTangentTrCurveVctMaster);
            normTangentTrCurveVctSlave = EMPIRE::MathLibrary::computeDenseDotProduct(noCoord,tangentTrCurveVctSlave,tangentTrCurveVctSlave);
            normTangentTrCurveVctSlave = sqrt(normTangentTrCurveVctSlave);
            cosPhiTangents = EMPIRE::MathLibrary::computeDenseDotProduct(noCoord,tangentTrCurveVctMaster,tangentTrCurveVctSlave);
            cosPhiTangents = cosPhiTangents/(normTangentTrCurveVctMaster*normTangentTrCurveVctSlave);

            // 3xiii.13. Check if the tangent vectors at the coupled trimming curves are zero and if yes go to the next Gauss point
            if(normTangentTrCurveVctMaster < tolVct && normTangentTrCurveVctSlave < tolVct)
                continue;
            else if((normTangentTrCurveVctMaster < tolVct && normTangentTrCurveVctSlave > tolVct) || (normTangentTrCurveVctMaster > tolVct && normTangentTrCurveVctSlave < tolVct))
                assert(false);

            // 3xiii.14. Check if the tangent vectors are aligned and if not assert error
            condAligned = cosPhiTangents*cosPhiTangents - 1;
            if (abs(condAligned) > tolAngle) {
                INFO_OUT() << "Found boundaries for which the tangent vectors are not aligned with angle = " << cosPhiTangents << endl;
            }
            assert(abs(condAligned) < tolAngle);

            // 3xiii.15. Check if the angle between the tangent vectors is 0° (have the same direction cos(phi) = 1) or 180° (have opposite directions cos(phi) = - 1) to formulate the interface constraint
            if (abs(cosPhiTangents - 1) < sqrt(tolAngle)) {
                factorTangent = - 1.0;
            } else if (abs(cosPhiTangents + 1) < sqrt(tolAngle)) {
                factorTangent = + 1.0;
            } else {
                assert(false);
            }

            // 3xiii.16. Determine the alignment of the normal vectors from both patches at their common interface
            normNormalTrCurveVctMaster = EMPIRE::MathLibrary::computeDenseDotProduct(noCoord,normalTrCurveVctMaster,normalTrCurveVctMaster);
            normNormalTrCurveVctMaster = sqrt(normNormalTrCurveVctMaster);
            normNormalTrCurveVctSlave = EMPIRE::MathLibrary::computeDenseDotProduct(noCoord,normalTrCurveVctSlave,normalTrCurveVctSlave);
            normNormalTrCurveVctSlave = sqrt(normNormalTrCurveVctSlave);
            cosPhiNormals = EMPIRE::MathLibrary::computeDenseDotProduct(noCoord,normalTrCurveVctMaster,normalTrCurveVctSlave);
            cosPhiNormals = cosPhiTangents/(normNormalTrCurveVctMaster*normNormalTrCurveVctSlave);

            // 3xiii.17. Check if the normal vectors at the coupled trimming curves are zero and if yes go to the next Gauss point
            if(normNormalTrCurveVctMaster < tolVct && normNormalTrCurveVctSlave < tolVct)
                continue;
            else if((normNormalTrCurveVctMaster < tolVct && normNormalTrCurveVctSlave > tolVct) || (normNormalTrCurveVctMaster > tolVct && normNormalTrCurveVctSlave < tolVct))
                assert(false);

            // 3xiii.18. Check if the boundary normal vectors have an angle of [0,90]U[270,360] meaning the the vectors are in the positive quadrant or [90,270] meaning that the vectors are in the negative quadrant
            if (cosPhiNormals >= 0)
                factorNormal = + 1.0;
            else
                factorNormal = - 1.0;
            factorNormal *= - 1.0;

            // Check whether the surface normal vectors have an angle of [0,90]U[270,360] meaning the the patches have the same normal orientation or [90,270] meaning that the patches have opossite normal orientation
//            if (cosPhiSurfaceNormals < 0.0) {
//                factorTangent = (-1.0)*factorTangent;
//                factorNormal = (-1.0)*factorNormal;
//            }

            // 3xiii.19. Compute the dual product matrices for the displacements
            EMPIRE::MathLibrary::computeTransposeMatrixProduct(noCoord,noDOFsLocMaster,noDOFsLocMaster,BDisplacementsGCMaster,BDisplacementsGCMaster,KPenaltyDisplacementMaster);
            EMPIRE::MathLibrary::computeTransposeMatrixProduct(noCoord,noDOFsLocSlave,noDOFsLocSlave,BDisplacementsGCSlave,BDisplacementsGCSlave,KPenaltyDisplacementSlave);
            EMPIRE::MathLibrary::computeTransposeMatrixProduct(noCoord,noDOFsLocMaster,noDOFsLocSlave,BDisplacementsGCMaster,BDisplacementsGCSlave,CPenaltyDisplacement);

            // 3xiii.20. Compute the dual product matrices for the bending rotations
            EMPIRE::MathLibrary::computeTransposeMatrixProduct(1, noDOFsLocMaster, noDOFsLocMaster ,BOperatorOmegaTMaster, BOperatorOmegaTMaster, KPenaltyBendingRotationMaster);
            EMPIRE::MathLibrary::computeTransposeMatrixProduct(1, noDOFsLocSlave, noDOFsLocSlave, BOperatorOmegaTSlave, BOperatorOmegaTSlave, KPenaltyBendingRotationSlave);
            EMPIRE::MathLibrary::computeTransposeMatrixProduct(1, noDOFsLocMaster, noDOFsLocSlave, BOperatorOmegaTMaster, BOperatorOmegaTSlave,CPenaltyBendingRotation);

            // 3xiii.21. Compute the dual product matrices for the twisting rotations
            EMPIRE::MathLibrary::computeTransposeMatrixProduct(1, noDOFsLocMaster, noDOFsLocMaster ,BOperatorOmegaNMaster, BOperatorOmegaNMaster, KPenaltyTwistingRotationMaster);
            EMPIRE::MathLibrary::computeTransposeMatrixProduct(1, noDOFsLocSlave, noDOFsLocSlave, BOperatorOmegaNSlave, BOperatorOmegaNSlave, KPenaltyTwistingRotationSlave);
            EMPIRE::MathLibrary::computeTransposeMatrixProduct(1, noDOFsLocMaster, noDOFsLocSlave, BOperatorOmegaNMaster, BOperatorOmegaNSlave,CPenaltyTwistingRotation);

            // 3xiii.22. Compute the element index tables for the master and slave patch
            int CPIndexMaster[noLocalBasisFctsMaster];
            int CPIndexSlave[noLocalBasisFctsSlave];
            patchMaster->getIGABasis()->getBasisFunctionsIndex(uKnotSpanMaster, vKnotSpanMaster, CPIndexMaster);
            patchSlave->getIGABasis()->getBasisFunctionsIndex(uKnotSpanSlave, vKnotSpanSlave, CPIndexSlave);

            // 3xiii.23. Compute the element freedom tables for the master and the slave patch
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

            // 3xiii.24. Loop over all DOFs of the master patch to assemble KPenaltyDisplacementMaster to the global coupling matrix CNN
            for(int i = 0; i < noDOFsLocMaster; i++){
                for(int j = 0; j < noDOFsLocMaster; j++){
                    // 3xiii.24i. Assemble the displacement coupling entries
                    if (propWeakPatchContinuityConditions.isPrimCoupled)
                        couplingMatrices->addCNNValue(EFTMaster[i], EFTMaster[j], alphaPrimary*KPenaltyDisplacementMaster[i*noDOFsLocMaster + j]*elementLengthOnGP);

                    // 3xiii.24ii. Assemble the bending rotation coupling entries
                    if (propWeakPatchContinuityConditions.isSecBendingCoupled)
                        couplingMatrices->addCNNValue(EFTMaster[i], EFTMaster[j], alphaSecondaryBending*KPenaltyBendingRotationMaster[i*noDOFsLocMaster + j]*elementLengthOnGP);

                    // 3xiii.24iii. Assemble the twisting rotation coupling entries
                    if (propWeakPatchContinuityConditions.isSecTwistingCoupled)
                        couplingMatrices->addCNNValue(EFTMaster[i], EFTMaster[j], alphaSecondaryTwisting*KPenaltyTwistingRotationMaster[i*noDOFsLocMaster + j]*elementLengthOnGP);
                }
            }

            // 3xiii.25. Loop over all DOFs of the slave patch to assemble KPenaltyDisplacementSlave to the global coupling matrix CNN
            for(int i = 0; i < noDOFsLocSlave; i++){
                for(int j = 0; j < noDOFsLocSlave; j++) {
                    // 3xiii.25i. Assemble the displacement coupling entries
                    if (propWeakPatchContinuityConditions.isPrimCoupled)
                        couplingMatrices->addCNNValue(EFTSlave[i], EFTSlave[j], alphaPrimary*KPenaltyDisplacementSlave[i*noDOFsLocSlave + j]*elementLengthOnGP);

                    // 3xiii.25ii. Assemble the bending rotation coupling entries
                    if (propWeakPatchContinuityConditions.isSecBendingCoupled)
                        couplingMatrices->addCNNValue(EFTSlave[i], EFTSlave[j], alphaSecondaryBending*KPenaltyBendingRotationSlave[i*noDOFsLocSlave + j]*elementLengthOnGP);

                    // 3xiii.25iii. Assemble the twisting rotation coupling entries
                    if (propWeakPatchContinuityConditions.isSecTwistingCoupled)
                        couplingMatrices->addCNNValue(EFTSlave[i], EFTSlave[j], alphaSecondaryTwisting*KPenaltyTwistingRotationSlave[i*noDOFsLocSlave + j]*elementLengthOnGP);
                }
            }

            // 3xiii.26. Loop over all DOFs of the master/slave patches to assemble CPenaltyDisplacement to the global coupling matrix CNN
            for(int i = 0; i < noDOFsLocMaster; i++){
                for(int j = 0; j < noDOFsLocSlave; j++){
                    // 3xiii.26i. Assemble the displacement coupling entries
                    if (propWeakPatchContinuityConditions.isPrimCoupled) {
                        couplingMatrices->addCNNValue(EFTMaster[i], EFTSlave[j], alphaPrimary*(-1.0)*CPenaltyDisplacement[i*noDOFsLocSlave + j]*elementLengthOnGP);
                        couplingMatrices->addCNNValue(EFTSlave[j], EFTMaster[i], alphaPrimary*(-1.0)*CPenaltyDisplacement[i*noDOFsLocSlave + j]*elementLengthOnGP);
                    }

                    // 3xiii.26ii. Assemble the bending rotation coupling entries
                    if (propWeakPatchContinuityConditions.isSecBendingCoupled) {
                        couplingMatrices->addCNNValue(EFTMaster[i], EFTSlave[j], alphaSecondaryBending*factorTangent*CPenaltyBendingRotation[i*noDOFsLocSlave + j]*elementLengthOnGP);
                        couplingMatrices->addCNNValue(EFTSlave[j], EFTMaster[i], alphaSecondaryBending*factorTangent*CPenaltyBendingRotation[i*noDOFsLocSlave + j]*elementLengthOnGP);
                    }

                    // 3xiii.26iii. Assemble the twisting rotation coupling entries
                    if (propWeakPatchContinuityConditions.isSecTwistingCoupled) {
                        couplingMatrices->addCNNValue(EFTMaster[i], EFTSlave[j], alphaSecondaryTwisting*factorNormal*CPenaltyTwistingRotation[i*noDOFsLocSlave + j]*elementLengthOnGP);
                        couplingMatrices->addCNNValue(EFTSlave[j], EFTMaster[i], alphaSecondaryTwisting*factorNormal*CPenaltyTwistingRotation[i*noDOFsLocSlave + j]*elementLengthOnGP);
                    }
                }
            }

            // 3xiii.27. Save the Gauss point data for further computation of the error in the L2 norm
            if(propErrorComputation.isInterfaceError){
                // Initialize variable storing the Gauss Point data
                std::vector<double> streamInterfaceGP;

                // elementLengthOnGP + noBasisFuncsI + (#indexCP, basisFuncValueI,...) + (#indexDOF, BtValueI, BnValueI,...) + noBasisFuncsJ + (#indexCP, basisFuncValueJ,...) + (#indexDOF, BtValueJ, BnValueJ,...) + factorTangent + factorNormal
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

                // Push back the Gauss Point values into the member variable
                streamInterfaceGPs.push_back(streamInterfaceGP);
            }
        } // End of Gauss Point loop

        // 3xiv. Delete pointers
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
                                                                      double* _surfaceNormalVct, IGAPatchSurface* _patch,
                                                                      double* _tangentTrCurveVct, double _u, double _v,
                                                                      int _uKnotSpan, int _vKnotSpan) {
    /*
     * Returns the B-operator matrix for the displacement field, the bending and the twisting rotations in the global Cartesian space.
     * Additionally the normal to the boundary vector which is tangent to the surface patch and the surace normal vector are returned.
     *
     * Hints on the size of the [in/out] arrays:
     *
     * double* _BDisplacementsGC = new double[noCoord*noDOFsLoc]
     * double* _BOperatorOmegaT = new double[noDOFsLoc]
     * double* _BOperatorOmegaN = new double[noDOFsLoc]
     * double _normalTrCurveVct[noCoord]
     *
     * Function layout:
     *
     * 1. Initialize auxiliary arrays
     *
     * 2. Compute the basis functions and their derivatives
     *
     * 3. Compute the base vectors and their derivatives
     *
     * 4. Compute the derivatives of the surface normal vector
     *
     * 5. Compute the normal to the boundary vector which is tangent to the surface
     *
     * 6. Compute the covariant metric tensor
     *
     * 7. Compute the contravariant base vectors of the master patch
     *
     * 8. Initialize the B-operator matrix for the displacement field
     *
     * 9. Compute the B-operator matrix for the displacement field and its parametric derivatives
     *
     * 10. Transform the normal and the tangent vectors to the covariant base for the master patch
     *
     * 11. Compute the curvature tensor in the contravariant basis
     *
     * 12. Compute the B-operator matrices for the rotation vector
     *
     * 13. Delete pointers
     */

    // 1. Initialize auxiliary arrays
    const int noCoord = 3;
    const int noParametricCoord = 2;
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
    int p = _patch->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree();
    int q = _patch->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree();
    int noLocalBasisFcts = (p + 1)*(q + 1);
    int noDOFsLoc = noCoord*noLocalBasisFcts;
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

    // 2. Compute the basis functions and their derivatives
    _patch->getIGABasis()->computeLocalBasisFunctionsAndDerivatives(basisFctsAndDerivs, derivDegreeBasis, _u, _uKnotSpan,
                                                                    _v, _vKnotSpan);

    // 3. Compute the base vectors and their derivatives
    _patch->computeBaseVectorsAndDerivatives(baseVctsAndDerivs, basisFctsAndDerivs, derivDegreeBaseVec,
                                             _uKnotSpan, _vKnotSpan);

    // 4. Compute the derivatives of the surface normal vector
    _patch->computeSurfaceNormalVectorAndDerivatives(surfNormalVctAndDervs, baseVctsAndDerivs, derivDegreeBaseVec);

    // 5. Compute the normal to the boundary vector which is tangent to the surface
    for(int i = 0; i < noCoord; i++){
        _surfaceNormalVct[i] = surfNormalVctAndDervs[i];
    }
    EMPIRE::MathLibrary::computeVectorCrossProduct(_surfaceNormalVct, _tangentTrCurveVct, _normalTrCurveVct);

    // 6. Compute the covariant metric tensor
    _patch->computeCovariantMetricTensor(covariantMetricTensor, baseVctsAndDerivs, derivDegreeBaseVec);

    // 7. Compute the contravariant base vectors of the master patch
    _patch->computeContravariantBaseVectors(contravariantBaseVcts, covariantMetricTensor, baseVctsAndDerivs, derivDegreeBaseVec);

    // 8. Initialize the B-operator matrix for the displacement field
    for(int iBF = 0; iBF < noCoord*noDOFsLoc; iBF++){
        _BDisplacementsGC[iBF] = 0.0;
        BdDisplacementsdUGC[iBF] = 0.0;
        BdDisplacementsdVGC[iBF] = 0.0;
    }

    // 9. Compute the B-operator matrix for the displacement field and its parametric derivatives
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

    // 10. Transform the normal and the tangent vectors to the covariant base for the master patch
    for(int iCov = 0; iCov < noParametricCoord; iCov++){
        for (int iCoord = 0; iCoord < noCoord; iCoord++)
            contravariantBaseVct[iCoord] = contravariantBaseVcts[noCoord*iCov + iCoord];
        tangentTrCurveVctCov[iCov] = EMPIRE::MathLibrary::computeDenseDotProduct(noCoord,contravariantBaseVct,_tangentTrCurveVct);
        normalTrCurveVctCov[iCov] = EMPIRE::MathLibrary::computeDenseDotProduct(noCoord,contravariantBaseVct,_normalTrCurveVct);
    }

    // 11. Compute the curvature tensor in the contravariant basis
    _patch->computeContravariantCurvatureTensor(contravariantCurvatureTensor, _surfaceNormalVct, baseVctsAndDerivs, derivDegreeBaseVec);

    // 12. Compute the B-operator matrices for the rotation vector
    for(int iCovCoord = 0; iCovCoord < noParametricCoord; iCovCoord++)
        for(int iCoord = 0; iCoord < noCoord; iCoord++)
            dT3Cov2GC[iCovCoord*noCoord + iCoord] = surfNormalVctAndDervs[noCoord*(iCovCoord + 1) + iCoord]; // Parametric derivatives of the transformation matrix from the covariant to the global Cartesian basis
    EMPIRE::MathLibrary::computeMatrixProduct(noParametricCoord, noCoord, noDOFsLoc, dT3Cov2GC, _BDisplacementsGC, commonBOperator1);
    EMPIRE::MathLibrary::computeMatrixProduct(1, noCoord, noDOFsLoc, _surfaceNormalVct, BdDisplacementsdUGC, commonBOperator2Part1);
    EMPIRE::MathLibrary::computeMatrixProduct(1, noCoord, noDOFsLoc, _surfaceNormalVct, BdDisplacementsdVGC, commonBOperator2Part2);
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

    // 13. Delete pointers
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
     *
     * Function layout:
     *
     * 1. Initialize a file stream to write out the Penalty parameters
     *
     * 2. Initialize auxiliary arrays
     *
     * 3. Get the weak patch continuity conditions
     *
     * 4. Loop over all the conditions for the application of weak continuity across patch interfaces
     * ->
     *    4i. Check if penalty factors are to be assigned manually in the xml file
     *   4ii. Get the index of the patch
     *  4iii. Get the number of Gauss Points for the given condition
     *   4iv. Get the parametric coordinates of the Gauss Points
     *    4v. Get the corresponding Gauss weights
     *   4vi. Get the tangent vectors at the curve of the given condition in the Cartesian space
     *  4vii. Get the product of the Jacobian transformations
     * 4viii. Get the master and the slave patch
     *   4ix. Get the polynomial orders of the patch
     *    4x. Get the number of local basis functions for the patch
     *   4xi. Get the number of the local DOFs for the patch
     *  4xii. Initialize flags on whether an element has been changed while looping over the Gauss Points
     * 4xiii. Initialize the element edge sizes
     *  4xiv. Loop over all the Gauss Points of the given condition
     *  ->
     *        4xiv.1. Get the parametric coordinates of the Gauss Point on the patch
     *        4xiv.2. Find the knot span indices of the Gauss point locations in the parameter space of the patch
     *        4xiv.3. Initialize the saved knot span indices
     *        4xiv.4. Initialize element edge sizes if an element has been crossed
     *        4xiv.5. Add the contribution from the Gauss Point to the element edge sizes
     *        4xiv.6. Save the knot span indices
     *  <-
     *   4xv. Check the element sizes for the last elements
     *  4xvi. Compute correspondingly the penalty factors
     * 4xvii. Write the Penalty parameters into a file
     * <-
     *
     * 5. Close file
     */

    // 1. Initialize a file stream to write out the Penalty parameters
    std::ofstream ofs;
    if (!_filename.empty()) {
        ofs.open(_filename.c_str(), std::ofstream::out);
        ofs << std::scientific;
        ofs << std::endl;
    }

    // 2. Initialize auxiliary arrays
    const int noCoord = 3;
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
    double* curveGPs;
    double* curveGPWeights;
    double* curveGPTangents;
    double* curveGPJacobianProducts;
    IGAPatchSurface* thePatch;

    // 3. Get the weak patch continuity conditions
    std::vector<WeakIGADirichletCurveCondition*> weakIGADirichletCurveConditions = meshIGA->getWeakIGADirichletCurveConditions();

    // 4. Loop over all the conditions for the application of weak continuity across patch interfaces
    for (int iWDCC = 0; iWDCC < weakIGADirichletCurveConditions.size(); iWDCC++){
        // 4i. Check if penalty factors are to be assigned manually in the xml file
        if(!propWeakCurveDirichletConditions.isAutomaticPenaltyParameters){
            weakDirichletCCAlphaPrimary[iWDCC] = propWeakCurveDirichletConditions.alphaPrim;
            weakDirichletCCAlphaSecondaryBending[iWDCC] = propWeakCurveDirichletConditions.alphaSecBending;
            weakDirichletCCAlphaSecondaryTwisting[iWDCC] = propWeakCurveDirichletConditions.alphaSecTwisting;
            continue;
        }

        // 4ii. Get the index of the patch
        patchIndex = weakIGADirichletCurveConditions[iWDCC]->getPatchIndex();

        // 4iii. Get the number of Gauss Points for the given condition
        noGPsOnCond = weakIGADirichletCurveConditions[iWDCC]->getCurveNumGP();

        // 4iv. Get the parametric coordinates of the Gauss Points
        curveGPs = weakIGADirichletCurveConditions[iWDCC]->getCurveGPs();

        // 4v. Get the corresponding Gauss weights
        curveGPWeights = weakIGADirichletCurveConditions[iWDCC]->getCurveGPWeights();

        // 4vi. Get the tangent vectors at the curve of the given condition in the Cartesian space
        curveGPTangents = weakIGADirichletCurveConditions[iWDCC]->getCurveGPTangents();

        // 4vii. Get the product of the Jacobian transformations
        curveGPJacobianProducts = weakIGADirichletCurveConditions[iWDCC]->getCurveGPJacobianProducts();

        // 4viii. Get the master and the slave patch
        thePatch = meshIGA->getSurfacePatch(patchIndex);

        // 4ix. Get the polynomial orders of the patch
        p = thePatch->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree();
        q = thePatch->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree();
        pMax = std::max(p, q);

        // 4x. Get the number of local basis functions for the patch
        noLocalBasisFcts = (p + 1)*(q + 1);

        // 4xi. Get the number of the local DOFs for the patch
        noDOFsLoc = noCoord*noLocalBasisFcts;

        // 4xii. Initialize flags on whether an element has been changed while looping over the Gauss Points
        isElementChanged = false;

        // 4xiii. Initialize the element edge sizes
        elEdgeSize = 0.0;
        minElEdgeSize = std::numeric_limits<double>::max();

        // 4xiv. Loop over all the Gauss Points of the given condition
        for(int iGP = 0; iGP < noGPsOnCond; iGP++){
            // 4xiv.1. Get the parametric coordinates of the Gauss Point on the patch
            uGP = curveGPs[2*iGP];
            vGP = curveGPs[2*iGP + 1];

            // 4xiv.2. Find the knot span indices of the Gauss point locations in the parameter space of the patch
            uKnotSpan = thePatch->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(uGP);
            vKnotSpan = thePatch->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(vGP);

            // 4xiv.3. Initialize the saved knot span indices
            if(iGP == 0){
                uKnotSpanSaved = uKnotSpan;
                vKnotSpanSaved = vKnotSpan;
            }

            // 4xiv.4. Initialize element edge sizes if an element has been crossed
            if(uKnotSpan != uKnotSpanSaved || vKnotSpan != vKnotSpanSaved){
                if(elEdgeSize < minElEdgeSize)
                    minElEdgeSize = elEdgeSize;
                elEdgeSize = 0.0;
            }

            // 4xiv.5. Add the contribution from the Gauss Point to the element edge sizes
            elEdgeSize += curveGPJacobianProducts[iGP];

            // 4xiv.6. Save the knot span indices
            uKnotSpanSaved = uKnotSpan;
            vKnotSpanSaved = vKnotSpan;

        } // End of Gauss Point loop

        // 4xv. Check the element sizes for the last elements
        if(elEdgeSize < minElEdgeSize)
            minElEdgeSize = elEdgeSize;

        // 4xvi. Compute correspondingly the penalty factors
        alphaPrim = pMax/minElEdgeSize;
        alphaSec = pMax/sqrt(minElEdgeSize);
        weakDirichletCCAlphaPrimary[iWDCC] = alphaPrim;
        weakDirichletCCAlphaSecondaryBending[iWDCC] = alphaPrim;

        // 4xvii. Write the Penalty parameters into a file
        if (!_filename.empty()) {
            ofs << "Weak Dirichlet curve condition on patch[" << patchIndex << "] :" << std::endl;
            ofs << "weakDirichletCCAlpha[" << iWDCC << "] = " << scientific << setprecision(15) << alphaPrim << std::endl;
            ofs << std::endl;
        }
    } // End of weak Dirichlet curve condition loop

    // 5. Close file
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
            DEBUG_OUT() << std::endl;
        }
    }
}

void IGAMortarMapper::computePenaltyParametersForPatchContinuityConditions(std::string _filename){
    /*
     * Compute the penalty factors related to the application of weak patch continuity conditions
     * as a function of the minimum element edge size across each interface.
     *
     * 1. Initialize a file stream to write out the Penalty parameters
     *
     * 2. Initialize auxiliary variables
     *
     * 3. Get the weak patch continuity conditions
     *
     * 4. Loop over all the conditions for the application of weak continuity across patch interfaces
     * ->
     *    4i. Check if penalty factors are to be assigned manually in the xml file
     *   4ii. Get the index of the master and slave patches
     *  4iii. Get the number of Gauss Points for the given condition
     *   4iv. Get the parametric coordinates of the Gauss Points
     *    4v. Get the corresponding Gauss weights
     *   4vi. Get the tangent vectors at the trimming curve of the given condition in the Cartesian space
     *  4vii. Get the product of the Jacobian transformations
     * 4viii. Get the master and the slave patch
     *   4ix. Get the polynomial orders of the master and the slave patch
     *    4x. Get the number of local basis functions for master and slave patch
     *   4xi. Get the number of the local DOFs for the master and slave patch
     *  4xii. Initialize flags on whether an element has been changed while looping over the Gauss Points
     * 4xiii. Initialize the element edge sizes
     *  4xiv. Loop over all the Gauss Points of the given condition
     *  ->
     *        4xiv.1. Get the parametric coordinates of the Gauss Point on the master patch
     *        4xiv.2. Get the parametric coordinates of the Gauss Point on the slave patch
     *        4xiv.3. Find the knot span indices of the Gauss point locations in the parameter space of the master patch
     *        4xiv.4. Find the knot span indices of the Gauss point locations in the parameter space of the slave patch
     *        4xiv.5. Initialize the saved knot span indices
     *        4xiv.6. Initialize element edge sizes if an element has been crossed
     *        4xiv.7. Add the contribution from the Gauss Point to the element edge sizes
     *        4xiv.8. Save the knot span indices
     *  <-
     *   4xv. Check the element sizes for the last elements
     *  4xvi. Get the minimum of the minimum element edge sizes between both patches
     * 4xvii. Compute correspondingly the penalty factors
     *4xviii. Write the Penalty parameters into a file
     * <-
     *
     * 5. Close file
     */

    // 1. Initialize a file stream to write out the Penalty parameters
    std::ofstream ofs;
    if (!_filename.empty()) {
        ofs.open(_filename.c_str(), std::ofstream::out);
        ofs << std::scientific;
        ofs << std::endl;
    }

    // 2. Initialize auxiliary variables
    const int noCoord = 3;
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
    double alphaBar;
    double* trCurveMasterGPs;
    double* trCurveSlaveGPs;
    double* trCurveGPWeights;
    double* trCurveMasterGPTangents;
    double* trCurveSlaveGPTangents;
    double* trCurveGPJacobianProducts;
    IGAPatchSurface* patchMaster;
    IGAPatchSurface* patchSlave;

    // 3. Get the weak patch continuity conditions
    std::vector<WeakIGAPatchContinuityCondition*> weakIGAPatchContinuityConditions = meshIGA->getWeakIGAPatchContinuityConditions();

    // 4. Loop over all the conditions for the application of weak continuity across patch interfaces
    for (int iWCC = 0; iWCC < weakIGAPatchContinuityConditions.size(); iWCC++){
        // 4i. Check if penalty factors are to be assigned manually in the xml file
        if(!propWeakPatchContinuityConditions.isAutomaticPenaltyParameters){
            weakPatchContinuityAlphaPrimaryIJ[iWCC] = propWeakPatchContinuityConditions.alphaPrim;
            weakPatchContinuityAlphaSecondaryBendingIJ[iWCC] = propWeakPatchContinuityConditions.alphaSecBending;
            weakPatchContinuityAlphaSecondaryTwistingIJ[iWCC] = propWeakPatchContinuityConditions.alphaSecTwisting;
            continue;
        }

        // 4ii. Get the index of the master and slave patches
        indexMaster = weakIGAPatchContinuityConditions[iWCC]->getMasterPatchIndex();
        indexSlave = weakIGAPatchContinuityConditions[iWCC]->getSlavePatchIndex();

        // 4iii. Get the number of Gauss Points for the given condition
        noGPsOnContCond = weakIGAPatchContinuityConditions[iWCC]->getTrCurveNumGP();

        // 4iv. Get the parametric coordinates of the Gauss Points
        trCurveMasterGPs = weakIGAPatchContinuityConditions[iWCC]->getTrCurveMasterGPs();
        trCurveSlaveGPs = weakIGAPatchContinuityConditions[iWCC]->getTrCurveSlaveGPs();

        // 4v. Get the corresponding Gauss weights
        trCurveGPWeights = weakIGAPatchContinuityConditions[iWCC]->getTrCurveGPWeights();

        // 4vi. Get the tangent vectors at the trimming curve of the given condition in the Cartesian space
        trCurveMasterGPTangents = weakIGAPatchContinuityConditions[iWCC]->getTrCurveMasterGPTangents();
        trCurveSlaveGPTangents = weakIGAPatchContinuityConditions[iWCC]->getTrCurveSlaveGPTangents();

        // 4vii. Get the product of the Jacobian transformations
        trCurveGPJacobianProducts = weakIGAPatchContinuityConditions[iWCC]->getTrCurveGPJacobianProducts();

        // 4viii. Get the master and the slave patch
        patchMaster = meshIGA->getSurfacePatch(indexMaster);
        patchSlave = meshIGA->getSurfacePatch(indexSlave);

        // 4ix. Get the polynomial orders of the master and the slave patch
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

        // 4x. Get the number of local basis functions for master and slave patch
        noLocalBasisFctsMaster = (pMaster + 1)*(qMaster + 1);
        noLocalBasisFctsSlave = (pSlave + 1)*(qSlave + 1);

        // 4xi. Get the number of the local DOFs for the master and slave patch
        noDOFsLocMaster = noCoord*noLocalBasisFctsMaster;
        noDOFsLocSlave = noCoord*noLocalBasisFctsSlave;

        // 4xii. Initialize flags on whether an element has been changed while looping over the Gauss Points
        isElementMasterChanged = false;
        isElementSlaveChanged = false;

        // 4xiii. Initialize the element edge sizes
        elEdgeSizeMaster = 0.0;
        elEdgeSizeSlave = 0.0;
        minElEdgeSizeMaster = std::numeric_limits<double>::max();
        minElEdgeSizeSlave = std::numeric_limits<double>::max();

        // 4xiv. Loop over all the Gauss Points of the given condition
        for(int iGP = 0; iGP < noGPsOnContCond; iGP++){
            // 4xiv.1. Get the parametric coordinates of the Gauss Point on the master patch
            uGPMaster = trCurveMasterGPs[2*iGP];
            vGPMaster = trCurveMasterGPs[2*iGP + 1];

            // 4xiv.2. Get the parametric coordinates of the Gauss Point on the slave patch
            uGPSlave = trCurveSlaveGPs[2*iGP];
            vGPSlave = trCurveSlaveGPs[2*iGP + 1];

            // 4xiv.3. Find the knot span indices of the Gauss point locations in the parameter space of the master patch
            uKnotSpanMaster = patchMaster->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(uGPMaster);
            vKnotSpanMaster = patchMaster->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(vGPMaster);

            // 4xiv.4. Find the knot span indices of the Gauss point locations in the parameter space of the slave patch
            uKnotSpanSlave = patchSlave->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(uGPSlave);
            vKnotSpanSlave = patchSlave->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(vGPSlave);

            // 4xiv.5. Initialize the saved knot span indices
            if(iGP == 0){
                uKnotSpanMasterSaved = uKnotSpanMaster;
                vKnotSpanMasterSaved = vKnotSpanMaster;
                uKnotSpanSlaveSaved = uKnotSpanSlave;
                vKnotSpanSlaveSaved = vKnotSpanSlave;
            }

            // 4xiv.6. Initialize element edge sizes if an element has been crossed
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

            // 4xiv.7. Add the contribution from the Gauss Point to the element edge sizes
            elEdgeSizeMaster += trCurveGPJacobianProducts[iGP];
            elEdgeSizeSlave += trCurveGPJacobianProducts[iGP];

            // 4xiv.8. Save the knot span indices
            uKnotSpanMasterSaved = uKnotSpanMaster;
            vKnotSpanMasterSaved = vKnotSpanMaster;
            uKnotSpanSlaveSaved = uKnotSpanSlave;
            vKnotSpanSlaveSaved = vKnotSpanSlave;
        } // End of Gauss Point loop

        // 4xv. Check the element sizes for the last elements
        if(elEdgeSizeMaster < minElEdgeSizeMaster)
            minElEdgeSizeMaster = elEdgeSizeMaster;
        if(elEdgeSizeSlave < minElEdgeSizeSlave)
            minElEdgeSizeSlave = elEdgeSizeSlave;

        // 4xvi. Get the minimum of the minimum element edge sizes between both patches
        minElEdgeSize = minElEdgeSizeMaster;
        if(minElEdgeSizeSlave < minElEdgeSize){
            minElEdgeSize = minElEdgeSizeSlave;
        }

        // 4xvii. Compute correspondingly the penalty factors
        alphaBar = pMax/minElEdgeSize;
        weakPatchContinuityAlphaPrimaryIJ[iWCC] = alphaBar;
        weakPatchContinuityAlphaSecondaryBendingIJ[iWCC] = alphaBar;
        weakPatchContinuityAlphaSecondaryTwistingIJ[iWCC] = alphaBar;

        // 4xviii. Write the Penalty parameters into a file
        if (!_filename.empty()) {
            ofs << "Coupling between patch[" << indexMaster << "] and patch[" << indexSlave << "]:" << std::endl;
            ofs << "weakPatchContinuityAlphaPrimaryIJ[" << iWCC << "] = " << alphaBar << std::endl;
            ofs << std::endl;
        }
    } // End of weak continuity condition loop

    // 5. Close file
    if (!_filename.empty())
        ofs.close();
}

void IGAMortarMapper::consistentMapping(const double* _slaveField, double *_masterField) {
    /*
     * Computes the master field given the slave field using the isogeometric mortar based method,
     * Cnn * x_master = Cnr * x_slave
     *
     * Function layout:
     *
     * 1. Initialize auxiliary variables
     *
     * 2. Compute the right hand side of the isogeometric mortar-based mapping method, C_NR * x_slave = tmpVec
     *
     * 3. Solve for the master field using the isogeometric mortar-based mapping method, Cnn * x_master = tmpVec
     *
     * 4. Delete pointers
     */

    // 1. Initialize auxiliary variables
    int size_N = couplingMatrices->getSizeN();
    double* tmpVec = new double[size_N]();

    // 2. Compute the right hand side of the isogeometric mortar-based mapping method, C_NR * x_slave = tmpVec
    couplingMatrices->getCnr()->mulitplyVec(false,const_cast<double *>(_slaveField), tmpVec, size_N);

    // 3. Solve for the master field using the isogeometric mortar-based mapping method, Cnn * x_master = tmpVec
    couplingMatrices->getCnn()->solve(_masterField, tmpVec);

    // 4. Delete pointers
    delete[] tmpVec;
}

void IGAMortarMapper::conservativeMapping(const double* _masterField, double *_slaveField) {
    /*
     * Computes the slave field given the master field using the transpose transformation matrix
     * from the isogeometric mortar-based method as,
     *
     * f_slave = (Cnn^(-1) * Cnr)^T * f_master
     *
     * This transformation is only valid when mapping forces (conservative mapping)
     *
     * Function layout:
     *
     * 1. Initialize auxiliary arrays
     *
     * 2. Compute the transformation matrix corresponding to the isogeometric mortar-based mapping
     *
     * 3. Tranpose multiply the right-hand side with the mortar tranformation matrix
     *
     * 4. Delete pointers
     */

    // 1. Initialize auxiliary arrays
    int size_N = couplingMatrices->getSizeN();
    double* tmpVec = new double[size_N];

    // 2. Compute the transformation matrix corresponding to the isogeometric mortar-based mapping
    couplingMatrices->getCnn()->solve(tmpVec, const_cast<double *>(_masterField));

    // 3. Tranpose multiply the right-hand side with the mortar tranformation matrix
    couplingMatrices->getCnr()->transposeMulitplyVec(tmpVec, _slaveField, size_N);

    // 4. Delete pointers
    delete[] tmpVec;
}

void IGAMortarMapper::computeErrorsConsistentMapping(const double* _slaveField, const double *_masterField) {
    /*
     * Computes the domain, interface and boundary errors corresponding to the isogeometric mortar-based mapping
     *
     * Function layout:
     *
     * 1. Initialize auxiliary arrays
     *
     * 2. Compute the relative error in terms of the L2 norm of the error in the domain
     *
     * 3. Compute the L2 norm of the error along the Dirichlet boundary
     *
     * 4. Compute the L2 norm of the error along the patch interfaces
     *
     * 5. Print messages on the errors
     */

    // 1. Initialize auxiliary arrays
    double errorL2Domain;
    double errorL2Interface[2];
    double errorL2Curve[2];

    // 2. Compute the relative error in terms of the L2 norm of the error in the domain
    if(propErrorComputation.isDomainError)
        errorL2Domain = computeDomainErrorInL2Norm4ConsistentMapping(_slaveField, _masterField);

    // 3. Compute the L2 norm of the error along the Dirichlet boundary
    if(propErrorComputation.isCurveError)
        if(isMappingIGA2FEM){
            computeIGADirichletCurveErrorInL2Norm(errorL2Curve, _slaveField);
        }else{
            computeIGADirichletCurveErrorInL2Norm(errorL2Curve, _masterField);
        }

    // 4. Compute the L2 norm of the error along the patch interfaces
    if(propErrorComputation.isInterfaceError)
        if(isMappingIGA2FEM){
            computeIGAPatchInterfaceErrorInL2Norm(errorL2Interface, _slaveField);
        }else{
            computeIGAPatchInterfaceErrorInL2Norm(errorL2Interface, _masterField);
        }

    // 5. Print messages on the errors
    printErrorMessage(infoOut, errorL2Domain, errorL2Curve, errorL2Interface);
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
     * Function layout:
     *
     * 1. Initialize auxiliary arrays
     *
     * 2. Loop over all the Gauss Points
     * ->
     *    2i. Get the Gauss Point Weight
     *   2ii. Get the product of the Jacobian transformations at the Gauss point
     *  2iii. Get the number of the basis functions of the finite element
     *   2iv. Initialize the field on the finite element mesh at the Gauss point
     *    2v. Loop over the nodes of the finite element
     *    ->
     *        2v.1. Get the value of the basis function
     *        2v.2. Get the index of the node
     *    <-
     *   2vi. Get the number of basis functions of the isogeometric discretization
     *  2vii. Initialize the field on the isogeometric discretization at the Gauss point
     * 2viii. Loop over the Control Points of the isogeometric discretization
     * ->
     *        2viii.1. Get the value of the basis function
     *        2viii.2. Get the index of the CP
     * <-
     *   2ix. Compute the difference of the vectors on the Gauss Point
     *    2x. Compute the norm of the difference of the fields at the Gauss Point
     *   2xi. Compute the norm of the reference field at the Gauss Point
     *  2xii. Add the contributions from the Gauss Point
     * <-
     *
     * 3. Compute the relative L2 norm of the mapping error
     *
     * 4. Return the relative L2 norm of the mapping error
     */

    // 1. Initialize auxiliary arrays
    double errorL2Domain = 0.0;
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
    int indexCP;
    int noCoord = 3;
    double tolNormSlaveField = 1e-6;

    // 2. Loop over all the Gauss Points
    for(int iGP = 0; iGP < streamGPs.size(); iGP++){
        // 2i. Get the Gauss Point Weight
        GW = streamGPs[iGP][0];

        // 2ii. Get the product of the Jacobian transformations at the Gauss point
        JacobianProducts = streamGPs[iGP][1];

        // 2iii. Get the number of the basis functions of the finite element
        noNodesFE = streamGPs[iGP][2];

        // 2iv. Initialize the field on the finite element mesh at the Gauss point
        for(int iCoord = 0; iCoord < noCoord; iCoord++)
            fieldFEM[iCoord] = 0.0;

        // 2v. Loop over the nodes of the finite element
        for(int iNodesFE = 0; iNodesFE < noNodesFE; iNodesFE++){
            // 2v.1. Get the value of the basis function
            basisFctFEM = streamGPs[iGP][3 + 2*iNodesFE + 1];

            // 2v.2. Get the index of the node
            indexNode = streamGPs[iGP][3 + 2*iNodesFE];
            for(int iCoord = 0; iCoord < noCoord; iCoord++)
                if(!isMappingIGA2FEM) {
                    fieldFEM[iCoord] += basisFctFEM*_slaveField[noCoord*indexNode + iCoord];
                } else
                    fieldFEM[iCoord] += basisFctFEM*_masterField[noCoord*indexNode + iCoord];
        }

        // 2vi. Get the number of basis functions of the isogeometric discretization
        noCPsIGA = streamGPs[iGP][3 + 2*noNodesFE];

        // 2vii. Initialize the field on the isogeometric discretization at the Gauss point
        for(int iCoord = 0; iCoord < noCoord; iCoord++)
            fieldIGA[iCoord] = 0.0;

        // 2viii. Loop over the Control Points of the isogeometric discretization
        for(int iCPsIGA = 0; iCPsIGA < noCPsIGA; iCPsIGA++){
            // 2viii.1. Get the value of the basis function
            basisFctIGA = streamGPs[iGP][3 + 2*noNodesFE + 2*iCPsIGA + 2];

            // 2viii.2. Get the index of the CP
            indexCP = streamGPs[iGP][3 + 2*noNodesFE + 2*iCPsIGA + 1];
            for(int iCoord = 0; iCoord < noCoord; iCoord++)
                if(!isMappingIGA2FEM)
                    fieldIGA[iCoord] += basisFctIGA*_masterField[noCoord*indexCP + iCoord];
                else
                    fieldIGA[iCoord] += basisFctIGA*_slaveField[noCoord*indexCP + iCoord];
        }

        // 2ix. Compute the difference of the vectors on the Gauss Point
        for(int iCoord = 0; iCoord < noCoord; iCoord++)
            errorVct[iCoord] = fieldFEM[iCoord] - fieldIGA[iCoord];

        // 2x. Compute the norm of the difference of the fields at the Gauss Point
        errorGPSquare = EMPIRE::MathLibrary::computeDenseDotProduct(noCoord, errorVct, errorVct);

        // 2xi. Compute the norm of the reference field at the Gauss Point
        if(!isMappingIGA2FEM)
            slaveFieldNorm = EMPIRE::MathLibrary::computeDenseDotProduct(noCoord, fieldFEM, fieldFEM);
        else
            slaveFieldNorm = EMPIRE::MathLibrary::computeDenseDotProduct(noCoord, fieldIGA, fieldIGA);

        // 2xii. Add the contributions from the Gauss Point
        errorL2Domain += errorGPSquare*JacobianProducts*GW;
        slaveFieldL2Domain += slaveFieldNorm*JacobianProducts*GW;
    }

    // 3. Compute the relative L2 norm of the mapping error
    errorL2Domain = sqrt(errorL2Domain);
    slaveFieldL2Domain = sqrt(slaveFieldL2Domain);
    if (slaveFieldL2Domain > tolNormSlaveField)
        errorL2Domain /= slaveFieldL2Domain;
    else
        WARNING_OUT() << "The norm of the slave field is smaller than the tolerance, no division of the mapping error is made" << std::endl;

    // 4. Return the relative L2 norm of the mapping error
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
     *
     * Function layout :
     *
     * 1. Initialize auxiliary arrays
     *
     * 2. Loop over all the Gauss Points
     * ->
     *    2i. Get the element length on the Gauss Point
     *   2ii. Get the number of basis functions of patch
     *  2iii. Initialize the field and its rotations on the isogeometric discretization at the Gauss point
     *   2iv. Loop over all basis functions of the patch and compute the field at the Gauss Point
     *   ->
     *        2iv.1. Get the index of the Control Point
     *        2iv.2. Get the value of the basis function
     *        2iv.3. Loop over all the Cartesian coordinates
     *   <-
     *    2v. Loop over all the DOFs of the patch and compute the rotations of the field at the Gauss point
     *    ->
     *        2v.1. Get the index of the DOF
     *        2v.2. Get the value of the B-operator for the bending rotation
     *        2v.3. Get the value of the B-operator for the twisting rotation
     *        2v.4. Compute the tangent and the bending rotations
     *    <-
     *   2vi. Compute the error vector for the displacements
     *  2vii. Compute the error in terms of the rotations
     * <-
     *
     * 3. Take the necessary for the norm square roots
     */

    // 1. Initialize auxiliary arrays
    for(int i = 0; i < 2; i++)
        _errorL2Curve[i] = 0.0;
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
    double errorTwistingRotation = 0.0;
    double normRotationSquare;
    double normErrorFieldSquare;
    double errorField[3];

    // 2. Loop over all the Gauss Points
    for(int iGP = 0; iGP < streamCurveGPs.size(); iGP++){
        // 2i. Get the element length on the Gauss Point
        elementLengthOnGP = streamCurveGPs[iGP][0];

        // 2ii. Get the number of basis functions of patch
        noCPs = streamCurveGPs[iGP][1];
        noDOFs = 3*noCPs;

        // 2iii. Initialize the field and its rotations on the isogeometric discretization at the Gauss point
        for(int iCoord = 0; iCoord < noCoord; iCoord++)
            field[iCoord] = 0.0;
        omegaT = 0.0;
        omegaN = 0.0;

        // 2iv. Loop over all basis functions of the patch and compute the field at the Gauss Point
        for(int iBFs = 0; iBFs < noCPs; iBFs++){
            // 2iv.1. Get the index of the Control Point
            indexCP = streamCurveGPs[iGP][1 + 1 + 2*iBFs];

            // 2iv.2. Get the value of the basis function
            basisFct = streamCurveGPs[iGP][1 + 1 + 2*iBFs + 1];

            // 2iv.3. Loop over all the Cartesian coordinates
            for(int iCoord = 0; iCoord < noCoord; iCoord++)
                field[iCoord] += basisFct*_fieldIGA[noCoord*indexCP + iCoord];
        }

        // 2v. Loop over all the DOFs of the patch and compute the rotations of the field at the Gauss point
        for(int iDOFs = 0; iDOFs < noDOFs; iDOFs++){
            // 2v.1. Get the index of the DOF
            indexDOF = streamCurveGPs[iGP][1 + 1 + 2*noCPs + 3*iDOFs];

            // 2v.2. Get the value of the B-operator for the bending rotation
            BoperatorT = streamCurveGPs[iGP][1 + 1 + 2*noCPs + 3*iDOFs + 1];

            // 2v.3. Get the value of the B-operator for the twisting rotation
            BoperatorN = streamCurveGPs[iGP][1 + 1 + 2*noCPs + 3*iDOFs + 2];

            // 2v.4. Compute the tangent and the bending rotations
            omegaT += BoperatorT*_fieldIGA[indexDOF];
            if (propWeakCurveDirichletConditions.isSecTwistingPrescribed)
                omegaN += BoperatorN*_fieldIGA[indexDOF];
        }

        // 2vi. Compute the error vector for the displacements
        for(int iCoord = 0; iCoord < noCoord; iCoord++){
            errorField[iCoord] = field[iCoord] - 0.0;
        }
        normErrorFieldSquare = EMPIRE::MathLibrary::computeDenseDotProduct(noCoord, errorField, errorField);
        _errorL2Curve[0] += normErrorFieldSquare*elementLengthOnGP;

        // 2vii. Compute the error in terms of the rotations
        errorBendingRotation = omegaT + 0.0;
        errorTwistingRotation = omegaN + 0.0;
        if (propWeakCurveDirichletConditions.isSecTwistingPrescribed)
            normRotationSquare = errorBendingRotation*errorBendingRotation + errorTwistingRotation*errorTwistingRotation;
        else
            normRotationSquare = errorBendingRotation*errorBendingRotation;
        _errorL2Curve[1] += normRotationSquare*elementLengthOnGP;
    }

    // 3. Take the necessary for the norm square roots
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
     *
     * Functions layout:
     *
     * 1. Initialize auxiliary arrays
     *
     * 2. Loop over all the Gauss Points
     * ->
     *    2i. Get the element length on the Gauss Point
     *   2ii. Get the number of basis functions of patch I
     *  2iii. Initialize the field and its rotations on the isogeometric discretization at the Gauss point
     *   2iv. Loop over all basis functions of patch I and compute the field at the Gauss Point
     *   ->
     *        2iv.1. Get the index of the Control Point
     *        2iv.2. Get the value of the basis function
     *        2iv.3. Loop over all the Cartesian coordinates
     *   <-
     *    2v. Loop over all the DOFs of patch I and compute the rotations of the field at the Gauss point
     *    ->
     *        2v.1. Get the index of the DOF
     *        2v.2. Get the value of the B-operator for the bending rotation
     *        2v.3. Get the value of the B-operator for the twisting rotation
     *        2v.4. Compute the tangent and the bending rotations
     *    <-
     *   2vi. Get the number of basis functions of patch J
     *  2vii. Loop over all basis functions of patch J and compute the field at the Gauss Point
     *  ->
     *        2vii.1. Get the index of the Control Point
     *        2vii.2. Get the value of the basis function
     *        2vii.3. Loop over all the Cartesian coordinates
     *  <-
     * 2viii. Loop over all the DOFs of patch J and compute the rotations of the field at the Gauss point
     * ->
     *        2viii.1. Get the index of the DOF
     *        2viii.2. Get the value of the B-operator for the bending rotation
     *        2viii.3. Get the value of the B-operator for the twisting rotation
     *        2viii.4. Compute the tangent and the bending rotations
     * <-
     *   2ix. Get the factors for the bending and the twisting rotation
     *    2x. Compute the error vector for the displacements
     *   2xi. Compute the error in terms of the rotations
     * <-
     *
     * 3. Take the necessary for the norm square roots
     */

    // 1. Initialize auxiliary arrays
    for(int i = 0; i < 2; i++)
        _errorL2Interface[i] = 0.0;
    int noCPsI;
    int noCPsJ;
    int noDOFsI;
    int noDOFsJ;
    int indexCP;
    int indexDOF;
    int noCoord = 3;
    double factorTangent;
    double factorNormal;
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
    double errorTwistingRotation = 0.0;
    double normRotationSquare;
    double normErrorFieldSquare;
    double errorField[3];

    // 2. Loop over all the Gauss Points
    for(int iGP = 0; iGP < streamInterfaceGPs.size(); iGP++){
        // 2i. Get the element length on the Gauss Point
        elementLengthOnGP = streamInterfaceGPs[iGP][0];

        // 2ii. Get the number of basis functions of patch I
        noCPsI = streamInterfaceGPs[iGP][1];
        noDOFsI = 3*noCPsI;

        // 2iii. Initialize the field and its rotations on the isogeometric discretization at the Gauss point
        for(int iCoord = 0; iCoord < noCoord; iCoord++){
            fieldI[iCoord] = 0.0;
            fieldJ[iCoord] = 0.0;
        }
        omegaTI = 0.0;
        omegaTJ = 0.0;
        omegaNI = 0.0;
        omegaNJ = 0.0;

        // 2iv. Loop over all basis functions of patch I and compute the field at the Gauss Point
        for(int iBFs = 0; iBFs < noCPsI; iBFs++){
            // 2iv.1. Get the index of the Control Point
            indexCP = streamInterfaceGPs[iGP][1 + 1 + 2*iBFs];

            // 2iv.2. Get the value of the basis function
            basisFct = streamInterfaceGPs[iGP][1 + 1 + 2*iBFs + 1];

            // 2iv.3. Loop over all the Cartesian coordinates
            for(int iCoord = 0; iCoord < noCoord; iCoord++)
                fieldI[iCoord] += basisFct*_fieldIGA[noCoord*indexCP + iCoord];
        }

        // 2v. Loop over all the DOFs of patch I and compute the rotations of the field at the Gauss point
        for(int iDOFs = 0; iDOFs < noDOFsI; iDOFs++){
            // 2v.1. Get the index of the DOF
            indexDOF = streamInterfaceGPs[iGP][1 + 1 + 2*noCPsI + 3*iDOFs];

            // 2v.2. Get the value of the B-operator for the bending rotation
            BoperatorT = streamInterfaceGPs[iGP][1 + 1 + 2*noCPsI + 3*iDOFs + 1];

            // 2v.3. Get the value of the B-operator for the twisting rotation
            BoperatorN = streamInterfaceGPs[iGP][1 + 1 + 2*noCPsI + 3*iDOFs + 2];

            // 2v.4. Compute the tangent and the bending rotations
            omegaTI += BoperatorT*_fieldIGA[indexDOF];
            omegaNI += BoperatorN*_fieldIGA[indexDOF];
        }

        // 2vi. Get the number of basis functions of patch J
        noCPsJ = streamInterfaceGPs[iGP][1 + 1 + 2*noCPsI + 3*noDOFsI];
        noDOFsJ = 3*noCPsJ;

        // 2vii. Loop over all basis functions of patch J and compute the field at the Gauss Point
        for(int iBFs = 0; iBFs < noCPsJ; iBFs++){
            // 2vii.1. Get the index of the Control Point
            indexCP = streamInterfaceGPs[iGP][1 + 1 + 2*noCPsI + 3*noDOFsI + 1 + 2*iBFs];

            // 2vii.2. Get the value of the basis function
            basisFct = streamInterfaceGPs[iGP][1 + 1 + 2*noCPsI + 3*noDOFsI + 1 + 2*iBFs + 1];

            // 2vii.3. Loop over all the Cartesian coordinates
            for(int iCoord = 0; iCoord < noCoord; iCoord++)
                fieldJ[iCoord] += basisFct*_fieldIGA[noCoord*indexCP + iCoord];
        }

        // 2viii. Loop over all the DOFs of patch J and compute the rotations of the field at the Gauss point
        for(int iDOFs = 0; iDOFs < noDOFsJ; iDOFs++){
            // 2viii.1. Get the index of the DOF
            indexDOF = streamInterfaceGPs[iGP][1 + 1 + 2*noCPsI + 3*noDOFsI + 1 + 2*noCPsJ + 3*iDOFs];

            // 2viii.2. Get the value of the B-operator for the bending rotation
            BoperatorT = streamInterfaceGPs[iGP][1 + 1 + 2*noCPsI + 3*noDOFsI + 1 + 2*noCPsJ + 3*iDOFs + 1];

            // 2viii.3. Get the value of the B-operator for the twisting rotation
            BoperatorN = streamInterfaceGPs[iGP][1 + 1 + 2*noCPsI + 3*noDOFsI + 1 + 2*noCPsJ + 3*iDOFs + 2];

            // 2viii.4. Compute the tangent and the bending rotations
            omegaTJ += BoperatorT*_fieldIGA[indexDOF];
            omegaNJ += BoperatorN*_fieldIGA[indexDOF];
        }

        // 2ix. Get the factors for the bending and the twisting rotation
        factorTangent = streamInterfaceGPs[iGP][1 + 1 + 2*noCPsI + 3*noDOFsI + 1 + 2*noCPsJ + 3*noDOFsJ];
        factorNormal = streamInterfaceGPs[iGP][1 + 1 + 2*noCPsI + 3*noDOFsI + 1 + 2*noCPsJ + 3*noDOFsJ + 1];

        // 2x. Compute the error vector for the displacements
        for(int iCoord = 0; iCoord < noCoord; iCoord++){
            errorField[iCoord] = fieldI[iCoord] - fieldJ[iCoord];
        }
        normErrorFieldSquare = EMPIRE::MathLibrary::computeDenseDotProduct(noCoord, errorField, errorField);
        _errorL2Interface[0] += normErrorFieldSquare*elementLengthOnGP;

        // 2xi. Compute the error in terms of the rotations
        errorBendingRotation = omegaTI + factorTangent*omegaTJ;
        if (propWeakPatchContinuityConditions.isSecTwistingCoupled)
            errorTwistingRotation = omegaNI + factorNormal*omegaNJ;
        normRotationSquare = errorBendingRotation*errorBendingRotation + errorTwistingRotation*errorTwistingRotation;
        _errorL2Interface[1] += normRotationSquare*elementLengthOnGP;
    }

    // 3. Take the necessary for the norm square roots
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
     *
     * Function layout:
     *
     * 1. Print message
     *
     * 2. Initialize auxiliary arrays
     *
     * 3. Create a field of ones to be mapped with the isogeometric mortar-based method
     *
     * 4. Compute the mapped field using the isogeometric mortar-based method
     *
     * 5. Compute the norm of the mapped field
     *
     * 6. Replace each row which's sum is not equal to 1.0 +- epsilon by sum value of Cnr
     * ### If the array of the inconsistent DOFs is not zero ###
     *    6i. Print info message
     *   6ii. Factorize Cnn
     *  6iii. Perform the consistent mapping using the unit field after enforcing consistency
     *   6iv. Compute the norm of the mapped field
     * ### If the array of the inconsistent DOFs not zero ###
     *
     * 7. Compute the norm of the mapped field
     *
     * 8. Print info message on the consistency check
     *
     * 9. Delete pointers
     */

    // 1. Print message
    INFO_OUT() << "Checking consistency" << std::endl;

    // 2. Initialize auxiliary arrays
    int denom;
    int size_N = couplingMatrices->getSizeN();
    int size_R = couplingMatrices->getSizeR();
    double norm;

    // 3. Create a field of ones to be mapped with the isogeometric mortar-based method
    double* ones = new double[size_R];
    for(int i = 0; i < size_R; i++) {
        ones[i] = 1.0;
    }

    // 4. Compute the mapped field using the isogeometric mortar-based method
    double* output = new double[size_N];
    this->consistentMapping(ones, output);

    // 5. Compute the norm of the mapped field
    norm = 0;
    vector<int> inconsistentDoF;
    for(int i = 0; i < size_N; i++) {
        if(fabs(output[i] - 1) > propConsistency.tolConsistency && output[i] != 0)
            inconsistentDoF.push_back(i);
        norm += output[i]*output[i];
    }

    // 6. Replace each row which's sum is not equal to 1.0 +- epsilon by sum value of Cnr
    if(!inconsistentDoF.empty()) { // ### If the array of the inconsistent DOFs is not zero ###
        // 6i. Print info message
        INFO_OUT() << "Mapping found inconsistent up to specified tolerance" << std::endl;
        INFO_OUT() << "inconsistent DOF size = " << inconsistentDoF.size() << std::endl;
        INFO_OUT() << "Enforcing consistency" << std::endl;
        for(vector<int>::iterator it = inconsistentDoF.begin(); it != inconsistentDoF.end(); it++) {
            couplingMatrices->deleterow(*it);       // deleterow might not be working for the new coupling matrix datastructure
            couplingMatrices->addCNNValue(*it ,*it , couplingMatrices->getCnr()->getRowSum(*it));
        }

        // 6ii. Factorize Cnn
        couplingMatrices->factorizeCnn();

        // 6iii. Perform the consistent mapping using the unit field after enforcing consistency
        this->consistentMapping(ones, output);

        // 6iv. Compute the norm of the mapped field
        norm = 0;
        for(int i = 0; i < size_N; i++) {
            norm += output[i]*output[i];
        }
    } else // ### If the array of the inconsistent DOFs not zero ###
        INFO_OUT() << "Mapping found consistent up to specified tolerance" << std::endl;

    // 7. Compute the norm of the mapped field
    denom = size_N - couplingMatrices->getIndexEmptyRowCnn().size();
    norm = sqrt(norm/denom);

    // 8. Print info message on the consistency check
    DEBUG_OUT() << "Checking consistency after enforcement"<<endl;
    DEBUG_OUT() << "Deviation from unit field : "<< fabs(norm - 1.0) << endl;
    if(fabs(norm - 1.0) > propConsistency.tolConsistency) {
        ERROR_OUT() << "Coupling matrices not consistent up to " << propConsistency.tolConsistency << "!"<<endl;
        exit(-1);
    }

    // 9. Delete pointers
    delete ones;
    delete output;
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
    /*
     * Print message containing the errors in the L2 norm
     */
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
