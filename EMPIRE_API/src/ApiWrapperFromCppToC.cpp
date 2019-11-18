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
/***********************************************************************************************//**
 * \file ApiWrapperFromCppToC.cpp
 * This file wraps the C++ API in a C API
 * \date 2/22/2012
 **************************************************************************************************/
#include "Empire.h"
#include "EMPIRE_API.h"
#include <assert.h>
#include <string.h>
#include <string>
#include <vector>
#include <cfloat>
#include <iostream>

using namespace EMPIRE;
using namespace std;

/// Lets make a EMPIRE_API object in the global scope
Empire *empire;

void EMPIRE_API_Connect(char* inputFileName) {
    empire = new Empire();
    empire->initEnvironment(inputFileName);
    empire->connect();
}

const char *EMPIRE_API_getUserDefinedText(char *elementName) {
    static string text;
    text = empire->getUserDefinedText(elementName);
    return text.data();
}

void EMPIRE_API_sendMesh(char *name, int numNodes, int numElems, double *nodes, int *nodeIDs,
        int *numNodesPerElem, int *elems) {
    empire->sendMesh(numNodes, numElems, nodes, nodeIDs, numNodesPerElem, elems);
}

void EMPIRE_API_sendSectionMesh(char *name, int numNodes, int numElems, double *nodes, int *nodeIDs,
        int *numNodesPerElem, int *elems, int numSections, int numRootSectionNodes,
        int numNormalSectionNodes, int numTipSectionNodes, double *rotationGlobal2Root,
        double *translationGlobal2Root) {
    empire->sendSectionMesh(numNodes, numElems, nodes, nodeIDs, numNodesPerElem, elems, numSections,
            numRootSectionNodes, numNormalSectionNodes, numTipSectionNodes, rotationGlobal2Root,
            translationGlobal2Root);
}

void EMPIRE_API_recvMesh(char *name, int *numNodes, int *numElems, double **nodes, int **nodeIDs,
        int **numNodesPerElem, int **elem) {
    empire->recvMesh(numNodes, numElems, nodes, nodeIDs, numNodesPerElem, elem);
}

void EMPIRE_API_sendIGAPatch(int _pDegree, int _uNumKnots, double* _uKnotVector, int _qDegree,
        int _vNumKnots, double* _vKnotVector, int _uNumControlPoints, int _vNumControlPoints,
        double* _cpNet, int* _nodeNet) {
    empire->sendIGAPatch(_pDegree, _uNumKnots, _uKnotVector, _qDegree, _vNumKnots, _vKnotVector,
            _uNumControlPoints, _vNumControlPoints, _cpNet, _nodeNet);
}

void EMPIRE_API_sendIGAMesh(char *_name, int _numPatches, int _numNodes) {
    empire->sendIGAMesh(_numPatches, _numNodes);
}

void EMPIRE_API_sendIGATrimmingInfo(int _isTrimmed, int _numLoops) {
    empire->sendIGATrimmingInfo(_isTrimmed, _numLoops);
}

void EMPIRE_API_sendIGATrimmingLoopInfo(int _inner, int _numCurves) {
    empire->sendIGATrimmingLoopInfo(_inner, _numCurves);
}

void EMPIRE_API_sendIGATrimmingCurve(int _direction, int _pDegree, int _uNumKnots,
        double* _uKnotVector, int _uNumControlPoints, double* _cpNet) {
    empire->sendIGATrimmingCurve(_direction, _pDegree, _uNumKnots, _uKnotVector, _uNumControlPoints,
            _cpNet);
}

void EMPIRE_API_sendIGANumDirichletConditions(int _numDirichletConditions){
    empire->sendIGANumDirichletConditions(_numDirichletConditions);
}

void EMPIRE_API_sendIGADirichletConditionInfo(int _patchCtr, int _patchBLCtr, int _patchBLTrCurveCtr,
                                           int _isGPProvided) {
    empire->sendIGADirichletConditionInfo(_patchCtr, _patchBLCtr, _patchBLTrCurveCtr,
                                       _isGPProvided);
}

void EMPIRE_API_sendIGADirichletConditionData(int _trCurveNumGP,
                                           double *_trCurveGPs, double *_trCurveGPWeights,
                                           double *_trCurveGPTangents,
                                           double *_trCurveGPJacobianProducts) {
    empire->sendIGADirichletConditionData(_trCurveNumGP,
                                       _trCurveGPs, _trCurveGPWeights,
                                       _trCurveGPTangents,
                                       _trCurveGPJacobianProducts);
}

void EMPIRE_API_sendIGANumPatchConnections(int _numConnections){
    empire->sendIGANumPatchConnections(_numConnections);
}

void EMPIRE_API_sendIGAPatchConnectionInfo(int _masterPatchCtr, int _masterPatchBLCtr, int _masterPatchBLTrCurveCtr,
                                           int _slavePatchCtr, int _slavePatchBLCtr, int _slavePatchBLTrCurveCtr,
                                           int _isGPProvided) {
    empire->sendIGAPatchConnectionInfo(_masterPatchCtr, _masterPatchBLCtr, _masterPatchBLTrCurveCtr,
                                       _slavePatchCtr, _slavePatchBLCtr, _slavePatchBLTrCurveCtr,
                                       _isGPProvided);
}

void EMPIRE_API_sendIGAPatchConnectionData(int _trCurveNumGP,
                                           double *_trCurveMasterGPs, double *_trCurveSlaveGPs, double *_trCurveGPWeights,
                                           double *_trCurveMasterGPTangents, double *_trCurveSlaveGPTangents,
                                           double *_trCurveGPJacobianProducts) {
    empire->sendIGAPatchConnectionData(_trCurveNumGP,
                                       _trCurveMasterGPs, _trCurveSlaveGPs, _trCurveGPWeights,
                                       _trCurveMasterGPTangents, _trCurveSlaveGPTangents,
                                       _trCurveGPJacobianProducts);
}

void EMPIRE_API_sendIGAPatchCouplingGaussPointsTest(int numPatches, int* numBRepsPerPatch, int totalNumGP, int totalNumBRePs,
        int* allID_slave, int* allNumElemsPerBRep, int* allNumGPsPerElem,
        double* allGPOfBRep_master, double* allGPOfBRep_slave, double* allGPOfBRep_weight,
        double* allTangents_master, double* allTangents_slave, double* allMappings) {

    empire->sendIGAPatchCouplingGaussPointsTest(numPatches, numBRepsPerPatch, totalNumGP, totalNumBRePs,
                allID_slave, allNumElemsPerBRep, allNumGPsPerElem,
                allGPOfBRep_master, allGPOfBRep_slave, allGPOfBRep_weight,
                allTangents_master, allTangents_slave, allMappings);
}

void EMPIRE_API_sendIGADirichletDofs(int numberOfClampedDofs, int* clampedDofs, int clampedDirections) {
    empire->sendIGADirichletDofs(numberOfClampedDofs, clampedDofs, clampedDirections);
}

void EMPIRE_API_sendDataField(char *name, int sizeOfArray, double *dataField) {
    empire->sendDataField(sizeOfArray, dataField);
}

void EMPIRE_API_recvDataField(char *name, int sizeOfArray, double *dataField) {
    empire->recvDataField(sizeOfArray, dataField);
}

void EMPIRE_API_sendSignal_double(char *name, int sizeOfArray, double *signal) {
    empire->sendSignal_double(name, sizeOfArray, signal);
}

void EMPIRE_API_recvSignal_double(char *name, int sizeOfArray, double *signal) {
    empire->recvSignal_double(name, sizeOfArray, signal);
}

void EMPIRE_API_sendConvergenceSignal(int signal) {
    empire->sendConvergenceSignal(signal);
}

int EMPIRE_API_recvConvergenceSignal() {
    return empire->recvConvergenceSignal();
}

void EMPIRE_API_printDataField(char *name, int sizeOfArray, double *dataField) {
    empire->printDataField(name, sizeOfArray, dataField);
}

void EMPIRE_API_Disconnect() {
    empire->disconnect();
    delete empire;
}
