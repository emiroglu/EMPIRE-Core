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
#include <assert.h>
#include <string.h>
#include "ClientCode.h"
#include "ServerCommunication.h"
#include "DataField.h"
#include "AbstractMesh.h"
#include "FEMesh.h"
#include "SectionMesh.h"
#include "IGAMesh.h"	 	
#include "IGAPatchSurface.h"
#include "WeakIGADirichletCurveCondition.h"
#include "WeakIGAPatchContinuityCondition.h"
#include "Signal.h"
#include "Message.h"
#include <stdlib.h>
#include <algorithm>

namespace EMPIRE {

using namespace std;

ClientCode::ClientCode(string _name) :
        name(_name), serverComm(NULL), nameToMeshMap(), nameToSignalMap() {
}

ClientCode::~ClientCode() {
    for (map<string, AbstractMesh*>::iterator it = nameToMeshMap.begin(); it != nameToMeshMap.end();
            it++) {
        delete it->second;
    }
    for (map<string, Signal*>::iterator it = nameToSignalMap.begin(); it != nameToSignalMap.end();
            it++) {
        delete it->second;
    }
}

void ClientCode::setServerCommunication(ServerCommunication *_serverComm) {
    assert(_serverComm != NULL);
    serverComm = _serverComm;
}

void ClientCode::recvFEMesh(std::string meshName, bool triangulateAll) {
    assert(serverComm != NULL);
    assert(nameToMeshMap.find(meshName) == nameToMeshMap.end());

    const int BUFFER_SIZE = 2;
    int meshInfo[BUFFER_SIZE]; // number of nodes, number of elements, number of nodes per element
    { // output to shell
        HEADING_OUT(3, "ClientCode", "receiving mesh (" + meshName + ") from [" + name + "]...",
                infoOut);
    }
    serverComm->receiveFromClientBlocking<int>(name, BUFFER_SIZE, meshInfo);
    int numNodes = meshInfo[0];
    int numElems = meshInfo[1];

    FEMesh *mesh = new FEMesh(meshName, numNodes, numElems, triangulateAll);
    serverComm->receiveFromClientBlocking<double>(name, numNodes * 3, mesh->nodes);
	serverComm->receiveFromClientBlocking<int>(name, numNodes, mesh->nodeIDs);
    serverComm->receiveFromClientBlocking<int>(name, numElems, mesh->numNodesPerElem);
    mesh->initElems();
    serverComm->receiveFromClientBlocking<int>(name, mesh->elemsArraySize, mesh->elems);
	nameToMeshMap.insert(pair<string, AbstractMesh*>(meshName, mesh));
    { // output to shell
        DEBUG_OUT() << (*mesh) << endl;
        mesh->computeBoundingBox();
        INFO_OUT() << mesh->boundingBox << endl;
    }
}

void ClientCode::recvSectionMesh(std::string meshName, bool triangulateAll) {
    assert(serverComm != NULL);
    assert(nameToMeshMap.find(meshName) == nameToMeshMap.end());

    const int BUFFER_SIZE = 2;
    int meshInfo[BUFFER_SIZE]; // number of nodes, number of elements, number of nodes per element
    { // output to shell
        HEADING_OUT(3, "ClientCode", "receiving mesh (" + meshName + ") from [" + name + "]...",
                infoOut);
    }
    serverComm->receiveFromClientBlocking<int>(name, BUFFER_SIZE, meshInfo);
    int numNodes = meshInfo[0];
    int numElems = meshInfo[1];

    SectionMesh *mesh = new SectionMesh(meshName, numNodes, numElems, triangulateAll);
    serverComm->receiveFromClientBlocking<double>(name, numNodes * 3, mesh->nodes);
    serverComm->receiveFromClientBlocking<int>(name, numNodes, mesh->nodeIDs);
    serverComm->receiveFromClientBlocking<int>(name, numElems, mesh->numNodesPerElem);
    mesh->initElems();
    serverComm->receiveFromClientBlocking<int>(name, mesh->elemsArraySize, mesh->elems);

    const int SECTION_INFO_SIZE = 4;
    int sectionInfo[4]; // number of sections, number of root section nodes, number of normal section nodes, number of tip section nodes
    serverComm->receiveFromClientBlocking<int>(name, SECTION_INFO_SIZE, sectionInfo);
    mesh->setNumSections(sectionInfo[0]);
    mesh->setNumRootSectionNodes(sectionInfo[1]);
    mesh->setNumNormalSectionNodes(sectionInfo[2]);
    mesh->setNumTipSectionNodes(sectionInfo[3]);

    double rotationGlobal2Root[9];
    double translationGlobal2Root[3];
    serverComm->receiveFromClientBlocking<double>(name, 9, rotationGlobal2Root);
    serverComm->receiveFromClientBlocking<double>(name, 3, translationGlobal2Root);
    mesh->setRotationGlobal2Root(rotationGlobal2Root);
    mesh->setTranslationGlobal2Root(translationGlobal2Root);

    nameToMeshMap.insert(pair<string, AbstractMesh*>(meshName, mesh));
    { // output to shell
        DEBUG_OUT() << (*mesh) << endl;
        mesh->computeBoundingBox();
        INFO_OUT() << mesh->boundingBox << endl;
    }
}

void ClientCode::sendMesh(std::string meshName){
    assert(serverComm != NULL);
//    assert(nameToMeshMap.find(meshName) == nameToMeshMap.end());

    AbstractMesh *mesh = this->getMeshByName(meshName);
    if (mesh->type == EMPIRE_Mesh_FEMesh){
        FEMesh* femesh = dynamic_cast<FEMesh*>( mesh );
        int numNodes2send = femesh->numNodes;
        int *numNodes2sendPtr = &numNodes2send;
        int numElems2send = femesh->numElems;
        int *numElems2sendPtr = &numElems2send;
        serverComm->sendToClientBlocking<int>(name, 1, numNodes2sendPtr);
        serverComm->sendToClientBlocking<int>(name, 1, numElems2sendPtr);
        serverComm->sendToClientBlocking<double>(name, femesh->numNodes * 3, femesh->nodes);
        serverComm->sendToClientBlocking<int>(name, femesh->numNodes, femesh->nodeIDs);
        serverComm->sendToClientBlocking<int>(name, femesh->numElems, femesh->numNodesPerElem);
        int count = 0;
        for (int i = 0; i < numElems2send; i++)
            count += femesh->numNodesPerElem[i];
        serverComm->sendToClientBlocking<int>(name, count, femesh->elems);
     } else{
        DEBUG_OUT() << "Not implemented for IGA yet" << endl;
        assert(false);
    }
}

void ClientCode::copyMesh(std::string meshName, AbstractMesh *meshToCopyFrom) {
    assert(nameToMeshMap.find(meshName) == nameToMeshMap.end());
    if (meshToCopyFrom->type == EMPIRE_Mesh_FEMesh){
        { // output to shell
            HEADING_OUT(3, "ClientCode", "copying mesh (" + meshToCopyFrom->name + ") to (" + meshName + ") of ["+name+"]...",
                    infoOut);
        }
        FEMesh* femesh = dynamic_cast<FEMesh*>( meshToCopyFrom );
        int numNodes = femesh->numNodes;
        int numElems = femesh->numElems;
        FEMesh *copyMesh = new FEMesh(meshName, numNodes, numElems, false); //triangulateAll=false
        copyMesh->nodes = femesh->nodes;
        copyMesh->nodeIDs = femesh->nodeIDs;
        copyMesh->numNodesPerElem = femesh->numNodesPerElem;
        copyMesh->initElems();
        copyMesh->elems = femesh->elems;
        nameToMeshMap.insert(pair<string, AbstractMesh*>(meshName, copyMesh));
        { // output to shell
            DEBUG_OUT() << (*copyMesh) << endl;
            copyMesh->computeBoundingBox();
            INFO_OUT() << copyMesh->boundingBox << endl;
        }
    } else if (meshToCopyFrom->type == EMPIRE_Mesh_IGAMesh) {
        { // output to shell
            HEADING_OUT(3, "ClientCode", "copying mesh (" + meshToCopyFrom->name + ") to (" + meshName + ") of ["+name+"]...",
                    infoOut);
        }
        IGAMesh* igaMesh = dynamic_cast<IGAMesh*>( meshToCopyFrom );
        int numPatches = igaMesh->getNumPatches();
        int numNodes = igaMesh->getNumNodes();
        IGAMesh* copyMesh = new IGAMesh(meshName, numNodes);

        for (int patchCount = 0; patchCount < numPatches; patchCount++) {
            int pDegree = igaMesh->getSurfacePatch(patchCount)->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree();
            int uNoKnots = igaMesh->getSurfacePatch(patchCount)->getIGABasis()->getUBSplineBasis1D()->getNoKnots();
            int qDegree = igaMesh->getSurfacePatch(patchCount)->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree();
            int vNoKnots = igaMesh->getSurfacePatch(patchCount)->getIGABasis()->getVBSplineBasis1D()->getNoKnots();
            int uNoControlPoints = igaMesh->getSurfacePatch(patchCount)->getUNoControlPoints();
            int vNoControlPoints = igaMesh->getSurfacePatch(patchCount)->getVNoControlPoints();

            double* uKnotVector = new double[uNoKnots];
            double* vKnotVector = new double[vNoKnots];

            for (int iUNoKnots = 0;  iUNoKnots < uNoKnots; iUNoKnots++)
                uKnotVector[iUNoKnots] = igaMesh->getSurfacePatch(patchCount)->getIGABasis()->getUBSplineBasis1D()->getKnotVector()[iUNoKnots];
            for (int iVNoKnots = 0;  iVNoKnots < vNoKnots; iVNoKnots++)
                vKnotVector[iVNoKnots] = igaMesh->getSurfacePatch(patchCount)->getIGABasis()->getVBSplineBasis1D()->getKnotVector()[iVNoKnots];

            double* controlPointNet = new double[uNoControlPoints * vNoControlPoints * 4];
            int* dofIndexNet = new int[uNoControlPoints * vNoControlPoints];

            IGAControlPoint** IGAcontrolPointNet = igaMesh->getSurfacePatch(patchCount)->getControlPointNet();
            for (int iCP = 0; iCP < uNoControlPoints * vNoControlPoints; iCP++) {
                dofIndexNet[iCP] = IGAcontrolPointNet[iCP]->getDofIndex();
                controlPointNet[iCP*4] = IGAcontrolPointNet[iCP]->getX();
                controlPointNet[iCP*4+1] = IGAcontrolPointNet[iCP]->getY();
                controlPointNet[iCP*4+2] = IGAcontrolPointNet[iCP]->getZ();
                controlPointNet[iCP*4+3] = IGAcontrolPointNet[iCP]->getW();
            }
            IGAPatchSurface* thePatch = copyMesh->addPatch(pDegree, uNoKnots, uKnotVector, qDegree, vNoKnots, vKnotVector,
                    uNoControlPoints, vNoControlPoints, controlPointNet, dofIndexNet);

            // deleting the pointers after object creation
            delete uKnotVector;
            delete vKnotVector;
            delete controlPointNet;
            delete dofIndexNet;

            // Add and linearize the trimming curves
            int isTrimmed = igaMesh->getSurfacePatch(patchCount)->getTrimming().isTrimmed();
            int numLoops = igaMesh->getSurfacePatch(patchCount)->getTrimming().getNumOfLoops();
            if(isTrimmed) {
                assert(numLoops>0);
                //Get every loop
                const std::vector<int> outerLoopIndices = igaMesh->getSurfacePatch(patchCount)->getTrimming().getOutterLoopIndex();
                for(int loopCount = 0; loopCount < numLoops; loopCount++) {

                    int inner = 1;
                    if ( std::find(outerLoopIndices.begin(), outerLoopIndices.end(), loopCount) != outerLoopIndices.end() )
                        inner = 0;
                    int numCurves = igaMesh->getSurfacePatch(patchCount)->getTrimming().getLoop(loopCount).getNoCurves();

                    thePatch->addTrimLoop(inner, numCurves);

                    //Get every curve
                    for(int curveCount = 0; curveCount < numCurves; curveCount++) {

                        int direction = igaMesh->getSurfacePatch(patchCount)->getTrimming().getLoop(loopCount).getDirection(curveCount);
                        int pDegree = igaMesh->getSurfacePatch(patchCount)->getTrimming().getLoop(loopCount).getIGACurve(curveCount).getIGABasis()->getPolynomialDegree();
                        int uNoKnots = igaMesh->getSurfacePatch(patchCount)->getTrimming().getLoop(loopCount).getIGACurve(curveCount).getIGABasis()->getNoKnots();
                        int uNoControlPoints = igaMesh->getSurfacePatch(patchCount)->getTrimming().getLoop(loopCount).getIGACurve(curveCount).getNoControlPoints();

                        double uKnotVector[uNoKnots];
                        for (int iKnot = 0; iKnot < uNoKnots; iKnot++){
                            uKnotVector[iKnot] = igaMesh->getSurfacePatch(patchCount)->getTrimming().getLoop(loopCount).getIGACurve(curveCount).getIGABasis()->getKnotVector()[iKnot];
                        }
                        double controlPointNet[uNoControlPoints*4];
                        for (int iCP = 0; iCP < uNoControlPoints; iCP++){
                            controlPointNet[iCP*4] = igaMesh->getSurfacePatch(patchCount)->getTrimming().getLoop(loopCount).getIGACurve(curveCount).getControlPointNet()[iCP].getX();
                            controlPointNet[iCP*4+1] = igaMesh->getSurfacePatch(patchCount)->getTrimming().getLoop(loopCount).getIGACurve(curveCount).getControlPointNet()[iCP].getY();
                            controlPointNet[iCP*4+2] = igaMesh->getSurfacePatch(patchCount)->getTrimming().getLoop(loopCount).getIGACurve(curveCount).getControlPointNet()[iCP].getZ();
                            controlPointNet[iCP*4+3] = igaMesh->getSurfacePatch(patchCount)->getTrimming().getLoop(loopCount).getIGACurve(curveCount).getControlPointNet()[iCP].getW();
                        }
                        thePatch->addTrimCurve(direction, pDegree, uNoKnots, uKnotVector,
                                uNoControlPoints, controlPointNet);

                    } // end curve
                } // end trimming loops
                thePatch->linearizeTrimming();
            } // end isTrimmed
        } // end patch
        nameToMeshMap.insert(pair<string, AbstractMesh*>(meshName, copyMesh));
        { // output to shell
            DEBUG_OUT() << (*copyMesh) << endl;
            copyMesh->computeBoundingBox();
            INFO_OUT() << copyMesh->boundingBox << endl;
        }
    } else {
        ERROR_BLOCK_OUT("ClientCode","copyMesh","Not implemented for IGA yet");
    }
}

void ClientCode::recvIGAMesh(std::string meshName) {
    assert(serverComm != NULL);
    assert(nameToMeshMap.find(meshName) == nameToMeshMap.end());

    const int BUFFER_SIZE_MESH = 2;
    int meshInfo[BUFFER_SIZE_MESH];
    { // output to shell
        HEADING_OUT(3, "ClientCode", "receiving IGA Mesh (" + meshName + ") from [" + name + "]...",
                infoOut);
    }
    serverComm->receiveFromClientBlocking<int>(name, BUFFER_SIZE_MESH, meshInfo);

    int numPatches = meshInfo[0];
    int numNodes = meshInfo[1];

    IGAMesh* theIGAMesh = new IGAMesh(meshName, numNodes);

    const int BUFFER_SIZE_PATCH = 6;
    int patchInfo[BUFFER_SIZE_PATCH];
    for (int patchCount = 0; patchCount < numPatches; patchCount++) {

        serverComm->receiveFromClientBlocking<int>(name, BUFFER_SIZE_PATCH, patchInfo);

        int pDegree = patchInfo[0];
        int uNoKnots = patchInfo[1];
        int qDegree = patchInfo[2];
        int vNoKnots = patchInfo[3];
        int uNoControlPoints = patchInfo[4];
        int vNoControlPoints = patchInfo[5];

        double* uKnotVector = new double[uNoKnots];
        double* vKnotVector = new double[vNoKnots];
        double* controlPointNet = new double[uNoControlPoints * vNoControlPoints * 4];
        int* dofIndexNet = new int[uNoControlPoints * vNoControlPoints];

        serverComm->receiveFromClientBlocking<double>(name, uNoKnots, uKnotVector);
        serverComm->receiveFromClientBlocking<double>(name, vNoKnots, vKnotVector);
        serverComm->receiveFromClientBlocking<double>(name, uNoControlPoints * vNoControlPoints * 4,
                controlPointNet);
        serverComm->receiveFromClientBlocking<int>(name, uNoControlPoints * vNoControlPoints,
                dofIndexNet);

        IGAPatchSurface* thePatch = theIGAMesh->addPatch(pDegree, uNoKnots, uKnotVector, qDegree, vNoKnots, vKnotVector,
                uNoControlPoints, vNoControlPoints, controlPointNet, dofIndexNet);

        // deleting the pointers after object creation
        delete uKnotVector;
        delete vKnotVector;
        delete controlPointNet;
        delete dofIndexNet;

        // Add and linearize the trimming curves
        int trimInfo[2];
        serverComm->receiveFromClientBlocking<int>(name, 2, trimInfo);
        int isTrimmed = trimInfo[0];
        int numLoops = trimInfo[1];
        if(isTrimmed) {
            assert(numLoops>0);
            //Get every loop
            for(int loopCount = 0; loopCount < numLoops; loopCount++) {
                const int BUFFER_SIZE_TRIM = 2;
                int loopInfo[BUFFER_SIZE_TRIM];
                serverComm->receiveFromClientBlocking<int>(name, BUFFER_SIZE_TRIM, loopInfo);

                int inner = loopInfo[0];
                int numCurves = loopInfo[1];
                thePatch->addTrimLoop(inner, numCurves);
                
                const int BUFFER_SIZE_CURVE = 4;
                int curveInfo[BUFFER_SIZE_CURVE];
                //Get every curve
                for(int curveCount = 0; curveCount < numCurves; curveCount++) {
                    serverComm->receiveFromClientBlocking<int>(name, BUFFER_SIZE_CURVE, curveInfo);
                    
                    int direction = curveInfo[0];
                    int pDegree = curveInfo[1];
                    int uNoKnots = curveInfo[2];
                    int uNoControlPoints = curveInfo[3];
                    
                    double* uKnotVector = new double[uNoKnots];
                    double controlPointNet[uNoControlPoints * 4];

                    serverComm->receiveFromClientBlocking<double>(name, uNoKnots, uKnotVector);
                    serverComm->receiveFromClientBlocking<double>(name, uNoControlPoints * 4, controlPointNet);

                    thePatch->addTrimCurve(direction, pDegree, uNoKnots, uKnotVector,
                    		uNoControlPoints, controlPointNet);

                    // deleting the pointers after object creation
                    delete uKnotVector;

                } // end curve
            } // end trimming loops
            thePatch->linearizeTrimming();
        } // end isTrimmed
    } // end patch

    // weak dirichlet condition data
    int numWeakDirichletCond;
    serverComm->receiveFromClientBlocking<int>(name, 1, &numWeakDirichletCond);

    for(int wDCCtr = 0; wDCCtr<numWeakDirichletCond; wDCCtr++){
        const int BUFFER_SIZE_CONNECTION_INFO = 3;
        int connectionInfo[BUFFER_SIZE_CONNECTION_INFO];
        serverComm->receiveFromClientBlocking<int>(name, BUFFER_SIZE_CONNECTION_INFO, connectionInfo);

        int patchIndex = connectionInfo[0];
        int patchBLIndex = connectionInfo[1];
        int patchBLTrCurveIndex = connectionInfo[2];
        int isGPProvided = 0;

        serverComm->receiveFromClientBlocking<int>(name, 1, &isGPProvided);
        if (isGPProvided != 1 && isGPProvided != 0) assert(false);

        WeakIGADirichletCurveCondition* theWeakDirichletCurveCond = theIGAMesh->addWeakDirichletCurveCondition(wDCCtr,
                                                            patchIndex, patchBLIndex, patchBLTrCurveIndex);

        if(isGPProvided){
            int trCurveNumGP;
            serverComm->receiveFromClientBlocking<int>(name, 1, &trCurveNumGP);

            double* trCurveGPs = new double[trCurveNumGP*2];
            double* trCurveGPWeights = new double[trCurveNumGP];
            double* trCurveGPTangents = new double[trCurveNumGP*3];
            double* trCurveGPJacobianProducts = new double[trCurveNumGP];

            serverComm->receiveFromClientBlocking<double>(name,trCurveNumGP*2, trCurveGPs);
            serverComm->receiveFromClientBlocking<double>(name,trCurveNumGP, trCurveGPWeights);
            serverComm->receiveFromClientBlocking<double>(name,trCurveNumGP*3, trCurveGPTangents);
            serverComm->receiveFromClientBlocking<double>(name,trCurveNumGP, trCurveGPJacobianProducts);

            theWeakDirichletCurveCond->addWeakDirichletCurveConditionGPData(trCurveNumGP,
                                               trCurveGPs, trCurveGPWeights,
                                               trCurveGPTangents,
                                               trCurveGPJacobianProducts);

            // deleting the pointers after object creation
            delete trCurveGPs;
            delete trCurveGPWeights;
            delete trCurveGPTangents;
            delete trCurveGPJacobianProducts;

        } else {
            theIGAMesh->createWeakDirichletCurveConditionGPData(wDCCtr);
        }
    }

    // patch coupling data
    int numWeakContCond;
    serverComm->receiveFromClientBlocking<int>(name, 1, &numWeakContCond);

    for(int wPCCCtr = 0; wPCCCtr<numWeakContCond; wPCCCtr++){
        const int BUFFER_SIZE_CONNECTION_INFO = 6;
        int connectionInfo[BUFFER_SIZE_CONNECTION_INFO];
        serverComm->receiveFromClientBlocking<int>(name, BUFFER_SIZE_CONNECTION_INFO, connectionInfo);

        int masterPatchIndex = connectionInfo[0];
        int masterPatchBLIndex = connectionInfo[1];
        int masterPatchBLTrCurveIndex = connectionInfo[2];
        int slavePatchIndex = connectionInfo[3];
        int slavePatchBLIndex = connectionInfo[4];
        int slavePatchBLTrCurveIndex = connectionInfo[5];
        int isGPProvided = 0;

        serverComm->receiveFromClientBlocking<int>(name, 1, &isGPProvided);
        if (isGPProvided != 1 && isGPProvided != 0) assert(false);

        WeakIGAPatchContinuityCondition* theWeakContCond = theIGAMesh->addWeakContinuityCondition(wPCCCtr,
                                                            masterPatchIndex, masterPatchBLIndex, masterPatchBLTrCurveIndex,
                                                            slavePatchIndex,  slavePatchBLIndex,  slavePatchBLTrCurveIndex);

        if(isGPProvided){
            int trCurveNumGP;
            serverComm->receiveFromClientBlocking<int>(name, 1, &trCurveNumGP);

            double* trCurveMasterGPs = new double[trCurveNumGP*2];
            double* trCurveSlaveGPs = new double[trCurveNumGP*2];
            double* trCurveGPWeights = new double[trCurveNumGP];
            double* trCurveMasterGPTangents = new double[trCurveNumGP*3];
            double* trCurveSlaveGPTangents = new double[trCurveNumGP*3];
            double* trCurveGPJacobianProducts = new double[trCurveNumGP];

            serverComm->receiveFromClientBlocking<double>(name,trCurveNumGP*2, trCurveMasterGPs);
            serverComm->receiveFromClientBlocking<double>(name,trCurveNumGP*2, trCurveSlaveGPs);
            serverComm->receiveFromClientBlocking<double>(name,trCurveNumGP, trCurveGPWeights);
            serverComm->receiveFromClientBlocking<double>(name,trCurveNumGP*3, trCurveMasterGPTangents);
            serverComm->receiveFromClientBlocking<double>(name,trCurveNumGP*3, trCurveSlaveGPTangents);
            serverComm->receiveFromClientBlocking<double>(name,trCurveNumGP, trCurveGPJacobianProducts);

            theWeakContCond->addWeakContinuityConditionGPData(trCurveNumGP,
                                               trCurveMasterGPs, trCurveSlaveGPs, trCurveGPWeights,
                                               trCurveMasterGPTangents, trCurveSlaveGPTangents,
                                               trCurveGPJacobianProducts);

            // deleting the pointers after object creation
            delete trCurveMasterGPs;
            delete trCurveSlaveGPs;
            delete trCurveGPWeights;
            delete trCurveMasterGPTangents;
            delete trCurveSlaveGPTangents;
            delete trCurveGPJacobianProducts;

        } else {
            theIGAMesh->createWeakContinuityConditionGPData(wPCCCtr);
        }

    }

//    // get the dirichlet boundary conditions
//    int dirichletBCInfo[2];
//    serverComm->receiveFromClientBlocking<int>(name, 2, dirichletBCInfo);
//    int numberOfClampedDofs = dirichletBCInfo[0];
//    int clampedDirections = dirichletBCInfo[1];

//    int* clampedDofs = new int[numberOfClampedDofs];
//    serverComm->receiveFromClientBlocking<int>(name, numberOfClampedDofs, clampedDofs);

//    theIGAMesh->setClampedDofs(numberOfClampedDofs, clampedDofs);
//    theIGAMesh->setClampedDirections(clampedDirections);

    { // output to shell
        DEBUG_OUT() << (*theIGAMesh) << endl;
        theIGAMesh->computeBoundingBox();
        INFO_OUT() << "\t+Number of Patches: " << theIGAMesh->getSurfacePatches().size() << endl;
        INFO_OUT() << "\t+---------------------------------" << endl << endl;
        INFO_OUT() << theIGAMesh->boundingBox << endl;
    }

    nameToMeshMap.insert(pair<string, AbstractMesh*>(meshName, theIGAMesh));
}

void ClientCode::recvDataField(std::string meshName, std::string dataFieldName) {
    assert(serverComm != NULL);
    assert(nameToMeshMap.find(meshName) != nameToMeshMap.end());
    AbstractMesh *mesh = nameToMeshMap[meshName];
    DataField *df = mesh->getDataFieldByName(dataFieldName);

    { // output to shell
        string info = "Emperor is receiving (" + meshName + ": " + dataFieldName + ") from [" + name
                + "] ...";
        INDENT_OUT(1, info, infoOut);
    }

    int sizeClientSay = -1;
    serverComm->receiveFromClientBlocking<int>(name, 1, &sizeClientSay);
    int bufferSize = df->numLocations * df->dimension;
    assert(bufferSize == sizeClientSay);

    serverComm->receiveFromClientBlocking<double>(name, bufferSize, df->data);
    DEBUG_OUT() << (*df) << endl;
}

void ClientCode::sendDataField(std::string meshName, std::string dataFieldName) {
    assert(serverComm != NULL);
    assert(nameToMeshMap.find(meshName) != nameToMeshMap.end());
    AbstractMesh *mesh = nameToMeshMap[meshName];
    DataField *df = mesh->getDataFieldByName(dataFieldName);

    { // output to shell
        string info = "Emperor is sending (" + meshName + ": " + dataFieldName + ") to [" + name
                + "] ...";
        INDENT_OUT(1, info, infoOut);
    }

    int bufferSize = df->numLocations * df->dimension;
    serverComm->sendToClientBlocking<int>(name, 1, &bufferSize);
    serverComm->sendToClientBlocking<double>(name, bufferSize, df->data);
    // debug
    // (*df).writeToFile("sentData.dat","sent data from empire", "end");
    DEBUG_OUT() << (*df) << endl;
}

void ClientCode::sendConvergenceSignal(bool convergent) {
    assert(serverComm != NULL);
    int tmpInt = (int) convergent;

    { // output to shell
        string converg;
        if (convergent)
            converg = "true";
        else
            converg = "false";
        string info = "Emperor is sending convergence signal (" + converg + ") to [" + name
                + "] ...";
        INDENT_OUT(1, info, infoOut);
    }

    serverComm->sendToClientBlocking<int>(name, 1, &tmpInt);
}

bool ClientCode::recvConvergenceSignal() {
    assert(serverComm != NULL);
    int tmpInt;
    serverComm->receiveFromClientBlocking<int>(name, 1, &tmpInt);
    bool convergent = (bool) tmpInt;
    { // output to shell
        string converg;
        if (convergent)
            converg = "true";
        else
            converg = "false";
        string info = "Emperor is receiving convergence signal (" + converg + ") from [" + name
                + "] ...";
        INDENT_OUT(1, info, infoOut);
    }
    return convergent;
}

AbstractMesh *ClientCode::getMeshByName(std::string meshName) {
    assert(nameToMeshMap.find(meshName) != nameToMeshMap.end());
    return nameToMeshMap[meshName];
}

void ClientCode::addSignal(std::string signalName, int size0, int size1, int size2) {

    if (nameToSignalMap.find(signalName) != nameToSignalMap.end()) {
        ERROR_OUT() << "Signal name: " << signalName << " already used!" << endl;
        assert(false);
    }
    nameToSignalMap.insert(
            pair<string, Signal*>(signalName, new Signal(signalName, size0, size1, size2)));
}

void ClientCode::recvSignal(std::string signalName) {
    assert(nameToSignalMap.find(signalName) != nameToSignalMap.end());
    { // output to shell
        string info = "Emperor is receiving (" + signalName + ") from [" + name + "] ...";
        INDENT_OUT(1, info, infoOut);
    }
    Signal *signal = nameToSignalMap[signalName];
    int signalSizeReceive = -1;
    const int NAME_STRING_LENGTH = ServerCommunication::NAME_STRING_LENGTH;
    char signalNameReceive[NAME_STRING_LENGTH];
    serverComm->receiveFromClientBlocking<char>(name, NAME_STRING_LENGTH, signalNameReceive);
    string tmp(signalNameReceive);
    { // output to shell
        string info = "Signal name received is \"" + tmp + "\"";
        INDENT_OUT(2, info, infoOut);
    }
    if (signalName.compare(tmp) != 0) {
        ERROR_OUT() << "Signal name received is " << tmp << ", however " << signalName
                << " expected!" << endl;
        assert(false);
    }

    serverComm->receiveFromClientBlocking<int>(name, 1, &signalSizeReceive);
    if (signalSizeReceive != signal->size) {
        ERROR_OUT() << "Signal " << signalSizeReceive << " size received is " << signalSizeReceive
                << ", however " << signal->size << " expected!" << endl;
        assert(false);
    }
    serverComm->receiveFromClientBlocking<double>(name, signal->size, signal->array);
    DEBUG_OUT() << (*signal) << endl;
}

void ClientCode::sendSignal(std::string signalName) {
    assert(nameToSignalMap.find(signalName) != nameToSignalMap.end());
    { // output to shell
        string info = "Emperor is sending (" + signalName + ") to [" + name + "] ...";
        INDENT_OUT(1, info, infoOut);
    }
    const Signal *signal = nameToSignalMap[signalName];
    const int NAME_STRING_LENGTH = ServerCommunication::NAME_STRING_LENGTH;
    char signalNameSend[NAME_STRING_LENGTH];
    strcpy(signalNameSend, signalName.c_str());
    serverComm->sendToClientBlocking<char>(name, NAME_STRING_LENGTH, signalNameSend);
    assert(signal->size != 0);
    int signalSize = signal->size;
    serverComm->sendToClientBlocking<int>(name, 1, &(signalSize));
    serverComm->sendToClientBlocking<double>(name, signal->size, signal->array);
    DEBUG_OUT() << *(const_cast<Signal*>(signal)) << endl;
}

Signal *ClientCode::getSignalByName(std::string signalName) {
    if (nameToSignalMap.find(signalName) == nameToSignalMap.end()) {
        ERROR_OUT("Signal name: "+signalName+" not found!");
        for (map<std::string, Signal*>::iterator it = nameToSignalMap.begin();
                it != nameToSignalMap.end(); it++) {
            ERROR_OUT("Possible signal names for client "+name+" are : "+it->first);
        }

        assert(nameToSignalMap.find(signalName) != nameToSignalMap.end());
    }
    return nameToSignalMap[signalName];
}

} /* namespace EMPIRE */

