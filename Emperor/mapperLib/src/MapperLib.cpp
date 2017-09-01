/*  Copyright &copy; 2013, TU Muenchen, Chair of Structural Analysis,
 *  Stefan Sicklinger, Tianyang Wang, Andreas Apostolatos, Altug Emiroglu Munich
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
 * \file MapperLib.cpp
 * This file wraps the C++ API in a C API
 * \date 5/12/2014
 **************************************************************************************************/
#include <assert.h>
#include <string>

#include "AbstractMesh.h"
#include "FEMesh.h"
#include "IGAMesh.h"
#include "IGAPatchSurface.h"
#include "WeakIGADirichletCurveCondition.h"
#include "WeakIGAPatchContinuityCondition.h"

#include "AbstractMapper.h"
#include "NearestNeighborMapper.h"
#include "NearestElementMapper.h"
#include "BarycentricInterpolationMapper.h"
#include "MortarMapper.h"
#include "IGAMortarMapper.h"
#include "MapperLib.h"

using namespace EMPIRE;
using namespace std;

/// AbstractMapper object map in the global scope
std::map <std::string, AbstractMapper*> mapperList;
/// AbstractMesh object map in the global scope
std::map <std::string, AbstractMesh*> meshList;

void init_FE_NearestNeighborMapper(char* mapperName, 
                                   int AnumNodes, const double *Anodes, int BnumNodes, const double *Bnodes){

    WARNING_OUT("In \"init_FE_NearestNeighborMapper\": This interface function of the mapper library is going to be removed!");
    WARNING_OUT("Please check \"initFEMNearestNeighborMapper\" for the new declaration.");

    std::string mapperNameToMap = std::string(mapperName);

    if (mapperList.count( mapperNameToMap )){
        ERROR_OUT("A mapper with name : \"" + mapperNameToMap + "\" has already been initialized!");
        exit(EXIT_FAILURE);
    }
    else{
        mapperList[mapperNameToMap] = new NearestNeighborMapper(AnumNodes, Anodes, BnumNodes, Bnodes);
        mapperList[mapperNameToMap]->buildCouplingMatrices();
        INFO_OUT("Mapper \"" + mapperNameToMap + "\" is generated");
    }
}

void init_FE_NearestElementMapper(char* mapperName,
                                  int AnumNodes, int AnumElems, const int *AnumNodesPerElem, const double *Anodes, const int *AnodeIDs, const int *Aelems,
                                  int BnumNodes, int BnumElems, const int *BnumNodesPerElem, const double *Bnodes, const int *BnodeIDs, const int *Belems){
    
    WARNING_OUT("In \"init_FE_NearestElementMapper\": This interface function of the mapper library is going to be removed!");
    WARNING_OUT("Please check \"initFEMNearestElementMapper\" for the new declaration.");

    std::string mapperNameToMap = std::string(mapperName);

    if (mapperList.count( mapperNameToMap )){
        ERROR_OUT("A mapper with name : \"" + mapperNameToMap + "\" has already been initialized!");
        exit(EXIT_FAILURE);
    }
    else{
        mapperList[mapperNameToMap] = new NearestElementMapper(AnumNodes, AnumElems, AnumNodesPerElem, Anodes, AnodeIDs, Aelems,
                                                               BnumNodes, BnumElems, BnumNodesPerElem, Bnodes, BnodeIDs, Belems);
        mapperList[mapperNameToMap]->buildCouplingMatrices();
        INFO_OUT("Mapper \"" + mapperNameToMap + "\" is generated");
    }
}

void init_FE_BarycentricInterpolationMapper(char* mapperName, 
                                            int AnumNodes, const double *Anodes, int BnumNodes, const double *Bnodes){

    WARNING_OUT("In \"init_FE_BarycentricInterpolationMapper\": This interface function of the mapper library is going to be removed!");
    WARNING_OUT("Please check \"initFEMBarycentricInterpolationMapper\" for the new declaration.");

    std::string mapperNameToMap = std::string(mapperName);

    if (mapperList.count( mapperNameToMap )){
        ERROR_OUT("A mapper with name : " + mapperNameToMap + "has already been initialized!");
        exit(EXIT_FAILURE);
    } else {
        mapperList[mapperNameToMap] = new BarycentricInterpolationMapper(AnumNodes, Anodes, BnumNodes, Bnodes);
        mapperList[mapperNameToMap]->buildCouplingMatrices();
        INFO_OUT("Mapper \"" + mapperNameToMap + "\" is generated");
    }
}

void init_FE_MortarMapper(char* mapperName,
                          int AnumNodes, int AnumElems, const int* AnumNodesPerElem, const double* Anodes, const int* AnodeIDs, const int* Aelems,
                          int BnumNodes, int BnumElems, const int* BnumNodesPerElem, const double* Bnodes, const int* BnodeIDs, const int* Belems,
                          int oppositeSurfaceNormal, int dual, int enforceConsistency){

    WARNING_OUT("In \"init_FE_MortarMapper\": This interface function of the mapper library is going to be removed!");
    WARNING_OUT("Please check \"initFEMMortarMapper\" for the new declaration.");

    bool _oppositeSurfaceNormal = false;
    bool _dual = false;
    bool _enforceConsistency = false;

    if (oppositeSurfaceNormal != 0) _oppositeSurfaceNormal = true;
    if (dual != 0) _dual = true;
    if (enforceConsistency != 0) _enforceConsistency = true;

    std::string mapperNameToMap = std::string(mapperName);

    if (mapperList.count( mapperNameToMap )){
        ERROR_OUT("A mapper with name : " + mapperNameToMap + "has already been initialized!");
        exit(EXIT_FAILURE);
    } else {
        mapperList[mapperNameToMap] = new MortarMapper(AnumNodes, AnumElems, AnumNodesPerElem, Anodes, AnodeIDs, Aelems,
                                                       BnumNodes, BnumElems, BnumNodesPerElem, Bnodes, BnodeIDs, Belems,
                                                       _oppositeSurfaceNormal, _dual, _enforceConsistency);
        mapperList[mapperNameToMap]->buildCouplingMatrices();
        INFO_OUT("Mapper \"" + mapperNameToMap + "\" is generated");
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////

void initFEMesh(char* meshName, int numNodes, int numElems, bool triangulateAll = false){

    std::string meshNameToMap = std::string(meshName);

    if (meshList.count( meshNameToMap )){
        ERROR_OUT("A mesh with name : " + meshNameToMap + " has already been initialized!");
        ERROR_OUT("Mesh not generated!");
        return;
    } else {
        FEMesh* tmpFEMesh = new FEMesh(meshNameToMap, numNodes, numElems, triangulateAll);
        meshList[meshNameToMap] = tmpFEMesh;
        INFO_OUT("FEMesh \"" + meshNameToMap +  "\" generated");
    }
}

void setNodesToFEMesh(char* meshName, int* nodeIDs, double* nodes){

    std::string meshNameInMap = std::string(meshName);

    if (!meshList.count( meshNameInMap )){
        ERROR_OUT("A mesh with name : " + meshNameInMap + " does not exist!");
        ERROR_OUT("Did nothing!");
        return;
    } else if (meshList[meshNameInMap]->type == EMPIRE_Mesh_FEMesh){
        FEMesh *tmpFEMesh = dynamic_cast<FEMesh *>(meshList[meshNameInMap]);
        for(int i=0; i<tmpFEMesh->numNodes;i++) tmpFEMesh->nodeIDs[i] = nodeIDs[i];
        for(int i=0; i<(tmpFEMesh->numNodes)*3;i++) tmpFEMesh->nodes[i] = nodes[i];
        INFO_OUT("Nodes set to \"" + meshNameInMap );
    }
}

void setElementsToFEMesh(char* meshName, int* numNodesPerElem, int* elems){

    std::string meshNameInMap = std::string(meshName);

    if (!meshList.count( meshNameInMap )){
        ERROR_OUT("A mesh with name : " + meshNameInMap + " does not exist!");
        ERROR_OUT("Did nothing!");
        return;
    } else if (meshList[meshNameInMap]->type == EMPIRE_Mesh_FEMesh){
        FEMesh *tmpFEMesh = dynamic_cast<FEMesh *>(meshList[meshNameInMap]);
        tmpFEMesh->elems = NULL;
        for(int i=0; i<tmpFEMesh->numElems;i++) tmpFEMesh->numNodesPerElem[i] = numNodesPerElem[i];
        tmpFEMesh->initElems();
        for(int i=0; i<tmpFEMesh->elemsArraySize;i++) tmpFEMesh->elems[i] = elems[i];
        INFO_OUT("Elements set to \"" + meshNameInMap );
    }
}

void initIGAMesh(char* meshName){

    std::string meshNameToMap = std::string(meshName);

    if (meshList.count( meshNameToMap )){
        ERROR_OUT("A mesh with name : " + meshNameToMap + " has already been initialized!");
        ERROR_OUT("Mesh not generated!");
        return;
    } else {
        meshList[meshNameToMap] = new IGAMesh(meshNameToMap);
        INFO_OUT("IGAMesh \"" + meshNameToMap +  "\" generated");
    }
}

void addPatchToIGAMesh(char* meshName,
                       int pDegree, int uNoKnots, double* uKnotVector,
                       int qDegree, int vNoKnots, double* vKnotVector,
                       int uNoControlPoints, int vNoControlPoints,
                       double* controlPointNet, int* dofIndexNet){

    std::string meshNameInMap = std::string(meshName);

    if (!meshList.count( meshNameInMap )){
        ERROR_OUT("A mesh with name : " + meshNameInMap + " does not exist!");
        ERROR_OUT("Did nothing!");
        return;
    } else if (meshList[meshNameInMap]->type == EMPIRE_Mesh_IGAMesh){
        IGAMesh *tmpIGAMesh = dynamic_cast<IGAMesh *>(meshList[meshNameInMap]);
        tmpIGAMesh->addPatch(pDegree, uNoKnots, uKnotVector,
                             qDegree, vNoKnots, vKnotVector,
                             uNoControlPoints, vNoControlPoints,
                             controlPointNet, dofIndexNet);
        INFO_OUT("Added patch to \"" + meshNameInMap );
    }
}

void addTrimmingLoopToPatch(char* meshName, int patchIndex,
                            int inner, int numCurves){

    std::string meshNameInMap = std::string(meshName);

    if (!meshList.count( meshNameInMap )){
        ERROR_OUT("A mesh with name : " + meshNameInMap + "does not exist!");
        ERROR_OUT("Did nothing!");
        return;
    } else if (meshList[meshNameInMap]->type == EMPIRE_Mesh_IGAMesh){
        IGAMesh *tmpIGAMesh = dynamic_cast<IGAMesh *>(meshList[meshNameInMap]);
        tmpIGAMesh->getSurfacePatch(patchIndex)->addTrimLoop(inner,numCurves);
        INFO_OUT()<<"Added trimming loop to \"" << meshNameInMap << "\" patch index " << patchIndex <<  std::endl;
    }
}

void addTrimmingCurveToTrimmingLoop(char* meshName, int patchIndex,
                                    int direction, int pDegree, int uNoKnots, double* uKnotVector,
                                    int uNoControlPoints, double* controlPoints){

    std::string meshNameInMap = std::string(meshName);

    if (!meshList.count( meshNameInMap )){
        ERROR_OUT("A mesh with name : " + meshNameInMap + "does not exist!");
        ERROR_OUT("Did nothing!");
        return;
    } else if (meshList[meshNameInMap]->type == EMPIRE_Mesh_IGAMesh){
        IGAMesh *tmpIGAMesh = dynamic_cast<IGAMesh *>(meshList[meshNameInMap]);
        tmpIGAMesh->getSurfacePatch(patchIndex)->addTrimCurve(direction, pDegree, uNoKnots, uKnotVector,
                                                              uNoControlPoints, controlPoints);
        INFO_OUT()<<"Added trimming curve to \"" << meshNameInMap << "\" patch index " << patchIndex <<  std::endl;
    }
}

void linearizeTrimmingLoops(char* meshName, int patchIndex){

    std::string meshNameInMap = std::string(meshName);

    if (!meshList.count( meshNameInMap )){
        ERROR_OUT("A mesh with name : " + meshNameInMap + " does not exist!");
        ERROR_OUT("Did nothing!");
        return;
    } else if (meshList[meshNameInMap]->type == EMPIRE_Mesh_IGAMesh){
        IGAMesh *tmpIGAMesh = dynamic_cast<IGAMesh *>(meshList[meshNameInMap]);
        tmpIGAMesh->getSurfacePatch(patchIndex)->linearizeTrimming();
        INFO_OUT()<<"Linearized trimming curves of \"" << meshNameInMap << "\" patch index " << patchIndex <<  std::endl;
    }
}

void addDirichletCurveConditionToIGAMesh(char* meshName,
                                         int conditionID, int patchIndex,
                                         int pDegree, int uNoKnots, double* uKnotVector,
                                         int uNoControlPoints, double* controlPointNet) {

    std::string meshNameInMap = std::string(meshName);

    if (!meshList.count( meshNameInMap )){
        ERROR_OUT("A mesh with name : " + meshNameInMap + " does not exist!");
        ERROR_OUT("Did nothing!");
        return;
    } else if (meshList[meshNameInMap]->type == EMPIRE_Mesh_IGAMesh){
        IGAMesh *tmpIGAMesh = dynamic_cast<IGAMesh *>(meshList[meshNameInMap]);
        WeakIGADirichletCurveCondition* tmpCond = tmpIGAMesh->addWeakDirichletCurveCondition(conditionID, patchIndex, pDegree, uNoKnots, uKnotVector, uNoControlPoints, controlPointNet);
        tmpCond->createGPData(tmpIGAMesh->getSurfacePatch(patchIndex));
        INFO_OUT()<<"Added a Dirichlet curve condition to \"" << meshNameInMap << "\" patch index " << patchIndex <<  std::endl;
    }
}

void addDirichletBoundaryConditionToIGAMesh(char* meshName,
                                            int conditionID,
                                            int patchIndex, int patchBLIndex, int patchBLTrCurveIndex) {

    std::string meshNameInMap = std::string(meshName);

    if (!meshList.count( meshNameInMap )){
        ERROR_OUT("A mesh with name : " + meshNameInMap + " does not exist!");
        ERROR_OUT("Did nothing!");
        return;
    } else if (meshList[meshNameInMap]->type == EMPIRE_Mesh_IGAMesh){
        IGAMesh *tmpIGAMesh = dynamic_cast<IGAMesh *>(meshList[meshNameInMap]);
        WeakIGADirichletCurveCondition* tmpCond = tmpIGAMesh->addWeakDirichletCurveCondition(conditionID, patchIndex, patchBLIndex, patchBLTrCurveIndex);
        tmpCond->createGPData(tmpIGAMesh->getSurfacePatch(patchIndex));
        INFO_OUT()<<"Added a Dirichlet boundary condition to \"" << meshNameInMap << "\" patch index " << patchIndex <<  std::endl;
    }
}

void addPatchContinuityConditionToIGAMesh(char* meshName,
                                          int connectionID,
                                          int masterPatchIndex, int masterPatchBLIndex, int masterPatchBLTrCurveIndex,
                                          int slavePatchIndex,  int slavePatchBLIndex,  int slavePatchBLTrCurveIndex) {

    std::string meshNameInMap = std::string(meshName);

    if (!meshList.count( meshNameInMap )){
        ERROR_OUT("A mesh with name : " + meshNameInMap + " does not exist!");
        ERROR_OUT("Did nothing!");
        return;
    } else if (meshList[meshNameInMap]->type == EMPIRE_Mesh_IGAMesh){
        IGAMesh *tmpIGAMesh = dynamic_cast<IGAMesh *>(meshList[meshNameInMap]);
        WeakIGAPatchContinuityCondition* tmpContCond = tmpIGAMesh->addWeakContinuityCondition(connectionID,
                                                                                              masterPatchIndex, masterPatchBLIndex, masterPatchBLTrCurveIndex,
                                                                                              slavePatchIndex, slavePatchBLIndex, slavePatchBLTrCurveIndex);
        tmpContCond->createGPData(tmpIGAMesh->getSurfacePatch(masterPatchIndex), tmpIGAMesh->getSurfacePatch(slavePatchIndex));
        INFO_OUT()<<"Added a continuity condition to \"" << meshNameInMap << "\" between patch indices " << masterPatchIndex << " and " << slavePatchIndex << std::endl;
    }
}

void addPatchContinuityConditionOnCurvesToIGAMesh(char* meshName,
                                                  int connectionID,
                                                  int masterPatchIndex, int pMaster, int uNoKnotsMaster, double* uKnotVectorMaster, int uNoControlPointsMaster, double* controlPointNetMaster,
                                                  int slavePatchIndex,  int pSlave, int uNoKnotsSlave, double* uKnotVectorSlave, int uNoControlPointsSlave, double* controlPointNetSlave) {

    std::string meshNameInMap = std::string(meshName);

    if (!meshList.count( meshNameInMap )){
        ERROR_OUT("A mesh with name : " + meshNameInMap + " does not exist!");
        ERROR_OUT("Did nothing!");
        return;
    } else if (meshList[meshNameInMap]->type == EMPIRE_Mesh_IGAMesh){
        IGAMesh *tmpIGAMesh = dynamic_cast<IGAMesh *>(meshList[meshNameInMap]);
        WeakIGAPatchContinuityCondition* tmpContCond = tmpIGAMesh->addWeakContinuityCondition(connectionID,
                                                                                              masterPatchIndex, pMaster, uNoKnotsMaster, uKnotVectorMaster, uNoControlPointsMaster, controlPointNetMaster,
                                                                                              slavePatchIndex, pSlave, uNoKnotsSlave, uKnotVectorSlave, uNoControlPointsSlave, controlPointNetSlave);
        tmpContCond->createGPData(tmpIGAMesh->getSurfacePatch(masterPatchIndex), tmpIGAMesh->getSurfacePatch(slavePatchIndex));
        INFO_OUT()<<"Added a continuity condition to \"" << meshNameInMap << "\" between patch indices " << masterPatchIndex << " and " << slavePatchIndex << std::endl;
    }
}


void initFEMMortarMapper(char* mapperName, char* AmeshName, char* BmeshName,
                         int oppositeSurfaceNormal, int dual, int enforceConsistency){

    std::string mapperNameToMap = std::string(mapperName);
    std::string aFEMeshNameInMap = std::string(AmeshName);
    std::string bFEMeshNameInMap = std::string(BmeshName);

    FEMesh *tmpaFEMesh;
    FEMesh *tmpbFEMesh;

    // check if the mesh with the given name is generated and is of correct type
    if (!meshList.count( aFEMeshNameInMap )){
        ERROR_OUT("A mesh with name : " + aFEMeshNameInMap + " does not exist!");
        ERROR_OUT("Mapper not generated!");
        return;
    } else if (meshList[aFEMeshNameInMap]->type != EMPIRE_Mesh_FEMesh){
        ERROR_OUT(aFEMeshNameInMap + " is not a type of FEMesh");
        ERROR_OUT("Mapper not generated!");
        return;
    } else {
        tmpaFEMesh = dynamic_cast<FEMesh *>(meshList[bFEMeshNameInMap]);
    }

    // check if the mesh with the given name is generated and is of correct type
    if (!meshList.count( bFEMeshNameInMap )){
        ERROR_OUT("A mesh with name : " + bFEMeshNameInMap + " does not exist!");
        ERROR_OUT("Mapper not generated!");
        return;
    } else if (meshList[bFEMeshNameInMap]->type != EMPIRE_Mesh_FEMesh){
        ERROR_OUT(bFEMeshNameInMap + " is not a type of FEMesh");
        ERROR_OUT("Mapper not generated!");
        return;
    } else {
        tmpbFEMesh = dynamic_cast<FEMesh *>(meshList[bFEMeshNameInMap]);
    }

    bool _oppositeSurfaceNormal = false;
    bool _dual = false;
    bool _enforceConsistency = false;

    if (oppositeSurfaceNormal != 0) _oppositeSurfaceNormal = true;
    if (dual != 0) _dual = true;
    if (enforceConsistency != 0) _enforceConsistency = true;

    if (mapperList.count( mapperNameToMap )){
        ERROR_OUT("A mapper with name : " + mapperNameToMap + " has already been initialized!");
        ERROR_OUT("Mapper not generated!");
        return;
    } else {

        if (!tmpaFEMesh->boundingBox.isComputed()) tmpaFEMesh->computeBoundingBox();
        if (!tmpbFEMesh->boundingBox.isComputed()) tmpbFEMesh->computeBoundingBox();

        mapperList[mapperNameToMap] = new MortarMapper(tmpaFEMesh->numNodes, tmpaFEMesh->numElems, tmpaFEMesh->numNodesPerElem, tmpaFEMesh->nodes, tmpaFEMesh->nodeIDs, tmpaFEMesh->elems,
                                                       tmpbFEMesh->numNodes, tmpbFEMesh->numElems, tmpbFEMesh->numNodesPerElem, tmpbFEMesh->nodes, tmpbFEMesh->nodeIDs, tmpbFEMesh->elems,
                                                       _oppositeSurfaceNormal, _dual, _enforceConsistency);
        INFO_OUT("Generated \"" +  mapperNameToMap);
    }

}

void initFEMNearestNeighborMapper(char* mapperName, char* AmeshName, char* BmeshName){

    std::string mapperNameToMap = std::string(mapperName);
    std::string aFEMeshNameInMap = std::string(AmeshName);
    std::string bFEMeshNameInMap = std::string(BmeshName);

    FEMesh *tmpaFEMesh;
    FEMesh *tmpbFEMesh;

    // check if the mesh with the given name is generated and is of correct type
    if (!meshList.count( aFEMeshNameInMap )){
        ERROR_OUT("A mesh with name : " + aFEMeshNameInMap + " does not exist!");
        ERROR_OUT("Mapper not generated!");
        return;
    } else if (meshList[aFEMeshNameInMap]->type != EMPIRE_Mesh_FEMesh){
        ERROR_OUT(aFEMeshNameInMap + " is not a type of FEMesh");
        ERROR_OUT("Mapper not generated!");
        return;
    } else {
        tmpaFEMesh = dynamic_cast<FEMesh *>(meshList[bFEMeshNameInMap]);
    }

    // check if the mesh with the given name is generated and is of correct type
    if (!meshList.count( bFEMeshNameInMap )){
        ERROR_OUT("A mesh with name : " + bFEMeshNameInMap + " does not exist!");
        ERROR_OUT("Mapper not generated!");
        return;
    } else if (meshList[bFEMeshNameInMap]->type != EMPIRE_Mesh_FEMesh){
        ERROR_OUT(bFEMeshNameInMap + " is not a type of FEMesh");
        ERROR_OUT("Mapper not generated!");
        return;
    } else {
        tmpbFEMesh = dynamic_cast<FEMesh *>(meshList[bFEMeshNameInMap]);
    }

    if (mapperList.count( mapperNameToMap )){
        ERROR_OUT("A mapper with name : " + mapperNameToMap + " has already been initialized!");
        ERROR_OUT("Mapper not generated!");
        return;
    } else {

        if (!tmpaFEMesh->boundingBox.isComputed()) tmpaFEMesh->computeBoundingBox();
        if (!tmpbFEMesh->boundingBox.isComputed()) tmpbFEMesh->computeBoundingBox();

        mapperList[mapperNameToMap] = new NearestNeighborMapper(tmpaFEMesh->numNodes, tmpaFEMesh->nodes,
                                                                tmpbFEMesh->numNodes, tmpbFEMesh->nodes);
        INFO_OUT("Generated \"" +  mapperNameToMap);
    }

}

void initFEMNearestElementMapper(char* mapperName, char* AmeshName, char* BmeshName){

    std::string mapperNameToMap = std::string(mapperName);
    std::string aFEMeshNameInMap = std::string(AmeshName);
    std::string bFEMeshNameInMap = std::string(BmeshName);

    FEMesh *tmpaFEMesh;
    FEMesh *tmpbFEMesh;

    // check if the mesh with the given name is generated and is of correct type
    if (!meshList.count( aFEMeshNameInMap )){
        ERROR_OUT("A mesh with name : " + aFEMeshNameInMap + " does not exist!");
        ERROR_OUT("Mapper not generated!");
        return;
    } else if (meshList[aFEMeshNameInMap]->type != EMPIRE_Mesh_FEMesh){
        ERROR_OUT(aFEMeshNameInMap + " is not a type of FEMesh");
        ERROR_OUT("Mapper not generated!");
        return;
    } else {
        tmpaFEMesh = dynamic_cast<FEMesh *>(meshList[bFEMeshNameInMap]);
    }

    // check if the mesh with the given name is generated and is of correct type
    if (!meshList.count( bFEMeshNameInMap )){
        ERROR_OUT("A mesh with name : " + bFEMeshNameInMap + " does not exist!");
        ERROR_OUT("Mapper not generated!");
        return;
    } else if (meshList[bFEMeshNameInMap]->type != EMPIRE_Mesh_FEMesh){
        ERROR_OUT(bFEMeshNameInMap + " is not a type of FEMesh");
        ERROR_OUT("Mapper not generated!");
        return;
    } else {
        tmpbFEMesh = dynamic_cast<FEMesh *>(meshList[bFEMeshNameInMap]);
    }

    if (mapperList.count( mapperNameToMap )){
        ERROR_OUT("A mapper with name : " + mapperNameToMap + " has already been initialized!");
        ERROR_OUT("Mapper not generated!");
        return;
    } else {

        if (!tmpaFEMesh->boundingBox.isComputed()) tmpaFEMesh->computeBoundingBox();
        if (!tmpbFEMesh->boundingBox.isComputed()) tmpbFEMesh->computeBoundingBox();

        mapperList[mapperNameToMap] = new NearestElementMapper(tmpaFEMesh->numNodes, tmpaFEMesh->numElems, tmpaFEMesh->numNodesPerElem, tmpaFEMesh->nodes, tmpaFEMesh->nodeIDs, tmpaFEMesh->elems,
                                                               tmpbFEMesh->numNodes, tmpbFEMesh->numElems, tmpbFEMesh->numNodesPerElem, tmpbFEMesh->nodes, tmpbFEMesh->nodeIDs, tmpbFEMesh->elems);
        INFO_OUT("Generated \"" +  mapperNameToMap);
    }
}

void initFEMBarycentricInterpolationMapper(char* mapperName, char* AmeshName, char* BmeshName){

    std::string mapperNameToMap = std::string(mapperName);
    std::string aFEMeshNameInMap = std::string(AmeshName);
    std::string bFEMeshNameInMap = std::string(BmeshName);

    FEMesh *tmpaFEMesh;
    FEMesh *tmpbFEMesh;

    // check if the mesh with the given name is generated and is of correct type
    if (!meshList.count( aFEMeshNameInMap )){
        ERROR_OUT("A mesh with name : " + aFEMeshNameInMap + " does not exist!");
        ERROR_OUT("Mapper not generated!");
        return;
    } else if (meshList[aFEMeshNameInMap]->type != EMPIRE_Mesh_FEMesh){
        ERROR_OUT(aFEMeshNameInMap + " is not a type of FEMesh");
        ERROR_OUT("Mapper not generated!");
        return;
    } else {
        tmpaFEMesh = dynamic_cast<FEMesh *>(meshList[bFEMeshNameInMap]);
    }

    // check if the mesh with the given name is generated and is of correct type
    if (!meshList.count( bFEMeshNameInMap )){
        ERROR_OUT("A mesh with name : " + bFEMeshNameInMap + " does not exist!");
        ERROR_OUT("Mapper not generated!");
        return;
    } else if (meshList[bFEMeshNameInMap]->type != EMPIRE_Mesh_FEMesh){
        ERROR_OUT(bFEMeshNameInMap + " is not a type of FEMesh");
        ERROR_OUT("Mapper not generated!");
        return;
    } else {
        tmpbFEMesh = dynamic_cast<FEMesh *>(meshList[bFEMeshNameInMap]);
    }

    if (mapperList.count( mapperNameToMap )){
        ERROR_OUT("A mapper with name : " + mapperNameToMap + " has already been initialized!");
        ERROR_OUT("Mapper not generated!");
        return;
    } else {

        if (!tmpaFEMesh->boundingBox.isComputed()) tmpaFEMesh->computeBoundingBox();
        if (!tmpbFEMesh->boundingBox.isComputed()) tmpbFEMesh->computeBoundingBox();

        mapperList[mapperNameToMap] = new BarycentricInterpolationMapper(tmpaFEMesh->numNodes, tmpaFEMesh->nodes,
                                                                         tmpbFEMesh->numNodes, tmpbFEMesh->nodes);
        INFO_OUT("Generated \"" +  mapperNameToMap);
    }
}

void initIGAMortarMapper(char* mapperName, char* IGAMeshName, char* FEMeshName, bool isWeakConditions, bool isMappingIGA2FEM){

    std::string mapperNameToMap = std::string(mapperName);
    std::string IGAMeshNameInMap = std::string(IGAMeshName);
    std::string FEMeshNameInMap = std::string(FEMeshName);

    IGAMesh *tmpIGAMesh;
    FEMesh *tmpFEMesh;

    // check if the mesh with the given name is generated and is of correct type
    if (!meshList.count( IGAMeshNameInMap )){
        ERROR_OUT("A mesh with name : " + IGAMeshNameInMap + " does not exist!");
        ERROR_OUT("Mapper not generated!");
        return;
    } else if (meshList[IGAMeshNameInMap]->type != EMPIRE_Mesh_IGAMesh){
        ERROR_OUT(IGAMeshNameInMap + " is not a type of IGAMesh");
        ERROR_OUT("Mapper not generated!");
        return;
    } else {
        tmpIGAMesh = dynamic_cast<IGAMesh *>(meshList[IGAMeshNameInMap]);
    }

    // check if the mesh with the given name is generated and is of correct type
    if (!meshList.count( FEMeshNameInMap )){
        ERROR_OUT("A mesh with name : " + FEMeshNameInMap + " does not exist!");
        ERROR_OUT("Mapper not generated!");
        return;
    } else if (meshList[FEMeshNameInMap]->type != EMPIRE_Mesh_FEMesh){
        ERROR_OUT(FEMeshNameInMap + " is not a type of FEMesh");
        ERROR_OUT("Mapper not generated!");
        return;
    } else {
        tmpFEMesh = dynamic_cast<FEMesh *>(meshList[FEMeshNameInMap]);
    }

    if (mapperList.count( mapperNameToMap )){
        ERROR_OUT("A mapper with name : " + mapperNameToMap + " has already been initialized!");
        ERROR_OUT("Mapper not generated!");
        return;
    } else {

        if (!tmpIGAMesh->boundingBox.isComputed()) tmpIGAMesh->computeBoundingBox();
        if (!tmpFEMesh->boundingBox.isComputed()) tmpFEMesh->computeBoundingBox();

        mapperList[mapperNameToMap] = new IGAMortarMapper(mapperNameToMap, tmpIGAMesh, tmpFEMesh, isWeakConditions, isMappingIGA2FEM);
        INFO_OUT("Generated \"" +  mapperNameToMap + "\"");
    }
}

void setParametersProjection(char* mapperName,
                             double maxProjectionDistance, int numRefinementForIntialGuess,
                             double maxDistanceForProjectedPointsOnDifferentPatches){

    std::string mapperNameInMap = std::string(mapperName);

    IGAMortarMapper *tmpIGAMortarMapper;

    // check if the mapper with the given name is generated and is of correct type
    if (!mapperList.count( mapperNameInMap )){
        ERROR_OUT("A mapper with name : " + mapperNameInMap + " does not exist!");
        ERROR_OUT("Did nothing!");
        return;
    } else if (mapperList[mapperNameInMap]->mapperType != EMPIRE_IGAMortarMapper){
        ERROR_OUT(mapperNameInMap + " is not a type of IGAMortarMapper");
        ERROR_OUT("Did nothing!");
        return;
    }
    else{
        tmpIGAMortarMapper = dynamic_cast<IGAMortarMapper *>(mapperList[mapperNameInMap]);
        tmpIGAMortarMapper->setParametersProjection(maxProjectionDistance, numRefinementForIntialGuess, maxDistanceForProjectedPointsOnDifferentPatches);
        INFO_OUT("Point projection parameters are set for \"" +  mapperNameInMap + "\"");
    }

}

void setParametersNewtonRaphson(char* mapperName,
                                int maxNumOfIterations, double tolerance){

    std::string mapperNameInMap = std::string(mapperName);

    IGAMortarMapper *tmpIGAMortarMapper;

    // check if the mapper with the given name is generated and is of correct type
    if (!mapperList.count( mapperNameInMap )){
        ERROR_OUT("A mapper with name : " + mapperNameInMap + " does not exist!");
        ERROR_OUT("Did nothing!");
        return;
    } else if (mapperList[mapperNameInMap]->mapperType != EMPIRE_IGAMortarMapper){
        ERROR_OUT(mapperNameInMap + " is not a type of IGAMortarMapper");
        ERROR_OUT("Did nothing!");
        return;
    }
    else{
        tmpIGAMortarMapper = dynamic_cast<IGAMortarMapper *>(mapperList[mapperNameInMap]);
        tmpIGAMortarMapper->setParametersNewtonRaphson(maxNumOfIterations, tolerance);
        INFO_OUT("Point projection Newton Raphson parameters are set for \"" +  mapperNameInMap + "\"");
    }
}

void setParametersNewtonRaphsonBoundary(char* mapperName,
                                        int maxNumOfIterations, double tolerance){

    std::string mapperNameInMap = std::string(mapperName);

    IGAMortarMapper *tmpIGAMortarMapper;

    // check if the mapper with the given name is generated and is of correct type
    if (!mapperList.count( mapperNameInMap )){
        ERROR_OUT("A mapper with name : " + mapperNameInMap + " does not exist!");
        ERROR_OUT("Did nothing!");
        return;
    } else if (mapperList[mapperNameInMap]->mapperType != EMPIRE_IGAMortarMapper){
        ERROR_OUT(mapperNameInMap + " is not a type of IGAMortarMapper");
        ERROR_OUT("Did nothing!");
        return;
    }
    else{
        tmpIGAMortarMapper = dynamic_cast<IGAMortarMapper *>(mapperList[mapperNameInMap]);
        tmpIGAMortarMapper->setParametersNewtonRaphsonBoundary(maxNumOfIterations, tolerance);
        INFO_OUT("Line projection parameters using Newton-Raphson are set for \"" +  mapperNameInMap + "\"");
    }
}

void setParametersBisection(char* mapperName,
                            int maxNumOfIterations, double tolerance){

    std::string mapperNameInMap = std::string(mapperName);

    IGAMortarMapper *tmpIGAMortarMapper;

    // check if the mapper with the given name is generated and is of correct type
    if (!mapperList.count( mapperNameInMap )){
        ERROR_OUT("A mapper with name : " + mapperNameInMap + " does not exist!");
        ERROR_OUT("Did nothing!");
        return;
    } else if (mapperList[mapperNameInMap]->mapperType != EMPIRE_IGAMortarMapper){
        ERROR_OUT(mapperNameInMap + " is not a type of IGAMortarMapper");
        ERROR_OUT("Did nothing!");
        return;
    } else {
        tmpIGAMortarMapper = dynamic_cast<IGAMortarMapper *>(mapperList[mapperNameInMap]);
        tmpIGAMortarMapper->setParametersBisection(maxNumOfIterations, tolerance);
        INFO_OUT("Line projection parameters using Bisection are set for \"" +  mapperNameInMap + "\"");
    }
}

void setParametersIntegration(char* mapperName,
                              int numGPTriangle, int numGPQuad){

    std::string mapperNameInMap = std::string(mapperName);

    IGAMortarMapper *tmpIGAMortarMapper;

    // check if the mapper with the given name is generated and is of correct type
    if (!mapperList.count( mapperNameInMap )){
        ERROR_OUT("A mapper with name : " + mapperNameInMap + " does not exist!");
        ERROR_OUT("Did nothing!");
        return;
    } else if (mapperList[mapperNameInMap]->mapperType != EMPIRE_IGAMortarMapper){
        ERROR_OUT(mapperNameInMap + " is not a type of IGAMortarMapper");
        ERROR_OUT("Did nothing!");
        return;
    }
    else{
        tmpIGAMortarMapper = dynamic_cast<IGAMortarMapper *>(mapperList[mapperNameInMap]);
        tmpIGAMortarMapper->setParametersIntegration(numGPTriangle, numGPQuad);
        INFO_OUT("Integration parameters are set for \"" +  mapperNameInMap + "\"");
    }
}

void setParametersConsistency(char* mapperName,
                              bool _enforceConsistency, double _tolConsistency){

    std::string mapperNameInMap = std::string(mapperName);

    IGAMortarMapper *tmpIGAMortarMapper;

    // check if the mapper with the given name is generated and is of correct type
    if (!mapperList.count( mapperNameInMap )){
        ERROR_OUT("A mapper with name : " + mapperNameInMap + " does not exist!");
        ERROR_OUT("Did nothing!");
        return;
    } else if (mapperList[mapperNameInMap]->mapperType != EMPIRE_IGAMortarMapper){
        ERROR_OUT(mapperNameInMap + " is not a type of IGAMortarMapper");
        ERROR_OUT("Did nothing!");
        return;
    } else {
        tmpIGAMortarMapper = dynamic_cast<IGAMortarMapper *>(mapperList[mapperNameInMap]);
        tmpIGAMortarMapper->setParametersConsistency(_enforceConsistency, _tolConsistency);
        INFO_OUT("Enforce consistency parameter is set for \"" +  mapperNameInMap + "\"");
    }
}

void setParametersIgaWeakDirichletConditions(char* mapperName,
                                             bool _isCurveConditions, bool _isSurfaceConditions,
                                             double _alphaPrim, double _alphaSecBending, double _alphaSecTwisting,
                                             int _isAutomaticPenaltyFactors) {
    std::string mapperNameInMap = std::string(mapperName);

    IGAMortarMapper *tmpIGAMortarMapper;

    // check if the mapper with the given name is generated and is of correct type
    if (!mapperList.count( mapperNameInMap )){
        ERROR_OUT("A mapper with name : " + mapperNameInMap + " does not exist!");
        ERROR_OUT("Did nothing!");
        return;
    } else if (mapperList[mapperNameInMap]->mapperType != EMPIRE_IGAMortarMapper){
        ERROR_OUT(mapperNameInMap + " is not a type of IGAMortarMapper");
        ERROR_OUT("Did nothing!");
        return;
    } else {
        tmpIGAMortarMapper = dynamic_cast<IGAMortarMapper *>(mapperList[mapperNameInMap]);
        tmpIGAMortarMapper->setParametersWeakCurveDirichletConditions(_isCurveConditions, _isAutomaticPenaltyFactors, _alphaPrim, _alphaSecBending, _alphaSecTwisting);
        tmpIGAMortarMapper->setParametersWeakSurfaceDirichletConditions(_isSurfaceConditions, _isAutomaticPenaltyFactors, _alphaPrim);
        INFO_OUT("Dirichlet condition parameters are set for \"" +  mapperNameInMap + "\"");
    }
}

void setParametersIgaPatchCoupling(char* mapperName,
                                   bool _isWeakPatchContinuityConditions, double _alphaPrim,
                                   double _alphaSecBending, double _alphaSecTwisting, int isAutomaticPenaltyFactors) {
    std::string mapperNameInMap = std::string(mapperName);

    IGAMortarMapper *tmpIGAMortarMapper;

    // check if the mapper with the given name is generated and is of correct type
    if (!mapperList.count( mapperNameInMap )){
        ERROR_OUT("A mapper with name : " + mapperNameInMap + " does not exist!");
        ERROR_OUT("Did nothing!");
        return;
    } else if (mapperList[mapperNameInMap]->mapperType != EMPIRE_IGAMortarMapper){
        ERROR_OUT(mapperNameInMap + " is not a type of IGAMortarMapper");
        ERROR_OUT("Did nothing!");
        return;
    } else {
        tmpIGAMortarMapper = dynamic_cast<IGAMortarMapper *>(mapperList[mapperNameInMap]);
        tmpIGAMortarMapper->setParametersWeakPatchContinuityConditions(_isWeakPatchContinuityConditions, isAutomaticPenaltyFactors, _alphaPrim, _alphaSecBending, _alphaSecTwisting);
        INFO_OUT("Patch coupling parameters are set for \"" +  mapperNameInMap + "\"");
    }
}

void setParametersErrorComputation(char* mapperName,
                                   bool _isErrorComputation, bool _isDomainError, bool _isInterfaceError, bool _isCurveError) {
    std::string mapperNameInMap = std::string(mapperName);

    IGAMortarMapper *tmpIGAMortarMapper;

    // check if the mapper with the given name is generated and is of correct type
    if (!mapperList.count( mapperNameInMap )){
        ERROR_OUT("A mapper with name : " + mapperNameInMap + " does not exist!");
        ERROR_OUT("Did nothing!");
        return;
    } else if (mapperList[mapperNameInMap]->mapperType != EMPIRE_IGAMortarMapper){
        ERROR_OUT(mapperNameInMap + " is not a type of IGAMortarMapper");
        ERROR_OUT("Did nothing!");
        return;
    } else {
        tmpIGAMortarMapper = dynamic_cast<IGAMortarMapper *>(mapperList[mapperNameInMap]);
        tmpIGAMortarMapper->setParametersErrorComputation(_isErrorComputation, _isDomainError, _isInterfaceError, _isCurveError);
        INFO_OUT("Error computation parameters are set for \"" +  mapperNameInMap + "\"");
    }
}

void buildCouplingMatrices(char *mapperName){

    std::string mapperNameInMap = std::string(mapperName);

    // check if the mapper with the given name is generated
    if (!mapperList.count( mapperNameInMap )){
        ERROR_OUT("A mapper with name : " + mapperNameInMap + " does not exist!");
        ERROR_OUT("Did nothing!");
        return;
    } else{
        mapperList[mapperNameInMap]->buildCouplingMatrices();
        INFO_OUT("Generated coupling matrices for \"" +  mapperNameInMap );
    }

}

void doConsistentMapping(char* mapperName, int dimension, int dataSizeA, const double* dataA, int dataSizeB, double* dataB){
    assert(dimension == 1 || dimension == 3);

    std::string mapperNameToMap = std::string(mapperName);

    if (!mapperList.count( mapperNameToMap )){
        ERROR_OUT("This mapper does not exist : " + mapperNameToMap);
        ERROR_OUT("Mapping not performed!");
        return;
    } else {
        // if a vector field is to be mapped x, y, z components are extracted
        if (dimension == 3){
            int sizeDataToMap = dataSizeA/dimension;
            int sizeDataToWrite = dataSizeB/dimension;

            double** dataAtoMap = new double*[dimension];
            double** dataBtoWrite = new double*[dimension];

            for (int i=0; i<dimension ; i++){
                dataAtoMap[i]=new double[sizeDataToMap];
                dataBtoWrite[i]=new double[sizeDataToWrite];
            }
            for (int i=0 ; i<dimension ; i++){
                for (int j=0 ; j<sizeDataToMap; j++){
                    dataAtoMap[i][j] = dataA[j*dimension+i];
                }

                mapperList[mapperNameToMap]->consistentMapping(dataAtoMap[i], dataBtoWrite[i]);
                for (int j=0 ; j<sizeDataToWrite; j++){
                    dataB[j*dimension+i] = dataBtoWrite[i][j];
                }
            }
            for (int i = 0; i<dimension; i++){
                delete[] dataAtoMap[i];
                delete[] dataBtoWrite[i];
            }
            delete[] dataAtoMap;
            delete[] dataBtoWrite;
        }
        // else field is mapped as it is
        else {
            mapperList[mapperNameToMap]->consistentMapping(dataA, dataB);
        }
    }
}

void doConservativeMapping(char* mapperName, int dimension, int dataSizeB, const double* dataB, int dataSizeA, double* dataA){
    assert(dimension == 1 || dimension == 3);

    std::string mapperNameToMap = std::string(mapperName);
    if (!mapperList.count( mapperNameToMap )){
        ERROR_OUT("This mapper does not exist : " + mapperNameToMap);
        ERROR_OUT("Mapping not performed!");
        return;
    } else {
        // if a vector field is to be mapped x, y, z components are extracted
        if (dimension == 3){
            int sizeDataToMap = dataSizeB/dimension;
            int sizeDataToWrite = dataSizeA/dimension;

            double** dataBtoMap = new double*[dimension];
            double** dataAtoWrite = new double*[dimension];

            for (int i=0; i<dimension ; i++){
                dataBtoMap[i]=new double[sizeDataToMap];
                dataAtoWrite[i]=new double[sizeDataToWrite];
            }
            for (int i=0 ; i<dimension ; i++){
                for (int j=0 ; j<sizeDataToMap; j++){
                    dataBtoMap[i][j] = dataB[j*dimension+i];
                }
                mapperList[mapperNameToMap]->conservativeMapping(dataBtoMap[i], dataAtoWrite[i]);
                for (int j=0 ; j<sizeDataToWrite; j++){
                    dataA[j*dimension+i] = dataAtoWrite[i][j];
                }
            }
            for (int i = 0; i<dimension; i++){
                delete[] dataBtoMap[i];
                delete[] dataAtoWrite[i];
            }
            delete[] dataBtoMap;
            delete[] dataAtoWrite;
        }
        // else field is mapped as it is
        else {
            mapperList[mapperNameToMap]->conservativeMapping(dataB, dataA);
        }
    }
}

void printMesh(char* meshName){
    std::string meshNameInMap = std::string(meshName);
    if (!meshList.count( meshNameInMap )){
        ERROR_OUT("Mesh with name: \"" + meshNameInMap + "\" does not exist : ");
        ERROR_OUT("Did nothing!");
    } else {
        if (!((meshList[meshNameInMap])->boundingBox.isComputed())) (meshList[meshNameInMap])->computeBoundingBox();

        INFO_OUT(" ======== " + meshNameInMap + " ======== ");
        INFO_OUT() << (meshList[meshNameInMap])->boundingBox << endl;
    }
}

void deleteMesh(char* meshName){
    std::string meshNameInMap = std::string(meshName);
    if (!meshList.count( meshNameInMap )){
        ERROR_OUT("Mesh with name: \"" + meshNameInMap + "\" does not exist : ");
        ERROR_OUT("Did nothing!");
    } else {
        delete meshList[meshNameInMap];
    }
}

void deleteMapper(char* mapperName){
    std::string mapperNameInMap = std::string(mapperName);
    if (!mapperList.count( mapperNameInMap )){
        ERROR_OUT("Mapper with name: \"" + mapperNameInMap + "\" does not exist : ");
        ERROR_OUT("Did nothing!");
        return;
    } else {
        delete mapperList[mapperNameInMap];
    }
}

void deleteAllMappers(){

    std::map<std::string, AbstractMapper*>::iterator iter = mapperList.begin();

    if(iter!=mapperList.end()){
        delete iter->second;
        mapperList.erase(iter);
    }
}
