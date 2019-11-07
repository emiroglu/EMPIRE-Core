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
#ifdef USE_INTEL_MKL
#include <mkl.h>
#include <mkl_lapacke.h>
#include <mkl_spblas.h>
#endif

#ifndef USE_INTEL_MKL
#include <Dense>
#endif

#ifdef FLANN
#include "flann/flann.hpp"
#endif

#ifdef ANN
#include "ANN/ANN.h"
#endif

#include "VertexMorphingMapper.h"
#include "Message.h"
#include <iostream>
#include <stdlib.h>
#include <set>
#include <assert.h>
#include <omp.h>
#include "Message.h"
#include "IGAMesh.h"
#include "FEMesh.h"

using namespace std;
namespace EMPIRE {

// set default value for MKL threads
int VertexMorphingMapper::mklSetNumThreads = 1;
// set default value for mapper threads
int VertexMorphingMapper::mapperSetNumThreads = 1;

double VertexMorphingMapper::EPS_XSI = 1e-6;
double VertexMorphingMapper::Polygon::EPS_PolygonPoint = 1e-6;

VertexMorphingMapper::VertexMorphingMapper(std::string _name, AbstractMesh *_meshA, AbstractMesh *_meshB):
    name(_name){

    // Check input
    assert(_meshA != NULL);
    assert(_meshB != NULL);

    // Assign meshes
    bool isMeshAFEM = (_meshA->type == EMPIRE_Mesh_FEMesh || _meshA->type == EMPIRE_Mesh_copyFEMesh);
    bool isMeshBFEM = (_meshB->type == EMPIRE_Mesh_FEMesh || _meshB->type == EMPIRE_Mesh_copyFEMesh);

    if (isMeshAFEM && isMeshBFEM){
        if (dynamic_cast<FEMesh *>(_meshA)->triangulate() == NULL){
            meshA = dynamic_cast<FEMesh *>(_meshA);
        }
        else
            meshA = dynamic_cast<FEMesh *>(_meshA)->triangulate();

        if (dynamic_cast<FEMesh *>(_meshB)->triangulate() == NULL){
            meshB = dynamic_cast<FEMesh *>(_meshB);
        }
        else
            meshB = dynamic_cast<FEMesh *>(_meshB)->triangulate();
    } else
        ERROR_BLOCK_OUT("VertexMorphingMapper","VertexMorphingMapper","Wrong type of mesh!");

    // Assign the mapper type
    mapperType = EMPIRE_VertexMorphingMapper;

    /// Initializing the sparse matrices
    /// This is a rectangular matrix of size
    /// masterNumNodes X slaveNumNodes
    C_BA = new MathLibrary::SparseMatrix<double>((size_t)(meshB->numNodes),(size_t)(meshA->numNodes));
    C_AB = new MathLibrary::SparseMatrix<double>((size_t)(meshA->numNodes),(size_t)(meshB->numNodes));
    
    // C_BB = new MathLibrary::SparseMatrix<double>((size_t)(meshB->numNodes),(size_t)(meshB->numNodes));

}

VertexMorphingMapper::~VertexMorphingMapper(){

    delete C_BA;
    delete C_AB;
    // delete C_BB;
    delete filterFunction;

#ifdef ANN
    for (int i = 0; i < meshA->numNodes; i++) {
        delete[] ANNSlaveNodes[i];
    }
    delete[] ANNSlaveNodes;
    for (int i = 0; i < meshB->numNodes; i++) {
        delete[] ANNMasterNodes[i];
    }
    delete[] ANNMasterNodes;
#endif

}

void VertexMorphingMapper::deleteTables() {
    for (int i = 0; i < meshA->numElems; i++){
        slaveDirectElemTable[i]->clear();
        delete slaveDirectElemTable[i];
    }
    delete[] slaveDirectElemTable;
    for (int i = 0; i < meshB->numElems; i++){
        masterDirectElemTable[i]->clear();
        delete masterDirectElemTable[i];
    }
    delete[] masterDirectElemTable;

    for (int i = 0; i < meshA->numNodes; i++){
        slaveNodeToElemTable[i]->clear();
        delete slaveNodeToElemTable[i];
    }
    delete[] slaveNodeToElemTable;
    for (int i = 0; i < meshB->numNodes; i++){
        masterNodeToElemTable[i]->clear();
        delete masterNodeToElemTable[i];
    }
    delete[] masterNodeToElemTable;

    // delete[] slaveElemToMasterNodeTable;

    delete[] masterFilterFunctionIntegrationOnSlave;
    delete[] slaveFilterFunctionIntegrationOnMaster;

}

void VertexMorphingMapper::deleteANNTree() {
#ifdef ANN
    delete[] ANNSlaveNodes;
    delete slaveNodesTree;
    delete[] ANNMasterNodes;
    delete masterNodesTree;
    annClose();
#endif
#ifdef FLANN
    delete FLANNSlaveNodes;
    delete FLANNkd_slaveTree;
    delete FLANNMasterNodes;
    delete FLANNkd_masterTree;
#endif
}

void VertexMorphingMapper::setParameters(EMPIRE_VMM_FilterType _filterType, double _filterRadius){

    // Set the filter radius for the search radius
    filterRadius = _filterRadius;

    // Set the filter function pointer and the num GPs for integration
    if (_filterType == EMPIRE_VMM_HatFilter) {
        filterFunction = new HatFilterFunction(filterRadius);
        numGPsOnTri = 3;
        numGPsOnQuad = 4;
    } else if (_filterType == EMPIRE_VMM_GaussianFilter) {
        filterFunction = new GaussianFilterFunction(filterRadius);
        numGPsOnTri = 12;
        numGPsOnQuad = 9;
    } else assert(false); // shouldnt come here

}

void VertexMorphingMapper::initialize(){

    // Initialize data that could be used later
    initANNTree();
    initTables();

}

void VertexMorphingMapper::initTables() {

    // using the map to store the nodeNumbers
    // but here the "key" is the node number, and the value is the position in nodeNumbers
    // the map is sorted automatically, so it is efficient for searching

    { // 1. compute slaveDirectElemTable
        slaveDirectElemTable = new vector<int>*[meshA->numElems];
        for (int i = 0; i < meshA->numElems; i++)
            slaveDirectElemTable[i] = new vector<int>;
        map<int, int> *slaveNodesMap = new map<int, int>();
        for (int i = 0; i < meshA->numNodes; i++)
            slaveNodesMap->insert(slaveNodesMap->end(), pair<int, int>(meshA->nodeIDs[i], i));
        int count = 0;
        for (int i = 0; i < meshA->numElems; i++) {
            const int numNodesSlaveElem = meshA->numNodesPerElem[i];
            for (int j = 0; j < numNodesSlaveElem; j++) {
                if(slaveNodesMap->find(meshA->elems[count + j])== slaveNodesMap->end())
                    ERROR_OUT()<< "Slave Node Label " << meshA->elems[count + j] << " is not part of slaveNodesMap." << endl;

                slaveDirectElemTable[i]->push_back(slaveNodesMap->at(meshA->elems[count + j]));
            }
            count += numNodesSlaveElem;
        }
        delete slaveNodesMap;
    }

    { // 2. compute slaveNodeToElemTable
        slaveNodeToElemTable = new vector<int>*[meshA->numNodes];
        for (int i = 0; i < meshA->numNodes; i++)
            slaveNodeToElemTable[i] = new vector<int>;
        for (int i = 0; i < meshA->numElems; i++) {
            const int numNodesSlaveElem = meshA->numNodesPerElem[i];
            for (int j = 0; j < numNodesSlaveElem; j++) {
                int nodePos = slaveDirectElemTable[i]->at(j);
                slaveNodeToElemTable[nodePos]->push_back(i);
            }
        }
    }

    { // 3. compute masterDirectElemTable
        masterDirectElemTable = new vector<int>*[meshB->numElems];
        for (int i = 0; i < meshB->numElems; i++)
            masterDirectElemTable[i] = new vector<int>;
        map<int, int> *masterNodesMap = new map<int, int>();
        for (int i = 0; i < meshB->numNodes; i++)
            masterNodesMap->insert(masterNodesMap->end(), pair<int, int>(meshB->nodeIDs[i], i));
        int count = 0;
        for (int i = 0; i < meshB->numElems; i++) {
            const int numNodesMasterElem = meshB->numNodesPerElem[i];
            for (int j = 0; j < numNodesMasterElem; j++) {
                if(masterNodesMap->find(meshB->elems[count + j])== masterNodesMap->end())
                    ERROR_OUT()<< "Master Node Label " << meshB->elems[count + j] << " is not part of masterNodesMap." << endl;

                masterDirectElemTable[i]->push_back(masterNodesMap->at(meshB->elems[count + j]));
            }
            count += numNodesMasterElem;
        }
        delete masterNodesMap;
    }

    { // 4. compute masterNodeToElemTable
        masterNodeToElemTable = new vector<int>*[meshB->numNodes];
        for (int i = 0; i < meshB->numNodes; i++)
            masterNodeToElemTable[i] = new vector<int>;
        for (int i = 0; i < meshB->numElems; i++) {
            const int numNodesMasterElem = meshB->numNodesPerElem[i];
            for (int j = 0; j < numNodesMasterElem; j++) {
                int nodePos = masterDirectElemTable[i]->at(j);
                masterNodeToElemTable[nodePos]->push_back(i);
            }
        }
    }
    
    masterFilterFunctionIntegrationOnSlave = new double[meshB->numNodes];
    for (int iNode = 0; iNode < meshB->numNodes; iNode++)
        masterFilterFunctionIntegrationOnSlave[iNode] = 0.0;

    slaveFilterFunctionIntegrationOnMaster = new double[meshA->numNodes];
    for (int iNode = 0; iNode < meshA->numNodes; iNode++)
        slaveFilterFunctionIntegrationOnMaster[iNode] = 0.0;

    // This is necessary for C_BB
    // slaveElemToMasterNodeTable = new vector<int>[meshA->numElems];

}


void VertexMorphingMapper::initANNTree() {

#ifdef ANN
    ANNSlaveNodes = new double*[meshA->numNodes]; // ANN uses 2D array
    for (int i = 0; i < meshA->numNodes; i++) {
        ANNSlaveNodes[i] = new double[3];
        for (int j = 0; j<3; j++)
            ANNSlaveNodes[i][j] = meshA->nodes[i * 3 + j];
    }
    slaveNodesTree = new ANNkd_tree(ANNSlaveNodes, meshA->numNodes, 3);

    ANNMasterNodes = new double*[meshB->numNodes]; // ANN uses 2D array
    for (int i = 0; i < meshB->numNodes; i++) {
        ANNMasterNodes[i] = new double[3];
        for (int j = 0; j<3; j++)
            ANNMasterNodes[i][j] = meshB->nodes[i * 3 + j];
    }
    masterNodesTree = new ANNkd_tree(ANNMasterNodes, meshB->numNodes, 3);

#endif
#ifdef FLANN
    FLANNSlaveNodes = new flann::Matrix<double>(const_cast<double*>(meshA->nodes), meshA->numNodes, 3);
    FLANNkd_slaveTree = new flann::Index<flann::L2<double> >(*FLANNSlaveNodes, flann::KDTreeSingleIndexParams(1));
    FLANNkd_slaveTree->buildIndex(); // Build binary tree for searching

    FLANNMasterNodes = new flann::Matrix<double>(const_cast<double*>(meshB->nodes), meshB->numNodes, 3);
    FLANNkd_masterTree = new flann::Index<flann::L2<double> >(*FLANNMasterNodes, flann::KDTreeSingleIndexParams(1));
    FLANNkd_masterTree->buildIndex(); // Build binary tree for searching

#endif
}

void VertexMorphingMapper::findSlaveNodes(const int _masterNodeIdx, vector<int>& _slaveNodeIdxs){

    const int dim = 3;
    double* masterNodeCoords = &meshB->nodes[_masterNodeIdx*dim];

    findSlaveNodes(masterNodeCoords, _slaveNodeIdxs);

}

void VertexMorphingMapper::findMasterNodes(const int _slaveNodeIdx, vector<int>& _masterNodeIdxs){

    const int dim = 3;
    double* slaveNodeCoords = &meshA->nodes[_slaveNodeIdx*dim];

    findMasterNodes(slaveNodeCoords, _masterNodeIdxs);

}

void VertexMorphingMapper::findSlaveNodes(double* _masterNodeCoords, vector<int>& _slaveNodeIdxs){

    // 1.1 Find the node indices in the filter radius
    int dim = 3;
    int numInfNodes = 0;
    int *infNodeIndices;

#ifdef ANN
    numInfNodes = slaveNodesTree->annkFRSearch(_masterNodeCoords, filterRadius, 0); // get the number of neighbors in a radius
    infNodeIndices = new int[numInfNodes];
    double *dummy = new double[numInfNodes];
    slaveNodesTree->annkFRSearch(&_masterNodeCoords, filterRadius*filterRadius, numInfNodes, infNodeIndices, dummy);// get the real neighbors (ANN uses the square of the radius)
    delete[] dummy;
#endif

#ifdef FLANN
    flann::Matrix<double> masterNodeCopyFlann(_masterNodeCoords, 1, 3);
    vector<vector<int> > indices_tmp;
    vector<vector<double> > dists_tmp;
    FLANNkd_slaveTree->radiusSearch(masterNodeCopyFlann, indices_tmp, dists_tmp, filterRadius*filterRadius, flann::SearchParams(1));
    numInfNodes = indices_tmp[0].size();
    infNodeIndices = new int[indices_tmp[0].size()];
    for (int j=0; j<indices_tmp[0].size(); j++)
        infNodeIndices[j] = indices_tmp[0][j];
#endif

    for (int iNodeIdx = 0; iNodeIdx < numInfNodes; iNodeIdx++)
        _slaveNodeIdxs.push_back(infNodeIndices[iNodeIdx]);
    
    delete[] infNodeIndices;
}

void VertexMorphingMapper::findMasterNodes(double* _slaveNodeCoords, vector<int>& _masterNodeIdxs){

    // 1.1 Find the node indices in the filter radius
    int dim = 3;
    int numInfNodes = 0;
    int *infNodeIndices;

#ifdef ANN
    numInfNodes = masterNodesTree->annkFRSearch(_slaveNodeCoords, filterRadius, 0); // get the number of neighbors in a radius
    infNodeIndices = new int[numInfNodes];
    double *dummy = new double[numInfNodes];
    masterNodesTree->annkFRSearch(&_slaveNodeCoords, filterRadius*filterRadius, numInfNodes, infNodeIndices, dummy);// get the real neighbors (ANN uses the square of the radius)
    delete[] dummy;
#endif

#ifdef FLANN
    flann::Matrix<double> slaveNodeCopyFlann(_slaveNodeCoords, 1, 3);
    vector<vector<int> > indices_tmp;
    vector<vector<double> > dists_tmp;
    FLANNkd_masterTree->radiusSearch(slaveNodeCopyFlann, indices_tmp, dists_tmp, filterRadius*filterRadius, flann::SearchParams(1));
    numInfNodes = indices_tmp[0].size();
    infNodeIndices = new int[indices_tmp[0].size()];
    for (int j=0; j<indices_tmp[0].size(); j++)
        infNodeIndices[j] = indices_tmp[0][j];
#endif

    for (int iNodeIdx = 0; iNodeIdx < numInfNodes; iNodeIdx++)
        _masterNodeIdxs.push_back(infNodeIndices[iNodeIdx]);
    
    delete[] infNodeIndices;
}

void VertexMorphingMapper::findSlaveElements(const int _masterNodeIdx, vector<int>& _slaveElemIdxs){

    vector<int> slaveNodeIdxs;
    findSlaveNodes(_masterNodeIdx, slaveNodeIdxs);

    // Loop over the influenced slave nodes
    for (vector<int>::iterator itSlaveNodeIdx = slaveNodeIdxs.begin(); itSlaveNodeIdx != slaveNodeIdxs.end(); itSlaveNodeIdx++) {
        
        // Get the master elems containing this slave node
        vector<int>* slaveNodeElems = slaveNodeToElemTable[*itSlaveNodeIdx];

        // Loop over the influenced slave elems list and add this element if it is not in the list already
        for (vector<int>::iterator itSlaveElem = slaveNodeElems->begin(); itSlaveElem != slaveNodeElems->end(); itSlaveElem++){
            if (find(_slaveElemIdxs.begin(), _slaveElemIdxs.end(), *itSlaveElem) ==  _slaveElemIdxs.end()){
                _slaveElemIdxs.push_back(*itSlaveElem);
            }
        }
    }

    // Add the master nodes to the related slave element table
    // This is necessary for the computation of C_BB
    // for (vector<int>::iterator itSlaveElemIdxs = _slaveElemIdxs.begin(); itSlaveElemIdxs != _slaveElemIdxs.end(); itSlaveElemIdxs++){
    //     slaveElemToMasterNodeTable[*itSlaveElemIdxs].push_back(_masterNodeIdx);
    // }
}

void VertexMorphingMapper::findMasterElements(const int _slaveNodeIdx, vector<int>& _masterElemIdxs){

    vector<int> masterNodeIdxs;
    findMasterNodes(_slaveNodeIdx, masterNodeIdxs);

    // Loop over the influenced slave nodes
    for (vector<int>::iterator itMasterNodeIdx = masterNodeIdxs.begin(); itMasterNodeIdx != masterNodeIdxs.end(); itMasterNodeIdx++) {
        
        // Get the master elems containing this slave node
        vector<int>* masterNodeElems = masterNodeToElemTable[*itMasterNodeIdx];

        // Loop over the influenced slave elems list and add this element if it is not in the list already
        for (vector<int>::iterator itMasterElem = masterNodeElems->begin(); itMasterElem != masterNodeElems->end(); itMasterElem++){
            if (find(_masterElemIdxs.begin(), _masterElemIdxs.end(), *itMasterElem) ==  _masterElemIdxs.end()){
                _masterElemIdxs.push_back(*itMasterElem);
            }
        }
    }

    // Add the master nodes to the related slave element table
    // This is necessary for the computation of C_BB
    // for (vector<int>::iterator itSlaveElemIdxs = _slaveElemIdxs.begin(); itSlaveElemIdxs != _slaveElemIdxs.end(); itSlaveElemIdxs++){
    //     slaveElemToMasterNodeTable[*itSlaveElemIdxs].push_back(_masterNodeIdx);
    // }
}


void VertexMorphingMapper::build_C_BA(){

    int dim = 3;

    vector<int> slaveElemIdxs;

    // Loop over the master nodes
    for (int masterNodeIdx = 0; masterNodeIdx < meshB->numNodes; masterNodeIdx++){

        // cout << "Master Node Idx: " << masterNodeIdx << endl;
        // Find influenced master elements
        findSlaveElements(masterNodeIdx, slaveElemIdxs);
        // cout << "Slave Element Idxs: ";
        // for (int i =0; i<slaveElemIdxs.size(); i++)
        //     cout << slaveElemIdxs.at(i) << " ";
        // cout << endl;

        // Loop over the influenced elements and integrate
        for (vector<int>::iterator itSlaveElemIdx = slaveElemIdxs.begin(); itSlaveElemIdx != slaveElemIdxs.end(); itSlaveElemIdx++)
            doIntegration_C_BA(*itSlaveElemIdx, masterNodeIdx);

        slaveElemIdxs.clear();
    }

    adjustFilterFunctions_BA();
}

void VertexMorphingMapper::build_C_AB(){

    int dim = 3;

    vector<int> masterElemIdxs;

    // Loop over the master nodes
    for (int slaveNodeIdx = 0; slaveNodeIdx < meshA->numNodes; slaveNodeIdx++){

        // Find influenced master elements
        findMasterElements(slaveNodeIdx, masterElemIdxs);

        // Loop over the influenced elements and integrate
        for (vector<int>::iterator itMasterElemIdx = masterElemIdxs.begin(); itMasterElemIdx != masterElemIdxs.end(); itMasterElemIdx++)
            doIntegration_C_AB(*itMasterElemIdx, slaveNodeIdx);

        masterElemIdxs.clear();
    }

    adjustFilterFunctions_AB();
}

void VertexMorphingMapper::build_C_BB(){

    // int numTriaNodes = 3;
    // int dim = 3;

    // // Build C_BB with the filter functions -> F^T . F
    // // Loop over the slave elements
    // for (int slaveElemIdx = 0; slaveElemIdx < meshA->numElems; slaveElemIdx++){

    //     std::cout << "SlaveElem: " << slaveElemIdx+1 << std::endl;
        
    //     vector<int>* slaveElemMasterNodes = &slaveElemToMasterNodeTable[slaveElemIdx];

    //     int numSlaveElemNodes = meshA->numNodesPerElem[slaveElemIdx];

    //     // Get the element to be integrated
    //     double slaveElem[numSlaveElemNodes*dim];
    //     for (int iNode = 0; iNode<numSlaveElemNodes; iNode++){
    //         int nodeIdx = slaveDirectElemTable[slaveElemIdx]->at(iNode);
    //         for (int iXYZ = 0; iXYZ<dim; iXYZ++)
    //             slaveElem[iNode*dim+iXYZ] = meshA->nodes[nodeIdx*dim+iXYZ];
    //     }

    //     // Loop over the influencing master nodes
    //     for (int ctrISlaveElemMasterNode = 0; ctrISlaveElemMasterNode < slaveElemMasterNodes->size(); ctrISlaveElemMasterNode++){
    //     // for (vector<int>::iterator itISlaveElemMasterNode = slaveElemMasterNodes->begin(); itISlaveElemMasterNode != slaveElemMasterNodes->end(); itISlaveElemMasterNode++){

    //         int iSlaveElemMasterNodeIdx = slaveElemMasterNodes->at(ctrISlaveElemMasterNode);
    //         // std::cout << iSlaveElemMasterNodeIdx << std::endl;
            
    //         VertexMorphingMapper::Polygon intPolygon = VertexMorphingMapper::Polygon();
            
    //         // Clip element with iMasterNode
    //         clipElementWithFilterRadius(iSlaveElemMasterNodeIdx, numSlaveElemNodes, slaveElem, intPolygon);
    //         intPolygon.triangulate();

    //         // Loop over the subtriangles of the polygon
    //         for (std::vector<double*>::iterator itTria = intPolygon.triangles.begin(); itTria !=intPolygon.triangles.end(); itTria++){
                
    //             // Loop over the influencing master nodes again
    //             for (int ctrJSlaveElemMasterNode = ctrISlaveElemMasterNode; ctrJSlaveElemMasterNode < slaveElemMasterNodes->size(); ctrJSlaveElemMasterNode++){
    //             // for (vector<int>::iterator itJSlaveElemMasterNode = slaveElemMasterNodes->begin(); itJSlaveElemMasterNode != slaveElemMasterNodes->end(); itJSlaveElemMasterNode++){

    //                 int jSlaveElemMasterNodeIdx = slaveElemMasterNodes->at(ctrJSlaveElemMasterNode);
                    
    //                 if (iSlaveElemMasterNodeIdx == jSlaveElemMasterNodeIdx){
    //                     // std::cout << iSlaveElemMasterNodeIdx << "==" << jSlaveElemMasterNodeIdx << std::endl;
    //                     // Create Gauss rule on subtriangle
    //                     EMPIRE::MathLibrary::GaussQuadratureOnTriangle *gaussQuadratureOnTria =
    //                             new EMPIRE::MathLibrary::GaussQuadratureOnTriangle(*itTria, 12);
    //                     for (int iGP = 0; iGP < gaussQuadratureOnTria->getNumGaussPoints(); iGP++) {
    //                         (*C_BB)(iSlaveElemMasterNodeIdx, jSlaveElemMasterNodeIdx) += 
    //                             gaussQuadratureOnTria->getWeight(iGP)*
    //                             filterFunction->computeFunction(&meshB->nodes[iSlaveElemMasterNodeIdx*dim], &gaussQuadratureOnTria->getGaussPointsGlobal()[iGP*dim])*
    //                             filterFunction->computeFunction(&meshB->nodes[jSlaveElemMasterNodeIdx*dim], &gaussQuadratureOnTria->getGaussPointsGlobal()[iGP*dim])*
    //                             gaussQuadratureOnTria->getDetJ(iGP);
    //                             // *
    //                             // (1.0/masterFilterFunctionIntegrationOnSlave[iSlaveElemMasterNodeIdx])*
    //                             // (1.0/masterFilterFunctionIntegrationOnSlave[jSlaveElemMasterNodeIdx]);
    //                     }
    //                 } else {

    //                     // std::cout << iSlaveElemMasterNodeIdx << "!=" << jSlaveElemMasterNodeIdx << std::endl;

    //                     VertexMorphingMapper::Polygon intSubPolygon = VertexMorphingMapper::Polygon();

    //                     // Clip triangle with the jMasterNode
    //                     clipElementWithFilterRadius(jSlaveElemMasterNodeIdx, numTriaNodes, *itTria, intSubPolygon);
    //                     intSubPolygon.triangulate();

    //                     // Loop over the subPolygonSubTria
    //                     for (std::vector<double*>::iterator itSubTria = intSubPolygon.triangles.begin(); itSubTria !=intSubPolygon.triangles.end(); itSubTria++){

    //                         // Create Gauss rule on subtriangle
    //                         EMPIRE::MathLibrary::GaussQuadratureOnTriangle *gaussQuadratureOnTria =
    //                                 new EMPIRE::MathLibrary::GaussQuadratureOnTriangle(*itSubTria, 12);
    //                         for (int iGP = 0; iGP < gaussQuadratureOnTria->getNumGaussPoints(); iGP++) {
    //                             double entryIJ = gaussQuadratureOnTria->getWeight(iGP)*
    //                                 filterFunction->computeFunction(&meshB->nodes[iSlaveElemMasterNodeIdx*dim], &gaussQuadratureOnTria->getGaussPointsGlobal()[iGP*dim])*
    //                                 filterFunction->computeFunction(&meshB->nodes[jSlaveElemMasterNodeIdx*dim], &gaussQuadratureOnTria->getGaussPointsGlobal()[iGP*dim])*
    //                                 gaussQuadratureOnTria->getDetJ(iGP);
    //                                 // *
    //                                 // (1.0/masterFilterFunctionIntegrationOnSlave[iSlaveElemMasterNodeIdx])*
    //                                 // (1.0/masterFilterFunctionIntegrationOnSlave[jSlaveElemMasterNodeIdx]);
    //                             (*C_BB)(iSlaveElemMasterNodeIdx, jSlaveElemMasterNodeIdx) += entryIJ;
    //                             (*C_BB)(jSlaveElemMasterNodeIdx, iSlaveElemMasterNodeIdx) += entryIJ;
    //                         }

    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }

    // cout << "VertexMorphingMapper::build_C_BB()" << endl;

    // // This function is copied and modified following dataFieldIntegration.cpp

    // // compute the sparsity map
    // // sparsity map has the information of a, ia, ja in a CSR formated matrix
    // map<int, double> **sparsityMap = new map<int, double>*[meshB->numNodes];
    // for (int i = 0; i < meshB->numNodes; i++)
    //     sparsityMap[i] = new map<int, double>;

    // for (int i = 0; i < meshB->numElems; i++) {
    //     const int numNodesThisElem = meshB->numNodesPerElem[i];
    //     double elem[numNodesThisElem * 3]; // this element
    //     int pos[numNodesThisElem]; // the position of the node in the nodeCoors
    //     for (int j = 0; j < numNodesThisElem; j++) {
    //         pos[j] = masterDirectElemTable[i]->at(j);
    //         for (int k = 0; k < 3; k++)
    //             elem[j * 3 + k] = meshB->nodes[pos[j] * 3 + k];
    //     }

    //     // replace the master element by the projection of it on its "element plane"
    //     // we do it here because we have done the same in MortarMapper
    //     if (numNodesThisElem == 4) {
    //         double masterElemNormal[3];
    //         EMPIRE::MathLibrary::computeNormalOfQuad(elem, true, masterElemNormal);
    //         double masterQuadCenter[3];
    //         EMPIRE::MathLibrary::computePolygonCenter(elem, 4, masterQuadCenter);
    //         double masterQuadPrj[12];
    //         EMPIRE::MathLibrary::projectToPlane(masterQuadCenter, masterElemNormal, elem, 4, masterQuadPrj);
    //         for (int i = 0; i < 12; i++)
    //             elem[i] = masterQuadPrj[i];
    //     }

    //     // make use of the symmetry
    //     double massMatrixElem[numNodesThisElem * numNodesThisElem];
    //     if (numNodesThisElem == 4)
    //     	EMPIRE::MathLibrary::computeMassMatrixOfQuad(elem, numGPsOnQuad, false, massMatrixElem);
    //     else if (numNodesThisElem == 3)
    //     	EMPIRE::MathLibrary::computeMassMatrixOfTrianlge(elem, numGPsOnTri, false, massMatrixElem);
    //     else
    //         assert(false);
    //     for (int j = 0; j < numNodesThisElem; j++) {
    //         for (int k = j; k < numNodesThisElem; k++) {
    //             int smaller, larger;
    //             if (pos[j] > pos[k]) {
    //                 larger = pos[j];
    //                 smaller = pos[k];
    //             } else {
    //                 larger = pos[k];
    //                 smaller = pos[j];
    //             }
    //             double massMatrixJK = massMatrixElem[j * numNodesThisElem + k];
    //             bool inserted = sparsityMap[smaller]->insert(
    //                     pair<int, double>(larger, massMatrixJK)).second; // *.second is a bool indicating inserted or not
    //             if (!inserted)
    //                 (sparsityMap[smaller]->at(larger)) += massMatrixJK; // at() returns a reference, so using "+=" is correct
    //         }
    //     }
    // }

    // // 2. according to sparsity map, compute a, ia, ja of CSR formated massMatrix
    // int nnz = 1; // number of non-zero entries
    // for (int i = 0; i < meshB->numNodes; i++) {
    //     nnz += sparsityMap[i]->size();
    // }
    // nnz--;
    // int count2 = 0;
    // for (int i = 0; i < meshB->numNodes; i++) {
    //     for (map<int, double>::iterator it = sparsityMap[i]->begin(); it != sparsityMap[i]->end();
    //             it++) {
    //     	// Edit Aditya
    //     	(*C_BB)(i,it->first) = it->second;
    //     	(*C_BB)(it->first,i) = it->second;
    //         count2++;
    //     }
    // }assert(nnz == count2);

    // // delete spasity map
    // for (int i = 0; i < meshB->numNodes; i++)
    //     delete sparsityMap[i];
    // delete[] sparsityMap;

}

void VertexMorphingMapper::buildCouplingMatrices(){

    initialize();

    build_C_BA();

    // build_C_AB();

    // build_C_BB();

}

void VertexMorphingMapper::doIntegration_C_BA(const int _slaveElemIdx, const int _masterNodeIdx)
{

    int dim = 3;
    int numNodes = meshA->numNodesPerElem[_slaveElemIdx];

    // Get the element to be integrated
    double slaveElem[meshA->numNodesPerElem[_slaveElemIdx]*dim];
    for (int iNode = 0; iNode<meshA->numNodesPerElem[_slaveElemIdx]; iNode++){
        int nodeIdx = slaveDirectElemTable[_slaveElemIdx]->at(iNode);
        for (int iXYZ = 0; iXYZ<dim; iXYZ++)
            slaveElem[iNode*dim+iXYZ] = meshA->nodes[nodeIdx*dim+iXYZ];
    }

    // Make an integration polygon to triangulate afterwards
    VertexMorphingMapper::Polygon* intPolygon = new VertexMorphingMapper::Polygon();

    // Find the clipping of the element and store in the integration polygon
    clipSlaveElementWithFilterRadius(_masterNodeIdx, numNodes, slaveElem, intPolygon);

    // Triangulate the polygon
    intPolygon->triangulate();

    // Loop over the subtriangles of the polygon
    for (std::vector<double*>::iterator itTria = intPolygon->triangles.begin(); itTria !=intPolygon->triangles.end(); itTria++){
        // Create Gauss rule on subtriangle
        EMPIRE::MathLibrary::GaussQuadratureOnTriangle *gaussQuadratureOnTria =
                new EMPIRE::MathLibrary::GaussQuadratureOnTriangle(*itTria, numGPsOnTri);
        // Compute the integration of the filter function
        // It is done this way instead of doing it on the matrix directly because it is a very expensive operation on the matrix
        for (int iGP = 0; iGP < gaussQuadratureOnTria->getNumGaussPoints(); iGP++) {
            masterFilterFunctionIntegrationOnSlave[_masterNodeIdx] += gaussQuadratureOnTria->getWeight(iGP)
                    * filterFunction->computeFunction(&meshB->nodes[_masterNodeIdx*dim], &gaussQuadratureOnTria->getGaussPointsGlobal()[iGP*dim])
                    * gaussQuadratureOnTria->getDetJ(iGP);
        }

        // Create the integrand on the unclipped element (This can be generated outside and passed in)
        FilterFunctionProduct* integrand = new FilterFunctionProduct(meshA->numNodesPerElem[_slaveElemIdx], slaveElem, &meshB->nodes[_masterNodeIdx*dim]);

        integrand->setGaussPoints(gaussQuadratureOnTria->gaussPointsGlobal, gaussQuadratureOnTria->numGaussPoints);
        integrand->computeFunctionProducts(filterFunction);

        for (int iNode = 0; iNode < meshA->numNodesPerElem[_slaveElemIdx]; iNode++) {
            gaussQuadratureOnTria->setIntegrandFunc(integrand);
            integrand->setFunctionID(iNode);
            (*C_BA)(_masterNodeIdx,slaveDirectElemTable[_slaveElemIdx]->at(iNode)) += gaussQuadratureOnTria->computeIntegral();
        }

        delete integrand;
        delete gaussQuadratureOnTria;

    }

    delete intPolygon;

}

void VertexMorphingMapper::doIntegration_C_AB(const int _masterElemIdx, const int _slaveNodeIdx)
{

    int dim = 3;
    int numNodes = meshB->numNodesPerElem[_masterElemIdx];

    // Get the element to be integrated
    double masterElem[meshB->numNodesPerElem[_masterElemIdx]*dim];
    for (int iNode = 0; iNode<meshB->numNodesPerElem[_masterElemIdx]; iNode++){
        int nodeIdx = masterDirectElemTable[_masterElemIdx]->at(iNode);
        for (int iXYZ = 0; iXYZ<dim; iXYZ++)
            masterElem[iNode*dim+iXYZ] = meshB->nodes[nodeIdx*dim+iXYZ];
    }

    // Make an integration polygon to triangulate afterwards
    VertexMorphingMapper::Polygon* intPolygon = new VertexMorphingMapper::Polygon();

    // Find the clipping of the element and store in the integration polygon
    clipMasterElementWithFilterRadius(_slaveNodeIdx, numNodes, masterElem, intPolygon);

    // Triangulate the polygon
    intPolygon->triangulate();

    // Loop over the subtriangles of the polygon
    for (std::vector<double*>::iterator itTria = intPolygon->triangles.begin(); itTria !=intPolygon->triangles.end(); itTria++){
        // Create Gauss rule on subtriangle
        EMPIRE::MathLibrary::GaussQuadratureOnTriangle *gaussQuadratureOnTria =
                new EMPIRE::MathLibrary::GaussQuadratureOnTriangle(*itTria, numGPsOnTri);
        // Compute the integration of the filter function
        // It is done this way instead of doing it on the matrix directly because it is a very expensive operation on the matrix
        for (int iGP = 0; iGP < gaussQuadratureOnTria->getNumGaussPoints(); iGP++) {
            slaveFilterFunctionIntegrationOnMaster[_slaveNodeIdx] += gaussQuadratureOnTria->getWeight(iGP)
                    * filterFunction->computeFunction(&meshA->nodes[_slaveNodeIdx*dim], &gaussQuadratureOnTria->getGaussPointsGlobal()[iGP*dim])
                    * gaussQuadratureOnTria->getDetJ(iGP);
        }

        // Create the integrand on the unclipped element (This can be generated outside and passed in)
        FilterFunctionProduct* integrand = new FilterFunctionProduct(meshB->numNodesPerElem[_masterElemIdx], masterElem, &meshA->nodes[_slaveNodeIdx*dim]);

        integrand->setGaussPoints(gaussQuadratureOnTria->gaussPointsGlobal, gaussQuadratureOnTria->numGaussPoints);
        integrand->computeFunctionProducts(filterFunction);

        for (int iNode = 0; iNode < meshB->numNodesPerElem[_masterElemIdx]; iNode++) {
            gaussQuadratureOnTria->setIntegrandFunc(integrand);
            integrand->setFunctionID(iNode);
            (*C_AB)(_slaveNodeIdx,masterDirectElemTable[_masterElemIdx]->at(iNode)) += gaussQuadratureOnTria->computeIntegral();
        }

        delete integrand;
        delete gaussQuadratureOnTria;

    }

    delete intPolygon;

}

void VertexMorphingMapper::clipSlaveElementWithFilterRadius(const int _masterNodeIdx, const int _numNodes, double* _elem, VertexMorphingMapper::Polygon* _polygon){

    // This function considers the clipping cases one by one instead of trying to clip each edge by default.
    // This way turned out to be more efficient than the latter.

    int dim = 3;

    // Count the influenced nodes by the master node
    int numInNodes = 0;
    std::vector<int> inNodesPos;
    std::vector<int> outNodesPos;
    int nodePos = -1;

    // Find the inside and outside nodes of this element
    for (int iNode = 0; iNode < _numNodes; iNode++){
        nodePos++;
        double dist = EMPIRE::MathLibrary::computeDenseEuclideanNorm(dim, &meshB->nodes[_masterNodeIdx*dim], &_elem[iNode*dim]);
        if (dist < filterRadius){
            inNodesPos.push_back(nodePos);
            numInNodes++;
        } else
            outNodesPos.push_back(nodePos);
    }

    // triaCase = 1: one node inside
    // triaCase = 2: two nodes inside
    // triaCase = 3: all nodes inside
    // quadCase = 1: one node inside
    // quadCase = 2: two adjacent nodes inside
    // quadCase = 3: three nodes inside
    // quadCase = 4: all nodes inside
    // quadCase = 5: two nodes accross each other inside
    int triaCase = 0;
    int quadCase = 0;
    if (_numNodes == 3) {
        triaCase += numInNodes;
    }
    else if (_numNodes == 4){
        quadCase += numInNodes;
        if (numInNodes == 2 && (abs(inNodesPos[0]-inNodesPos[1]) == 2))
            quadCase += 3;
    } else assert(false);

    // triaCase == 1: one node inside
    if (triaCase == 1){
        // cout << "triaCase == 1" << endl;
        double* P0 = &_elem[inNodesPos.at(0)*dim];  // inside point
        double* P1 = &_elem[outNodesPos.at(0)*dim]; // next point
        double* P2 = &_elem[outNodesPos.at(1)*dim]; // prev point
        // cout << "got points" << endl;
        
        // xsi01 and xsi02 must have only 1 clipping
        vector<double> xsi01;
        vector<double> xsi02;
        // xsi12 might have 1, 2 or no clippings
        vector<double> xsi12;

        // Find clippings
        // cout << "finding clippings" << endl;
        bool isClipping01 = findClippingWithMaster(_masterNodeIdx, P0, P1, xsi01);  // P01
        bool isClipping02 = findClippingWithMaster(_masterNodeIdx, P0, P2, xsi02);  // P02
        bool isClipping12 = findClippingWithMaster(_masterNodeIdx, P1, P2, xsi12);  // P12
        // cout << "found clippings " << xsi01.size() << " " << xsi02.size() << " " << xsi12.size() << endl;

        // Compute the global cartesian coordinates of the clipping points
        double P01[3];
        double P02[3];
        if (isClipping01){
            for (int iXYZ = 0; iXYZ < dim; iXYZ++)
                P01[iXYZ] = (1.0-xsi01.at(0)) * P0[iXYZ] + xsi01.at(0) * P1[iXYZ];
        }
        if (isClipping02){
            for (int iXYZ = 0; iXYZ < dim; iXYZ++)
                P02[iXYZ] = (1.0-xsi02.at(0)) * P0[iXYZ] + xsi02.at(0) * P2[iXYZ];
        }
        xsi01.clear();
        xsi02.clear();

        // for (int iXYZ = 0; iXYZ < dim; iXYZ++){
        //     P01[iXYZ] = (1.0-xsi01.at(0)) * P0[iXYZ] + xsi01.at(0) * P1[iXYZ];
        //     P02[iXYZ] = (1.0-xsi02.at(0)) * P0[iXYZ] + xsi02.at(0) * P2[iXYZ];
        // }
        
        // Add P0
        _polygon->addPoint(P0);
        // if it doesnt clip the edge across add P01 and P02
        if (xsi12.size() == 0){
            if (isClipping01)
                _polygon->addPoint(P01);
            if (isClipping02)
                _polygon->addPoint(P02);
        } // if it clips once add P01, P12, P02
        else if (xsi12.size() == 1){
            if (isClipping01)
                _polygon->addPoint(P01);
            double P12[3];
            for (int iXYZ = 0; iXYZ < dim; iXYZ++)
                P12[iXYZ] = (1.0-xsi12.at(0)) * P1[iXYZ] + xsi12.at(0) * P2[iXYZ];
            _polygon->addPoint(P12);
            if (isClipping02)
                _polygon->addPoint(P02);
        }
        // if it clips twice
        else if (xsi12.size() == 2){
            if (isClipping01)
                _polygon->addPoint(P01);
            double P12_1[3];
            double P12_2[3];
            // Sort xsi12 so that the order on the edge is P1->P12_1->P12_2->P2
            sort(xsi12.begin(), xsi12.end());
            for (int iXYZ = 0; iXYZ < dim; iXYZ++){
                P12_1[iXYZ] = (1.0-xsi12.at(0)) * P1[iXYZ] + xsi12.at(0) * P2[iXYZ];
                P12_2[iXYZ] = (1.0-xsi12.at(1)) * P1[iXYZ] + xsi12.at(1) * P2[iXYZ];
            }
            _polygon->addPoint(P12_1);
            _polygon->addPoint(P12_2);
            if (isClipping02)
                _polygon->addPoint(P02);
        }
        xsi12.clear();

    } // triaCase == 2: two nodes inside
    else if(triaCase == 2) {

        double* P0 = &_elem[outNodesPos.at(0)*dim];  // outside point
        double* P1 = &_elem[inNodesPos.at(0)*dim];   // next point
        double* P2 = &_elem[inNodesPos.at(1)*dim];   // prev point

        // Find clipping
        vector<double> xsi01;
        vector<double> xsi02;
        bool isClipping01 = findClippingWithMaster(_masterNodeIdx, P0, P1, xsi01);  // P01
        bool isClipping02 = findClippingWithMaster(_masterNodeIdx, P0, P2, xsi02);  // P02

        // Add points to the polygon
        _polygon->addPoint(P1);
        _polygon->addPoint(P2);

        // Compute the global Cartesian coordinates of the clippings if they are valid
        if (isClipping02){
            double P02[3];
            for (int iXYZ = 0; iXYZ < dim; iXYZ++)
                P02[iXYZ] = (1.0-xsi02.at(0)) * P0[iXYZ] + xsi02.at(0) * P2[iXYZ];
            xsi02.clear();
            _polygon->addPoint(P02);
        }        

        if (isClipping01){
            double P01[3];
            for (int iXYZ = 0; iXYZ < dim; iXYZ++)
                P01[iXYZ] = (1.0-xsi01.at(0)) * P0[iXYZ] + xsi01.at(0) * P1[iXYZ];
            xsi01.clear();
            _polygon->addPoint(P01);
        }        

    } // triaCase == 3: 
    else if (triaCase == 3) {

        _polygon->addPoint(&_elem[inNodesPos.at(0)*dim]);
        _polygon->addPoint(&_elem[inNodesPos.at(1)*dim]);
        _polygon->addPoint(&_elem[inNodesPos.at(2)*dim]);

    } // quadCase == 1: one node inside
    else if (quadCase == 1){

        double* P0 = &_elem[inNodesPos.at(0)*dim];

        // Add the inside node to the integration polygon
        _polygon->addPoint(P0);
        // Loop over the outside nodes and find a clipping between inside node and outside node
        for (int iOutNode = 1; iOutNode < _numNodes; iOutNode++){
            double P0n[3];
            double* Pn = &_elem[((inNodesPos.at(0)+iOutNode)%_numNodes)*dim];
            vector<double> xsi0n;
            bool isClipping = findClippingWithMaster(_masterNodeIdx, P0, Pn, xsi0n);
            // Compute the global cartesian coordinates of the clipping location
            for (int iXYZ = 0; iXYZ < dim; iXYZ++)
                P0n[iXYZ] = (1.0-xsi0n.at(0)) * P0[iXYZ] + xsi0n.at(0) * Pn[iXYZ];
            xsi0n.clear();
            // Add the point to the integration polygon
            _polygon->addPoint(P0n);
        }
    } // quadCase == 2: two adjacent nodes inside
    else if (quadCase == 2){

        // The nodes are always ordered in the counter-clockwise order
        double* P0;
        double* P1;
        double* P2;
        double* P3;

        // if first and the last nodes are inside
        if (abs(inNodesPos.at(0) - inNodesPos.at(1)) == 3){
            P0 = &_elem[inNodesPos.at(1)*dim];
            P1 = &_elem[inNodesPos.at(0)*dim];
            P2 = &_elem[((inNodesPos.at(0)+_numNodes+1)%_numNodes)*dim];
            P3 = &_elem[((inNodesPos.at(1)+_numNodes-1)%_numNodes)*dim];
        } // if two consecutive nodes are inside
        else {
            P0 = &_elem[inNodesPos.at(0)*dim];
            P1 = &_elem[inNodesPos.at(1)*dim];
            P2 = &_elem[((inNodesPos.at(1)+_numNodes+1)%_numNodes)*dim];
            P3 = &_elem[((inNodesPos.at(0)+_numNodes-1)%_numNodes)*dim];
        }

        // Find clipping
        vector<double> xsi03;
        vector<double> xsi12;
        bool isClipping03 = findClippingWithMaster(_masterNodeIdx, P0, P3, xsi03);  // P03
        bool isClipping12 = findClippingWithMaster(_masterNodeIdx, P1, P2, xsi12);  // P12

        // Compute the global Cartesian coordinates of the clippings
        double P03[3];
        double P12[3];
        for (int iXYZ = 0; iXYZ < dim; iXYZ++){
            P03[iXYZ] = (1.0-xsi03.at(0)) * P0[iXYZ] + xsi03.at(0) * P3[iXYZ];
            P12[iXYZ] = (1.0-xsi12.at(0)) * P1[iXYZ] + xsi12.at(0) * P2[iXYZ];
        }
        xsi03.clear();
        xsi12.clear();

        // Add the points to the polygon
        _polygon->addPoint(P0);
        _polygon->addPoint(P1);
        _polygon->addPoint(P12);
        _polygon->addPoint(P03);

    } // quadCase == 3: three nodes inside
    else if (quadCase == 3){

        // Add the inside nodes in counter-clockwise order wrt the outside node
        for (int iInNode = 1; iInNode < 4; iInNode++){
            double* Pn = &_elem[((outNodesPos.at(0)+iInNode)%_numNodes)*dim];
            _polygon->addPoint(Pn);
        }

        // Loop over the inside nodes in clockwise order wrt the outside node and add the clippings
        double* P0 = &_elem[outNodesPos.at(0)*dim];
        for (int iInNode = 3; iInNode > 0; iInNode--){
            // Find clipping
            vector<double> xsi0n;
            double* Pn = &_elem[((outNodesPos.at(0)+iInNode)%_numNodes)*dim];
            bool isClipping = findClippingWithMaster(_masterNodeIdx, P0, Pn, xsi0n);
            if (isClipping) {
                // Compute the global Cartesian coordinates of the clipping
                double P0n[3];
                for (int iXYZ = 0; iXYZ < dim; iXYZ++)
                    P0n[iXYZ] = (1.0-xsi0n.at(0)) * P0[iXYZ] + xsi0n.at(0) * Pn[iXYZ];
                // Add the point to the polygon
                _polygon->addPoint(P0n);
            }
        }

    } // quadCase == 4: 
    else if (quadCase == 4) {
        
        _polygon->addPoint(&_elem[inNodesPos.at(0)*dim]);
        _polygon->addPoint(&_elem[inNodesPos.at(1)*dim]);
        _polygon->addPoint(&_elem[inNodesPos.at(2)*dim]);
        _polygon->addPoint(&_elem[inNodesPos.at(3)*dim]);

    } // quadCase == 5: two nodes accross each other inside
    else if (quadCase == 5){
        if (Message::isDebugMode())
            WARNING_BLOCK_OUT("VertexMorphingMapper","clipElementWithFilterRadius", "quadCase == 5 is not implemented yet!");
        return;
    } // unknown case
    else {
        if (Message::isDebugMode()){
            cout << "Tria case: " << triaCase << " Quad case: " << quadCase << endl;
            WARNING_BLOCK_OUT("VertexMorphingMapper","clipElementWithFilterRadius", "Unknown case!");
        }
        return;
    }
    
    inNodesPos.clear();
    outNodesPos.clear();

}

void VertexMorphingMapper::clipMasterElementWithFilterRadius(const int _slaveNodeIdx, const int _numNodes, double* _elem, VertexMorphingMapper::Polygon* _polygon){

    // This function considers the clipping cases one by one instead of trying to clip each edge by default.
    // This way turned out to be more efficient than the latter.

    int dim = 3;

    // Count the influenced nodes by the master node
    int numInNodes = 0;
    std::vector<int> inNodesPos;
    std::vector<int> outNodesPos;
    int nodePos = -1;

    // Find the inside and outside nodes of this element
    for (int iNode = 0; iNode < _numNodes; iNode++){
        nodePos++;
        double dist = EMPIRE::MathLibrary::computeDenseEuclideanNorm(dim, &meshA->nodes[_slaveNodeIdx*dim], &_elem[iNode*dim]);
        if (dist < filterRadius){
            inNodesPos.push_back(nodePos);
            numInNodes++;
        } else
            outNodesPos.push_back(nodePos);
    }

    // triaCase = 1: one node inside
    // triaCase = 2: two nodes inside
    // triaCase = 3: all nodes inside
    // quadCase = 1: one node inside
    // quadCase = 2: two adjacent nodes inside
    // quadCase = 3: three nodes inside
    // quadCase = 4: all nodes inside
    // quadCase = 5: two nodes accross each other inside
    int triaCase = 0;
    int quadCase = 0;
    if (_numNodes == 3) {
        triaCase += numInNodes;
    }
    else if (_numNodes == 4){
        quadCase += numInNodes;
        if (numInNodes == 2 && (abs(inNodesPos[0]-inNodesPos[1]) == 2))
            quadCase += 3;
    } else assert(false);

    // triaCase == 1: one node inside
    if (triaCase == 1){
        // cout << "triaCase == 1" << endl;
        double* P0 = &_elem[inNodesPos.at(0)*dim];  // inside point
        double* P1 = &_elem[outNodesPos.at(0)*dim]; // next point
        double* P2 = &_elem[outNodesPos.at(1)*dim]; // prev point
        // cout << "got points" << endl;
        
        // xsi01 and xsi02 must have only 1 clipping
        vector<double> xsi01;
        vector<double> xsi02;
        // xsi12 might have 1, 2 or no clippings
        vector<double> xsi12;

        // Find clippings
        // cout << "finding clippings" << endl;
        bool isClipping01 = findClippingWithSlave(_slaveNodeIdx, P0, P1, xsi01);  // P01
        bool isClipping02 = findClippingWithSlave(_slaveNodeIdx, P0, P2, xsi02);  // P02
        bool isClipping12 = findClippingWithSlave(_slaveNodeIdx, P1, P2, xsi12);  // P12
        // cout << "found clippings " << xsi01.size() << " " << xsi02.size() << " " << xsi12.size() << endl;

        // Compute the global cartesian coordinates of the clipping points
        double P01[3];
        double P02[3];
        if (isClipping01){
            for (int iXYZ = 0; iXYZ < dim; iXYZ++)
                P01[iXYZ] = (1.0-xsi01.at(0)) * P0[iXYZ] + xsi01.at(0) * P1[iXYZ];
        }
        if (isClipping02){
            for (int iXYZ = 0; iXYZ < dim; iXYZ++)
                P02[iXYZ] = (1.0-xsi02.at(0)) * P0[iXYZ] + xsi02.at(0) * P2[iXYZ];
        }
        xsi01.clear();
        xsi02.clear();

        // for (int iXYZ = 0; iXYZ < dim; iXYZ++){
        //     P01[iXYZ] = (1.0-xsi01.at(0)) * P0[iXYZ] + xsi01.at(0) * P1[iXYZ];
        //     P02[iXYZ] = (1.0-xsi02.at(0)) * P0[iXYZ] + xsi02.at(0) * P2[iXYZ];
        // }
        
        // Add P0
        _polygon->addPoint(P0);
        // if it doesnt clip the edge across add P01 and P02
        if (xsi12.size() == 0){
            if (isClipping01)
                _polygon->addPoint(P01);
            if (isClipping02)
                _polygon->addPoint(P02);
        } // if it clips once add P01, P12, P02
        else if (xsi12.size() == 1){
            if (isClipping01)
                _polygon->addPoint(P01);
            double P12[3];
            for (int iXYZ = 0; iXYZ < dim; iXYZ++)
                P12[iXYZ] = (1.0-xsi12.at(0)) * P1[iXYZ] + xsi12.at(0) * P2[iXYZ];
            _polygon->addPoint(P12);
            if (isClipping02)
                _polygon->addPoint(P02);
        }
        // if it clips twice
        else if (xsi12.size() == 2){
            if (isClipping01)
                _polygon->addPoint(P01);
            double P12_1[3];
            double P12_2[3];
            // Sort xsi12 so that the order on the edge is P1->P12_1->P12_2->P2
            sort(xsi12.begin(), xsi12.end());
            for (int iXYZ = 0; iXYZ < dim; iXYZ++){
                P12_1[iXYZ] = (1.0-xsi12.at(0)) * P1[iXYZ] + xsi12.at(0) * P2[iXYZ];
                P12_2[iXYZ] = (1.0-xsi12.at(1)) * P1[iXYZ] + xsi12.at(1) * P2[iXYZ];
            }
            _polygon->addPoint(P12_1);
            _polygon->addPoint(P12_2);
            if (isClipping02)
                _polygon->addPoint(P02);
        }
        xsi12.clear();

    } // triaCase == 2: two nodes inside
    else if(triaCase == 2) {

        double* P0 = &_elem[outNodesPos.at(0)*dim];  // outside point
        double* P1 = &_elem[inNodesPos.at(0)*dim];   // next point
        double* P2 = &_elem[inNodesPos.at(1)*dim];   // prev point

        // Find clipping
        vector<double> xsi01;
        vector<double> xsi02;
        bool isClipping01 = findClippingWithSlave(_slaveNodeIdx, P0, P1, xsi01);  // P01
        bool isClipping02 = findClippingWithSlave(_slaveNodeIdx, P0, P2, xsi02);  // P02

        // Add points to the polygon
        _polygon->addPoint(P1);
        _polygon->addPoint(P2);

        // Compute the global Cartesian coordinates of the clippings if they are valid
        if (isClipping02){
            double P02[3];
            for (int iXYZ = 0; iXYZ < dim; iXYZ++)
                P02[iXYZ] = (1.0-xsi02.at(0)) * P0[iXYZ] + xsi02.at(0) * P2[iXYZ];
            xsi02.clear();
            _polygon->addPoint(P02);
        }        

        if (isClipping01){
            double P01[3];
            for (int iXYZ = 0; iXYZ < dim; iXYZ++)
                P01[iXYZ] = (1.0-xsi01.at(0)) * P0[iXYZ] + xsi01.at(0) * P1[iXYZ];
            xsi01.clear();
            _polygon->addPoint(P01);
        }        

    } // triaCase == 3: 
    else if (triaCase == 3) {

        _polygon->addPoint(&_elem[inNodesPos.at(0)*dim]);
        _polygon->addPoint(&_elem[inNodesPos.at(1)*dim]);
        _polygon->addPoint(&_elem[inNodesPos.at(2)*dim]);

    } // quadCase == 1: one node inside
    else if (quadCase == 1){

        double* P0 = &_elem[inNodesPos.at(0)*dim];

        // Add the inside node to the integration polygon
        _polygon->addPoint(P0);
        // Loop over the outside nodes and find a clipping between inside node and outside node
        for (int iOutNode = 1; iOutNode < _numNodes; iOutNode++){
            double P0n[3];
            double* Pn = &_elem[((inNodesPos.at(0)+iOutNode)%_numNodes)*dim];
            vector<double> xsi0n;
            bool isClipping = findClippingWithSlave(_slaveNodeIdx, P0, Pn, xsi0n);
            // Compute the global cartesian coordinates of the clipping location
            for (int iXYZ = 0; iXYZ < dim; iXYZ++)
                P0n[iXYZ] = (1.0-xsi0n.at(0)) * P0[iXYZ] + xsi0n.at(0) * Pn[iXYZ];
            xsi0n.clear();
            // Add the point to the integration polygon
            _polygon->addPoint(P0n);
        }
    } // quadCase == 2: two adjacent nodes inside
    else if (quadCase == 2){

        // The nodes are always ordered in the counter-clockwise order
        double* P0;
        double* P1;
        double* P2;
        double* P3;

        // if first and the last nodes are inside
        if (abs(inNodesPos.at(0) - inNodesPos.at(1)) == 3){
            P0 = &_elem[inNodesPos.at(1)*dim];
            P1 = &_elem[inNodesPos.at(0)*dim];
            P2 = &_elem[((inNodesPos.at(0)+_numNodes+1)%_numNodes)*dim];
            P3 = &_elem[((inNodesPos.at(1)+_numNodes-1)%_numNodes)*dim];
        } // if two consecutive nodes are inside
        else {
            P0 = &_elem[inNodesPos.at(0)*dim];
            P1 = &_elem[inNodesPos.at(1)*dim];
            P2 = &_elem[((inNodesPos.at(1)+_numNodes+1)%_numNodes)*dim];
            P3 = &_elem[((inNodesPos.at(0)+_numNodes-1)%_numNodes)*dim];
        }

        // Find clipping
        vector<double> xsi03;
        vector<double> xsi12;
        bool isClipping03 = findClippingWithSlave(_slaveNodeIdx, P0, P3, xsi03);  // P03
        bool isClipping12 = findClippingWithSlave(_slaveNodeIdx, P1, P2, xsi12);  // P12

        // Compute the global Cartesian coordinates of the clippings
        double P03[3];
        double P12[3];
        for (int iXYZ = 0; iXYZ < dim; iXYZ++){
            P03[iXYZ] = (1.0-xsi03.at(0)) * P0[iXYZ] + xsi03.at(0) * P3[iXYZ];
            P12[iXYZ] = (1.0-xsi12.at(0)) * P1[iXYZ] + xsi12.at(0) * P2[iXYZ];
        }
        xsi03.clear();
        xsi12.clear();

        // Add the points to the polygon
        _polygon->addPoint(P0);
        _polygon->addPoint(P1);
        _polygon->addPoint(P12);
        _polygon->addPoint(P03);

    } // quadCase == 3: three nodes inside
    else if (quadCase == 3){

        // Add the inside nodes in counter-clockwise order wrt the outside node
        for (int iInNode = 1; iInNode < 4; iInNode++){
            double* Pn = &_elem[((outNodesPos.at(0)+iInNode)%_numNodes)*dim];
            _polygon->addPoint(Pn);
        }

        // Loop over the inside nodes in clockwise order wrt the outside node and add the clippings
        double* P0 = &_elem[outNodesPos.at(0)*dim];
        for (int iInNode = 3; iInNode > 0; iInNode--){
            // Find clipping
            vector<double> xsi0n;
            double* Pn = &_elem[((outNodesPos.at(0)+iInNode)%_numNodes)*dim];
            bool isClipping = findClippingWithSlave(_slaveNodeIdx, P0, Pn, xsi0n);
            if (isClipping) {
                // Compute the global Cartesian coordinates of the clipping
                double P0n[3];
                for (int iXYZ = 0; iXYZ < dim; iXYZ++)
                    P0n[iXYZ] = (1.0-xsi0n.at(0)) * P0[iXYZ] + xsi0n.at(0) * Pn[iXYZ];
                // Add the point to the polygon
                _polygon->addPoint(P0n);
            }
        }

    } // quadCase == 4: 
    else if (quadCase == 4) {
        
        _polygon->addPoint(&_elem[inNodesPos.at(0)*dim]);
        _polygon->addPoint(&_elem[inNodesPos.at(1)*dim]);
        _polygon->addPoint(&_elem[inNodesPos.at(2)*dim]);
        _polygon->addPoint(&_elem[inNodesPos.at(3)*dim]);

    } // quadCase == 5: two nodes accross each other inside
    else if (quadCase == 5){
        if (Message::isDebugMode())
            WARNING_BLOCK_OUT("VertexMorphingMapper","clipElementWithFilterRadius", "quadCase == 5 is not implemented yet!");
        return;
    } // unknown case
    else {
        if (Message::isDebugMode()){
            cout << "Tria case: " << triaCase << " Quad case: " << quadCase << endl;
            WARNING_BLOCK_OUT("VertexMorphingMapper","clipElementWithFilterRadius", "Unknown case!");
        }
        return;
    }
    
    inNodesPos.clear();
    outNodesPos.clear();

}

bool VertexMorphingMapper::findClippingWithMaster(const int _masterNodeIdx, double* _P0, double* _Pn, vector<double>& _xsi){

    int dim = 3;

    // Get the node
    double* Pc = &meshB->nodes[_masterNodeIdx*dim];

    // Some precomputations
    double P0P0 = EMPIRE::MathLibrary::computeDenseDotProduct(dim, _P0, _P0);
    double PnPn = EMPIRE::MathLibrary::computeDenseDotProduct(dim, _Pn, _Pn);
    double PcPc = EMPIRE::MathLibrary::computeDenseDotProduct(dim, Pc, Pc);
    double P0Pn = EMPIRE::MathLibrary::computeDenseDotProduct(dim, _P0, _Pn);
    double PcP0 = EMPIRE::MathLibrary::computeDenseDotProduct(dim, Pc, _P0);
    double PcPn = EMPIRE::MathLibrary::computeDenseDotProduct(dim, Pc, _Pn);

    // compute the line parameter: P0n(xsi) = (1-xsi)*P0 + xsi*Pn
    // it is a root finding problem of a 2nd degree polynomial
    double xsi_pre;
    double a_xsi = P0P0 - 2.0*P0Pn + PnPn;
    double b_xsi = 2.0*(-P0P0 + P0Pn + PcP0 - PcPn);
    double c_xsi = P0P0 - 2.0*PcP0 + PcPc - filterRadius*filterRadius;
    double discriminant = pow(b_xsi,2)-4.0*a_xsi*c_xsi;
    // it cannot be negative
    if (discriminant < 0.0 )
        return false;

    // If there are multiple roots
    if (discriminant > 0.0) {

        // Check the validity of the roots and add them if they are valid
        xsi_pre = (-b_xsi + sqrt(discriminant)) / (2.0*a_xsi);
        if (clampXsi(xsi_pre)) {
            _xsi.push_back(xsi_pre);
        }
        xsi_pre = (-b_xsi - sqrt(discriminant)) / (2.0*a_xsi);
        if (clampXsi(xsi_pre)) {
            _xsi.push_back(xsi_pre);
        }

    } // if there is single root
    else if (discriminant == 0.0) {
        xsi_pre = -b_xsi / (2.0*a_xsi);
        if (clampXsi(xsi_pre)) {
            _xsi.push_back(xsi_pre);
        }
    }

    if (_xsi.size() > 0)
        return true;
    else
        return false;
}

bool VertexMorphingMapper::findClippingWithSlave(const int _slaveNodeIdx, double* _P0, double* _Pn, vector<double>& _xsi){

    int dim = 3;

    // Get the node
    double* Pc = &meshA->nodes[_slaveNodeIdx*dim];

    // Some precomputations
    double P0P0 = EMPIRE::MathLibrary::computeDenseDotProduct(dim, _P0, _P0);
    double PnPn = EMPIRE::MathLibrary::computeDenseDotProduct(dim, _Pn, _Pn);
    double PcPc = EMPIRE::MathLibrary::computeDenseDotProduct(dim, Pc, Pc);
    double P0Pn = EMPIRE::MathLibrary::computeDenseDotProduct(dim, _P0, _Pn);
    double PcP0 = EMPIRE::MathLibrary::computeDenseDotProduct(dim, Pc, _P0);
    double PcPn = EMPIRE::MathLibrary::computeDenseDotProduct(dim, Pc, _Pn);

    // compute the line parameter: P0n(xsi) = (1-xsi)*P0 + xsi*Pn
    // it is a root finding problem of a 2nd degree polynomial
    double xsi_pre;
    double a_xsi = P0P0 - 2.0*P0Pn + PnPn;
    double b_xsi = 2.0*(-P0P0 + P0Pn + PcP0 - PcPn);
    double c_xsi = P0P0 - 2.0*PcP0 + PcPc - filterRadius*filterRadius;
    double discriminant = pow(b_xsi,2)-4.0*a_xsi*c_xsi;
    // it cannot be negative
    if (discriminant < 0.0 )
        return false;

    // If there are multiple roots
    if (discriminant > 0.0) {

        // Check the validity of the roots and add them if they are valid
        xsi_pre = (-b_xsi + sqrt(discriminant)) / (2.0*a_xsi);
        if (clampXsi(xsi_pre)) {
            _xsi.push_back(xsi_pre);
        }
        xsi_pre = (-b_xsi - sqrt(discriminant)) / (2.0*a_xsi);
        if (clampXsi(xsi_pre)) {
            _xsi.push_back(xsi_pre);
        }

    } // if there is single root
    else if (discriminant == 0.0) {
        xsi_pre = -b_xsi / (2.0*a_xsi);
        if (clampXsi(xsi_pre)) {
            _xsi.push_back(xsi_pre);
        }
    }

    if (_xsi.size() > 0)
        return true;
    else
        return false;
}

bool VertexMorphingMapper::clampXsi(double& _xsi){

    // if the parameter is outside the admissible range then invalid -> returns false

    // if the parameter is within the tolerances
    // clamp it to the range 0 <= xsi <= 1
    if (_xsi > -EPS_XSI && _xsi < 0.0) {
        _xsi = 0.0;
    }
    if (_xsi > 1.0 && _xsi < 1.0 + EPS_XSI) {
        _xsi = 1.0;
    }
    if (_xsi>=0.0 && _xsi <= 1.0)
        return true;

    return false;

}

void VertexMorphingMapper::adjustFilterFunctions_BA(){

    for (int iRow = 0; iRow < C_BA->getNumberOfRows(); iRow++)
        (*C_BA).multiplyRowWith(iRow, 1.0/masterFilterFunctionIntegrationOnSlave[iRow]);

    // filterFunctionCorrectionFactors = new double[meshB->numNodes];
    // double* RHS = new double[meshB->numNodes];
    // for (int i = 0; i<meshB->numNodes; i++)
    //     RHS[i] = 1.0;

    // cout << "Solving for filter function correction factors" << endl;
    // (*C_corr).solve(filterFunctionCorrectionFactors, RHS);
    // cout << "Finished!" << endl;

    // cout << "Adjusting filter sizes" << endl;
    // for (int iRow = 0; iRow < C_BA->getNumberOfRows(); iRow++)
    //     (*C_BA).multiplyRowWith(iRow, filterFunctionCorrectionFactors[iRow]);
    // cout << "Finished!" << endl;

    // delete[] filterFunctionCorrectionFactors;
    // delete[] RHS;

}

void VertexMorphingMapper::adjustFilterFunctions_AB(){

    for (int iRow = 0; iRow < C_AB->getNumberOfRows(); iRow++)
        (*C_AB).multiplyRowWith(iRow, 1.0/slaveFilterFunctionIntegrationOnMaster[iRow]);

}

void VertexMorphingMapper::consistentMapping(const double *slaveField, double *masterField) {

    //    INFO_OUT() << "VertexMorphingMapper::consistentMapping ->  Norm of the Slave field :: "  << EMPIRE::MathLibrary::computeVectorLength(slaveField) << endl;
    //    INFO_OUT() << "VertexMorphingMapper::consistentMapping ->  Norm of the Master field :: " << EMPIRE::MathLibrary::computeVectorLength(masterField) << endl;

    double *slaveFieldCopy = new double[meshA->numNodes];
    for (int i = 0; i < meshA->numNodes; i++)
        slaveFieldCopy[i] = slaveField[i];

    // matrix vector product (W_tmp = C_BA * W_A)
    (*C_BA).mulitplyVec(false,slaveFieldCopy,masterField,meshB->numNodes);

    delete[] slaveFieldCopy;

}

void VertexMorphingMapper::conservativeMapping(const double *masterField, double *slaveField) {

    //    INFO_OUT() << "VertexMorphingMapper::consistentMapping ->  Norm of the Slave field :: "  << EMPIRE::MathLibrary::computeVectorLength(slaveField) << endl;
    //    INFO_OUT() << "VertexMorphingMapper::consistentMapping ->  Norm of the Master field :: " << EMPIRE::MathLibrary::computeVectorLength(masterField) << endl;

    double *masterFieldCopy = new double[meshB->numNodes];
    for (int i = 0; i < meshB->numNodes; i++)
        masterFieldCopy[i] = masterField[i];

    // // matrix vector product (W_tmp = C_BA * W_A)
    // (*C_AB).mulitplyVec(false,masterFieldCopy,slaveField,meshA->numNodes);

    (*C_BA).transposeMulitplyVec(masterFieldCopy, slaveField, meshA->numNodes);


    delete[] masterFieldCopy;

}

void VertexMorphingMapper::computeErrorsConsistentMapping(const double *slaveField, const double *masterField){

    ERROR_BLOCK_OUT("VertexMorphingMapper","computeErrorsConsistentMapping", "Not Implemented!");
    assert(false);

}

// FilterFunction classes
double VertexMorphingMapper::HatFilterFunction::computeFunction(double* _supportCenter, double* _globalCoor){

    // compute distance to the support center
    double dist = EMPIRE::MathLibrary::computeDenseEuclideanNorm(3, _supportCenter, _globalCoor);

    // if the distance is smaller equal to filterRadius
    if (dist <= filterRadius) // return function value
        return 1.0 - dist/filterRadius;
    else // the filter function is 0 out of its effective radius
        return 0.0;

}

double VertexMorphingMapper::GaussianFilterFunction::computeFunction(double* _supportCenter, double* _globalCoor){

    // compute distance to the support center
    double dist = EMPIRE::MathLibrary::computeDenseEuclideanNorm(3, _supportCenter, _globalCoor);

    // if the distance is smaller equal to filterRadius
    if (dist <= filterRadius) // return function value
        return exp(-pow(dist,2) / (2 * pow(filterRadius,2) / 9.0));
    else // the filter function is 0 out of its effective radius
        return 0.0;
   
}

// Integrand classes
VertexMorphingMapper::FilterFunctionProduct::FilterFunctionProduct(int _numNodes, double* _elem, double* _controlNode)
{

    int dim=3;

    gaussPoints = NULL;
    functionProducts = NULL;

    numNodes = _numNodes;
    elem = new double[_numNodes*dim];
    for (int iNode = 0; iNode < numNodes; iNode++){
        for (int iXYZ = 0; iXYZ < dim; iXYZ++)
            elem[iNode*dim + iXYZ] = _elem[iNode*dim + iXYZ];
    }

    controlNode = new double[dim];
    for (int iXYZ = 0; iXYZ < dim; iXYZ++)
        controlNode[iXYZ] = _controlNode[iXYZ];

}

VertexMorphingMapper::FilterFunctionProduct::~FilterFunctionProduct(){

    delete[] controlNode;
    delete[] elem;

    if (gaussPoints != NULL)
        delete[] gaussPoints;

    if (functionProducts != NULL) {
        for (int i=0; i < numNodes; i++)
            delete[] functionProducts[i];
        delete[] functionProducts;
    }

}

void VertexMorphingMapper::FilterFunctionProduct::setGaussPoints(const double *_gaussPoints, int _numGaussPoints) {
    assert(gaussPoints == NULL);
    int dim = 3;
    numGaussPoints = _numGaussPoints;
    gaussPoints = new double[numGaussPoints * dim];
    for (int i = 0; i < numGaussPoints * dim; i++)
        gaussPoints[i] = _gaussPoints[i];
}

void VertexMorphingMapper::FilterFunctionProduct::computeFunctionProducts(AbstractFilterFunction *_filterFunction){
    assert(functionProducts == NULL);
    assert(gaussPoints != NULL);

    int dim = 3;

    // initiate function products container
    functionProducts = new double*[numNodes];
    for (int i = 0; i < numNodes; i++)
        functionProducts[i] = new double[numGaussPoints];

    // loop over gauss points
    for (int iGP = 0; iGP < numGaussPoints; iGP++) {
        double *gaussPoint = &gaussPoints[iGP * dim];

        // *. compute shape functions on element (on the Gauss point)
        double shapeFuncValues[numNodes];
        double normal[3];
        bool inside;
        if (numNodes == 3) {
            EMPIRE::MathLibrary::computeNormalOfTriangle(elem, true, normal);
            int planeToProject = EMPIRE::MathLibrary::computePlaneToProject(normal);
            inside = EMPIRE::MathLibrary::computeLocalCoorInTriangle(elem, planeToProject,
                                                                     gaussPoint, shapeFuncValues);
        } else if (numNodes == 4) {
            double localCoords[2];
            EMPIRE::MathLibrary::computeNormalOfQuad(elem, true, normal);
            int planeToProject = EMPIRE::MathLibrary::computePlaneToProject(normal);
            inside = EMPIRE::MathLibrary::computeLocalCoorInQuad(elem, planeToProject,
                                                                 gaussPoint, localCoords);
            EMPIRE::MathLibrary::computeShapeFuncOfQuad(localCoords, shapeFuncValues);
        } else assert(false);

        //        // debug
        //        if (!inside){
        //            cout<<"Error in computing local coordinates in tria master element"<<endl;
        //            cout<<"GP coordinates: "<<gaussPoint[0]<<"\t"<<gaussPoint[1]<<"\t"<<gaussPoint[2]<<endl;
        //            cout<<"Element nodes:"<<endl;
        //            for(int ctr=0;ctr<numNodes;ctr++){
        //                cout<<ctr+1<<": "<<elem[ctr*dim]<<"\t"<<elem[ctr*dim+1]<<"\t"<<elem[ctr*dim+2]<<endl;
        //            }
        //        }
        //        // debug end
        //        assert(inside);

        // *. compute shape function products (on the Gauss point)
        for (int iNode = 0; iNode < numNodes; iNode++)
            functionProducts[iNode][iGP] = shapeFuncValues[iNode] * _filterFunction->computeFunction(controlNode, gaussPoint);
    }
}

double VertexMorphingMapper::FilterFunctionProduct::operator ()(double *gaussPoint) {
    assert(gaussPoints!=NULL);
    assert(functionProducts!=NULL);
    for (int i = 0; i < numGaussPoints; i++) {
        if (gaussPoints[i * 3 + 0] == gaussPoint[0] && gaussPoints[i * 3 + 1] == gaussPoint[1]
                && gaussPoints[i * 3 + 2] == gaussPoint[2]) {
            return functionProducts[funcID][i];
        }
    }
    // should not come here
    assert(false);
    return 0.0;
}

void VertexMorphingMapper::FilterFunctionProduct::setFunctionID(int _funcID) {
    funcID = _funcID;
}

VertexMorphingMapper::Polygon::~Polygon(){

    if (triangles.size() > 0){
        for (vector<double*>::iterator itTriangle = triangles.begin(); itTriangle != triangles.end(); itTriangle++)
            delete[] *itTriangle;
        triangles.clear();
    }

    if (points.size() > 0){
        for (vector<double*>::iterator itPoint = points.begin(); itPoint != points.end(); itPoint++)
            delete[] *itPoint;
        points.clear();
    }

}

void VertexMorphingMapper::Polygon::addPoint(double* _point){
    
    int dim = 3;
    double* Padd = new double[dim];
    for (int iXYZ = 0; iXYZ < dim; iXYZ++)
        Padd[iXYZ] = _point[iXYZ];

    if ( points.size() > 0 ){
        double pointDist = EMPIRE::MathLibrary::computeDenseEuclideanNorm(3, Padd, points.back());
        if (pointDist > EPS_PolygonPoint)
            points.push_back(Padd);
    } else
        points.push_back(Padd);

}

void VertexMorphingMapper::Polygon::printPolygon(){

    cout << "Polygon Points:" << endl;
    for (int iPoint=0; iPoint < points.size(); iPoint++)
        cout << points.at(iPoint)[0] << " " << points.at(iPoint)[1] << " " << points.at(iPoint)[2] << endl;

}

void VertexMorphingMapper::Polygon::triangulate(){

    if (Message::isDebugMode()){
        if (points.size() == 0)
            WARNING_BLOCK_OUT("VertexMorphingMapper::Polygon","triangulate()","Polygon does not contain points. Skipping triangulation!");
        else if (points.size() == 1)
            WARNING_BLOCK_OUT("VertexMorphingMapper::Polygon","triangulate()","Polygon contains only one point. Skipping triangulation!");
        else if (points.size() == 2)
            WARNING_BLOCK_OUT("VertexMorphingMapper::Polygon","triangulate()","Polygon contains only two points. Skipping triangulation!");
    }

    if (points.size() < 3)
        return;

    int dim = 3;
    int numNodes = 3;

    for (int iTriangle = 0; iTriangle < points.size()-2; iTriangle++){
        double* tmpTriangle = new double[numNodes*dim];
        for (int iNode = 0; iNode < numNodes; iNode++){
            for (int iXYZ = 0; iXYZ < dim; iXYZ++){
                tmpTriangle[iNode*dim+iXYZ] = points.at(iTriangle+iNode)[iXYZ];
            }
        }
        triangles.push_back(tmpTriangle);
    }
}

void VertexMorphingMapper::writeCartesianPolygons(string _fileName, std::map<int, std::vector<std::vector<double*> > >& _polygons){
    
    
    for(map<int, std::vector<std::vector<double*> > >::iterator itControlNode=_polygons.begin(); itControlNode!=_polygons.end(); itControlNode++) {
        ofstream out;
        string iControlNode_string;
        stringstream tmp_ss;
        tmp_ss << itControlNode->first;
        string outName = name + "_" +_fileName + tmp_ss.str() + ".vtk";
        out.open(outName.c_str(), ofstream::out);
        out << "# vtk DataFile Version 2.0\n";
        out << "clipped FE elements\n";
        out << "ASCII\nDATASET POLYDATA\n";
        // 	int controlNodeIdx = itControlNode->first();
        string points, pointsHeader, polygons, polygonsHeader, nodeColor, patchColorHeader;
        int pointsNumber=0, polygonsNumber=0, polygonsEntriesNumber=0;
        for (std::vector<std::vector<double*> >::iterator itPolygon=itControlNode->second.begin(); itPolygon!=itControlNode->second.end(); itPolygon++){
            int numEdge = 0;
            for (std::vector<double*>::iterator itNode=itPolygon->begin(); itNode!=itPolygon->end(); itNode++){
                stringstream pointStream;
                // 		cout << (*itNode)[0] << " " << (*itNode)[1] << " " << (*itNode)[2] << endl;
                pointStream << (*itNode)[0];
                pointStream << " " << (*itNode)[1];
                pointStream << " " << (*itNode)[2];
                pointStream << "\n";
                points += pointStream.str();
                pointsNumber++;
                numEdge++;
            }
            stringstream colorStream;
            /// Concatenate new polygon color
            colorStream << itControlNode->first << "\n";
            nodeColor += colorStream.str();
            polygonsNumber++;
            polygonsEntriesNumber += numEdge + 1;
            stringstream polygonStream;
            polygonStream << numEdge;
            for(int i=numEdge;i>0;i--) {
                polygonStream << " " << pointsNumber - i;
            }
            polygonStream << "\n";
            polygons += polygonStream.str();
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
        header << "CELL_DATA " << polygonsNumber << "\nSCALARS node_belonging int 1\nLOOKUP_TABLE default\n";
        patchColorHeader = header.str();
        out << pointsHeader << points;
        out << polygonsHeader << polygons;
        out << patchColorHeader << nodeColor;
        out.close();
    }

}

} /* namespace EMPIRE */

