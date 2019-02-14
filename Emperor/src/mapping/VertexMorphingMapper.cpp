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

const int VertexMorphingMapper::numGPsMassMatrixTri = 6;
const int VertexMorphingMapper::numGPsMassMatrixQuad = 4;
const int VertexMorphingMapper::numGPsOnClipTri = 6;
const int VertexMorphingMapper::numGPsOnClipQuad = 12;

VertexMorphingMapper::VertexMorphingMapper(std::string _name, FEMesh *_mesh, EMPIRE_VMM_FilterType _filterType, double _filterRadius):
    name(_name), filterType(_filterType), filterRadius(_filterRadius) {

    // Check input
    assert(_mesh != NULL);

    bool isMeshFEM = (_mesh->type == EMPIRE_Mesh_FEMesh || _mesh->type == EMPIRE_Mesh_copyFEMesh);

    if (isMeshFEM) {
        if (dynamic_cast<FEMesh *>(_mesh)->triangulate() == NULL)
            mesh = dynamic_cast<FEMesh *>(_mesh);
        else
            mesh = dynamic_cast<FEMesh *>(_mesh)->triangulate();
    } else {
        ERROR_BLOCK_OUT("VertexMorphingMapper","VertexMorphingMapper","Wrong type of mesh!");
    }

    // Assign the mapper type
    mapperType = EMPIRE_VertexMorphingMapper;

}

void VertexMorphingMapper::initTables() {
    // using the map to store the nodeNumbers
    // but here the "key" is the node number, and the value is the position in nodeNumbers
    // the map is sorted automatically, so it is efficient for searching

    // 1. compute directElemTable

    const int numElems = mesh->numElems;
    const int numNodes = mesh->numNodes;

    directElemTable = new vector<int>*[numElems];
    for (int i = 0; i < numElems; i++)
        directElemTable[i] = new vector<int>;
    map<int, int> *nodesMap = new map<int, int>();
    for (int i = 0; i < numNodes; i++)
        nodesMap->insert(nodesMap->end(), pair<int, int>(mesh->nodeIDs[i], i));
    int count = 0;
    for (int i = 0; i < numElems; i++) {
        const int numNodesElem = mesh->numNodesPerElem[i];
        for (int j = 0; j < numNodesElem; j++) {
            if(nodesMap->find(mesh->elems[count + j])== nodesMap->end() ){
                ERROR_OUT()<< "Slave Node Label " << mesh->elems[count + j] << " is not part of slaveNodesMap." << endl;
            }
            directElemTable[i]->push_back(nodesMap->at(mesh->elems[count + j]));
        }
        count += numNodesElem;
    }
    delete nodesMap;

    // 2. compute slaveNodeToElemTable
    nodeToElemTable = new vector<int>*[numNodes];
    for (int i = 0; i < numNodes; i++)
        nodeToElemTable[i] = new vector<int>;
    for (int i = 0; i < numElems; i++) {
        const int numNodesElem = mesh->numNodesPerElem[i];
        for (int j = 0; j < numNodesElem; j++) {
            int nodePos = directElemTable[i]->at(j);
            nodeToElemTable[nodePos]->push_back(i);
        }
    }

}

void VertexMorphingMapper::initANNTree() {
#ifdef ANN
    ANNNodes = new double*[mesh->numNodes]; // ANN uses 2D array
    for (int i = 0; i < mesh->numNodes; i++) {
        ANNNodes[i] = new double[3];
        for (int j = 0; j<3; j++)
            ANNNodes[i][j] = mesh->nodes[i * 3 + j];
    }
    nodesTree = new ANNkd_tree(ANNNodes, mesh->numNodes, 3);
#endif
#ifdef FLANN
    FLANNNodes = new flann::Matrix<double>(const_cast<double*>(mesh->nodes), mesh->numNodes, 3);
    FLANNkd_tree = new flann::Index<flann::L2<double> >(*FLANNNodes, flann::KDTreeSingleIndexParams(1));
    FLANNkd_tree->buildIndex(); // Build binary tree for searching
#endif
}

void VertexMorphingMapper::deleteANNTree() {
#ifdef ANN
    delete[] ANNNodes;
    delete nodesTree;
    annClose();
#endif
#ifdef FLANN
    delete FLANNNodes;
    delete FLANNkd_tree;
#endif
}

void VertexMorphingMapper::findCandidates(double* controlNode, std::set<int> *infNodeIdxs, std::set<int> *infElemIdxs)
{
    // Use fixed radius search on end points of the element,
    // all elements containing these points are the overlapped candidates
    // 1. find all the influenced elements

    // OpenMP parallelize this loop
    // Ann is not thread safe

    // 1.1 Find the node indices in the filter radius
    int numInfNodes = 0;
    int *infNodeIndices;
#ifdef ANN
    numInfNodes = nodesTree->annkFRSearch(controlNode, filterRadius, 0); // get the number of neighbors in a radius
    infNodeIndices = new int[numInfNodes];
    double *dummy = new double[numInfNodes];
    nodesTree->annkFRSearch(&controlNode, filterRadius, numInfNodes, infNodeIndices, dummy);// get the real neighbors (ANN uses the square of the radius)
    delete dummy;
#endif

#ifdef FLANN
    flann::Matrix<double> controlNodeCopyFlann(controlNode, 1, 3);
    vector<vector<int> > indices_tmp;
    vector<vector<double> > dists_tmp;
    FLANNkd_tree->radiusSearch(controlNodeCopyFlann, indices_tmp, dists_tmp, filterRadius, flann::SearchParams(1));
    numInfNodes = indices_tmp[0].size();
    infNodeIndices = new int[indices_tmp[0].size()];
    for (int j=0; j<indices_tmp[0].size(); j++)
        infNodeIndices[j] = indices_tmp[0][j];
#endif

    // 1.2 Fill up the influenced nodes and the influenced elements sets
    for (int i = 0; i < numInfNodes; i++) {

        // Influenced nodes
        infNodeIdxs->insert(infNodeIndices[i]);
        vector<int> * vecInfElemIndices = nodeToElemTable[infNodeIndices[i]];

        // Influenced elements
        for (vector<int>::iterator it = vecInfElemIndices->begin(); it < vecInfElemIndices->end(); it++)
            infElemIdxs->insert(*it);

    }
    delete[] infNodeIndices;
}

void VertexMorphingMapper::buildCouplingMatrices(){

    const int numNodes = mesh->numNodes;

    // Loop over nodes to find influenced nodes and elements
    for (int iNode = 0; iNode < numNodes; iNode++) {

        // Find candidates
        double* controlNodePtr = &mesh->nodes[iNode*3];
        std::set<int> tmp_infNodes;
        std::set<int> tmp_infElems;
        findCandidates(controlNodePtr, &tmp_infNodes, &tmp_infElems);

        // Loop over the influenced nodes
        for (std::set<int>::iterator itNode = tmp_infNodes.begin(); itNode != tmp_infNodes.end(); itNode++){

            // Loop over the elements that contain the influenced node
            for (std::set<int>::iterator itElem = tmp_infElems.begin(); itElem != tmp_infElems.end(); itElem++){

                // Loop over the element's nodes and check if they are all inside the filter radius
                for (std::vector<int>::iterator itElemNode = directElemTable[*itElem]->begin(); itElemNode != directElemTable[*itElem]->end(); itElemNode++){

                    bool allInfluenced = tmp_infNodes.find(*itElemNode) != tmp_infNodes.end();

                }
            }
        }
    }
}

}
