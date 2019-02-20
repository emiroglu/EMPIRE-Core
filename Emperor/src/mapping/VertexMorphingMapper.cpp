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

const int VertexMorphingMapper::numGPsOnTri = 12;
const int VertexMorphingMapper::numGPsMassMatrixTri = 12;

VertexMorphingMapper::VertexMorphingMapper(std::string _name, AbstractMesh *_meshA, AbstractMesh *_meshB, EMPIRE_VMM_FilterType _filterType, double _filterRadius, bool _mortar):
    name(_name), filterRadius(_filterRadius), mortar(_mortar){

    // Check input
    assert(_meshA != NULL);
    assert(_meshB != NULL);

    // Assign meshes
    bool isMeshAFEM = (_meshA->type == EMPIRE_Mesh_FEMesh || _meshA->type == EMPIRE_Mesh_copyFEMesh);
    bool isMeshBFEM = (_meshB->type == EMPIRE_Mesh_FEMesh || _meshB->type == EMPIRE_Mesh_copyFEMesh);

    if (isMeshAFEM && isMeshBFEM){
        if (dynamic_cast<FEMesh *>(_meshA)->triangulate() == NULL){
            ERROR_BLOCK_OUT("VertexMorphingMapper","VertexMorphingMapper","Mesh A has to be triangulated!");
        }
        else
            meshA = dynamic_cast<FEMesh *>(_meshA)->triangulate();

        if (dynamic_cast<FEMesh *>(_meshB)->triangulate() == NULL){
            ERROR_BLOCK_OUT("VertexMorphingMapper","VertexMorphingMapper","Mesh B has to be triangulated!");
        }
        else
            meshB = dynamic_cast<FEMesh *>(_meshB)->triangulate();
    } else
        ERROR_BLOCK_OUT("VertexMorphingMapper","VertexMorphingMapper","Wrong type of mesh!");

    // Assign the mapper type
    mapperType = EMPIRE_VertexMorphingMapper;

    /// Initializing the sparse matrices
    /// This is symmetric square matrix of size "masterNumNodes". For first attempt
    /// we store it as full matrix
    if (mortar)
        C_BB = new MathLibrary::SparseMatrix<double>(meshB->numNodes, false);

    /// This is a rectangular matrix of size
    /// masterNumNodes X slaveNumNodes
    C_BA = new MathLibrary::SparseMatrix<double>((size_t)(meshB->numNodes),(size_t)(meshA->numNodes));

    // Initialize data that could be used later
    initTables();
    initANNTree();

    // Set the filter function pointer
    if (_filterType == EMPIRE_VMM_HatFilter)
        filterFunction = new HatFilterFunction(filterRadius);
    else if (_filterType == EMPIRE_VMM_GaussianFilter)
        filterFunction = new GaussianFilterFunction(filterRadius);
    else
        assert(false); // shouldnt come here

}

VertexMorphingMapper::~VertexMorphingMapper(){

    if (mortar)
        delete C_BB;
    delete C_BA;
    delete filterFunction;

#ifdef ANN
    for (int i = 0; i < slaveNumNodes; i++) {
        delete[] ANNSlaveNodes[i];
    }
    delete[] ANNSlaveNodes;
#endif

}

void VertexMorphingMapper::deleteTables() {
    for (int i = 0; i < meshA->numElems; i++)
        delete slaveDirectElemTable[i];
    delete[] slaveDirectElemTable;

    for (int i = 0; i < meshB->numElems; i++)
        delete masterDirectElemTable[i];
    delete[] masterDirectElemTable;

    for (int i = 0; i < meshA->numNodes; i++)
        delete slaveNodeToElemTable[i];
    delete[] slaveNodeToElemTable;

    for (int i = 0; i < meshB->numNodes; i++)
        delete masterNodeToElemTable[i];
    delete[] masterNodeToElemTable;

}

void VertexMorphingMapper::deleteANNTree() {
#ifdef ANN
    delete[] ANNSlaveNodes;
    delete slaveNodesTree;
    annClose();
#endif
#ifdef FLANN
    delete FLANNSlaveNodes;
    delete FLANNkd_tree;
#endif
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
                if(slaveNodesMap->find(meshA->elems[count + j])== slaveNodesMap->end()){
                ERROR_OUT()<< "Slave Node Label " << meshA->elems[count + j] << " is not part of slaveNodesMap." << endl;
                }
                slaveDirectElemTable[i]->push_back(slaveNodesMap->at(meshA->elems[count + j]));
            }
            count += numNodesSlaveElem;
        }
        delete slaveNodesMap;
    }

    { // 2. compute slaveDirectElemTable
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
                if(masterNodesMap->find(meshB->elems[count + j])== masterNodesMap->end() ){
                ERROR_OUT()<< "Master Node Label " << meshB->elems[count + j] << " is not part of masterNodesMap." << endl;
                }
                masterDirectElemTable[i]->push_back(masterNodesMap->at(meshB->elems[count + j]));
            }
            count += numNodesMasterElem;
        }
        delete masterNodesMap;
    }

    { // 3. compute slaveNodeToElemTable
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

    { // 4. compute masterNodeToElemTable
        masterNodeToElemTable = new vector<int>*[meshA->numNodes];
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
}

void VertexMorphingMapper::initANNTree() {
#ifdef ANN
    ANNNodes = new double*[meshA->numNodes]; // ANN uses 2D array
    for (int i = 0; i < meshA->numNodes; i++) {
        ANNNodes[i] = new double[3];
        for (int j = 0; j<3; j++)
            ANNNodes[i][j] = meshA->nodes[i * 3 + j];
    }
    slaveNodesTree = new ANNkd_tree(ANNNodes, meshA->numNodes, 3);
#endif
#ifdef FLANN
    FLANNSlaveNodes = new flann::Matrix<double>(const_cast<double*>(meshA->nodes), meshA->numNodes, 3);
    FLANNkd_tree = new flann::Index<flann::L2<double> >(*FLANNSlaveNodes, flann::KDTreeSingleIndexParams(1));
    FLANNkd_tree->buildIndex(); // Build binary tree for searching
#endif
}

void VertexMorphingMapper::findCandidates(double* _controlNode, std::set<int>* _infNodeIdxs, std::set<int>* _infElemIdxs)
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
    numInfNodes = slaveNodesTree->annkFRSearch(_controlNode, filterRadius, 0); // get the number of neighbors in a radius
    infNodeIndices = new int[numInfNodes];
    double *dummy = new double[numInfNodes];
    slaveNodesTree->annkFRSearch(&_controlNode, filterRadius*filterRadius, numInfNodes, infNodeIndices, dummy);// get the real neighbors (ANN uses the square of the radius)
    delete dummy;
#endif

#ifdef FLANN
    flann::Matrix<double> controlNodeCopyFlann(_controlNode, 1, 3);
    vector<vector<int> > indices_tmp;
    vector<vector<double> > dists_tmp;
    FLANNkd_tree->radiusSearch(controlNodeCopyFlann, indices_tmp, dists_tmp, filterRadius*filterRadius, flann::SearchParams(1));
    numInfNodes = indices_tmp[0].size();
    infNodeIndices = new int[indices_tmp[0].size()];
    for (int j=0; j<indices_tmp[0].size(); j++)
        infNodeIndices[j] = indices_tmp[0][j];
#endif

    // 1.2 Fill up the influenced nodes and the influenced elements sets
    for (int i = 0; i < numInfNodes; i++) {
        // Influenced nodes
        _infNodeIdxs->insert(infNodeIndices[i]);
        vector<int> * vecInfElemIndices = slaveNodeToElemTable[infNodeIndices[i]];

        // Fill up unique influenced elements
        for (vector<int>::iterator it = vecInfElemIndices->begin(); it != vecInfElemIndices->end(); it++){
            if (_infElemIdxs->find(*it) == _infElemIdxs->end())
                _infElemIdxs->insert(*it);
        }
    }
    delete[] infNodeIndices;
}

void VertexMorphingMapper::buildCouplingMatrices(){

    // Loop over nodes to find influenced nodes and elements
    const int masterNumNodes = meshB->numNodes;
    for (int iNode = 0; iNode < masterNumNodes; iNode++) {

        // Find candidates
        double* controlNode = &meshB->nodes[iNode*3];
        std::set<int> tmp_infNodes;
        std::set<int> tmp_infElems;
        findCandidates(controlNode, &tmp_infNodes, &tmp_infElems);

        // Loop over the influenced elements
        for (std::set<int>::iterator itElem = tmp_infElems.begin(); itElem != tmp_infElems.end(); itElem++){

            // Loop over the element's nodes and check if they are all inside the filter radius by checking against the influenced nodes list
            std::vector<bool> elemNodeInside;
            for (std::vector<int>::iterator itElemNode = slaveDirectElemTable[*itElem]->begin(); itElemNode != slaveDirectElemTable[*itElem]->end(); itElemNode++)
                elemNodeInside.push_back(tmp_infNodes.find(*itElemNode) != tmp_infNodes.end());

            double contributions[3] = {0.0};
            // if all the nodes are inside do full integration
            if (find(elemNodeInside.begin(), elemNodeInside.end(), false) == elemNodeInside.end()){
                doFullIntegration(controlNode, *itElem, contributions);
            } else { // If some nodes are outside do partial integration
                doPartialIntegration(controlNode, *itElem, &elemNodeInside, contributions);
            }

            // Assemble the contributions from this element (index j) to the control node (index i)
            assemble_C_BA(iNode, *itElem, contributions);
        }
    }

    // Make sure that the integration of the filter function in its effective radius is 1 by C_BA manipulation
    adjustFilter();

    // Delete unnecessary variables
    deleteANNTree();
    deleteTables();
}

void VertexMorphingMapper::doFullIntegration(double* _controlNode, int _elemIdx, double* _contributions){

    int dim = 3;
    int numNodes = 3;

    // Get the triangle
    double tria[9];
    int nodeIdx = -1;
    for (int iNode = 0; iNode<numNodes; iNode++){
        nodeIdx = slaveDirectElemTable[_elemIdx]->at(iNode);
        for (int iXYZ = 0; iXYZ<dim; iXYZ++)
            tria[iNode*dim+iXYZ] = meshA->nodes[nodeIdx*dim+iXYZ];
    }

    // Create Gauss rule on the triangle
    EMPIRE::MathLibrary::GaussQuadratureOnTriangle *gaussQuadratureOnTria =
            new EMPIRE::MathLibrary::GaussQuadratureOnTriangle(tria, numGPsOnTri);

//    // Compute the integration of the filter function
//    This is taken care by C_BA manipulation now (see adjustFilter())
//    for (int iGP = 0; iGP < gaussQuadratureOnTria->numGaussPoints; iGP++)
//        _intFilter += gaussQuadratureOnTria->getWeight(iGP)
//                * filterFunction->computeFunction(_controlNode, &gaussQuadratureOnTria->gaussPointsGlobal[iGP*dim])
//                * gaussQuadratureOnTria->area;

    // Create the integrand on the triangle
    FilterFunctionProduct* integrand = new FilterFunctionProduct(tria, _controlNode);

    integrand->setGaussPoints(gaussQuadratureOnTria->gaussPointsGlobal, gaussQuadratureOnTria->numGaussPoints);
    integrand->computeFunctionProducts(filterFunction);

    for (int iNode = 0; iNode < 3; iNode++) {
        gaussQuadratureOnTria->setIntegrandFunc(integrand);
        integrand->setFunctionID(iNode);
        _contributions[iNode] += gaussQuadratureOnTria->computeIntegral();
    }

    delete gaussQuadratureOnTria;
    delete integrand;

}

void VertexMorphingMapper::doPartialIntegration(double* _controlNode, int _elemIdx, std::vector<bool>* _elemNodeInside, double* _contributions)
{

    int dim = 3;
    int numNodes = 3;

    // Get the triangle
    double tria[9];
    int nodeIdx = -1;
    for (int iNode = 0; iNode<numNodes; iNode++){
        nodeIdx = slaveDirectElemTable[_elemIdx]->at(iNode);
        for (int iXYZ = 0; iXYZ<dim; iXYZ++)
            tria[iNode*dim+iXYZ] = meshA->nodes[nodeIdx*dim+iXYZ];
    }

    // Count the influenced nodes in this element
    int numInfNodesInElem = std::count(_elemNodeInside->begin(), _elemNodeInside->end(), true);
    if (numInfNodesInElem != 1 && numInfNodesInElem != 2 )
        assert(false); // it should not come here

    // find the singled out node
    int nodePos = 0;

    // if there is one node inside
    if (numInfNodesInElem == 1) {
        for (std::vector<bool>::iterator itNode=_elemNodeInside->begin(); !(*itNode); itNode++)
            nodePos++;
    } // else if there is one node outside
    else {
        for (std::vector<bool>::iterator itNode=_elemNodeInside->begin(); *itNode; itNode++)
            nodePos++;
    }

    // assign this node and other nodes coordinates
    double p0[3] = {tria[nodePos*dim],tria[nodePos*dim+1], tria[nodePos*dim+2]};
    double p1[3] = {tria[((nodePos+1)*dim)%9],tria[((nodePos+1)*dim+1)%9], tria[((nodePos+1)*dim+2)%9]};
    double p2[3] = {tria[((nodePos+2)*dim)%9],tria[((nodePos+2)*dim+1)%9], tria[((nodePos+2)*dim+2)%9]};

    // Some precomputations
    double n0n0 = EMPIRE::MathLibrary::computeDenseDotProduct(dim, p0, p0);
    double n1n1 = EMPIRE::MathLibrary::computeDenseDotProduct(dim, p1, p1);
    double n2n2 = EMPIRE::MathLibrary::computeDenseDotProduct(dim, p2, p2);
    double ncnc = EMPIRE::MathLibrary::computeDenseDotProduct(dim, _controlNode, _controlNode);
    double n0n1 = EMPIRE::MathLibrary::computeDenseDotProduct(dim, p0, p1);
    double n0n2 = EMPIRE::MathLibrary::computeDenseDotProduct(dim, p0, p2);
    double n1n2 = EMPIRE::MathLibrary::computeDenseDotProduct(dim, p1, p2);
    double ncn0 = EMPIRE::MathLibrary::computeDenseDotProduct(dim, _controlNode, p0);
    double ncn1 = EMPIRE::MathLibrary::computeDenseDotProduct(dim, _controlNode, p1);
    double ncn2 = EMPIRE::MathLibrary::computeDenseDotProduct(dim, _controlNode, p2);

    /// find clippings
    double xsi1 = 0.0;
    double xsi2 = 0.0;
    // these are initiated as -1 so that they can be distinguished later on
    double xsi31 = -1.0;
    double xsi32 = -1.0;

    /// find clipping between this and next
    // local coordinate
    {
        double a_xsi1 = n0n0 - 2.0*n0n1 + n1n1;
        double b_xsi1 = 2.0*(-n0n0 + n0n1 + ncn0 - ncn1);
        double c_xsi1 = n0n0 - 2.0*ncn0 + ncnc - filterRadius*filterRadius;
        double discriminant = pow(b_xsi1,2)-4.0*a_xsi1*c_xsi1;
        assert (discriminant >= 0); // it cannot be negative
        xsi1 = (-b_xsi1 + sqrt(discriminant)) / (2.0*a_xsi1);
        xsi1 = (xsi1 >= 0.0 && xsi1 <= 1.0) ? xsi1 : (-b_xsi1 - sqrt(discriminant)) / (2.0*a_xsi1);
        assert(xsi1 >= 0.0 && xsi1 <= 1.0); // 0 <= xsi1 <= 1
    }
    // global coordinates
    double p01[3];
    for (int iDim = 0; iDim < dim; iDim++)
        p01[iDim] = (1.0-xsi1) * p0[iDim] + xsi1*p1[iDim];

    /// find clipping between this and prev
    // local coordinate
    {
        double a_xsi2 = n0n0 - 2.0*n0n2 + n2n2;
        double b_xsi2 = 2.0*(-n0n0 + n0n2 + ncn0 - ncn2);
        double c_xsi2 = n0n0 - 2.0*ncn0 + ncnc - filterRadius*filterRadius;

        double discriminant = pow(b_xsi2,2)-4.0*a_xsi2*c_xsi2;
        assert (discriminant >= 0); // it cannot be negative
        xsi2 = (-b_xsi2 + sqrt(discriminant)) / (2.0*a_xsi2);
        xsi2 = (xsi2 >= 0.0 && xsi2 <= 1.0) ? xsi2 : (-b_xsi2 - sqrt(discriminant)) / (2.0*a_xsi2);
        assert(xsi2 >= 0.0 && xsi2 <= 1.0); // 0 <= xsi2 <= 1
    }
    // global coordinates
    double p02[3];
    for (int iDim = 0; iDim < dim; iDim++)
        p02[iDim] = (1.0-xsi2) * p0[iDim] + xsi2*p2[iDim];

    /// try to find clipping between next and prev
    // local coordinate
    double p12_1[3];
    double p12_2[3];
    if (numInfNodesInElem == 1) {

        double a_xsi3 = n1n1 - 2.0*n1n2 + n2n2;
        double b_xsi3 = 2.0*(-n1n1 + n1n2 + ncn1 - ncn2);
        double c_xsi3 = n1n1 - 2.0*ncn1 + ncnc - filterRadius*filterRadius;
        double discriminant = pow(b_xsi3,2)-4.0*a_xsi3*c_xsi3;

        // if there are multiple roots
        if (discriminant > 0) {

            double xsi31_pre = (-b_xsi3 + sqrt(discriminant)) / (2.0*a_xsi3);
            double xsi32_pre = (-b_xsi3 - sqrt(discriminant)) / (2.0*a_xsi3);

            // if both of them are inside set both
            if ((xsi31_pre >= 0.0 && xsi31_pre <= 1.0) && (xsi32_pre >= 0.0 && xsi32_pre<= 1.0)){
                xsi31 = xsi31_pre;
                xsi32 = xsi32_pre;
                for (int iDim = 0; iDim < dim; iDim++)
                    p12_2[iDim] = (1.0-xsi32) * p1[iDim] + xsi32*p2[iDim];

            } // if one of them is outside the limits set only the inside one to xsi31
            else if ((xsi31_pre >= 0.0 && xsi31_pre <= 1.0) && !(xsi32_pre >= 0.0 && xsi32_pre<= 1.0)) {
                xsi31 = xsi31_pre;
            } else if (!(xsi31_pre >= 0.0 && xsi31_pre <= 1.0) && (xsi32_pre >= 0.0 && xsi32_pre<= 1.0)) {
                xsi31 = xsi32_pre;
            }

        } // if there is single root
        else if (discriminant == 0) {
            xsi31 = -b_xsi3 / (2.0*a_xsi3);
        }
        for (int iDim = 0; iDim < dim; iDim++)
            p12_1[iDim] = (1.0-xsi31) * p1[iDim] + xsi31*p2[iDim];
    }

    // make subtriangles
    std::vector<double*> subtrias;
    if (numInfNodesInElem == 1){

        double subtria[9];
        for (int iDim = 0; iDim<dim; iDim++) {
            subtria[iDim] = p0[iDim];
            subtria[dim+iDim] = p01[iDim];
            subtria[2*dim+iDim] = p02[iDim];
        }
        subtrias.push_back(subtria);

        // if the edge across the node has two intersections
        if (xsi31 != -1.0 && xsi32 != -1.0){
            double subtria1[9];
            double subtria2[9];
            for (int iDim = 0; iDim<dim; iDim++) {
                subtria1[iDim] = p01[iDim];
                subtria1[dim+iDim] = p12_1[iDim];
                subtria1[2*dim+iDim] = p02[iDim];

                subtria2[iDim] = p12_1[iDim];
                subtria2[dim+iDim] = p12_2[iDim];
                subtria2[2*dim+iDim] = p02[iDim];
            }
            subtrias.push_back(subtria1);
            subtrias.push_back(subtria2);
        } // if there is only one valid intersection
        else if (xsi31 != -1.0 && xsi32 == -1.0) {
            double subtria1[9];
            for (int iDim = 0; iDim<dim; iDim++) {
                subtria1[iDim] = p01[iDim];
                subtria1[dim+iDim] = p12_1[iDim];
                subtria1[2*dim+iDim] = p02[iDim];
            }
            subtrias.push_back(subtria1);
        }

    } else {

        double subtria1[9];
        double subtria2[9];

        for (int iDim = 0; iDim<dim; iDim++) {
            subtria1[iDim] = p01[iDim];
            subtria1[dim+iDim] = p1[iDim];
            subtria1[2*dim+iDim] = p02[iDim];

            subtria2[iDim] = p1[iDim];
            subtria2[dim+iDim] = p02[iDim];
            subtria2[2*dim+iDim] = p2[iDim];
        }

        subtrias.push_back(subtria1);
        subtrias.push_back(subtria2);

    }

    // Loop over the subtriangles and integrate
    for (std::vector<double*>::iterator itSubtria = subtrias.begin(); itSubtria != subtrias.end(); itSubtria++){

        // Create Gauss rule on subtriangle
        EMPIRE::MathLibrary::GaussQuadratureOnTriangle *gaussQuadratureOnTria =
                new EMPIRE::MathLibrary::GaussQuadratureOnTriangle(*itSubtria, numGPsOnTri);
//        // Compute the integration of the filter function
//        This is taken care by C_BA manipulation now (see adjustFilter())
//        for (int iGP = 0; iGP < gaussQuadratureOnTria->numGaussPoints; iGP++)
//            _intFilter += gaussQuadratureOnTria->getWeight(iGP)
//                    * filterFunction->computeFunction(_controlNode, &gaussQuadratureOnTria->gaussPointsGlobal[iGP*dim])
//                    * gaussQuadratureOnTria->area;

        // Create the integrand on the unclipped triangle
        FilterFunctionProduct* integrand = new FilterFunctionProduct(tria, _controlNode);

        integrand->setGaussPoints(gaussQuadratureOnTria->gaussPointsGlobal, gaussQuadratureOnTria->numGaussPoints);
        integrand->computeFunctionProducts(filterFunction);

        for (int iNode = 0; iNode < 3; iNode++) {
            gaussQuadratureOnTria->setIntegrandFunc(integrand);
            integrand->setFunctionID(iNode);
            _contributions[iNode] += gaussQuadratureOnTria->computeIntegral();
        }

        delete gaussQuadratureOnTria;
        delete integrand;
    }

    subtrias.clear();

}

void VertexMorphingMapper::assemble_C_BA(int _controlNodeIdx, int _elemIdx, double* _contributions){

    int ctr = 0;
    for (std::vector<int>::iterator itElemNode = slaveDirectElemTable[_elemIdx]->begin(); itElemNode != slaveDirectElemTable[_elemIdx]->end(); itElemNode++){
        (*C_BA)(_controlNodeIdx, *itElemNode) += _contributions[ctr];
        ctr += 1;
    }
}

void VertexMorphingMapper::adjustFilter(){

    int ctr = 0;
    double fac = 0.0;
    // Loop over the rows
    for (int iRow = 0; iRow < C_BA->getNumberOfRows(); iRow++){
        // Count the nonzero columns
        ctr = 0;
        for (int iColumn = 0; iColumn < C_BA->getNumberOfColumns(); iColumn++){
            if ((*C_BA)(iRow,iColumn) != 0.0);
                ctr++;
        }

        // Compute the factor for making the C_BA consistent and multiply the row
        fac = static_cast<double>(ctr)/(C_BA->getRowSum(iRow)*C_BA->getNumberOfColumns());
        C_BA->multiplyRowWith(iRow, fac);
    }

}

void VertexMorphingMapper::consistentMapping(const double *slaveField, double *masterField) {
    double *slaveFieldCopy = new double[meshA->numNodes];

//    INFO_OUT() << "MortarMapper::consistentMapping ->  Norm of the Slave field :: "  << EMPIRE::MathLibrary::computeVectorLength(slaveField) << endl;
//    INFO_OUT() << "MortarMapper::consistentMapping ->  Norm of the Master field :: " << EMPIRE::MathLibrary::computeVectorLength(masterField) << endl;

    for (int i = 0; i < meshA->numNodes; i++)
        slaveFieldCopy[i] = slaveField[i];

    // matrix vector product (W_tmp = C_BA * W_A)
    (*C_BA).mulitplyVec(false,slaveFieldCopy,masterField,meshB->numNodes);

    delete[] slaveFieldCopy;

}

// FilterFunction classes
double VertexMorphingMapper::HatFilterFunction::computeFunction(double* _supportCenter, double* _globalCoor){

    // compute distance to the support center
    double dist = EMPIRE::MathLibrary::computeDenseEuclideanNorm(3, _supportCenter, _globalCoor);

    // return function value
    return 1.0 - dist/filterRadius;
}

double VertexMorphingMapper::GaussianFilterFunction::computeFunction(double* _supportCenter, double* _globalCoor){

    // compute distance to the support center
    double dist = EMPIRE::MathLibrary::computeDenseEuclideanNorm(3, _supportCenter, _globalCoor);

    // return function value
    return exp(-pow(dist,2) / (2 * pow(filterRadius,2) / 9.0));
}

// Integrand classes
VertexMorphingMapper::FilterFunctionProduct::FilterFunctionProduct(double* _elem, double* _controlNode)
{

    int dim = 3;
    int numNodes = 3;

    gaussPoints = NULL;
    functionProducts = NULL;

    elem = new double[9];
    for (int iNode = 0; iNode < numNodes; iNode++){
        for (int iXYZ = 0; iXYZ < dim; iXYZ++)
           elem[iNode*dim + iXYZ] = _elem[iNode*dim + iXYZ];
    }

    controlNode = new double[3];
    for (int iXYZ = 0; iXYZ < 3; iXYZ++)
        controlNode[iXYZ] = _controlNode[iXYZ];

}

VertexMorphingMapper::FilterFunctionProduct::~FilterFunctionProduct(){

    delete controlNode;
    delete gaussPoints;
    delete elem;

    for (int i=0; i < 3; i++)
        delete functionProducts[i];
    delete functionProducts;

}

void VertexMorphingMapper::FilterFunctionProduct::setGaussPoints(const double *_gaussPoints, int _numGaussPoints) {
    assert(gaussPoints == NULL);
    numGaussPoints = _numGaussPoints;
    gaussPoints = new double[numGaussPoints * 3];
    for (int i = 0; i < numGaussPoints * 3; i++)
        gaussPoints[i] = _gaussPoints[i];
}

void VertexMorphingMapper::FilterFunctionProduct::computeFunctionProducts(AbstractFilterFunction *_filterFunction){
    assert(functionProducts == NULL);
    assert(gaussPoints != NULL);

    int dim = 3;
    int numNodes = 3;

    // initiate function products container
    functionProducts = new double*[numNodes];
    for (int i = 0; i < numNodes; i++)
        functionProducts[i] = new double[numGaussPoints];

    // loop over gauss points
    for (int iGP = 0; iGP < numGaussPoints; iGP++) {
        double *gaussPoint = &gaussPoints[iGP * dim];

        // *. compute shape functions on element (on the Gauss point)
        double shapeFuncValues[3];
        double normal[3];
        EMPIRE::MathLibrary::computeNormalOfTriangle(elem, true, normal);
        int planeToProject = EMPIRE::MathLibrary::computePlaneToProject(normal);
        bool inside = EMPIRE::MathLibrary::computeLocalCoorInTriangle(elem, planeToProject,
                                                                      gaussPoint, shapeFuncValues);

        // debug
        if (!inside){
            cout<<"Error in computing local coordinates in tria master element"<<endl;
            cout<<"GP coordinates: "<<gaussPoint[0]<<"\t"<<gaussPoint[1]<<"\t"<<gaussPoint[2]<<endl;
            cout<<"Element nodes:"<<endl;
            for(int ctr=0;ctr<numNodes;ctr++){
                cout<<ctr+1<<": "<<elem[ctr*dim]<<"\t"<<elem[ctr*dim+1]<<"\t"<<elem[ctr*dim+2]<<endl;
            }
        }
        // debug end
        assert(inside);

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


} /* namespace EMPIRE */

