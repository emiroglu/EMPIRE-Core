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
#include "DataFieldIntegration.h"
#include "MathLibrary.h"
#include "Message.h"
#include <map>
#include <vector>
#include <assert.h>
#include <iostream>
#include <stdlib.h>


using namespace std;

namespace EMPIRE {

int DataFieldIntegration::mklSetNumThreads = 1; // set default value for MKL threads
const int DataFieldIntegration::numGPsMassMatrixTri = 6;
const int DataFieldIntegration::numGPsMassMatrixQuad = 4;

DataFieldIntegration::DataFieldIntegration(int _numNodes, int _numElems,
        const int *_numNodesPerElem, const double *_nodes, const int *_nodeIDs, const int *_elems) :
        numNodes(_numNodes), numElems(_numElems), numNodesPerElem(_numNodesPerElem), nodes(_nodes), nodeIDs(
                _nodeIDs), elems(_elems) {
    // Edit Aditya
    massMatrix = new EMPIRE::MathLibrary::SparseMatrix<double>(numNodes,true);

    // compute directElemTable which link the element to the all its nodes' positions (instead of IDs)
    vector<int> **directElemTable = new vector<int>*[numElems];
    for (int i = 0; i < numElems; i++)
        directElemTable[i] = new vector<int>;
    map<int, int> *nodesMap = new map<int, int>();
    for (int i = 0; i < numNodes; i++)
        nodesMap->insert(nodesMap->end(), pair<int, int>(nodeIDs[i], i));
    int count = 0;
    for (int i = 0; i < numElems; i++) {
        const int numNodesThisElem = numNodesPerElem[i];
        for (int j = 0; j < numNodesThisElem; j++) {
            directElemTable[i]->push_back(nodesMap->at(elems[count + j]));
        }
        count += numNodesThisElem;
    }
    delete nodesMap;

    // compute the sparsity map
    // sparsity map has the information of a, ia, ja in a CSR formated matrix
    map<int, double> **sparsityMap = new map<int, double>*[numNodes];
    for (int i = 0; i < numNodes; i++)
        sparsityMap[i] = new map<int, double>;

    for (int i = 0; i < numElems; i++) {
        const int numNodesThisElem = numNodesPerElem[i];
        double elem[numNodesThisElem * 3]; // this element
        int pos[numNodesThisElem]; // the position of the node in the nodeCoors
        for (int j = 0; j < numNodesThisElem; j++) {
            pos[j] = directElemTable[i]->at(j);
            for (int k = 0; k < 3; k++)
                elem[j * 3 + k] = nodes[pos[j] * 3 + k];
        }

        // replace the master element by the projection of it on its "element plane"
        // we do it here because we have done the same in MortarMapper
        if (numNodesThisElem == 4) {
            double masterElemNormal[3];
            EMPIRE::MathLibrary::computeNormalOfQuad(elem, true, masterElemNormal);
            double masterQuadCenter[3];
            EMPIRE::MathLibrary::computePolygonCenter(elem, 4, masterQuadCenter);
            double masterQuadPrj[12];
            EMPIRE::MathLibrary::projectToPlane(masterQuadCenter, masterElemNormal, elem, 4, masterQuadPrj);
            for (int i = 0; i < 12; i++)
                elem[i] = masterQuadPrj[i];
        }

        // make use of the symmetry
        double massMatrixElem[numNodesThisElem * numNodesThisElem];
        if (numNodesThisElem == 4)
        	EMPIRE::MathLibrary::computeMassMatrixOfQuad(elem, numGPsMassMatrixQuad, false, massMatrixElem);
        else if (numNodesThisElem == 3)
        	EMPIRE::MathLibrary::computeMassMatrixOfTrianlge(elem, numGPsMassMatrixTri, false,
                    massMatrixElem);
        else
            assert(false);
        for (int j = 0; j < numNodesThisElem; j++) {
            for (int k = j; k < numNodesThisElem; k++) {
                int smaller, larger;
                if (pos[j] > pos[k]) {
                    larger = pos[j];
                    smaller = pos[k];
                } else {
                    larger = pos[k];
                    smaller = pos[j];
                }
                double massMatrixJK = massMatrixElem[j * numNodesThisElem + k];
                bool inserted = sparsityMap[smaller]->insert(
                        pair<int, double>(larger, massMatrixJK)).second; // *.second is a bool indicating inserted or not
                if (!inserted)
                    (sparsityMap[smaller]->at(larger)) += massMatrixJK; // at() returns a reference, so using "+=" is correct
            }
        }
    }

    // 2. according to sparsity map, compute a, ia, ja of CSR formated massMatrix
    int nnz = 1; // number of non-zero entries
    for (int i = 0; i < numNodes; i++) {
        nnz += sparsityMap[i]->size();
    }
    nnz--;
    int count2 = 0;
    for (int i = 0; i < numNodes; i++) {
        for (map<int, double>::iterator it = sparsityMap[i]->begin(); it != sparsityMap[i]->end();
                it++) {
        	// Edit Aditya
        	(*massMatrix)(i,count2) = it->second;
            count2++;
        }
    }assert(nnz == count2);

    // delete spasity map
    for (int i = 0; i < numNodes; i++)
        delete sparsityMap[i];
    delete[] sparsityMap;

    for (int i = 0; i < numElems; i++)
        delete directElemTable[i];
    delete[] directElemTable;

    //MortarMath::printCSRMatrixSymmetric(massMatrix_A, massMatrix_IA, massMatrix_JA, numNodes);
}

DataFieldIntegration::~DataFieldIntegration() {
    delete[] massMatrix;
}

void DataFieldIntegration::integrate(const double *tractions, double *forces) {
    int n = numNodes;
    char up = 'u';
    // This routine supports only one-based indexing of the input arrays.
    // Edit Aditya
    massMatrix->determineCSR();
    massMatrix->mulitplyVec(false,const_cast<double *>(tractions),forces,n);

}

void DataFieldIntegration::deIntegrate(const double *forces, double *tractions) {
    int n = numNodes;
    char up = 'u';
    // use pardiso solver
    // pardiso forward and backward substitution
    double *ddum = new double[n]; // dummy but the memory is asked for
    // Edit Aditya
    massMatrix->solve(tractions,const_cast<double*>(forces));

    delete[] ddum;
}

} /* namespace EMPIRE */
