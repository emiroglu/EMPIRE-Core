/*  Copyright &copy; 2013, TU Muenchen, Chair of Structural Analysis,
 *  Fabien Pean, Andreas Apostolatos, Chenshen Wu,
 *  Ragnar Björnsson, Munich
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

#include "IGAMortarCouplingMatrices.h"
#include "MathLibrary.h"
#include "IGAMortarMapper.h"
#include "IGAPatchSurface.h"
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

using namespace std;

namespace EMPIRE {

IGAMortarCouplingMatrices::IGAMortarCouplingMatrices(int _size_N , int _size_R, bool _isExpanded)
{
    // Initialize sizes
    size_N = _size_N;
    size_R = _size_R;

    // Initialize coupling matrices
    Cnn = new MathLibrary::SparseMatrix<double>(size_N, false);
    Cnr = new MathLibrary::SparseMatrix<double>(size_N, size_R);

    // Flag on whether the expanded version of the coupling matrices is assumed
    isExpanded = _isExpanded;
}


IGAMortarCouplingMatrices::~IGAMortarCouplingMatrices() {
    delete Cnn;
    delete Cnr;
}

void IGAMortarCouplingMatrices::factorizeCnn() {
    Cnn->factorize();
}

void IGAMortarCouplingMatrices::enforceCnn() {
    /*
     * Checks if a row is empty and if yes adds 1.0 in the diagonal
     */
    indexEmptyRowCnn.reserve(size_N);
    for(int i = 0; i < size_N; i++) {
        if(Cnn->isRowEmpty(i)) {
            (*Cnn)(i,i) = 1;
            indexEmptyRowCnn.push_back(i);
        }
    }
}

}

//void IGAMortarCouplingMatrices::applyDirichletBCs(std::vector<int> clampedIds) {

//    if(clampedIds.size() > 0) {
//        if(isClampedDofs || isIGAWeakDirichletConditions || isIGAPatchContinuityConditions ) {
//            for(int i = 0 ; i < 3*size_N ; i++) {
//                if ( std::find(clampedIds.begin(), clampedIds.end(), i)!=clampedIds.end() ) {
//                    for(int j = 0 ; j < i ; j++) {
//                        (*C_NN_BCs)(j,i) = 0;
//                        (*C_NN_BCs)(i,j) = 0;
//                    }
//                    (*C_NN_BCs)(i,i) = 1;
//                }
//                else {
//                    for(int j = 0 ; j < 3*size_N ; j++) {
//                        if(!(std::find(clampedIds.begin(), clampedIds.end(), j)!=clampedIds.end()))
//                            (*C_NN_BCs)(i,j) = (*C_NN_expanded)(i,j);
//                    }
//                    for(int j = 0 ; j < 3*size_R ; j++) {
//                        (*C_NR_BCs)(i,j) = (*C_NR_expanded)(i,j);
//                    }
//                }
//            }
//        }
//        else {

//            std::vector<int> clampedNodes;
//            for(int i = 0 ; i < clampedIds.size() ; i++) {
//                clampedNodes.push_back((int)clampedIds[i]/3);
//                std::cout<<"clampedNodes = "<<clampedNodes[i]<<std::endl;
//            }

//            for(int i = 0 ; i < size_N ; i++) {
//                if ( std::find(clampedNodes.begin(), clampedNodes.end(), i)!=clampedNodes.end() ) {
//                    for(int j = 0 ; j < i ; j++) {
//                        (*C_NN_BCs)(j,i) = 0;
//                        (*C_NN_BCs)(i,j) = 0;
//                    }
//                    (*C_NN_BCs)(i,i) = 1;
//                }
//                else {
//                    for(int j = 0 ; j < size_N ; j++) {
//                        if(!(std::find(clampedNodes.begin(), clampedNodes.end(), j)!=clampedNodes.end()))
//                            (*C_NN_BCs)(i,j) = (*C_NN)(i,j);
//                    }
//                    for(int j = 0 ; j < size_R ; j++) {
//                        (*C_NR_BCs)(i,j) = (*C_NR)(i,j);
//                    }
//                }
//            }
//        }
//    }
//    else
//        WARNING_OUT()<<"Nothing to clamp"<<std::endl;

//    INFO_OUT()<<"Applying Dirichlet BCs finished"<<std::endl;
//}
