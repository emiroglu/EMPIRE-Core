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
#include "WeakIGADirichletCondition.h"
#include "Message.h"
#include <iostream>
#include <assert.h>

using namespace std;

namespace EMPIRE {

WeakIGADirichletCondition::WeakIGADirichletCondition(int _ID,
                                 int _patchIndex, int _patchBLIndex, int _patchBLTrCurveIndex,
                                 int _isGPProvided) :
        AbstractCondition(_ID),
        patchIndex(_patchIndex), patchBLIndex(_patchBLIndex), patchBLTrCurveIndex(_patchBLTrCurveIndex),
        isGPProvided(_isGPProvided)
        {

    type = EMPIRE_WeakIGADirichletCondition;

}

void WeakIGADirichletCondition::addWeakDirichletConditionGPData(int _trCurveNumGP,
                            double* _trCurveGPs, double* _trCurveGPWeights,
                            double* _trCurveGPTangents,
                            double* _trCurveGPJacobianProducts){

    trCurveNumGP = _trCurveNumGP;

    trCurveGPs = new double[2*_trCurveNumGP];
    trCurveGPWeights = new double[_trCurveNumGP];
    trCurveGPTangents = new double[3*_trCurveNumGP];
    trCurveGPJacobianProducts = new double[_trCurveNumGP];

    for (int i = 0; i < trCurveNumGP; i++){
        for (int j = 0; j < 2; j++){
            trCurveGPs[2*i+j]=_trCurveGPs[2*i+j];
        }
        trCurveGPWeights[i] = _trCurveGPWeights[i];
        for (int j = 0; j < 3; j++){
            trCurveGPTangents[3*i+j]=_trCurveGPTangents[3*i+j];
        }
        trCurveGPJacobianProducts[i] = _trCurveGPJacobianProducts[i];
    }
}

WeakIGADirichletCondition::~WeakIGADirichletCondition() {
    delete trCurveGPs;
    delete trCurveGPWeights;
    delete trCurveGPTangents;
    delete trCurveGPJacobianProducts;
}


} /* namespace EMPIRE */
