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
#include "WeakIGAPatchContinuityCondition.h"
#include "Message.h"
#include <iostream>
#include <assert.h>

using namespace std;

namespace EMPIRE {

WeakIGAPatchContinuityCondition::WeakIGAPatchContinuityCondition(int _ID, 
								 int _masterPatchIndex, int _masterPatchBLIndex, int _masterPatchBLTrCurveIndex,
								 int _slavePatchIndex, int _slavePatchBLIndex, int _slavePatchBLTrCurveIndex,
                                 int _isGPProvided) :
        AbstractCondition(_ID),
        masterPatchIndex(_masterPatchIndex), masterPatchBLIndex(_masterPatchBLIndex), masterPatchBLTrCurveIndex(_masterPatchBLTrCurveIndex),
        slavePatchIndex(_slavePatchIndex), slavePatchBLIndex(_slavePatchBLIndex), slavePatchBLTrCurveIndex(_slavePatchBLTrCurveIndex),
        isGPProvided(_isGPProvided)
        {

    type = EMPIRE_WeakIGAPatchContinuityCondition;

}

void WeakIGAPatchContinuityCondition::addWeakContinuityConditionGPData(int _trCurveNumGP,
							double* _trCurveMasterGPs, double* _trCurveSlaveGPs, double* _trCurveGPWeights,
							double* _trCurveMasterGPTangents, double* _trCurveSlaveGPTangents,
							double* _trCurveGPJacobianProducts){
  
    trCurveNumGP = _trCurveNumGP;

    trCurveMasterGPs = new double[2*_trCurveNumGP];
    trCurveSlaveGPs = new double[2*_trCurveNumGP];
    trCurveGPWeights = new double[_trCurveNumGP];
    trCurveMasterGPTangents = new double[3*_trCurveNumGP];
    trCurveSlaveGPTangents = new double[3*_trCurveNumGP];
    trCurveGPJacobianProducts = new double[_trCurveNumGP];

    for (int i = 0; i < trCurveNumGP; i++){
        for (int j = 0; j < 2; j++){
            trCurveMasterGPs[2*i+j]=_trCurveMasterGPs[2*i+j];
            trCurveSlaveGPs[2*i+j]=_trCurveSlaveGPs[2*i+j];
        }
        trCurveGPWeights[i] = _trCurveGPWeights[i];
        for (int j = 0; j < 3; j++){
            trCurveMasterGPTangents[3*i+j]=_trCurveMasterGPTangents[3*i+j];
            trCurveSlaveGPTangents[3*i+j]=_trCurveSlaveGPTangents[3*i+j];
        }
        trCurveGPJacobianProducts[i] = _trCurveGPJacobianProducts[i];
    }
}

WeakIGAPatchContinuityCondition::~WeakIGAPatchContinuityCondition() {
    delete trCurveMasterGPs;
    delete trCurveSlaveGPs;
    delete trCurveGPWeights;
    delete trCurveMasterGPTangents;
    delete trCurveSlaveGPTangents;
    delete trCurveGPJacobianProducts;
}


} /* namespace EMPIRE */
