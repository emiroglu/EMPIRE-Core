/*  Copyright &copy; 2013, TU Muenchen, Chair of Structural Analysis,
 *  Stefan Sicklinger, Tianyang Wang, Andreas Apostolatos, Munich
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
// Inclusion of standard libraries
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

// Inclusion of user defined libraries
#include "IGAControlPoint.h"
#include "IGAPatchSurface.h"
#include "IGAMesh.h"
#include "WeakIGADirichletCurveCondition.h"
#include "WeakIGADirichletSurfaceCondition.h"
#include "WeakIGAPatchContinuityCondition.h"
#include "DataField.h"
#include "Message.h"

using namespace std;

namespace EMPIRE {

IGAMesh::IGAMesh(std::string _name, int _numNodes) :
    AbstractMesh(_name), numNodes(_numNodes) {
    type = EMPIRE_Mesh_IGAMesh;

}

IGAMesh::~IGAMesh() {

    for (int i = 0; i < surfacePatches.size(); i++)
        delete surfacePatches[i];

    for (int i = 0; i < weakIGADirichletCurveConditions.size(); i++)
        delete weakIGADirichletCurveConditions[i];

    for (int i = 0; i < weakIGADirichletSurfaceConditions.size(); i++)
        delete weakIGADirichletSurfaceConditions[i];

    for (int i = 0; i < weakIGAPatchContinuityConditions.size(); i++)
        delete weakIGAPatchContinuityConditions[i];
}

IGAPatchSurface* IGAMesh::addPatch(int _pDegree, int _uNoKnots, double* _uKnotVector, int _qDegree,
                                   int _vNoKnots, double* _vKnotVector, int _uNoControlPoints, int _vNoControlPoints,
                                   double* _controlPointNet, int* _dofIndexNet) {

    std::string patchName = name + " Patch";
    int IDBasis = 0;

    int numCPs = _uNoControlPoints * _vNoControlPoints;
    IGAControlPoint **cpNet;
    cpNet = new IGAControlPoint*[numCPs];

    for (int i = 0; i < numCPs; i++) {
        if (_dofIndexNet[i] < numNodes && _dofIndexNet[i] >= 0)
            cpNet[i] = new IGAControlPoint(_dofIndexNet[i], &_controlPointNet[i * 4]);
        else {
            ERROR_OUT() << "DOF " << _dofIndexNet[i] << " has not been defined" << endl;
            exit(-1);
        }
    }

    surfacePatches.push_back(
                new IGAPatchSurface(IDBasis, _pDegree, _uNoKnots, _uKnotVector, _qDegree, _vNoKnots,
                                    _vKnotVector, _uNoControlPoints, _vNoControlPoints, cpNet));
    return surfacePatches.back();
}

void IGAMesh::computeBoundingBox() {
    if (boundingBox.isComputed())
        return;

    boundingBox[0] = surfacePatches[0]->getControlPointNet()[0]->getX();
    boundingBox[1] = surfacePatches[0]->getControlPointNet()[0]->getX();
    boundingBox[2] = surfacePatches[0]->getControlPointNet()[0]->getY();
    boundingBox[3] = surfacePatches[0]->getControlPointNet()[0]->getY();
    boundingBox[4] = surfacePatches[0]->getControlPointNet()[0]->getZ();
    boundingBox[5] = surfacePatches[0]->getControlPointNet()[0]->getZ();

    for (int patchCount = 0; patchCount < surfacePatches.size(); patchCount++) {
        IGAPatchSurface* patch = surfacePatches[patchCount];
        patch->computeBoundingBox();
    }
    for (int patchCount = 0; patchCount < surfacePatches.size(); patchCount++) {
        IGAPatchSurface* patch = surfacePatches[patchCount];
        double xmin = patch->getBoundingBox(0);
        double xmax = patch->getBoundingBox(1);
        double ymin = patch->getBoundingBox(2);
        double ymax = patch->getBoundingBox(3);
        double zmin = patch->getBoundingBox(4);
        double zmax = patch->getBoundingBox(5);
        if (xmin < boundingBox[0])
            boundingBox[0] = xmin;
        else if ( xmax > boundingBox[1])
            boundingBox[1] = xmax;
        if (ymin < boundingBox[2])
            boundingBox[2] = ymin;
        else if (ymax > boundingBox[3])
            boundingBox[3] = ymax;
        if (zmin < boundingBox[4])
            boundingBox[4] = zmin;
        else if (zmax > boundingBox[5])
            boundingBox[5] = zmax;
    }
    boundingBox.isComputed(true);
}

void IGAMesh::addDataField(string _dataFieldName, EMPIRE_DataField_location _location,
                           EMPIRE_DataField_dimension _dimension, EMPIRE_DataField_typeOfQuantity _typeOfQuantity) {
    int numLocations = -1;
    if (_location == EMPIRE_DataField_atNode)
        numLocations = numNodes;
    else
        assert(false);

    assert(nameToDataFieldMap.find(_dataFieldName) == nameToDataFieldMap.end());
    DataField *dataField = new DataField(_dataFieldName, _location, numLocations, _dimension,
                                         _typeOfQuantity);
    nameToDataFieldMap.insert(pair<string, DataField*>(_dataFieldName, dataField));

}

WeakIGADirichletCurveCondition* IGAMesh::addWeakDirichletCurveCondition(int _conditionID,
                                                                        int _patchIndex,
                                                                        int _pDegree, int _uNoKnots, double* _uKnotVector,
                                                                        int _uNoControlPoints, double* _controlPointNet){

    weakIGADirichletCurveConditions.push_back(new WeakIGADirichletCurveCondition(_conditionID, _patchIndex,
                                                                                 new IGAPatchCurve(0, _pDegree, _uNoKnots, _uKnotVector, _uNoControlPoints, _controlPointNet)));

    return weakIGADirichletCurveConditions.back();

}

WeakIGADirichletCurveCondition* IGAMesh::addWeakDirichletCurveCondition(int _conditionID,
                                                                        int _patchIndex, IGAPatchCurve* _dirichletCurve){

    weakIGADirichletCurveConditions.push_back(new WeakIGADirichletCurveCondition(_conditionID, _patchIndex, _dirichletCurve));

    return weakIGADirichletCurveConditions.back();

}

WeakIGADirichletCurveCondition* IGAMesh::addWeakDirichletCurveCondition(int _conditionID,
                                                                        int _patchIndex, int _patchBLIndex, int _patchBLTrCurveIndex){

    weakIGADirichletCurveConditions.push_back(new WeakIGADirichletCurveCondition(_conditionID,
                                                                                 _patchIndex, _patchBLIndex, _patchBLTrCurveIndex));
    return weakIGADirichletCurveConditions.back();

}

void IGAMesh::createWeakDirichletCurveConditionGPData(int _conditionIndex){

    /*
     * This function calls the corresponding function of the weak Dirichlet condition
     */

    int patchIndex = weakIGADirichletCurveConditions[_conditionIndex]->getPatchIndex();

    // Get the patch
    IGAPatchSurface* thePatch = getSurfacePatch(patchIndex);

    // Create GP Data for the patch
    weakIGADirichletCurveConditions[_conditionIndex]->createGPData(thePatch);

}

WeakIGADirichletSurfaceCondition* IGAMesh::addWeakDirichletSurfaceCondition(int _conditionID, int _patchIndex, IGAPatchSurfaceTrimmingLoop* theLoop) {

    weakIGADirichletSurfaceConditions.push_back(new WeakIGADirichletSurfaceCondition(_conditionID, _patchIndex, theLoop));

    return weakIGADirichletSurfaceConditions.back();
}

WeakIGADirichletSurfaceCondition* IGAMesh::addWeakDirichletSurfaceCondition(int _conditionID, int _patchIndex, int _patchBLIndex) {

    weakIGADirichletSurfaceConditions.push_back(new WeakIGADirichletSurfaceCondition(_conditionID, _patchIndex, _patchBLIndex));

    return weakIGADirichletSurfaceConditions.back();
}

void IGAMesh::createWeakDirichletSurfaceConditionGPData(int _connectionIndex) {

    /*
     * This function calls the corresponding function of the weak Dirichlet condition
     */

    int patchIndex = weakIGADirichletSurfaceConditions[_connectionIndex]->getPatchIndex();

    // Get the patch
    IGAPatchSurface* thePatch = getSurfacePatch(patchIndex);

    // Create GP Data for the patch
    weakIGADirichletSurfaceConditions[_connectionIndex]->createGPData(thePatch);
}

WeakIGAPatchContinuityCondition* IGAMesh::addWeakContinuityCondition(int _connectionID,
                                                                     int _masterPatchIndex, int _masterPatchBLIndex, int _masterPatchBLTrCurveIndex,
                                                                     int _slavePatchIndex, int _slavePatchBLIndex, int _slavePatchBLTrCurveIndex){

    weakIGAPatchContinuityConditions.push_back(new WeakIGAPatchContinuityCondition(_connectionID,
                                                                                   _masterPatchIndex, _masterPatchBLIndex, _masterPatchBLTrCurveIndex,
                                                                                   _slavePatchIndex, _slavePatchBLIndex, _slavePatchBLTrCurveIndex));
    return weakIGAPatchContinuityConditions.back();

}

WeakIGAPatchContinuityCondition* IGAMesh::addWeakContinuityCondition(int _connectionID,
                                                                     int _masterPatchIndex, IGAPatchCurve* _masterCurve,
                                                                     int _slavePatchIndex,  IGAPatchCurve* _slaveCurve) {

    weakIGAPatchContinuityConditions.push_back(new WeakIGAPatchContinuityCondition(_connectionID,
                                                                                   _masterPatchIndex, _masterCurve,
                                                                                   _slavePatchIndex,  _slaveCurve));
    return weakIGAPatchContinuityConditions.back();

}

WeakIGAPatchContinuityCondition* IGAMesh::addWeakContinuityCondition(int _connectionID,
                                                                     int _masterPatchIndex, int _pMaster, int _uNoKnotsMaster, double* _uKnotVectorMaster, int _uNoControlPointsMaster, double* _controlPointNetMaster,
                                                                     int _slavePatchIndex,  int _pSlave, int _uNoKnotsSlave, double* _uKnotVectorSlave, int _uNoControlPointsSlave, double* _controlPointNetSlave) {

    weakIGAPatchContinuityConditions.push_back(new WeakIGAPatchContinuityCondition(_connectionID,
                                                                                   _masterPatchIndex, _pMaster, _uNoKnotsMaster, _uKnotVectorMaster, _uNoControlPointsMaster, _controlPointNetMaster,
                                                                                   _slavePatchIndex,  _pSlave, _uNoKnotsSlave, _uKnotVectorSlave, _uNoControlPointsSlave, _controlPointNetSlave));
    return weakIGAPatchContinuityConditions.back();

}

void IGAMesh::createWeakContinuityConditionGPData(int _connectionIndex){

    /*
     * This function calls the corresponding function of the weak continuity condition
     */

    int masterPatchIndex = weakIGAPatchContinuityConditions[_connectionIndex]->getMasterPatchIndex();
    int slavePatchIndex = weakIGAPatchContinuityConditions[_connectionIndex]->getSlavePatchIndex();

    // Get the master and slave patches
    IGAPatchSurface* masterPatch = getSurfacePatch(masterPatchIndex);
    IGAPatchSurface* slavePatch = getSurfacePatch(slavePatchIndex);

    // Create GP Data for the patch pair
    weakIGAPatchContinuityConditions[_connectionIndex]->createGPData(masterPatch, slavePatch);

}

Message &operator<<(Message & _message, const IGAMesh & _mesh) {
    _message <<endl;
    _message << "\t" << "---------------------------------Start Mesh" << endl;
    _message << "\t" << "IGA Mesh name: " << _mesh.name << endl;
    _message << "\t\tNumber of Patches: " << _mesh.getSurfacePatches().size() << endl;
    _message << "\t\tNumber of Nodes: " << _mesh.getNumNodes() << endl;
    //_message << "\t" << "---------------------------------" << endl;
    for (int k = 0; k < _mesh.getSurfacePatches().size(); k++) {
        _message << "\t" << "---------------------------------Start Patch" << endl;
        _message << "\tPatch[" << k << "]:" << endl;
        _message << *(_mesh.getSurfacePatches()[k]);
    }
    _message << "\t" << "---------------------------------End Mesh" << endl;
    return _message;
}

}/* namespace EMPIRE */

