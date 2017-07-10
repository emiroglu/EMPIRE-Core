/*  Copyright &copy; 2013, TU Muenchen, Chair of Structural Analysis,
 *  Fabien Pean, Andreas Apostolatos, Chenshen Wu,
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
/***********************************************************************************************//**
 * \file IGAMesh.h
 * This file holds the class IGAMesh.h
 * \date 6/8/2013
 **************************************************************************************************/

#ifndef IGAMesh_H_
#define IGAMesh_H_

#include <string>
#include <cfloat>
// Inclusion of user defined libraries
#include "AbstractMesh.h"
#include "IGAPatchCouplingCaratData.h"
#include "WeakIGAPatchContinuityCondition.h"
#include "WeakIGADirichletCondition.h"

namespace EMPIRE {
class DataField;
class Message;
class IGAPatchSurface;
class IGAControlPoint;

class IGAPatchCouplingCaratData;

/********//**
 * \brief class IGAMesh is a specialization of the class AbstractMesh used for IGA Mesh containing number of IGA surface patches
 ***********/

class IGAMesh: public AbstractMesh {

protected:

    /// Array of IGA Surface Patches
    std::vector<IGAPatchSurface*> surfacePatches;

    /// Vector of all clampedDofs
    std::vector<int> clampedDofs;

    /// the lowest number of clamped directions
    int clampedDirections;

    /// Object which has all the patch coupling data
    IGAPatchCouplingCaratData* couplingData;

    /// Vector of all weak Dirichlet conditions
    std::vector<WeakIGADirichletCondition*> weakIGADirichletConditions;

    /// Vector of all weak patch continuity conditions
    std::vector<WeakIGAPatchContinuityCondition*> weakIGAPatchContinuityConditions;

    /// The number of the Control Points in the IGAMesh
    int numNodes;

    /// The constructor, the destructor and the copy constructor
public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _name The name of the IGA mesh
     * \param[in] _numControlPoints The number of the Control Points
     * \author Chenshen Wu
     ***********/
    IGAMesh(std::string _name, int _numNodes);

    /***********************************************************************************************
     * \brief Destructor
     * \author Chenshen Wu
     ***********/
    ~IGAMesh();

    /***********************************************************************************************
     * brief Add a new surface patch to the IGA mesh
     * \param[in] _pDegree The polynomial degree of the IGA 2D patch in the u-direction
     * \param[in] _uNoKnots The number of knots for the knot vector in the u-direction
     * \param[in] _uKnotVector The underlying knot vector of the IGA 2D patch in the u-direction
     * \param[in] _qDegree The polynomial degree of the IGA 2D patch in the v-direction
     * \param[in] _vNoKnots The number of knots for the knot vector in the v-direction
     * \param[in] _vKnotVector The underlying knot vector of the IGA 2D patch in the v-direction
     * \param[in] _uNoControlPoints The number of the Control Points for the 2D NURBS patch in the u-direction
     * \param[in] _vNoControlPoints The number of the Control Points for the 2D NURBS patch in the v-direction
     * \param[in] _controlPointNet The set of the Control Points related to the 2D NURBS patch
     * \param[in] _dofIndexNet The index of the dof of the each Control Points related to
     * \return The pointer to the patch just created
     * \author Chenshen Wu
     ***********/
    IGAPatchSurface* addPatch(int _pDegree, int _uNoKnots, double* _uKnotVector, int _qDegree, int _vNoKnots,
                  double* _vKnotVector, int _uNoControlPoints, int _vNoControlPoints,
                  double* controlPointNet, int* _dofIndexNet);

    /// Specializing abstract functions from AbstractMesh class
public:
    /***********************************************************************************************
     * \brief Add a new data field to this mesh
     * \param[in] _dataFieldName name of the data field
     * \param[in] _location at node or at element centroid
     * \param[in] _dimension vector or scalar
     * \param[in] _typeOfQuantity field or field integral
     * \author Chenshen Wu
     ***********/
    void addDataField(std::string _dataFieldName, EMPIRE_DataField_location _location,
            EMPIRE_DataField_dimension _dimension, EMPIRE_DataField_typeOfQuantity _typeOfQuantity);

    /***********************************************************************************************
     * \brief Compute the bounding box of the mesh
     * \author Chenshen Wu
     ***********/
    void computeBoundingBox();

    /***********************************************************************************************
     * brief Add a new weak condition to the IGA mesh
     * \param[in] _connectionID The ID of the condition
     * \param[in] _patchIndex The index of the patch in the EMPIRE data structure
     * \param[in] _patchBLIndex The index of the patch boundary loop in the EMPIRE data structure
     * \param[in] _patchBLTrCurveIndex The index of the patch trimming curve in the current boundary loop in the EMPIRE data structure
     * \param[in] _isGPprovided Flag if the GP data is provided
     * \return The pointer to the weak condition just created
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
     WeakIGADirichletCondition* addWeakDirichletCondition(int _connectionID,
                          int _patchIndex, int _patchBLIndex, int _patchBLTrCurveIndex,
                          bool _isGPProvided);

     /***********************************************************************************************
      * brief Create the GP data for the weak dirichlet condition if not provided
      * \author Andreas Apostolatos, Altug Emiroglu
      ***********/
     void createWeakDirichletConditionGPData();

    /***********************************************************************************************
     * brief Add a new weak condition to the IGA mesh
     * \param[in] _connectionID The ID of the condition
     * \param[in] _masterPatchIndex The index of the master patch in the EMPIRE data structure
     * \param[in] _masterPatchBLIndex The index of the master patch boundary loop in the EMPIRE data structure
     * \param[in] _masterPatchBLTrCurveIndex The index of the master patch trimming curve in the current boundary loop in the EMPIRE data structure
     * \param[in] _slavePatchIndex The index of the slave patch in the EMPIRE data structure
     * \param[in] _slavePatchBLIndex The index of the slave patch boundary loop in the EMPIRE data structure
     * \param[in] _slavePatchBLTrCurveIndex The index of the slave patch trimming curve in the current boundary loop in the EMPIRE data structure
     * \param[in] _isGPprovided Flag if the GP data is provided
     * \return The pointer to the weak condition just created
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
     WeakIGAPatchContinuityCondition* addWeakContinuityCondition(int _connectionID,
                          int _masterPatchIndex, int _masterPatchBLIndex, int _masterPatchBLTrCurveIndex,
                          int _slavePatchIndex,  int _slavePatchBLIndex,  int _slavePatchBLTrCurveIndex,
                          bool _isGPProvided);

     /***********************************************************************************************
      * brief Create the GP data for the weak continuity condition if not provided
      * \author Andreas Apostolatos, Altug Emiroglu
      ***********/
     void createWeakContinuityConditionGPData();

    /***********************************************************************************************
     * \brief Add coupling data of a patch
     * \param[in] _patchCounter number of the patch
     * \param[in] _BRepCounter number of the BReP
     * \param[in] _GP_m gauss point coordinates on master side
     * \param[in] _GP_s gauss point coordinates on slave side
     * \param[in] _GP_w gauss point weights
     * \param[in] _tang_m tangent on master side
     * \param[in] _tang_s tangent on slave side
     * \param[in] _map mapping from element space to physical space
     * \param[in] _ID_s ID of slave patch
     * \param[in] _NumElemsOfBRep number of linear elements on the BReP
     * \param[in] _NumGPsOfElem number of gauss points on each linear element
     * \author Ragnar Björnsson
     ***********/
    void addCouplingData(int _patchCounter, int _BRepCounter, double* _GP_m, double* _GP_s, double* _GP_w,
                         double* _tang_m, double* _tang_s, double* _map,
                         int _ID_s, int _NumElemsOfBRep, int _NumGPsOfElem);

    /// Get and set functions
public:
    /***********************************************************************************************
     * \brief Get the surface patches
     * \return A container vector of type std::vector<IGAPatchSurface*>
     * \author Fabien Pean, Chenshen Wu
     ***********/
    inline std::vector<IGAPatchSurface*> getSurfacePatches() {
		return surfacePatches;
    }
    inline const std::vector<IGAPatchSurface*>& getSurfacePatches() const {
        return surfacePatches;
    }
    /***********************************************************************************************
     * \brief Get a specific patch
     * \param[in] The id of the patch
     * \return The pointer to the patch
     * \author Fabien Pean
     ***********/
    inline IGAPatchSurface* getSurfacePatch(const unsigned int i) {
		return surfacePatches.at(i);
    }
    inline IGAPatchSurface* getSurfacePatch(const unsigned int i) const {
        return surfacePatches.at(i);
    }
    inline IGAPatchSurface* operator[](const unsigned int i) {
    	return surfacePatches.at(i);
    }
    inline const IGAPatchSurface* operator[](const unsigned int i) const {
    	return surfacePatches.at(i);
    }

    /***********************************************************************************************
     * \brief Get the number of patches
     * \author Fabien Pean
     ***********/
    inline int getNumPatches() const {
    	return surfacePatches.size();
    }

    /***********************************************************************************************
     * \brief Get the number of the Nodes
     * \param[out] The number of the Nodes
     * \author Chenshen Wu
     ***********/
    inline int getNumNodes() const {
        return numNodes;
    }

    /***********************************************************************************************
     * \brief Returns the array of all weak IGA Dirichlet conditions
     * \param[out] The array of all weak IGA Dirichlet conditions
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    std::vector<WeakIGADirichletCondition*> getWeakIGADirichletConditions() const{
        return weakIGADirichletConditions;
    }

    /***********************************************************************************************
     * \brief Returns the array of all weak IGA patch continuity conditions
     * \param[out] The array of all weak IGA patch continuity conditions
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    std::vector<WeakIGAPatchContinuityCondition*> getWeakIGAPatchContinuityConditions() const{
        return weakIGAPatchContinuityConditions;
    }

    /***********************************************************************************************
     * \brief Initilize a object of IGAPatchCouplingCaratData
     * \param[in] numPatches the number of patches
     * \param[in] numBrepsPerPatch number of couplings for each patch
     * \author Ragnar Björnsson
     ***********/
    void initializePatchCouplingData(int numPatches, int* numBrepsPerPatch) {
        couplingData = new IGAPatchCouplingCaratData(numPatches,numBrepsPerPatch);
    }

    /***********************************************************************************************
     * \brief get IGAPatchCouplingCaratData object
     * \author Ragnar Björnsson
     ***********/
    IGAPatchCouplingCaratData* getIGAPatchCouplingData() {
        return couplingData;
    }

    /***********************************************************************************************
     * \brief set clamped dofs
     * \param[in] numberOfClampedDofs the number of clamped dofs
     * \param[in] _clampedDofs the clamped dofs
     * \author Ragnar Björnsson
     ***********/
    void setClampedDofs(int numberOfClampedDofs, int* _clampedDofs) {
        for(int i = 0 ; i < numberOfClampedDofs ; i++)
            clampedDofs.push_back(_clampedDofs[i]);
    }

    /***********************************************************************************************
     * \brief set clamped dircetions
     * \param[in] _clampedDirections the lowest number of clamped directions
     * \author Ragnar Björnsson
     ***********/
    void setClampedDirections(int _clampedDirections){
        clampedDirections = _clampedDirections;
    }

    /***********************************************************************************************
     * \brief get clamped dofs
     * \author Ragnar Björnsson
     ***********/
    std::vector<int> getClampedDofs() {
        return clampedDofs;
    }

    /***********************************************************************************************
     * \brief get clamped directions
     * \author Ragnar Björnsson
     ***********/
    int getClampedDirections() {
        return clampedDirections;
    }

};

/***********************************************************************************************
 * \brief Allows for nice debug output
 * \author Fabien Pean, Chenshen Wu
 ***********/
Message &operator<<(Message &message, const IGAMesh &mesh);

}/* namespace EMPIRE */

#endif /* IGAPatchSurface_H_ */
