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
    IGAPatchCouplingCaratData* couplingData = NULL;

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
