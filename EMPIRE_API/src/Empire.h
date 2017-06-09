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
/***********************************************************************************************//**
 * \file Empire.h
 * This file holds the class of Empire for the API
 * \date 3/18/2012
 **************************************************************************************************/

#ifndef EMPIRE_H_
#define EMPIRE_H_

#include <string>
#include <vector>
// how do I include cfloat, #include <cfloat> doesnt seem to work
// check ManagePortability.h in carat

namespace EMPIRE {
/********//**
 * \brief This class is the main class for EMPIRE (Client = API)
 *
 * \author Stefan Sicklinger
 ***********/
class Empire {
public:
    /***********************************************************************************************
     * \brief Constructor
     *
     * \author Stefan Sicklinger
     ***********/
    Empire();
    /***********************************************************************************************
     * \brief Destructor
     *
     * \author Stefan Sicklinger
     ***********/
    virtual ~Empire();
    /***********************************************************************************************
     * \brief connects API with Emperor
     *
     * \author Stefan Sicklinger
     ***********/
    void connect();
    /***********************************************************************************************
     * \brief disconnects API from Emperor
     *
     * \author Stefan Sicklinger
     ***********/
    void disconnect();
    /***********************************************************************************************
     * \brief Initializes Meta-database (parsing done) and ClientCommunication
     *
     * \author Stefan Sicklinger
     ***********/
    void initEnvironment(char *inputFileName);
    /***********************************************************************************************
     * \brief Get user defined text by the element name in the XML input file
     * \param[in] elementName name of the XML element
     * \return user defined text
     * \author Tianyang Wang
     ***********/
    std::string getUserDefinedText(std::string elementName);
    /***********************************************************************************************
     * \brief Send the mesh to the server
     * \param[in] numNodes number of nodes
     * \param[in] numElems number of elements
     * \param[in] nodes coordinates of all nodes
     * \param[in] nodeIDs IDs of all nodes
     * \param[in] numNodesPerElem number of nodes per element
     * \param[in] elems connectivity table of all elements
     * \author Tianyang Wang
     ***********/
    void sendMesh(int numNodes, int numElems, double *nodes, int *nodeIDs, int *numNodesPerElem,
            int *elems);

    /***********************************************************************************************
     * \brief Send the section mesh to the server (for mapping with beam elements)
     * \param[in] name name of the mesh
     * \param[in] numNodes number of nodes
     * \param[in] numElems number of elements
     * \param[in] nodes coordinates of all nodes
     * \param[in] nodeIDs IDs of all nodes
     * \param[in] numNodesPerElem number of nodes per element
     * \param[in] elems connectivity table of all elements
     * \param[in] numSections number of sections
     * \param[in] numRootSectionNodes number of nodes of the root section
     * \param[in] numNormalSectionNodes number of nodes of every normal section
     * \param[in] numTipSectionNodes number of nodes of the tip section
     * \param[in] rotationGlobal2Root rotation matrix from the global coordinate system to the root section system
     * \param[in] translationGlobal2Root translation vector from the global coordinate system to the root section system
     ***********/
    void sendSectionMesh(int numNodes, int numElems, double *nodes, int *nodeIDs,
            int *numNodesPerElem, int *elems, int numSections, int numRootSectionNodes,
            int numNormalSectionNodes, int numTipSectionNodes, double *rotationGlobal2Root,
            double *translationGlobal2Root);

    /***********************************************************************************************
     * \brief Receive mesh from the server
     * \param[in] numNodes number of nodes
     * \param[in] numElems number of elements
     * \param[in] nodes coordinates of all nodes
     * \param[in] nodeIDs IDs of all nodes
     * \param[in] numNodesPerElem number of nodes per element
     * \param[in] elems connectivity table of all elements
     * \author Altug Emiroglu
     ***********/
    void recvMesh(int *numNodes, int *numElems, double **nodes, int **nodeIDs,
            int **numNodesPerElem, int **elems);

    /***********************************************************************************************
     * \brief Send the IGA patch to the server
     * \param[in] _pDegree The polynomial degree of the IGA 2D patch in the u-direction
     * \param[in] _uNumKnots The number of knots for the knot vector in the u-direction
     * \param[in] _uKnotVector The underlying knot vector of the IGA 2D patch in the u-direction
     * \param[in] _qDegree The polynomial degree of the IGA 2D patch in the v-direction
     * \param[in] _vNumKnots The number of knots for the knot vector in the v-direction
     * \param[in] _vKnotVector The underlying knot vector of the IGA 2D patch in the v-direction
     * \param[in] _uNumControlPoints The number of the Control Points for the 2D NURBS patch in the u-direction
     * \param[in] _vNumControlPoints The number of the Control Points for the 2D NURBS patch in the v-direction
     * \param[in] _cpNet The set of the Control Points related to the 2D NURBS patch
     * \param[in] _nodeNet The set of the dof index Control Points related to the 2D NURBS patch
     * \author Chenshen Wu 
     ***********/
    void sendIGAPatch(int _pDegree, int _uNumKnots, double* _uKnotVector, int _qDegree,
            int _vNumKnots, double* _vKnotVector, int _uNumControlPoints, int _vNumControlPoints,
            double* _cpNet, int* _nodeNet);

    /***********************************************************************************************
     * \brief Send the IGA mesh to the server
     * \param[in] _numPatches The number of the patches out of which the IGA mesh consists
     * \param[in] _numNodes The number of the Control Points which are needed for the computation of the coupling matrices
     * \author Chenshen Wu
     ***********/
    void sendIGAMesh(int _numPatches, int _numNodes);

    /***********************************************************************************************
     * \brief Send the IGA trimming information to the server
     * \param[in] _isTrimmed Whether the current considered patch is trimmed
     * \param[in] _numLoops The number of loops defining boundary
     * \author Fabien Pean
     ***********/
    void sendIGATrimmingInfo(int _isTrimmed, int _numLoops);

    /***********************************************************************************************
     * \brief Send the IGA trimming information about the loop to the server
     * \param[in] _inner whether loop is outter boundary loop or inner
     * \param[in] _numCurves The number of curves defining the loop
     * \author Fabien Pean
     ***********/
    void sendIGATrimmingLoopInfo(int _inner, int _numCurves);

    /***********************************************************************************************
     * \brief Send the IGA trimming curve to the server
     * \param[in] direction The direction of the curve if is following standard or not
     * \param[in] _pDegree The polynomial degree of the IGA 1D curve in the u-direction
     * \param[in] _uNumKnots The number of knots for the knot vector in the u-direction
     * \param[in] _uKnotVector The underlying knot vector of the IGA 1D curve in the u-direction
     * \param[in] _uNoControlPoints The number of the Control Points for the 1D NURBS patch in the u-direction
     * \param[in] _controlPointNet The set of the Control Points related to the 1D NURBS patch
     * \author Fabien Pean
     ***********/
    void sendIGATrimmingCurve(int _direction, int _pDegree, int _uNumKnots, double* _uKnotVector,
            int _uNumControlPoints, double* _cpNet);

    /***********************************************************************************************
     * \brief Send the coupling information of the current connection to the server
     * \param[in] _numConnections Number of connections between patches in the multipatch geometry
     * \author Andreas Apostolatos, Altug Emiroglu
     ************/
    void sendIGANumPatchConnections(int _numConnections);

    /***********************************************************************************************
     * \brief Send the coupling information of the current connection to the server
     * \param[in] _masterPatchIndex The index of the master patch in the EMPIRE data structure
     * \param[in] _masterPatchBLIndex The index of the master patch boundary loop in the EMPIRE data structure
     * \param[in] _masterPatchBLTrCurveIndex The index of the master patch trimming curve in the current boundary loop in the EMPIRE data structure
     * \param[in] _slavePatchIndex The index of the slave patch in the EMPIRE data structure
     * \param[in] _slavePatchBLIndex The index of the slave patch boundary loop in the EMPIRE data structure
     * \param[in] _slavePatchBLTrCurveIndex The index of the slave patch trimming curve in the current boundary loop in the EMPIRE data structure
     * \param[in] _isGPprovided Flag if the GP data is provided
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    void sendIGAPatchConnectionInfo(int _masterPatchIndex, int _masterPatchBLIndex, int _masterPatchBLTrCurveIndex,
                                    int _slavePatchIndex, int _slavePatchBLIndex, int _slavePatchBLTrCurveIndex,
                                    int _isGPProvided);

    /***********************************************************************************************
     * \brief Send the coupling information of the current connection to the server
     * \param[in] _trCurveNumGP The total number of GPs on the trimming curve
     * \param[in] _trCurveMasterGPs The parametric coordinates of the GPs on the master trimming curve in the master patch parameter space
     * \param[in] _trCurveSlaveGPs The parametric coordinates of the GPs on the slave trimming curve in the slave patch parameter space
     * \param[in] _trCurveGPWeights The GP weights
     * \param[in] _trCurveMasterGPTangents The tangent to the trimming curve vector in the parameter space of the master patch (third coordinate is 0 or unused)
     * \param[in] _trCurveSlaveGPTangents The tangent to the trimming curve vector in the parameter space of the slave patch (third coordinate is 0 or unused)
     * \param[in] _trCurveGPJacobianProducts The Jacobian products for the transformation of the integrals
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    void sendIGAPatchConnectionData(int _trCurveNumGP,
                                    double* _trCurveMasterGPs, double* _trCurveSlaveGPs, double* _trCurveGPWeights,
                                    double* _trCurveMasterGPTangents, double* _trCurveSlaveGPTangents,
                                    double* _trCurveGPJacobianProducts);

    /***********************************************************************************************
     * \brief Send IGAPatchCoupling to the server
     * \param[in] numPatches number of patches
     * \param[in] numBRepsPerPatch number of BReps each patch has
     * \param[in] totalNumGP total number of gausspoints on all patch coupling interfaces
     * \param[in] totalNumBRePs total number of BReps
     * \param[in] allID_slave all ids of the slave patches
     * \param[in] allNumElemsPerBRep number of elements for all BReps
     * \param[in] allNumGPsPerElem number of gauss points for each elem for all BReps
     * \param[in] allGPOfBRep_master the parametric coordinates of the gauss points on the master patch
     * \param[in] allGPOfBRep_slave the parametric coordinates of the gauss points on the slave patch
     * \param[in] allGPOfBRep_weight the weight of the gauss points
     * \param[in] allTangents_master all the tangents on the cartesian space on the master patch
     * \param[in] allTangents_slave all the tangents on the cartesian space on the slave patch
     * \param[in] allMappings the mapping of each gauss point from the parent element space to the cartesian space
     * \author Ragnar Björnsson
     ***********/
    void sendIGAPatchCouplingGaussPointsTest(int numPatches, int* numBRepsPerPatch, int totalNumGP, int totalNumBRePs,
            int* allID_slave, int* allNumElemsPerBRep, int* allNumGPsPerElem,
            double* allGPOfBRep_master, double* allGPOfBRep_slave, double* allGPOfBRep_weight,
            double* allTangents_master, double* allTangents_slave, double* allMappings);

    /***********************************************************************************************
     * \brief Send IGA dirichlet boundary conditions to the server
     * \param[in] numberOfClampedDofs total number of clamped IGA dofs
     * \param[in] clampedDofs clamped IGA dofs
     * \param[in] clampedDirections lowest number of clamped directions for all clamped nodes
     * \author Ragnar Björnsson
     ***********/
    void sendIGADirichletDofs(int numberOfClampedDofs, int* clampedDofs, int clampedDirections);

    /***********************************************************************************************
     * \brief Send data field to the server
     * \param[in] sizeOfArray size of the array (data field)
     * \param[in] dataField the data field to be sent
     * \author Tianyang Wang
     ***********/
    void sendDataField(int sizeOfArray, double *dataField);
    /***********************************************************************************************
     * \brief Receive data field from the server
     * \param[in] sizeOfArray size of the array (data field)
     * \param[out] dataField the data field to be received
     * \author Tianyang Wang
     ***********/
    void recvDataField(int sizeOfArray, double *dataField);
    /***********************************************************************************************
     * \brief Send signal to the server
     * \param[in] name name of the signal
     * \param[in] sizeOfArray size of the array (signal)
     * \param[in] signal the signal
     ***********/
    void sendSignal_double(char *name, int sizeOfArray, double *signal);
    /***********************************************************************************************
     * \brief Receive signal from the server
     * \param[in] name name of the signal
     * \param[in] sizeOfArray size of the array (signal)
     * \param[in] signal the signal
     ***********/
    void recvSignal_double(char *name, int sizeOfArray, double *signal);
    /***********************************************************************************************
     * \brief Send the convergence signal of an loop
     * \param[in] signal 1 means convergence, 0 means non-convergence
     ***********/
    void sendConvergenceSignal(int signal);
    /***********************************************************************************************
     * \brief Receive the convergence signal of an loop
     * \return 1 means convergence, 0 means non-convergence
     * \author Tianyang Wang
     ***********/
    int recvConvergenceSignal();
    /***********************************************************************************************
     * \brief A simple debug function showing the content of the data field
     * \param[in] name name of the data field
     * \param[in] sizeOfArray size of the array (data field)
     * \param[in] dataField the data field to be printed
     * \author Tianyang Wang
     ***********/
    void printDataField(char *name, int sizeOfArray, double *dataField);
};

}/* namespace EMPIRE */

#endif /* EMPIRE_H_ */
