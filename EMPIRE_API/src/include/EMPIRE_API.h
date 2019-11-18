/*           This file has been prepared for Doxygen automatic documentation generation.          */
/***********************************************************************************************//**
 * \mainpage
 *
 * \section LICENSE
 *  Copyright &copy; 2013, TU Muenchen, Chair of Structural Analysis, Stefan Sicklinger, Tianyang Wang, Munich \n
 *  All rights reserved. \n
 *
 *  This file is part of EMPIRE.
 *
 *  EMPIRE is free software: you can redistribute it and/or modify \n
 *  it under the terms of the GNU General Public License as published by \n
 *  the Free Software Foundation, either version 3 of the License, or \n
 *  (at your option) any later version. \n
 *
 *  EMPIRE is distributed in the hope that it will be useful, \n
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of \n
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the \n
 *  GNU General Public License for more details. \n
 *
 *  You should have received a copy of the GNU General Public License \n
 *  along with EMPIRE.  If not, see <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses/</a>.
 *
 *
 * \section DESCRIPTION
 * This is the API of EMPIRE. It consists of dynamic link library and a header file called EMPIRE_API.h
 *
 *
 * \section COMPILATION
 *  export CC=icc
 *  export CXX=icpc
 *  cd build
 *  cmake ..
 * There are the following make targets available:
 * - make (compilation and linking)
 * - make clean (remove object files and executable including all folders)
 * - make doc (generates documentation) html main file is  /EMPEROR/doc/html/index.html
 * - make cleandoc (removes documentation)
 *
 *
 * \section HOWTO
 * Please find all further information on
 * <a href="http://empire.st.bv.tum.de">EMPIRE Project</a>
 *
 *
 * <EM> Note: The Makefile suppresses per default all compile and linking command output to the terminal.
 *       You may enable this information by make VEREBOSE=1</EM>
 *
 **************************************************************************************************/

/***********************************************************************************************//**
 * \file EMPIRE_API.h
 * This file defines the EMPIRE API for Co-Simulation
 * \date 2/22/2012
 **************************************************************************************************/

#ifndef EMPIRE_API_H_
#define EMPIRE_API_H_

#ifdef __cplusplus
extern "C" { ///Define extern C if C++ compiler is used
#endif
/***********************************************************************************************
 * \brief Establishes the necessary connection with the Emperor
 ***********/
void EMPIRE_API_Connect(char* inputFileName);

/***********************************************************************************************
 * \brief Get user defined text by the element name in the XML input file
 * \param[in] elementName name of the XML element
 * \return user defined text
 ***********/
const char *EMPIRE_API_getUserDefinedText(char *elementName);

/***********************************************************************************************
 * \brief Send the mesh to the server
 * \param[in] name name of the mesh
 * \param[in] numNodes number of nodes
 * \param[in] numElems number of elements
 * \param[in] nodes coordinates of all nodes
 * \param[in] nodeIDs IDs of all nodes
 * \param[in] numNodesPerElem number of nodes per element
 * \param[in] elems connectivity table of all elements
 ***********/
void EMPIRE_API_sendMesh(char *name, int numNodes, int numElems, double *nodes, int *nodeIDs,
        int *numNodesPerElem, int *elems);

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
void EMPIRE_API_sendSectionMesh(char *name, int numNodes, int numElems, double *nodes, int *nodeIDs,
        int *numNodesPerElem, int *elems, int numSections, int numRootSectionNodes,
        int numNormalSectionNodes, int numTipSectionNodes, double *rotationGlobal2Root,
        double *translationGlobal2Root);

/***********************************************************************************************
 * \brief Recieve mesh from the server
 * \param[in] name name of the mesh
 * \param[in] numNodes number of nodes
 * \param[in] numElems number of elements
 * \param[in] nodes coordinates of all nodes
 * \param[in] nodeIDs IDs of all nodes
 * \param[in] numNodesPerElem number of nodes per element
 * \param[in] elems connectivity table of all elements
 ***********/
void EMPIRE_API_recvMesh(char *name, int *numNodes, int *numElems, double **nodes, int **nodeIDs,
        int **numNodesPerElem, int **elem);

/***********************************************************************************************
 * \brief Send the IGA patch to the server
 * \param[in] _name name of the field
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
 ***********/
void EMPIRE_API_sendIGAPatch(int _pDegree, int _uNumKnots, double* _uKnotVector, int _qDegree,
        int _vNumKnots, double* _vKnotVector, int _uNumControlPoints, int _vNumControlPoints,
        double* _cpNet, int* _nodeNet);

/***********************************************************************************************
 * \brief Send the IGA patch to the server
 * \param[in] _name name of the field
 * \param[in] _numPatches The number of the patches contained in the IGA mesh
 * \param[in] _numNodes The number of nodes of the analysis model, i.e. merged Control Points are seen as one node
 ***********/
void EMPIRE_API_sendIGAMesh(char *_name, int _numPatches, int _numNodes);

/***********************************************************************************************
 * \brief Send the IGA trimming information to the server
 * \param[in] _isTrimmed Whether the current considered patch is trimmed
 * \param[in] _numLoops The number of loops defining boundary
 * \author Fabien Pean
 ***********/
void EMPIRE_API_sendIGATrimmingInfo(int _isTrimmed, int _numLoops);

/***********************************************************************************************
 * \brief Send the IGA trimming information about the loop to the server
 * \param[in] _inner whether loop is outter boundary loop or inner
 * \param[in] _numCurves The number of curves defining the loop
 * \author Fabien Pean
 ***********/
void EMPIRE_API_sendIGATrimmingLoopInfo(int _inner, int _numCurves);

/***********************************************************************************************
 * \brief Send a IGA trimming curve to the server
 * \param[in] direction The direction of the curve if is following standard or not
 * \param[in] _pDegree The polynomial degree of the IGA 1D curve in the u-direction
 * \param[in] _uNumKnots The number of knots for the knot vector in the u-direction
 * \param[in] _uKnotVector The underlying knot vector of the IGA 1D curve in the u-direction
 * \param[in] _uNoControlPoints The number of the Control Points for the 1D NURBS patch in the u-direction
 * \param[in] _controlPointNet The set of the Control Points related to the 1D NURBS patch
 * \author Fabien Pean
 ***********/
void EMPIRE_API_sendIGATrimmingCurve(int _direction, int _pDegree, int _uNumKnots,
        double* _uKnotVector, int _uNumControlPoints, double* _cpNet);

/***********************************************************************************************
 * \brief Send the information of the current Dirichlet condition to the server
 * \param[in] _numDirichletConditions Number of connections between patches in the multipatch geometry
 * \author Andreas Apostolatos, Altug Emiroglu
 ************/
void EMPIRE_API_sendIGANumDirichletConditions(int _numDirichletConditions);

/***********************************************************************************************
 * \brief Send the coupling information of the current connection to the server
 * \param[in] _patchCtr The index of the patch in the EMPIRE data structure
 * \param[in] _patchBLCtr The index of the patch boundary loop in the EMPIRE data structure
 * \param[in] _patchBLTrCurveCtr The index of the patch trimming curve in the current boundary loop in the EMPIRE data structure
 * \param[in] _isGPprovided Flag if the GP data is provided
 * \author Andreas Apostolatos, Altug Emiroglu
 ***********/
void EMPIRE_API_sendIGADirichletConditionInfo(int _patchCtr, int _patchBLCtr, int _patchBLTrCurveCtr,
                                           int _isGPProvided);

/***********************************************************************************************
 * \brief Send the coupling information of the current connection to the server
 * \param[in] _trCurveNumGP The total number of GPs on the trimming curve
 * \param[in] _trCurveGPs The parametric coordinates of the GPs on the master trimming curve in the patch parameter space
 * \param[in] _trCurveGPWeights The GP weights
 * \param[in] _trCurveGPTangents The tangent to the trimming curve vector in the parameter space of the patch
 * \param[in] _trCurveGPJacobianProducts The Jacobian products for the transformation of the integrals
 * \author Andreas Apostolatos, Altug Emiroglu
 ***********/
void EMPIRE_API_sendIGADirichletConditionData(int _trCurveNumGP,
                                           double* _trCurveGPs, double* _trCurveGPWeights,
                                           double* _trCurveGPTangents,
                                           double* _trCurveGPJacobianProducts);

/***********************************************************************************************
 * \brief Send the coupling information of the current connection to the server
 * \param[in] _numConnections Number of connections between patches in the multipatch geometry
 * \author Andreas Apostolatos, Altug Emiroglu
 ************/
void EMPIRE_API_sendIGANumPatchConnections(int _numConnections);

/***********************************************************************************************
 * \brief Send the coupling information of the current connection to the server
 * \param[in] _masterPatchCtr The index of the master patch in the EMPIRE data structure
 * \param[in] _masterPatchBLCtr The index of the master patch boundary loop in the EMPIRE data structure
 * \param[in] _masterPatchBLTrCurveCtr The index of the master patch trimming curve in the current boundary loop in the EMPIRE data structure
 * \param[in] _slavePatchCtr The index of the slave patch in the EMPIRE data structure
 * \param[in] _slavePatchBLCtr The index of the slave patch boundary loop in the EMPIRE data structure
 * \param[in] _slavePatchBLTrCurveCtr The index of the slave patch trimming curve in the current boundary loop in the EMPIRE data structure
 * \param[in] _isGPprovided Flag if the GP data is provided
 * \author Andreas Apostolatos, Altug Emiroglu
 ***********/
void EMPIRE_API_sendIGAPatchConnectionInfo(int _masterPatchCtr, int _masterPatchBLCtr, int _masterPatchBLTrCurveCtr,
                                           int _slavePatchCtr, int _slavePatchBLCtr, int _slavePatchBLTrCurveCtr,
                                           int _isGPProvided);

/***********************************************************************************************
 * \brief Send the coupling information of the current connection to the server
 * \param[in] _trCurveNumGP The total number of GPs on the trimming curve
 * \param[in] _trCurveMasterGPs The parametric coordinates of the GPs on the master trimming curve in the master patch parameter space
 * \param[in] _trCurveSlaveGPs The parametric coordinates of the GPs on the slave trimming curve in the slave patch parameter space
 * \param[in] _trCurveGPWeights The GP weights
 * \param[in] _trCurveMasterGPTangents The tangent to the trimming curve vector in the parameter space of the master patch
 * \param[in] _trCurveSlaveGPTangents The tangent to the trimming curve vector in the parameter space of the slave patch
 * \param[in] _trCurveGPJacobianProducts The Jacobian products for the transformation of the integrals
 * \author Andreas Apostolatos, Altug Emiroglu
 ***********/
void EMPIRE_API_sendIGAPatchConnectionData(int _trCurveNumGP,
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
void EMPIRE_API_sendIGAPatchCouplingGaussPointsTest(int numPatches, int* numBRepsPerPatch, int totalNumGP, int totalNumBRePs,
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
void EMPIRE_API_sendIGADirichletDofs(int numberOfClampedDofs, int* clampedDofs, int clampedDirections);

/***********************************************************************************************
 * \brief Send data field to the server
 * \param[in] name name of the field
 * \param[in] sizeOfArray size of the array (data field)
 * \param[in] dataField the data field to be sent
 ***********/
void EMPIRE_API_sendDataField(char *name, int sizeOfArray, double *dataField);

/***********************************************************************************************
 * \brief Receive data field from the server
 * \param[in] name name of the field
 * \param[in] sizeOfArray size of the array (data field)
 * \param[out] dataField the data field to be received
 ***********/
void EMPIRE_API_recvDataField(char *name, int sizeOfArray, double *dataField);

/***********************************************************************************************
 * \brief Send signal to the server
 * \param[in] name name of the signal
 * \param[in] sizeOfArray size of the array (signal)
 * \param[in] signal the signal
 ***********/
void EMPIRE_API_sendSignal_double(char *name, int sizeOfArray, double *signal);

/***********************************************************************************************
 * \brief Receive signal from the server
 * \param[in] name name of the signal
 * \param[in] sizeOfArray size of the array (signal)
 * \param[in] signal the signal
 ***********/
void EMPIRE_API_recvSignal_double(char *name, int sizeOfArray, double *signal);

/***********************************************************************************************
 * \brief Receive the convergence signal of an loop
 * \return 1 means convergence, 0 means non-convergence
 ***********/
int EMPIRE_API_recvConvergenceSignal();

/***********************************************************************************************
 * \brief Send the convergence signal of an loop
 * \param[in] signal 1 means convergence, 0 means non-convergence
 ***********/
void EMPIRE_API_sendConvergenceSignal(int signal);

/***********************************************************************************************
 * \brief A simple debug function showing the content of the data field
 * \param[in] name name of the data field
 * \param[in] sizeOfArray size of the array (data field)
 * \param[in] dataField the data field to be printed
 ***********/
void EMPIRE_API_printDataField(char *name, int sizeOfArray, double *dataField);

/***********************************************************************************************
 * \brief Performs disconnection and finalization operations to the Emperor
 ***********/
void EMPIRE_API_Disconnect(void);

/// maximum length of a name string
static const int EMPIRE_API_NAME_STRING_LENGTH = 80;

#ifdef __cplusplus
}
#endif

#endif /* EMPIRE_API_H_ */
