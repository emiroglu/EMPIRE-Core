/*  Copyright &copy; 2013, TU Muenchen, Chair of Structural Analysis,
 *  Stefan Sicklinger, Tianyang Wang, Andreas Apostolatos, Altug Emiroglu Munich
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
 * \file MapperLib.h
 * This file defines the EMPIRE Mapper Library API
 * \date 5/12/2014
 **************************************************************************************************/

#ifndef MAPPERLIB_H_
#define MAPPERLIB_H_

#ifdef __cplusplus
extern "C" { ///Define extern C if C++ compiler is used
#endif

/***********************************************************************************************
 * \brief Initializes and inserts a MortarMapper to the mapper list
 * \param[in] mapperName name of the mapper
 * \param[in] AnumNodes number of nodes on slave side
 * \param[in] AnumElems number of elements on slave side
 * \param[in] AnumNodesPerElem number of nodes per element on slave side
 * \param[in] Anodes coordinates of all nodes on slave side
 * \param[in] AnodeIDs IDs of all nodes on slave side
 * \param[in] Aelems connectivity table of all elements on slave side
 * \param[in] BnumNodes number of nodes on master side
 * \param[in] BnumElems number of elements on master side
 * \param[in] BnumNodesPerElem number of nodes per element on master side
 * \param[in] Bnodes coordinates of all nodes on master side
 * \param[in] BnodeIDs IDs of all nodes on master side
 * \param[in] Belems connectivity table of all elements on master side
 * \param[in] oppositeSurfaceNormal whether the interface of master side and of master side have opposite normals  or not (true or false)
 * \param[in] dual whether or not to use dual mortar (true or false)
 * \param[in] enforceConsistency whether or not to enforce consistency
 * \author Altug Emiroglu
 ***********/
void init_FE_MortarMapper(char* mapperName,
                          int AnumNodes, int AnumElems, const int* AnumNodesPerElem, const double* Anodes, const int* AnodeIDs, const int* Aelems,
                          int BnumNodes, int BnumElems, const int* BnumNodesPerElem, const double* Bnodes, const int* BnodeIDs, const int* Belems,
                          int oppositeSurfaceNormal, int dual, int enforceConsistency);

/***********************************************************************************************
 * \brief Initializes and inserts a NearestNeighborMapper to the mapper list
 * \param[in] mapperName name of the mapper
 * \param[in] AnumNodes number of nodes on slave side
 * \param[in] Anodes coordinates of all nodes on slave side
 * \param[in] BnumNodes number of nodes on master side
 * \param[in] Bnodes coordinates of all nodes on master side
 * \author Altug Emiroglu
 ***********/
void init_FE_NearestNeighborMapper(char* mapperName, 
                                   int AnumNodes, const double *Anodes,
                                   int BnumNodes, const double *Bnodes);

/***********************************************************************************************
 * \brief Initializes and inserts a NearestElementMapper to the mapper list
 * \param[in] mapperName name of the mapper
 * \param[in] AnumNodes number of nodes on slave side
 * \param[in] AnumElems number of elements on slave side
 * \param[in] AnumNodesPerElem number of nodes per element on slave side
 * \param[in] Anodes coordinates of all nodes on slave side
 * \param[in] AnodeIDs IDs of all nodes on slave side
 * \param[in] Aelems connectivity table of all elements on slave side
 * \param[in] BnumNodes number of nodes on master side
 * \param[in] BnumElems number of elements on master side
 * \param[in] BnumNodesPerElem number of nodes per element on master side
 * \param[in] Bnodes coordinates of all nodes on master side
 * \param[in] BnodeIDs IDs of all nodes on master side
 * \param[in] Belems connectivity table of all elements on master side
 * \author Altug Emiroglu
***********/
void init_FE_NearestElementMapper(char* mapperName,
                                  int AnumNodes, int AnumElems, const int *AnumNodesPerElem, const double *Anodes, const int *AnodeIDs, const int *Aelems,
                                  int BnumNodes, int BnumElems, const int *BnumNodesPerElem, const double *Bnodes, const int *BnodeIDs, const int *Belems);

/***********************************************************************************************
 * \brief Initializes and inserts a BarycentricInterpolationMapper to the mapper list
 * \param[in] mapperName name of the mapper
 * \param[in] AnumNodes number of nodes on slave side
 * \param[in] AnumElems number of elements on slave side
 * \param[in] BnumNodes number of nodes on master side
 * \param[in] BnumElems number of elements on master side
 * \author Altug Emiroglu
***********/
void init_FE_BarycentricInterpolationMapper(char* mapperName, 
                                            int AnumNodes, const double *Anodes,
                                            int BnumNodes, const double *Bnodes);

/***********************************************************************************************
 * \brief Initializes and inserts a FEMesh to the mesh list
 * \param[in] meshName name of the mesh
 * \param[in] numNodes number of nodes
 * \param[in] numElems number of elements
 * \param[in] triangulateAll flag to determinf if elements should be triangulated
 * \author Altug Emiroglu
***********/
void initFEMesh(char* meshName, int numNodes, int numElems, bool triangulateAll);

/***********************************************************************************************
 * \brief Sets the nodes of a previously initialized FEMesh
 * \param[in] meshName name of the mesh
 * \param[in] nodeIDs node IDs
 * \param[in] nodes XYZ coordinates of the nodes with the ordering of the nodeIDs
 * \author Altug Emiroglu
***********/
void setNodesToFEMesh(char* meshName, int* nodeIDs, double* nodes);

/***********************************************************************************************
 * \brief Sets the elements of a previously initialized FEMesh
 * \param[in] meshName name of the mesh
 * \param[in] numNodesPerElem number of nodes per element
 * \param[in] elems ordered nodeIDs of each element collected into an array
 * \author Altug Emiroglu
***********/
void setElementsToFEMesh(char* meshName, int* numNodesPerElem, int* elems);

/***********************************************************************************************
 * \brief Initializes and inserts a IGAMesh to the mesh list
 * \param[in] meshName name of the mesh
 * \param[in] numNodes total number of CPs in a multipatch geometry
 * \author Altug Emiroglu
***********/
void initIGAMesh(char* meshName, int numNodes);

/***********************************************************************************************
 * brief Adds a new patch patch to a previously initialized IGA mesh
 * \param[in] meshName name of the mesh
 * \param[in] pDegree The polynomial degree of the IGA 2D patch in the u-direction
 * \param[in] uNoKnots The number of knots for the knot vector in the u-direction
 * \param[in] uKnotVector The underlying knot vector of the IGA 2D patch in the u-direction
 * \param[in] qDegree The polynomial degree of the IGA 2D patch in the v-direction
 * \param[in] vNoKnots The number of knots for the knot vector in the v-direction
 * \param[in] vKnotVector The underlying knot vector of the IGA 2D patch in the v-direction
 * \param[in] uNoControlPoints The number of the Control Points for the 2D NURBS patch in the u-direction
 * \param[in] vNoControlPoints The number of the Control Points for the 2D NURBS patch in the v-direction
 * \param[in] controlPointNet The set of the Control Points related to the 2D NURBS patch
 * \param[in] dofIndexNet The index of the dof of the each Control Points related to. Equivalent to equation ID
 * \author Altug Emiroglu
 ***********/
void addPatchToIGAMesh(char* meshName,
                       int pDegree, int uNoKnots, double* uKnotVector,
                       int qDegree, int vNoKnots, double* vKnotVector,
                       int uNoControlPoints, int vNoControlPoints,
                       double* controlPointNet, int* dofIndexNet);

/***********************************************************************************************
 * \brief Adds trimming loop information to a previously added patch of a mesh
 * \param[in] meshName name of the mesh
 * \param[in] idxSurfacePatch index of the surface patch that is previously added to the mesh
 * \param[in] inner 0 for outter and 1 for inner
 * \param[in] numCurves Number of curves to be received for this loop
 * \author Altug Emiroglu
 ***********/
void addTrimmingLoopToPatch(char* meshName, int idxSurfacePatch,
                            int inner, int numCurves);

/***********************************************************************************************
 * \brief Adds a NURBS curve for the current loop and its attached information.
 *        It always adds the curve to the last initialized trimming loop
 * \param[in] direction The direction of the curve according to its CP and knot vector ordering
 * \param[in] pDegree The polynomial degree of the IGA 1D curve in the u-direction
 * \param[in] uNoKnots The number of knots for the knot vector in the u-direction
 * \param[in] uKnotVector The underlying knot vector of the IGA 1D curve in the u-direction
 * \param[in] uNoControlPoints The number of the Control Points for the 1D NURBS patch in the u-direction
 * \param[in] controlPoints The set of the Control Points related to the 1D NURBS curve
 * \author Altug Emiroglu
 ***********/
void addTrimmingCurveToTrimmingLoop(char* meshName, int idxSurfacePatch,
                                    int direction, int pDegree, int uNoKnots, double* uKnotVector,
                                    int uNoControlPoints, double* controlPoints);

/***********************************************************************************************
 * \brief Initializes and inserts a MortarMapper to the mapper list
 * \param[in] mapperName name of the mapper
 * \param[in] AmeshName a previously initialized slave FEMesh name
 * \param[in] BmeshName a previously initialized master FEMesh name
 * \author Altug Emiroglu
 ***********/
void initFEMNearestNeighborMapper(char* mapperName,
                                  char* AmeshName, char* BmeshName);

/***********************************************************************************************
 * \brief Initializes and inserts a MortarMapper to the mapper list
 * \param[in] mapperName name of the mapper
 * \param[in] AmeshName a previously initialized slave FEMesh name
 * \param[in] BmeshName a previously initialized master FEMesh name
 * \author Altug Emiroglu
 ***********/
void initFEMNearestElementMapper(char* mapperName,
                                 char* AmeshName, char* BmeshName);

/***********************************************************************************************
 * \brief Initializes and inserts a MortarMapper to the mapper list
 * \param[in] mapperName name of the mapper
 * \param[in] AmeshName a previously initialized slave FEMesh name
 * \param[in] BmeshName a previously initialized master FEMesh name
 * \author Altug Emiroglu
 ***********/
void initFEMBarycentricInterpolationMapper(char* mapperName,
                                           char* AmeshName, char* BmeshName);

/***********************************************************************************************
 * \brief Initializes and inserts a MortarMapper to the mapper list
 * \param[in] mapperName name of the mapper
 * \param[in] AmeshName a previously initialized slave FEMesh name
 * \param[in] BmeshName a previously initialized master FEMesh name
 * \param[in] oppositeSurfaceNormal whether the master side and slave side have opposite normals(true or false)
 * \param[in] dual whether or not to use dual mortar (true or false)
 * \param[in] enforceConsistency whether or not to enforce consistency
 * \author Altug Emiroglu
 ***********/
void initFEMMortarMapper(char* mapperName,
                         char* AmeshName, char* BmeshName,
                         int oppositeSurfaceNormal, int dual, int enforceConsistency);

/***********************************************************************************************
 * \brief Initializes and inserts a MortarMapper to the mapper list
 * \param[in] mapperName name of the mapper
 * \param[in] IGAMeshName a previously initialized IGAMeshName name
 * \param[in] FEMeshName a previously initialized FEMeshName name
 * \param[in] isMappingIGA2FEM direction of the mapping to be generated (true=IGA->FEM or false=FEM->IGA)
 * \author Altug Emiroglu
 ***********/
void initIGAMortarMapper(char* mapperName, char* IGAMeshName, char* FEMeshName, bool isMappingIGA2FEM);

/***********************************************************************************************
 * \brief Set parameter for the projection of mesh onto the NURBS surface
 * \param[in] mapperName name of the mapper
 * \param[in] maxProjectionDistance The max distance allowed between FE mesh and NURBS surface
 * \param[in] numRefinementForIntialGuess The number of test point to find initial guess for Newton-Raphson scheme
 * \param[in] maxDistanceForProjectedPointsOnDifferentPatches The max authorized distance between two projected points from a same physical node
 * \author Altug Emiroglu
 ***********/
void setParametersProjection(char* mapperName,
                             double maxProjectionDistance, int numRefinementForIntialGuess,
                             double maxDistanceForProjectedPointsOnDifferentPatches);

/***********************************************************************************************
 * \brief Set parameter for Newton-Raphson scheme of projection on NURBS patch
 * \param[in] mapperName name of the mapper
 * \param[in] newtonRaphsonMaxIt The number of iteration for Newton-Raphson scheme of projecting a node on a NURBS patch
 * \param[in] newtonRaphsonTol The tolerance for Newton-Raphson scheme of projecting a node on a NURBS patch
* \author Altug Emiroglu
 ***********/
void setParametersNewtonRaphson(char* mapperName,
                                int maxNumOfIterations, double tolerance);

/***********************************************************************************************
 * \brief Set parameter for Newton-Raphson scheme of projection on NURBS patch boundary
 * \param[in] mapperName name of the mapper
 * \param[in] newtonRaphsonBoundaryMaxIt The number of iteration for Newton-Raphson scheme of projecting a node on a NURBS patch boundary
 * \param[in] newtonRaphsonBoundaryTol The tolerance for Newton-Raphson scheme of projecting a node on a NURBS patch boundary
 * \author Altug Emiroglu
 ***********/
void setParametersNewtonRaphsonBoundary(char* mapperName,
                                        int maxNumOfIterations, double tolerance);

/***********************************************************************************************
 * \brief Set parameter for bisection scheme of projection on NURBS patch boundary
 * \param[in] mapperName name of the mapper
 * \param[in] bisectionMaxIt The number of iteration for bisection scheme of projecting a node on a NURBS patch boundary
 * \param[in] bisectionTol The tolerance for bisection scheme of projecting a node on a NURBS patch boundary
 * \author Altug Emiroglu
 ***********/
void setParametersBisection(char* mapperName,
                            int maxNumOfIterations, double tolerance);

/***********************************************************************************************
 * \brief Set parameter for integration
 * \param[in] mapperName name of the mapper
 * \param[in] numGPsTriangle The number of Gauss points when performs integration on triangle
 * \param[in] numGPsQuad The number of Gauss points when performs integration on quadrilateral
 * \author Altug Emiroglu
 ***********/
void setParametersIntegration(char* mapperName,
                              int numGPTriangle, int numGPQuad);

/***********************************************************************************************
 * \brief Build Coupling Matrices
 * \param[in] mapperName name of the mapper
 * \author Altug Emiroglu
 ***********/
void buildCouplingMatrices(char* mapperName);

/***********************************************************************************************
 * \brief Performs consistent mapping on fields (e.g. displacements or tractions) with the previously initialized mapper with name mapperName
 * \param[in] mapperName name of the mapper
 *
 * \param[in] dimension 1 or 3 dimensional data
 * \param[in] dataSizeA size of data for fieldA
 * \param[in] dataA the field of mesh A
 * \param[in] dataSizeB size of data for fieldB
 * \param[out] dataB the field of mesh B
 * \author Altug Emiroglu
***********/
void doConsistentMapping(char* mapperName, int dimension, int dataSizeA, const double* dataA, int dataSizeB, double* dataB);

/***********************************************************************************************
 * \brief Performs conservative mapping on fields (e.g. displacements or tractions) with the previously initialized mapper with name mapperName
 * \param[in] mapperName name of the mapper
 *
 * \param[in] dimension 1 or 3 dimensional data
 * \param[in] dataSizeB size of data for fieldB
 * \param[in] dataB the field of mesh B
 * \param[in] dataSizeA size of data for fieldA
 * \param[out] dataA the field of mesh A
 * \author Altug Emiroglu
***********/
void doConservativeMapping(char* mapperName, int dimension, int dataSizeB, const double* dataB, int dataSizeA, double* dataA);

/***********************************************************************************************
 * \brief Deletes a previously initialized mapper
 * \param[in] mapperName name of the mapper to be deleted from the mapper list
 * \author Altug Emiroglu
 ***********/
void deleteMapper(char* mapperName);

/***********************************************************************************************
 * \brief Delets all the initialized mappers in the mapper list
 ***********/
void deleteAllMappers();

#ifdef __cplusplus
}
#endif

#endif // MAPPERLIB_H_
