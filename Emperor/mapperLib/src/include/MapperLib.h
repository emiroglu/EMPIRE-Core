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
 * \date 2/22/2012
 **************************************************************************************************/

#ifndef MAPPERLIB_H_
#define MAPPERLIB_H_

#ifdef __cplusplus
extern "C" { ///Define extern C if C++ compiler is used
#endif

/***********************************************************************************************
 * \brief Initializes and inserts a MortarMapper to the mapper list
 * \param[in] mapperName name of the mapper
 *
 * \param[in] AnumNodes number of nodes on slave side
 * \param[in] AnumElems number of elements on slave side
 * \param[in] AnumNodesPerElem number of nodes per element on slave side
 * \param[in] Anodes coordinates of all nodes on slave side
 * \param[in] AnodeIDs IDs of all nodes on slave side
 * \param[in] Aelems connectivity table of all elements on slave side
 *
 * \param[in] BnumNodes number of nodes on master side
 * \param[in] BnumElems number of elements on master side
 * \param[in] BnumNodesPerElem number of nodes per element on master side
 * \param[in] Bnodes coordinates of all nodes on master side
 * \param[in] BnodeIDs IDs of all nodes on master side
 * \param[in] Belems connectivity table of all elements on master side
 *
 * \param[in] oppositeSurfaceNormal whether the interface of master side and of master side have opposite normals  or not (true or false)
 * \param[in] dual whether or not to use dual mortar (true or false)
 * \param[in] enforceConsistency whether or not to enforce consistency
 *
 * \author Altug Emiroglu
 ***********/
void init_FE_MortarMapper(char* mapperName,
                          int AnumNodes, int AnumElems, const int* AnumNodesPerElem, const double* Anodes, const int* AnodeIDs, const int* Aelems,
                          int BnumNodes, int BnumElems, const int* BnumNodesPerElem, const double* Bnodes, const int* BnodeIDs, const int* Belems,
                          int oppositeSurfaceNormal, int dual, int enforceConsistency);

/***********************************************************************************************
 * \brief Initializes and inserts a NearestNeighborMapper to the mapper list
 * \param[in] mapperName name of the mapper
 *
 * \param[in] AnumNodes number of nodes on slave side
 * \param[in] Anodes coordinates of all nodes on slave side
 *
 * \param[in] BnumNodes number of nodes on master side
 * \param[in] Bnodes coordinates of all nodes on master side
 *
 * \author Altug Emiroglu
 ***********/
void init_FE_NearestNeighborMapper(char* mapperName, 
                                   int AnumNodes, const double *Anodes,
                                   int BnumNodes, const double *Bnodes);

/***********************************************************************************************
 * \brief Initializes and inserts a NearestElementMapper to the mapper list
 * \param[in] mapperName name of the mapper
 *
 * \param[in] AnumNodes number of nodes on slave side
 * \param[in] AnumElems number of elements on slave side
 * \param[in] AnumNodesPerElem number of nodes per element on slave side
 * \param[in] Anodes coordinates of all nodes on slave side
 * \param[in] AnodeIDs IDs of all nodes on slave side
 * \param[in] Aelems connectivity table of all elements on slave side
 *
 * \param[in] BnumNodes number of nodes on master side
 * \param[in] BnumElems number of elements on master side
 * \param[in] BnumNodesPerElem number of nodes per element on master side
 * \param[in] Bnodes coordinates of all nodes on master side
 * \param[in] BnodeIDs IDs of all nodes on master side
 * \param[in] Belems connectivity table of all elements on master side
 *
 * \author Altug Emiroglu
***********/
void init_FE_NearestElementMapper(char* mapperName,
                                  int AnumNodes, int AnumElems, const int *AnumNodesPerElem, const double *Anodes, const int *AnodeIDs, const int *Aelems,
                                  int BnumNodes, int BnumElems, const int *BnumNodesPerElem, const double *Bnodes, const int *BnodeIDs, const int *Belems);

/***********************************************************************************************
 * \brief Initializes and inserts a BarycentricInterpolationMapper to the mapper list
 * \param[in] mapperName name of the mapper
 *
 * \param[in] AnumNodes number of nodes on slave side
 * \param[in] AnumElems number of elements on slave side

 * \param[in] BnumNodes number of nodes on master side
 * \param[in] BnumElems number of elements on master side
 *
 * \author Altug Emiroglu
***********/
void init_FE_BarycentricInterpolationMapper(char* mapperName, 
                                            int AnumNodes, const double *Anodes,
                                            int BnumNodes, const double *Bnodes);

/***********************************************************************************************
 * \brief Performs consistent mapping on fields (e.g. displacements or tractions) with the previously initialized mapper with name mapperName
 * \param[in] mapperName name of the mapper
 *
 * \param[in] dimension 1 or 3 dimensional data
 * \param[in] dataSize size of data for fieldA
 * \param[in] dataA the field of mesh A
 * \param[out] dataB the field of mesh B
 * \author Altug Emiroglu
***********/
void doConsistentMapping(char* mapperName, int dimension, int dataSize, const double* dataA, double* dataB);

/***********************************************************************************************
 * \brief Performs conservative mapping on fields (e.g. displacements or tractions) with the previously initialized mapper with name mapperName
 * \param[in] mapperName name of the mapper
 *
 * \param[in] dimension 1 or 3 dimensional data
 * \param[in] dataSize size of data for fieldB
 * \param[in] dataB the field of mesh B
 * \param[out] dataA the field of mesh A
 * \author Altug Emiroglu
***********/
void doConservativeMapping(char* mapperName, int dimension, int dataSize, const double* dataB, double* dataA);

/***********************************************************************************************
 * \brief Deletes a previously initialized mapper
 * \param[in] mapperName name of the mapper to be deleted from the mapper list
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
