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
 * \file MapperLib.cpp
 * This file wraps the C++ API in a C API
 * \date 5/12/2014
 **************************************************************************************************/
#include "AbstractMapper.h"
#include "NearestNeighborMapper.h"
#include "NearestElementMapper.h"
#include "BarycentricInterpolationMapper.h"
#include "MortarMapper.h"
#include "MapperLib.h"
#include <assert.h>
#include <string>

using namespace EMPIRE;
using namespace std;

/// Lets make an AbstractMapper object map in the global scope
std::map <std::string, AbstractMapper*> mapperList;

void init_FE_NearestNeighborMapper(char* mapperName, 
                                   int AnumNodes, const double *Anodes, int BnumNodes, const double *Bnodes){
    std::string mapperNameToMap = std::string(mapperName);

    if (mapperList.count( mapperNameToMap )){
        ERROR_OUT("A mapper with name : " + mapperNameToMap + "has already been initialized!");
        exit(EXIT_FAILURE);
    }
    else{
        mapperList[mapperNameToMap] = new NearestNeighborMapper(AnumNodes, Anodes, BnumNodes, Bnodes);
    }
}

void init_FE_NearestElementMapper(char* mapperName,
                                  int AnumNodes, int AnumElems, const int *AnumNodesPerElem, const double *Anodes, const int *AnodeIDs, const int *Aelems,
                                  int BnumNodes, int BnumElems, const int *BnumNodesPerElem, const double *Bnodes, const int *BnodeIDs, const int *Belems){
    
    std::string mapperNameToMap = std::string(mapperName);
    
    if (mapperList.count( mapperNameToMap )){
        ERROR_OUT("A mapper with name : " + mapperNameToMap + "has already been initialized!");
        exit(EXIT_FAILURE);
    }
    else{
        mapperList[mapperNameToMap] = new NearestElementMapper(AnumNodes, AnumElems, AnumNodesPerElem, Anodes, AnodeIDs, Aelems,
                                                               BnumNodes, BnumElems, BnumNodesPerElem, Bnodes, BnodeIDs, Belems);
    }
}

void init_FE_BarycentricInterpolationMapper(char* mapperName, 
                                            int AnumNodes, const double *Anodes, int BnumNodes, const double *Bnodes){

    std::string mapperNameToMap = std::string(mapperName);
    
    if (mapperList.count( mapperNameToMap )){
        ERROR_OUT("A mapper with name : " + mapperNameToMap + "has already been initialized!");
        exit(EXIT_FAILURE);
    } else {
        mapperList[mapperNameToMap] = new BarycentricInterpolationMapper(AnumNodes, Anodes, BnumNodes, Bnodes);
    }
}

void init_FE_MortarMapper(char* mapperName,
                          int AnumNodes, int AnumElems, const int* AnumNodesPerElem, const double* Anodes, const int* AnodeIDs, const int* Aelems,
                          int BnumNodes, int BnumElems, const int* BnumNodesPerElem, const double* Bnodes, const int* BnodeIDs, const int* Belems,
                          int oppositeSurfaceNormal, int dual, int enforceConsistency){

    bool _oppositeSurfaceNormal = false;
    bool _dual = false;
    bool _enforceConsistency = false;

    if (oppositeSurfaceNormal != 0) _oppositeSurfaceNormal = true;
    if (dual != 0) _dual = true;
    if (enforceConsistency != 0) _enforceConsistency = true;

    std::string mapperNameToMap = std::string(mapperName);

    if (mapperList.count( mapperNameToMap )){
        ERROR_OUT("A mapper with name : " + mapperNameToMap + "has already been initialized!");
        exit(EXIT_FAILURE);
    } else {
        mapperList[mapperNameToMap] = new MortarMapper(AnumNodes, AnumElems, AnumNodesPerElem, Anodes, AnodeIDs, Aelems,
                                                       BnumNodes, BnumElems, BnumNodesPerElem, Bnodes, BnodeIDs, Belems,
                                                       _oppositeSurfaceNormal, _dual, _enforceConsistency);
    }

}

void doConsistentMapping(char* mapperName, int dimension, int dataSizeA, const double* dataA, int dataSizeB, double* dataB){
    assert(dimension == 1 || dimension == 3);
    
    std::string mapperNameToMap = std::string(mapperName);

    if (!mapperList.count( mapperNameToMap )){
        ERROR_OUT("This mapper does not exist : " + mapperNameToMap);
        exit(EXIT_FAILURE);
    } else {
        // if a vector field is to be mapped x, y, z components are extracted
        if (dimension == 3){
            int sizeDataToMap = dataSizeA/dimension;
            int sizeDataToWrite = dataSizeB/dimension;

            double** dataAtoMap = new double*[dimension];
            double** dataBtoWrite = new double*[dimension];

            for (int i=0; i<dimension ; i++){
                dataAtoMap[i]=new double[sizeDataToMap];
                dataBtoWrite[i]=new double[sizeDataToWrite];
            }
            for (int i=0 ; i<dimension ; i++){
                for (int j=0 ; j<sizeDataToWrite; j++){
                    dataAtoMap[i][j] = dataA[j*dimension+i];
                }

                mapperList[mapperNameToMap]->consistentMapping(dataAtoMap[i], dataBtoWrite[i]);
                for (int j=0 ; j<sizeDataToWrite; j++){
                    dataB[j*dimension+i] = dataBtoWrite[i][j];
                }
            }
            for (int i = 0; i<dimension; i++){
                delete[] dataAtoMap[i];
                delete[] dataBtoWrite[i];
            }
            delete[] dataAtoMap;
            delete[] dataBtoWrite;
        }
        // else field is mapped as it is
        else {
            mapperList[mapperNameToMap]->consistentMapping(dataA, dataB);
        }
    }
}

void doConservativeMapping(char* mapperName, int dimension, int dataSizeB, const double* dataB, int dataSizeA, double* dataA){
    assert(dimension == 1 || dimension == 3);

    std::string mapperNameToMap = std::string(mapperName);
    if (!mapperList.count( mapperNameToMap )){
        ERROR_OUT("This mapper does not exist : " + mapperNameToMap);
        exit(EXIT_FAILURE);
    } else {
        // if a vector field is to be mapped x, y, z components are extracted
        if (dimension == 3){
            int sizeDataToMap = dataSizeB/dimension;
            int sizeDataToWrite = dataSizeA/dimension;

            double** dataBtoMap = new double*[dimension];
            double** dataAtoWrite = new double*[dimension];

            for (int i=0; i<dimension ; i++){
                dataBtoMap[i]=new double[sizeDataToMap];
                dataAtoWrite[i]=new double[sizeDataToWrite];
            }
            for (int i=0 ; i<dimension ; i++){
                for (int j=0 ; j<sizeDataToMap; j++){
                    dataBtoMap[i][j] = dataB[j*dimension+i];
                }
                mapperList[mapperNameToMap]->conservativeMapping(dataBtoMap[i], dataAtoWrite[i]);
                for (int j=0 ; j<sizeDataToWrite; j++){
                    dataA[j*dimension+i] = dataAtoWrite[i][j];
                }
            }
            for (int i = 0; i<dimension; i++){
                delete[] dataBtoMap[i];
                delete[] dataAtoWrite[i];
            }
            delete[] dataBtoMap;
            delete[] dataAtoWrite;
        }
        // else field is mapped as it is
        else {
            mapperList[mapperNameToMap]->conservativeMapping(dataB, dataA);
        }
    }
}   

void deleteMapper(char* mapperName){
    std::string mapperNameToMap = std::string(mapperName);
    if (!mapperList.count( mapperNameToMap )){
        ERROR_OUT("This mapper does not exist : " + mapperNameToMap);
        exit(EXIT_FAILURE);
    } else {
        delete mapperList[mapperNameToMap];
    }
}

void deleteAllMappers(){
    
    std::map<std::string, AbstractMapper*>::iterator iter = mapperList.begin();
    
    if(iter!=mapperList.end()){
        delete iter->second;
        mapperList.erase(iter);
    }
}
