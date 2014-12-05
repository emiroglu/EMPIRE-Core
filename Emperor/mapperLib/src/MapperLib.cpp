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

//TODO : add comments

//TODO : before initialization of a mapper check if another one with the same name exists

void init_FE_NearestNeighborMapper(char* mapperName, 
				int AnumNodes, const double *Anodes, int BnumNodes, const double *Bnodes){
    
    std::string mapperNameToMap = std::string(mapperName); 
    
    mapperList[mapperNameToMap] = new NearestNeighborMapper(AnumNodes, Anodes, BnumNodes, Bnodes);
}

void init_FE_NearestElementMapper(char* mapperName,
			       int AnumNodes, int AnumElems, const int *AnumNodesPerElem, const double *Anodes, const int *AnodeIDs, const int *Aelems, 
			       int BnumNodes, int BnumElems, const int *BnumNodesPerElem, const double *Bnodes, const int *BnodeIDs, const int *Belems){
    
    std::string mapperNameToMap = std::string(mapperName); 
    
    mapperList[mapperNameToMap] = new NearestElementMapper(AnumNodes, AnumElems, AnumNodesPerElem, Anodes, AnodeIDs, Aelems, 
							   BnumNodes, BnumElems, BnumNodesPerElem, Bnodes, BnodeIDs, Belems);
}

void init_FE_BarycentricInterpolationMapper(char* mapperName, 
					 int AnumNodes, const double *Anodes, int BnumNodes, const double *Bnodes){

    std::string mapperNameToMap = std::string(mapperName);
    
    mapperList[mapperNameToMap] = new BarycentricInterpolationMapper(AnumNodes, Anodes, BnumNodes, Bnodes);
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
    
    mapperList[mapperNameToMap] = new MortarMapper(AnumNodes, AnumElems, AnumNodesPerElem, Anodes, AnodeIDs, Aelems,
						   BnumNodes, BnumElems, BnumNodesPerElem, Bnodes, BnodeIDs, Belems,
						   _oppositeSurfaceNormal, _dual, _enforceConsistency);

}

//TODO : before doing mapping check if the mapper with the called name exists
void doConsistentMapping(char* mapperName, int dimension, int dataSize, const double* dataA, double* dataB){
    assert(dimension == 1 || dimension == 3);
    
    std::string mapperNameToMap = std::string(mapperName);
    
    // if a vector field is to be mapped x, y, z components are extracted
    if (dimension == 3){
        int sizeDataToMap = dataSize/dimension;

        double** dataAtoMap = new double*[dimension];
        double** dataBtoWrite = new double*[dimension];

        for (int i=0; i<dimension ; i++){
            dataAtoMap[i]=new double[sizeDataToMap];
            dataBtoWrite[i]=new double[sizeDataToMap];
        }
        for (int i=0 ; i<dimension ; i++){
            for (int j=0 ; j<sizeDataToMap; j++){
                dataAtoMap[i][j] = dataA[j*dimension+i];
            }
            
            mapperList[mapperNameToMap]->consistentMapping(dataAtoMap[i], dataBtoWrite[i]);
            for (int j=0 ; j<sizeDataToMap; j++){
                dataB[j*dimension+i] = dataBtoWrite[i][j];
            }
        }
    }
    // else field is mapped as it is
    else {
        mapperList[mapperNameToMap]->consistentMapping(dataA, dataB);
    }
}

void doConservativeMapping(char* mapperName, int dimension, int dataSize, const double* dataB, double* dataA){
    assert(dimension == 1 || dimension == 3);

    std::string mapperNameToMap = std::string(mapperName);
    
    // if a vector field is to be mapped x, y, z components are extracted
    if (dimension == 3){
        int sizeDataToMap = dataSize/dimension;

        double** dataBtoMap = new double*[dimension];
        double** dataAtoWrite = new double*[dimension];

        for (int i=0; i<dimension ; i++){
            dataBtoMap[i]=new double[sizeDataToMap];
            dataAtoWrite[i]=new double[sizeDataToMap];
        }
        for (int i=0 ; i<dimension ; i++){
            for (int j=0 ; j<sizeDataToMap; j++){
                dataBtoMap[i][j] = dataB[j*dimension+i];
            }
            mapperList[mapperNameToMap]->conservativeMapping(dataBtoMap[i], dataAtoWrite[i]);
            for (int j=0 ; j<sizeDataToMap; j++){
                dataA[j*dimension+i] = dataAtoWrite[i][j];
            }
        }
    }
    // else field is mapped as it is
    else {
        mapperList[mapperNameToMap]->conservativeMapping(dataB, dataA);
    }
}   

//TODO : before deletion of a mapper check if a mapper with given name exists
void deleteMapper(char* mapperName){
     std::string mapperNameToMap = std::string(mapperName);
     delete mapperList[mapperNameToMap];
}

void deleteAllMappers(){
    
    std::map<std::string, AbstractMapper*>::iterator iter = mapperList.begin();
    
    if(iter!=mapperList.end()){
      delete iter->second;
      mapperList.erase(iter);
    }
    
}
