#ifndef MAPPERLIB_H_
#define MAPPERLIB_H_

#ifdef __cplusplus
extern "C" { ///Define extern C if C++ compiler is used
#endif
//TODO : add comments
void init_FE_MortarMapper(char* mapperName, int AnumNodes, int AnumElems, const int* AnumNodesPerElem, const double* Anodes, const int* AnodeIDs, const int* Aelems, int BnumNodes, int BnumElems, const int* BnumNodesPerElem, const double* Bnodes, const int* BnodeIDs, const int* Belems, int oppositeSurfaceNormal, int dual, int enforceConsistency);

void init_FE_NearestNeighborMapper(char* mapperName, 
				int AnumNodes, const double *Anodes, int BnumNodes, const double *Bnodes);

void init_FE_NearestElementMapper(char* mapperName,
			       int AnumNodes, int AnumElems, const int *AnumNodesPerElem, const double *Anodes, const int *AnodeIDs, const int *Aelems, 
			       int BnumNodes, int BnumElems, const int *BnumNodesPerElem, const double *Bnodes, const int *BnodeIDs, const int *Belems);

void init_FE_BarycentricInterpolationMapper(char* mapperName, 
					 int AnumNodes, const double *Anodes, int BnumNodes, const double *Bnodes);

void doConsistentMapping(char* mapperName, int dimension, int dataSize, const double* dataA, double* dataB);

void doConservativeMapping(char* mapperName, int dimension, int dataSize, const double* dataB, double* dataA);

void deleteMapper(char* mapperName);

void deleteAllMappers();

#ifdef __cplusplus
}
#endif

#endif // MAPPERLIB_H_
