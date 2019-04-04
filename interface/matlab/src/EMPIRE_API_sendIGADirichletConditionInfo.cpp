#include "matrix.h"
#include "mex.h"
#include <assert.h>
#include "EMPIRE_API.h"
#include "HelperFunctions.h"
#include <iostream>

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    assert(nrhs==4);
    assert(nlhs==0);

#define PATCHCTR_IN                prhs[0]	
#define PATCHBLCTR_IN              prhs[1]
#define PATCHBLTCCTR_IN            prhs[2]
#define ISGPDATAPROVIDED_IN        prhs[3]
	
	// Counter of the patch
	assert(mxIsDouble(PATCHCTR_IN));
    assert(mxGetNumberOfElements(PATCHCTR_IN) == 1);
    int patchCtr = (int) mxGetPr(PATCHCTR_IN)[0]; // cast from double to int

    // Counter of the boundary loop of the patch
    assert(mxIsDouble(PATCHBLCTR_IN));
    assert(mxGetNumberOfElements(PATCHBLCTR_IN) == 1);
    int patchBLCtr = (int) mxGetPr(PATCHBLCTR_IN)[0]; // cast from double to int

    // Counter of the trimming curve of the boundary loop of the patch
    assert(mxIsDouble(PATCHBLTCCTR_IN));
    assert(mxGetNumberOfElements(PATCHBLTCCTR_IN) == 1);
    int patchBLTrCurveCtr = (int) mxGetPr(PATCHBLTCCTR_IN)[0]; // cast from double to int

    // Flag on whether the Gauss Point data for the application of weak Dirichlet boundary conditions are provided
    assert(mxIsDouble(ISGPDATAPROVIDED_IN));
    assert(mxGetNumberOfElements(ISGPDATAPROVIDED_IN) == 1);
    int isDirichletGPProvided = (int) mxGetPr(ISGPDATAPROVIDED_IN)[0];

	// Send the information for the application of weak Dirichlet boundary conditions to Empire through API
    EMPIRE_API_sendIGADirichletConditionInfo(patchCtr, patchBLCtr,  patchBLTrCurveCtr, isDirichletGPProvided);

#undef PATCHCTR_IN
#undef PATCHBLCTR_IN
#undef PATCHBLTCCTR_IN
#undef ISGPDATAPROVIDED_IN
}