#include "matrix.h"
#include "mex.h"
#include <assert.h>
#include "EMPIRE_API.h"
#include "HelperFunctions.h"
#include <iostream>

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    assert(nrhs==7);
    assert(nlhs==0);

#define MASTERPATCHCTR_IN          prhs[0]
#define MASTERPATCHBLCTR_IN        prhs[1]
#define MASTERPATCHBLTCCTR_IN      prhs[2]
#define SLAVEPATCHCTR_IN           prhs[3]
#define SLAVEPATCHBLCTR_IN         prhs[4]
#define SLAVEPATCHBLTCCTR_IN       prhs[5]
#define ISCOUPLINGGPPROVIDED_IN    prhs[6]
	
	// Counter of the master patch
	assert(mxIsDouble(MASTERPATCHCTR_IN));
    assert(mxGetNumberOfElements(MASTERPATCHCTR_IN) == 1);
    int masterPatchCtr = (int) mxGetPr(MASTERPATCHCTR_IN)[0]; // cast from double to int

    // Counter of the boundary loop of the master patch
    assert(mxIsDouble(MASTERPATCHBLCTR_IN));
    assert(mxGetNumberOfElements(MASTERPATCHBLCTR_IN) == 1);
    int masterPatchBLCtr = (int) mxGetPr(MASTERPATCHBLCTR_IN)[0]; // cast from double to int

    // Counter of the trimming curve of the boundary loop of the master patch
    assert(mxIsDouble(MASTERPATCHBLTCCTR_IN));
    assert(mxGetNumberOfElements(MASTERPATCHBLTCCTR_IN) == 1);
    int masterPatchBLTrCurveCtr = (int) mxGetPr(MASTERPATCHBLTCCTR_IN)[0]; // cast from double to int

    // Counter of the slave patch
    assert(mxIsDouble(SLAVEPATCHCTR_IN));
    assert(mxGetNumberOfElements(SLAVEPATCHCTR_IN) == 1);
    int slavePatchCtr = (int) mxGetPr(SLAVEPATCHCTR_IN)[0];
	
	// Counter of the boundary loop of the slave patch
    assert(mxIsDouble(SLAVEPATCHBLCTR_IN));
    assert(mxGetNumberOfElements(SLAVEPATCHBLCTR_IN) == 1);
    int slavePatchBLCtr = (int) mxGetPr(SLAVEPATCHBLCTR_IN)[0];
	
	// Counter of the trimming curve of the boundary loop of the slave patch
    assert(mxIsDouble(SLAVEPATCHBLTCCTR_IN));
    assert(mxGetNumberOfElements(SLAVEPATCHBLTCCTR_IN) == 1);
    int slavePatchBLTrCurveCtr = (int) mxGetPr(SLAVEPATCHBLTCCTR_IN)[0];
	
	// Flag on whether the Gauss Point data for the application of the interface continuity conditions are provided
    assert(mxIsDouble(ISCOUPLINGGPPROVIDED_IN));
    assert(mxGetNumberOfElements(ISCOUPLINGGPPROVIDED_IN) == 1);
    int isCouplingGPProvided = (int) mxGetPr(ISCOUPLINGGPPROVIDED_IN)[0];

	// Send the information for the application of weak Dirichlet boundary conditions to Empire through API
    EMPIRE_API_sendIGAPatchConnectionInfo(masterPatchCtr, masterPatchBLCtr,  masterPatchBLTrCurveCtr, slavePatchCtr, slavePatchBLCtr, slavePatchBLTrCurveCtr, isCouplingGPProvided);

#undef MASTERPATCHCTR_IN
#undef MASTERPATCHBLCTR_IN
#undef MASTERPATCHBLTCCTR_IN
#undef SLAVEPATCHCTR_IN
#undef SLAVEPATCHBLCTR_IN
#undef SLAVEPATCHBLTCCTR_IN
#undef ISCOUPLINGGPPROVIDED_IN
}