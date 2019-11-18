#include "matrix.h"
#include "mex.h"
#include <assert.h>
#include "EMPIRE_API.h"
#include <iostream>

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    assert(nrhs==2);
    assert(nlhs==0);

#define INNER_IN           prhs[0]
#define NUMCURVES_IN       prhs[1]

    // Flag on trimming
    assert(mxIsDouble(INNER_IN));
    assert(mxGetNumberOfElements(INNER_IN) == 1);
    int isInner = (int) mxGetPr(INNER_IN)[0]; // cast from double to int

    // Number of trimming loops
    assert(mxIsDouble(NUMCURVES_IN));
    assert(mxGetNumberOfElements(NUMCURVES_IN) == 1);
    int numCurves = (int) mxGetPr(NUMCURVES_IN)[0]; // cast from double to int

    // Send the IGA trimming 
	EMPIRE_API_sendIGATrimmingLoopInfo(isInner, numCurves);

#undef INNER_IN
#undef NUMCURVES_IN
}
