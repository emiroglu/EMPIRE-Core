#include "matrix.h"
#include "mex.h"
#include <assert.h>
#include "EMPIRE_API.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    assert(nrhs==2);
    assert(nlhs==0);

#define ISTRIMMED_IN       prhs[0]
#define NUMLOOPS_IN       prhs[1]

    // Flag on trimming
    assert(mxIsDouble(ISTRIMMED_IN));
    assert(mxGetNumberOfElements(ISTRIMMED_IN) == 1);
    int isTrimmed = (int) mxGetPr(ISTRIMMED_IN)[0]; // cast from double to int

    // Number of trimming loops
    assert(mxIsDouble(NUMLOOPS_IN));
    assert(mxGetNumberOfElements(NUMLOOPS_IN) == 1);
    int numLoops = (int) mxGetPr(NUMLOOPS_IN)[0]; // cast from double to int

	EMPIRE_API_sendIGATrimmingInfo(isTrimmed, numLoops);

#undef ISTRIMMED_IN
#undef NUMLOOPS_IN
}
