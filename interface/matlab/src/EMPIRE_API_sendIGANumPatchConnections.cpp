#include "matrix.h"
#include "mex.h"
#include <assert.h>
#include "EMPIRE_API.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    assert(nrhs==1);
    assert(nlhs==0);

#define NUMCONNECTIONS_IN       prhs[0]

    // size of the datafield array
    assert(mxIsDouble(NUMCONNECTIONS_IN));
    assert(mxGetNumberOfElements(NUMCONNECTIONS_IN) == 1);
    int numConnections = (int) mxGetPr(NUMCONNECTIONS_IN)[0]; // cast from double to int

	EMPIRE_API_sendIGANumPatchConnections(numConnections);

#undef NUMCONNECTIONS_IN
}
