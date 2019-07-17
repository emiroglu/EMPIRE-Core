#include "matrix.h"
#include "mex.h"
#include <assert.h>
#include "EMPIRE_API.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    assert(nrhs==1);
    assert(nlhs==0);

#define NUMDIRICHLETCONDITIONS_IN       prhs[0]

    // size of the datafield array
    assert(mxIsDouble(NUMDIRICHLETCONDITIONS_IN));
    assert(mxGetNumberOfElements(NUMDIRICHLETCONDITIONS_IN) == 1);
    int numDirichletConditions = (int) mxGetPr(NUMDIRICHLETCONDITIONS_IN)[0]; // cast from double to int

	EMPIRE_API_sendIGANumDirichletConditions(numDirichletConditions);

#undef NUMDIRICHLETCONDITIONS_IN
}
