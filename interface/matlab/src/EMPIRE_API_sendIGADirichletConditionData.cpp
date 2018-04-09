#include "matrix.h"
#include "mex.h"
#include <assert.h>
#include "EMPIRE_API.h"
#include "HelperFunctions.h"
#include <iostream>

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    assert(nrhs==5);
    assert(nlhs==0);

#define TRCURVENUMGP_IN                    prhs[0]
#define TRCURVEGPS_IN                      prhs[1]
#define TRCURVEGPWEIGHTS_IN                prhs[2]
#define TRCURVEGPTANGENTS_IN               prhs[3]
#define TRCURVEGPJACOBIANPRODUCTS_IN       prhs[4]
	
	// Number of Gauss Points on the trimming curve
	assert(mxIsDouble(TRCURVENUMGP_IN));
	assert(mxGetNumberOfElements(TRCURVENUMGP_IN) == 1);
	int trCurveNumGP = (int) mxGetPr(TRCURVENUMGP_IN)[0]; // cast from double to int
	
	// Parametric location of the Gauss Points on the trimming curve
	assert(mxIsDouble(TRCURVEGPS_IN));
	assert(mxGetNumberOfElements(TRCURVEGPS_IN) == 2*trCurveNumGP);
	double *trCurveGPs = mxGetPr(TRCURVEGPS_IN);
	
	// Gauss Weights
	assert(mxIsDouble(TRCURVEGPWEIGHTS_IN));
	assert(mxGetNumberOfElements(TRCURVEGPWEIGHTS_IN) == trCurveNumGP);
	double *trCurveGPWeights = mxGetPr(TRCURVEGPWEIGHTS_IN); // cast from double to int
	
	// Tangent vectors along the trimming curve in the physical space
	assert(mxIsDouble(TRCURVEGPTANGENTS_IN));
	assert(mxGetNumberOfElements(TRCURVEGPTANGENTS_IN) == 3*trCurveNumGP);
	double *trCurveGPTangents = mxGetPr(TRCURVEGPTANGENTS_IN);
	
	// Jacobian products (including the Gauss Weights)
	assert(mxIsDouble(TRCURVEGPJACOBIANPRODUCTS_IN));
	assert(mxGetNumberOfElements(TRCURVEGPJACOBIANPRODUCTS_IN) == trCurveNumGP);
	double *trCurveGPJacobianProducts = mxGetPr(TRCURVEGPJACOBIANPRODUCTS_IN);
	
	// Send the data for the application of weak Dirichlet boundary conditions to Empire through API
    EMPIRE_API_sendIGADirichletConditionData(trCurveNumGP, trCurveGPs,  trCurveGPWeights, trCurveGPTangents, trCurveGPJacobianProducts);

#undef TRCURVENUMGP_IN
#undef TRCURVEGPS_IN
#undef TRCURVEGPWEIGHTS_IN
#undef TRCURVEGPTANGENTS_IN
#undef TRCURVEGPJACOBIANPRODUCTS_IN
}