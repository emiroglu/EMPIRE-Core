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

#define COUPLINGCURVENUMGP_IN              prhs[0]
#define TRCURVEMASTERGPS_IN                prhs[1]
#define TRCURVESLAVEGPS_IN                 prhs[2]
#define TRCURVEGPWEIGHTS_IN                prhs[3]
#define TRCURVEMASTERGPTANGENTS_IN         prhs[4]
#define TRCURVESLAVEGPTANGENTS_IN          prhs[5]
#define TRCURVEGPJACOBIANPRODUCTS_IN       prhs[6]
	
	// Number of Gauss Points on the trimming curve
	assert(mxIsDouble(COUPLINGCURVENUMGP_IN));
	assert(mxGetNumberOfElements(COUPLINGCURVENUMGP_IN) == 1);
	int couplingCurveNumGP = (int) mxGetPr(COUPLINGCURVENUMGP_IN)[0]; // cast from double to int
	
	// Parametric location of the Gauss Points on the master trimming curve
	assert(mxIsDouble(TRCURVEMASTERGPS_IN));
	assert(mxGetNumberOfElements(TRCURVEMASTERGPS_IN) == 2*couplingCurveNumGP);
	double *trCurveMasterGPs = mxGetPr(TRCURVEMASTERGPS_IN);
	
	// Parametric location of the Gauss Points on the slave trimming curve
	assert(mxIsDouble(TRCURVESLAVEGPS_IN));
	assert(mxGetNumberOfElements(TRCURVESLAVEGPS_IN) == 2*couplingCurveNumGP);
	double *trCurveSlaveGPs = mxGetPr(TRCURVESLAVEGPS_IN);
	
	// Gauss Weights
	assert(mxIsDouble(TRCURVEGPWEIGHTS_IN));
	assert(mxGetNumberOfElements(TRCURVEGPWEIGHTS_IN) == couplingCurveNumGP);
	double *trCurveGPWeights = mxGetPr(TRCURVEGPWEIGHTS_IN); // cast from double to int
	
	// Tangent vectors along the master trimming curve in the physical space
	assert(mxIsDouble(TRCURVEMASTERGPTANGENTS_IN));
	assert(mxGetNumberOfElements(TRCURVEMASTERGPTANGENTS_IN) == 3*couplingCurveNumGP);
	double *trCurveMasterGPTangents = mxGetPr(TRCURVEMASTERGPTANGENTS_IN);
	
	// Tangent vectors along the slave trimming curve in the physical space
	assert(mxIsDouble(TRCURVESLAVEGPTANGENTS_IN));
	assert(mxGetNumberOfElements(TRCURVESLAVEGPTANGENTS_IN) == 3*couplingCurveNumGP);
	double *trCurveSlaveGPTangents = mxGetPr(TRCURVESLAVEGPTANGENTS_IN);
	
	// Jacobian products (including the Gauss Weights)
	assert(mxIsDouble(TRCURVEGPJACOBIANPRODUCTS_IN));
	assert(mxGetNumberOfElements(TRCURVEGPJACOBIANPRODUCTS_IN) == couplingCurveNumGP);
	double *trCurveGPJacobianProducts = mxGetPr(TRCURVEGPJACOBIANPRODUCTS_IN);
	
	// Send the data for the application of weak Dirichlet boundary conditions to Empire through API
    EMPIRE_API_sendIGAPatchConnectionData(couplingCurveNumGP, trCurveMasterGPs, trCurveSlaveGPs,trCurveGPWeights, trCurveMasterGPTangents, trCurveSlaveGPTangents, trCurveGPJacobianProducts);

#undef COUPLINGCURVENUMGP_IN
#undef TRCURVEMASTERGPS_IN
#undef TRCURVESLAVEGPS_IN
#undef TRCURVEGPWEIGHTS_IN
#undef TRCURVEMASTERGPTANGENTS_IN
#undef TRCURVESLAVEGPTANGENTS_IN
#undef TRCURVEGPJACOBIANPRODUCTS_IN
}