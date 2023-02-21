#include "matrix.h"
#include "mex.h"
#include <assert.h>
#include "EMPIRE_API.h"
#include "HelperFunctions.h"

#include <iostream>
using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    assert(nrhs==6);
    assert(nlhs==0);

#define DIRECTION_IN                prhs[0]	
#define P_DEGREE_IN                 prhs[1]
#define U_NUM_KNOTS_IN              prhs[2]
#define U_KNOT_VECTOR_IN            prhs[3]
#define U_NUM_CONTROL_PTS_IN        prhs[4]
#define CP_NET_IN                   prhs[5]
	
	// direction of the trimming curve
	assert(mxIsDouble(DIRECTION_IN));
    assert(mxGetNumberOfElements(DIRECTION_IN) == 1);
    int direction = (int) mxGetPr(DIRECTION_IN)[0]; // cast from double to int

    // p degree
    assert(mxIsDouble(P_DEGREE_IN));
    assert(mxGetNumberOfElements(P_DEGREE_IN) == 1);
    int pDegree = (int) mxGetPr(P_DEGREE_IN)[0]; // cast from double to int

    // u number of knots
    assert(mxIsDouble(U_NUM_KNOTS_IN));
    assert(mxGetNumberOfElements(U_NUM_KNOTS_IN) == 1);
    int uNumKnots = (int) mxGetPr(U_NUM_KNOTS_IN)[0]; // cast from double to int

    // u knot vector
    assert(mxIsDouble(U_KNOT_VECTOR_IN));
    assert(mxGetNumberOfElements(U_KNOT_VECTOR_IN) == uNumKnots);
    double *uKnotVector = mxGetPr(U_KNOT_VECTOR_IN);

    // u number of control points
    assert(mxIsDouble(U_NUM_CONTROL_PTS_IN));
    assert(mxGetNumberOfElements(U_NUM_CONTROL_PTS_IN) == 1);
    int uNumControlPoints = (int) mxGetPr(U_NUM_CONTROL_PTS_IN)[0]; // cast from double to int

    // coordinates of control points
    assert(mxIsDouble(CP_NET_IN));
    assert(mxGetNumberOfElements(CP_NET_IN) == uNumControlPoints * 4);
    double *cpNet = mxGetPr(CP_NET_IN);

	// Send the trimming curve information through API to Empire
    EMPIRE_API_sendIGATrimmingCurve(direction, pDegree,  uNumKnots, uKnotVector, uNumControlPoints, cpNet);

#undef DIRECTION_IN
#undef P_DEGREE_IN
#undef U_NUM_KNOTS_IN
#undef U_KNOT_VECTOR_IN
#undef U_NUM_CONTROL_PTS_IN
#undef CP_NET_IN
}