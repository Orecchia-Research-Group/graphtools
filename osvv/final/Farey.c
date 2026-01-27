/*
C MATLAB function: Farey
 
USAGE: [num,den] = Farey(fractnum, fracden, p)

PURPOSE: compute the best fractional approximation larger than fracnum/fracden with
numerator <= p.
 
NOTES:  p < 100000

*/

#include <math.h>
#include <stdint.h>
#include "mex.h"
#include "matrix.h"
#include "farey.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    int64_t num_in;
    int64_t den_in;
    int64_t p;
    int64_t num_out;
    int64_t den_out;

    mwSize dims[] = {1, 1};
    int64_t *temp;

    /*%%%%%%%%%%%%%%%%% ARGUMENT LOADING &  CHECKING %%%%%%%%%%%%%%%%%%%%%*/
    /* CHECK CORRECT NUMBER OF INPUT/OUTPUTS */
    if (nrhs != 3 || nlhs != 2)
        mexErrMsgTxt("Error in Farey. Incorrect usage.\n");

    /* CHECK TYPES */
    if (mxGetClassID(prhs[0]) != mxINT64_CLASS)
        mexErrMsgTxt("Error in Farey. Numerator must be of class int64.\n");

    if (mxGetClassID(prhs[1]) != mxINT64_CLASS)
        mexErrMsgTxt("Error in Farey. Denominator  must be of class int64.\n");

    if (mxGetClassID(prhs[2]) != mxINT64_CLASS)
        mexErrMsgTxt("Error in Farey. Precision must be of class int64.\n");


    num_in = ((int64_t *) mxGetPr(prhs[0]))[0];
    den_in = ((int64_t *) mxGetPr(prhs[1]))[0];
    p = ((int64_t *) mxGetPr(prhs[2]))[0];

    if (num_in < 0 || den_in < 0)
        mexErrMsgTxt("Error in Farey. Numerators and denominators must be positive.\n");

    if (p > 100000 || p < 0)
        mexErrMsgTxt("Error in Farey. Precision out of range.\n");

    /*%%%%%%%%%%%%%%%%%%%%%%%%% MAIN BODY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
    farey(num_in, den_in, p, &num_out, &den_out);

    /*%%%%%%%%%%%%%%%%%%%%% TERMINATION AND CLEANING %%%%%%%%%%%%%%%%%%%%%%%%*/

    plhs[0] = mxCreateNumericArray(2, dims, mxINT64_CLASS, mxREAL);
    plhs[1] = mxCreateNumericArray(2, dims, mxINT64_CLASS, mxREAL);

    temp = (int64_t *) mxGetPr(plhs[0]);
    *temp = num_out;

    temp = (int64_t *) mxGetPr(plhs[1]);
    *temp = den_out;

}
