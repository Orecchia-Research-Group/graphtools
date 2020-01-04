/*
C MATLAB function: cutweird
 
USAGE: [(double) wr_num, (int64) wr_den,  (double) wr] = cutexp(sparse matrix (double) G, vector (int64) cut, vector (int64) bisec);

PURPOSE: compute the weird ratio of the cut with respect to the near bisection bisec.

NOTES: 
   - if weird ratio is negative, return abs.value (i.e. weirdratio of complement cut)
   - G is assumed to be undirected, no check for that
   - G is assumed to be sparse, program will check this
   - ASSUMING C LONG TYPE IS 64 BITS


mexFunction INPUTS;
   nrhs = 3
   nlhs = 3
   
*/
#include "mex.h"
#include "matrix.h"

#define abs(x) (x) > 0 ? (x) : -(x)

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    const mxArray *G;
    mwIndex *col;
    mwIndex *row;

    double *array_G;
    long *cut;
    long *reciprocal_cut;
    long size_cut;
    long reciprocal_size_cut;
    long *bisec;
    long size_bisec;
    long i, j;
    long n;
    double cutedges = 0;
    long size_intersect;
    long reciprocal_size_intersect;
    long denominator;
    long lamda;

    int *mask_cut;
    int *reciprocal_mask_cut;
    mxArray *exp;

    mwSize dims[] = {1, 1};
    long *temp;

    /*%%%%%%%%%%%%%%%%% ARGUMENT LOADING &  CHECKING %%%%%%%%%%%%%%%%%%%%%*/

    /* CHECK CORRECT NUMBER OF INPUT/OUTPUTS */
    if (nrhs != 5 || nlhs != 3)
        mexErrMsgTxt("Error in cutweird. Incorrect usage.\n");

    /* CHECK TYPES */
    if (!mxIsSparse(prhs[0]))
        mexErrMsgTxt("Error in cutweird. Graph must be sparse.\n");

    if (mxGetClassID(prhs[1]) != mxINT64_CLASS)
        mexErrMsgTxt("Error in cutweird. Cut must be of class int64.\n");

    if (mxGetClassID(prhs[2]) != mxINT64_CLASS)
        mexErrMsgTxt("Error in cutweird. Cut must be of class int64.\n");

    if (mxGetClassID(prhs[3]) != mxINT64_CLASS)
        mexErrMsgTxt("Error in cutweird. Bisec must be of class int64.\n");


    /* LOAD ARGUMENTS */
    G = prhs[0];
    n = mxGetM(G);
    col = mxGetJc(G);
    row = mxGetIr(G);
    array_G = mxGetPr(G);

    cut = (long *) mxGetPr(prhs[1]);
    size_cut = mxGetM(prhs[1]);
    reciprocal_cut = (long *) mxGetPr(prhs[2]);
    reciprocal_size_cut = mxGetM(prhs[2]);

    bisec = (long *) mxGetPr(prhs[3]);
    size_bisec = mxGetM(prhs[3]);

    lamda = mxGetScalar(prhs[4]);


    /*%%%%%%%%%%%%%%%%%%%%%%%%% MAIN BODY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

    /* PREPARE CUT MASK  AND COMPUTE INTERSECT*/
    mask_cut = calloc(sizeof(*mask_cut), n + 1);
    if (mask_cut == NULL)
        mexErrMsgTxt("Error allocating memory in cutweird.");

    for (i = 0; i < size_cut; i++)
        mask_cut[cut[i]] = 1;

    reciprocal_mask_cut = calloc(sizeof(*reciprocal_mask_cut), n + 1);
    if (reciprocal_mask_cut == NULL)
        mexErrMsgTxt("Error allocating memory in cutweird.");

    for (i = 0; i < reciprocal_size_cut; i++)
        reciprocal_mask_cut[reciprocal_cut[i]] = 1;

    size_intersect = 0;
    for (i = 0; i < size_bisec; i++)
        if (mask_cut[bisec[i]] == 1)
            size_intersect++;

    reciprocal_size_intersect = 0;
    for (i = 0; i < size_bisec; i++)
        if (reciprocal_mask_cut[bisec[i]] == 1)
            reciprocal_size_intersect++;


    /* COMPUTE WEIRD RATIO*/
    /* COMPUTE EDGES CUT */
    for (i = 0; i < size_cut; i++) {
        if (reciprocal_mask_cut[cut[i]] == 1)
            continue;
        for (j = col[cut[i] - 1]; j < col[cut[i]]; j++)
            if (mask_cut[row[j] + 1] == 0)
                cutedges += array_G[j];
    }

    printf("cutedges: %lf\n", cutedges);
    if (lamda > 0) cutedges += (size_cut + reciprocal_size_cut - n);
    /* COMPUTE DENOMINATOR */
    if (abs(2 * size_intersect - size_cut) < abs(2 * reciprocal_size_intersect - reciprocal_size_cut))
        denominator = 2 * size_intersect - size_cut;
    else
        denominator = 2 * reciprocal_size_intersect - reciprocal_size_cut;

    if (denominator < 0)
        denominator = denominator * (-1);

    /*%%%%%%%%%%%%%%%%%%%%% TERMINATION AND CLEANING %%%%%%%%%%%%%%%%%%%%%%%%*/

    plhs[0] = mxCreateDoubleScalar(cutedges);
    plhs[1] = mxCreateNumericArray(2, dims, mxINT64_CLASS, mxREAL);
    plhs[2] = mxCreateDoubleScalar(cutedges / (double) denominator);
    temp = (long *) mxGetPr(plhs[1]);
    *temp = denominator;

    free(mask_cut);

} 



  


  
