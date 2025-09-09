/*
C MATLAB function: cutexp
 
USAGE: [(double) ex_num, (int64) ex_den,  (double) ex] = cutexp(sparse matrix (double) G, vector (int64) cut);

PURPOSE: compute the expansion of cut in G

NOTES: 
   - cut will be complemented if larger than half the vertices;
   - G is assumed to be undirected, no check for that
   - G is assumed to be sparse, program will check this
   - ASSUMING C int64_t TYPE IS 64 BITS


mexFunction INPUTS;
   nrhs = 2
   nlhs = 3
   
*/
#include <math.h>
#include <stdint.h>
#include "mex.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    const mxArray *G;
    mwIndex *col;
    mwIndex *row;

    double *array_G;
    int64_t *cut;
    int64_t *reciprocal_cut;
    int64_t size_cut;
    int64_t reciprocal_size_cut;
    int64_t i, j;
    int64_t n;
    int64_t lamda_num;
    int64_t lamda_den;
    int64_t cut_vol = 0;
    int64_t reciprocal_cut_vol = 0;
    double cutedges = 0;

    int *mask_cut;
    int *reciprocal_mask_cut;
    int64_t *weight;

    mwSize dims[] = {1, 1};
    int64_t *temp;

    /*%%%%%%%%%%%%%%%%% ARGUMENT LOADING &  CHECKING %%%%%%%%%%%%%%%%%%%%%*/

    /* CHECK CORRECT NUMBER OF INPUT/OUTPUTS */
    if (nrhs != 6 || nlhs != 3)
        mexErrMsgTxt("Error in cutexp. Incorrect usage.\n");

    /* CHECK TYPES */
    if (!mxIsSparse(prhs[0]))
        mexErrMsgTxt("Error in cutexp. Graph must be sparse.\n");

    if (mxGetClassID(prhs[1]) != mxINT64_CLASS)
        mexErrMsgTxt("Error in cutexp. lamda_num must be of class int64.\n");

    if (mxGetClassID(prhs[2]) != mxINT64_CLASS)
        mexErrMsgTxt("Error in cutexp. lamda_den must be of class int64.\n");

    if (mxGetClassID(prhs[3]) != mxINT64_CLASS)
        mexErrMsgTxt("Error in cutexp. Weight must be of class int64.\n");

    if (mxGetClassID(prhs[4]) != mxINT64_CLASS)
        mexErrMsgTxt("Error in cutexp. Cut must be of class int64.\n");

    if (mxGetClassID(prhs[5]) != mxINT64_CLASS)
        mexErrMsgTxt("Error in cutexp. Reciprocal cut must be of class int64.\n");

    /* LOAD ARGUMENTS */
    G = prhs[0];
    n = mxGetM(G);
    col = mxGetJc(G);
    row = mxGetIr(G);
    array_G = mxGetPr(G);
    lamda_num = ((int64_t *) mxGetPr(prhs[1]))[0];
    lamda_den = ((int64_t *) mxGetPr(prhs[2]))[0];
    if (lamda_num < 0) {
        lamda_den = 1l;
    }

    weight = (int64_t *) mxGetPr(prhs[3]);

    cut = (int64_t *) mxGetPr(prhs[4]);
    reciprocal_cut = (int64_t *) mxGetPr(prhs[5]);
    size_cut = mxGetM(prhs[4]);
    reciprocal_size_cut = mxGetM(prhs[5]);

    /* PREPARE MASK - COMPLEMENT CUT IF NECESSARY */
    mask_cut = calloc(sizeof(*mask_cut), n + 1);
    if (mask_cut == NULL)
        mexErrMsgTxt("Error allocating memory in cutexp.");

    reciprocal_mask_cut = calloc(sizeof(*reciprocal_mask_cut), n + 1);
    if (reciprocal_mask_cut == NULL)
        mexErrMsgTxt("Error allocating memory in cutexp.");


    for (i = 0; i < size_cut; i++)
        mask_cut[cut[i] - 1] = 1;

    for (i = 0; i < reciprocal_size_cut; i++)
        reciprocal_mask_cut[reciprocal_cut[i] - 1] = 1;

    for (i = 0; i < n; i++) {
        if (mask_cut[i]) {
            cut_vol += weight[i] * lamda_den;
            if (!reciprocal_mask_cut[i]) {
                for (j = col[i]; j < col[i + 1]; j++)
                    if (!mask_cut[row[j]])
                        cutedges += array_G[j] * lamda_den;
            }
            else if (lamda_num > 0){
                cutedges += weight[i] * lamda_num;
            }
        }
        if (reciprocal_mask_cut[i])
            reciprocal_cut_vol += weight[i] * lamda_den;
    }

    int64_t denominator = cut_vol < reciprocal_cut_vol ? cut_vol : reciprocal_cut_vol;
    /*%%%%%%%%%%%%%%%%%%%%% TERMINATION AND CLEANING %%%%%%%%%%%%%%%%%%%%%%%%*/
    plhs[0] = mxCreateDoubleScalar(cutedges);
    plhs[1] = mxCreateNumericArray(2, dims, mxINT64_CLASS, mxREAL);
    plhs[2] = mxCreateDoubleScalar(cutedges / denominator);
    temp = (int64_t *) mxGetPr(plhs[1]);
    temp[0] = denominator;

    free(mask_cut);
    free(reciprocal_mask_cut);

} 
