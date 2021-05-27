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
#include <stdio.h>
#include <math.h>
#include "mex.h"
#include "matrix.h"
#include "farey.h"

#define abs(x) (x) > 0 ? (x) : -(x)

long gcd(long a, long b)
{
    // Everything divides 0
    if ((!a) || (!b))
        return a + b;
    // a is greater
    if (a > b)
        return gcd(a - b, b);
    return gcd(a, b - a);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {       // TODO: Add weight vector argument
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
    long *weight;
    long weight_bisec = 0;
    long weight_recip = 0;
    long w_bisec;
    long w_recip;
    long p = 10000;
    long i, j;
    long n;
    double cutedges = 0;
    long denominator = 0;
    long lamda_num;
    long lamda_den;

    int *mask_cut;
    int *reciprocal_mask_cut;
    int *mask_bisec;

    mwSize dims[] = {1, 1};
    long *temp;

    /*%%%%%%%%%%%%%%%%% ARGUMENT LOADING &  CHECKING %%%%%%%%%%%%%%%%%%%%%*/

    //fprintf(stderr, "cutweird: nrhs = %d nlhs = %d\n", nrhs, nlhs);

    /* CHECK CORRECT NUMBER OF INPUT/OUTPUTS */
    if (nrhs != 7 || nlhs != 3)
        mexErrMsgTxt("Error in cutweird. Incorrect usage.\n");

    /* CHECK TYPES */
    if (!mxIsSparse(prhs[0]))
        mexErrMsgTxt("Error in cutweird. Graph must be sparse.\n");

    if (mxGetClassID(prhs[1]) != mxINT64_CLASS)
        mexErrMsgTxt("Error in cutweird. Cut must be of class int64.\n");

    if (mxGetClassID(prhs[2]) != mxINT64_CLASS)
        mexErrMsgTxt("Error in cutweird. Reciprocal cut must be of class int64.\n");

    if (mxGetClassID(prhs[3]) != mxINT64_CLASS)
        mexErrMsgTxt("Error in cutweird. Bisec must be of class int64.\n");

    if (mxGetClassID(prhs[4]) != mxINT64_CLASS)
        mexErrMsgTxt("Error in cutweird. Weight must be of class int64.\n");
    
    if (mxGetClassID(prhs[5]) != mxINT64_CLASS)
        mexErrMsgTxt("Error in cutweird. lamda_num must be of class int64.\n");
    
    if (mxGetClassID(prhs[6]) != mxINT64_CLASS)
        mexErrMsgTxt("Error in cutweird. lamda_den must be of class int64.\n");




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
    weight = (long *) mxGetPr(prhs[4]);
    lamda_num = ((long *) mxGetPr(prhs[5]))[0];
    lamda_den = ((long *) mxGetPr(prhs[6]))[0];
    if (lamda_num < 0) {
        lamda_den = 1l;
    }

    // fprintf(stderr, "lamda_num = %ld, lamda_den = %ld\n", lamda_num, lamda_den);


    /*%%%%%%%%%%%%%%%%%%%%%%%%% MAIN BODY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

    /* PREPARE CUT MASK  AND COMPUTE INTERSECT*/
    if (!(mask_cut = calloc(sizeof(*mask_cut), n + 1))) mexErrMsgTxt("Error allocating memory in cutweird.");
    else for (i = 0; i < size_cut; i++) mask_cut[cut[i] - 1] = 1;

    if (!(reciprocal_mask_cut = calloc(sizeof(*reciprocal_mask_cut), n + 1))) mexErrMsgTxt("Error allocating memory in cutweird.");
    else for (i = 0; i < reciprocal_size_cut; i++) reciprocal_mask_cut[reciprocal_cut[i] - 1] = 1;

    if (!(mask_bisec = calloc(sizeof(*mask_bisec), n + 1))) mexErrMsgTxt("Error allocating memory in cutweird.");
    else for (i = 0; i < size_bisec; i++) mask_bisec[bisec[i] - 1] = 1;

    for (i = 0; i < n; i++) {
        if (mask_bisec[i]) weight_bisec += weight[i];
        else weight_recip += weight[i];
    }

    w_bisec = weight_bisec;
    w_recip = weight_recip;
    
    if (abs(((double)weight_bisec) /  weight_recip - 1) < 1e-4) {
        w_bisec = w_recip = 1;
    } else {
        farey(weight_recip, weight_bisec, p, &w_recip, &w_bisec);
    }
#ifdef DEBUG
    fprintf(stderr, "Initial vol(L)=%ld vol(R)=%ld ratio=%lf\n", weight_bisec, weight_recip, abs(((double)weight_bisec) /  weight_recip - 1));
    fprintf(stderr, "vol(L)=%ld vol(R)=%ld\n", w_bisec, w_recip);
#endif

    /* COMPUTE WEIRD RATIO*/
    /* COMPUTE EDGES CUT */

    for (i = 0; i < n; i++) {
        if (mask_cut[i] && mask_bisec[i]) {                 // π(S && A)
            denominator += weight[i] * lamda_den * w_recip;
        }

        if (mask_cut[i] && (!mask_bisec[i]) && (!reciprocal_mask_cut[i])) {     // - π(L && !A), L = S \ T
            denominator -= weight[i] * lamda_den * w_bisec;
        }

        if (mask_cut[i] && reciprocal_mask_cut[i] && (lamda_num > 0)) {     // λ π(C)
            cutedges += weight[i] * w_recip * lamda_num;
        }

        if (mask_cut[i] && !reciprocal_mask_cut[i]) {                   // w(E(L, R))
            for (j = col[i]; j < col[i + 1]; j++) {
                if (!mask_cut[row[j]]) {
                    cutedges += array_G[j] * lamda_den * w_recip;
                }
            }
        }
        if (cutedges < 0) {
            printf("Overflow of cutedges detected\n");
        }
    }

#ifdef DEBUG
    fprintf(stderr, "cutedges=%lf denominator=%ld\n", cutedges, denominator);
#endif

    // denominator = e2 * size_intersect - size_cut + (size_cut + reciprocal_size_cut - n - size_overlap_intersect);
    if (denominator < 0)
        denominator = denominator * (-1);

    /*%%%%%%%%%%%%%%%%%%%%% TERMINATION AND CLEANING %%%%%%%%%%%%%%%%%%%%%%%%*/

    long g = gcd(abs(cutedges), abs(denominator));
    if (g > 0) {
        cutedges /= g;
        denominator /= g;
    }
    plhs[0] = mxCreateDoubleScalar(cutedges);
    plhs[1] = mxCreateNumericArray(2, dims, mxINT64_CLASS, mxREAL);
    plhs[2] = mxCreateDoubleScalar(cutedges / (double) denominator);
    temp = (long *) mxGetPr(plhs[1]);
    *temp = denominator;

    free(mask_cut);
    free(reciprocal_mask_cut);
    free(mask_bisec);

}







