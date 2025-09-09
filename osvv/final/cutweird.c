/*
C MATLAB function: cutweird

USAGE: [(double) wr_num, (int64) wr_den,  (double) wr] = cutexp(sparse matrix (double) G, vector (int64) cut, vector (int64) bisec);

PURPOSE: compute the weird ratio of the cut with respect to the near bisection bisec.

NOTES:
   - if weird ratio is negative, return abs.value (i.e. weirdratio of complement cut)
   - G is assumed to be undirected, no check for that
   - G is assumed to be sparse, program will check this
   - ASSUMING C int64_t TYPE IS 64 BITS


mexFunction INPUTS;
   nrhs = 3
   nlhs = 3

*/
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include "mex.h"
#include "matrix.h"
#include "farey.h"

#define abs(x) (x) > 0 ? (x) : -(x)

int64_t gcd(int64_t a, int64_t b)
{
    // Everything divides 0
    if ((!a) || (!b))
        return a + b;
    // a is greater
    if (a > b)
        return gcd(a % b, b);
    return gcd(a, b % a);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {       // TODO: Add weight vector argument
    const mxArray *G;
    mwIndex *col;
    mwIndex *row;

    double *array_G;
    int64_t *cut;
    int64_t *reciprocal_cut;
    int64_t size_cut;
    int64_t reciprocal_size_cut;
    int64_t *source_set;
    int64_t *sink_set;
    int64_t source_set_size;
    int64_t sink_set_size;
    int64_t *weight;
    int64_t source_set_volume = 0;
    int64_t sink_set_volume = 0;
    int64_t w_bisec;
    int64_t w_recip;
    int64_t p = 10000;
    int64_t i, j;
    int64_t n;
    double cutedges = 0;
    int64_t denominator = 0;
    int64_t lamda_num;
    int64_t lamda_den;

    int *mask_cut;
    int *reciprocal_mask_cut;
    int *source_set_mask;
    int *sink_set_mask;

    mwSize dims[] = {1, 1};
    int64_t *temp;

    /*%%%%%%%%%%%%%%%%% ARGUMENT LOADING &  CHECKING %%%%%%%%%%%%%%%%%%%%%*/

    //fprintf(stderr, "cutweird: nrhs = %d nlhs = %d\n", nrhs, nlhs);

    /* CHECK CORRECT NUMBER OF INPUT/OUTPUTS */
    if (nrhs != 8 || nlhs != 3)
        mexErrMsgTxt("Error in cutweird. Incorrect usage.\n");

    /* CHECK TYPES */
    if (!mxIsSparse(prhs[0]))
        mexErrMsgTxt("Error in cutweird. Graph must be sparse.\n");

    if (mxGetClassID(prhs[1]) != mxINT64_CLASS)
        mexErrMsgTxt("Error in cutweird. Cut must be of class int64.\n");

    if (mxGetClassID(prhs[2]) != mxINT64_CLASS)
        mexErrMsgTxt("Error in cutweird. Reciprocal cut must be of class int64.\n");

    if (mxGetClassID(prhs[3]) != mxINT64_CLASS)
        mexErrMsgTxt("Error in cutweird. sink_set must be of class int64.\n");

    if (mxGetClassID(prhs[4]) != mxINT64_CLASS)
        mexErrMsgTxt("Error in cutweird. sink_set must be of class int64.\n");

    if (mxGetClassID(prhs[5]) != mxINT64_CLASS)
        mexErrMsgTxt("Error in cutweird. Weight must be of class int64.\n");
    
    if (mxGetClassID(prhs[6]) != mxINT64_CLASS)
        mexErrMsgTxt("Error in cutweird. lamda_num must be of class int64.\n");
    
    if (mxGetClassID(prhs[7]) != mxINT64_CLASS)
        mexErrMsgTxt("Error in cutweird. lamda_den must be of class int64.\n");




    /* LOAD ARGUMENTS */
    G = prhs[0];
    n = mxGetM(G);
    col = mxGetJc(G);
    row = mxGetIr(G);
    array_G = mxGetPr(G);

    cut = (int64_t *) mxGetPr(prhs[1]);
    size_cut = mxGetM(prhs[1]);
    reciprocal_cut = (int64_t *) mxGetPr(prhs[2]);
    reciprocal_size_cut = mxGetM(prhs[2]);

    source_set = (int64_t *) mxGetPr(prhs[3]);
    source_set_size = mxGetM(prhs[3]);
    sink_set = (int64_t *) mxGetPr(prhs[4]);
    sink_set_size = mxGetM(prhs[4]);
    weight = (int64_t *) mxGetPr(prhs[5]);
    lamda_num = ((int64_t *) mxGetPr(prhs[6]))[0];
    lamda_den = ((int64_t *) mxGetPr(prhs[7]))[0];
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

    if (!(source_set_mask = calloc(sizeof(*source_set_mask), n + 1))) mexErrMsgTxt("Error allocating memory in cutweird.");
    else for (i = 0; i < source_set_size; i++) source_set_mask[source_set[i] - 1] = 1;

    if (!(sink_set_mask = calloc(sizeof(*sink_set_mask), n + 1))) mexErrMsgTxt("Error allocating memory in cutweird.");
    else for (i = 0; i < sink_set_size; i++) sink_set_mask[sink_set[i] - 1] = 1;

    for (i = 0; i < n; i++) {
        if (source_set_mask[i]) source_set_volume += weight[i];
        if (sink_set_mask[i]) sink_set_volume += weight[i];
    }

    w_bisec = source_set_volume;
    w_recip = sink_set_volume;
    
    if ((sink_set_volume > source_set_volume) && (sink_set_volume / (double) source_set_volume - 1 < 0.0001)) {
        w_bisec = w_recip = 1;
#ifdef DEBUG
        fprintf(stderr, "Source and sink sets have almost the same weight: fabs(((double)source_set_volume) /  sink_set_volume - 1) = %lf < 1e-4\n", fabs(((double)source_set_volume) / sink_set_volume - 1));
#endif
    } else {
        cfarey(source_set_volume, sink_set_volume, p, &w_bisec, &w_recip);
#ifdef DEBUG
        fprintf(stderr, "Calling farey(%ld, %ld, %ld, %ld, %ld)\n", source_set_volume, sink_set_volume, p, w_bisec, w_recip);
#endif
    }
#ifdef DEBUG
    fprintf(stderr, "Initial vol(L)=%ld vol(R)=%ld ratio=%lf\n", source_set_volume, sink_set_volume, fabs(((double)source_set_volume) /  sink_set_volume - 1));
    fprintf(stderr, "vol(L)=%ld vol(R)=%ld\n", w_bisec, w_recip);
#endif

    /* COMPUTE WEIRD RATIO*/
    /* COMPUTE EDGES CUT */

    for (i = 0; i < n; i++) {
        if (mask_cut[i] && source_set_mask[i]) {                 // π(S && A)
            denominator += weight[i] * lamda_den * w_recip;
        }

        if (mask_cut[i] && sink_set_mask[i] && !reciprocal_mask_cut[i]) {     // - π(L && B), L = S \ T
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

    // denominator = e2 * size_intersect - size_cut + (size_cut + reciprocal_size_cut - n - size_overlap_intersect);
    if (denominator < 0)
        denominator = denominator * (-1);

    /*%%%%%%%%%%%%%%%%%%%%% TERMINATION AND CLEANING %%%%%%%%%%%%%%%%%%%%%%%%*/

    int64_t g = gcd(abs(cutedges), abs(denominator));
    if (g > 0) {
        cutedges /= g;
        denominator /= g;
    }

#ifdef DEBUG
    fprintf(stderr, "cutedges=%lf denominator=%ld\n", cutedges, denominator);
#endif

    cfarey(cutedges, denominator, p, &w_bisec, &w_recip);

    plhs[0] = mxCreateDoubleScalar(w_bisec);
    plhs[1] = mxCreateDoubleScalar(w_recip);
    plhs[2] = mxCreateDoubleScalar(w_bisec / (double) w_recip);

    free(mask_cut);
    free(reciprocal_mask_cut);
    free(source_set_mask);
    free(sink_set_mask);

}







