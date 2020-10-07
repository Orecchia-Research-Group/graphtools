/* C MATLAB function: Pairing

PURPOSE: prepares the flow problem and calls hipr which performs it.
         Obtains hipr output and converts it to MATLAB objects.

USAGE: 
function [flow, cut, matching]= Pairing(G, bisec, source_modifier, sink_modifier, original_modifier);

INPUTS: Note: vertex indices start at 1
 -G: a sparse graph
 -bisec: an array of 64bits integers, each representing an index of a node
 in the desired bisection
 - source_modifier: a 64 bit integer representing the capacity to put for edges between source
 and bisec
 - sink_modifier: a 64 bit integer representing the capacity to put for edges between sink and complement of bisec
 - original_modifier: capacity to multiply edges in G
 - lambda: percentage of degree flow that can pass through the node

OUTPUTS:
 - flow: value of flow routed
 - cut: mincut (list of indices), smaller side of mincut is returned
 - matching: demand flow routed between bisec and complement. Note that
 if no matching is required by MATLAB the flow computation does not waste time computing it.
*/

#include <math.h>
#include "mex.h"
#include "matrix.h"
#include "timer.h"

/* PROTOTYPE
function [flow, cut, matching]= Pairing(G, bisec, source_modifier, sink_modifier, original_modifier);
*/

void hipr 
( 
   long ninput, 
   long minput, 
   long *tails, 
   long *heads, 
   long *weights, 
   long s, 
   long t,
   int **output_set,
   long **mheads,
   long **mtails,
   long **mweights,
   long* nedges,
   long* fflow,
   int route_flag
);



mxArray* Sparse(long* heads, long* tails, long* weights, long m, long n );

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    const mxArray *G;
    long *bisec;
    long source_modifier;
    long sink_modifier;
    long original_modifier;
    long internal_modifier;

    long N;
    long M;
    long size_bisec;
    long *col_G;
    long *row_G;
    double *pr_G;
    long *volume;

    long *tails;
    long *heads;
    long *weights;
    long *degrees;
    long sideWeight[2] = {0, 0};
    long sideNom = 1;
    long sideDen = 1;

    long i;
    long j;
    long k;
    long h;
    long reciprocalOffset;
    int *mask;
    mxArray *temp;

    long n;
    long m;
    int *output_set;
    long *mheads = NULL;
    long *mtails = NULL;
    long *mweights = NULL;
    long nedges;
    long fflow;
    long size_cut;
    long reciprocal_size_cut;

    mxArray *matching;
    mxArray *cut, *reciprocalCut;
    mxArray *flow;
    long *cut_pr, *reciprocalCut_pr;

    mwSize dims[] = {1, 1};
    long *temp_a;

    float t1, t2;
    int route_flag;

    /*  t1 =timer();*/

    if (nrhs > 7 || nrhs < 6 || nlhs > 4 || nlhs < 3)
        mexErrMsgTxt("Error in usage of Pairing.\n");


    /* EXTRACT DATA FROM MATLAB */
    G = prhs[0];
    bisec = (long *) mxGetPr(prhs[1]);
    size_bisec = mxGetM(prhs[1]);
    volume = (long *) mxGetPr(prhs[2]);
    source_modifier = ((long *) mxGetPr(prhs[3]))[0];
    sink_modifier = ((long *) mxGetPr(prhs[4]))[0];
    original_modifier = ((long *) mxGetPr(prhs[5]))[0];

    if (nrhs > 6) internal_modifier = ((long *) mxGetPr(prhs[6]))[0];

    N = mxGetM(G);
    mexCallMATLAB(1, &temp, 1, &G, "nnz");
    M = mxGetScalar(temp);

    col_G = (long *) mxGetJc(G);
    row_G = (long *) mxGetIr(G);
    pr_G = (double *) mxGetPr(G);

    reciprocalOffset = N * (nrhs > 6);

    /* CONSTRUCT THE FLOW PROBLEM IN THE REPRESENTATION
       NEED TO ADD SOURCE / SINK AND RELATIVE EDGES
    */

    mask = calloc(sizeof(*mask), N + 1);
    if (!mask) {
        fprintf(stderr, "Error allocating mask\n");
    }

    for (h = 0; h < size_bisec; h++)
        mask[bisec[h]] = 1;

    k = 0;

    // Mixed cut. Nodes 1-N have incoming edges. Nodes N+1-2N have outgoing edges.
    // There is an edge from node k to node k+N with capacity lambda.

    tails = calloc(sizeof(*tails), M + N + reciprocalOffset);
    heads = calloc(sizeof(*heads), M + N + reciprocalOffset);
    weights = calloc(sizeof(*weights), M + N + reciprocalOffset);
    degrees = calloc(sizeof(*degrees), N + 1);

    if (!(tails && heads && weights && degrees)) {
        fprintf(stderr, "Error allocating memory for edge information\n");
    }

    long zero_degree = 0;
    for (i = 0; i < N; i++) {
        for (j = col_G[i]; j < col_G[i + 1]; j++) {
            heads[k] = i + 1;
            tails[k] = row_G[j] + reciprocalOffset + 1;
            weights[k] = pr_G[j] * original_modifier;
            degrees[i + 1] += weights[k];
            k++;
        }
        if (degrees[i + 1] == 0) {
            // fprintf(stderr, "Node %6ld has degree 0.\n", i + 1);
            zero_degree++;
        }
        if (i == 0 && degrees[i + 1] == 0) fprintf(stderr, "u: %ld. j: %ld. w: %lf. index: %ld. original_modifier: %ld\n",
                i + 1, row_G[col_G[i]] + 1, pr_G[col_G[i]], col_G[i], original_modifier);
    }
    if (zero_degree > 0) fprintf(stderr, "There are %ld with zero degree\n", zero_degree);

    long zero_internal = 0;
    for (h = 0; h < reciprocalOffset; h++) {
        heads[k] = h + reciprocalOffset + 1;
        tails[k] = h + 1;
        weights[k] = volume[h] * internal_modifier;
        if (weights[k] <= 0) {
            // fprintf(stderr, "Internal edge for node %6ld was reduced to 0.\n", h + 1);
            zero_internal++;
        }
        k++;
    }
    if (zero_internal > 0) fprintf(stderr, "There are %ld whose internal edges were reduced to zero\n", zero_internal);

    for (h = 0; h < N; h++) {
        if (mask[h + 1] == 1) {
            heads[k] = h + 1;
            tails[k] = N + reciprocalOffset + 1;
            weights[k] = source_modifier * volume[h];
        } else {
            heads[k] = N + reciprocalOffset + 2;
            tails[k] = h + reciprocalOffset + 1;
            weights[k] = sink_modifier * volume[h];
        }

        k++;
    }

    /* CALL HI_PR - modified to output flow - would prefer for hipr to allocate this memory*/

    /* m-arrays are initialized within hi_pr */
    n = N + 2 + reciprocalOffset;
    m = M + N + reciprocalOffset;

    if (nlhs == 3) {
        route_flag = 0;
    } else
        route_flag = 1;
    /*  t1 = timer() - t1;*/
    hipr(n, m, tails, heads, weights, N + reciprocalOffset + 1, N + reciprocalOffset + 2, &output_set, &mheads, &mtails, &mweights, &nedges, &fflow,
         route_flag);
    /*  t2 = timer();*/




    /* INITIALIZE MATCHING */
    if (route_flag == 1) {
        // all the heads in matching are in [N+1, 2N]. To get them to the initial space, subtract N from every element in mheads.
        if (nrhs > 4)
            for (k = 0; k < nedges; k++) {
                if (mheads[k] > N) mheads[k] -= N;
                if (mtails[k] > N) mtails[k] -= N;
            }

        matching = Sparse(mheads, mtails, mweights, nedges, N);
    }

    /* INITIALIZE FLOW */
    flow = mxCreateNumericArray(2, dims, mxINT64_CLASS, mxREAL);
    temp_a = (long *) mxGetPr(flow);
    *temp_a = fflow;


    /* INITIALIZE CUT - recall cut returned is sink side! so need to complement*/
    size_cut = 0;
    for (i = 0; i < N; i++)
        if (output_set[i] == 0)
            size_cut++;

    reciprocal_size_cut = 0;
    for (i = reciprocalOffset; i < N + reciprocalOffset; i++)
        if (output_set[i] != 0)
            reciprocal_size_cut++;

    j = 0;
    k = 0;
    dims[0] = size_cut;
    cut = mxCreateNumericArray(2, dims, mxINT64_CLASS, mxREAL);
    cut_pr = (long *) mxGetPr(cut);
    dims[0] = reciprocal_size_cut;
    reciprocalCut = mxCreateNumericArray(2, dims, mxINT64_CLASS, mxREAL);
    reciprocalCut_pr = (long *) mxGetPr(reciprocalCut);

    for (i = 0; i < N; i++) {
        if (output_set[i] == 0) {
            cut_pr[k] = (long) i + 1;
            k++;
        }
        if (output_set[i+reciprocalOffset] != 0) {
            reciprocalCut_pr[j] = (long) i + 1;
            j++;
        }
    }

   plhs[0] = flow;
   plhs[1] = cut;
   plhs[2] = reciprocalCut;
  
   if(route_flag == 1)
     plhs[3] = matching;

   /*   t2 = timer() -t2;
	fprintf(stderr, "Oth tm: %f", t2 + t1);*/
   free(mask);
   free(heads);
   free(tails);
   free(weights);
   free(degrees);
   if(mheads) free(mheads);
   if(mtails) free(mtails);
   if(mweights) free(mweights);
   free(output_set);
}
 
