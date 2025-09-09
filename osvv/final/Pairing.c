/* C MATLAB function: Pairing

PURPOSE: prepares the flow problem and calls hipr which performs it.
         Obtains hipr output and converts it to MATLAB objects.

USAGE: 
function [flow, cut, matching]= Pairing(G, source_set, sink_set, source_modifier, sink_modifier, original_modifier[, internal_modifier]);

INPUTS: Note: vertex indices start at 1
 -G: a sparse graph
 -source_set: an array of 64bits integers, each representing an index of a node in the left side of the partition
 -sink_set: an array of 64bits integers, each representing an index of a node in the right side of the partition
 - source_modifier: a 64 bit integer representing the capacity to put for edges between source
 and source_set
 - sink_modifier: a 64 bit integer representing the capacity to put for edges between sink and sink_set
 - original_modifier: capacity to multiply edges in G
 - lambda: percentage of degree flow that can pass through the node

OUTPUTS:
 - flow: value of flow routed
 - cut: mincut (list of indices), smaller side of mincut is returned
 - matching: demand flow routed between source_set and sink_set. Note that
 if no matching is required by MATLAB the flow computation does not waste time computing it.
*/

#include <math.h>
#include <string.h>
#include "mex.h"
#include "matrix.h"
#include "timer.h"
#include "flow.h"

/* PROTOTYPE
function [flow, cut, matching]= Pairing(G, source_set, sink_set, source_modifier, sink_modifier, original_modifier[, internal_modifier]);
*/

mxArray* Sparse(int64_t* heads, int64_t* tails, int64_t* weights, int64_t m, int64_t n );

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    const mxArray *G;
    int64_t *source_set;
    int64_t *sink_set;
    int64_t source_modifier;
    int64_t sink_modifier;
    int64_t original_modifier;
    int64_t internal_modifier;

    int64_t N;
    int64_t M;
    int64_t source_set_size;
    int64_t sink_set_size;
    int64_t *col_G;
    int64_t *row_G;
    double *pr_G;
    int64_t *volume;

    int64_t *tails;
    int64_t *heads;
    int64_t *weights;
    int64_t *degrees;

    int64_t i;
    int64_t j;
    int64_t k;
    int64_t h;
    int64_t reciprocalOffset;
    int *source_set_mask;
    int *sink_set_mask;
    mxArray *temp;

    int64_t n;
    int64_t m;
    int64_t *output_set;
    int64_t *mheads = NULL;
    int64_t *mtails = NULL;
    int64_t *mweights = NULL;
    int64_t nedges;
    int64_t fflow;
    int64_t size_cut;
    int64_t reciprocal_size_cut;
    int64_t internal_edge_count;
    int64_t internalNodes;
    int64_t matching_index = 0;

    mxArray *matching;
    mxArray *cut, *reciprocalCut;
    mxArray *flow;
    int64_t *cut_pr, *reciprocalCut_pr;

    mwSize dims[] = {1, 1};
    int64_t *temp_a;

    // float t1, t2;
    int route_flag;

    /*  t1 =timer();*/

    if (nrhs > 9 || nrhs < 8 || nlhs > 4 || nlhs < 3)
        mexErrMsgTxt("Error in usage of Pairing.\n");


    /* EXTRACT DATA FROM MATLAB */
    G = prhs[0];
    source_set = (int64_t *) mxGetPr(prhs[1]);
    source_set_size = mxGetM(prhs[1]);
    sink_set = (int64_t *) mxGetPr(prhs[2]);
    sink_set_size = mxGetM(prhs[2]);
    volume = (int64_t *) mxGetPr(prhs[3]);
    char *matching_algorithm = mxArrayToString(prhs[4]);
    if (!strcmp(matching_algorithm, "dinic"))
        matching_index = 0;
    else if (!strcmp(matching_algorithm, "dynamic"))
        matching_index = 1;
    else
        mexErrMsgTxt("Error in recognizing the matching algorithm");
    source_modifier = ((int64_t *) mxGetPr(prhs[5]))[0];
    sink_modifier = ((int64_t *) mxGetPr(prhs[6]))[0];
    original_modifier = ((int64_t *) mxGetPr(prhs[7]))[0];

    if (nrhs > 8) internal_modifier = ((int64_t *) mxGetPr(prhs[8]))[0];
    else internal_modifier = 1;

#ifdef DEBUG
    fprintf(stderr, "Pairing: source_modifier=%ld\n", source_modifier);
    fprintf(stderr, "Pairing: sink_modifier=%ld\n", sink_modifier);
    fprintf(stderr, "Pairing: original_modifier=%ld\n", original_modifier);
    fprintf(stderr, "Pairing: internal_modifier=%ld\n", internal_modifier);
#endif

    N = mxGetM(G);
    mexCallMATLAB(1, &temp, 1, &G, "nnz");
    M = mxGetScalar(temp);

    col_G = (int64_t *) mxGetJc(G);
    row_G = (int64_t *) mxGetIr(G);
    pr_G = (double *) mxGetPr(G);

    internal_edge_count = 0;
    for (int64_t i = 0; i < N; i++) {
        internal_edge_count += (volume[i] > 0);
    }
    reciprocalOffset = N * (nrhs > 6);
    internalNodes = internal_edge_count * (reciprocalOffset > 0);

    /* CONSTRUCT THE FLOW PROBLEM IN THE REPRESENTATION
       NEED TO ADD SOURCE / SINK AND RELATIVE EDGES
    */

    source_set_mask = calloc(sizeof(*source_set_mask), N + 1);
    if (!source_set_mask) {
        fprintf(stderr, "Error allocating mask\n");
    }

    for (h = 0; h < source_set_size; h++)
        source_set_mask[source_set[h]] = 1;

    sink_set_mask = calloc(sizeof(*sink_set_mask), N + 1);
    if (!sink_set_mask) {
        fprintf(stderr, "Error allocating mask\n");
    }

    for (h = 0; h < sink_set_size; h++)
        sink_set_mask[sink_set[h]] = 1;

    k = 0;

    // Mixed cut. Nodes 1-N have incoming edges. Nodes N+1-2N have outgoing edges.
    // There is an edge from node k to node k+N with capacity lambda.

    tails = calloc(sizeof(*tails), M + internal_edge_count + internalNodes);
    heads = calloc(sizeof(*heads), M + internal_edge_count + internalNodes);
    weights = calloc(sizeof(*weights), M + internal_edge_count + internalNodes);
    degrees = calloc(sizeof(*degrees), N + 1);

    if (!(tails && heads && weights && degrees)) {
        fprintf(stderr, "Error allocating memory for edge information\n");
    }

    int64_t zero_degree = 0;
    for (i = 0; i < N; i++) {
        for (j = col_G[i]; j < col_G[i + 1]; j++) {
            if (i == row_G[j]) continue;
            heads[k] = i + 1;
            tails[k] = row_G[j] + reciprocalOffset * (volume[row_G[j]] > 0) + 1;
            weights[k] = ((int64_t) pr_G[j]) * original_modifier;
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

    int64_t zero_internal = 0;
    for (h = 0; h < reciprocalOffset; h++) {
        if (volume[h] == 0) continue;
        heads[k] = h + reciprocalOffset + 1;
        tails[k] = h + 1;
        weights[k] = volume[h] * internal_modifier;
        if (weights[k] <= 0) {
            // fprintf(stderr, "Internal edge for node %6ld was reduced to 0.\n", h + 1);
            zero_internal++;
        }
        k++;
    }
    if (zero_internal > 0) fprintf(stderr, "There are %ld nodes whose internal edges were reduced to zero\n", zero_internal);

    for (h = 0; h < N; h++) {
        if (volume[h] == 0) continue;
        if (source_set_mask[h + 1] == 1) {
            heads[k] = h + 1;
            tails[k] = N + internalNodes + 1;
            weights[k] = source_modifier * volume[h];
            k++;
        }
        if (sink_set_mask[h + 1] == 1) {
            heads[k] = N + internalNodes + 2;
            tails[k] = h + reciprocalOffset + 1;
            weights[k] = sink_modifier * volume[h];
            k++;
        }
    }

    /* CALL HI_PR - modified to output flow - would prefer for hipr to allocate this memory*/

    /* m-arrays are initialized within hi_pr */
    n = N + 2 + internalNodes;
    m = k;

    if (nlhs == 3) {
        route_flag = 0;
    } else {
        route_flag = 1;
    }
    /*  t1 = timer() - t1;*/
    hipr(n, m, tails, heads, weights, N + internal_edge_count + 1, N + internal_edge_count + 2, &output_set, &mheads, &mtails, &mweights, &nedges, &fflow,
         route_flag, matching_index);
    /*  t2 = timer();*/

#ifdef DEBUG
    fprintf(stderr, "Original nodes: %ld\n", N);
    fprintf(stderr, "Original edges: %ld\n", M);
    fprintf(stderr, "Flow nodes: %ld\n", n);
    fprintf(stderr, "Flow arcs: %ld\n", m);
    fprintf(stderr, "Flow returned: %ld\n", fflow);
#endif

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
    temp_a = (int64_t *) mxGetPr(flow);
    *temp_a = fflow;


    /* INITIALIZE CUT - recall cut returned is sink side! so need to complement*/
    size_cut = 0;
    for (i = 0; i < N; i++)
        if (output_set[i] == 0)
            size_cut++;

    reciprocal_size_cut = 0;
    for (i = internalNodes; i < N + internalNodes; i++)
        if (output_set[i] != 0)
            reciprocal_size_cut++;

    j = 0;
    k = 0;
    dims[0] = size_cut;
    cut = mxCreateNumericArray(2, dims, mxINT64_CLASS, mxREAL);
    cut_pr = (int64_t *) mxGetPr(cut);
    dims[0] = reciprocal_size_cut;
    reciprocalCut = mxCreateNumericArray(2, dims, mxINT64_CLASS, mxREAL);
    reciprocalCut_pr = (int64_t *) mxGetPr(reciprocalCut);
    
    #ifdef DEBUG
        fprintf(stderr, "Preparing to create cut = %ld and reciprocal cut = %ld\n", size_cut, reciprocal_size_cut);
    #endif
    
    for (i = 0; i < N; i++) {
        if (output_set[i] == 0) {
            cut_pr[k] = (int64_t) i + 1;
            k++;
        }
        if (output_set[i + reciprocalOffset * (volume[i] > 0)] != 0) {
            reciprocalCut_pr[j] = (int64_t) i + 1;
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
   free(source_set_mask);
   free(sink_set_mask);
   free(heads);
   free(tails);
   free(weights);
   free(degrees);
   if(mheads) free(mheads);
   if(mtails) free(mtails);
   if(mweights) free(mweights);
   free(output_set);
}
 
