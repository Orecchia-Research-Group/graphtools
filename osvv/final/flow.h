/*
C function: hipr.

INPUTS: Note that vertex indices go from 1 to n.
 - n: number of vertices
 - m: number of arcs
 - tails: pointer to array of tails of the arcs. tails[i] is the tail of the (i-1)th arc.
 - heads: pointer to array of heads of the arcs.
 - weights: pointer to array of weights of the arcs.
 - s: index of source. (should be n-1)
 - t: index of sink . (should be n)
 - output_set: pointer to array of size n-2 (no sink and source) filled in by hipr to be mask for mincut.
 - mheads: pointer to array of heads of arcs of routed matching.
 - mtails: pointer to array of tails of arcs of routed matching.
 - mweights: pointer to array of weights of arcs of routed matching.
 - fflow: pointer to long which becomes equal to flow routed.
*/

#ifndef FLOW_SEEN
#define FLOW_SEEN


void hipr(
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
        long *nedges,
        long *fflow,
        int route_flag
);

#endif