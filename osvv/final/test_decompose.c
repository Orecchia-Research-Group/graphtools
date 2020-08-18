#include "stdio.h"
#include "stdlib.h"
#include "assert.h"
#include "flow.c"


int main() {
    //greater sized tree with n=10 (including source and sink) and m = 14
    node *a1, *a2, *a3, *a4, *a5, *a6, *a7, *a8, *source, *sink;
    
    node *tails[14] = {a1, a2, a3, a4, a5, a6, a6, a8, a8, a7, sink, sink, sink, sink};
    node *heads[14] = {source, source, source, a1, a1, a1, a2, a2, a2, a3, a4, a5, a6, a7};
    node *flows[14] = {2, 1, 2, 0, 2, 1, 1, 2, 0, 1, 2, 2, 1, 2};
    
    long s = 9;
    long t = 10;
    int **output_set[8];
    long **mheads, **mtails, **mweights;
    long *nedges, *fflow;
    int route_flag;
    
    //call of function, hipr that is being tested, :
    hipr(10, 14, tails, heads, flows, s, t, output_set, mheads, mtails, mweights, nedges, fflow, route_flag);
    
    //testing the mheads, mtails, mweights which is what was changed by the decomposition of the paths
    
    //checking if mheads, mweights, mtails are of correct length
    assert (sizeOf(mheads)==12);
    assert (sizeOf(mtails)==12);
    assert (sizeOf(mweights)==12);
    
    
    //the kth path is at the 2kth index
    
    assert(mheads[0]==a_1);
    assert(mtails[0]==a_5);
    assert(mweights[0]==2)
    
    assert(mheads[2]==a_1);
    assert(mtails[2]==a_6);
    assert(mweights[2]==1);
    
    assert(mheads[4]==a_2);
    assert(mtails[4]==a_6);
    assert(mweights[4]==1);
    
    assert(mheads[6]==a_2);
    assert(mtails[6]==a_7);
    assert(mweights[6]==1);
    
    assert(mheads[8]==a_2);
    assert(mtails[8]==a_8);
    assert(mweights[8]==0);
    
    assert(mheads[10]==a_3);
    assert(mtails[10]==a_8);
    assert(mweights[10]==1);
    
}