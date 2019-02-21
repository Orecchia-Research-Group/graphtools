/* C MATLAB function: Pairing

PURPOSE: prepares the flow problem and calls hipr which performs it.
         Obtains hipr output and converts it to MATLAB objects.

USAGE: 
function [flow, cut, matching]= Pairing(G, bisec, cap_add, cap_orig); 

INPUTS: Note: vertex indices start at 1
 -G: a sparse graph
 -bisec: an array of 64bits integers, each representing an index of a node
 in the desired bisection
 - cap_add: a 64 bit integer representing the capacity to put for edges between source
 and bisec and between sink and complement of bisec
 - cap-orig: capacity to assign to edges in G

OUTPUTS:
 - flow: value of flow routed
 - cut: mincut (list of indices), smaller side of mincut is returned
 - matching: demand flow routed between bisec and complement. Note that
 if no matching is required by MATLAB the flow computation does not waste time computing it.
*/

#include "mex.h"
#include "matrix.h"
#include "timer.h"

/* PROTOTYPE
function [flow, cut, matching]= Pairing(G, bisec, cap_add, cap_orig); 
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

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  const mxArray* G;
  long* bisec;
  long cap_add;
  long cap_orig;
  
  long N;
  long M;
  long size_bisec;
  long* col_G;
  long* row_G;
  double* pr_G;
  
  
  long *tails;
  long *heads;
  long *weights;
  

  long i;
  long j;
  long k;
  long h;
  int* mask;
  mxArray* temp;

  long n;
  long m;
  int* output_set;
  long* mheads;
  long* mtails;
  long* mweights;
  long nedges;
  long fflow;
  long size_cut;

  mxArray* matching;
  mxArray* cut;
  mxArray* flow;
  long* cut_pr;

  mwSize dims[] = {1,1};
  long* temp_a;

  float t1, t2;
  int route_flag;
 
  /*  t1 =timer();*/

  if(nrhs != 4 || nlhs > 3 || nlhs < 2)
    mexErrMsgTxt("Error in usage of Pairing.\n");

  /* EXTRACT DATA FROM MATLAB */
  G = prhs[0];
  bisec = (long*) mxGetPr(prhs[1]);
  size_bisec = mxGetM(prhs[1]);
  cap_add = ((long*) mxGetPr(prhs[2]))[0];
  cap_orig = ((long*) mxGetPr(prhs[3]))[0];
  
  N = mxGetM(G);
  mexCallMATLAB(1, &temp, 1, &G, "nnz");
  M = mxGetScalar(temp);
 
  col_G = (long*) mxGetJc(G);
  row_G = (long*) mxGetIr(G);
  pr_G = (double*) mxGetPr(G);

  
  /* CONSTRUCT THE FLOW PROBLEM IN THW REPRESENTATION 
     NEED TO ADD SOURCE / SINK AND RELATIVE EDGES
  */
  
  tails = calloc(sizeof(*tails), M + N);
  heads = calloc(sizeof(*heads), M + N);
  weights = calloc(sizeof(*weights), M + N);
  mask = calloc(sizeof(*mask), N+1);

  for(h=0; h<size_bisec; h++)
    mask[bisec[h]]= 1;
  
  k = 0;


  for(i=0; i < N; i++)
    for(j=col_G[i]; j < col_G[i+1]; j++)
      {
	heads[k] = i+1;
	tails[k] = row_G[j]+1;
	weights[k] = pr_G[j] * cap_orig;

	k++;
      }

  
  for(h=0; h < N; h++) 
    {
      if(mask[h+1] == 1) 
	{
	  heads[k] = h+1;
	  tails[k] = N+1;
	}
      else
	{
	  heads[k] = N+2;
	  tails[k] = h+1;
	}
      weights[k] = cap_add;
      
      k++;
    }
  

  /* CALL HI_PR - modified to output flow - would prefer for hipr to allocate this memory*/


  /* m-arrays are initialized within hi_pr */
  n = N+2;
  m = M+N;

  if(nlhs == 2)
    {
      route_flag = 0;
    }
  else
    route_flag = 1;


  /*  t1 = timer() - t1;*/
  hipr(n, m, tails, heads, weights, N+1, N+2, &output_set, &mheads, &mtails, &mweights, &nedges, &fflow, route_flag);
  /*  t2 = timer();*/
  


  
  /* INITIALIZE MATCHING */
  if(route_flag == 1)
    matching = Sparse(mheads, mtails, mweights, nedges, N);
   
   /* INITIALIZE FLOW */
  flow = mxCreateNumericArray(2, dims, mxINT64_CLASS, mxREAL);
  temp_a = (long*) mxGetPr(flow);
  *temp_a = fflow;  
  

   /* INITIALIZE CUT - recall cut returned is sink side! so need to complement*/
   size_cut = 0;
   for(i=0; i < N; i++)
     if(output_set[i] != 0)
       size_cut++;
   
   size_cut = N - size_cut;

   k = 0;
   dims[0]=size_cut;
   cut = mxCreateNumericArray(2,dims,mxINT64_CLASS, mxREAL);
   cut_pr = (long*) mxGetPr(cut);

   for(i=0; i< N; i++)
     if(output_set[i] == 0) 
       {
	 cut_pr[k]= (long) i+1;
	 k++;
       }

   plhs[0] = flow;
   plhs[1] = cut;
  
   if(route_flag == 1)
     plhs[2] = matching;

   /*   t2 = timer() -t2;
	fprintf(stderr, "Oth tm: %f", t2 + t1);*/
}
 
