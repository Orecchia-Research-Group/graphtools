/*
C MATLAB function: cutexp

USAGE: [(double) ex_num, (int64) ex_den,  (double) ex] = cutexp(sparse matrix (double) G, vector (int64) cut);

PURPOSE: compute the expansion of cut in G

NOTES:
   - cut will be complemented if larger than half the vertices;
   - G is assumed to be undirected, no check for that
   - G is assumed to be sparse, program will check this
   - ASSUMING C LONG TYPE IS 64 BITS


mexFunction INPUTS;
   nrhs = 2
   nlhs = 3

*/
#include "mex.h"
#include "matrix.h"

#define min(a, b) ( ( (a) < (b) ) ? a : b )

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  const mxArray* G;
  mwIndex* col;
  mwIndex* row;

  double* array_G;
  long* cut;
  long size_cut;
  long i, j;
  long n;
  double cutedges = 0;

  int* mask_cut;
  mxArray* exp;

  mwSize dims[]={1,1};
  long* temp;

  /*%%%%%%%%%%%%%%%%% ARGUMENT LOADING &  CHECKING %%%%%%%%%%%%%%%%%%%%%*/

  /* CHECK CORRECT NUMBER OF INPUT/OUTPUTS */
  if(nrhs != 2 || nlhs != 3)
    mexErrMsgTxt("Error in cutexp. Incorrect usage.\n");

  /* CHECK TYPES */
  if(!mxIsSparse(prhs[0]))
    mexErrMsgTxt("Error in cutexp. Graph must be sparse.\n");

  if(mxGetClassID(prhs[1]) != mxINT64_CLASS)
    mexErrMsgTxt("Error in cutexp. Cut must be of class int64.\n");


  /* LOAD ARGUMENTS */
  G = prhs[0];
  n = mxGetM(G);
  col = mxGetJc(G);
  row = mxGetIr(G);
  array_G = mxGetPr(G);

  cut = (long*) mxGetPr(prhs[1]);
  size_cut = mxGetM(prhs[1]);

  /* PREPARE MASK - COMPLEMENT CUT IF NECESSARY */
  mask_cut = calloc(sizeof(*mask_cut), n+1);
  if (mask_cut == NULL)
    mexErrMsgTxt("Error allocating memory in cutexp.");


  for(i=0; i < size_cut; i++)
      mask_cut[cut[i]] = 1;


  // if(size_cut > n/2)    /* COMPLEMENTATION */
  //   {
  //   for(i=1; i < n+1 ; i++)
  //     mask_cut[i] = !mask_cut[i];
  //
  //   size_cut = n - size_cut;
  //
  //   j = 0;
  //   for(i=1; i < n+1; i++)
  //     {
	// if(mask_cut[i]) {
	//   cut[j] = i;
	//   j++;
	// }
  //     }
  //   }


  /* COMPUTE EXPANSION */
  for(i=0; i < size_cut; i++)
      for(j = col[cut[i]-1]; j < col[cut[i]]; j++)
	  if(mask_cut[row[j] + 1] == 0)
	    cutedges += array_G[j];

  /*%%%%%%%%%%%%%%%%%%%%% TERMINATION AND CLEANING %%%%%%%%%%%%%%%%%%%%%%%%*/

  plhs[0] = mxCreateDoubleScalar(cutedges);
  plhs[1] = mxCreateNumericArray(2, dims, mxINT64_CLASS, mxREAL);
  plhs[2] = mxCreateDoubleScalar(cutedges/min(size_cut, n - sizecut));
  temp = (long*) mxGetPr(plhs[1]);
  temp[0] = size_cut;

  free(mask_cut);

}
