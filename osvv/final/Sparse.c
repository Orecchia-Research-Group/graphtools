/*
C MATLAB function: Sparse

PURPOSE: this is a C function to create a MATLAB sparse square matrix.
USAGE: same as matlab sparse command.

NOTES: the inputs are long but the sparse array holds them in double
format as this is the only type supproted by MATLAB with the sparse format.

INPUTS:
 - heads: array of arc heads.
 - tails: array of arc tails.
 - weights: array of arc weights.
 - m: number of nnz entries
 - n: matrix dimension

OUTPUTS:
 - mxArray*: a pointer to the sparse Array
*/

#include <stdint.h>
#include "mex.h"
#include "matrix.h"

mxArray* Sparse(int64_t* heads, int64_t* tails, int64_t* weights, int64_t m, int64_t n )
{
  mxArray* H;
  mxArray* T;
  mxArray* W;

  mxArray* plhs[1];
  mxArray* prhs[5];

  double* hh;
  double* tt;
  double* ww;

  mxArray* temp;
  int i;

  H = mxCreateDoubleMatrix(m,1,mxREAL);
  T = mxCreateDoubleMatrix(m,1,mxREAL);
  W = mxCreateDoubleMatrix(m,1,mxREAL);
 
  hh = mxGetPr(H); 
  tt = mxGetPr(T); 
  ww = mxGetPr(W);
 
  for(i=0; i < m; i++) 
    {
      hh[i] = (double) heads[i];
      tt[i] = (double) tails[i];
      ww[i] = (double) weights[i];
    }

  prhs[0] = H; 
  prhs[1] = T; 
  prhs[2] = W;
  prhs[3] = mxCreateDoubleScalar(n);
  prhs[4] = mxCreateDoubleScalar(n);

  mexCallMATLAB(1, plhs,5, prhs, "sparse");

  mxDestroyArray(H);
  mxDestroyArray(T); 
  mxDestroyArray(W); 
  mxDestroyArray(prhs[3]); 
  mxDestroyArray(prhs[4]);
  
  temp = plhs[0];
  
  return  temp;
}
