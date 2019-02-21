/*
C MATLAB function: sweepcut
USAGE: [(double) ex_num, (int64) ex_den,  (double) ex, (double) side] = sweepcut(sparse matrix (double) G, vector (int64) ordering);

PURPOSE: compute the expansion of the best sweep cut of G

NOTES: 
   - side will be 0 if cut is small side of ordering. Different double ow.
   - G is assumed to be sparse, program will check this
   - ASSUMING C LONG TYPE IS 64 BITS


mexFunction INPUTS;
   nrhs = 2
   nlhs = 4
   
*/
#include "mex.h"
#include "matrix.h"
#include <algorithm>
#include <tr1/unordered_set>
#include <cfloat>
using namespace std;

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  const mxArray* G;
  mwIndex* col;
  mwIndex* row;
  
  double* array_G;
  long* ordering;
  long size_ordering;
  long i, j, v, u;
  long n;
  long half;
  long den;
  
  double cutedges = 0;  

  double min_cutedges;
  long min_den;
  double min_exp;
  double min_side;
  double side;
  double exp;

  mwSize dims[]={1,1}; 
  long* temp;

  /*%%%%%%%%%%%%%%%%% ARGUMENT LOADING &  CHECKING %%%%%%%%%%%%%%%%%%%%%*/

  /* CHECK CORRECT NUMBER OF INPUT/OUTPUTS */
  if(nrhs != 2 || nlhs != 4) 
    mexErrMsgTxt("Error in sweepcut. Incorrect usage.\n");
 
  /* CHECK TYPES */
  if(!mxIsSparse(prhs[0]))
    mexErrMsgTxt("Error in sweepcut. Graph must be sparse.\n");

  if(mxGetClassID(prhs[1]) != mxINT64_CLASS)
    mexErrMsgTxt("Error in sweepcut. Ordering must be of class int64.\n");

  /* SHOULDN'T  I CHECK MATRIX IS SQUARED? AND  POSITIVE ENTRIES */
  
  /* LOAD ARGUMENTS */
  G = prhs[0];
  n = mxGetM(G);
  col = mxGetJc(G);
  row = mxGetIr(G);
  array_G = mxGetPr(G);


  ordering = (long*) mxGetPr(prhs[1]);
  size_ordering = mxGetM(prhs[1]);
  if (size_ordering != n) {
    mexErrMsgTxt("Error in sweepcut. Vector  must have same size as graph.\n"); 
  }



  /*%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN BODY %%%%%%%%%%%%%%%%%%%%%%%%%*/

  min_exp = DBL_MAX;
  half = n/2;

  tr1::unordered_set<long> cutset;

  /* COMPUTE EXPANSION */
  for(i=0; i < n-1; i++) {
    v = ordering[i] - 1;
    cutset.insert(v);
    for(j = col[v]; j < col[v+1]; j++) {  /* loop over u=row[j]'s adjacent to v*/
      u = row[j];
      if(cutset.find(u) == cutset.end()) /* u is not in current sweep, i.e. u is ahead */
	cutedges += array_G[j];
      else
	cutedges -= array_G[j]; /* u is in current sweep , i.e. u is behind*/
    }
    
    if(cutedges < 0) /* for numerical stability */
      cutedges = 0;
    
    if( i+1  > half) { /* SHOULD ONLY HAVE TO DO THIS ONCE */
      den = n-i-1;
      side = 1;
    }
    else {
      den = i+1;
      side = 0;
    }
    
    exp = cutedges/den;

    if(exp < min_exp) {
      min_exp = exp;
      min_den = den;
      min_cutedges = cutedges;
      min_side = side;
    }

    /*    mexPrintf("Exp: %lf, Size: %d. Edges: %lf.\n", exp, den, cutedges); */


  }

  /*%%%%%%%%%%%%%%%%%%%%% TERMINATION AND CLEANING %%%%%%%%%%%%%%%%%%%%%%%%*/
  
  plhs[0] = mxCreateDoubleScalar(min_cutedges); 
  plhs[1] = mxCreateNumericArray(2, dims, mxINT64_CLASS, mxREAL);
  plhs[2] = mxCreateDoubleScalar(min_exp);
  plhs[3] = mxCreateDoubleScalar(side);
  
  temp = (long*) mxGetPr(plhs[1]);
  temp[0] = min_den; 
} 



  


  
