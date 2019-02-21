/*
C MATLAB function: sweepcutquot
USAGE: [(double) cut_num, (int64) cut_den,  (double) cutquot, (double) lrange, (double) urange] = sweepcut(sparse matrix (double) G, vector (int64) ordering, vector (int64) quotweight);

PURPOSE: compute the cut quotient of the best sweep cut of G. works for weighted graphs too.

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
  long* qw;
  long weight;
  long cw;

  long i, j, v, u;
  long n;
  long half;
  long den;
  
  double cutedges = 0;  

  double min_cutedges;
  long min_den;
  double min_score;
  double min_side;
  double min_index;
  double score;
  double side;
  double lrange;
  double urange;

  mwSize dims[]={1,1}; 
  long* temp;

  /*%%%%%%%%%%%%%%%%% ARGUMENT LOADING &  CHECKING %%%%%%%%%%%%%%%%%%%%%*/

  /* CHECK CORRECT NUMBER OF INPUT/OUTPUTS */
  if(nrhs != 3 || nlhs != 5) 
    mexErrMsgTxt("Error in sweepcut. Incorrect usage.\n");
 
  /* CHECK TYPES */
  if(!mxIsSparse(prhs[0]))
    mexErrMsgTxt("Error in sweepcut. Graph must be sparse.\n");

  if(mxGetClassID(prhs[1]) != mxINT64_CLASS)
    mexErrMsgTxt("Error in sweepcut. Ordering must be of class int64.\n");

  /* SHOULDN'T  I CHECK MATRIX IS SQUARED? AND  POSITIVE ENTRIES */


  if(mxGetClassID(prhs[2]) != mxINT64_CLASS)
    mexErrMsgTxt("Error in cutquot. quotweight must be of class int64.\n");

  if(mxGetM(prhs[2]) != mxGetM(prhs[0]))
     mexErrMsgTxt("Error in cutquot. quotweight must be of size equal to the number of vertices in the graph.");

  
  /* LOAD ARGUMENTS */
  G = prhs[0];
  n = mxGetM(G);
  col = mxGetJc(G);
  row = mxGetIr(G);
  array_G = mxGetPr(G);



  ordering = (long*) mxGetPr(prhs[1]);
  size_ordering = mxGetM(prhs[1]);
  
  /*  if (size_ordering != n) {
    mexErrMsgTxt("Error in sweepcut. ordering  must have same size as graph.\n"); 
    }
  */

  qw = (long*) mxGetPr(prhs[2]);

  /*  COMPUTE TOTAL WEIGHT */
  weight = 0;
  for(i=0; i < n; i++) {
    if(qw[i] < 0) 
      mexErrMsgTxt("Error in quotweight. Weights must be nonnegative");
    weight += qw[i];
  }





  /*%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN BODY %%%%%%%%%%%%%%%%%%%%%%%%%*/

  min_score = DBL_MAX;
  half = weight/2;
  cw = 0;

  tr1::unordered_set<long> cutset;

  /* COMPUTE SCORE */
  for(i=0; i < size_ordering-1; i++) {
    v = ordering[i] - 1;
    cutset.insert(v);
    cw += qw[v];

    for(j = col[v]; j < col[v+1]; j++) {  /* loop over u=row[j]'s adjacent to v*/
      u = row[j];
      if(cutset.find(u) == cutset.end()) /* u is not in current sweep, i.e. u is ahead */
	cutedges += array_G[j];
      else
	cutedges -= array_G[j]; /* u is in current sweep , i.e. u is behind*/
    }
    
    if(cutedges < 0) /* for numerical stability */
      cutedges = 0;
    
    if( cw  > half) { /* SHOULD ONLY HAVE TO DO THIS ONCE */
      den = weight - cw;
      side = 1;
    }
    else {
      den = cw;
      side = 0;
    }
    
    score = cutedges/den;

    if(score < min_score) {
      min_score = score;
      min_den = den;
      min_cutedges = cutedges;
      min_side = side;
      min_index = i;
    }

    /*    mexPrintf("Exp: %lf, Size: %d. Edges: %lf.\n", score, den, cutedges); */


  }

  if(min_side == 0) {
    lrange = 1;
    urange = min_index + 1;
  }
  else {
    lrange = min_index + 2;
    urange = n;
  }

  /*%%%%%%%%%%%%%%%%%%%%% TERMINATION AND CLEANING %%%%%%%%%%%%%%%%%%%%%%%%*/
  
  plhs[0] = mxCreateDoubleScalar(min_cutedges); 
  plhs[1] = mxCreateNumericArray(2, dims, mxINT64_CLASS, mxREAL);
  plhs[2] = mxCreateDoubleScalar(min_score);
  plhs[3] = mxCreateDoubleScalar(lrange);
  plhs[4] = mxCreateDoubleScalar(urange);

  
  temp = (long*) mxGetPr(plhs[1]);
  temp[0] = min_den; 
} 



  


  
