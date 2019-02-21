/*
C MATLAB function: Farey
 
USAGE: [num,den] = Farey(fractnum, fracden, p)

PURPOSE: compute the best fractional approximation larger than fracnum/fracden with
numerator <= p.
 
NOTES:  p < 100000

*/


#include "mex.h"
#include "matrix.h"
#include <math.h>

void mexFunction (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
  long num_in;
  long den_in;
  long p;
  long num_out;
  long den_out;

  mwSize dims[] = {1,1};
  long* temp;

  long error;
  long a,b,c,d;
  long h, k;

  /*%%%%%%%%%%%%%%%%% ARGUMENT LOADING &  CHECKING %%%%%%%%%%%%%%%%%%%%%*/


  /* CHECK CORRECT NUMBER OF INPUT/OUTPUTS */
  if(nrhs != 3 || nlhs != 2) 
    mexErrMsgTxt("Error in Farey. Incorrect usage.\n");
 
  /* CHECK TYPES */
  if(mxGetClassID(prhs[0]) != mxINT64_CLASS)
    mexErrMsgTxt("Error in Farey. Numerator must be of class int64.\n");

  if(mxGetClassID(prhs[1]) != mxINT64_CLASS)
    mexErrMsgTxt("Error in Farey. Denominator  must be of class int64.\n");
 
  if(mxGetClassID(prhs[2]) != mxINT64_CLASS)
    mexErrMsgTxt("Error in Farey. Precision must be of class int64.\n");

  
  num_in =  ((long*) mxGetPr(prhs[0]))[0];
  den_in =  ((long*) mxGetPr(prhs[1]))[0]; 
  p = ((long*) mxGetPr(prhs[2]))[0];

  if(num_in < 0 || den_in < 0)
    mexErrMsgTxt("Error in Farey. Numerators and denominators must be positive.\n");

  if( p > 100000 || p < 0)
    mexErrMsgTxt("Error in Farey. Precision out of range.\n");
      

  /*%%%%%%%%%%%%%%%%%%%%%%%%% MAIN BODY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  
  

  a = 0;
  b = 1;
  c = ceil( ((double) num_in)/den_in);
  d = 1;

  for( ; ; )
    {
      h = a + c;
      k = b + d;
      
      if( h > p) 
	{
	  h = c;
	  k = d;
	  break;
	}
   
      error = h * den_in - k * num_in; /* error is positive if h/k > num/den */

      if(error == 0)
	break;

      if(error < 0)
	{
	  a = h;
	  b = k;
	}
      
      if(error > 0)
	{
	  c = h;
	  d = k;
	}
    }


  /*%%%%%%%%%%%%%%%%%%%%% TERMINATION AND CLEANING %%%%%%%%%%%%%%%%%%%%%%%%*/
  
  plhs[0] = mxCreateNumericArray(2, dims, mxINT64_CLASS, mxREAL);
  plhs[1] = mxCreateNumericArray(2, dims, mxINT64_CLASS, mxREAL);
  
  temp = (long*) mxGetPr(plhs[0]);
  *temp = h;
  
  temp = (long*) mxGetPr(plhs[1]);
  *temp = k;
  
}



