#include "mex.h"

/* function d=nrm_ai2(A,sigma)
calculate the 2-norm (squared) of each column of A-sigma*I 
A should be real, but sigma may be complex.
*/

void mexFunction(int nlhs,
                 mxArray *plhs[],
                 int nrhs,
                 const mxArray *prhs[])
{
  int i,j;
  int64_t nrows, ncols;

  int64_t *irA, *jcA;
  double *prA, *nrm_ai2;
  double *prSigma,*piSigma;
  double sigma_r, sigma_i;
  
  nrows = mxGetM(prhs[0]);
  ncols = mxGetN(prhs[0]);

  /* copy x to the output arg */
  plhs[0]=mxCreateDoubleMatrix(nrows,1,false);

  irA = mxGetIr(prhs[0]);
  jcA = mxGetJc(prhs[0]);
  prA = mxGetPr(prhs[0]);
  nrm_ai2 = mxGetPr(plhs[0]);

  prSigma = mxGetPr(prhs[1]);
  piSigma = mxGetPi(prhs[1]);
  sigma_r=prSigma[0];
  if (piSigma)
  {
    sigma_i=piSigma[0];
  }
  else
  {
    sigma_i=0.0;
  }

  for (i=0; i<nrows; i++)
  {
    double sum=sigma_r*sigma_r+sigma_i*sigma_i;
    for (j=jcA[i]; j<jcA[i+1]; j++)
    {
      if (irA[j]==i)
      {
        sum-=2.0*prA[j]*sigma_r;
      }
    }
    for (j=jcA[i]; j<jcA[i+1]; j++)
    {
      sum+=prA[j]*prA[j];
    }
    nrm_ai2[i]=sum;
  }

return;
}
