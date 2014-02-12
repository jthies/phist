#include "mex.h"

/* function x=dkswp(A,sigma,b,x,omega,nrm_ai2)                    */
/* Kaczmarz forward/backward sweep on A'-sigma*I                        */
/* note: you have to pass in A^T to this function!                      */
/* A should be real, but sigma may be complex. b is assumed             */
/* to be real, but x may be complex.                                    */
/* nrm_ai2(i) should be ||A(i,:)-sigma*I||_2^2.                         */
void mexFunction(int nlhs,
                 mxArray *plhs[],
                 int nrhs,
                 const mxArray *prhs[])
{
  int i,j;
  int64_t nrows, nvecs;

  int64_t *irA, *jcA;
  double *prA, *xr,*xi, *br, *bi, *nrm_ai2, *pr_omega, *diagA;
  double *xr_in, *xi_in;
  double *prSigma, *piSigma;
  double omega, sigma_r, sigma_i;
  
  nrows = mxGetM(prhs[0]);
  nvecs = mxGetN(prhs[3]);
  
  if (nvecs!=1)
  {
    fprintf(stderr, "WARNING - handling only the first column of X correctly\n"
                    "(in file %s, line %d)\n",__FILE__,__LINE__);
  }

  /* do we have a complex shift or complex vectors? */
  int is_cmplx=mxIsComplex(prhs[1]);
  is_cmplx|=mxIsComplex(prhs[2]);
  is_cmplx|=mxIsComplex(prhs[3]);

  /* copy x to the output arg, make it complex if necessary */
  if (is_cmplx && !mxIsComplex(prhs[3]))
  {
    plhs[0]=mxCreateDoubleMatrix(nrows, nvecs, mxCOMPLEX);
    xr = mxGetPr(plhs[0]);
    xi = mxGetPi(plhs[0]);
    xr_in = mxGetPr(prhs[3]);
    for (i=0;i<nrows;i++)
    {
      xr[i]=xr_in[i];
      xi[i]=0.0;
    }
  }
  else
  {
    plhs[0]=mxDuplicateArray(prhs[3]);
    xr = mxGetPr(plhs[0]);
    xi = mxGetPi(plhs[0]);
  }

  if (xi==NULL)
  {
    xi=(double*)malloc(nrows*sizeof(double));
    for (i=0;i<nrows;i++)
    {
      xi[i]=0.0;
    }
  }

  irA = mxGetIr(prhs[0]);
  jcA = mxGetJc(prhs[0]);
  prA = mxGetPr(prhs[0]);
  prSigma=mxGetPr(prhs[1]);
  sigma_r=prSigma[0];
  piSigma=mxGetPi(prhs[1]);
  if (piSigma)
  {
    sigma_i=piSigma[0];
  }
  else
  {
    sigma_i=0.0;
  }
  br = mxGetPr(prhs[2]);
  bi = mxGetPi(prhs[2]);
  pr_omega = mxGetPr(prhs[4]);
  nrm_ai2 = mxGetPr(prhs[5]);
  omega=pr_omega[0];

  /* real formulation of complex system: 
             | A-sr I     si I | |xr|   |br|
        Ax = |  -si I   A-sr I | |xi| = | 0|
  */
  for (i=0; i<nrows; i++)
  {
    /* contribution from rhs (b) and shifted diagonal */
    double scal_r=-br[i]-sigma_r*xr[i]+sigma_i*xi[i];
    double scal_i=-sigma_r*xi[i]-sigma_i*xr[i];
    if (bi) scal_i-=bi[i];
    /* contributions from the unshifted matrix A */
    for (j=jcA[i]; j<jcA[i+1]; j++)
    {
      scal_r+=prA[j]*xr[irA[j]];
      scal_i+=prA[j]*xi[irA[j]];
    }
    /* note: we require the nrm_ai2 array to already include the shift */
    scal_r*=omega/nrm_ai2[i];
    scal_i*=omega/nrm_ai2[i];

    /* diagonal (cross-)terms */
    xr[i]+=scal_r*sigma_r+scal_i*sigma_i;
    xi[i]+=scal_i*sigma_r - scal_r*sigma_i;
    for (j=jcA[i]; j<jcA[i+1]; j++)
    {
      xr[irA[j]] -= scal_r*prA[j];
      xi[irA[j]] -= scal_i*prA[j];
    }
  }

  for (i=nrows-1; i>=0; i--)
  {
    double scal_r=-br[i]-sigma_r*xr[i]+sigma_i*xi[i];
    double scal_i=-sigma_r*xi[i]-sigma_i*xr[i];
    if (bi) scal_i-=bi[i];
    for (j=jcA[i]; j<jcA[i+1]; j++)
    {
      scal_r+=prA[j]*xr[irA[j]];
      scal_i+=prA[j]*xi[irA[j]];
    }
    scal_r*=omega/nrm_ai2[i];
    scal_i*=omega/nrm_ai2[i];
    xr[i]+=scal_r*sigma_r+scal_i*sigma_i;
    xi[i]+=scal_i*sigma_r - scal_r*sigma_i;
    for (j=jcA[i]; j<jcA[i+1]; j++)
    {
      xr[irA[j]] -= scal_r*prA[j];
      xi[irA[j]] -= scal_i*prA[j];
    }
  }

if (mxGetPi(plhs[0])==NULL)
{
  free(xi);
}

return;
}
