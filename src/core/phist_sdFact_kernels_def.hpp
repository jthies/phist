/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
/*! \file phist_sdFact_kernels_def.hpp
 * Implementation of some serial LAPaCK like functions with standard precision
 * \author "Melven Roehrig-Zoellner" <Melven.Roehrig-Zoellner@DLR.de>
 * \author "Jonas Thies" <Jonas.Thies@DLR.de>
 */
#include "phist_config.h"
#include "phist_macros.h"
#include "phist_typedefs.h"

//TODO: this file was dapted from the high-precision variant and the implementation is *not finished*

#ifdef PHIST_SDMATS_ROWMAJOR
#error "functionality not implemented for row-major sdMats"
#endif

// calculates a possibly low rank approximation of a lower cholesky factor of an spd matrix
// higher-precision + pivoting + stable low rank approximation
extern "C" void SUBR(cholesky)(_ST_ *__restrict__ a, phist_lidx n, phist_lidx lda, phist_lidx *perm, int *rank, 
        _MT_ rankTol, int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_FCN(__FUNCTION__);

  // permutation
  int p[n];
  
  for(int i = 0; i < n; i++)
    p[i] = i;
  // constructed L
  _ST_ l[n*n], lC[n*n];
  for(int i = 0; i < n*n; i++)
  {
    l[i] = st::zero();
  }
  // diagonal entries
  _ST_ d[n];
  _MT_ diagNorm = mt::zero();
  for(int i = 0; i < n; i++)
  {
    d[i]=a[i*lda+i];
    diagNorm = std::max(diagNorm,st::abs(d[i]));
  }
  if( diagNorm <= rankTol )
  {
    PHIST_OUT(PHIST_WARNING,"all diagonal entries tiny in %s\n", __FUNCTION__);
    // use absolute criterion for small diagonals (will give rank 0 below)
    diagNorm = mt::one();
  }

  *rank = 0;
  while(*rank < n)
  {
    // check rank
    _MT_ err = mt::zero();
    for(int i = *rank; i < n; i++)
    {
      err = std::max(err,st::abs(d[p[i]]));
    }
//printf("step %d, err %e\n", *rank, err);
    if( err <= rankTol*diagNorm ) break;

    int m = *rank;
    *rank = *rank + 1;
    // find next pivot
    {
      int i = m;
      for(int j = m+1; j < n; j++)
      {
        if( st::abs(d[p[j]]) > st::abs(d[p[i]]) )
        {
          i = j;
        }
      }
      // swap p[i] p[m]
      int tmp = p[i];
      p[i] = p[m];
      p[m] = tmp;
//printf("pivot %d, perm", i);
//for(int j = 0; j < n; j++)
//  printf(" %d",p[j]);
    }
//printf("\n");

    // l_m,p[m] = sqrt(d_p[m])
    l[p[m]*n+m]=std::sqrt(d[p[m]]);
    _ST_ div_lmm=st::one()/l[p[m]*n+m];
//printf("m=%d,p[m]=%d: d_pm = %e, l_m,pm = %e, div_lmm = %e\n", m, p[m], l[p[m]*n+m], div_lmm);

    for(int i = m+1; i < n; i++)
    {
      // l_m,p[i] = 1/l_m,p[m] * ( a_p[m],p[i] - sum_j=0^m-2 l_j,p[m]*l_j,p[i] )
      {
        _ST_ s = a[p[i]*lda+p[m]];
        for(int j = 0; j < m; j++)
        {
          _ST_ lj;
          s-=st::conj(l[p[m]*n+j])*l[p[i]*n+j];
        }
        s*=div_lmm;
        l[p[i]*n+m] = s;
      }
      // d_p[i] = d_p[i]-l_m,p[i]^2
      {
        d[p[i]] -= st::conj(l[p[i]*n+m]) * l[p[i]*n+m];
//printf("d[p[i=%d]] <- %e\n", i,-s);
      }
    }
  }

  // store result in a
  for(int i = 0; i < n; i++)
  {
    perm[i] = p[i];
    for(int j = 0; j < n; j++)
    {
      a[i*lda+j] = l[i*n+j];
    }
  }
  *iflag=0;
}


// apply backward substitution with permuted upper triangular matrix to k vectors in col-major storage
extern "C" void SUBR(backwardSubst)(const _ST_ *__restrict__ r, phist_lidx n, phist_lidx ldr, phist_lidx *p, int rank,
        _ST_ *__restrict__ x, phist_lidx k, phist_lidx ldx, int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;

  for(int l = 0; l < k; l++)
  {
    _ST_ newXl[n];
    for (int i= n-1; i>=rank; i--)
    {
      newXl[i]=st::zero();
    }
    for(int i = rank-1; i >= 0; i--)
    {
      _ST_ rii_inv=st::one()/r[p[i]*ldr+i];
      // for j = i+1..n
      // x_i,l <- x_i,l - r_i,p[j]*x_p[j],l
      for(int j = i+1; j < n; j++)
      {
        x[l*ldx+i] -= r[p[j]*ldr+i]*newXl[j];
      }
//printf("x_p[i=%d],l=%d-sum...: %e\n", i, l, -s);

      // x_p[i] = x_p[i]/r_i,p[i]
      newXl[i] =x[l*ldx+i]*rii_inv;
//printf("new x_p[i=%d],l=%d: %e\n", i, l, -s);
    }

    for(int i = 0; i < rank; i++)
    {
      x[l*ldx+p[i]] = newXl[i];
    }
    for (int i=rank; i<n; i++)
    {
      x[l*ldx+p[i]] = st::zero();
    }
  }
}

// apply forward substitution with permuted transposed upper triangular matrix
extern "C" void SUBR(forwardSubst)(const _ST_ *__restrict__ r, phist_lidx n, phist_lidx ldr, phist_lidx *p, int rank,
        _ST_ *__restrict__ x, phist_lidx k, phist_lidx ldx, int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;

  for(int l = 0; l < k; l++)
  {
    _ST_ newXl[n];
    for(int i = 0; i < rank; i++)
    {
      _ST_ s=x[l*ldx+p[i]];
      // for j = 1,i-1
      // x_p[i],l <- x_p[i],l - r_j,p[i]*x_p[j],l
      for(int j = 0; j < i; j++)
      {
        s -= r[p[i]*ldr+j]*newXl[j];
      }

      // x_p[i] = x_p[i]/r_i,p[i]
      newXl[i] =s/r[p[i]*ldr+i];
    }

    // unpermute result
    for(int i = 0; i < n; i++)
    {
      x[l*ldx+i] = newXl[i];
    }
    for(int i = rank; i < n; i++)
    {
      x[l*ldx+p[i]] = st::zero();
    }
  }
}

extern "C" void SUBR(qb)(_ST_ *__restrict__ a,
                    _ST_ *__restrict__ bi,
                    _MT_ *D,
                    phist_lidx n, phist_lidx lda, int *rank, _MT_ rankTol, int* iflag)
{
#include "phist_std_typedefs.hpp"
  // compute sqrt(diag(A)) and its inverse
  _ST_ d[n], di[n];
  for (int i=0; i<n; i++)
  {
    // prevent nans
    if( st::real(a[i*lda+i]) <= rankTol )
    {
      d[i] = 0.;
      di[i] = 0.;
    }
    else
    {
      d[i]=std::sqrt(a[i*lda+i]);
      di[i]=_ST_(1)/d[i];
    }
    if (D!=NULL) D[i]=st::real(d[i]);
  }
  
  // scale input matrix from left and right => diagonal elements 1
  for (int j=0; j<n; j++)
  {
    for (int i=0; i<n; i++)
    {
      a[j*lda+i]*=di[i]*di[j];
    }
  }

  // compute eigenpairs of the scaled matrix A^
  _MT_ w[n];
  #ifdef IS_COMPLEX
    PHIST_CHK_IERR(*iflag=PHIST_LAPACKE(heevd)
        (SDMAT_FLAG, 'V' , 'U', n, (mt::blas_cmplx_t*)a, lda, w),*iflag);
#else
    PHIST_CHK_IERR(*iflag=PHIST_LAPACKE(syevd)
        (SDMAT_FLAG, 'V' , 'U', n, a, lda, w),*iflag);
#endif


  // determine rank of input matrix and 1/sqrt(w)
  _MT_ wi[n];
  *rank=(int)n;
  
  // set w=sqrt(w) and wi=1/sqrt(w)
  for(int i=0; i<n; i++)
  {
    PHIST_SOUT(PHIST_INFO,"w[%d]=%e\n",i,w[i]);
    if (w[i]<rankTol)
    {
      (*rank)--;
      w[i]=_MT_(0);
      wi[i]=_MT_(0);
      PHIST_SOUT(PHIST_INFO,"<-0\n");
    }
    else
    {
      w[i] = std::sqrt(w[i]);
      wi[i]= _MT_(1)/w[i];
    }
  }

  // scale the eigenvector matrix, B^ <- Dinv * B * Einv
  // As B is orthonormal, the inverse of B^ is given by
  // E*B'*D (returned in bi, biC)
  for(int i=0; i<n; i++)
  {
    for(int j=0;j<n;j++)
    {
      if (bi!=NULL)
      {
        bi[i*lda+j]=a[j*lda+i]*d[i]*w[j];
      }
      a[j*lda+i]*=di[i]*wi[j];
    }
  }
  
  // finally we reverse the order of the scaled eigenvectors because
  // this way those associated with 0 eigenvalues will appear last instead
  // of first. This is important for our orthogonalization schemes, which
  // may randomize the last few columns in case of (near) singularity (see
  // flag PHIST_ORTHOG_RANDOMIZE_NULLSPACE)
  for (int i=0; i<n; i++)
  {
    for (int j=0; j<n/2; j++)
    {
      int k=n-j-1;
      std::swap(a[j*lda+i],a[k*lda+i]);
      if (bi!=NULL)
      {
        std::swap(bi[i*lda+j],bi[i*lda+k]);
      }
    }
  }
  
  *iflag=0;
}

