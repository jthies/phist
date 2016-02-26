/*! \file phist_chol_def.hpp
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

// threshold at which to call a matrix rank deficient
#ifdef SINGTOL
#undef SINGTOL
#endif
#define SINGTOL 10*mt::eps()

// calculates a possibly low rank approximation of a lower cholesky factor of an spd matrix
// higher-precision + pivoting + stable low rank approximation
extern "C" void SUBR(cholesky)(_ST_ *__restrict__ a, lidx_t n, lidx_t lda, lidx_t *perm, int *rank, int* iflag)
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
    diagNorm += st::real(st::conj(d[i])*d[i]);
  }
  if( diagNorm == mt::zero() )
  {
    PHIST_OUT(PHIST_WARNING,"zero diagonal in %s\n", __FUNCTION__);
    diagNorm = mt::eps();
  }

  *rank = 0;
  while(*rank < n)
  {
    // check rank
    _MT_ err = 0;
    for(int i = *rank; i < n; i++)
      err = err + st::abs(d[p[i]]);
//printf("step %d, err %e\n", *rank, err);
    if( err < SINGTOL*diagNorm )
      break;

    int m = *rank;
    *rank = *rank + 1;
    // find next pivot
    {
      int i = m;
      for(int j = m+1; j < n; j++)
        if( st::abs(d[p[j]]) > st::abs(d[p[i]]) )
          i = j;
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
          s-=l[p[m]*n+j]*l[p[i]*n+j];
        }
        s*=div_lmm;
        l[p[i]*n+m] = s;
      }
      // d_p[i] = d_p[i]-l_m,p[i]^2
      {
        d[p[i]] -= l[p[i]*n+m] * l[p[i]*n+m];
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
extern "C" void SUBR(backwardSubst)(const _ST_ *__restrict__ r, lidx_t n, lidx_t ldr, lidx_t *p, int rank,
        _ST_ *__restrict__ x, lidx_t k, lidx_t ldx, int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;

  for(int l = 0; l < k; l++)
  {
    _ST_ newXl[n];
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

    for(int i = 0; i < n; i++)
    {
      x[l*ldx+p[i]] = newXl[i];
    }
  }
}

// apply forward substitution with permuted transposed upper triangular matrix
extern "C" void SUBR(forwardSubst)(const _ST_ *__restrict__ r, lidx_t n, lidx_t ldr, lidx_t *p, int rank,
        _ST_ *__restrict__ x, lidx_t k, lidx_t ldx, int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;

  for(int l = 0; l < k; l++)
  {
    _ST_ newXl[n];
    for(int i = 0; i < rank; i++)
    {
      _ST_ rii_inv=st::one()/r[p[i]*ldr+i];
      // for j = 1,i-1
      // x_p[i],l <- x_p[i],l - r_j,p[i]*x_p[j],l
      for(int j = 0; j < i; j++)
      {
        x[l*ldx+i] -= r[p[j]*ldr+i]*newXl[j];
      }

      // x_p[i] = x_p[i]/r_i,p[i]
      newXl[i] =x[l*ldx+i]*rii_inv;
    }

    // unpermute result
    for(int i = 0; i < n; i++)
    {
      x[l*ldx+i] = newXl[i];
    }
  }
}

extern "C" void SUBR(qb)(_ST_ *__restrict__ a,
                    _ST_ *__restrict__ bi,
                    lidx_t n, lidx_t lda, int *rank, int* iflag)
{
  *iflag=-99;
}

