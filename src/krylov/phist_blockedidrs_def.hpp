/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/

#ifdef M_ELEM
#undef M_ELEM
#endif

// helper function: trace_dot(x, y) = trace(x'*y)
void SUBR(trace_dot)(TYPE(const_mvec_ptr) x, TYPE(const_mvec_ptr) y, _ST_* result, int* iflag)
{
  int nvec;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(x,&nvec,iflag),*iflag);
  _ST_ dot[nvec];
  PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(x, y, dot, iflag), *iflag);
  *result=0;
  for (int i=0; i<nvec; i++) *result+=dot[i];
  return;
}

void SUBR(frob_norm)(TYPE(const_mvec_ptr) x, _MT_* result, int *iflag)
{
#include "phist_std_typedefs.hpp"
  _ST_ dot;
  PHIST_CHK_IERR(SUBR(trace_dot)(x,x,&dot,iflag),*iflag);
  *result = std::sqrt(st::real(dot));
}

// helper function: pdots[i] = trace_dot(P[i],w), for i=0...s-1
void SUBR(P_dot)(TYPE(mvec_ptr) P[], TYPE(const_mvec_ptr) w, int s, _ST_* pdots, int* iflag)
{
  for (int i=0; i<s; i++) PHIST_CHK_IERR(SUBR(trace_dot)(P[i], w, pdots+i, iflag),*iflag);
}

// modified Gram-Schmidt orthogonalization sweep to make P[i], i=1..s mutually orthogonal
// all P[i] must have the same number of columns, namely numSys.
void SUBR(mgs)(TYPE(mvec_ptr) P[], int s, int numSys, int* iflag)
{
  _ST_ alpha=0;
  int nvec;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(P[0],&nvec,iflag),*iflag);
  double nrms[nvec];
  for (int j=0; j<s; j++)
  {
    for (int k=0; k<j; k++)
    {
      PHIST_CHK_IERR(SUBR(trace_dot)(P[k], P[j], &alpha, iflag), *iflag);
      PHIST_CHK_IERR(SUBR(mvec_add_mvec)(-alpha, P[k], 1, P[j], iflag), *iflag);
    }
    PHIST_CHK_IERR(SUBR(mvec_normalize)(P[j], nrms,iflag), *iflag);
  }
}


// implementation of IDR(s) on several systems simultaneously.
// The implementation follows Martin van Gijzens Fortran code found here: http://homepage.tudelft.nl/1w5b5/idrs-software.html
//
extern "C" void SUBR(blockedIDRs_iterate)(TYPE(const_linearOp_ptr) Aop, TYPE(const_linearOp_ptr) Pop,
        TYPE(const_mvec_ptr) rhs, TYPE(mvec_ptr) sol, TYPE(const_mvec_ptr) V,
        int numSys, int* nIter, _MT_ tol, int s, int* iflag)
{
#include "phist_std_typedefs.hpp"
  using mvec_owner=phist::MvecOwner<ST>;
  using sdMat_owner=phist::SdMatOwner<ST>;

  *iflag = 0;
  if (numSys==0) return; // do not appear in timing stats
  PHIST_ENTER_FCN(__FUNCTION__);

  int maxIter = (*nIter)>0 ? *nIter: 1000;

  bool inispace=false; //TODO: Martins implementation allows the user to provide an initial space,
                       //      this sounds like a very good idea in the context of Jacobi-Davidson

  int n_omega=0;
  ST *user_omega = nullptr; //TODO: original implementation allows passing in an array of relaxation parameters,
                               //      we just use the default omega=0.7
  ST omega = 1.0;
  *iflag=0;

  if( Aop == nullptr || rhs == nullptr || sol == nullptr || s<1 )
  {
    *iflag = PHIST_INVALID_INPUT;
    return;
  }

  PHIST_CHK_IERR(*iflag=(Pop==nullptr)?0:PHIST_NOT_IMPLEMENTED,*iflag);

  // check dimensions
  {
    int ncrhs,ncsol;
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(rhs,&ncrhs,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(sol,&ncsol,iflag),*iflag);
    // note: we actually don't care about the dimensions of sol as long as it is large 
    // enough since we just work with views in this function. However, we issue a warning 
    // (positive iflag) if the array dimensions are larger than expected so that the user 
    // doesn't unexpectedly end up with empty rows/cols
    if (ncsol<numSys || ncrhs<numSys) {*iflag = -1; return;}
    if (ncsol>numSys || ncrhs>numSys)
    {
      PHIST_SOUT(PHIST_VERBOSE,"REMARK: input vectors/matrix to IDRs are larger than necessary.\n");
      PHIST_SOUT(PHIST_VERBOSE,"        sol has %d columns (expecting %d)\n",                                        ncsol, ncrhs);
    }
  }

  // Allocate multivectors (mvecs) and arrays of s multivectors.
  // We do this by creating C++ objects that own the pointers (MvecOwners),
  // and then retrieve their raw pointers to conveniently pass them to phist
  // functions. The C++ objects (marked by a leading underscore) are deleted
  // at the end of the scope (function), cleaning up the memory allocated by
  // the kernel library.

  // create local variables:
  phist_const_map_ptr map=nullptr;
  phist_const_comm_ptr comm=nullptr;
  PHIST_CHK_IERR(SUBR(mvec_get_map)(rhs, &map, iflag), *iflag);
  PHIST_CHK_IERR(phist_map_get_comm(map,&comm,iflag),*iflag);
  mvec_owner _x, _b;
  mvec_owner _v(map, numSys,iflag); PHIST_CHK_IERR("allocate v",*iflag);
  mvec_owner _r(map, numSys,iflag); PHIST_CHK_IERR("allocate r",*iflag);
  mvec_owner _t(map, numSys,iflag); PHIST_CHK_IERR("allocate t",*iflag);
  mvec_owner _P[s], _U[s], _G[s];

  // to be able to call phist functions conveniently, get raw pointers
  mvec_ptr r=_r.get(), t=_t.get(), v=_v.get();

  mvec_ptr x=nullptr;
  mvec_ptr b=nullptr;
  PHIST_CHK_IERR(SUBR(mvec_view_block)(sol, &x, 0, numSys-1, iflag),*iflag);
  _x.set(x);
  PHIST_CHK_IERR(SUBR(mvec_view_block)((mvec_ptr)rhs, &b, 0, numSys-1, iflag),*iflag);
  _b.set(b);

  mvec_ptr P[s], G[s], U[s];

  for (int i=0; i<s; i++)
  {
     PHIST_CHK_IERR(SUBR(mvec_create)(&U[i], map, numSys, iflag), *iflag);
     _U[i].set(U[i]);
     PHIST_CHK_IERR(SUBR(mvec_create)(&G[i], map, numSys, iflag), *iflag);
     _G[i].set(G[i]);
     PHIST_CHK_IERR(SUBR(mvec_create)(&P[i], map, numSys, iflag), *iflag);
     _P[i].set(P[i]);
  }

  *nIter = 0;

  // compute initial residual
  MT normR, normR0, normB, tolB, tolR0;
  ST alpha[s], beta[s], gamma[s];
  PHIST_CHK_IERR(SUBR(mvec_add_mvec)(1,b,0,r,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(linearOp_apply)(-1, Aop, x, 1, r, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(frob_norm)(b, &normB, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(frob_norm)(r, &normR, iflag), *iflag);
  normR0 = normR;
  tolB  = tol * normB;
  tolR0 = tol * normR0;

  PHIST_SOUT(PHIST_VERBOSE, "IDR iter %d:\t%8.4e\t%8.4e\n", *nIter, normR/normR0, normR);

  //check if the initial solution is not already a solution within the prescribed
  //tolerance for all of the systems. If it is for one or a few, we do at least one
  //iteration so that the outer JD loop progresses.
  if (normR < tolB)
  {
    *iflag=0;
    return;
  }


  // setting the first column of P to r reproduces BiCGStab for IDR(1)
  PHIST_CHK_IERR(SUBR(mvec_add_mvec)(1,r,0,P[0], iflag), *iflag);
  for (int i=1; i<s; i++)
  {
    PHIST_CHK_IERR(SUBR(mvec_random)(P[i], iflag), *iflag);
  }

  // orthogonalize the s P vectors mutually
  PHIST_CHK_IERR(SUBR(mgs)(P, s, numSys, iflag), *iflag);
  ST kappa = 0.7;

  sdMat_owner _M(s,s,comm,iflag); PHIST_REPORT_IERR("memory allocation failed",*iflag);
  sdMat_ptr M = _M.get();
  PHIST_CHK_IERR(SUBR(sdMat_put_value)(M,0,iflag),*iflag);
  ST* M_raw=nullptr;
  int ldM;
  PHIST_CHK_IERR(SUBR(sdMat_extract_view)(M,&M_raw,&ldM,iflag),*iflag);
#ifdef PHIST_SDMATS_ROW_MAJOR
#define M_ELEM(II,JJ) M_raw[ (II)*ldM+(JJ) ]
#else
#define M_ELEM(II,JJ) M_raw[ (JJ)*ldM+(II) ]
#endif

  // This concludes the initialisation phase

  int ii=0, jj=0; // inner iteration counters

  // Main iteration loop, build G-spaces:

  while (  *nIter < maxIter )  // start of iteration loop
  {
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Generate s vectors in G_j
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++

    // New right-hand side for small system: f[i] = trace(dot(P[i],r))
    ST f[s];
    PHIST_CHK_IERR(SUBR(P_dot)(P, r, s, f, iflag), *iflag);

    for (int k=0; k<s; k++)
    {
      ii++;

      // Compute new v
      PHIST_CHK_IERR(SUBR(mvec_add_mvec)(1, r, 0, v, iflag), *iflag);
      if ( jj > 0 )
      {
        // Solve small system (Note: M is lower triangular) and make v orthogonal to P:
        for (int i=k; i<s; i++)
        {
          gamma[i] = f[i];
          for (int j = k; j<i; j++)
          {
            gamma[i] -= M_ELEM(i,j)*gamma[j];
          }
          gamma[i] = gamma[i]/M_ELEM(i,i);
          //v = v - gamma[i]*G[i]
          PHIST_CHK_IERR(SUBR(mvec_add_mvec)(-gamma[i], G[i], 1, v, iflag), *iflag);
        }

        // Compute new U[k]
        if (Pop!=nullptr)
        {
          PHIST_CHK_IERR(SUBR(linearOp_apply)(omega,Pop,v,0,t,iflag),*iflag);
        }
        else
        {
          PHIST_CHK_IERR(SUBR(mvec_add_mvec)(omega,v,0,t,iflag),*iflag);
        }
        for (int i=k; i<s; i++)
        {
          //t = t + gamma[i]*U[i]
          PHIST_CHK_IERR(SUBR(mvec_add_mvec)(gamma[i],U[i],1,t,iflag),*iflag);
        }
        //U[k] = t
        std::swap(U[k], t);
      }
      else if (!inispace)
      {
        // Updates for the first s iterations (in G_0):
        //U[k] = v/M1
        if (Pop!=nullptr)
        {
          PHIST_CHK_IERR(SUBR(linearOp_apply)(1,Pop,v,0,U[k],iflag),*iflag);
        }
        else
        {
          PHIST_CHK_IERR(SUBR(mvec_add_mvec)(1,v,0,U[k],iflag),*iflag);
        }
      }

      // Compute new G[k], G[k] is in space G_j
      PHIST_CHK_IERR(SUBR(linearOp_apply)(1,Aop,U[k],0,G[k],iflag),*iflag);

      // Bi-Orthogonalise the new basis vectors: 
      ST mu[s];
      PHIST_CHK_IERR(SUBR(P_dot)( P, G[k], s, mu, iflag), *iflag);
      for (int i = 0; i<k; i++)
      {
        alpha[i] = mu[i];
        for (int j = 0; j<i; j++)
        {
          alpha[i] -= M_ELEM(i,j)*alpha[j];
        }
        alpha[i] /= M_ELEM(i,i);
        //G[k] = G[k] - G[i]*alpha[i]
        PHIST_CHK_IERR(SUBR(mvec_add_mvec)(-alpha[i],G[i],1,G[k],iflag),*iflag);
        //U[k] = U[k] - U[i]*alpha[i]
        PHIST_CHK_IERR(SUBR(mvec_add_mvec)(-alpha[i],U[i],1,U[k],iflag),*iflag);
        //mu(k:s)  = mu(k:s)  - M(k:s,i)*alpha(i)
        for (int j=k; j<s; j++) mu[j] -= M_ELEM(j,i)*alpha[i];
      }
      for (int j=k; j<s; j++) M_ELEM(j,k) = mu[j];

      // Compute Hessenberg matrix?
      /*
      if ( out_H .and. ii <= nritz .and. k  > 1 ) &
          H(ii-k+1:ii-1,ii) =  alpha(1:k-1)/beta(1:k-1)
       */
      // Break down?
      if ( std::abs(M_ELEM(k,k)) <= mt::eps() )
      {
        *iflag=-3;
        PHIST_REPORT_IERR("breakdown encountered in IDR(s)",*iflag);
      }

      // Make r orthogonal to p_i, i = 1..k, update solution and residual 
      beta[k] = f[k]/M_ELEM(k,k);
      //r = r - beta(k)*G(:,:,k)
      PHIST_CHK_IERR(SUBR(mvec_add_mvec)(-beta[k],G[k],1,r,iflag),*iflag);
      ///x = x + beta(k)*U(:,:,k)
      PHIST_CHK_IERR(SUBR(mvec_add_mvec)(beta[k],U[k],1,x,iflag),*iflag);

      // New f = P'*r (first k  components are zero)
      if ( k < s )
      {
        for (int j=k+1; j<s; j++) f[j]   = f[j] - beta[k]*M_ELEM(j,k);
      }

        // Compute Hessenberg matrix?
        /*
        if ( out_H .and. ii <= nritz ) then
            H(ii,ii) = 1./beta(k)
            l = max(1,ii-s)
            H(l+1:ii+1,ii) = (H(l+1:ii+1,ii) - H(l:ii,ii))
            H(l:ii+1,ii)   = H(l:ii+1,ii)/om
         end if
         */
      // Check for convergence
      PHIST_CHK_IERR(SUBR(frob_norm)(r, &normR, iflag), *iflag);
      (*nIter)++;
      PHIST_SOUT(PHIST_VERBOSE, "IDR iter %d:\t%8.4e\t%8.4e\n", *nIter, normR/normR0, normR);
      if ( normR < tolB || normR < tolR0 )
      {
        *iflag = 0;
        return;
      }
      else if ( *nIter == maxIter )
      {
        *iflag = 1;
        return;
      }
    }// k-loop: Now we have computed s+1 vectors in G_j

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Compute first residual in G_j+1
    //!+++++++++++++++++++++++++++++++++++++++++++++++++++++++

    // Update G-space counter
    jj++;

    // Compute first residual in G_j+1
    // Note: r is already perpendicular to P so v = r

    // Preconditioning:
    if (Pop!=nullptr)
    {
      PHIST_CHK_IERR(SUBR(linearOp_apply)(1,Pop,r,0,v,iflag),*iflag);
    }
    else
    {
      PHIST_CHK_IERR(SUBR(mvec_add_mvec)(1,r,0,v,iflag),*iflag);
    }
    //t = A*v
    PHIST_CHK_IERR(SUBR(linearOp_apply)(1,Aop,v,0,t,iflag),*iflag);

    // Computation of a new omega
    if ( user_omega )
    {
      omega = user_omega[jj%n_omega];
    }
    else
    {
      // 'Maintaining the convergence':
      MT nr, nt;
      PHIST_CHK_IERR(SUBR(frob_norm)(r,&nr, iflag), *iflag);
      PHIST_CHK_IERR(SUBR(frob_norm)(t,&nt, iflag), *iflag);
      ST tr;
      PHIST_CHK_IERR(SUBR(trace_dot)(t,r,&tr,iflag),*iflag);
      MT rho = st::abs(tr/(nt*nr));
      omega=tr/(nt*nt);
      if ( rho < kappa )
      {
        omega = omega*kappa/rho;
      }
    }
    if ( std::abs(omega) <= mt::eps() )
    {
      *iflag = 3;
      PHIST_REPORT_IERR("stagnation detected in IDR(s)",*iflag);
    }

    // Update solution and residual
    // r = r - om*t 
    PHIST_CHK_IERR(SUBR(mvec_add_mvec)(-omega, t, 1, r, iflag), *iflag);
    // x = x + om*v
    PHIST_CHK_IERR(SUBR(mvec_add_mvec)(omega, v, 1, x, iflag), *iflag);

    // Check for convergence
    PHIST_CHK_IERR(SUBR(frob_norm)(r,&normR,iflag),*iflag);
    (*nIter)++;
    PHIST_SOUT(PHIST_VERBOSE, "IDR iter %d:\t%8.4e\t%8.4e\n", *nIter, normR/normR0, normR);

    if ( normR < tolB || normR < tolR0)
    {
         *iflag = 0;
    }
    else if ( *nIter == maxIter )
    {
         *iflag = 1;
    }
  } // end of while loop

  // Set output parameters
  /* TODO
   r = b - A*x
   normr = FROB_NORM(r)

   if ( info == 0 .and. normr > tolb ) info = 2
   if ( out_iterations ) iterations = iter
   if ( out_relres )     relres=normr/normb
   if ( out_flag )       flag = info
*/
}

extern "C" void SUBR( IDRs ) (TYPE(const_linearOp_ptr) Aop, TYPE(const_linearOp_ptr) Pop,
		TYPE(const_mvec_ptr) rhs, TYPE(mvec_ptr) sol, int* nIter, _MT_ tol, int s, int* iflag)
{
#include "phist_std_typedefs.hpp"
  *iflag = 0;
  PHIST_ENTER_FCN(__FUNCTION__);
  int num_sol,num_rhs;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(rhs,&num_rhs,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(sol,&num_sol,iflag),*iflag);
  PHIST_CHK_IERR(num_sol==num_rhs?0: PHIST_INVALID_INPUT, *iflag);
  PHIST_CHK_IERR(SUBR(blockedIDRs_iterate)(Aop, Pop, rhs, sol, nullptr, num_rhs, nIter, tol, s, iflag),*iflag);
}


