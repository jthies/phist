#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <complex.h>

#include "scamac_config.h"
#ifdef SCAMAC_EXTERNAL_CBLAS
#include SCAMAC_EXTERNAL_CBLAS
#else
#include "scamac_cblas.h"
#endif
#ifdef SCAMAC_EXTERNAL_CLAPACK
#include SCAMAC_EXTERNAL_CLAPACK
#else
#include "scamac_clapack.h"
#endif

#include "scamac_lanczos.h"


ScamacErrorCode scamac_lanczos_ev_mat_real(const scamac_sparsemat_st *sm, double tol, double *ev1, double *ev2, double *eps1, double *eps2) {
  ScamacErrorCode err;

  if (!sm || !ev1 || !ev2 || !eps1 || !eps2) {
    return SCAMAC_ENULL;
  }
  if (sm->valtype != SCAMAC_VAL_REAL) {
    return SCAMAC_EINVALID;
  }
  if (sm->nr != sm->nc) {
    return SCAMAC_EINVALID;
  }

  //maximal number of Lanczos iterations (not const .. because it's passed to LAPACK routine)
  int nl=200;

  int nv = sm->nc;

  // two vectors (long)
  double *x, *y;
  //tridiagonal matrix elements (short/few)
  double *d,*e,*z,*work,*dhp,*ehp;

  x=malloc(nv * sizeof *x);
  if (!x) {
    return SCAMAC_EMALLOCFAIL;
  }
  y=malloc(nv * sizeof *y);
  if (!y) {
    return SCAMAC_EMALLOCFAIL;
  }
  d=malloc(nl * sizeof *d);
  if (!d) {
    return SCAMAC_EMALLOCFAIL;
  }
  e=malloc(nl * sizeof *e);
  if (!e) {
    return SCAMAC_EMALLOCFAIL;
  }
  dhp=malloc(nl * sizeof *dhp);
  if (!dhp) {
    return SCAMAC_EMALLOCFAIL;
  }
  ehp=malloc(nl * sizeof *ehp);
  if (!ehp) {
    return SCAMAC_EMALLOCFAIL;
  }
  z=malloc(nl*nl * sizeof *z);
  if (!z) {
    return SCAMAC_EMALLOCFAIL;
  }
  work=malloc((2*nl-2) * sizeof *work);
  if (!work) {
    return SCAMAC_EMALLOCFAIL;
  }

  //random vector
  //dm_random_vector(nv,x);
  int idist = 2;
  int iseed[4] = {4,5,6,7};
  dlarnv_(&idist,iseed,&nv,x);
  //normalize
  double dotx = cblas_dnrm2(nv,x,1);
  cblas_dscal(nv,1.0/dotx,x,1);

  int il;
  int happybreak = 0;

  for (il=0; il<nl; il++) {
    if (il==0) {
      err=scamac_sparsemat_mvm(sm,x,y, 1.0,0.0,0.0);
      if (err) {
        return err|SCAMAC_EINTERNAL;
      }
      d[0]=cblas_ddot(nv,x,1,y,1);     //  d[0]=dot_product(nv,x,y);
      cblas_daxpy(nv,-d[0],x,1,y,1); // Y=Y-d(1)*X
      e[0]=cblas_dnrm2(nv,y,1); // dot_product(Y,Y))
      if ( fabs(e[0]) < DBL_EPSILON ) {
        happybreak=1;
      } else {
        cblas_dscal(nv,1.0/e[0],y,1); // Y=Y/e(1)
      }
    } else if ( il%2 == 0 ) {
      err=scamac_sparsemat_mvm(sm,x,y, 1.0, -e[il-1],0.0);
      if (err) {
        return err|SCAMAC_EINTERNAL;
      }
      d[il] = cblas_ddot(nv,x,1,y,1); // d(IL)=dot_product(X,Y)
      cblas_daxpy(nv,-d[il],x,1,y,1); // Y=Y-d(IL)*X
      e[il]=cblas_dnrm2(nv,y,1); // sqrt(ddot_dot_product(Y,Y))
      if ( fabs(e[il])< DBL_EPSILON ) {
        happybreak=1;
      } else {
        cblas_dscal(nv,1.0/e[il],y,1);
      }
    } else {
      err=scamac_sparsemat_mvm(sm,y,x, 1.0, -e[il-1],0.0);
      if (err) {
        return err|SCAMAC_EINTERNAL;
      }
      d[il] = cblas_ddot(nv,y,1,x,1);      //     d(IL)=dot_product(Y,X)
      cblas_daxpy(nv,-d[il],y,1,x,1);      //     X=X-d(IL)*Y
      e[il]=cblas_dnrm2(nv,x,1); //     e(IL)=sqrt(dot_product(X,X))
      if ( fabs(e[il]) < DBL_EPSILON ) {
        happybreak=1;
      } else {
        cblas_dscal(nv,1.0/e[il],x,1);  // X=X/e(IL)
      }
    }

    // determine actual eigenvalue estimates = extremal Ritz values
    if (il==0) {
      *ev1=d[0];
      *ev2=d[0];
      *eps1=0.0;
      *eps2=0.0;
    } else {
      // ! --- needs Z be set to 1 initially? Check.
      //dm_identity(nl,1.0,z);
      cblas_dcopy(il+1,d,1,dhp,1);
      cblas_dcopy(il+1,e,1,ehp,1);
      char compz='I';
      int ilplus1 = il+1;
      int INFO;
      dsteqr_(&compz, &ilplus1,dhp,ehp,z,&nl,work,&INFO);
      if (INFO) {
        return SCAMAC_EFAIL|SCAMAC_EINTERNAL;
      }
      *ev1=dhp[0];
      *ev2=dhp[il];
      //error bounds from vector   (note: e(IL) is not touched in DSTEQR)
      *eps1=fabs(z[il+nl*0])*fabs(e[il]);
      *eps2=fabs(z[il+nl*il])*fabs(e[il]);
      if (*ev2-*ev1>sqrt(DBL_EPSILON)) { // accuracy achieved?
        // relative error
        if ( (*eps1/(*ev2-*ev1)<tol) && (*eps2/(*ev2-*ev1)<tol) ) {
          break;
        }
      } else {
        //absolute error
        if ( (*eps1<tol) && (*eps2<tol) ) {
          break;
        }
      }
    }

    if (happybreak) break;

  }

  free(x);
  free(y);
  free(d);
  free(e);
  free(dhp);
  free(ehp);
  free(z);
  free(work);


  if (il >= nl) {
    return SCAMAC_ENOTCONVERGED;   // do loop used up all iterations
  } else {
    return SCAMAC_EOK;
  }

}



ScamacErrorCode scamac_lanczos_ev_mat_cplx(const scamac_sparsemat_st *sm, double tol, double *ev1, double *ev2, double *eps1, double *eps2) {
  ScamacErrorCode err;

  if (!sm || !ev1 || !ev2 || !eps1 || !eps2) {
    return SCAMAC_ENULL;
  }
  if (sm->valtype != SCAMAC_VAL_COMPLEX) {
    return SCAMAC_EINVALID;
  }
  if (sm->nr != sm->nc) {
    return SCAMAC_EINVALID;
  }


  //maximal number of Lanczos iterations (not const .. because it's passed to LAPACK routine)
  int nl=200;

  int nv = sm->nc;

  // two vectors (long)
  double complex *x, *y;
  //tridiagonal matrix elements (short/few)
  double *d,*e,*z,*work,*dhp,*ehp;

  x=malloc(nv * sizeof *x);
  if (!x) {
    return SCAMAC_EMALLOCFAIL;
  }
  y=malloc(nv * sizeof *y);
  if (!y) {
    return SCAMAC_EMALLOCFAIL;
  }
  d=malloc(nl * sizeof *d);
  if (!d) {
    return SCAMAC_EMALLOCFAIL;
  }
  e=malloc(nl * sizeof *e);
  if (!e) {
    return SCAMAC_EMALLOCFAIL;
  }
  dhp=malloc(nl * sizeof *dhp);
  if (!dhp) {
    return SCAMAC_EMALLOCFAIL;
  }
  ehp=malloc(nl * sizeof *ehp);
  if (!ehp) {
    return SCAMAC_EMALLOCFAIL;
  }
  z=malloc(nl*nl * sizeof *z);
  if (!z) {
    return SCAMAC_EMALLOCFAIL;
  }
  work=malloc((2*nl-2) * sizeof *work);
  if (!work) {
    return SCAMAC_EMALLOCFAIL;
  }

  //random vector
  //dm_random_vector(nv,x);
  int idist = 2;
  int iseed[4] = {4,5,6,7};
  // here we go. Instead of zlarnv_, we use "dlarnv" on vector of length 2nv, to avoid problems with complex (C)LAPACK types
  // zlarnv_(&idist,iseed,&nv,x);
  int twonv = 2*nv;
  dlarnv_(&idist,iseed,&twonv,(double *) x);
  //normalize
  double dotx = cblas_dznrm2(nv,x,1);
  cblas_zdscal(nv,1.0/dotx,x,1);

  int il;
  int happybreak = 0;

  for (il=0; il<nl; il++) {
    if (il==0) {
      err=scamac_sparsemat_mvm_cplx(sm,x,y, 1.0,0.0,0.0);
      if (err) {
        return err|SCAMAC_EINTERNAL;
      }
      d[0]=cblas_ddot(2*nv, (double *) x,1, (double *) y,1);
      cblas_daxpy(2*nv,-d[0],(double *) x,1, (double *) y,1); // Y=Y-d(1)*X
      e[0]=cblas_dnrm2(2*nv, (double *) y,1); // dot_product(Y,Y))
      if ( fabs(e[0]) < DBL_EPSILON ) {
        happybreak=1;
      } else {
        cblas_dscal(2*nv, 1.0/e[0],(double *) y,1); // Y=Y/e(1)
      }
    } else if ( il%2 == 0 ) {
      err=scamac_sparsemat_mvm_cplx(sm,x,y, 1.0, -e[il-1],0.0);
      if (err) {
        return err|SCAMAC_EINTERNAL;
      }
      //dd=cblas_zdotc(nv,x,1,y,1);
      d[il]=cblas_ddot(2*nv,  (double *) x,1, (double *) y,1); // d(IL)=dot_product(X,Y)
      cblas_daxpy(2*nv,-d[il],(double *) x,1,(double *) y,1); // Y=Y-d(IL)*X
      e[il]=cblas_dnrm2(2*nv, (double *) y,1); // sqrt(ddot_dot_product(Y,Y))
      if ( fabs(e[il])< DBL_EPSILON ) {
        happybreak=1;
      } else {
        cblas_dscal(2*nv, 1.0/e[il], (double *) y,1);
      }
    } else {
      err=scamac_sparsemat_mvm_cplx(sm,y,x, 1.0, -e[il-1],0.0);
      if (err) {
        return err|SCAMAC_EINTERNAL;
      }
      d[il]=cblas_ddot(2*nv, (double *) y,1, (double *) x,1);
      cblas_daxpy(2*nv,-d[il],(double *) y,1,(double *) x,1);      //     X=X-d(IL)*Y
      e[il]=cblas_dnrm2(2*nv,(double *) x,1); //     e(IL)=sqrt(dot_product(X,X))
      if ( fabs(e[il]) < DBL_EPSILON ) {
        happybreak=1;
      } else {
        cblas_dscal(2*nv, 1.0/e[il],(double *) x,1);  // X=X/e(IL)
      }
    }

    // determine actual eigenvalue estimates = extremal Ritz values
    if (il==0) {
      *ev1=d[0];
      *ev2=d[0];
      *eps1=0.0;
      *eps2=0.0;
    } else {
      // ! --- needs Z be set to 1 initially? Check.
      //dm_identity(nl,1.0,z);
      cblas_dcopy(il+1,d,1,dhp,1);
      cblas_dcopy(il+1,e,1,ehp,1);
      char compz='I';
      int ilplus1 = il+1;
      int INFO;
      dsteqr_(&compz, &ilplus1,dhp,ehp,z,&nl,work,&INFO);
      if (INFO) {
        return SCAMAC_EFAIL|SCAMAC_EINTERNAL;
      }
      *ev1=dhp[0];
      *ev2=dhp[il];
      //error bounds from vector   (note: e(IL) is not touched in DSTEQR)
      *eps1=fabs(z[il+nl*0])*fabs(e[il]);
      *eps2=fabs(z[il+nl*il])*fabs(e[il]);
      if (*ev2-*ev1>sqrt(DBL_EPSILON)) { // accuracy achieved?
        // relative error
        if ( (*eps1/(*ev2-*ev1)<tol) && (*eps2/(*ev2-*ev1)<tol) ) {
          break;
        }
      } else {
        //absolute error
        if ( (*eps1<tol) && (*eps2<tol) ) {
          break;
        }
      }
    }

    if (happybreak) break;

  }

  free(x);
  free(y);
  free(d);
  free(e);
  free(dhp);
  free(ehp);
  free(z);
  free(work);


  if (il >= nl) {
    return SCAMAC_ENOTCONVERGED;   // do loop used up all iterations
  } else {
    return SCAMAC_EOK;
  }

}



ScamacErrorCode scamac_lanczos_ev_mat(const scamac_sparsemat_st *sm, double tol, double *ev1, double *ev2, double *eps1, double *eps2) {
  ScamacErrorCode err;
  if (sm->valtype == SCAMAC_VAL_REAL) {
    err = scamac_lanczos_ev_mat_real(sm, tol, ev1, ev2, eps1, eps2);
  } else if (sm->valtype == SCAMAC_VAL_COMPLEX) {
    err = scamac_lanczos_ev_mat_cplx(sm, tol, ev1, ev2, eps1, eps2);
  } else {
    err = SCAMAC_ECORRUPTED;
  }
  return err;
}

ScamacErrorCode scamac_lanczos_ev(const ScamacGenerator *gen, double tol, double *ev1, double *ev2, double *eps1, double *eps2) {
  ScamacErrorCode err;
  scamac_sparsemat_st * sm;
  err = scamac_sparsemat_from_generator(gen, &sm);
  if (err) {
    return err;
  }
  err = scamac_lanczos_ev_mat(sm, tol, ev1, ev2, eps1, eps2);
  scamac_sparsemat_free(sm);
  return err;
}

