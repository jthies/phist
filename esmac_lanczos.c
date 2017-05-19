#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <gsl/gsl_cblas.h>

#include "esmac_lanczos.h"



// LAPACK
extern void dsteqr_(char *compz, int *n, double *d, double *e, double *z, int *ldz, double *work, int *info);
extern void dlarnv_(int *idist, int *iseed, int *n, double *x);

int esmac_lanczos_ev_mat(const esmac_sparsemat_t *sm, double tol, double *ev1, double *ev2, double *eps1, double *eps2) {
  
  //maximal number of Lanczos iterations (not const .. because it's passed to LAPACK routine)
  int nl=200;
  
  if (sm->nr != sm->nc) {
    printf("%s: need square matrix\n",__func__);
    exit(EXIT_FAILURE);
  }
  
  int nv = sm->nc;

  // two vectors (long)
  double *x, *y;
  //tridiagonal matrix elements (short/few)
  double *d,*e,*z,*work,*dhp,*ehp;

  x=malloc(nv * sizeof *x);
  y=malloc(nv * sizeof *y);
  d=malloc(nl * sizeof *d);
  e=malloc(nl * sizeof *e);
  dhp=malloc(nl * sizeof *dhp);
  ehp=malloc(nl * sizeof *ehp);
  z=malloc(nl*nl * sizeof *z);
  work=malloc((2*nl-2) * sizeof *work);

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
  
  for (il=0;il<nl;il++) {
    if (il==0) {
      esmac_sparsemat_mvm(sm,x,y, 1.0,0.0,0.0);
      d[0]=cblas_ddot(nv,x,1,y,1);     //  d[0]=dot_product(nv,x,y);
      cblas_daxpy(nv,-d[0],x,1,y,1); // Y=Y-d(1)*X
      e[0]=cblas_dnrm2(nv,y,1); // dot_product(Y,Y))
      if ( fabs(e[0]) < DBL_EPSILON ) {
	happybreak=1;
      } else {
	cblas_dscal(nv,1.0/e[0],y,1);	// Y=Y/e(1)
      }
    } else {
      if ( il%2 == 0 ) {
	esmac_sparsemat_mvm(sm,x,y, 1.0, -e[il-1],0.0);
	d[il] = cblas_ddot(nv,x,1,y,1); // d(IL)=dot_product(X,Y)
	cblas_daxpy(nv,-d[il],x,1,y,1); // Y=Y-d(IL)*X
	e[il]=cblas_dnrm2(nv,y,1); // sqrt(ddot_dot_product(Y,Y))
	if ( fabs(e[il])< DBL_EPSILON ) {
	  happybreak=1;
	} else {
	  cblas_dscal(nv,1.0/e[il],y,1); 
	}
      } else {
	esmac_sparsemat_mvm(sm,y,x, 1.0, -e[il-1],0.0);
	d[il] = cblas_ddot(nv,y,1,x,1);      //     d(IL)=dot_product(Y,X)
  	cblas_daxpy(nv,-d[il],y,1,x,1);      //     X=X-d(IL)*Y
	e[il]=cblas_dnrm2(nv,x,1); //     e(IL)=sqrt(dot_product(X,X))
	if ( fabs(e[il]) < DBL_EPSILON ) {
	  happybreak=1;
	} else {
	  cblas_dscal(nv,1.0/e[il],x,1);  // X=X/e(IL)
	}
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
      //	 printf("LAPACK\n");
      dsteqr_(&compz, &ilplus1 ,dhp,ehp,z,&nl,work,&INFO);
      //	 printf("done\n");
      if (INFO != 0) {
	printf("%s: Error in DSTEQR : %d",__func__,INFO);
	exit(EXIT_FAILURE);
      }
      *ev1=dhp[0];
      *ev2=dhp[il];
      //error bounds from vector   (note: e(IL) is not touched in DSTEQR)
      *eps1=fabs(z[il+nl*0])*fabs(e[il]);
      *eps2=fabs(z[il+nl*il])*fabs(e[il]);
      if (*ev2-*ev1>sqrt(DBL_EPSILON)) { // accuracy achieved?
	// relative error
	if ( (*eps1/(*ev2-*ev1)<tol) && (*eps2/(*ev2-*ev1)<tol) ) {break;}
      } else {
	//absolute error
	if ( (*eps1<tol) && (*eps2<tol) ) {break;}
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
    return nl;   // do loop used up all iterations
  } else {
    return -il;
  }
  
}


int esmac_lanczos_ev(const esmac_generator_t *gen, double tol, double *ev1, double *ev2, double *eps1, double *eps2) {
  esmac_sparsemat_t *sm = esmac_sparsemat_from_generator(gen);
  int info;
  info = esmac_lanczos_ev_mat(sm, tol, ev1, ev2, eps1, eps2);
  esmac_sparsemat_free(sm);
  return info;
}

