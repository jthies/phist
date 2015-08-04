#include "phist_config.h"

#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#ifdef PHIST_HAVE_LIKWID
#include "likwid.h"
#endif

#include "phist_macros.h"
#include "phist_enums.h"
#include "phist_kernels.h"
#include "phist_operator.h"
#include "phist_svrr.h"

#include "phist_ScalarTraits.hpp"

#include "phist_gen_d.h"
#include "phist_driver_utils.h"

typedef phist::ScalarTraits<ST> st;
typedef phist::ScalarTraits<MT> mt;

// workaround for #131 (put_func fails with Intel 15.0up02)
void SUBR(mvec_put_func_substitute)(TYPE(mvec_ptr) V, int (*func)(gidx_t,lidx_t,void*),int* ierr)
{
  _ST_ *val;
  lidx_t nloc, lda, ncols;
  gidx_t ilower, iupper;
  const_map_ptr_t map;
  PHIST_CHK_IERR(SUBR(mvec_get_map)(V,&map,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&ncols,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvec_my_length)(V,&nloc,ierr),*ierr);
  PHIST_CHK_IERR(phist_map_get_ilower(map,&ilower,ierr),*ierr);
  PHIST_CHK_IERR(phist_map_get_iupper(map,&iupper,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvec_from_device)(V,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvec_extract_view)(V,&val,&lda,ierr),*ierr);
  for (lidx_t i=0; i<nloc; i++)
  {
    gidx_t ii=i+ilower;
    for (lidx_t jj=0; jj<ncols; jj++)
    {
#ifdef PHIST_MVECS_ROW_MAJOR
      func(ii,jj,&val[i*lda+jj]);
#else
      func(ii,jj,&val[jj*lda+i]);
#endif
    }
  }
  PHIST_CHK_IERR(SUBR(mvec_to_device)(V,ierr),*ierr);
}

// prototypes and mvec initialization functions

// function that performs numerical tests
static int run_tests1(Dmvec_ptr_t X, bool high_prec);

extern "C" {

// Hilbert matrix generating function
int hilbert(ghost_gidx_t i, ghost_lidx_t j, void* val)
{
  ST* dval=(ST*)val;
  *dval = st::one()/((ST)i+(ST)j+st::one());
  std::cout<<"TROET "<<i<<" "<<j<<" "<<*dval<<std::endl;
  return 0;
}

// synthetic matrix from Yamasaki et al
int synth1(ghost_gidx_t i, ghost_lidx_t j, void* val)
{
  ST* dval=(ST*)val;
  if (i==0)
  {
    *dval = st::one();
  }
  else if (i==(gidx_t)(j+1))
  {
    ST e=st::eps();
    *dval = st::rand()*e*e*e;
  }
  else
  {
    *dval=st::zero();
  }
  return 0;
}

}//extern "C"

//! benchmark the function svrr, given k, n, m
//! do k times: randomize an n x m mvec and orthogonalize
//! it using svrr.
int main(int argc, char** argv)
{
  int ierr;

  comm_ptr_t comm = NULL;
  const_map_ptr_t cmap = NULL;
  map_ptr_t map = NULL;
  mvec_ptr_t X = NULL, Q=NULL;
  sparseMat_ptr_t A = NULL;
  
  PHIST_ICHK_IERR(phist_kernels_init(&argc,&argv,&ierr),ierr);

PHIST_MAIN_TASK_BEGIN

  gidx_t n;
  int m;
  
  PHIST_ICHK_IERR(phist_comm_create(&comm,&ierr),ierr);
  PHIST_ICHK_IERR(SUBR(create_matrix)(&A,comm,"fem_lap32.mm",&ierr),ierr);
  PHIST_ICHK_IERR(SUBR(sparseMat_get_domain_map)(A,&cmap,&ierr),ierr);

  // TEST CASE 1: K_20(A,1)
  for (m=20;m<=30;m+=10)
  {
    PHIST_ICHK_IERR(SUBR(mvec_create)(&X,cmap,m,&ierr),ierr);
  
    TYPE(mvec_ptr) vi=NULL,vj=NULL;
    PHIST_ICHK_IERR(SUBR(mvec_put_value)(X,st::one(),&ierr),ierr);
  
    PHIST_ICHK_IERR(SUBR(mvec_view_block)(X,&vi,0,0,&ierr),ierr);
    MT dum;
    PHIST_ICHK_IERR(SUBR(mvec_normalize)(vi,&dum,&ierr),ierr);
  
    for (int i=1; i<m; i++)
    {
      PHIST_ICHK_IERR(SUBR(mvec_view_block)(X,&vj,i,i,&ierr),ierr);
      PHIST_ICHK_IERR(SUBR(sparseMat_times_mvec)(st::one(),A,vi,st::zero(),vj,&ierr),ierr);
      PHIST_ICHK_IERR(SUBR(mvec_normalize)(vj,&dum,&ierr),ierr);
      vi=vj;
      vj=NULL;
    }

    PHIST_SOUT(PHIST_INFO,"PROBLEM: Lapl%d matrix, standard\n",m);
    run_tests1(X,false);

    PHIST_SOUT(PHIST_INFO,"PROBLEM: Lapl%d matrix, high_prec\n",m);
    run_tests1(X,true);
    
    PHIST_ICHK_IERR(SUBR(mvec_delete)(X,&ierr),ierr);
  }

  n=20;
  m=20;
  PHIST_ICHK_IERR(phist_map_create(&map,comm,n,&ierr),ierr);
  PHIST_ICHK_IERR(SUBR(mvec_create)(&X,map,m,&ierr),ierr);
  //PHIST_ICHK_IERR(SUBR(mvec_put_func)(X,&hilbert,&ierr),ierr);
  PHIST_ICHK_IERR(SUBR(mvec_put_func_substitute)(X,&hilbert,&ierr),ierr);
    PHIST_SOUT(PHIST_INFO,"PROBLEM: Hilbert%dx%d matrix, standard\n",(int)n,(int)m);
    run_tests1(X,false);
    PHIST_SOUT(PHIST_INFO,"PROBLEM: Hilbert%dx%d matrix, high_prec\n",(int)n,(int)m);
    run_tests1(X,true);
    PHIST_ICHK_IERR(SUBR(mvec_delete)(X,&ierr),ierr);
    PHIST_ICHK_IERR(phist_map_delete(map,&ierr),ierr);

  n=21;
  m=20;
  PHIST_ICHK_IERR(phist_map_create(&map,comm,n,&ierr),ierr);
  PHIST_ICHK_IERR(SUBR(mvec_create)(&X,map,m,&ierr),ierr);
  PHIST_ICHK_IERR(SUBR(mvec_put_func_substitute)(X,&synth1,&ierr),ierr);
  //PHIST_ICHK_IERR(SUBR(mvec_put_func)(X,&synth1,&ierr),ierr);
    PHIST_SOUT(PHIST_INFO,"PROBLEM: Synth%dx%d matrix, standard\n",(int)n,(int)m);
    run_tests1(X,false);
    PHIST_SOUT(PHIST_INFO,"PROBLEM: Synth%dx%d matrix, high_prec\n",(int)n,(int)m);
    run_tests1(X,true);
    PHIST_ICHK_IERR(SUBR(mvec_delete)(X,&ierr),ierr);
    PHIST_ICHK_IERR(phist_map_delete(map,&ierr),ierr);
  
PHIST_MAIN_TASK_END

  return ierr;
}







int run_tests1(Dmvec_ptr_t X, bool high_prec)
{
#ifndef PHIST_HIGH_PRECISION_KERNELS
  if (high_prec) 
  {
    PHIST_SOUT(PHIST_ERROR,"high precision kernels not available\n");
    return -99;
  }
#endif
  int ierr=0;
  const_comm_ptr_t comm=NULL;
  const_map_ptr_t map=NULL;
  mvec_ptr_t Q=NULL;
  sdMat_ptr_t R = NULL, QtQ=NULL;
  
  lidx_t n;
  int m;

  PHIST_ICHK_IERR(SUBR(mvec_get_map)(X,&map,&ierr),ierr);
  PHIST_ICHK_IERR(phist_map_get_comm(map,&comm,&ierr),ierr);
  PHIST_ICHK_IERR(SUBR(mvec_my_length)(X,&n,&ierr),ierr);
  PHIST_ICHK_IERR(SUBR(mvec_num_vectors)(X,&m,&ierr),ierr);
  PHIST_ICHK_IERR(SUBR(mvec_create)(&Q,map,m,&ierr),ierr);
  PHIST_ICHK_IERR(SUBR(sdMat_create)(&R,m,m,comm,&ierr),ierr);
  PHIST_ICHK_IERR(SUBR(sdMat_create)(&QtQ,m,m,comm,&ierr),ierr);
  
  PHIST_ICHK_IERR(SUBR(mvec_add_mvec)(st::one(),X,st::zero(),Q,&ierr),ierr);
  
  if (!high_prec && n<200)
  {
    PHIST_SOUT(PHIST_INFO,"input matrix:\n");
    PHIST_ICHK_IERR(SUBR(mvec_print)(X,&ierr),ierr);
  }
  
  for (int it=0; it<5; it++)
  {

    if( high_prec )
    {
      ierr = PHIST_ROBUST_REDUCTIONS;
    }
    PHIST_ICHK_NEG_IERR(SUBR(svrr)(Q,R,&ierr),ierr);
    if (ierr) PHIST_SOUT(PHIST_INFO,"found rank(V)=%d\n",m-ierr);

    // check the result
    if( high_prec )
    {
      ierr = PHIST_ROBUST_REDUCTIONS;
    }
    PHIST_ICHK_IERR(SUBR(mvecT_times_mvec)(st::one(),Q,Q,st::zero(),QtQ,&ierr),ierr);
    PHIST_ICHK_IERR(SUBR(sdMat_from_device)(QtQ,&ierr),ierr);
    ST* q=NULL;
    lidx_t ldq;
    PHIST_ICHK_IERR(SUBR(sdMat_extract_view)(QtQ,&q,&ldq,&ierr),ierr);
    MT err=mt::zero();
    for (int i=0;i<m;i++)
      for (int j=0;j<m;j++)
      {
        MT val=(i==j)? mt::one() : mt::zero();
        err+= (q[i*ldq+j] - val)*(q[i*ldq+j]-val);
      }

    PHIST_SOUT(PHIST_INFO,"STEP %d, ||Q^TQ-I||_2: %8.4e\n",it+1,err);
    PHIST_SOUT(PHIST_INFO,"Example for Q^T Q:\n");
    PHIST_ICHK_IERR(SUBR(sdMat_print)(QtQ,&ierr),ierr);

  }

  PHIST_ICHK_IERR(SUBR(mvec_delete)(Q,&ierr),ierr);
  PHIST_ICHK_IERR(SUBR(sdMat_delete)(R,&ierr),ierr);
  PHIST_ICHK_IERR(SUBR(sdMat_delete)(QtQ,&ierr),ierr);

  return ierr;
}

