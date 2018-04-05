/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
namespace phist
{

  namespace testing
  {

void SUBR(read_mat)(const char* filebase,phist_const_comm_ptr comm, 
        int nglob, int mglob, TYPE(sparseMat_ptr) *ptr, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  int iflag_in=*iflag;
  std::stringstream dimstr;
  dimstr<<nglob;
  if (mglob!=nglob) dimstr<<"x"<<mglob;
  iflag_in |= PHIST_OUTLEV>=PHIST_DEBUG ? 0 : PHIST_SPARSEMAT_QUIET;
  *ptr = NULL;
  char tpc = ::phist::ScalarTraits< _ST_ >::type_char();
  char mmfile[256],hbfile[256],binfile[256],crsfile[256];
  sprintf(mmfile,"%c%s%d.mm",tpc,filebase,nglob);
  if (::phist::ScalarTraits< _ST_ >::is_complex())
  {
    sprintf(hbfile,"%c%s%s.cua",tpc,filebase,dimstr.str().c_str());
  } 
  else
  {
    sprintf(hbfile,"%c%s%s.rua",tpc,filebase,dimstr.str().c_str());
  }

  sprintf(binfile,"%c%s%s.bin",tpc,filebase,dimstr.str().c_str());
  sprintf(crsfile,"%c%s%s.crs",tpc,filebase,dimstr.str().c_str());
  
  PHIST_SOUT(PHIST_DEBUG, "Looking for matrix \'%s\'..\n", filebase);

  PHIST_SOUT(PHIST_DEBUG, "... try \'%s\'\n", binfile);
  *iflag = iflag_in;
  SUBR(sparseMat_read_bin)(ptr,comm,binfile,iflag);
  if (*iflag==PHIST_SUCCESS) return;

  // try same format, different extension
  PHIST_SOUT(PHIST_DEBUG, "... try \'%s\'\n", crsfile);
  *iflag = PHIST_OUTLEV>=PHIST_DEBUG ? 0 : PHIST_SPARSEMAT_QUIET;
  SUBR(sparseMat_read_bin)(ptr,comm,crsfile,iflag);
  if (*iflag==PHIST_SUCCESS) return;
  
  PHIST_SOUT(PHIST_DEBUG, "... try \'%s\'\n", mmfile);
  *iflag = iflag_in;
  SUBR(sparseMat_read_mm)(ptr,comm,mmfile,iflag);
  if (*iflag==PHIST_SUCCESS) return;

  PHIST_SOUT(PHIST_DEBUG, "... try \'%s\'\n", hbfile);
  *iflag = iflag_in;
  SUBR(sparseMat_read_hb)(ptr,comm,hbfile,iflag);
}//read_mat

// some functions for initializing matrices and mvecs
int PHIST_TG_PREFIX(idfunc)(ghost_gidx row, ghost_lidx *len, ghost_gidx* cols, void* vval, void *arg)
{
  static phist_gidx gnrows=-1;
  static phist_gidx gncols=-1;
  if (row==-1)
  {
    gnrows=cols[0];
    gncols=cols[1];
    return 0;
  }
  else if (row==-2)
  {
    gnrows=-1;
    gncols=-1;
  }
  if (gnrows>0 && row>=gnrows) return -1;
  if (row<0) return -1;
  if (gncols>0 && row>=gncols)
  {
    *len=0;
    return 0;
  }
  *len=1;
  _ST_* val = (_ST_*)vval;
  _ST_* diagp= (_ST_*)arg;
  _ST_ diag=(diagp==NULL)? (_ST_)1.0: (*diagp);
  val[0]=diag;
  cols[0]=row;
  return 0;
}

  // some row function that uses a workspace. Will only create the identity matrix if
  // initialized using idfunc_init_workspace is called. Requires arg to be a struct with
  // a pointer to idfunc_with_workspace_arg.
  int PHIST_TG_PREFIX(idfunc_with_workspace)(ghost_gidx row, ghost_lidx *len, ghost_gidx* cols, void* vval, void *v_arg)
  {
    PHIST_TG_PREFIX(idfunc_workspace)* wsp = 
        (PHIST_TG_PREFIX(idfunc_workspace)*)v_arg;
    PHIST_TG_PREFIX(idfunc_with_workspace_arg)* arg = wsp->arg;
    _ST_* val = (_ST_*)vval;
    if (row<0 || row>=arg->gnrows || row>=arg->gncols)
    {
      return -1; // array out of bounds
    }
    else if (wsp->data == NULL)
    {
      return -2; // not initialized
    }
    *len=1;
    cols[0] = row;
    val[0] = wsp->data[0] * arg->scale;
    return 0;
  }

  // 'constructor'
  int PHIST_TG_PREFIX(idfunc_init_workspace)(void *v_arg, void** work)
  {
    PHIST_TG_PREFIX(idfunc_with_workspace_arg)* arg =
        (PHIST_TG_PREFIX(idfunc_with_workspace_arg)*)v_arg;
    if (!arg) return -1;
    if (*work!=nullptr)
    {
      PHIST_TG_PREFIX(idfunc_workspace)* wsp =(PHIST_TG_PREFIX(idfunc_workspace)*)(*work);
      if (wsp->arg == arg)
      {
        // second call (destroy)
        delete [] wsp->data;
        delete wsp;
        *work=NULL;
      }
      else
      {
        //probably first call, but *work not set to NULL
        return PHIST_INVALID_INPUT;
      }
    }
    else
    {
      // first call (create)
      PHIST_TG_PREFIX(idfunc_workspace)* my_work=new PHIST_TG_PREFIX(idfunc_workspace);
      my_work->data=new _ST_[1];
      my_work->data[0]=_ST_(1.0);
      my_work->arg=arg;
      *work=(void*)my_work;
    }
    return 0;
  }



/*
int PHIST_TG_PREFIX(some_rowFunc)(ghost_gidx row, ghost_lidx *len, ghost_gidx* cols, void* vval, void *arg)
{
#include "phist_std_typedefs.hpp"
  _ST_* val = (_ST_*)vval;

  *len=5;
  for (int i=0; i<*len; i++)
  {
    cols[i]=(ghost_gidx)(((row+i-2)*3)%_N_);
    if (cols[i]<0) cols[i]+=_N_;
    val[i]=(ST)(i+1)/(ST)(row+1) + st::cmplx_I()*(ST)(row-cols[i]);
  }
  return 0;
}
*/



// a big macro to define the content of our (tridiagonal and block diagonal matrices quickly further down.
// The matrix will be tridiag([ALOWER, ADIAG, AUPPER],-1:1,n,n) if BLOCKDIAG==false, and
// kron(I,[ADIAG,AUPPER; ALOWER, ADIAG]) if BLOCKDIAG==true (in MATLAB notation). The macro also takes care
// of the initialization of the problem size so that this code doesn't have to be replicated.
#ifndef TRIDIAG
#define TRIDIAG(ADIAG,AUPPER,ALOWER,BLOCKDIAG) \
  static ghost_gidx gnrows=-1;\
  _ST_ *vals=(_ST_*)vval;\
  \
  if (vals) vals[0]=(ADIAG);\
  if (cols && row>=0) cols[0]=row;\
\
  if (row==-1)\
  {\
    gnrows=cols[0];\
    return 0;\
  }\
  else if (row==-2)\
  {\
    gnrows=-2;\
  }\
  else if (gnrows<0)\
  {\
    PHIST_SOUT(PHIST_ERROR,"%s not correctly initialized, call with row=-1 and cols[0]=gnrows first!",__FUNCTION__);\
    return -1;\
  }\
  else if (row==0 || ((BLOCKDIAG) && (row%2==0)))\
  {\
    *len=2;\
    cols[1]=row+1;\
    vals[1]=(AUPPER);\
  }\
  else if (row==gnrows-1 || ((BLOCKDIAG) && (row%2==1)))\
  {\
    *len=2;\
    cols[1]=row-1;\
    vals[1]=(ALOWER);\
  }\
  else\
  {\
    *len=3;\
    cols[1]=row-1;\
    cols[2]=row+1;\
    vals[1]=(ALOWER);\
    vals[2]=(AUPPER);\
}
#endif

int PHIST_TG_PREFIX(hpd_tridiag)(ghost_gidx row, ghost_lidx *len, ghost_gidx* cols, void* vval, void *arg)
{
#include "phist_std_typedefs.hpp"
#ifdef IS_COMPLEX
    TRIDIAG(st::one(),ST(-0.5)*st::cmplx_I(),ST(+0.5)*st::cmplx_I(),false);
#else
    TRIDIAG(st::one(),ST(-0.5),ST(-0.5),false);
#endif
  return 0;
}

// 1D laplacian, tridiag([-1 2 -1])
int PHIST_TG_PREFIX(lapl_tridiag)(ghost_gidx row, ghost_lidx *len, ghost_gidx* cols, void* vval, void *arg)
{
#include "phist_std_typedefs.hpp"
  TRIDIAG((_ST_)2.0,ST(-1.0),ST(-1.0),false);
  return 0;
}

//! creates a simple tridiagonal non-Hermitian but positive definite matrix. For usage info, see hpd_tridiag.
int PHIST_TG_PREFIX(nhpd_tridiag)(ghost_gidx row, ghost_lidx *len, ghost_gidx* cols, void* vval, void *arg)
{
#include "phist_std_typedefs.hpp"
  _ST_ a=(_ST_)2.0, b=(_ST_)-0.9, c=(_ST_)-1.1;
  TRIDIAG(a,b,c,false);
  return 0; 
}

//! creates a simple tridiagonal non-Hermitian and indefinite matrix. For usage info, see hpd_tridiag.
int PHIST_TG_PREFIX(hid_tridiag)(ghost_gidx row, ghost_lidx *len, ghost_gidx* cols, void* vval, void *arg)
{
#include "phist_std_typedefs.hpp"
  // construct laplacian but shift it so that a few eigenvalues are negative

  static double L=-1;
  if (row==-1) L=(cols[0]+1); // L=gnrows+1
  int k=10; // number of eigenvalues to shift over the axis
  static const double pi=4.0*std::atan(1.0);
  double pi_div_L=pi/L;
  double ev_k  =(k    *pi_div_L)*(k    *pi_div_L);
  double ev_kp1=((k+1)*pi_div_L)*((k+1)*pi_div_L);
  double shift=(0.5*(ev_k+ev_kp1));

  TRIDIAG(ST(2.0-shift),ST(-1.0),ST(-1.0),false);
  return 0;
}

//! creates an approximate inverse of hpd_tridiag (the inverse of the 2x2 block diagonal approximation of A)
int PHIST_TG_PREFIX(hpd_tridiag_ainv)(ghost_gidx row, ghost_lidx *len, ghost_gidx* cols, void* vval, void *arg)
{
#include "phist_std_typedefs.hpp"
#ifdef IS_COMPLEX
  TRIDIAG((_ST_)0.8,ST(0.4)*st::cmplx_I(),ST(0.4)*st::cmplx_I(),true);
#else
  TRIDIAG(ST(4./3.),ST(2./3.),ST(2./3.),true);
#endif
  return 0;
}

//! creates an approximate inverse of lapl_tridiag (the inverse of the 2x2 block diagonal approximation of A)
int PHIST_TG_PREFIX(lapl_tridiag_ainv)(ghost_gidx row, ghost_lidx *len, ghost_gidx* cols, void* vval, void *arg)
{
#include "phist_std_typedefs.hpp"
  TRIDIAG(ST(2./3.),ST(1./3.),ST(1./3.),true);
  return 0;
}

//! creates an approximate inverse of nhpd_tridiag (the inverse of the 2x2 block diagonal approximation of A)
int PHIST_TG_PREFIX(nhpd_tridiag_ainv)(ghost_gidx row, ghost_lidx *len, ghost_gidx* cols, void* vval, void *arg)
{
#include "phist_std_typedefs.hpp"

  _ST_ a=(_ST_)2.0, b=(_ST_)-0.9, c=(_ST_)-1.1;
  _ST_ s=st::one()/(a*a-b*c);
  TRIDIAG(a*s,-c*s,-b*s,true);
  return 0;
}

//! creates an approximate inverse of nhid_tridiag (the inverse of the 2x2 block diagonal approximation of A)
int PHIST_TG_PREFIX(hid_tridiag_ainv)(ghost_gidx row, ghost_lidx *len, ghost_gidx* cols, void* vval, void *arg)
{
#include "phist_std_typedefs.hpp"
  static double L=-1;
  if (row==-1) L=(cols[0]+1); // L=gnrows+1
  int k=10; // number of eigenvalues to shift over the axis
  static const double pi=4.0*std::atan(1.0);
  double pi_div_L=pi/L;
  double ev_k  =(k    *pi_div_L)*(k    *pi_div_L);
  double ev_kp1=((k+1)*pi_div_L)*((k+1)*pi_div_L);
  double shift=(0.5*(ev_k+ev_kp1));

  _ST_ a=ST(2.0-shift), b=(_ST_)-1.0, c=(_ST_)-1.0;
  _ST_ s=st::one()/(a*a-b*c);
  TRIDIAG(a*s,-c*s,-b*s,true);
  return 0;
}
//! creates an approximate inverse of nhid_tridiag (the inverse of the 2x2 block diagonal approximation of A)
int PHIST_TG_PREFIX(nhid_tridiag_ainv)(ghost_gidx row, ghost_lidx *len, ghost_gidx* cols, void* vval, void *arg)
{
  return -99; // not implemented
}
  // defines a matrix with only a constant subdiagonal and an entry (1,N)
  // that defines a periodic "right shift" of vector elements.
  int PHIST_TG_PREFIX(right_shift_perio)(ghost_gidx row, ghost_lidx *len, ghost_gidx* cols, void* vval, void *arg)
  {
#include "phist_std_typedefs.hpp"
    static ghost_gidx gnrows=-1;
    _ST_ *vals=(_ST_*)vval;
  
    if (row==-1)
    {
      gnrows=cols[0];
      return 0;
    }
    else if (row==-2)
    {
      gnrows=-1;
    }
    else if (gnrows<0)
    {
      PHIST_SOUT(PHIST_ERROR,"%s not correctly initialized, call with row=-1 and cols[0]=gnrows first!",__FUNCTION__);
      return -1;
    }
    else
    {
      *len=1;
      cols[0]=row>0?row-1:gnrows-1;
      vals[0]=st::one();
    }
    return 0;
  }

  // defines a matrix with only a constant subdiagonal that defines a non-periodic "right shift" of vector elements.
  int PHIST_TG_PREFIX(right_shift)(ghost_gidx row, ghost_lidx *len, ghost_gidx* cols, void* vval, void *arg)
  {
#include "phist_std_typedefs.hpp"
    _ST_ *vals=(_ST_*)vval;
    
    if (row<0) return 0;
  
    *len=0;
    if (row>0)
    {
      *len=1;
      cols[0]=row-1;
      vals[0]=st::one();
    }
    return 0;
  }

  // defines a matrix with only a constant superdiagonal and an entry (N,1)
  // that defines a periodic "left shift" of vector elements.
  int PHIST_TG_PREFIX(left_shift_perio)(ghost_gidx row, ghost_lidx *len, ghost_gidx* cols, void* vval, void *arg)
  {
#include "phist_std_typedefs.hpp"
    static ghost_gidx gnrows=-1;
    _ST_ *vals=(_ST_*)vval;
  
    if (row==-1)
    {
      gnrows=cols[0];
      return 0;
    }
    else if (row==-2)
    {
      gnrows=-1;
      return 0;
    }
    else if (gnrows<0)
    {
      PHIST_SOUT(PHIST_ERROR,"%s not correctly initialized, call with row=-1 and cols[0]=gnrows first!",__FUNCTION__);
      return -1;
    }
    else
    {
      *len=1;
      cols[0]=row<gnrows-1?row+1:0;
      vals[0]=st::one();
    }
    return 0;
  }


  // defines a matrix with only a constant superdiagonal that defines a non-periodic "left shift" of vector elements.
  int PHIST_TG_PREFIX(left_shift)(ghost_gidx row, ghost_lidx *len, ghost_gidx* cols, void* vval, void *arg)
  {
#include "phist_std_typedefs.hpp"
    static ghost_gidx gnrows=-1;
    _ST_ *vals=(_ST_*)vval;
    if (row==-1)
    {
      gnrows=cols[0];
      return 0;
    }
    else if (row==-2)
    {
      gnrows=-1;
      return 0;
    }
    else if (gnrows<0)
    {
      PHIST_SOUT(PHIST_ERROR,"%s not correctly initialized, call with row=-1 and cols[0]=gnrows first!",__FUNCTION__);
      return -1;
    }
    else
    {
      *len=0;
      if (row<gnrows-1)
      {
        *len=1;
        cols[0]=row+1;
        vals[0]=st::one();
      }
    }
    return 0;
  }




  int PHIST_TG_PREFIX(mvec123func)(ghost_gidx i, ghost_lidx j, void* val, void* last_arg)
  {
#include "phist_std_typedefs.hpp"
    _ST_* v= (_ST_*)val;
    int *int_arg=(int*)last_arg;
    int N  = int_arg[0];
    int NV = int_arg[1];
    *v = ST(i+1 + N*j);
    return 0;
  }

  int PHIST_TG_PREFIX(mvec321func)(ghost_gidx i, ghost_lidx j, void* val, void* last_arg)
  {
#include "phist_std_typedefs.hpp"
    _ST_* v= (_ST_*)val;
    int *int_arg=(int*)last_arg;
    int N  = int_arg[0];
    int NV = int_arg[1];
    *v = ST((N-i) + N*(NV-(j+1)));
    return 0;
  }


  }//testing
}//phist
