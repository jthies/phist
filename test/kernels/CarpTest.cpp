/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
#include "phist_config.h"

#include "phist_tools.h"
#include "kernels/phist_kernels.h"
#include "matfuncs/matfuncs.h"
#include "gtest/phist_gtest.h"


#include "kernels/phist_kernel_flags.h"
#include "KernelTestWithSparseMat.h"
#include "KernelTestWithVectors.h"
#include "KernelTestWithSdMats.h"

#include "phist_normF.h"

#include <complex>

#include "phist_gen_d.h"
#include "krylov/phist_carp_cg_kernels_decl.hpp"
#include "phist_gen_clean.h"

using namespace ::testing;

/* tests involving carp sweeps are only executed if the carp kernel is implemented,
if the setup routine returns -99, further tests will not be run
Other tests compare "RC" variants with complex arithmetic implementations, these are
only run if the kernel lib supports the Z type.
*/
#ifdef PHIST_HAVE_CMPLX

double MvecsEqualZD(phist_Zmvec* zvec, phist_Dmvec* dvec_r, phist_Dmvec* dvec_i, double relTo=0.0)
{

  int iflag=0;
  phist_Dmvec_from_device(dvec_r,&iflag);
  if (iflag!=PHIST_SUCCESS) return -1.0e11;
  phist_Dmvec_from_device(dvec_i,&iflag);
  if (iflag!=PHIST_SUCCESS) return -1.0e12;
  phist_Zmvec_from_device(zvec,&iflag);
  if (iflag!=PHIST_SUCCESS) return -1.0e13;

  double *dval_r, *dval_i;
  phist_d_complex *zval;
  phist_lidx lda_z, lda_r, lda_i;
  int nvec_r, nvec_i, nvec_z;
  phist_lidx nloc_r, nloc_i, nloc_z;
  
  phist_Dmvec_extract_view(dvec_r,&dval_r,&lda_r,&iflag);
  if (iflag!=PHIST_SUCCESS) return -1.0e21;

  phist_Dmvec_extract_view(dvec_i,&dval_i,&lda_i,&iflag);
  if (iflag!=PHIST_SUCCESS) return -1.0e22;

  phist_Zmvec_extract_view(zvec,&zval,&lda_z,&iflag);
  if (iflag!=PHIST_SUCCESS) return -1.0e23;

  phist_Dmvec_num_vectors(dvec_r,&nvec_r,&iflag);
  if (iflag!=PHIST_SUCCESS) return -1.0e24;

  phist_Dmvec_num_vectors(dvec_i,&nvec_i,&iflag);
  if (iflag!=PHIST_SUCCESS) return -1.0e25;

  phist_Zmvec_num_vectors(zvec,&nvec_z,&iflag);
  if (iflag!=PHIST_SUCCESS) return -1.0e26;

  phist_Dmvec_my_length(dvec_r,&nloc_r,&iflag);
  if (iflag!=PHIST_SUCCESS) return -1.0e27;

  phist_Dmvec_my_length(dvec_i,&nloc_i,&iflag);
  if (iflag!=PHIST_SUCCESS) return -1.0e28;

  phist_Zmvec_my_length(zvec,&nloc_z,&iflag);
  if (iflag!=PHIST_SUCCESS) return -1.0e29;
  
  if (nvec_r!=nvec_i || nvec_r!=nvec_z) return -1.0e30;
  if (nloc_r!=nloc_i || nloc_r!=nloc_z) return -1.0e31;
#ifdef PHIST_MVECS_ROW_MAJOR
  double zval_r[nloc_z*lda_r], zval_i[nloc_z*lda_i];
#else
  double zval_r[nvec_z*lda_r], zval_i[nvec_z*lda_i];
#endif
  for (phist_lidx i=0; i<nloc_z; i++)
  {
    for (int j=0; j<nvec_z; j++)
    {
      zval_r[VIDX(i,j,lda_r)] = std::real(zval[VIDX(i,j,lda_z)]);
      zval_i[VIDX(i,j,lda_i)] = std::imag(zval[VIDX(i,j,lda_z)]);
    }
//    PHIST_SOUT(PHIST_DEBUG,"%d\t%16.8e%+16.8ei  <-> %16.8e%+16.8ei\n", 
//        i,zval_r[VIDX(i,0,lda_r)],zval_i[VIDX(i,0,lda_r)],
//        dval_r[VIDX(i,0,lda_r)],dval_i[VIDX(i,0,lda_r)]);
  }
  double result_r = TestWithType<double>::ArraysEqual(dval_r,zval_r,nloc_z,nvec_z,lda_r,1,KernelTest::vflag_,relTo);
  double result_i = TestWithType<double>::ArraysEqual(dval_i,zval_i,nloc_z,nvec_z,lda_r,1,KernelTest::vflag_,relTo);
  if (std::abs(result_r-1.0)>=std::abs(result_i-1.0))
  { 
    return result_r; 
  }
  else
  { 
    return result_i; 
  }
}

void MvecCopyX2Z(phist_Dx_mvec *xvec, phist_Zmvec *zvec, int* iflag)
{
  *iflag=0;
  PHIST_CHK_IERR(phist_Dmvec_from_device(xvec->v_,iflag),*iflag);
  PHIST_CHK_IERR(phist_Dmvec_from_device(xvec->vi_,iflag),*iflag);
  PHIST_CHK_IERR(phist_Zmvec_from_device(zvec,iflag),*iflag);

  double *dval_r, *dval_i;
  phist_d_complex *zval;
  phist_lidx lda_z, lda_r, lda_i;
  int nvec_r, nvec_i, nvec_z;
  phist_lidx nloc_r, nloc_i, nloc_z;
  
  PHIST_CHK_IERR(phist_Dmvec_extract_view(xvec->v_,&dval_r,&lda_r,iflag),*iflag);
  PHIST_CHK_IERR(phist_Dmvec_extract_view(xvec->vi_,&dval_i,&lda_i,iflag),*iflag);
  PHIST_CHK_IERR(phist_Zmvec_extract_view(zvec,&zval,&lda_z,iflag),*iflag);

  PHIST_CHK_IERR(phist_Dmvec_num_vectors(xvec->v_,&nvec_r,iflag),*iflag);
  PHIST_CHK_IERR(phist_Dmvec_num_vectors(xvec->vi_,&nvec_i,iflag),*iflag);
  PHIST_CHK_IERR(phist_Zmvec_num_vectors(zvec,&nvec_z,iflag),*iflag);

  PHIST_CHK_IERR(phist_Dmvec_my_length(xvec->v_,&nloc_r,iflag),*iflag);
  PHIST_CHK_IERR(phist_Dmvec_my_length(xvec->vi_,&nloc_i,iflag),*iflag);
  PHIST_CHK_IERR(phist_Zmvec_my_length(zvec,&nloc_z,iflag),*iflag);
  if (nvec_r!=nvec_i || nvec_r!=nvec_z) {*iflag=-1; return;}
  if (nloc_r!=nloc_i || nloc_r!=nloc_z) {*iflag=-2; return;}

  for (phist_lidx i=0; i<nloc_z; i++)
  {
    for (int j=0; j<nvec_z; j++)
    {
      zval[VIDX(i,j,lda_z)] =                           dval_r[VIDX(i,j,lda_r)] + 
        phist::ScalarTraits<phist_d_complex>::cmplx_I()*dval_i[VIDX(i,j,lda_i)];
    }
  }
  PHIST_CHK_IERR(phist_Zmvec_to_device(zvec,iflag),*iflag);
}
#endif

int zshift_bench3D(ghost_gidx row, ghost_lidx *nnz, ghost_gidx* cols, void* vals, void* data)
{
  int iflag=0;
  phist_d_complex* shift = (phist_d_complex*)data;
  phist_d_complex* zvals = (phist_d_complex*)vals;
  double dvals[7];
  PHIST_ICHK_IERR(iflag=MATPDE3D_rowFunc(row,nnz,cols,dvals,NULL),iflag);
  for (int i=0; i<*nnz; i++) 
  {
    zvals[i]=(phist_d_complex)dvals[i];
    if (cols[i]==row && shift!=NULL) zvals[i] -= (*shift);
  }
  return iflag;
}

int zshift_idfunc(ghost_gidx row, ghost_lidx *nnz, ghost_gidx* cols, void* vals, void* data)
{
  phist_d_complex* shift = (phist_d_complex*)data;
  phist_d_complex* zvals = (phist_d_complex*)vals;
  *nnz=1;
  cols[0]=row;
  zvals[0]= (phist_d_complex)1.0;
  if (shift!=NULL) zvals[0] -= (*shift);
  return 0;
}


#define CLASSFILE_DEF "CarpTest_def.hpp"

#define MATNAME MATNAME_BENCH3D_8_A1
#define ZMATFUNC zshift_bench3D
#define _BASENAME_ CarpTest_A

#define _N_ 512
#define _M_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 512
#define _M_ 4
#include "../phist_typed_test_gen.h"

#define _N_ 512
#define _M_ 7
#include "../phist_typed_test_gen.h"

#undef MATNAME
#define MATNAME MATNAME_IDFUNC
#undef ZMATFUNC
#define ZMATFUNC zshift_idfunc
#undef _BASENAME_
#define _BASENAME_ CarpTest_I

#define _N_ 99
#define _M_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 99
#define _M_ 4
#include "../phist_typed_test_gen.h"

#define _N_ 99
#define _M_ 7
#include "../phist_typed_test_gen.h"
