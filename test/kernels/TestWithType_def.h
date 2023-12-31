/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
/** specialization of TestWithType for known types
 */
template<>
class TestWithType< _ST_ > : public virtual testing::Test
{
public:

// including these typedefs allows us to write
// things like st::sqrt in derived classes
#include "phist_std_typedefs.hpp"

static bool typeImplemented_;

static void SetUpTestCase()
{
  int iflag=-1;
  SUBR(type_avail)(&iflag);
  typeImplemented_=(iflag==0);

  phist_random_init();
}

static void TearDownTestCase()
{
}
static _MT_ MvecEqual(TYPE(mvec_ptr) V, _ST_ value)
{
  int iflag;
  _ST_* val;
  phist_lidx lda,n;
  int m;
  
  SUBR(mvec_from_device)(V,&iflag);
  if (iflag!=PHIST_SUCCESS) return (_MT_)(-mt::one());
  SUBR(mvec_extract_view)(V,&val,&lda,&iflag);
  if (iflag!=PHIST_SUCCESS) return (_MT_)(-2*mt::one());
  SUBR(mvec_my_length)(V,&n,&iflag);
  if (iflag!=PHIST_SUCCESS) return (_MT_)(-2*mt::one());
  SUBR(mvec_num_vectors)(V,&m,&iflag);
  if (iflag!=PHIST_SUCCESS) return (_MT_)(-2*mt::one());
  
  _MT_ return_value;
  if (n*m>100000)
  {
    return_value=ArrayEqualOMP(val,n,m,lda,1,value,KernelTest::vflag_);
  }
  else
  {
    return_value=ArrayEqual(val,n,m,lda,1,value,KernelTest::vflag_);
  }
#if PHIST_OUTLEV>=PHIST_DEBUG
  int print_vec_loc,print_vec;
  print_vec_loc=(std::abs(return_value-mt::one())>std::sqrt(mt::eps())&&(n<=100))?1:0;
  MPI_Allreduce(&print_vec_loc,&print_vec,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  if (print_vec)
  {
    PHIST_SOUT(PHIST_DEBUG,"MvecEqual, val=%8.4e%+8.4ei, vec:\n",st::real(value),st::imag(value));
    SUBR(mvec_print)(V,&iflag);
  }
#endif
  return return_value;
}

static int MvecContainsInfOrNaN(TYPE(mvec_ptr) V)
{
  int iflag;
  _ST_* val;
  phist_lidx lda,n;
  int m;
  
  PHIST_ICHK_IERR(SUBR(mvec_from_device)(V,&iflag),iflag);
  PHIST_ICHK_IERR(SUBR(mvec_extract_view)(V,&val,&lda,&iflag),iflag);
  PHIST_ICHK_IERR(SUBR(mvec_my_length)(V,&n,&iflag),iflag);
  PHIST_ICHK_IERR(SUBR(mvec_num_vectors)(V,&m,&iflag),iflag);
  
  return ArrayContainsInfOrNaN(val,n,m,lda,1,KernelTest::vflag_);
}

static _MT_ MvecsEqual(TYPE(mvec_ptr) V1, TYPE(mvec_ptr) V2, _MT_ relTo = mt::zero())
{
  int iflag;
  _ST_ *val,*val2;
  phist_lidx lda,n,lda2,n2;
  int m,m2;
  
  SUBR(mvec_from_device)(V1,&iflag);
  if (iflag!=PHIST_SUCCESS) return (_MT_)(-mt::one());
  SUBR(mvec_from_device)(V2,&iflag);
  if (iflag!=PHIST_SUCCESS) return (_MT_)(-mt::one());
  
  SUBR(mvec_extract_view)(V1,&val,&lda,&iflag);
  if (iflag!=PHIST_SUCCESS) return (_MT_)(-2*mt::one());
  SUBR(mvec_extract_view)(V2,&val2,&lda2,&iflag);
  if (iflag!=PHIST_SUCCESS) return (_MT_)(-2*mt::one());
  
  SUBR(mvec_my_length)(V1,&n,&iflag);
  if (iflag!=PHIST_SUCCESS) return (_MT_)(-2*mt::one());
  SUBR(mvec_my_length)(V2,&n2,&iflag);
  if (iflag!=PHIST_SUCCESS) return (_MT_)(-2*mt::one());

  SUBR(mvec_num_vectors)(V1,&m,&iflag);
  if (iflag!=PHIST_SUCCESS) return (_MT_)(-2*mt::one());
  SUBR(mvec_num_vectors)(V2,&m2,&iflag);
  if (iflag!=PHIST_SUCCESS) return (_MT_)(-2*mt::one());
 
  // vectors not equal: dimensions mismatch
  if (n!=n2||m!=m2) return (_MT_)(-mt::one());
  MT return_value=mt::one();
  if (n*m>100000)
  {
    if (lda!=lda2)
    {
      return_value=ArraysEqualWithDifferentLDA_OMP(val,val2,n,m,lda,lda2,1,KernelTest::vflag_,relTo);
    }
    else
    {
      return_value=ArraysEqualOMP(val,val2,n,m,lda,1,KernelTest::vflag_,relTo);
    }
  }
  else
  {
    if (lda!=lda2)
    {
      return_value=ArraysEqualWithDifferentLDA(val,val2,n,m,lda,lda2,1,KernelTest::vflag_,relTo);
    }
    else
    {
      return_value=ArraysEqual(val,val2,n,m,lda,1,KernelTest::vflag_,relTo);
    }
  }
#if PHIST_OUTLEV>=PHIST_DEBUG
  int print_vecs_loc,print_vecs;
  print_vecs_loc=(std::abs(return_value-mt::one())>std::sqrt(mt::eps())&&(n<=100))?1:0;
  MPI_Allreduce(&print_vecs_loc,&print_vecs,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  if (print_vecs)
  {
    PHIST_SOUT(PHIST_DEBUG,"MvecsEqual vec1:\n");
    SUBR(mvec_print)(V1,&iflag);
    PHIST_SOUT(PHIST_DEBUG,"MvecsEqual vec2:\n");
    SUBR(mvec_print)(V2,&iflag);
  }
#endif
  return return_value;
}

// compare the *HOST SIDE* of an sdMat with a scalar constant
static _MT_ SdMatEqual(TYPE(sdMat_ptr) M, _ST_ value)
{
  int iflag;
  _ST_* val;
  phist_lidx lda;
  int n,m;
  
  SUBR(sdMat_extract_view)(M,&val,&lda,&iflag);
  if (iflag!=PHIST_SUCCESS) return (_MT_)(-2*mt::one());
  SUBR(sdMat_get_nrows)(M,&n,&iflag);
  if (iflag!=PHIST_SUCCESS) return (_MT_)(-2*mt::one());
  SUBR(sdMat_get_ncols)(M,&m,&iflag);
  if (iflag!=PHIST_SUCCESS) return (_MT_)(-2*mt::one());
  
  MT return_value=ArrayEqual(val,n,m,lda,1,value,KernelTest::mflag_);
#if PHIST_OUTLEV>=PHIST_DEBUG
  if (std::abs(return_value-mt::one())>std::sqrt(mt::eps()))
  {
    PHIST_SOUT(PHIST_DEBUG,"SdMatEqual, val=%8.4e%+8.4ei,  mat:\n",st::real(value),st::imag(value));
    SUBR(sdMat_print)(M,&iflag);
  }
#endif
  return return_value;
}

static int SdMatContainsInfOrNaN(TYPE(sdMat_ptr) M)
{
  int iflag;
  _ST_* val;
  phist_lidx lda;
  int n,m;
  
  PHIST_ICHK_IERR(SUBR(sdMat_extract_view)(M,&val,&lda,&iflag),iflag);
  PHIST_ICHK_IERR(SUBR(sdMat_get_nrows)(M,&n,&iflag),iflag);
  PHIST_ICHK_IERR(SUBR(sdMat_get_ncols)(M,&m,&iflag),iflag);
  
  return ArrayContainsInfOrNaN(val,n,m,lda,1,KernelTest::mflag_);
}

// compare the *HOST SIDE* of two sdMats
static _MT_ SdMatsEqual(TYPE(sdMat_ptr) M1, TYPE(sdMat_ptr) M2, _MT_ relTo = mt::zero())
{
  int iflag;
  _ST_ *val, *val2;
  phist_lidx lda,lda2;
  int n,m,n2,m2;
  
  SUBR(sdMat_extract_view)(M1,&val,&lda,&iflag);
  if (iflag!=PHIST_SUCCESS) return (_MT_)(-2*mt::one());
  SUBR(sdMat_extract_view)(M2,&val2,&lda2,&iflag);
  if (iflag!=PHIST_SUCCESS) return (_MT_)(-2*mt::one());
  
  SUBR(sdMat_get_nrows)(M1,&n,&iflag);
  if (iflag!=PHIST_SUCCESS) return (_MT_)(-2*mt::one());
  SUBR(sdMat_get_nrows)(M2,&n2,&iflag);
  if (iflag!=PHIST_SUCCESS) return (_MT_)(-2*mt::one());

  SUBR(sdMat_get_ncols)(M1,&m,&iflag);
  if (iflag!=PHIST_SUCCESS) return (_MT_)(-2*mt::one());
  SUBR(sdMat_get_ncols)(M2,&m2,&iflag);
  if (iflag!=PHIST_SUCCESS) return (_MT_)(-2*mt::one());
 
  // vectors not equal: dimensions mismatch
  if (n!=n2||m!=m2) return (_MT_)(-mt::one());

  MT return_value;

  if (lda!=lda2)
  {
    return_value=ArraysEqualWithDifferentLDA(val,val2,n,m,lda,lda2,1,KernelTest::mflag_,relTo);
  }
  else
  {
    return_value=ArraysEqual(val,val2,n,m,lda,1,KernelTest::mflag_,relTo);
  }
#if PHIST_OUTLEV>=PHIST_DEBUG
  if (std::abs(return_value-mt::one())>std::sqrt(mt::eps()))
  {
    PHIST_SOUT(PHIST_DEBUG,"SdMatsEqual mat1:\n");
    SUBR(sdMat_print)(M1,&iflag);
    PHIST_SOUT(PHIST_DEBUG,"SdMatsEqual mat2:\n");
    SUBR(sdMat_print)(M2,&iflag);
  }
#endif

  return return_value;
}

static int ArrayContainsInfOrNaN(const _ST_* array, int n, int m, phist_lidx lda, phist_lidx stride, bool swap_nm=false)
{
  MT *A=(MT*)array;
  int N = swap_nm? m: n;
  int M = swap_nm? n: m;
  int dt_size = (int)(sizeof(ST)/sizeof(MT));
  for (int i=0;i<N*stride;i+=stride)
  {
    for (int j=0; j<M;j++)
    {
      for (int k=0; k<dt_size; k++)
      {
        if (std::isnan(A[(j*lda+i)*dt_size+k])) return 1;
        if (std::isinf(A[(j*lda+i)*dt_size+k])) return 2;
      }
    }
  }
  return 0;
}

static _MT_ ArrayEqual(const _ST_* array, int n, int m, phist_lidx lda, phist_lidx stride, _ST_ value, bool swap_nm=false)
{
  MT maxval=mt::zero();
  MT scal= st::abs(value);
  int N = swap_nm? m: n;
  int M = swap_nm? n: m;
  if (scal==mt::zero()) scal=mt::one();
  for (int i=0;i<N*stride;i+=stride)
    {
    for (int j=0; j<M;j++)
      {
      maxval=mt::max(st::abs(array[j*lda+i]-value)/scal,maxval);
      }
    }
  return (MT)1.0+maxval;
}

// a faster variant for multi core CPUs that allows us to test larger cases
static inline _MT_ ArrayEqualOMP(const _ST_* array, int n, int m, phist_lidx lda, phist_lidx stride, _ST_ value, bool swap_nm=false)
{
  MT maxval=mt::zero();
  MT scal= st::abs(value);
  int N = swap_nm? m: n;
  int M = swap_nm? n: m;
  if (scal==mt::zero()) scal=mt::one();
// Intel produces "multiple definitions of ..." errors for OpenMP clauses here (in a header included multiple times)
#ifndef __INTEL_COMPILER
#pragma omp parallel for reduction(max:maxval) schedule(static)
#endif
  for (int i=0;i<N*stride;i+=stride)
  {
    for (int j=0; j<M;j++)
    {
      maxval=mt::max(st::abs(array[j*lda+i]-value)/scal,maxval);
    }
  }
  return (MT)1.0+maxval;
}


static _MT_ ArrayParallelReplicated(const _ST_* array, int n, int m, phist_lidx lda, phist_lidx stride, bool swap_nm=false)
{
  _ST_ buff[n*m];
  MT maxval=mt::zero();
  int N = swap_nm? m: n;
  int M = swap_nm? n: m;
  int rank = 0;
#ifdef PHIST_HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  {
    int nglobal_min, nglobal_max, mglobal_min, mglobal_max;
    MPI_Allreduce(&n, &nglobal_min, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&n, &nglobal_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&m, &mglobal_min, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&m, &mglobal_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if( nglobal_min != nglobal_max || mglobal_min != mglobal_max )
      return (MT)-1.0;
  }
#endif

  if( rank == 0 )
  {
    for (int i=0;i<N;i++)
      for (int j=0; j<M;j++)
        buff[i*M+j] = array[j*lda+i*stride];
  }
#ifdef PHIST_HAVE_MPI
  MPI_Bcast(buff, sizeof(_ST_)*n*m, MPI_BYTE, 0, MPI_COMM_WORLD);
#endif
  for (int i=0;i<N;i++)
  {
    for (int j=0; j<M;j++)
    {
      MT scal= st::abs(buff[i*M+j]);
      if (scal==mt::zero()) scal=mt::one();
      maxval=mt::max(st::abs(array[j*lda+i*stride]-buff[i*M+j])/scal,maxval);
    }
  }
  return (MT)1.0+maxval;
}

static _MT_ ArraysEqual(const _ST_* arr1,const _ST_* arr2, int n, int m, phist_lidx lda, phist_lidx stride, bool swap_n_m=false, _MT_ relTo=mt::zero())
{
    return ArraysEqualWithDifferentLDA(arr1,arr2,n,m,lda,lda,stride,swap_n_m,relTo);
}

static _MT_ ArraysEqualOMP(const _ST_* arr1,const _ST_* arr2, int n, int m, phist_lidx lda, phist_lidx stride, bool swap_n_m=false, _MT_ relTo=mt::zero())
{
    return ArraysEqualWithDifferentLDA_OMP(arr1,arr2,n,m,lda,lda,stride,swap_n_m,relTo);
}
  
static _MT_ ArraysEqualWithDifferentLDA(const _ST_* arr1,const _ST_* arr2, int n, int m, 
phist_lidx lda1, phist_lidx lda2, phist_lidx stride, bool swap_n_m=false, _MT_ relTo = mt::zero())
{
  int N = swap_n_m? m: n;
  int M = swap_n_m? n: m;
  _MT_ maxval=mt::zero();
#ifdef PHIST_TESTING
  int max_i=0,max_j=0;
#endif
  for (int i=0;i<N*stride;i+=stride)
  {
    for (int j=0;j<M;j++)
    {
      MT mn = st::abs(arr1[j*lda1+i]-arr2[j*lda2+i]);
      MT pl = relTo;
      if (pl==mt::zero()) pl = (st::abs(arr1[j*lda1+i])+st::abs(arr2[j*lda2+i]))*(MT)0.5;
      if (pl==mt::zero()) pl=mt::one();
      maxval=mt::max(mn/pl,maxval);
#ifdef PHIST_TESTING
      if (maxval==mn/pl) 
      {
        max_i=i;
        max_j=j;
      }
#endif
    }
  }
#ifdef PHIST_TESTING
  if (maxval>std::sqrt(mt::eps()))
  {
    PHIST_OUT(PHIST_INFO,"max relative difference is %e in location (i=%d,j=%d)\n",maxval,max_i,max_j);
#ifdef IS_COMPLEX
    PHIST_OUT(PHIST_INFO,"a1(%d,%d)=%16.8e%+16.8ei, a2(%d,%d)=%16.8e%+16.8ei\n",
        max_i,max_j, st::real(arr1[max_j*lda1+max_i]), st::imag(arr1[max_j*lda1+max_i]),
        max_i,max_j, st::real(arr2[max_j*lda2+max_i]), st::imag(arr2[max_j*lda2+max_i]));
#else
    PHIST_OUT(PHIST_INFO,"a1(%d,%d)=%16.8e, a2(%d,%d)=%16.8e\n",max_i,max_j,
        arr1[max_j*lda1+max_i],max_i,max_j,arr2[max_j*lda2+max_i]);
#endif
    PHIST_OUT(PHIST_INFO,"With lda1=%d, lda2=%d, N=%d, M=%d\n",lda1,lda2,n,m);
  }
#endif
  return mt::one()+maxval;
}

static inline _MT_ ArraysEqualWithDifferentLDA_OMP(const _ST_* arr1,const _ST_* arr2, int n, int m, 
phist_lidx lda1, phist_lidx lda2, phist_lidx stride, bool swap_n_m=false, _MT_ relTo = mt::zero())
{
  int N = swap_n_m? m: n;
  int M = swap_n_m? n: m;
  _MT_ maxval=mt::zero();
  std::stringstream ss_deb;
// Intel produces "multiple definitions of ..." errors for OpenMP clauses here (in a header included multiple times)
#ifndef __INTEL_COMPILER
#pragma omp parallel for reduction(max:maxval) schedule(static)
#endif
  for (int i=0;i<N*stride;i+=stride)
  {
    for (int j=0;j<M;j++)
    {
      MT mn = st::abs(arr1[j*lda1+i]-arr2[j*lda2+i]);
      MT pl = relTo;
      if (pl==mt::zero()) pl = (st::abs(arr1[j*lda1+i])+st::abs(arr2[j*lda2+i]))*(MT)0.5;
      if (pl==mt::zero()) pl=mt::one();
      maxval=mt::max(mn/pl,maxval);
    }
  }
  //std::cout << "MAX VAL: "<<maxval<<std::endl;
  return mt::one()+maxval;
}

static inline _ST_ random_number() 
{
  return (2*((MT)std::rand()/(MT)RAND_MAX)-1) +
         (2*((MT)std::rand()/(MT)RAND_MAX)-1) * st::cmplx_I();
}

};


