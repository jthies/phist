
/** 
 */
template<>
class KernelTestWithType< _ST_ >
{
public:

// including these typedefs allows us to write
// things like st::sqrt in derived classes
#include "phist_std_typedefs.hpp"

bool typeImplemented_;
bool verbose_;

/** Set up method.
 * Fills internal data vector with values 1.0, 2.0 and 3.0.
 */
virtual void SetUp() {

int ierr;
SUBR(type_avail)(&ierr);
typeImplemented_=(ierr==0);

int rank=0;
#ifdef PHIST_HAVE_MPI
MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#endif
verbose_=rank==0;

if (verbose_ && false)
  {
  std::cout << "data type: ";
#ifdef IS_COMPLEX
  std::cout << " complex ";
#else
  std::cout << " real ";
#endif
#ifdef _IS_DOUBLE_
  std::cout << "double" << std::endl;
#else
  std::cout << "float" << std::endl;
#endif
  }
}

virtual void TearDown() {
}

static _MT_ MvecEqual(TYPE(mvec_ptr) V, _ST_ value)
  {
  int ierr;
  _ST_* val;
  lidx_t lda,n,m;
  
  SUBR(mvec_from_device)(V,&ierr);
  if (ierr!=PHIST_SUCCESS) return (_MT_)(-mt::one());
  SUBR(mvec_extract_view)(V,&val,&lda,&ierr);
  if (ierr!=PHIST_SUCCESS) return (_MT_)(-2*mt::one());
  SUBR(mvec_my_length)(V,&n,&ierr);
  if (ierr!=PHIST_SUCCESS) return (_MT_)(-2*mt::one());
  SUBR(mvec_num_vectors)(V,&m,&ierr);
  if (ierr!=PHIST_SUCCESS) return (_MT_)(-2*mt::one());
  
  return ArrayEqual(val,n,m,lda,1,value,KernelTest::vflag_);
  }

static _MT_ MvecsEqual(TYPE(mvec_ptr) V1, TYPE(mvec_ptr) V2)
  {
  int ierr;
  _ST_ *val,*val2;
  lidx_t lda,n,m,lda2,n2,m2;
  
  SUBR(mvec_from_device)(V1,&ierr);
  if (ierr!=PHIST_SUCCESS) return (_MT_)(-mt::one());
  SUBR(mvec_from_device)(V2,&ierr);
  if (ierr!=PHIST_SUCCESS) return (_MT_)(-mt::one());
  
  SUBR(mvec_extract_view)(V1,&val,&lda,&ierr);
  if (ierr!=PHIST_SUCCESS) return (_MT_)(-2*mt::one());
  SUBR(mvec_extract_view)(V1,&val2,&lda2,&ierr);
  if (ierr!=PHIST_SUCCESS) return (_MT_)(-2*mt::one());
  
  SUBR(mvec_my_length)(V1,&n,&ierr);
  if (ierr!=PHIST_SUCCESS) return (_MT_)(-2*mt::one());
  SUBR(mvec_my_length)(V2,&n2,&ierr);
  if (ierr!=PHIST_SUCCESS) return (_MT_)(-2*mt::one());

  SUBR(mvec_num_vectors)(V1,&m,&ierr);
  if (ierr!=PHIST_SUCCESS) return (_MT_)(-2*mt::one());
  SUBR(mvec_num_vectors)(V2,&m2,&ierr);
  if (ierr!=PHIST_SUCCESS) return (_MT_)(-2*mt::one());
 
  // vectors not equal: dimensions mismatch
  if (n!=n2||m!=m2) return (_MT_)(mt::one());
  if (lda!=lda2) return (_MT_)(-99*mt::one()); // test not implemented
  return ArraysEqual(val,val2,n,m,lda,1,KernelTest::vflag_);
  }

static _MT_ SdMatEqual(TYPE(mvec_ptr) V, _ST_ value)
  {
  int ierr;
  _ST_* val;
  lidx_t lda,n,m;
  
  SUBR(sdMat_from_device)(V,&ierr);
  if (ierr!=PHIST_SUCCESS) return (_MT_)(-mt::one());
  SUBR(sdMat_extract_view)(V,&val,&lda,&ierr);
  if (ierr!=PHIST_SUCCESS) return (_MT_)(-2*mt::one());
  SUBR(sdMat_get_nrows)(V,&n,&ierr);
  if (ierr!=PHIST_SUCCESS) return (_MT_)(-2*mt::one());
  SUBR(sdMat_get_ncols)(V,&m,&ierr);
  if (ierr!=PHIST_SUCCESS) return (_MT_)(-2*mt::one());
  
  return ArrayEqual(val,n,m,lda,1,value,KernelTest::mflag_);
  }

static _MT_ SdMatsEqual(TYPE(sdMat_ptr) V1, TYPE(sdMat_ptr) V2)
  {
  int ierr;
  _ST_ *val, *val2;
  lidx_t lda,n,m,lda2,n2,m2;
  
  SUBR(sdMat_from_device)(V1,&ierr);
  if (ierr!=PHIST_SUCCESS) return (_MT_)(-mt::one());
  SUBR(sdMat_from_device)(V2,&ierr);
  if (ierr!=PHIST_SUCCESS) return (_MT_)(-mt::one());
  
  SUBR(sdMat_extract_view)(V1,&val,&lda,&ierr);
  if (ierr!=PHIST_SUCCESS) return (_MT_)(-2*mt::one());
  SUBR(sdMat_extract_view)(V1,&val2,&lda2,&ierr);
  if (ierr!=PHIST_SUCCESS) return (_MT_)(-2*mt::one());
  
  SUBR(sdMat_get_nrows)(V1,&n,&ierr);
  if (ierr!=PHIST_SUCCESS) return (_MT_)(-2*mt::one());
  SUBR(sdMat_get_nrows)(V2,&n2,&ierr);
  if (ierr!=PHIST_SUCCESS) return (_MT_)(-2*mt::one());

  SUBR(sdMat_get_ncols)(V1,&m,&ierr);
  if (ierr!=PHIST_SUCCESS) return (_MT_)(-2*mt::one());
  SUBR(sdMat_get_ncols)(V2,&m2,&ierr);
  if (ierr!=PHIST_SUCCESS) return (_MT_)(-2*mt::one());
 
  // vectors not equal: dimensions mismatch
  if (n!=n2||m!=m2) return (_MT_)(mt::one());
  if (lda!=lda2) return (_MT_)(-99*mt::one()); // test not implemented
  return ArraysEqual(val,val2,n,m,lda,1,KernelTest::mflag_);
  }

static _MT_ ArrayEqual(const _ST_* array, int n, int m, lidx_t lda, lidx_t stride, _ST_ value, bool swap_nm=false)
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
      maxval=std::max(st::abs(array[j*lda+i]-value)/scal,maxval);
      }
    }
  return (MT)1.0+maxval;
  }


static _MT_ ArraysEqual(const _ST_* arr1,const _ST_* arr2, int n, int m, lidx_t lda, lidx_t stride, bool swap_n_m=false)
  {
  int N = swap_n_m? m: n;
  int M = swap_n_m? n: m;
  _MT_ maxval=mt::zero();
  for (int i=0;i<N*stride;i+=stride)
    {
    for (int j=0;j<M;j++)
      {
      MT mn = st::abs(arr1[j*lda+i]-arr2[j*lda+i]);
      MT pl = (st::abs(arr1[j*lda+i])+st::abs(arr2[j*lda+i]))*(MT)0.5;
      if (pl==mt::zero()) pl=mt::one();
      maxval=std::max(mn/pl,maxval);
      //std::cout << arr1[j*lda+i]<< "\t SAME ?? AS \t"<<arr2[j*lda+i]<<std::endl;
      }
//    std::cout << std::endl;
    }
  //std::cout << "MAX VAL: "<<maxval<<std::endl;
  return mt::one()+maxval;
  }

static inline _ST_ random_number() 
  {
  return (MT)(2*std::rand()-RAND_MAX)/(MT)RAND_MAX + 
         (MT)(2*std::rand()-RAND_MAX)/(MT)RAND_MAX * st::cmplx_I();
  }

};


