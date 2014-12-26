
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

int iflag;
SUBR(type_avail)(&iflag);
typeImplemented_=(iflag==0);

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
  int iflag;
  _ST_* val;
  lidx_t lda,n;
  int m;
  
  SUBR(mvec_from_device)(V,&iflag);
  if (iflag!=PHIST_SUCCESS) return (_MT_)(-mt::one());
  SUBR(mvec_extract_view)(V,&val,&lda,&iflag);
  if (iflag!=PHIST_SUCCESS) return (_MT_)(-2*mt::one());
  SUBR(mvec_my_length)(V,&n,&iflag);
  if (iflag!=PHIST_SUCCESS) return (_MT_)(-2*mt::one());
  SUBR(mvec_num_vectors)(V,&m,&iflag);
  if (iflag!=PHIST_SUCCESS) return (_MT_)(-2*mt::one());
  
  return ArrayEqual(val,n,m,lda,1,value,KernelTest::vflag_);
  }

static _MT_ MvecsEqual(TYPE(mvec_ptr) V1, TYPE(mvec_ptr) V2)
  {
  int iflag;
  _ST_ *val,*val2;
  lidx_t lda,n,lda2,n2;
  int m,m2;
  
  SUBR(mvec_from_device)(V1,&iflag);
  if (iflag!=PHIST_SUCCESS) return (_MT_)(-mt::one());
  SUBR(mvec_from_device)(V2,&iflag);
  if (iflag!=PHIST_SUCCESS) return (_MT_)(-mt::one());
  
  SUBR(mvec_extract_view)(V1,&val,&lda,&iflag);
  if (iflag!=PHIST_SUCCESS) return (_MT_)(-2*mt::one());
  SUBR(mvec_extract_view)(V1,&val2,&lda2,&iflag);
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
  if (n!=n2||m!=m2) return (_MT_)(mt::one());
  if (lda!=lda2) return (_MT_)(-99*mt::one()); // test not implemented
  return ArraysEqual(val,val2,n,m,lda,1,KernelTest::vflag_);
  }

static _MT_ SdMatEqual(TYPE(mvec_ptr) V, _ST_ value)
  {
  int iflag;
  _ST_* val;
  lidx_t lda;
  int n,m;
  
  SUBR(sdMat_from_device)(V,&iflag);
  if (iflag!=PHIST_SUCCESS) return (_MT_)(-mt::one());
  SUBR(sdMat_extract_view)(V,&val,&lda,&iflag);
  if (iflag!=PHIST_SUCCESS) return (_MT_)(-2*mt::one());
  SUBR(sdMat_get_nrows)(V,&n,&iflag);
  if (iflag!=PHIST_SUCCESS) return (_MT_)(-2*mt::one());
  SUBR(sdMat_get_ncols)(V,&m,&iflag);
  if (iflag!=PHIST_SUCCESS) return (_MT_)(-2*mt::one());
  
  return ArrayEqual(val,n,m,lda,1,value,KernelTest::mflag_);
  }

static _MT_ SdMatsEqual(TYPE(sdMat_ptr) V1, TYPE(sdMat_ptr) V2)
  {
  int iflag;
  _ST_ *val, *val2;
  lidx_t lda,lda2;
  int n,m,n2,m2;
  
  SUBR(sdMat_from_device)(V1,&iflag);
  if (iflag!=PHIST_SUCCESS) return (_MT_)(-mt::one());
  SUBR(sdMat_from_device)(V2,&iflag);
  if (iflag!=PHIST_SUCCESS) return (_MT_)(-mt::one());
  
  SUBR(sdMat_extract_view)(V1,&val,&lda,&iflag);
  if (iflag!=PHIST_SUCCESS) return (_MT_)(-2*mt::one());
  SUBR(sdMat_extract_view)(V1,&val2,&lda2,&iflag);
  if (iflag!=PHIST_SUCCESS) return (_MT_)(-2*mt::one());
  
  SUBR(sdMat_get_nrows)(V1,&n,&iflag);
  if (iflag!=PHIST_SUCCESS) return (_MT_)(-2*mt::one());
  SUBR(sdMat_get_nrows)(V2,&n2,&iflag);
  if (iflag!=PHIST_SUCCESS) return (_MT_)(-2*mt::one());

  SUBR(sdMat_get_ncols)(V1,&m,&iflag);
  if (iflag!=PHIST_SUCCESS) return (_MT_)(-2*mt::one());
  SUBR(sdMat_get_ncols)(V2,&m2,&iflag);
  if (iflag!=PHIST_SUCCESS) return (_MT_)(-2*mt::one());
 
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

static _MT_ ArrayParallelReplicated(const _ST_* array, int n, int m, lidx_t lda, lidx_t stride, bool swap_nm=false)
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
      maxval=std::max(st::abs(array[j*lda+i*stride]-buff[i*M+j])/scal,maxval);
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
  return (2*((MT)std::rand()/(MT)RAND_MAX)-1) +
         (2*((MT)std::rand()/(MT)RAND_MAX)-1) * st::cmplx_I();
  }

};


