
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


