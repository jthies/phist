#include "../tools/TestHelpers.h"
/*! Test fixture. */
template<phist_gidx _Nglob, int _Nvec, int _useViews, int _numberOfVectorsInitialized, int _multipleDefinitionCounter>
class KernelTestWithVectors<_ST_,_Nglob,_Nvec, _useViews, _numberOfVectorsInitialized, _multipleDefinitionCounter> : 
        public virtual KernelTestWithMap<_Nglob> ,
        public virtual TestWithType< _ST_ >
  {

public:
  // make some stuff known to the compiler to prevent errors about undeclared specifiers
  using TestWithType<_ST_>::typeImplemented_;
  using KernelTestWithMap<_Nglob>::map_;
  using KernelTestWithMap<_Nglob>::nloc_;
  using KernelTestWithMap<_Nglob>::problemTooSmall_;
  using KernelTestWithMap<_Nglob>::iflag_;
  using KernelTestWithMap<_Nglob>::vflag_;

static void SetUpTestCase()
{
  KernelTestWithMap<_Nglob>::SetUpTestCase();
  TestWithType<_ST_>::SetUpTestCase();


  if (typeImplemented_ && !problemTooSmall_)
  {
    bool align_p2=(useViews_==2);
    pad_pre_=0, pad_post_=0;
    nvecPadded_=nvec_;
    if (useViews_)
    {
      pad_pre_=align_p2?8:3;
      pad_post_=7;
      // padding to power of 2
      int pow_p2=(int)ceil(log((double)(nvec_+pad_pre_+pad_post_))/log(2.0))+1;
      nvecPadded_=align_p2?(int)pow(2.0,(double)pow_p2): nvec_+pad_pre_+pad_post_;
      pad_post_ = nvecPadded_-pad_pre_-nvec_;
    }
    if( _numberOfVectorsInitialized >= 1 )
    {
      PHISTTEST_MVEC_CREATE(&mem1_,map_,nvecPadded_,&iflag_);
      ASSERT_EQ(0,iflag_);
    }

    if( _numberOfVectorsInitialized >= 2 )
    {
      PHISTTEST_MVEC_CREATE(&mem2_,map_,nvecPadded_,&iflag_);
      ASSERT_EQ(0,iflag_);
    }

    if( _numberOfVectorsInitialized >= 3 )
    {
      PHISTTEST_MVEC_CREATE(&mem3_,map_,nvecPadded_,&iflag_);
      ASSERT_EQ(0,iflag_);
    }

    // if requested, set vecX to views of memX
    if (useViews_)
    {
      vec1_=NULL; vec2_=NULL; vec3_=NULL;
      if( _numberOfVectorsInitialized >= 1 )
      {
        SUBR(mvec_view_block)(mem1_,&vec1_,pad_pre_,pad_pre_+nvec_-1,&iflag_);
        ASSERT_EQ(0,iflag_);
      }
      if( _numberOfVectorsInitialized >= 2 )
      {
        SUBR(mvec_view_block)(mem2_,&vec2_,pad_pre_,pad_pre_+nvec_-1,&iflag_);
        ASSERT_EQ(0,iflag_);
      }
      if( _numberOfVectorsInitialized >= 3 )
      {
        SUBR(mvec_view_block)(mem3_,&vec3_,pad_pre_,pad_pre_+nvec_-1,&iflag_);
        ASSERT_EQ(0,iflag_);
      }
    }
    else
    {
      vec1_=mem1_;
      vec2_=mem2_;
      vec3_=mem3_;
    }
  }
}


  /*! Set up routine.
   */
virtual void SetUp()
{
  KernelTestWithMap<_Nglob>::SetUp();

  if (typeImplemented_ && !problemTooSmall_)
  {
    if( useViews_ )
    {
      if( _numberOfVectorsInitialized >= 1 )
      {
        SUBR(mvec_put_value)(mem1_,(_ST_)-101.,&iflag_);
        ASSERT_EQ(0,iflag_);
      }
      if( _numberOfVectorsInitialized >= 2 )
      {
        SUBR(mvec_put_value)(mem2_,(_ST_)-102.,&iflag_);
        ASSERT_EQ(0,iflag_);
      }
      if( _numberOfVectorsInitialized >= 3 )
      {
        SUBR(mvec_put_value)(mem3_,(_ST_)-103.,&iflag_);
        ASSERT_EQ(0,iflag_);
      }
    }

    if( _numberOfVectorsInitialized >= 1 )
    {
      SUBR(mvec_put_value)(vec1_,st::zero(),&iflag_);
      ASSERT_EQ(0,iflag_);
    }
    if( _numberOfVectorsInitialized >= 2 )
    {
      SUBR(mvec_put_value)(vec2_,st::zero(),&iflag_);
      ASSERT_EQ(0,iflag_);
    }
    if( _numberOfVectorsInitialized >= 3 )
    {
      SUBR(mvec_put_value)(vec3_,st::zero(),&iflag_);
      ASSERT_EQ(0,iflag_);
    }

    // vectors created with the same function should get the same stride (lda)
    phist_lidx lda;

    // extract raw data of complete memory block
    if( _numberOfVectorsInitialized >= 1 )
    {
      SUBR(mvec_extract_view)(mem1_,&mem1_vp_,&lda,&iflag_);
      ASSERT_EQ(0,iflag_);
      lda_=lda;
    }
    if( _numberOfVectorsInitialized >= 2 )
    {
      SUBR(mvec_extract_view)(mem2_,&mem2_vp_,&lda,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_EQ(lda_,lda);
    }
    if( _numberOfVectorsInitialized >= 3 )
    {
      SUBR(mvec_extract_view)(mem3_,&mem3_vp_,&lda,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_EQ(lda_,lda);
    }


    // extract raw data of viewed block
    phist_lidx lda2;
    if( _numberOfVectorsInitialized >= 1 )
    {
      SUBR(mvec_extract_view)(vec1_,&vec1_vp_,&lda2,&iflag_);
      ASSERT_EQ(0,iflag_);
      lda = lda2;
    }
    if( _numberOfVectorsInitialized >= 2 )
    {
      SUBR(mvec_extract_view)(vec2_,&vec2_vp_,&lda2,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_EQ(lda,lda2);
    }
    if( _numberOfVectorsInitialized >= 3 )
    {
      SUBR(mvec_extract_view)(vec3_,&vec3_vp_,&lda2,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_EQ(lda,lda2);
    }

    // lda might be irrelevant for small views and not set by the kernel library
#ifdef PHIST_MVECS_ROW_MAJOR
    if( nloc_ > 1 ) {
      ASSERT_EQ(lda_,lda2);
    } else {
      ASSERT_EQ(lda_,lda_);
    }
#else
    if( nvec_ > 1 ) {
      ASSERT_EQ(lda_,lda2);
    } else {
      ASSERT_EQ(lda_,lda_);
    }
#endif

    stride_=1;
    if (useViews_)
    {
      PHIST_SOUT(PHIST_DEBUG,"Setting up the views with pad_pre %d and pad_post %d (complete padding %d, lda %d)\n", 
          pad_pre_, pad_post_, nvecPadded_, lda_);
    }
  }
}


  /*! Clean up.
   */
virtual void TearDown() 
{
  if (typeImplemented_ && !problemTooSmall_)
  {
    // verify nobody touched the unviewed parts!
    if (useViews_)
    {
      // download memory from device 
      if( _numberOfVectorsInitialized >= 1 )
      {
        SUBR(mvec_from_device)((void*)mem1_,&iflag_);
      }
      if( _numberOfVectorsInitialized >= 2 )
      {
        SUBR(mvec_from_device)((void*)mem2_,&iflag_);
      }
      if( _numberOfVectorsInitialized >= 3 )
      {
        SUBR(mvec_from_device)((void*)mem3_,&iflag_);
      }

      // check pre padding is still the same
      if( _numberOfVectorsInitialized >= 1 )
      {
        ASSERT_REAL_EQ(mt::one(), ArrayEqual(mem1_vp_,nloc_,pad_pre_,lda_,stride_,(_ST_)-101.,vflag_));
      }
      if( _numberOfVectorsInitialized >= 2 )
      {
        ASSERT_REAL_EQ(mt::one(), ArrayEqual(mem2_vp_,nloc_,pad_pre_,lda_,stride_,(_ST_)-102.,vflag_));
      }
      if( _numberOfVectorsInitialized >= 3 )
      {
        ASSERT_REAL_EQ(mt::one(), ArrayEqual(mem3_vp_,nloc_,pad_pre_,lda_,stride_,(_ST_)-103.,vflag_));
      }

      // check post padding is still the same
#ifdef PHIST_MVECS_ROW_MAJOR
      if( _numberOfVectorsInitialized >= 1 )
      {
        ASSERT_REAL_EQ(mt::one(), ArrayEqual(mem1_vp_+pad_pre_+nvec_,nloc_,pad_post_,lda_,stride_,(_ST_)-101.,vflag_));
      }
      if( _numberOfVectorsInitialized >= 2 )
      {
        ASSERT_REAL_EQ(mt::one(), ArrayEqual(mem2_vp_+pad_pre_+nvec_,nloc_,pad_post_,lda_,stride_,(_ST_)-102.,vflag_));
      }
      if( _numberOfVectorsInitialized >= 3 )
      {
        ASSERT_REAL_EQ(mt::one(), ArrayEqual(mem3_vp_+pad_pre_+nvec_,nloc_,pad_post_,lda_,stride_,(_ST_)-103.,vflag_));
      }
#else
      if( _numberOfVectorsInitialized >= 1 )
      {
        ASSERT_REAL_EQ(mt::one(), ArrayEqual(mem1_vp_+(pad_pre_+nvec_)*lda_,nloc_,pad_post_,lda_,stride_,(_ST_)-101.,vflag_));
      }
      if( _numberOfVectorsInitialized >= 2 )
      {
        ASSERT_REAL_EQ(mt::one(), ArrayEqual(mem2_vp_+(pad_pre_+nvec_)*lda_,nloc_,pad_post_,lda_,stride_,(_ST_)-102.,vflag_));
      }
      if( _numberOfVectorsInitialized >= 3 )
      {
        ASSERT_REAL_EQ(mt::one(), ArrayEqual(mem3_vp_+(pad_pre_+nvec_)*lda_,nloc_,pad_post_,lda_,stride_,(_ST_)-103.,vflag_));
      }
#endif
    }
  }
  KernelTestWithMap<_Nglob>::TearDown();
}

static void TearDownTestCase()
{
  if (typeImplemented_)
  {
    if( _numberOfVectorsInitialized >= 1 )
    {
      SUBR(mvec_delete)(vec1_,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
    if( _numberOfVectorsInitialized >= 2 )
    {
      SUBR(mvec_delete)(vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
    if( _numberOfVectorsInitialized >= 3 )
    {
      SUBR(mvec_delete)(vec3_,&iflag_);
      ASSERT_EQ(0,iflag_);
    }

    // delete memory blocks if vecX are views
    if (useViews_)
    {
      if( _numberOfVectorsInitialized >= 1 )
      {
        SUBR(mvec_delete)(mem1_,&iflag_);
        ASSERT_EQ(0,iflag_);
      }
      if( _numberOfVectorsInitialized >= 2 )
      {
        SUBR(mvec_delete)(mem2_,&iflag_);
        ASSERT_EQ(0,iflag_);
      }
      if( _numberOfVectorsInitialized >= 3 )
      {
        SUBR(mvec_delete)(mem3_,&iflag_);
        ASSERT_EQ(0,iflag_);
      }
    }

    vec1_ = vec2_ = vec3_ = NULL;
    mem1_ = mem2_ = mem3_ = NULL;
  }

  KernelTestWithMap<_Nglob>::TearDownTestCase();
}


// tolerance for tests depending on the vector length
inline static MT releps(TYPE(const_mvec_ptr) V=NULL)
{
  if (V==NULL) return std::sqrt((MT)_Nglob)*mt::eps();
  int nvec,iflag;
  SUBR(mvec_num_vectors)(V,&nvec,&iflag);
  MT *nrms = new MT[nvec];
  MT max_nrm=0;
  SUBR(mvec_norm2)(V,nrms,&iflag);
  for (int i=0;i<nvec;i++)
  {
    max_nrm = std::max(max_nrm,nrms[i]);
  }
  if (max_nrm<4*mt::eps()) max_nrm=mt::sqrt((MT)_Nglob);
  delete[] nrms;
  return max_nrm*mt::eps();
}

//! in-place reduction operation on scalar data type (for testing with MPI)
static int global_sum(ST* value, int count, MPI_Comm mpi_comm)
{
  int iflag=0;
#ifdef PHIST_HAVE_MPI
        ST* gvalue = new ST[count];
        for (int i=0;i<count;i++) gvalue[i]=st::zero();
        iflag=MPI_Allreduce(value,gvalue,count,
                st::mpi_type(), MPI_SUM, mpi_comm);
        for (int i=0;i<count;i++) value[i]=gvalue[i];
        delete [] gvalue;
#endif
  return iflag;
}

#if defined(PHIST_HIGH_PRECISION_KERNELS) && defined(IS_DOUBLE) && !defined(IS_COMPLEX)

//! in-place reduction operation on scalar data type (for testing with MPI)
static int global_prec_sum(ST* value, ST* err, int count, MPI_Comm mpi_comm)
{
  std::cout << "global_prec_sum local value[0]: "<<value[0] << " + " << err[0]<<std::endl;
# ifndef PHIST_HAVE_MPI
      // not implemented
      return -99;
# endif
  int iflag=0;
  int nproc;
  MPI_Comm_size(mpi_comm,&nproc);
      ST gvalue[count*nproc];
      ST gerr[count*nproc];

        iflag=MPI_Allgather(value,count,st::mpi_type(),gvalue,count,st::mpi_type(),mpi_comm);
        if (iflag) return iflag;
        iflag=MPI_Allgather(err,count,st::mpi_type(),gerr,count,st::mpi_type(),mpi_comm);
        if (iflag) return iflag;
        
        for (int i=0; i<count; i++) 
        {
          value[i]=st::zero();
          err[i]=st::zero();
        }
        prec_reduction_k(nproc,count,gvalue,gerr,value,err);
  std::cout << "global_prec_sum global value[0]: "<<value[0] << " " << gvalue[0]<<std::endl;
  return iflag;
}

#endif
//! in-place reduction operation on scalar data type (for testing with MPI)
static int global_msum(MT* value, int count, MPI_Comm mpi_comm)
{
  int iflag=0;
#ifdef PHIST_HAVE_MPI
        ST* gvalue = new MT[count];        
        iflag=MPI_Allreduce(value,gvalue,count,
                mt::mpi_type(), MPI_SUM, mpi_comm);
        for (int i=0;i<count;i++) value[i]=gvalue[i];
        delete [] gvalue;
#endif
  return iflag;
}

  //! tests if each column of an mv is normalized in the B-norm
  static MT ColsAreBNormalized(const ST* V_vp, const ST* BV_vp, phist_lidx nloc, 
        phist_lidx ldV, phist_lidx ldBV, phist_lidx stride,
        MPI_Comm mpi_comm)
  {
    MT res=1.0;
    // see if all columns in vec2 have B-norm 1
    ST norms[nvec_],errs[nvec_];
    for (int j=0;j<nvec_;j++)
    {
#if defined(PHIST_HIGH_PRECISION_KERNELS) && defined(IS_DOUBLE) && !defined(IS_COMPLEX)
      ST sum=st::zero(),err=st::zero();
      for (int i=0;i<stride*nloc;i+=stride)
      {
        ST val=V_vp[VIDX(i,j,ldV)];
        ST Bval=BV_vp[VIDX(i,j,ldBV)];
        DOUBLE_4DOTADD(val,Bval,sum,err);
      }
      norms[j]=sum;
      errs[j]=err;
#else
      ST sum=st::zero();
      for (int i=0;i<stride*nloc;i+=stride)
      {
        ST val=V_vp[VIDX(i,j,ldV)];
        ST Bval=BV_vp[VIDX(i,j,ldBV)];
        sum+=st::conj(val)*Bval; 
      }
      norms[j]=sum;
#endif
    }
#if defined(PHIST_HIGH_PRECISION_KERNELS) && defined(IS_DOUBLE) && !defined(IS_COMPLEX)
    global_prec_sum(norms,errs,nvec_,mpi_comm);
    for (int j=0;j<nvec_;j++)
    {
      double sqrt_a,sqrt_aC,divsqrt_a,divsqrt_aC;
      DOUBLE_4SQRT_NEWTONRAPHSON_FMA(norms[j],errs[j],sqrt_a,sqrt_aC,divsqrt_a,divsqrt_aC);
      norms[j]=sqrt_a;
    }
#else
    global_sum(norms,nvec_,mpi_comm);
    for (int j=0;j<nvec_;j++) norms[j]=st::sqrt(norms[j]);
#endif
    res=ArrayEqual(norms,nvec_,1,nvec_,1,st::one());
    return res;
  }

  // check if vectors are mutually B-orthogonal, V'BV=I
  static MT ColsAreBOrthogonal(ST* V_vp, ST* BV_vp, phist_lidx nloc, 
        phist_lidx ldV, phist_lidx ldBV, phist_lidx stride,MPI_Comm mpi_comm) 
  {
    MT res=mt::one();
    int nsums=(nvec_*nvec_-nvec_)/2;
    if(nsums == 0)
    {
      return mt::one();
    }
    else
    {
      ST sums[nsums+1],errs[nsums+1];
      int k=0;
      for (int j1=0;j1<nvec_;j1++)
        for (int j2=j1+1;j2<nvec_;j2++)
        {
#if defined(PHIST_HIGH_PRECISION_KERNELS) && defined(IS_DOUBLE) && !defined(IS_COMPLEX)
          // TODO: we need high precision arithmetic for other data types, for instance in the ScalarTraits
          ST sum=st::zero(),err=st::zero();
          for (int i=0;i<stride*nloc;i+=stride)
          {
            double val1=V_vp[VIDX(i,j1,ldV)];
            double val2=BV_vp[VIDX(i,j2,ldBV)];
            DOUBLE_4DOTADD(val1,val2,sum,err);
          }
          errs[k]=err;
          sums[k++]=sum;
#else
          ST sum=st::zero();
          for (int i=0;i<stride*nloc;i+=stride)
          {
            ST val1=V_vp[VIDX(i,j1,ldV)];
            ST val2=BV_vp[VIDX(i,j2,ldBV)];
            sum+=val1*st::conj(val2);
          }
          errs[k]=st::zero();
          sums[k++]=sum;
#endif
        }
#if defined(PHIST_HIGH_PRECISION_KERNELS) && defined(IS_DOUBLE) && !defined(IS_COMPLEX)
      global_prec_sum(sums,errs,nsums,mpi_comm);
#else
      global_sum(sums,nsums,mpi_comm);
#endif
      res=ArrayEqual(sums,nsums,1,nsums,1,st::zero());
      return res;
    }
  }

  //! tests if each column of an mv is normalized in the 2-norm, v_i'*v_i=1
  static MT ColsAreNormalized(const ST* V_vp, phist_lidx nloc, 
        phist_lidx ldV, phist_lidx stride,
        MPI_Comm mpi_comm)
    {
      return ColsAreBNormalized(V_vp,V_vp,nloc,ldV,ldV,stride,mpi_comm);
    }

  // check if vectors are mutually orthogonal in the 2-norm, V'V=I
  static MT ColsAreOrthogonal(ST* V_vp, phist_lidx nloc, phist_lidx ldV, phist_lidx stride,MPI_Comm mpi_comm) 
  {
    return ColsAreBOrthogonal(V_vp,V_vp,nloc,ldV,ldV,stride,mpi_comm);
  }

  // check if vectors are mutually orthogonal after QR factorization
  static void PrintVector(int outlev, std::string label, 
        ST* vec_vp, phist_lidx nloc, phist_lidx lda, phist_lidx stride,MPI_Comm mpi_comm) 
  {
    if (outlev<PHIST_OUTLEV) return;
    int rank=0, np=1;
#ifdef PHIST_HAVE_MPI
    MPI_Comm_rank(mpi_comm,&rank);
    MPI_Comm_size(mpi_comm,&np);
#endif    
    std::ostringstream oss;
    if (rank==0)
    {
      oss << std::endl<<label <<"="<<std::endl;
      oss << "nproc  "<<np<<std::endl;
      oss << "nglob  "<<_Nglob<<std::endl;
      oss << "nloc   "<<nloc<<std::endl;
      oss << "lda    "<<lda<<std::endl;
      oss << "stride "<<stride<<std::endl;
#ifdef PHIST_MVECS_ROW_MAJOR
      oss << "row-major storage"<<std::endl;
#else
      oss << "col-major storage"<<std::endl;
#endif
    }
    if (_Nglob>1000)
    {
      if (rank==0)
      {
        oss << "(values not printed because vector too big)"<<std::endl;
      }
    }
    else
    {
      for (int i=0;i<stride*nloc;i+=stride)
      {
        for (int j=0;j<nvec_;j++)
        {
          oss << vec_vp[VIDX(i,j,lda)]<<"  ";
        }//j
        oss << std::endl;
      }//i
    }//if small enough print values
    PHIST_ORDERED_OUT(outlev,mpi_comm,oss.str().c_str());
    return;
  }
    
  static bool pointerUnchanged(TYPE(mvec_ptr) V, ST* expected_location, int expected_lda)
  { 
    int iflag;
    phist_lidx lda;
    ST* ptr;
    SUBR(mvec_extract_view)(V,&ptr,&lda,&iflag);
    return ( (iflag==0)&&(lda==expected_lda)&&(ptr==expected_location) );
  }

  //! if _useViews=true, these are larger memory blocks
  //! holding vecX_ as an inner view. Otherwise they are NULL.
private:
  static TYPE(mvec_ptr) mem1_, mem2_, mem3_;
  ST *mem1_vp_ = NULL, *mem2_vp_ = NULL, *mem3_vp_ = NULL;
  static int nvecPadded_, pad_pre_, pad_post_;

public:

  static TYPE(mvec_ptr) vec1_, vec2_, vec3_;
  ST *vec1_vp_ = NULL, *vec2_vp_ = NULL, *vec3_vp_ = NULL;
  static const int nvec_=_Nvec;
  static const int useViews_=_useViews;
  phist_lidx lda_, stride_;
};

template<phist_gidx _Nglob, int _Nvec, int _useViews, int _numberOfVectorsInitialized, int _multipleDefinitionCounter>
TYPE(mvec_ptr) KernelTestWithVectors<_ST_,_Nglob,_Nvec, _useViews, _numberOfVectorsInitialized, _multipleDefinitionCounter>::mem1_ = NULL;

template<phist_gidx _Nglob, int _Nvec, int _useViews, int _numberOfVectorsInitialized, int _multipleDefinitionCounter>
TYPE(mvec_ptr) KernelTestWithVectors<_ST_,_Nglob,_Nvec, _useViews, _numberOfVectorsInitialized, _multipleDefinitionCounter>::mem2_ = NULL;

template<phist_gidx _Nglob, int _Nvec, int _useViews, int _numberOfVectorsInitialized, int _multipleDefinitionCounter>
TYPE(mvec_ptr) KernelTestWithVectors<_ST_,_Nglob,_Nvec, _useViews, _numberOfVectorsInitialized, _multipleDefinitionCounter>::mem3_ = NULL;


template<phist_gidx _Nglob, int _Nvec, int _useViews, int _numberOfVectorsInitialized, int _multipleDefinitionCounter>
int KernelTestWithVectors<_ST_,_Nglob,_Nvec, _useViews, _numberOfVectorsInitialized, _multipleDefinitionCounter>::nvecPadded_;

template<phist_gidx _Nglob, int _Nvec, int _useViews, int _numberOfVectorsInitialized, int _multipleDefinitionCounter>
int KernelTestWithVectors<_ST_,_Nglob,_Nvec, _useViews, _numberOfVectorsInitialized, _multipleDefinitionCounter>::pad_pre_;

template<phist_gidx _Nglob, int _Nvec, int _useViews, int _numberOfVectorsInitialized, int _multipleDefinitionCounter>
int KernelTestWithVectors<_ST_,_Nglob,_Nvec, _useViews, _numberOfVectorsInitialized, _multipleDefinitionCounter>::pad_post_;


template<phist_gidx _Nglob, int _Nvec, int _useViews, int _numberOfVectorsInitialized, int _multipleDefinitionCounter>
TYPE(mvec_ptr) KernelTestWithVectors<_ST_,_Nglob,_Nvec, _useViews, _numberOfVectorsInitialized, _multipleDefinitionCounter>::vec1_ = NULL;

template<phist_gidx _Nglob, int _Nvec, int _useViews, int _numberOfVectorsInitialized, int _multipleDefinitionCounter>
TYPE(mvec_ptr) KernelTestWithVectors<_ST_,_Nglob,_Nvec, _useViews, _numberOfVectorsInitialized, _multipleDefinitionCounter>::vec2_ = NULL;

template<phist_gidx _Nglob, int _Nvec, int _useViews, int _numberOfVectorsInitialized, int _multipleDefinitionCounter>
TYPE(mvec_ptr) KernelTestWithVectors<_ST_,_Nglob,_Nvec, _useViews, _numberOfVectorsInitialized, _multipleDefinitionCounter>::vec3_ = NULL;


template<phist_gidx _Nglob, int _Nvec, int _useViews, int _numberOfVectorsInitialized, int _multipleDefinitionCounter>
const int KernelTestWithVectors<_ST_,_Nglob,_Nvec, _useViews, _numberOfVectorsInitialized, _multipleDefinitionCounter>::nvec_;

template<phist_gidx _Nglob, int _Nvec, int _useViews, int _numberOfVectorsInitialized, int _multipleDefinitionCounter>
const int KernelTestWithVectors<_ST_,_Nglob,_Nvec, _useViews, _numberOfVectorsInitialized, _multipleDefinitionCounter>::useViews_;

