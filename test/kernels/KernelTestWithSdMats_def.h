
/*! Test fixture. */
template<int _Nrows, int _Ncols, int _useViews, int _multipleDefinitionCounter>
class KernelTestWithSdMats<_ST_,_Nrows,_Ncols, _useViews, _multipleDefinitionCounter> : 
        public virtual TestWithType< _ST_ >,
        public virtual KernelTest
  {

public:

static void SetUpTestCase()
{
  KernelTest::SetUpTestCase();
  TestWithType::SetUpTestCase();

  if (typeImplemented_)
  {
    bool align_p2=(useViews_==2);
    pad_nr_pre_=0, pad_nc_pre_=0;
    pad_nr_post_=0, pad_nc_post_=0;
    nr_padded_=nrows_, nc_padded_=ncols_;
    if (useViews_)
    {
      pad_nr_pre_=align_p2?8:3;
      pad_nr_post_=1;
      pad_nc_pre_=align_p2?8:5;
      pad_nc_post_=1;
    // padding to power of 2
    int nr_pow_p2=(int)ceil(log((double)(nrows_+pad_nr_pre_+pad_nr_post_))/log(2.0))+1;
    int nc_pow_p2=(int)ceil(log((double)(ncols_+pad_nc_pre_+pad_nc_post_))/log(2.0))+1;
    int nr_padded_p2=(int)pow(2.0,(double)nr_pow_p2);
    int nc_padded_p2=(int)pow(2.0,(double)nc_pow_p2);
    nr_padded_=align_p2?nr_padded_p2:nrows_+pad_nr_pre_+pad_nr_post_;
    nc_padded_=align_p2?nc_padded_p2:ncols_+pad_nc_pre_+pad_nc_post_;
    pad_nr_post_ = nr_padded_-nrows_-pad_nr_pre_;
    pad_nc_post_ = nc_padded_-ncols_-pad_nc_pre_;
  }

    SUBR(sdMat_create)(&mem1_,nr_padded_,nc_padded_,comm_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(sdMat_create)(&mem2_,nr_padded_,nc_padded_,comm_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(sdMat_create)(&mem3_,nr_padded_,nc_padded_,comm_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(sdMat_create)(&mem4_,nr_padded_,nc_padded_,comm_,&iflag_);
    ASSERT_EQ(0,iflag_);
    
    if (useViews_)
    {
      mat1_=NULL; mat2_=NULL; mat3_=NULL; mat4_=NULL;
      SUBR(sdMat_view_block)(mem1_,&mat1_,pad_nr_pre_,pad_nr_pre_+nrows_-1,
                             pad_nc_pre_,pad_nc_pre_+ncols_-1,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_view_block)(mem2_,&mat2_,pad_nr_pre_,pad_nr_pre_+nrows_-1,
                             pad_nc_pre_,pad_nc_pre_+ncols_-1,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_view_block)(mem3_,&mat3_,pad_nr_pre_,pad_nr_pre_+nrows_-1,
                             pad_nc_pre_,pad_nc_pre_+ncols_-1,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_view_block)(mem4_,&mat4_,pad_nr_pre_,pad_nr_pre_+nrows_-1,
                             pad_nc_pre_,pad_nc_pre_+ncols_-1,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
    else
    {
      mat1_=mem1_;
      mat2_=mem2_;
      mat3_=mem3_;
      mat4_=mem4_;
    }
    
  }
}

  /*! Set up routine.
   */
virtual void SetUp()
{
  KernelTest::SetUp();

  if (typeImplemented_)
  {
    if (useViews_)
    {
      // fill memory with checkable data
      SUBR(sdMat_put_value)(mem1_,(_ST_)-1001.,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_put_value)(mem2_,(_ST_)-1002.,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_put_value)(mem3_,(_ST_)-1003.,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_put_value)(mem4_,(_ST_)-1004.,&iflag_);
      ASSERT_EQ(0,iflag_);
    }

    // zero usable block
    SUBR(sdMat_put_value)(mat1_,st::zero(),&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(sdMat_put_value)(mat2_,st::zero(),&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(sdMat_put_value)(mat3_,st::zero(),&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(sdMat_put_value)(mat4_,st::zero(),&iflag_);
    ASSERT_EQ(0,iflag_);

    // get pointers to the whole memory block
    lidx_t lda;
    SUBR(sdMat_extract_view)(mem1_,&mem1_vp_,&lda,&iflag_);
    ASSERT_EQ(0,iflag_);
    this->m_lda_ = lda;
    SUBR(sdMat_extract_view)(mem2_,&mem2_vp_,&lda,&iflag_);
    ASSERT_EQ(0,iflag_);
    ASSERT_EQ(this->m_lda_,lda);
    SUBR(sdMat_extract_view)(mem3_,&mem3_vp_,&lda,&iflag_);
    ASSERT_EQ(0,iflag_);
    ASSERT_EQ(this->m_lda_,lda);
    SUBR(sdMat_extract_view)(mem4_,&mem4_vp_,&lda,&iflag_);
    ASSERT_EQ(0,iflag_);
    ASSERT_EQ(this->m_lda_,lda);

    // get raw access to the matrix storage
    lidx_t lda2;
    SUBR(sdMat_extract_view)(mat1_,&mat1_vp_,&lda,&iflag_);
    ASSERT_EQ(0,iflag_);
    lda2 = lda;
    SUBR(sdMat_extract_view)(mat2_,&mat2_vp_,&lda,&iflag_);
    ASSERT_EQ(0,iflag_);
    ASSERT_EQ(lda2,lda);
    SUBR(sdMat_extract_view)(mat3_,&mat3_vp_,&lda,&iflag_);
    ASSERT_EQ(0,iflag_);
    ASSERT_EQ(lda2,lda);
    SUBR(sdMat_extract_view)(mat4_,&mat4_vp_,&lda,&iflag_);
    ASSERT_EQ(0,iflag_);
    ASSERT_EQ(lda2,lda);

    // lda might be irrelevant for small views and thus not set by the kernel library
#ifdef PHIST_SDMATS_ROW_MAJOR
    if( nrows_ > 1 ) {
      ASSERT_EQ(this->m_lda_,lda2);
    }
#else
    if( ncols_ > 1 ) {
      ASSERT_EQ(this->m_lda_,lda2);
    }
#endif
  }
}

  /*! Clean up.
   */
virtual void TearDown() 
{
  if (typeImplemented_)
  {
    if (useViews_)
    {
   
      // download memory from device 
      SUBR(sdMat_from_device)((void*)mem1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_from_device)((void*)mem2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_from_device)((void*)mem3_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_from_device)((void*)mem4_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // check pre col padding is still the same
      ASSERT_REAL_EQ(mt::one(), ArrayEqual(mem1_vp_,nrows_,pad_nc_pre_,this->m_lda_,1,(_ST_)-1001.,this->mflag_));
      ASSERT_REAL_EQ(mt::one(), ArrayEqual(mem2_vp_,nrows_,pad_nc_pre_,this->m_lda_,1,(_ST_)-1002.,this->mflag_));
      ASSERT_REAL_EQ(mt::one(), ArrayEqual(mem3_vp_,nrows_,pad_nc_pre_,this->m_lda_,1,(_ST_)-1003.,this->mflag_));
      ASSERT_REAL_EQ(mt::one(), ArrayEqual(mem4_vp_,nrows_,pad_nc_pre_,this->m_lda_,1,(_ST_)-1004.,this->mflag_));

      // check pre row padding is still the same
      ASSERT_REAL_EQ(mt::one(), ArrayEqual(mem1_vp_,pad_nr_pre_,ncols_,this->m_lda_,1,(_ST_)-1001.,this->mflag_));
      ASSERT_REAL_EQ(mt::one(), ArrayEqual(mem2_vp_,pad_nr_pre_,ncols_,this->m_lda_,1,(_ST_)-1002.,this->mflag_));
      ASSERT_REAL_EQ(mt::one(), ArrayEqual(mem3_vp_,pad_nr_pre_,ncols_,this->m_lda_,1,(_ST_)-1003.,this->mflag_));
      ASSERT_REAL_EQ(mt::one(), ArrayEqual(mem4_vp_,pad_nr_pre_,ncols_,this->m_lda_,1,(_ST_)-1004.,this->mflag_));

      // check post col padding is still the same
#ifdef PHIST_SDMATS_ROW_MAJOR
      ASSERT_REAL_EQ(mt::one(), ArrayEqual(mem1_vp_+pad_nc_pre_+ncols_,nrows_,pad_nc_post_,this->m_lda_,1,(_ST_)-1001.,this->mflag_));
      ASSERT_REAL_EQ(mt::one(), ArrayEqual(mem2_vp_+pad_nc_pre_+ncols_,nrows_,pad_nc_post_,this->m_lda_,1,(_ST_)-1002.,this->mflag_));
      ASSERT_REAL_EQ(mt::one(), ArrayEqual(mem3_vp_+pad_nc_pre_+ncols_,nrows_,pad_nc_post_,this->m_lda_,1,(_ST_)-1003.,this->mflag_));
      ASSERT_REAL_EQ(mt::one(), ArrayEqual(mem4_vp_+pad_nc_pre_+ncols_,nrows_,pad_nc_post_,this->m_lda_,1,(_ST_)-1004.,this->mflag_));
#else
      ASSERT_REAL_EQ(mt::one(), ArrayEqual(mem1_vp_+(pad_nc_pre_+ncols_)*this->m_lda_,nrows_,pad_nc_post_,this->m_lda_,1,(_ST_)-1001.,this->mflag_));
      ASSERT_REAL_EQ(mt::one(), ArrayEqual(mem2_vp_+(pad_nc_pre_+ncols_)*this->m_lda_,nrows_,pad_nc_post_,this->m_lda_,1,(_ST_)-1002.,this->mflag_));
      ASSERT_REAL_EQ(mt::one(), ArrayEqual(mem3_vp_+(pad_nc_pre_+ncols_)*this->m_lda_,nrows_,pad_nc_post_,this->m_lda_,1,(_ST_)-1003.,this->mflag_));
      ASSERT_REAL_EQ(mt::one(), ArrayEqual(mem4_vp_+(pad_nc_pre_+ncols_)*this->m_lda_,nrows_,pad_nc_post_,this->m_lda_,1,(_ST_)-1004.,this->mflag_));
#endif

      // check post row padding is still the same
#ifdef PHIST_SDMATS_ROW_MAJOR
      ASSERT_REAL_EQ(mt::one(), ArrayEqual(mem1_vp_+(pad_nr_pre_+nrows_)*this->m_lda_,pad_nr_post_,ncols_,this->m_lda_,1,(_ST_)-1001.,this->mflag_));
      ASSERT_REAL_EQ(mt::one(), ArrayEqual(mem2_vp_+(pad_nr_pre_+nrows_)*this->m_lda_,pad_nr_post_,ncols_,this->m_lda_,1,(_ST_)-1002.,this->mflag_));
      ASSERT_REAL_EQ(mt::one(), ArrayEqual(mem3_vp_+(pad_nr_pre_+nrows_)*this->m_lda_,pad_nr_post_,ncols_,this->m_lda_,1,(_ST_)-1003.,this->mflag_));
      ASSERT_REAL_EQ(mt::one(), ArrayEqual(mem4_vp_+(pad_nr_pre_+nrows_)*this->m_lda_,pad_nr_post_,ncols_,this->m_lda_,1,(_ST_)-1004.,this->mflag_));
#else
      ASSERT_REAL_EQ(mt::one(), ArrayEqual(mem1_vp_+pad_nr_pre_+nrows_,pad_nr_post_,ncols_,this->m_lda_,1,(_ST_)-1001.,this->mflag_));
      ASSERT_REAL_EQ(mt::one(), ArrayEqual(mem2_vp_+pad_nr_pre_+nrows_,pad_nr_post_,ncols_,this->m_lda_,1,(_ST_)-1002.,this->mflag_));
      ASSERT_REAL_EQ(mt::one(), ArrayEqual(mem3_vp_+pad_nr_pre_+nrows_,pad_nr_post_,ncols_,this->m_lda_,1,(_ST_)-1003.,this->mflag_));
      ASSERT_REAL_EQ(mt::one(), ArrayEqual(mem4_vp_+pad_nr_pre_+nrows_,pad_nr_post_,ncols_,this->m_lda_,1,(_ST_)-1004.,this->mflag_));
#endif
    }
  }

  mat1_vp_ = mat2_vp_ = mat3_vp_ = mat4_vp_ = NULL;
  mem1_vp_ = mem2_vp_ = mem3_vp_ = mem4_vp_ = NULL;

  KernelTest::TearDown();
}

static void TearDownTestCase()
{
  if (typeImplemented_)
  {
    if (useViews_)
    {
      SUBR(sdMat_delete)(mat1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_delete)(mat2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_delete)(mat3_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_delete)(mat4_,&iflag_);    
      ASSERT_EQ(0,iflag_);
    }
    SUBR(sdMat_delete)(mem4_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(sdMat_delete)(mem3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(sdMat_delete)(mem2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(sdMat_delete)(mem1_,&iflag_);
    ASSERT_EQ(0,iflag_);
  }

  mat1_ = mat2_ = mat3_ = mat4_ = NULL;
  mem1_ = mem2_ = mem3_ = mem4_ = NULL;

  KernelTest::TearDownTestCase();
}

static void PrintSdMat(std::ostream& os, std::string label, 
        ST* mat_vp, lidx_t lda, lidx_t stride,MPI_Comm mpi_comm)
{
    int rank=0, np=1;
#ifdef PHIST_HAVE_MPI
    MPI_Comm_rank(mpi_comm,&rank);
    MPI_Comm_size(mpi_comm,&np);
#endif    
    if (rank==0)
    {
      os << std::endl<<label <<"="<<std::endl;
      os << "nproc  "<<np<<std::endl;
      os << "nrows  "<<nrows_<<std::endl;
      os << "ncols   "<<ncols_<<std::endl;
      os << "lda    "<<lda<<std::endl;
      os << "stride "<<stride<<std::endl;
#ifdef PHIST_SDMATS_ROW_MAJOR
      os << "row-major storage"<<std::endl;
#else
      os << "col-major storage"<<std::endl;
#endif
    }
    for (int p=0;p<np;p++)
    {
      if (p==rank)
      {
        os << " @ rank "<<p<<" @"<<std::endl;
        for (int i=0;i<stride*nrows_;i+=stride)
        {
          for (int j=0;j<ncols_;j++)
          {
            os << mat_vp[MIDX(i,j,lda)]<<"  ";
          }//j
          os << std::endl;
        }//i
        os << std::flush;
      }//rank==p
#ifdef PHIST_HAVE_MPI
        MPI_Barrier(mpi_comm);
        MPI_Barrier(mpi_comm);
        MPI_Barrier(mpi_comm);
        MPI_Barrier(mpi_comm);
        MPI_Barrier(mpi_comm);
#endif
    }//p
    return;
}

  static bool pointerUnchanged(TYPE(sdMat_ptr) V, ST* expected_location, int expected_lda)
  {
    int iflag = 0;
    lidx_t lda = 0;
    ST* ptr = NULL;
    SUBR(sdMat_extract_view)(V,&ptr,&lda,&iflag);
    if( iflag!= PHIST_SUCCESS )
      return false;

    lidx_t n = 0;
#ifdef PHIST_SDMATS_ROW_MAJOR
    SUBR(sdMat_get_nrows)(V,&n,&iflag);
#else
    SUBR(sdMat_get_ncols)(V,&n,&iflag);
#endif
    if( iflag!= PHIST_SUCCESS )
      return false;

    if( ptr!=expected_location )
      return false;

    // epetra doesn't set lda appropriately for single cols/rows (as it is not needed)
    if( n > 1 && lda!=expected_lda )
      return false;

    return true;
  }

  /*! compare sdMats on several procs
   */
  static void sdMat_parallel_check(TYPE(const_sdMat_ptr) mat, int* iflag)
  {
    *iflag = 0;
    // TODO: use correct communicator
    int n,m;
    PHIST_CHK_IERR(SUBR(sdMat_from_device)((void*)mat,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(mat, &m, iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(mat, &n, iflag),*iflag);
    _ST_* buff = new _ST_[m*n];
    _ST_* mat_raw;
    lidx_t lda;
    PHIST_CHK_IERR(SUBR(sdMat_extract_view)((TYPE(sdMat_ptr))mat, &mat_raw, &lda, iflag),*iflag);
    // copy data to buffer
    for(int j = 0; j < n; j++)
      for(int i = 0; i < m; i++)
        buff[j*m+i] = mat_raw[MIDX(i,j,lda)];
    // broadcast
    PHIST_CHK_IERR(*iflag = MPI_Bcast(buff,m*n,::phist::ScalarTraits<_ST_>::mpi_type(),0,MPI_COMM_WORLD),*iflag);
    // check
    int error = 0;
    for(int j = 0; j < n; j++)
      for(int i = 0; i < m; i++)
        if( buff[j*m+i] != mat_raw[MIDX(i,j,lda)] )
          error = 1;
    int globError = 0;
    PHIST_CHK_IERR(*iflag = MPI_Allreduce(&error,&globError,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD),*iflag);

    delete[] buff;

    if( globError )
    {
#if PHIST_OUTLEV>=PHIST_DEBUG
      PHIST_CHK_IERR(SUBR(sdMat_print)(mat,iflag),*iflag);
#endif
      *iflag = -1;
      return;
    }
  }
  
  static TYPE(sdMat_ptr) mem1_, mem2_, mem3_, mem4_;
  static TYPE(sdMat_ptr) mat1_, mat2_, mat3_, mat4_;
  _ST_ *mat1_vp_, *mat2_vp_, *mat3_vp_, *mat4_vp_;
  _ST_ *mem1_vp_, *mem2_vp_, *mem3_vp_, *mem4_vp_;
  static const int nrows_=_Nrows;
  static const int ncols_=_Ncols;
  static const int useViews_=_useViews;
  static int pad_nr_pre_, pad_nc_pre_, pad_nr_post_, pad_nc_post_;
  static int nr_padded_, nc_padded_;
  lidx_t m_lda_;
};


template<int _Nrows, int _Ncols, int _useViews, int _multipleDefinitionCounter>
TYPE(sdMat_ptr) KernelTestWithSdMats<_ST_,_Nrows,_Ncols,_useViews,_multipleDefinitionCounter>::mem1_;

template<int _Nrows, int _Ncols, int _useViews, int _multipleDefinitionCounter>
TYPE(sdMat_ptr) KernelTestWithSdMats<_ST_,_Nrows,_Ncols,_useViews,_multipleDefinitionCounter>::mem2_;

template<int _Nrows, int _Ncols, int _useViews, int _multipleDefinitionCounter>
TYPE(sdMat_ptr) KernelTestWithSdMats<_ST_,_Nrows,_Ncols,_useViews,_multipleDefinitionCounter>::mem3_;

template<int _Nrows, int _Ncols, int _useViews, int _multipleDefinitionCounter>
TYPE(sdMat_ptr) KernelTestWithSdMats<_ST_,_Nrows,_Ncols,_useViews,_multipleDefinitionCounter>::mem4_;


template<int _Nrows, int _Ncols, int _useViews, int _multipleDefinitionCounter>
TYPE(sdMat_ptr) KernelTestWithSdMats<_ST_,_Nrows,_Ncols,_useViews,_multipleDefinitionCounter>::mat1_;

template<int _Nrows, int _Ncols, int _useViews, int _multipleDefinitionCounter>
TYPE(sdMat_ptr) KernelTestWithSdMats<_ST_,_Nrows,_Ncols,_useViews,_multipleDefinitionCounter>::mat2_;

template<int _Nrows, int _Ncols, int _useViews, int _multipleDefinitionCounter>
TYPE(sdMat_ptr) KernelTestWithSdMats<_ST_,_Nrows,_Ncols,_useViews,_multipleDefinitionCounter>::mat3_;

template<int _Nrows, int _Ncols, int _useViews, int _multipleDefinitionCounter>
TYPE(sdMat_ptr) KernelTestWithSdMats<_ST_,_Nrows,_Ncols,_useViews,_multipleDefinitionCounter>::mat4_;


template<int _Nrows, int _Ncols, int _useViews, int _multipleDefinitionCounter>
const int KernelTestWithSdMats<_ST_,_Nrows,_Ncols,_useViews,_multipleDefinitionCounter>::nrows_;

template<int _Nrows, int _Ncols, int _useViews,int _multipleDefinitionCounter>
const int KernelTestWithSdMats<_ST_,_Nrows,_Ncols,_useViews,_multipleDefinitionCounter>::ncols_;

template<int _Nrows, int _Ncols, int _useViews, int _multipleDefinitionCounter>
const int KernelTestWithSdMats<_ST_,_Nrows,_Ncols,_useViews,_multipleDefinitionCounter>::useViews_;


template<int _Nrows, int _Ncols, int _useViews, int _multipleDefinitionCounter>
int KernelTestWithSdMats<_ST_,_Nrows,_Ncols,_useViews,_multipleDefinitionCounter>::pad_nr_pre_;

template<int _Nrows, int _Ncols, int _useViews, int _multipleDefinitionCounter>
int KernelTestWithSdMats<_ST_,_Nrows,_Ncols,_useViews,_multipleDefinitionCounter>::pad_nc_pre_;

template<int _Nrows, int _Ncols, int _useViews, int _multipleDefinitionCounter>
int KernelTestWithSdMats<_ST_,_Nrows,_Ncols,_useViews,_multipleDefinitionCounter>::pad_nr_post_;

template<int _Nrows, int _Ncols, int _useViews, int _multipleDefinitionCounter>
int KernelTestWithSdMats<_ST_,_Nrows,_Ncols,_useViews,_multipleDefinitionCounter>::pad_nc_post_;

template<int _Nrows, int _Ncols, int _useViews, int _multipleDefinitionCounter>
int KernelTestWithSdMats<_ST_,_Nrows,_Ncols,_useViews,_multipleDefinitionCounter>::nr_padded_;

template<int _Nrows, int _Ncols, int _useViews, int _multipleDefinitionCounter>
int KernelTestWithSdMats<_ST_,_Nrows,_Ncols,_useViews,_multipleDefinitionCounter>::nc_padded_;
