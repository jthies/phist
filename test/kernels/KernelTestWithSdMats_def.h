
/*! Test fixture. */
template<int _Nrows, int _Ncols, int _useViews>
class KernelTestWithSdMats<_ST_,_Nrows,_Ncols, _useViews> : 
        public virtual KernelTestWithType< _ST_ >,
        public virtual KernelTest
  {

public:

  /*! Set up routine.
   */
virtual void SetUp()
{
  KernelTest::SetUp();
  KernelTestWithType< _ST_ >::SetUp();
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

    SUBR(sdMat_create)(&mem1_,nr_padded_,nc_padded_,this->comm_,&this->iflag_);
    ASSERT_EQ(0,this->iflag_);
    SUBR(sdMat_create)(&mem2_,nr_padded_,nc_padded_,this->comm_,&this->iflag_);
    ASSERT_EQ(0,this->iflag_);
    SUBR(sdMat_create)(&mem3_,nr_padded_,nc_padded_,this->comm_,&this->iflag_);
    ASSERT_EQ(0,this->iflag_);
    SUBR(sdMat_create)(&mem4_,nr_padded_,nc_padded_,this->comm_,&this->iflag_);
    ASSERT_EQ(0,this->iflag_);
    
    if (useViews_)
    {
      mat1_=NULL; mat2_=NULL; mat3_=NULL; mat4_=NULL;
      SUBR(sdMat_view_block)(mem1_,&mat1_,pad_nr_pre_,pad_nr_pre_+nrows_-1,
                             pad_nc_pre_,pad_nc_pre_+ncols_-1,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
      SUBR(sdMat_view_block)(mem2_,&mat2_,pad_nr_pre_,pad_nr_pre_+nrows_-1,
                             pad_nc_pre_,pad_nc_pre_+ncols_-1,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
      SUBR(sdMat_view_block)(mem3_,&mat3_,pad_nr_pre_,pad_nr_pre_+nrows_-1,
                             pad_nc_pre_,pad_nc_pre_+ncols_-1,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
      SUBR(sdMat_view_block)(mem4_,&mat4_,pad_nr_pre_,pad_nr_pre_+nrows_-1,
                             pad_nc_pre_,pad_nc_pre_+ncols_-1,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);

      // fill memory with checkable data
      SUBR(sdMat_put_value)(mem1_,(_ST_)-1001.,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
      SUBR(sdMat_put_value)(mem2_,(_ST_)-1002.,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
      SUBR(sdMat_put_value)(mem3_,(_ST_)-1003.,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
      SUBR(sdMat_put_value)(mem4_,(_ST_)-1004.,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);

      // zero usable block
      SUBR(sdMat_put_value)(mat1_,st::zero(),&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
      SUBR(sdMat_put_value)(mat2_,st::zero(),&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
      SUBR(sdMat_put_value)(mat3_,st::zero(),&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
      SUBR(sdMat_put_value)(mat4_,st::zero(),&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
    }
    else
    {
      mat1_=mem1_;
      mat2_=mem2_;
      mat3_=mem3_;
      mat4_=mem4_;
    }
    
    // get pointers to the whole memory block
    lidx_t lda;
    SUBR(sdMat_extract_view)(mem1_,&mem1_vp_,&lda,&this->iflag_);
    ASSERT_EQ(0,this->iflag_);
    this->m_lda_ = lda;
    SUBR(sdMat_extract_view)(mem2_,&mem2_vp_,&lda,&this->iflag_);
    ASSERT_EQ(0,this->iflag_);
    ASSERT_EQ(this->m_lda_,lda);
    SUBR(sdMat_extract_view)(mem3_,&mem3_vp_,&lda,&this->iflag_);
    ASSERT_EQ(0,this->iflag_);
    ASSERT_EQ(this->m_lda_,lda);
    SUBR(sdMat_extract_view)(mem4_,&mem4_vp_,&lda,&this->iflag_);
    ASSERT_EQ(0,this->iflag_);
    ASSERT_EQ(this->m_lda_,lda);

    // get raw access to the matrix storage
    lidx_t lda2;
    SUBR(sdMat_extract_view)(mat1_,&mat1_vp_,&lda,&this->iflag_);
    ASSERT_EQ(0,this->iflag_);
    lda2 = lda;
    SUBR(sdMat_extract_view)(mat2_,&mat2_vp_,&lda,&this->iflag_);
    ASSERT_EQ(0,this->iflag_);
    ASSERT_EQ(lda2,lda);
    SUBR(sdMat_extract_view)(mat3_,&mat3_vp_,&lda,&this->iflag_);
    ASSERT_EQ(0,this->iflag_);
    ASSERT_EQ(lda2,lda);
    SUBR(sdMat_extract_view)(mat4_,&mat4_vp_,&lda,&this->iflag_);
    ASSERT_EQ(0,this->iflag_);
    ASSERT_EQ(lda2,lda);

    // lda might be irrelevant for small views and thus not set by the kernel library
#ifdef PHIST_SDMATS_ROW_MAJOR
    if( nrows_ > 1 ) {
#else
    if( ncols_ > 1 ) {
#endif
      ASSERT_EQ(this->m_lda_,lda2);
    }
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

      SUBR(sdMat_delete)(mat1_,&this->iflag_);
      SUBR(sdMat_delete)(mat2_,&this->iflag_);
      SUBR(sdMat_delete)(mat3_,&this->iflag_);
      SUBR(sdMat_delete)(mat4_,&this->iflag_);    
    }
    SUBR(sdMat_delete)(mem4_,&this->iflag_);
    SUBR(sdMat_delete)(mem3_,&this->iflag_);
    SUBR(sdMat_delete)(mem2_,&this->iflag_);
    SUBR(sdMat_delete)(mem1_,&this->iflag_);
  }
  KernelTestWithType< _ST_ >::TearDown();
  KernelTest::TearDown();
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

  
  TYPE(sdMat_ptr) mem1_, mem2_, mem3_, mem4_;
  TYPE(sdMat_ptr) mat1_, mat2_, mat3_, mat4_;
  _ST_ *mat1_vp_, *mat2_vp_, *mat3_vp_, *mat4_vp_;
  _ST_ *mem1_vp_, *mem2_vp_, *mem3_vp_, *mem4_vp_;
  static const int nrows_=_Nrows;
  static const int ncols_=_Ncols;
  static const int useViews_=_useViews;
  int pad_nr_pre_, pad_nc_pre_, pad_nr_post_, pad_nc_post_;
  int nr_padded_, nc_padded_;
  lidx_t m_lda_;
};


template<int _Nrows, int _Ncols, int _useViews>
const int KernelTestWithSdMats<_ST_,_Nrows,_Ncols,_useViews>::nrows_;

template<int _Nrows, int _Ncols, int _useViews>
const int KernelTestWithSdMats<_ST_,_Nrows,_Ncols,_useViews>::ncols_;

template<int _Nrows, int _Ncols, int _useViews>
const int KernelTestWithSdMats<_ST_,_Nrows,_Ncols,_useViews>::useViews_;

