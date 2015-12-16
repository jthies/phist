// basic operations for some extended matrix/vector types

  // constructor - does not allocate memory
  TYPE(x_mvec)::TYPE(x_mvec)()
  {
    v_=NULL;
    vp_=NULL;
    vi_=NULL;
    vpi_=NULL;
    own_mvecs_=true;
  }

  // constructor that views given mvecs or takes ownership
  TYPE(x_mvec)::TYPE(x_mvec)(TYPE(mvec_ptr) v, TYPE(mvec_ptr) vi, int naug, bool take_ownership, int* iflag)
  {
    v_=v;
    vp_=NULL;
    vi_=vi;
    vpi_=NULL;
    int nvec;
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(v_,&nvec,iflag),*iflag);
    if (naug>0)
    {
      const_map_ptr_t map=NULL;
      PHIST_CHK_IERR(mvec_get_map)(v_,&map,iflag),*iflag);
      const_comm_ptr_t comm=NULL;
      PHIST_CHK_IERR(map_get_comm(map,&comm,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(sdMat_create)(&vp_,naug,nvec,comm,iflag),*iflag);
      if (rc)
      {
        PHIST_CHK_IERR(SUBR(sdMat_create)(&vpi_,naug,nvec,comm,iflag),*iflag);
      }
    }
    own_mvecs_=take_ownership;
  }
  
  // imaginary part is allocated only if rc=true, augmented part only if naug>0.
  TYPE(x_mvec)::allocate(const_map_ptr_t map, int nvec, int naug, bool rc, int* iflag)
  {
    *iflag=0;
    if (v_!=NULL) 
    {
      // either allocate has already been called, or
      // the constructor with given mvecs has been used.
      *iflag=PHIST_INVALID_INPUT;
      return;
    }
    PHIST_CHK_IERR(SUBR(mvec_create)(&v_,map,nvec,iflag),*iflag);
    if (rc)
    {
      PHIST_CHK_IERR(SUBR(mvec_create)(&vi_,map,nvec,iflag),*iflag);
    }
    if (naug>0)
    {
      const_comm_ptr_t comm=NULL;
      PHIST_CHK_IERR(map_get_comm(map,&comm,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(sdMat_create)(&vp_,naug,nvec,comm,iflag),*iflag);
      if (rc)
      {
        PHIST_CHK_IERR(SUBR(sdMat_create)(&vpi_,naug,nvec,comm,iflag),*iflag);
      }
    }
  }

  // destructor
  TYPE(x_mvec)::~TYPE(x_mvec)()
  {
    deallocate();
  }

void TYPE(x_mvec)::deallocate()
{
    int iflag;
    if (!own_mvecs_)
    {
      if (v_) PHIST_CHK_IERR(SUBR(mvec_delete)(v_,&iflag),iflag); v_=NULL;
      if (vi_) PHIST_CHK_IERR(SUBR(mvec_delete)(vi_,&iflag),iflag); vi_=NULL;
    }
    if (vp_) PHIST_CHK_IERR(SUBR(sdMat_delete)(vp_,&iflag),iflag); vp_=NULL;
    if (vpi_) PHIST_CHK_IERR(SUBR(sdMat_delete)(vpi_,&iflag),iflag); vpi=NULL;
}
  


//! return true if both vectors have an allocated imaginary part,
//! false if none of them has, and throw an exception if only one of them has.
inline bool rc_variant(TYPE(x_mvec) const* v1, TYPE(x_mvec) const* v2)
{
#ifdef IS_COMPLEX
  return false;
#endif
  if (v1->vi_==NULL && v2->vi_==NULL) return false;
  if (v1->vi_!=NULL && v2->vi_!=NULL) return true;
  throw "either both or none of the vectors must have an imaginary part!"
#endif
}

//! returns true if matrix and vectors are all 'complex in real arithmetic'
inline bool rc_variant(TYPE(const_x_sparseMat_ptr) A, TYPE(x_mvec) const* v1, TYPE(x_mvec) const* v2)
{
  bool rc=rc_variant(v1,v2);
  return rc && (A->sigma_i_!=NULL);
}

//! return true if both vectors are augmented by additional rows,
//! false if none of them has, and throw an exception if only one of them has.
//! We do not check if the number of augmented rows is different, in that case
//! some later phist call will return an error.
inline bool aug_variant(TYPE(x_mvec) const* v1, TYPE(x_mvec) const* v2)
{
  if (v1->vp_ ==NULL && v2->vp_ ==NULL
    &&v1->vpi_==NULL && v2->vpi_==NULL) return false;

  if (v1->vp_ !=NULL && v2->vp_ !=NULL)
  {
    if (rc_variant(v1,v2)
    {
      if (v1->vpi_!=NULL && v2->vpi_!=NULL) return true;
    }
    return true;
  }
  throw "either both or none of the vectors must be augmented!"
}

// returns true if both vectors and matrix are augmented by additional rows (rows and cols)
inline bool aug_variant(TYPE(const_x_sparseMat_ptr) A, TYPE(x_mvec) const* v1, TYPE(x_mvec) const* v2)
{
  bool rc=rc_variant(v1,v2);
  return rc && A->Vproj_!=NULL;
}

//!
void SUBR(x_mvec_add_mvec)(_ST_ alpha, TYPE(x_mvec) const* V,
                            _ST_ beta, TYPE(x_mvec)* W, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  bool rc =  rc_variant(V,W);
  bool aug= aug_variant(V,W);
  
  PHIST_CHK_IERR(SUBR(mvec_add_mvec)(alpha,V->v_,beta,W->v_,iflag),*iflag);
  if (aug) PHIST_CHK_IERR(SUBR(sdMat_add_sdMat)(alpha,V->vp_,beta,W->vp_,iflag),*iflag);
  if (rc)
  { 
    PHIST_CHK_IERR(SUBR(mvec_add_mvec)(alpha,V->vi_,beta,W->vi_,iflag),*iflag);
    if (aug) PHIST_CHK_IERR(SUBR(sdMat_add_sdMat)(alpha,V->vpi_,beta,W->vpi_,iflag),*iflag);
  }
  return;
}

//! missing kernel function, I don't want to introduce it to the interface at this point,
//! we could move it to kernels/common if it is useful in other places
//! the _add signifies that dot is *not* initialized in this function but the result is added
//! to whatever is in there already, this is because of our special application in this file.
void SUBR(my_sdMat_dot_sdMat_add)(TYPE(const_sdMat_ptr) A, TYPE(const_sdMat_ptr) B, _ST_* dots, int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_FCN(__FUNCTION__);
  ST *a,*b;
  lidx_t lda,ldb;
  int nr, nc,nrb,ncb;

  // if this has a significant latency we could do it asynchronously while
  // computing the dot products for the actual vectors in x_mvec_dot_mvec
  PHIST_CHK_IERR(SUBR(sdMat_from_device)(A,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_from_device)(B,iflag),*iflag);

  PHIST_CHK_IERR(SUBR(sdMat_extract_view)(A,&a,&lda,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_extract_view)(B,&b,&ldb,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(A,&nr,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(A,&nc,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(B,&nrb,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(B,&ncb,iflag),*iflag);

  if (nr!=nrb || nc!=ncb)
  {
    *iflag=PHIST_INVALID_INPUT;
    return;
  }
  for (int j=0; j<nc; j++)
  {
    // note: the caller must make sure dots is initialized
    for (int i=0; i<nr; i++)
    {
#ifdef PHIST_SDMATS_ROW_MAJOR
      dots[j]+=st::conj(a[i*lda+j])*b[i*lda+j];
#else
      dots[j]+=st::conj(a[j*lda+i])*b[j*lda+i];
#endif
    }
  }
}


//! dot product of two possibly complex vectors, in the complex case the result is found in dots,
//! in the case of real arithmetic but complex vectors, dotsi is filled with the imaginary parts.
void SUBR(x_mvec_dot_mvec)(TYPE(x_mvec_ptr) v, TYPE(x_mvec_ptr) w,
                            _ST_   *dots, _MT_* dotsi, int *iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_FCN(__FUNCTION__);

  *iflag=0;
  
  bool aug=aug_variant(V,W);
  bool rc=rc_variant(V,W);
  
  PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(v->v_,w->v_,dots,iflag),*iflag);
  if (aug)
  {
    PHIST_CHK_IERR(SUBR(sdMat_dot_sdMat_add)(v->vp_,w->vp_,dots,iflag),*iflag);    
  }

  if (rc)
  {
    ST tmp[nvec];
    PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(v->vi_,w->vi_,tmp,iflag),*iflag);
    if (aug)
    {
      PHIST_CHK_IERR(SUBR(sdMat_dot_sdMat_add)(v->vpi_,w->vpi_,tmp,iflag),*iflag);    
    }
    for (int j=0;j<nvec;j++)
    {
      dots[j]+=tmp[j];
    }

    if ( (v!=w) && dotsi!=NULL)
    {
      PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(v->v_,w->vi_,dotsi,iflag),*iflag);      
      PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(v->vi_,w->v_,tmp,iflag),*iflag);      
      if (aug)
      {
        PHIST_CHK_IERR(SUBR(sdMat_dot_sdMat_add)(v->vp_,w->vpi_,dotsi,iflag),*iflag);    
        PHIST_CHK_IERR(SUBR(sdMat_dot_sdMat_add)(v->vpi_,w->vp_,tmp,iflag),*iflag);    
      }
      for (int j=0;j<nvec;j++)
      {
        dotsi[j]-=tmp[j];
      }
    }
    else if (dotsi!=NULL)
    {
      for (int j=0;j<nvec;j++)
      {
        dotsi[j]=mt::zero();
      }
    }
  }
  return;
}


//! rc matrix is
//! A-sigma_r[i]          sigma_i[i]
//!  -sigma_i[i]        A-sigma_r[i]
void SUBR(x_sparseMat_times_mvec)(_ST_ alpha, TYPE(const_x_sparseMat_ptr) A, TYPE(x_mvec) const* X,
                       _ST_ beta, TYPE(x_mvec)* Y, int *iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_FCN(__FUNCTION__);

  // these calls also check the sanity of the input (i.e. throw an exception
  // if, for instance, the vectors are not augmented but the matrix is)
  bool rc=rc_variant(A,X,Y);
  bool aug=aug_variant(A,X,Y);
  
  int nvec;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(X->v_,&nvec,iflag),*iflag);
  
  ST shifts[nvec];
  for (int i=0;i<nvec;i++)
  {
    shifts[i]=-(ST)A->sigma_r[i];
#ifdef IS_COMPLEX
    shifts[i]-=A->sigma_i[i]*st::cmplx_I();
#endif
  }
  PHIST_CHK_IERR(SUBR(sparseMat_times_mvec_vadd_mvec)(alpha,A->A_,shifts,X->v_,beta,Y->v_,iflag),*iflag);

  if (rc)
  {
    //yi=alpha(A-srI)xi+beta*yi
    PHIST_CHK_IERR(SUBR(sparseMat_times_mvec_vadd_mvec)
      (alpha,A->A_,shifts,X->vi_,beta,Y->vi_,iflag),*iflag);
    //y+=alpha*si*xi
    for (int i=0;i<nvec;i++) shifts[i]=alpha*A->sigma_i_[i];
    PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(shifts, X->vi_,st::one(),Y->v_,iflag),*iflag);
    // yi-=alpha*si*xr
    for (int i=0;i<nvec;i++) shifts[i]=-alpha*sigma_i[i];
    PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(shifts,X->v_,st::one(),Y->vi_,iflag),*iflag);
  }
  
  if (aug)
  {
    PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(alpha,A->Vproj_,X->vp_,beta,Y->v_);
    PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(alpha,A->Vproj_,X->v_,beta,Y->vp_);

    if (rc)
    {
      //yi=V*xpi+beta*yi
      PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(alpha,A->Vproj_,X->vpi_,beta,Y->vi_);
      //ypi=V'xi+beta*ypi
      PHIST_CHK_IERR(SUBR(mvecT_times_mvec)
        (alpha,A->Vproj,X->vi_,beta,Y->vpi_,iflag),*iflag);
    }
  }
  return;
}
