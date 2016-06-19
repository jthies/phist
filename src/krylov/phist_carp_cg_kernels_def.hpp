
#ifndef PHIST_KERNEL_LIB_BUILTIN
# ifndef VIEW_VDAT
# define VIEW_VDAT
# endif
#endif

// basic operations for some extended matrix/vector types

  // constructor - does not allocate memory
  TYPE(x_mvec)::TYPE(x_mvec)()
  {
    vdat_=NULL;
    v_=NULL;
    vp_=NULL;
    vi_=NULL;
    vpi_=NULL;
  }

  // constructor that allocates memory and copies given vector(s)
  TYPE(x_mvec)::TYPE(x_mvec)(TYPE(mvec_ptr) v, TYPE(mvec_ptr) vi, int naug, int* iflag)
  {

    vdat_=NULL;
    v_=NULL;
    vi_=NULL;
    vp_=NULL;
    vpi_=NULL;
    nvec_=-1;

    bool rc=vi!=NULL;
    phist_const_map_ptr map=NULL;
    PHIST_CHK_IERR(SUBR(mvec_get_map)(v,&map,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(v,&nvec_,iflag),*iflag);

    PHIST_CHK_IERR(allocate(map,nvec_,naug,rc,iflag),*iflag);

    PHIST_CHK_IERR(SUBR(mvec_set_block)(v_,v,0,nvec_-1,iflag),*iflag);
    if (rc)
    {
      PHIST_CHK_IERR(SUBR(mvec_set_block)(vi_,vi,0,nvec_-1,iflag),*iflag);
    }

    if (naug>0)
    {
      phist_const_comm_ptr comm=NULL;
      PHIST_CHK_IERR(phist_map_get_comm(map,&comm,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(sdMat_create)(&vp_,naug,nvec_,comm,iflag),*iflag);
      if (vi!=NULL)
      {
        PHIST_CHK_IERR(SUBR(sdMat_create)(&vpi_,naug,nvec_,comm,iflag),*iflag);
      }
    }
  }
  
  // imaginary part is allocated only if rc=true, augmented part only if naug>0.
  void TYPE(x_mvec)::allocate(phist_const_map_ptr map, int nvec, int naug, bool rc, int* iflag)
  {
    *iflag=0;
    if (vdat_!=NULL||v_!=NULL||vi_!=NULL) 
    {
      // either allocate has already been called, or
      // the constructor with given mvecs has been used.
      *iflag=PHIST_INVALID_INPUT;
      return;
    }
    nvec_=nvec;
#ifdef VIEW_VDAT
    // let v and vi view the first and last columns of vdat, this way the communication
    // in the CARP kernel can be done in one step
    int actual_nvec=rc?2*nvec:nvec;
    PHIST_CHK_IERR(SUBR(mvec_create)(&vdat_,map,actual_nvec,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(mvec_view_block)(vdat_,&v_,0,nvec_-1,iflag),*iflag);
    if (rc)
    {
      PHIST_CHK_IERR(SUBR(mvec_view_block)(vdat_,&vi_,nvec_,actual_nvec-1,iflag),*iflag);
    }
#else
    // the builtin implementation doesn't support this data layout yet, so allocate separate vectors
    vdat_=NULL;
    PHIST_CHK_IERR(SUBR(mvec_create)(&v_,map,nvec_,iflag),*iflag);
    if (rc)
    {
      PHIST_CHK_IERR(SUBR(mvec_create)(&vi_,map,nvec_,iflag),*iflag);
    }
#endif
    if (naug>0)
    {
      phist_const_comm_ptr comm=NULL;
      PHIST_CHK_IERR(phist_map_get_comm(map,&comm,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(sdMat_create)(&vp_,naug,nvec_,comm,iflag),*iflag);
      if (rc)
      {
        PHIST_CHK_IERR(SUBR(sdMat_create)(&vpi_,naug,nvec_,comm,iflag),*iflag);
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
    if (v_) PHIST_CHK_IERR(SUBR(mvec_delete)(v_,&iflag),iflag); v_=NULL;
    if (vi_) PHIST_CHK_IERR(SUBR(mvec_delete)(vi_,&iflag),iflag); vi_=NULL;
    if (vdat_) PHIST_CHK_IERR(SUBR(mvec_delete)(vdat_,&iflag),iflag); vdat_=NULL;

    if (vp_) PHIST_CHK_IERR(SUBR(sdMat_delete)(vp_,&iflag),iflag); vp_=NULL;
    if (vpi_) PHIST_CHK_IERR(SUBR(sdMat_delete)(vpi_,&iflag),iflag); vpi_=NULL;
}
  
//!
void SUBR(x_mvec_add_mvec)(_ST_ alpha, TYPE(x_mvec) const* V,
                            _ST_ beta, TYPE(x_mvec)* W, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  bool rc =  rc_variant(V,W);
  bool aug= aug_variant(V,W);

#ifdef VIEW_VDAT  
  PHIST_CHK_IERR(SUBR(mvec_add_mvec)(alpha,V->vdat_,beta,W->vdat_,iflag),*iflag);
#else
  PHIST_CHK_IERR(SUBR(mvec_add_mvec)(alpha,V->v_,beta,W->v_,iflag),*iflag);
  if (rc)
  {
    PHIST_CHK_IERR(SUBR(mvec_add_mvec)(alpha,V->vi_,beta,W->vi_,iflag),*iflag);
  }
#endif
  if (aug)
  { 
    PHIST_CHK_IERR(SUBR(sdMat_add_sdMat)(alpha,V->vp_,beta,W->vp_,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_add_sdMat)(alpha,V->vpi_,beta,W->vpi_,iflag),*iflag);
  }
  return;
}

//! missing kernel function, B(:,i)=alpha[i]*A(:,i)+beta*B(:,i)
//! note: can't make A const here because we need a "from_device" call to do it on the host
//! if this is a GPU process.
void SUBR(sdMat_vadd_sdMat)(_ST_ const alpha[], TYPE(sdMat_ptr) A,
                            _ST_       beta,    TYPE(sdMat_ptr)       B, int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_FCN(__FUNCTION__);
  
  ST *a,*b;
  phist_lidx lda,ldb;
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
    for (int i=0; i<nr; i++)
    {
#ifdef PHIST_SDMATS_ROW_MAJOR
      b[i*lda+j]=alpha[j]*a[i*lda+j] + beta*b[i*lda+j];
#else
      b[j*lda+i]=alpha[j]*a[j*lda+i] + beta*b[j*lda+i];
#endif
    }
  }
  PHIST_CHK_IERR(SUBR(sdMat_to_device)(B,iflag),*iflag);
}

//! missing kernel function, A(:,i)*=alpha[i]
void SUBR(sdMat_vscale)(TYPE(sdMat_ptr) A, _ST_ const alpha[], int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_FCN(__FUNCTION__);
  
  ST *a;
  phist_lidx lda;
  int nr, nc;

  // if this has a significant latency we could do it asynchronously while
  // computing the dot products for the actual vectors in x_mvec_dot_mvec
  PHIST_CHK_IERR(SUBR(sdMat_from_device)(A,iflag),*iflag);

  PHIST_CHK_IERR(SUBR(sdMat_extract_view)(A,&a,&lda,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(A,&nr,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(A,&nc,iflag),*iflag);

  for (int j=0; j<nc; j++)
  {
    for (int i=0; i<nr; i++)
    {
#ifdef PHIST_SDMATS_ROW_MAJOR
      a[i*lda+j]*=alpha[j];
#else
      a[j*lda+i]*=alpha[j];
#endif
    }
  }
  PHIST_CHK_IERR(SUBR(sdMat_to_device)(A,iflag),*iflag);
}

//! missing kernel function, I don't want to introduce it to the interface at this point,
//! we could move it to kernels/common if it is useful in other places
//! the _add signifies that dot is *not* initialized in this function but the result is added
//! to whatever is in there already, this is because of our special application in this file.
void SUBR(sdMat_dot_sdMat_add)(TYPE(sdMat_ptr) A, TYPE(sdMat_ptr) B, _ST_* dots, int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_FCN(__FUNCTION__);
  ST *a,*b;
  phist_lidx lda,ldb;
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
void SUBR(x_mvec_dot_mvec)(TYPE(x_mvec)* v, TYPE(x_mvec)* w,
                            _ST_   *dots, _MT_* dotsi, int *iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_FCN(__FUNCTION__);

  *iflag=0;
  
  bool aug=aug_variant(v,w);
  bool rc=rc_variant(v,w);

  int nvec=v->nvec_;
  int actual_nvec=rc?v->nvec_*2: nvec;
  ST tmp[actual_nvec];
  
#ifdef VIEW_VDAT
  PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(v->vdat_,w->vdat_,tmp,iflag),*iflag);
#else
  PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(v->v_,w->v_,tmp,iflag),*iflag);
  if (rc)
  {
    PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(v->vi_,w->vi_,tmp+nvec,iflag),*iflag);
  }
#endif
  if (aug)
  {
    PHIST_CHK_IERR(SUBR(sdMat_dot_sdMat_add)(v->vp_,w->vp_,tmp,iflag),*iflag);    
  }

  for (int j=0;j<nvec;j++)
  {
    dots[j]=tmp[j];
  }

#ifndef IS_COMPLEX
  // avoid uninitialized imaginary part
  if (dotsi!=NULL)
  {
    for (int j=0;j<nvec;j++)
    {
      dotsi[j]=mt::zero();
    }
  }
  if (rc)
  {
    if (aug)
    {
      PHIST_CHK_IERR(SUBR(sdMat_dot_sdMat_add)(v->vpi_,w->vpi_,tmp+nvec,iflag),*iflag);    
    }
    for (int j=0;j<nvec;j++)
    {
      dots[j]+=tmp[nvec+j];
    }

    if ( (v!=w) && dotsi!=NULL)
    {
      PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(v->v_,w->vi_,tmp,iflag),*iflag);      
      PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(v->vi_,w->v_,tmp+nvec,iflag),*iflag);      
      if (aug)
      {
        PHIST_CHK_IERR(SUBR(sdMat_dot_sdMat_add)(v->vp_,w->vpi_,tmp,iflag),*iflag);    
        PHIST_CHK_IERR(SUBR(sdMat_dot_sdMat_add)(v->vpi_,w->vp_,tmp+nvec,iflag),*iflag);    
      }
      for (int j=0;j<nvec;j++)
      {
        dotsi[j]=tmp[j]-tmp[nvec+j];
      }
    }
  }
#endif
  return;
}


//! rc matrix is
//! A-sigma_r[i]          sigma_i[i]
//!  -sigma_i[i]        A-sigma_r[i]
void SUBR(x_sparseMat_times_mvec)(_ST_ alpha, TYPE(x_sparseMat) const* A, TYPE(x_mvec) const* X,
                       _ST_ beta, TYPE(x_mvec)* Y, int *iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_FCN(__FUNCTION__);

  // these calls also check the sanity of the input (i.e. throw an exception
  // if, for instance, the vectors are not augmented but the matrix is)
  bool rc=rc_variant(A,X,Y);
  bool aug=aug_variant(A,X,Y);
  
  int nvec=X->nvec_;
  int actual_nvec=rc?2*nvec:nvec;
  
  ST shifts[actual_nvec];
  for (int i=0;i<nvec;i++)
  {
    shifts[i]=-(ST)A->sigma_r_[i];
    if (rc) shifts[nvec+i]=shifts[i];
#ifdef IS_COMPLEX
    shifts[i]-=A->sigma_i_[i]*st::cmplx_I();
#endif
  }
#ifdef VIEW_VDAT
  PHIST_CHK_IERR(SUBR(sparseMat_times_mvec_vadd_mvec)(alpha,A->A_,shifts,X->vdat_,beta,Y->vdat_,iflag),*iflag);
#else
  PHIST_CHK_IERR(SUBR(sparseMat_times_mvec_vadd_mvec)(alpha,A->A_,shifts,X->v_,beta,Y->v_,iflag),*iflag);
  if (rc)
  {
    PHIST_CHK_IERR(SUBR(sparseMat_times_mvec_vadd_mvec)(alpha,A->A_,shifts+nvec,X->vi_,beta,Y->vi_,iflag),*iflag);
  }
#endif
  // note: this is not optimal in row-major storage with the interleaved storage vdat=[v vi], if it is expensive
  // we could add a kernel function doing this directly on vdat
  if (rc)
  {
    //y+=alpha*si*xi
    for (int i=0;i<nvec;i++) shifts[i]=alpha*A->sigma_i_[i];
    PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(shifts, X->vi_,st::one(),Y->v_,iflag),*iflag);
    // yi-=alpha*si*xr
    for (int i=0;i<nvec;i++) shifts[i]=-alpha*A->sigma_i_[i];
    PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(shifts,X->v_,st::one(),Y->vi_,iflag),*iflag);
  }
  
  if (aug)
  {
    //y+=V*xp
    PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(alpha,A->Vproj_,X->vp_,st::one(),Y->v_,iflag),*iflag);
    //yp=V'x+beta*yp
    PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(alpha,A->Vproj_,X->v_,beta,Y->vp_,iflag),*iflag);

    if (rc)
    {
      //yi+=V*xpi
      PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(alpha,A->Vproj_,X->vpi_,st::one(),Y->vi_,iflag),*iflag);
      //ypi=V'xi+beta*ypi
      PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(alpha,A->Vproj_,X->vi_,beta,Y->vpi_,iflag),*iflag);
    }
  }
  return;
}

void SUBR(x_carp_sweep)(TYPE(x_sparseMat) const* A,TYPE(const_mvec_ptr) b,TYPE(x_mvec)* x,
        void* carp_data, _MT_ const omega[],int *iflag)
{
#include "phist_std_typedefs.hpp"
  //TODO Jonas: we should have an interface that makes use of vdat_=[v_ vi_] in the RC case,
  //     that may save 50% messages (but require the same communication volume)

  // double carp sweep in place, updates r=dkswp(sI-A,omega,r)
#ifndef IS_COMPLEX
  if (rc_variant(A,x,x))
  {
    if (aug_variant(A,x,x))
    {
       PHIST_CHK_IERR(SUBR(carp_sweep_aug_rc)(A->A_, A->sigma_r_, A->sigma_i_,A->Vproj_,b,x->v_,x->vi_,
            x->vp_,x->vpi_,carp_data,omega,iflag),*iflag);
    }
    else
    {
      PHIST_CHK_IERR(SUBR(carp_sweep_rc)(A->A_, A->sigma_r_, A->sigma_i_,b,x->v_,x->vi_,
            carp_data, omega, iflag),*iflag);
    }
  }
  else
#endif
  {
    _ST_* sigma=NULL;
#ifdef IS_COMPLEX
    sigma=new _ST_[x->nvec_];
    for (int i=0; i<x->nvec_;i++) sigma[i]=A->sigma_r_[i]+A->sigma_i_[i]*st::cmplx_I();
#else
    sigma=A->sigma_r_;
#endif
    if (aug_variant(A,x,x))
    {
      PHIST_CHK_IERR(SUBR(carp_sweep_aug)(A->A_, sigma,A->Vproj_,b,x->v_,
            x->vp_,carp_data,omega,iflag),*iflag);
    }
    else
    {
      PHIST_CHK_IERR(SUBR(carp_sweep)(A->A_,sigma,b,x->v_,
            carp_data, omega, iflag),*iflag);
      
    }
#ifdef IS_COMPLEX
    delete [] sigma;
#endif  
  }
}


//! Y = alpha X + beta Y
//!
//! yr = alpha_r xr + beta_r yr 
//!    - alpha_i xi 
//! yi = alpha_r xi + beta_r xi 
//!    + alpha_i xr 
void SUBR(x_mvec_vadd_mvec)(_ST_ const alpha[], _MT_ const alpha_i[], TYPE(x_mvec) const* X, _ST_ beta, TYPE(x_mvec)* Y, int* iflag)

{
#include "phist_std_typedefs.hpp"
  *iflag=0;
  bool rc = rc_variant(X,Y);
  bool aug=aug_variant(X,Y);
  int nvec=X->nvec_,actual_nvec;

  actual_nvec=rc?2*nvec:nvec;
  ST tmp_alpha[actual_nvec];
  
  for (int i=0; i<nvec; i++)
  {
    tmp_alpha[i]=alpha[i];
    if (rc) tmp_alpha[nvec+i]=alpha[i];
  }
#ifdef VIEW_VDAT
  PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(tmp_alpha,X->vdat_,beta,Y->vdat_,iflag),*iflag);
#else
  PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(tmp_alpha,X->v_,beta,Y->v_,iflag),*iflag);
  if (rc)
  {
    PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(tmp_alpha+nvec,X->vi_,beta,Y->vi_,iflag),*iflag);
  }
#endif
  if (aug)
  {
    PHIST_CHK_IERR(SUBR(sdMat_vadd_sdMat)(alpha,X->vp_,beta,Y->vp_,iflag),*iflag);
  }
#ifndef IS_COMPLEX
  if (rc)
  {
    for (int i=0; i<actual_nvec; i++)
    {
      tmp_alpha[i]=-(_ST_)alpha_i[i];
    }
    
    // y-=alpha_i*xi
    PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(tmp_alpha,X->vi_,st::one(),Y->v_,iflag),*iflag);
    // yi+=alpha_i*x
    PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(alpha_i,X->v_,st::one(),Y->vi_,iflag),*iflag);
    
    if (aug)
    {
      PHIST_CHK_IERR(SUBR(sdMat_vadd_sdMat)(alpha,X->vpi_,beta,Y->vpi_,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(sdMat_vadd_sdMat)(tmp_alpha,X->vpi_,st::one(),Y->vp_,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(sdMat_vadd_sdMat)(alpha_i,X->vp_,st::one(),Y->vpi_,iflag),*iflag);
    }
  }
#endif
  return;
}

//! scale columns i of v by real scalar alpha[i]
void SUBR(x_mvec_vscale)(TYPE(x_mvec)* v, _MT_ const alpha[], int* iflag)
{
  *iflag=0;
  bool rc=rc_variant(v,v);
  bool aug=aug_variant(v,v);

  int nvec=v->nvec_;
  _ST_ tmp_alpha[2*nvec];
  for (int i=0;i<nvec; i++) 
  {
    tmp_alpha[i]=(_ST_)alpha[i];
    tmp_alpha[nvec+i]=(_ST_)alpha[i];
  }
#ifdef VIEW_VDAT
  PHIST_CHK_IERR(SUBR(mvec_vscale)(v->vdat_,tmp_alpha,iflag),*iflag);
#else
  PHIST_CHK_IERR(SUBR(mvec_vscale)(v->v_,tmp_alpha,iflag),*iflag);
  if (rc)
  {
    PHIST_CHK_IERR(SUBR(mvec_vscale)(v->vi_,tmp_alpha+nvec,iflag),*iflag);
  }
#endif
  if (aug)
  {
    PHIST_CHK_IERR(SUBR(sdMat_vscale)(v->vp_,tmp_alpha,iflag),*iflag);    
    if (rc)
    {
      PHIST_CHK_IERR(SUBR(sdMat_vscale)(v->vpi_,tmp_alpha,iflag),*iflag);    
    }
  }
}
