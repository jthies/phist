/*! \file magma/kernels_def.hpp
 * included by magma/kernels.cpp
 * \author "Melven Roehrig-Zoellner <Melven.Roehrig-Zoellner@DLR.de>
 * \author "Jonas Thies <Jonas.Thies@DLR.de>
 *
*/

#ifdef MAGMA
#undef MAGMA
#undef MAGMABLAS
#undef MAGMA_HELPER
#undef MAGMA_HELPER2
#endif
#define MAGMABLAS(s) MAGMA_HELPER2(magmablas_,SPREFIX(s))
#define MAGMA(s) MAGMA_HELPER2(magma_,SPREFIX(s))
#define MAGMA_HELPER2(p,s) MAGMA_HELPER(p,s)
#define MAGMA_HELPER(p,s) p ## s

#ifdef MAGMA_ST
#undef MAGMA_ST
#undef MAGMAM
#endif
#ifdef IS_COMPLEX
# ifdef IS_DOUBLE
#  define MAGMA_ST magmaDoubleComplex
#  define MAGMAM(s) magma_dz ## s
# else
#  define MAGMA_ST magmaFloatComplex
#  define MAGMAM(s) magma_sc ## s
# endif
#else
#define MAGMA_ST _ST_
#define MAGMAM(s) MAGMA(s)
#endif

extern "C" void SUBR(type_avail)(int *iflag)
{
  *iflag=0;
}

extern "C" void SUBR(sparseMat_read_mm)(TYPE(sparseMat_ptr)* A, const_comm_ptr_t vcomm,
        const char* filename,int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMat_read_bin)(TYPE(sparseMat_ptr)* A, const_comm_ptr_t vcomm,
const char* filename,int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMat_read_hb)(TYPE(sparseMat_ptr)* A, const_comm_ptr_t vcomm,
const char* filename,int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMat_get_row_map)(TYPE(const_sparseMat_ptr) A, const_map_ptr_t* map, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMat_get_col_map)(TYPE(const_sparseMat_ptr) A, const_map_ptr_t* map, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMat_get_domain_map)(TYPE(const_sparseMat_ptr) A, const_map_ptr_t* map, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMat_get_range_map)(TYPE(const_sparseMat_ptr) A, const_map_ptr_t* map, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_create)(TYPE(mvec_ptr)* vV, 
    const_map_ptr_t map, lidx_t nvec, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t*,V,vV,*iflag);

  *V = new Traits<_ST_>::mvec_t;
  (*V)->map = map;
  PHIST_CHK_IERR(phist_map_get_local_length(map,&(*V)->n,iflag),*iflag);
  (*V)->nvec = nvec;
  (*V)->stride = (*V)->n;
  (*V)->is_view = false;
  MAGMA(malloc)((MAGMA_ST**)&(*V)->gpuData,(*V)->nvec*(*V)->n);
  PHIST_CHK_IERR(*iflag = (*V)->gpuData ? 0 : -1, *iflag);
  MAGMA(malloc_pinned)((MAGMA_ST**)&(*V)->cpuData,(*V)->nvec*(*V)->n);
  PHIST_CHK_IERR(*iflag = (*V)->cpuData ? 0 : -1, *iflag);
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_create_view)(TYPE(mvec_ptr)* vV, const_map_ptr_t map, 
    _ST_* values, lidx_t lda, int nvec,
    int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sdMat_create)(TYPE(sdMat_ptr)* vM, 
    int nrows, int ncols, const_comm_ptr_t comm, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t*,M,vM,*iflag);

  *M = new Traits<_ST_>::sdMat_t;
  (*M)->nrows = nrows;
  (*M)->ncols = ncols;
  (*M)->stride = nrows;
  (*M)->is_view = false;
  MAGMA(malloc)((MAGMA_ST**)&(*M)->gpuData,nrows*ncols);
  PHIST_CHK_IERR(*iflag = (*M)->gpuData ? 0 : -1, *iflag);
  MAGMA(malloc_cpu)((MAGMA_ST**)&(*M)->cpuData,nrows*ncols);
  PHIST_CHK_IERR(*iflag = (*M)->cpuData ? 0 : -1, *iflag);
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMat_create_view)(TYPE(sdMat_ptr)* vM, const_comm_ptr_t comm,
        _ST_* values, lidx_t lda, int nrows, int ncols,
        int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
  /*
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t*,M,vM,*iflag);

  (*M)->nrows = nrows;
  (*M)->ncols = ncols;
  (*M)->stride = lda;
  (*M)->is_view = true;
  (*M)->data = values;
  */
}
                  

extern "C" void SUBR(mvec_get_map)(TYPE(const_mvec_ptr) vV, const_map_ptr_t* map, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t,V,vV,*iflag);
  *map = V->map;
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_num_vectors)(TYPE(const_mvec_ptr) vV, int* nvec, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t,V,vV,*iflag);
  *nvec = V->nvec;
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMat_get_nrows)(TYPE(const_sdMat_ptr) vM, int* nrows, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t,M,vM,*iflag);
  *nrows = M->nrows;
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMat_get_ncols)(TYPE(const_sdMat_ptr) vM, int* ncols, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t,M,vM,*iflag);
  *ncols = M->ncols;
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_extract_view)(TYPE(mvec_ptr) vV, _ST_** val, lidx_t* lda, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t,V,vV,*iflag);
  *val = V->cpuData;
  *lda = V->stride;
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMat_extract_view)(TYPE(sdMat_ptr) vM, _ST_** val, lidx_t* lda, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t,M,vM,*iflag);
  *val = M->cpuData;
  *lda = M->stride;
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_to_mvec)(TYPE(const_mvec_ptr) vV, TYPE(mvec_ptr) vW, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t,W,vW,*iflag);

  MAGMABLAS(lacpy)(MagmaFull,W->n,W->nvec,
      (const MAGMA_ST*)V->gpuData,V->stride,
      (MAGMA_ST*)W->gpuData,W->stride);
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_view_block)(TYPE(mvec_ptr) vV,
    TYPE(mvec_ptr)* vVblock,
    int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);

  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t*,Vblock,vVblock,*iflag);

  if( (*Vblock) == NULL )
    *Vblock = new Traits<_ST_>::mvec_t;
  (*Vblock)->stride = V->stride;
  (*Vblock)->map = V->map;
  (*Vblock)->n = V->n;
  (*Vblock)->nvec = jmax-jmin+1;
  (*Vblock)->is_view = true;
  (*Vblock)->gpuData = V->gpuData + jmin*V->stride;
  (*Vblock)->cpuData = V->cpuData + jmin*V->stride;
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_get_block)(TYPE(const_mvec_ptr) vV,
    TYPE(mvec_ptr) vVblock,
    int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);

  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t,Vblock,vVblock,*iflag);

  MAGMABLAS(lacpy)(MagmaFull,Vblock->n,Vblock->nvec,(const MAGMA_ST*)V->gpuData+jmin*V->stride,V->stride,
      (MAGMA_ST*)Vblock->gpuData,Vblock->stride);
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_set_block)(TYPE(mvec_ptr) vV,
    TYPE(const_mvec_ptr) vVblock,
    int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);

  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t,Vblock,vVblock,*iflag);

  MAGMABLAS(lacpy)(MagmaFull,Vblock->n,Vblock->nvec,(const MAGMA_ST*)Vblock->gpuData,Vblock->stride,
      (MAGMA_ST*)V->gpuData+jmin*V->stride,V->stride);
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMat_view_block)(TYPE(mvec_ptr) vM, 
    TYPE(mvec_ptr)* vMblock,
    int imin, int imax, int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t,M,vM,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t*,Mblock,vMblock,*iflag);

  if( (*Mblock) == NULL )
    *Mblock = new Traits<_ST_>::sdMat_t;
  (*Mblock)->stride = M->stride;
  (*Mblock)->nrows = imax-imin+1;
  (*Mblock)->ncols = jmax-jmin+1;
  (*Mblock)->is_view = true;
  (*Mblock)->gpuData = M->gpuData + jmin*M->stride + imin;
  (*Mblock)->cpuData = M->cpuData + jmin*M->stride + imin;
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMat_get_block)(TYPE(const_mvec_ptr) vM, 
    TYPE(mvec_ptr) vMblock,
    int imin, int imax, int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t,M,vM,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t,Mblock,vMblock,*iflag);

  MAGMABLAS(lacpy)(MagmaFull,Mblock->nrows,Mblock->ncols,(const MAGMA_ST*)M->gpuData+jmin*M->stride+imin,M->stride,
      (MAGMA_ST*)Mblock->gpuData,Mblock->stride);
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMat_set_block)(TYPE(sdMat_ptr) vM, 
    TYPE(const_sdMat_ptr) vMblock,
    int imin, int imax, int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t,M,vM,*iflag);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t,Mblock,vMblock,*iflag);

  MAGMABLAS(lacpy)(MagmaFull,Mblock->nrows,Mblock->ncols,(const MAGMA_ST*)Mblock->gpuData,Mblock->stride,
      (MAGMA_ST*)M->gpuData+jmin*M->stride+imin,M->stride);
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sparseMat_delete)(TYPE(sparseMat_ptr) A, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_delete)(TYPE(mvec_ptr) vV, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t,V,vV,*iflag);
  if( !V->is_view )
  {
    magma_free(V->gpuData);
    magma_free_pinned(V->cpuData);
  }
  delete(V);
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMat_delete)(TYPE(sdMat_ptr) vM, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t,M,vM,*iflag);
  if( !M->is_view )
  {
    magma_free(M->gpuData);
    magma_free_cpu(M->cpuData);
  }
  delete(M);
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_put_value)(TYPE(mvec_ptr) vV, _ST_ value, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t,V,vV,*iflag);
  
  MAGMABLAS(laset)(MagmaFull,V->n,V->nvec,
      *(reinterpret_cast<const MAGMA_ST*>(&value)),
      *(reinterpret_cast<const MAGMA_ST*>(&value)),
      (MAGMA_ST*)V->gpuData,V->stride);
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_put_func)(TYPE(mvec_ptr) vV,
        phist_mvec_elemFunc funPtr,void* last_arg, int *iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t,V,vV,*iflag);
  for(int j = 0; j < V->nvec; j++)
    for(int i = 0; i < V->n; i++)
      PHIST_CHK_IERR(*iflag=funPtr(i,j,(void*)&(V->cpuData[j*V->stride+i]),last_arg),*iflag);

  MAGMA(setmatrix)(V->n,V->nvec,
      (const MAGMA_ST*)V->cpuData,V->stride,
      (MAGMA_ST*)V->gpuData,V->stride);
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMat_put_value)(TYPE(mvec_ptr) vM, _ST_ value, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t,M,vM,*iflag);

  MAGMABLAS(laset)(MagmaFull,M->nrows,M->ncols,
      *(reinterpret_cast<const MAGMA_ST*>(&value)),
      *(reinterpret_cast<const MAGMA_ST*>(&value)),
      (MAGMA_ST*)M->gpuData,M->stride);
  *iflag = PHIST_SUCCESS;
}

#ifndef PHIST_BUILTIN_RNG
extern "C" void SUBR(mvec_random)(TYPE(mvec_ptr) vV, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t,V,vV,*iflag);
  for(int j = 0; j < V->nvec; j++)
    for(int i = 0; i < V->n; i++)
      V->cpuData[j*V->stride+i] = st::rand();

  PHIST_CHK_IERR(SUBR(mvec_to_device)(vV,iflag),*iflag);
  *iflag = PHIST_SUCCESS;
}
#endif
extern "C" void SUBR(mvec_print)(TYPE(const_mvec_ptr) vV, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t,V,vV,*iflag);
  MAGMA(print)(V->n,V->nvec,(MAGMA_ST*)V->cpuData,V->stride);
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMat_print)(TYPE(const_sdMat_ptr) vM, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t,M,vM,*iflag);
  MAGMA(print)(M->nrows,M->ncols,(MAGMA_ST*)M->cpuData,M->stride);
  *iflag = PHIST_SUCCESS;
}

#ifndef PHIST_BUILTIN_RNG
extern "C" void SUBR(sdMat_random)(TYPE(sdMat_ptr) vM, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t,M,vM,*iflag);
  for(int j = 0; j < M->ncols; j++)
    for(int i = 0; i < M->nrows; i++)
      M->cpuData[j*M->stride+i] = st::prand();

  PHIST_CHK_IERR(SUBR(sdMat_to_device)(vM,iflag),*iflag);
  *iflag = PHIST_SUCCESS;
}
#endif

extern "C" void SUBR(sdMat_identity)(TYPE(sdMat_ptr) vM, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t,M,vM,*iflag);

  _ST_ one = st::one();
  _ST_ zero = st::zero();

  MAGMABLAS(laset)(MagmaFull,M->nrows,M->ncols,
      *(reinterpret_cast<const MAGMA_ST*>(&zero)),
      *(reinterpret_cast<const MAGMA_ST*>(&one)),
      (MAGMA_ST*)M->gpuData,M->stride);
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_norm2)(TYPE(const_mvec_ptr) vV,
    _MT_* vnrm, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t,V,vV,*iflag);

  for(int j = 0; j < V->nvec; j++)
    vnrm[j] = MAGMAM(nrm2)(V->n,(const MAGMA_ST*)V->gpuData+j*V->stride,1);
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_normalize)(TYPE(mvec_ptr) vV,
    _MT_* vnrm, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t,V,vV,*iflag);

  for(int j = 0; j < V->nvec; j++)
  {
    vnrm[j] = MAGMAM(nrm2)(V->n,(const MAGMA_ST*)V->gpuData+j*V->stride,1);
    _ST_ scal = vnrm[j] > mt::zero() ? st::one()/vnrm[j] : st::zero();
    MAGMA(scal)(V->n,
        *(reinterpret_cast<const MAGMA_ST*>(&scal)),
        (MAGMA_ST*)V->gpuData+j*V->stride,1);
  }
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_scale)(TYPE(mvec_ptr) vV, 
    _ST_ scalar, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t,V,vV,*iflag);

  MAGMA(scal)(V->n*V->nvec,
      *(reinterpret_cast<const MAGMA_ST*>(&scalar)),
      (MAGMA_ST*)V->gpuData,1);
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_vscale)(TYPE(mvec_ptr) vV, 
    const _ST_* scalar, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t,V,vV,*iflag);

  for(int j = 0; j < V->nvec; j++)
  {
    MAGMA(scal)(V->n,
        *(reinterpret_cast<const MAGMA_ST*>(&scalar[j])),
        (MAGMA_ST*)V->gpuData+j*V->stride,1);
  }
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_add_mvec)(_ST_ alpha, TYPE(const_mvec_ptr) vV,
                                    _ST_ beta,  TYPE(mvec_ptr)       vW,
                                    int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t,W,vW,*iflag);

  MAGMA(scal)(W->nvec*W->n,
      *(reinterpret_cast<const MAGMA_ST*>(&beta)),
      (MAGMA_ST*)W->gpuData,1);
  MAGMA(axpy)(W->nvec*W->n,
      *(reinterpret_cast<const MAGMA_ST*>(&alpha)),
      (const MAGMA_ST*)V->gpuData,1,
      (MAGMA_ST*)W->gpuData,1);
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_vadd_mvec)(const _ST_ alpha[], TYPE(const_mvec_ptr) vV,
                                    _ST_ beta,  TYPE(mvec_ptr)       vW,
                                    int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t,W,vW,*iflag);

  for(int j = 0; j < W->nvec; j++)
  {
    MAGMA(scal)(W->n,
        *(reinterpret_cast<const MAGMA_ST*>(&beta)),
        (MAGMA_ST*)W->gpuData+j*W->stride,1);
    MAGMA(axpy)(W->n,
      *(reinterpret_cast<const MAGMA_ST*>(&alpha[j])),
      (const MAGMA_ST*)V->gpuData+j*V->stride,1,
      (MAGMA_ST*)W->gpuData+j*W->stride,1);
  }
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMat_add_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) vA,
                                      _ST_ beta,  TYPE(sdMat_ptr)       vB, 
                                      int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t,A,vA,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t,B,vB,*iflag);

  PHIST_CHK_IERR(SUBR(sdMat_from_device)(const_cast<TYPE(sdMat_ptr)>(vA),iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_from_device)(vB,iflag),*iflag);
  for(int j = 0; j < B->ncols; j++)
    for(int i = 0; i < B->nrows; i++)
      B->cpuData[j*B->stride+i] = beta*B->cpuData[j*B->stride+i] + alpha*A->cpuData[j*A->stride+i];
  PHIST_CHK_IERR(SUBR(sdMat_to_device)(vB,iflag),*iflag);
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMatT_add_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) vA,
    _ST_ beta,  TYPE(sdMat_ptr)       vB, 
    int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t,A,vA,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t,B,vB,*iflag);

  PHIST_CHK_IERR(SUBR(sdMat_from_device)(const_cast<TYPE(sdMat_ptr)>(vA),iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_from_device)(vB,iflag),*iflag);
  for(int j = 0; j < B->ncols; j++)
    for(int i = 0; i < B->nrows; i++)
      B->cpuData[j*B->stride+i] = beta*B->cpuData[j*B->stride+i] + alpha*st::conj(A->cpuData[i*A->stride+j]);
  PHIST_CHK_IERR(SUBR(sdMat_to_device)(vB,iflag),*iflag);
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sparseMat_times_mvec_communicate)(TYPE(const_sparseMat_ptr) vA, TYPE(const_mvec_ptr) vx, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag = 0;
}

extern "C" void SUBR(sparseMat_times_mvec)(_ST_ alpha, TYPE(const_sparseMat_ptr) A, 
    TYPE(const_mvec_ptr) x, _ST_ beta, TYPE(mvec_ptr) y, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMatT_times_mvec)(_ST_ alpha, TYPE(const_sparseMat_ptr) A, 
    TYPE(const_mvec_ptr) x, _ST_ beta, TYPE(mvec_ptr) y, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

//! y[i]=alpha*(A*x[i]+shifts[i]*x[i]) + beta*y[i]
extern "C" void SUBR(sparseMat_times_mvec_vadd_mvec)(_ST_ alpha, TYPE(const_sparseMat_ptr) A,
        const _ST_ shifts[], TYPE(const_mvec_ptr) x, _ST_ beta, TYPE(mvec_ptr) y, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_dot_mvec)(TYPE(const_mvec_ptr) vV, 
                                    TYPE(const_mvec_ptr) vW, 
                                    _ST_* s, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t,W,vW,*iflag);

  for(int j = 0; j < V->nvec; j++)
  {
#ifdef IS_COMPLEX
    MAGMA_ST tmp = MAGMA(dotc)(V->n,(const MAGMA_ST*)V->gpuData+j*V->stride,1,(const MAGMA_ST*)W->gpuData+j*W->stride,1);
    s[j] = *(reinterpret_cast<const _ST_*>(&tmp));
#else
    s[j] = MAGMA(dot)(V->n,(const MAGMA_ST*)V->gpuData+j*V->stride,1,(const MAGMA_ST*)W->gpuData+j*W->stride,1);
#endif
  }
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_times_sdMat)(_ST_ alpha, TYPE(const_mvec_ptr) vV, 
    TYPE(const_sdMat_ptr) vM, 
    _ST_ beta, TYPE(mvec_ptr) vW, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t,M,vM,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t,W,vW,*iflag);

  MAGMA(gemm)(MagmaNoTrans,MagmaNoTrans,W->n,W->nvec,V->nvec,
      *(reinterpret_cast<const MAGMA_ST*>(&alpha)),
      (const MAGMA_ST*)V->gpuData,V->stride,
      (const MAGMA_ST*)M->gpuData,M->stride,
      *(reinterpret_cast<const MAGMA_ST*>(&beta)),
      (MAGMA_ST*)W->gpuData,W->stride);
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMat_times_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) vA, 
    TYPE(const_sdMat_ptr) vB, 
    _ST_ beta, TYPE(sdMat_ptr) vC, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t,A,vA,*iflag);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t,B,vB,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t,C,vC,*iflag);

  MAGMA(gemm)(MagmaNoTrans,MagmaNoTrans,C->nrows,C->ncols,A->ncols,
      *(reinterpret_cast<const MAGMA_ST*>(&alpha)),
      (const MAGMA_ST*)A->gpuData,A->stride,
      (const MAGMA_ST*)B->gpuData,B->stride,
      *(reinterpret_cast<const MAGMA_ST*>(&beta)),
      (MAGMA_ST*)C->gpuData,C->stride);
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMatT_times_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) vA, 
    TYPE(const_sdMat_ptr) vB, 
    _ST_ beta, TYPE(sdMat_ptr) vC, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t,A,vA,*iflag);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t,B,vB,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t,C,vC,*iflag);

#ifdef IS_COMPLEX
  magma_trans_t transA = MagmaConjTrans;
#else
  magma_trans_t transA = MagmaTrans;
#endif

  MAGMA(gemm)(transA,MagmaNoTrans,C->nrows,C->ncols,A->nrows,
      *(reinterpret_cast<const MAGMA_ST*>(&alpha)),
      (const MAGMA_ST*)A->gpuData,A->stride,
      (const MAGMA_ST*)B->gpuData,B->stride,
      *(reinterpret_cast<const MAGMA_ST*>(&beta)),
      (MAGMA_ST*)C->gpuData,C->stride);
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMat_times_sdMatT)(_ST_ alpha, TYPE(const_sdMat_ptr) vA,
    TYPE(const_sdMat_ptr) vB,
    _ST_ beta, TYPE(sdMat_ptr) vC, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t,A,vA,*iflag);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t,B,vB,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t,C,vC,*iflag);

#ifdef IS_COMPLEX
  magma_trans_t transB = MagmaConjTrans;
#else
  magma_trans_t transB = MagmaTrans;
#endif

  MAGMA(gemm)(MagmaNoTrans,transB,C->nrows,C->ncols,A->nrows,
      *(reinterpret_cast<const MAGMA_ST*>(&alpha)),
      (const MAGMA_ST*)A->gpuData,A->stride,
      (const MAGMA_ST*)B->gpuData,B->stride,
      *(reinterpret_cast<const MAGMA_ST*>(&beta)),
      (MAGMA_ST*)C->gpuData,C->stride);
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvecT_times_mvec)(_ST_ alpha, TYPE(const_mvec_ptr) vV, 
                                       TYPE(const_mvec_ptr) vW, 
                                       _ST_ beta, TYPE(sdMat_ptr) vM, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t,W,vW,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t,M,vM,*iflag);

#ifdef IS_COMPLEX
  magma_trans_t transV = MagmaConjTrans;
#else
  magma_trans_t transV = MagmaTrans;
#endif

  MAGMA(gemm)(transV,MagmaNoTrans,M->nrows,M->ncols,V->n,
      *(reinterpret_cast<const MAGMA_ST*>(&alpha)),
      (const MAGMA_ST*)V->gpuData,V->stride,
      (const MAGMA_ST*)W->gpuData,W->stride,
      *(reinterpret_cast<const MAGMA_ST*>(&beta)),
      (MAGMA_ST*)M->gpuData,M->stride);
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_QR)(TYPE(mvec_ptr) vV, TYPE(sdMat_ptr) vR, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t,R,vR,*iflag);

  // allocate required work buffer
  MAGMA_ST* gpuWork = NULL;
  MAGMA_ST* cpuWork = NULL;
  MAGMA(malloc)(&gpuWork,3*V->nvec*V->nvec);
  MAGMA(malloc_pinned)(&cpuWork,3*V->nvec*V->nvec);

  // check zero rank
  bool zeroRank = true;
  for(int j = 0; j < V->nvec; j++)
  {
    _MT_ nrm = MAGMAM(nrm2)(V->n,(const MAGMA_ST*)V->gpuData+j*V->stride,1);
    if( nrm > 1000*mt::eps() )
      zeroRank = false;
  }

  // run magmas QR
  if( zeroRank )
  {
    PHIST_CHK_IERR(SUBR(sdMat_put_value)(R,st::zero(),iflag),*iflag);
  }
  else
  {
    PHIST_CHK_IERR(MAGMA(gegqr_gpu)(3,V->n,V->nvec,(MAGMA_ST*)V->gpuData,V->stride,gpuWork,cpuWork,iflag),*iflag);
    // resulting R is in cpuWork
    MAGMA(setmatrix)(R->nrows,R->ncols,cpuWork,R->nrows,(MAGMA_ST*)R->gpuData,R->stride);
  }
  //MAGMA(print_gpu)(R->nrows,R->ncols,(const MAGMA_ST*)R->gpuData,R->stride);
  // get rank
  int rank = 0;
  if( !zeroRank )
  {
    _MT_ ref = st::abs(*(reinterpret_cast<const _ST_*>(&cpuWork[0])));
    for(int i = 0; i < V->nvec; i++)
      if( st::abs(*(reinterpret_cast<const _ST_*>(&cpuWork[i*V->nvec+i]))) > 1000*mt::eps()*ref )
        rank++;
  }

  // randomize null space
  if( rank < V->nvec )
  {
    for(int j = rank; j < V->nvec; j++)
      for(int i = 0; i < V->n; i++)
        V->cpuData[j*V->stride+i] = st::rand();
    MAGMA(setmatrix)(V->n,V->nvec-rank,
        (const MAGMA_ST*)V->cpuData+rank*V->stride,V->stride,
        (MAGMA_ST*)V->gpuData+rank*V->stride,V->stride);

    PHIST_CHK_IERR(MAGMA(gegqr_gpu)(3,V->n,V->nvec,(MAGMA_ST*)V->gpuData,V->stride,gpuWork,cpuWork,iflag),*iflag);
  }

  // free work buffer
  magma_free(gpuWork);
  magma_free_pinned(cpuWork);

  // return rank
  *iflag = V->nvec-rank;
}

extern "C" void SUBR(mvec_gather_mvecs)(TYPE(mvec_ptr) V, TYPE(const_mvec_ptr) W[], int nblocks, int *iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_MARK_AS_EXPERIMENTAL(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_scatter_mvecs)(TYPE(const_mvec_ptr) V, TYPE(mvec_ptr) W[], int nblocks, int *iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_MARK_AS_EXPERIMENTAL(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

// NOTE: see the description of sparseMat_read_mm on how we treat input flags for this function
extern "C" void SUBR(sparseMat_create_fromRowFunc)(TYPE(sparseMat_ptr) *A, const_comm_ptr_t vcomm,
        gidx_t nrows, gidx_t ncols, lidx_t maxnne,
                phist_sparseMat_rowFunc rowFunPtr, void* last_arg,
                int *iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

void SUBR(mvec_to_device)(TYPE(mvec_ptr) vV, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t,V,vV,*iflag);

  MAGMA(setmatrix)(V->n,V->nvec,
      (const MAGMA_ST*)V->cpuData,V->stride,
      (MAGMA_ST*)V->gpuData,V->stride);
  *iflag = PHIST_SUCCESS;
}

void SUBR(mvec_from_device)(TYPE(mvec_ptr) vV, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t,V,vV,*iflag);

  MAGMA(getmatrix)(V->n,V->nvec,
      (const MAGMA_ST*)V->gpuData,V->stride,
      (MAGMA_ST*)V->cpuData,V->stride);
  *iflag = PHIST_SUCCESS;
}

void SUBR(sdMat_to_device)(TYPE(sdMat_ptr) vM, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t,M,vM,*iflag);

  MAGMA(setmatrix)(M->nrows,M->ncols,
      (const MAGMA_ST*)M->cpuData,M->stride,
      (MAGMA_ST*)M->gpuData,M->stride);
  *iflag = PHIST_SUCCESS;
}

void SUBR(sdMat_from_device)(TYPE(sdMat_ptr) vM, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t,M,vM,*iflag);

  MAGMA(getmatrix)(M->nrows,M->ncols,
      (const MAGMA_ST*)M->gpuData,M->stride,
      (MAGMA_ST*)M->cpuData,M->stride);
  *iflag = PHIST_SUCCESS;
}

