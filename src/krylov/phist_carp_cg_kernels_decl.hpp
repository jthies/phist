//! \name internal data structures for mixed-precision CARP-CG
//! \note thse should not be used directly by the user!
//!@{

//! e(x)tended sparse matrix type with complex shift and augmented 
//! by additional (dense) rows and columns.
typedef struct TYPE(x_sparseMat) {
  TYPE(const_sparseMat_ptr) A_;
  _MT_ *sigma_r_;
  _MT_ *sigma_i_;
  TYPE(const_mvec_ptr)      Vproj_;
} TYPE(x_sparseMat);

//! represent augmented complex vector of the form
//! |v + i vi |
//! |vp+ i vpi|
class TYPE(x_mvec) 
{
  
  public:

    TYPE(mvec_ptr)      v_;
    TYPE(mvec_ptr)      vi_;
    TYPE(sdMat_ptr)     vp_;
    TYPE(sdMat_ptr)     vpi_;

  //! constructor - does not allocate memory
  TYPE(x_mvec)();

  //! constructor that views existing mvecs and optionally takes ownership
  TYPE(x_mvec)(TYPE(mvec_ptr) v, TYPE(mvec_ptr) vi, int naug, bool take_ownership, int* iflag);

  //! imaginary part is allocated only if rc=true, augmented part only if naug>0.
  void allocate(const_map_ptr_t map, int nvec, int naug, bool rc, int* iflag);

  //! destructor
  ~TYPE(x_mvec)();

  //!
  void deallocate();

  private:

  bool own_mvecs_;

};

//@}

//! \name some kernel functions for the special-purpose data types x_sparseMat and x_mvec
//@{

//!
void SUBR(x_mvec_add_mvec)(_ST_ alpha, TYPE(x_mvec) const* V,
        _ST_ beta, TYPE(x_mvec)* W, int* iflag);
                            
//! dot product of two possibly complex vectors, in the complex case the result is found in dots,
//! in the case of real arithmetic but complex vectors, dotsi is filled with the imaginary parts.
//! If the imaginary parts are not wanted dotsi may be NULL.
void SUBR(x_mvec_dot_mvec)(TYPE(x_mvec)* v, TYPE(x_mvec)* w,
                            _ST_   *dots, _MT_* dotsi, int *iflag);

//! Y = alpha X + beta Y
//!
//! yr = alpha_r xr + beta_r yr 
//!    - alpha_i xi 
//! yi = alpha_r xi + beta_r xi 
//!    + alpha_i xr 
void SUBR(x_mvec_vadd_mvec)(_ST_ const alphas[], _MT_ const alphas_i[], TYPE(x_mvec) const* X, _ST_ beta, TYPE(x_mvec)* Y, int* iflag);

//! scale columns i of v by real scalar alpha[i]
void SUBR(x_mvec_vscale)(TYPE(x_mvec)* v, _ST_ const alpha[], int* iflag);

//! there are two separate kernel functions, with and without augmented rows/cols, this wrapper calls the appropriate one for us.
void SUBR(x_carp_sweep)(TYPE(x_sparseMat) const* A,TYPE(const_mvec_ptr) b,TYPE(x_mvec)* x, void* carp_data, _MT_ const omega[],int *iflag);

//@}

#ifndef PHIST_CARP_CG_KERNEL_HELPERS
#define PHIST_CARP_CG_KERNEL_HELPERS

//! return true if both vectors have an allocated imaginary part,
//! false if none of them has, and throw an exception if only one of them has.
inline bool rc_variant(TYPE(x_mvec) const* v1, TYPE(x_mvec) const* v2)
{
#ifdef IS_COMPLEX
  return false;
#else
  if (v1->vi_==NULL && v2->vi_==NULL) return false;
  if (v1->vi_!=NULL && v2->vi_!=NULL) return true;
  throw "either both or none of the vectors must have an imaginary part!";
#endif
}

//! returns true if matrix and vectors are all 'complex in real arithmetic'
inline bool rc_variant(TYPE(x_sparseMat) const* A, TYPE(x_mvec) const* v1, TYPE(x_mvec) const* v2)
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
    if (rc_variant(v1,v2))
    {
      if (v1->vpi_!=NULL && v2->vpi_!=NULL) return true;
    }
    return true;
  }
  throw "either both or none of the vectors must be augmented!";
}

// returns true if both vectors and matrix are augmented by additional rows (rows and cols)
inline bool aug_variant(TYPE(x_sparseMat) const* A, TYPE(x_mvec) const* v1, TYPE(x_mvec) const* v2)
{
  bool rc=rc_variant(v1,v2);
  return rc && A->Vproj_!=NULL;
}

#endif // PHIST_CARP_CG_KERNEL_HELPERS
