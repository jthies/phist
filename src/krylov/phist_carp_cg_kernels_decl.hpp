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
