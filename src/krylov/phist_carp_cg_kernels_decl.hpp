//! \name internal data structures for mixed-precision CARP-CG
//! \note thse should not be used directly by the user!
//!@{

//! e(x)tended sparse matrix type with complex shift and augmented 
//! by additional (dense) rows and columns.
typedef struct TYPE(x_sparseMat) {
  TYPE(const_sparseMat_ptr) A_;
  _MT_ const *sigma_r_;
  _MT_ const *sigma_i_;
  TYPE(const_mvec_ptr)      Vproj_;
} TYPE(const_x_sparseMat);

typedef TYPE(const_x_sparseMat) const* TYPE(const_x_sparseMat_ptr);


//! represent augmented complex vector of the form
//! |v + i vi |
//! |vp+ i vpi|
class TYPE(x_mvec) {
  TYPE(mvec_ptr)      v_;
  TYPE(mvec_ptr)      vi_;
  TYPE(sdMat_ptr)     vp_;
  TYPE(sdMat_ptr)     vpi_;
  
  bool own_mvecs_;

  public:

  //! constructor - does not allocate memory
  TYPE(x_mvec)();

  //! constructor that views existing mvecs and optionally takes ownership
  TYPE(x_mvec)(TYPE(mvec_ptr) v, TYPE(mvec_ptr) vi, int naug, bool take_ownership, int* iflag);

  //! imaginary part is allocated only if rc=true, augmented part only if naug>0.
  allocate(const_map_ptr_t map, int nvec, int naug, bool rc, int* iflag);

  //! destructor
  ~TYPE(x_mvec)();

  //!
  void deallocate();

};

//@}
