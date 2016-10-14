#ifndef PHIST_MEM_OWNER_HPP
#define PHIST_MEM_OWNER_HPP

//! This file contains simple classes to wrap mvec, sdMat and sparseMat
//! pointers in C++ objects which will take care of their deletion when
//! destroyed themselves. There is no reference counting mechanism, but
//! the wrapper objects could be put into std::shared_ptr's to achieve 
//! automatic memory management.
//!
//! A typical use of these classes is:
//! 
//! {
//! TYPE(mvec_ptr) tmp_V;
//! PHIST_CHK_IERR(SUBR(mvec_create)(&tmp_V,map,nvec,iflag),*iflag);
//! mvecOwner<ST> ownV(tmpV);
//! PHIST_CHK_IERR(... some other calls ...)
//! ...
//! wherever the program leaves the scope, tmp_V is deleted without
//! the user having to call SUBR(mvec_delete)
//! }

//! mvec owner object
template<typename T> class MvecOwner
{
};

//! sdMat owner object
template<typename T> class SdMatOwner
{
};

//! sparseMat owner object
template<typename T> class SparseMatOwner
{
};

#include "phist_gen_d.h"
#include "phist_memOwner_decl.hpp"

#ifdef PHIST_HAVE_CMPLX
#include "phist_gen_z.h"
#include "phist_memOwner_decl.hpp"
#endif

#ifdef PHIST_HAVE_SP
#include "phist_gen_s.h"
#include "phist_memOwner_decl.hpp"

# ifdef PHIST_HAVE_CMPLX
# include "phist_gen_c.h"
# include "phist_memOwner_decl.hpp"
# endif
# include "phist_gen_clean.h"

#endif
#endif
