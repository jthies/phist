#ifndef PHIST_PRECON_TRAITS_HPP
#define PHIST_PRECON_TRAITS_HPP

#include "phist_config.h"
#include "phist_enums.h"
#include "phist_ScalarTraits.hpp"

namespace phist {

//! traits class for making preconditioners usable in phist programs. The user
//! should not use this class directly but the C wrapper defined in phist_precon.h

//! default implementation of the traits class that allows us to interface
//! with our own and third-party preconditioners. The default implementation
//! just gives error messages about missing template specialization.
template<typename ST,phist_Eprecon PT>
class PreconTraits
{

  typedef ScalarTraits<ST> st;
  typedef void mvec_t;

public:


  static void Usage()
  {
    PHIST_SOUT(PHIST_ERROR,"In order to use a certain preconditioner, you first have to implement class PreconTraits for it.\n"
                           "This is currently only possible from within PHIST as the preconditioner has to exist as an enum phist_Eprecon\n"
                           "entry. However, you can also construct any preconditioning object your self and pass it to the solvers as\n"
                           "TYPE(linearOp), in which case you have complete freedom of implementation even if PHIST is pre-installed.\n");
  }

  static void NotImplemented(int* iflag)
  {
    PHIST_SOUT(PHIST_ERROR,"class PreconTraits is missing a specialization for data type %c and preconditioner type %s\n",
                        st::type_char(), precon2str(PT));
    *iflag=PHIST_NOT_IMPLEMENTED;
  }

  static void Create(void** P, 
        const void* A, ST sigma, const void* B, 
        const void* Vkern, const void* BVkern,
        const char* options, int* iflag)
  {
    NotImplemented(iflag);
    return;
  }
  static void Delete(void* P, int* iflag)
  {
    NotImplemented(iflag);
    return;
  }
  static void Apply(ST alpha, void const* P, mvec_t const* X, ST beta, mvec_t* Y, int* iflag)
  {
    NotImplemented(iflag);
    return;
  }
  static void ApplyT(ST alpha, void const* P, mvec_t const* X, ST beta, mvec_t* Y, int* iflag)
  {
    NotImplemented(iflag);
    return;
  }
  static void ApplyShifted(ST alpha, const void* P, ST const * sigma,
          mvec_t const* X, ST beta,  mvec_t* Y, int* iflag)
  {
    NotImplemented(iflag);
    return;
  }
};


template<typename ST>
class PreconTraits<ST,phist_NO_PRECON>
{

  typedef ScalarTraits<ST> st;
  typedef void mvec_t;

public:


  static void Usage()
  {
    PHIST_SOUT(PHIST_INFO,"NO_PRECON is equivalent to P=I, this 'preconditioner' does not\n"
                          "need any options\n");
  }

  static void Create(void** P, 
        const void* A, ST sigma, const void* B, 
        const void* Vkern, const void* BVkern,
        const char* options, int* iflag)
  {
    *P=NULL;
    return;
  }

  static void Delete(void* P, int* iflag)
  {
    *iflag=0;
    return;
  }
  
  static void Apply(ST alpha, void const* P, mvec_t const* X, ST beta, mvec_t* Y, int* iflag)
  {
    *iflag=0;
    PHIST_CHK_IERR(SUBR(mvec_add_mvec)(alpha,X,beta,Y,iflag),*iflag);
    return;
  } 

  static void ApplyT(ST alpha, void const* P, mvec_t const* X, ST beta, mvec_t* Y, int* iflag)
  {
    Apply(alpha,P,X,beta,Y,iflag);
  }
  static void ApplyShifted(ST alpha, const void* P, ST const * sigma,
          mvec_t const* X, ST beta,  mvec_t* Y, int* iflag)
  {
    int nvec;
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(X,&nvec,iflag),*iflag);
    _ST_ alphas[nvec];
    for (int i=0; i<nvec;i++) alphas[i]=alpha*(st::one()-sigma[i]);
    PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(alphas,X,beta,Y,iflag),*iflag);
    return;
  }
};

}

// depending on the installation, include source files with template specializations
# if defined(PHIST_KERNEL_LIB_EPETRA)
# include "tpl/phist_Ifpack_def.hpp"
# include "tpl/phist_ML_def.hpp"
# include "tpl/phist_MueLU_def.hpp"
# elif defined(PHIST_KERNEL_LIB_TPETRA)
# include "tpl/phist_Ifpack2_def.hpp"
# include "tpl/phist_Amesos2_def.hpp"
# include "tpl/phist_MueLU_def.hpp"
# endif

#endif
