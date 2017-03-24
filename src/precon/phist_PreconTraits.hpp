/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
#ifndef PHIST_PRECON_TRAITS_HPP
#define PHIST_PRECON_TRAITS_HPP

#include "phist_config.h"
#include "phist_enums.h"
#include "phist_operator.h"
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
        const char* options, void* last_arg, int* iflag)
  {
    NotImplemented(iflag);
    return;
  }

  static void Update(void* P, const void* A, ST sigma, const void* B,
        const void* Vkern, const void* BVkern,
        int* iflag)
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


// specialization for phist_NO_PRECON: P=identity operator
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
        const char* options, void* last_arg, int* iflag)
  {
    *P=NULL;
    return;
  }
  static void Update(void* P, const void* A, ST sigma, const void* B,
        const void* Vkern, const void* BVkern,
        int* iflag)
  {
    return;
  }
  static void Delete(void* P, int* iflag)
  {
    *iflag=0;
    return;
  }
  
  static void Apply(ST alpha, void const* P, mvec_t const* X, ST beta, mvec_t* Y, int* iflag)
  {
    *iflag=PHIST_NOT_IMPLEMENTED;
    return;
  } 

  static void ApplyT(ST alpha, void const* P, mvec_t const* X, ST beta, mvec_t* Y, int* iflag)
  {
    Apply(alpha,P,X,beta,Y,iflag);
  }
  static void ApplyShifted(ST alpha, const void* P, ST const * sigma,
          mvec_t const* X, ST beta,  mvec_t* Y, int* iflag)
  {
    *iflag=PHIST_NOT_IMPLEMENTED;
    return;
  }
};

// specialization for phist_USER_PRECON: assume that the user has wrapped
// his preconditioner as a phist_XlinearOp and dispatch to the function  
// pointers here.
template<typename ST>
class PreconTraits<ST,phist_USER_PRECON>
{

  typedef ScalarTraits<ST> st;
  typedef typename st::linearOp_t linearOp;
  typedef typename st::mvec_t* mvec_ptr;
  typedef typename st::mvec_t const* const_mvec_ptr;

public:


  static void Usage()
  {
    PHIST_SOUT(PHIST_INFO,"USER_PRECON assumes that the preconditioner is provided by\n"
                          "the user in the form of a phist_XlinearOp.\n"
                          "The constructed preconditioner with all the linearOp functions set should be\n"
                          "passed to precon_create via the 'last_arg' argument.\n");
  }

  static void Create(void** P, 
        const void* A, ST sigma, const void* B, 
        const_mvec_ptr Vkern, const_mvec_ptr BVkern,
        const char* options, void* last_arg, int* iflag)
  {
    *P=last_arg;
    *iflag=0;
    return;
  }
  static void Update(void* P, const void* A, ST sigma, const void* B,
        const void* Vkern, const void* BVkern,
        int* iflag)
  {
    PHIST_CAST_PTR_FROM_VOID(linearOp,userOp,P,*iflag);
    if (userOp->update!=NULL)
    {
      // NOTE: can't update A or B via this interface unless the preconditioner
      // has his own pointers to/copies of the matrices
      PHIST_CHK_IERR(userOp->update(userOp->A,userOp->aux,sigma,Vkern,BVkern,iflag),*iflag);
    }
    else
    {
      static int first_time=1;
      if (first_time)
      {
        PHIST_SOUT(PHIST_WARNING,"your user-provided preconditioner does not have the 'update' function pointer\n"
                             "set, so either set it or disable preconUpdate in your options file.\n"
                             "Not updating the preconditioner.\n");
        first_time=0;
      }
    }
    return;
  }
  static void Delete(void* P, int* iflag)
  {
    PHIST_CAST_PTR_FROM_VOID(linearOp,userOp,P,*iflag);
    if (userOp->destroy!=NULL)
    {
      PHIST_CHK_IERR(userOp->destroy(userOp,iflag),*iflag);
    }
    *iflag=0;
    return;
  }
  
  static void Apply(ST alpha, void const* P, const_mvec_ptr X, ST beta, mvec_ptr Y, int* iflag)
  {
    PHIST_CAST_PTR_FROM_VOID(linearOp,userOp,P,*iflag);
    PHIST_CHK_IERR(*iflag = userOp->apply!=NULL?0:PHIST_INVALID_INPUT,*iflag);
    PHIST_CHK_IERR(userOp->apply(alpha,userOp->A,X,beta,Y,iflag),*iflag);
    return;
  } 

  static void ApplyT(ST alpha, void const* P, const_mvec_ptr X, ST beta, mvec_ptr Y, int* iflag)
  {
    PHIST_CAST_PTR_FROM_VOID(linearOp,userOp,P,*iflag);
    PHIST_CHK_IERR(*iflag = userOp->applyT!=NULL?0:PHIST_INVALID_INPUT,*iflag);
    PHIST_CHK_IERR(userOp->applyT(alpha,userOp->A,X,beta,Y,iflag),*iflag);
    return;
  }
  static void ApplyShifted(ST alpha, const void* P, ST const * sigma,
          const_mvec_ptr X, ST beta,  mvec_ptr Y, int* iflag)
  {
    PHIST_CAST_PTR_FROM_VOID(linearOp,userOp,P,*iflag);
    if (userOp->apply_shifted!=NULL)
    {
      PHIST_CHK_IERR(userOp->apply_shifted(alpha,userOp->A,sigma,X,beta,Y,iflag),*iflag);
    }
    else
    {
      Apply(alpha,P,X,beta,Y,iflag);
    }
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
