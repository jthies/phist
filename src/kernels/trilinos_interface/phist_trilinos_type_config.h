#ifndef _ST_
#error "this header should be included after a phist_gen_X.h header."
#endif

#ifndef PHIST_TRILINOS_TYPE_AVAIL
/*! \def PHIST_TRILINOS_TYPE_AVAIL

    after including e.g. phist_gen_c.h (which defines _ST_ etc.),
    this macro can be used to check wether the Trilinos installation
    used actually supports that data type. If not, PHIST_TRILIOS_TYPE_AVAIL
    will be #undefined. You still have to check if specific Trilinos packages
    are available, e.g. PHIST_HAVE_TEUCHOS, PHIST_HAVE_BELOS.

 */
#define PHIST_TRILINOS_TYPE_AVAIL 1
#endif

#ifdef PHIST_HAVE_TEUCHOS
# include "Teuchos_ConfigDefs.hpp"
# if defined(IS_COMPLEX) && !defined(HAVE_TEUCHOS_COMPLEX)
# warning "disabling Trilinos for complex types C, Z"
# undef PHIST_TRILINOS_TYPE_AVAIL
# endif
#endif

#ifdef PHIST_TRILINOS_TYPE_AVAIL

# ifdef PHIST_KERNEL_LIB_TPETRA
# include "TpetraCore_config.h"
#  if defined(IS_DOUBLE)&&defined(IS_COMPLEX)&&!defined(HAVE_TPETRA_INST_COMPLEX_DOUBLE)
#  warning "disabling Trilinos for double complex type Z"
#  undef PHIST_TRILINOS_TYPE_AVAIL
#  elif !defined(IS_DOUBLE)&&defined(IS_COMPLEX)&&!defined(HAVE_TPETRA_INST_COMPLEX_FLOAT)
#  warning "disabling Trilinos for float complex type C"
#  undef PHIST_TRILINOS_TYPE_AVAIL
#  elif defined(IS_DOUBLE)&&!defined(IS_COMPLEX)&&!defined(HAVE_TPETRA_INST_DOUBLE)
#  warning "disabling Trilinos for double real type D"
#  undef PHIST_TRILINOS_TYPE_AVAIL
#  elif !defined(IS_DOUBLE)&&!defined(IS_COMPLEX)&&!defined(HAVE_TPETRA_INST_FLOAT)
#  warning "disabling Trilinos for float real type S"
#  undef PHIST_TRILINOS_TYPE_AVAIL
#  endif
# endif
#endif
