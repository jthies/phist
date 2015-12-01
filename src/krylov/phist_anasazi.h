#ifndef PHIST_ANASAZI_H
#define PHIST_ANASAZI_H

#include "phist_config.h"

#ifndef DOXYGEN

/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif
#include "phist_operator.h"

#endif //DOXYGEN

#include "phist_enums.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {

  BKS=0, // block Krylov Schur
  TMD=1  // TraceMin Davidson
} phist_anasaziType;

#ifdef PHIST_HAVE_SP
#include "phist_gen_s.h"
#include "phist_anasazi_decl.h"
#include "phist_gen_c.h"
#include "phist_anasazi_decl.h"
#endif
#include "phist_gen_d.h"
#include "phist_anasazi_decl.h"
#include "phist_gen_z.h"
#include "phist_anasazi_decl.h"
#include "phist_gen_clean.h"

#ifdef __cplusplus
}
#endif

#endif
