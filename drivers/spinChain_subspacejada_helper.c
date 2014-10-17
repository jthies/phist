#include "phist_config.h"
#ifdef PHIST_KERNEL_LIB_GHOST

#include "ghost.h"

#ifdef PHIST_HAVE_ESSEX_PHYSICS
#include "essex-physics/matfuncs.h"
#else
#include "matfuncs.h"
#endif

GHOST_REGISTER_DT_D(my_datatype_)

void init_mtraits(ghost_sparsemat_traits_t* mtraits)
{
    ghost_sparsemat_traits_t newmtraits = GHOST_SPARSEMAT_TRAITS_INITIALIZER;
    newmtraits.format = GHOST_SPARSEMAT_CRS;
    newmtraits.flags = GHOST_SPARSEMAT_DEFAULT;
    newmtraits.datatype = my_datatype_;
    *mtraits = newmtraits;
}

#endif
