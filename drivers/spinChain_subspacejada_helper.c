#include "ghost.h"
#include "ghost/util.h"
#include "matfuncs.h"

GHOST_REGISTER_DT_D(my_datatype_)

void init_mtraits(ghost_sparsemat_traits_t* mtraits)
{
    ghost_sparsemat_traits_t newmtraits = GHOST_SPARSEMAT_TRAITS_INITIALIZER;
    newmtraits.format = GHOST_SPARSEMAT_CRS;
    newmtraits.flags = GHOST_SPARSEMAT_DEFAULT;
    newmtraits.datatype = my_datatype_;
    *mtraits = newmtraits;
}
