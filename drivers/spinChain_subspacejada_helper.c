#include "ghost.h"
#include "ghost/util.h"
#include "matfuncs.h"

GHOST_REGISTER_DT_D(my_datatype_)

void init_mtraits(ghost_mtraits_t* mtraits)
{
    ghost_mtraits_t newmtraits = GHOST_MTRAITS_INIT(.format = GHOST_SPM_FORMAT_CRS, .flags = GHOST_SPM_DEFAULT, .datatype = my_datatype_);
    *mtraits = newmtraits;
}
