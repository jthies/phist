# try to find out if METIS is installed with 64-bit integers
include(CheckCXXSourceCompiles)

set(CMAKE_REQUIRED_INCLUDES "${ParMETIS_INCLUDE_DIRS}")

CHECK_CXX_SOURCE_COMPILES("
" METIS_USABLE
)
#include "metis.h"
#ifequal IDXTYPEWIDTH == 32
#error "METIS is not compiled with 64-bit integers!"
#endif
