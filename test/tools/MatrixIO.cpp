
#include "MatrixIO.h"
#include <iostream>

#include "phist_kernels.h"
#include "phist_ScalarTraits.hpp"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef PHIST_HAVE_SP

#include "phist_gen_s.h"
#include "MatrixIO_def.hpp"

#include "phist_gen_c.h"
#include "MatrixIO_def.hpp"

#endif

#include "phist_gen_d.h"
#include "MatrixIO_def.hpp"

#include "phist_gen_z.h"
#include "MatrixIO_def.hpp"

#ifdef __cplusplus
} //extern "C"
#endif
