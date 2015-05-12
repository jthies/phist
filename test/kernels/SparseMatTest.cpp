#include "phist_config.h"
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "gtest/gtest.h"
//#include "gmock/gmock.h"

#include "kernels/phist_kernels.h"
#include "KernelTestWithVectors.h"
#include "KernelTestWithSdMats.h"
#include "../tools/MatrixIO.h"

#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

using namespace ::testing;

#define _N_ 25
#define _NV_ 1

#ifdef PHIST_HAVE_SP
# define CLASSNAME SSparseMatTest_25_1
# include "phist_gen_s.h"
# include "SparseMatTest_def.hpp"
# undef CLASSNAME

# define CLASSNAME CSparseMatTest_25_1
# include "phist_gen_c.h"
# include "SparseMatTest_def.hpp"
# undef CLASSNAME
#endif

#define CLASSNAME DSparseMatTest_25_1
#include "phist_gen_d.h"
#include "SparseMatTest_def.hpp"
#undef CLASSNAME

#define CLASSNAME ZSparseMatTest_25_1
#include "phist_gen_z.h"
#include "SparseMatTest_def.hpp"
#undef CLASSNAME

#undef _N_
#undef _NV_
#define _N_ 25
#define _NV_ 4

#ifdef PHIST_HAVE_SP

# define CLASSNAME SSparseMatTest_25_4
# include "phist_gen_s.h"
# include "SparseMatTest_def.hpp"
# undef CLASSNAME

# define CLASSNAME CSparseMatTest_25_4
# include "phist_gen_c.h"
# include "SparseMatTest_def.hpp"
#undef CLASSNAME

#endif

#define CLASSNAME DSparseMatTest_25_4
#include "phist_gen_d.h"
#include "SparseMatTest_def.hpp"
#undef CLASSNAME

#define CLASSNAME ZSparseMatTest_25_4
#include "phist_gen_z.h"
#include "SparseMatTest_def.hpp"

#undef CLASSNAME
#undef _N_
#undef _NV_
#define _N_ 25
#define _NV_ 7

#ifdef PHIST_HAVE_SP

# define CLASSNAME SSparseMatTest_25_7
# include "phist_gen_s.h"
# include "SparseMatTest_def.hpp"
# undef CLASSNAME

# define CLASSNAME CSparseMatTest_25_7
# include "phist_gen_c.h"
# include "SparseMatTest_def.hpp"
#undef CLASSNAME

#endif

#define CLASSNAME DSparseMatTest_25_7
#include "phist_gen_d.h"
#include "SparseMatTest_def.hpp"
#undef CLASSNAME

#define CLASSNAME ZSparseMatTest_25_7
#include "phist_gen_z.h"
#include "SparseMatTest_def.hpp"

//////////////////////////////////////////
// symmetric matrix tests               //
//////////////////////////////////////////

#undef CLASSNAME
#undef _N_
#undef _NV_
#define _N_ 20
#define _NV_ 8

#ifdef PHIST_HAVE_SP

# define CLASSNAME SSymmMatTest_20_8
# include "phist_gen_s.h"
# include "SymmMatTest_def.hpp"
# undef CLASSNAME

# define CLASSNAME CSymmMatTest_20_8
# include "phist_gen_c.h"
# include "SymmMatTest_def.hpp"
#undef CLASSNAME

#endif

# define CLASSNAME DSymmMatTest_20_8
#include "phist_gen_d.h"
#include "SymmMatTest_def.hpp"
#undef CLASSNAME

# define CLASSNAME ZSymmMatTest_20_8
#include "phist_gen_z.h"
#include "SymmMatTest_def.hpp"

#undef CLASSNAME
#undef _N_
#undef _NV_
#define _N_ 163
#define _NV_ 30

#ifdef PHIST_HAVE_SP

# define CLASSNAME SSymmMatTest_163_30
# include "phist_gen_s.h"
# include "SymmMatTest_def.hpp"
# undef CLASSNAME

# define CLASSNAME CSymmMatTest_163_30
# include "phist_gen_c.h"
# include "SymmMatTest_def.hpp"
#undef CLASSNAME

#endif

# define CLASSNAME DSymmMatTest_163_30
#include "phist_gen_d.h"
#include "SymmMatTest_def.hpp"
#undef CLASSNAME

# define CLASSNAME ZSymmMatTest_163_30
#include "phist_gen_z.h"
#include "SymmMatTest_def.hpp"
