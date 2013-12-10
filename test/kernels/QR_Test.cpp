
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "gtest/gtest.h"
//#include "gmock/gmock.h"

#include "KernelTestWithVectors.h"
#include "KernelTestWithSdMats.h"

#include "phist_kernels.h"

#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

using namespace testing;

#ifdef CLASSNAME
#undef CLASSNAME
#endif

#define _N_ 24
#define _NV_ 1

#define CLASSNAME SQR_Test_24_1
#include "phist_gen_s.h"
#include "QR_Test_def.hpp"

#undef CLASSNAME
#define CLASSNAME DQR_Test_24_1

#include "phist_gen_d.h"
#include "QR_Test_def.hpp"

#undef CLASSNAME
#define CLASSNAME CQR_Test_24_1

#include "phist_gen_c.h"
#include "QR_Test_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZQR_Test_24_1

#include "phist_gen_z.h"
#include "QR_Test_def.hpp"

#undef _N_
#define _N_ 59
#undef _NV_
#define _NV_ 5

#undef CLASSNAME
#define CLASSNAME SQR_Test_59_5

#include "phist_gen_s.h"
#include "QR_Test_def.hpp"

#undef CLASSNAME
#define CLASSNAME DQR_Test_59_5

#include "phist_gen_d.h"
#include "QR_Test_def.hpp"

#undef CLASSNAME
#define CLASSNAME CQR_Test_59_5

#include "phist_gen_c.h"
#include "QR_Test_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZQR_Test_59_5

#include "phist_gen_z.h"
#include "QR_Test_def.hpp"

#if 1
// let's try something bigger...
#undef _N_
#define _N_ 9999
#undef _NV_
#define _NV_ 65

#undef CLASSNAME
#define CLASSNAME SQR_Test_10k_65

#include "phist_gen_s.h"
#include "QR_Test_def.hpp"

#undef CLASSNAME
#define CLASSNAME DQR_Test_10k_65

#include "phist_gen_d.h"
#include "QR_Test_def.hpp"

#undef CLASSNAME
#define CLASSNAME CQR_Test_10k_65

#include "phist_gen_c.h"
#include "QR_Test_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZQR_Test_10k_65

#include "phist_gen_z.h"
#include "QR_Test_def.hpp"

#endif
