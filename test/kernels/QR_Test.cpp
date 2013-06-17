
#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "phist_kernels.h"

#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

using namespace testing;

#ifdef CLASSNAME
#undef CLASSNAME
#endif

#define _N_ 8
#define _NV_ 1

#define CLASSNAME SQR_Test_8_1
#include "phist_gen_s.h"
#include "QR_Test_def.hpp"

#undef CLASSNAME
#define CLASSNAME DQR_Test_8_1

#include "phist_gen_d.h"
#include "QR_Test_def.hpp"

#undef CLASSNAME
#define CLASSNAME CQR_Test_8_1

#include "phist_gen_c.h"
#include "QR_Test_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZQR_Test_8_1

#include "phist_gen_z.h"
#include "QR_Test_def.hpp"

#undef _N_
#define _N_ 12
#undef _NV_
#define _NV_ 5

#undef CLASSNAME
#define CLASSNAME SQR_Test_12_5

#include "phist_gen_s.h"
#include "QR_Test_def.hpp"

#undef CLASSNAME
#define CLASSNAME DQR_Test_12_5

#include "phist_gen_d.h"
#include "QR_Test_def.hpp"

#undef CLASSNAME
#define CLASSNAME CQR_Test_12_5

#include "phist_gen_c.h"
#include "QR_Test_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZQR_Test_12_5

#include "phist_gen_z.h"
#include "QR_Test_def.hpp"

