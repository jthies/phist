#include "TestWithType.h"

#ifdef PHIST_HAVE_SP

# include "phist_gen_s.h"
bool TestWithType<_ST_>::typeImplemented_ = false;

# include "phist_gen_c.h"
bool TestWithType<_ST_>::typeImplemented_ = false;

#endif

# include "phist_gen_d.h"
bool TestWithType<_ST_>::typeImplemented_ = false;

# include "phist_gen_z.h"
bool TestWithType<_ST_>::typeImplemented_ = false;
