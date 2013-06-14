#ifndef KERNEL_TEST_WITH_TYPE_H
#define KERNEL_TEST_WITH_TYPE_H

#include "KernelTest.h"

/** 
 */
template<typename ST>
class KernelTestWithType
{
public:

//! set to true if the kernel lib supports the type (S/D/C/Z), false otherwise
bool typeImplemented_;


/** Set up method.
 * Fills internal data vector with values 1.0, 2.0 and 3.0.
 */
virtual void SetUp() {
typeImplemented_=false;
}

virtual void TearDown() {
}

};

#include "phist_gen_s.h"
#include "KernelTestWithType_def.h"

#include "phist_gen_c.h"
#include "KernelTestWithType_def.h"

#include "phist_gen_d.h"
#include "KernelTestWithType_def.h"

#include "phist_gen_z.h"
#include "KernelTestWithType_def.h"

#endif
