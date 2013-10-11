#ifndef PHIST_JADA_H
#define PHIST_JADA_H

typedef enum eigSort_t {

LM=0,   // largest magnitude
SM=1,   // smallest magnitude
LR=2,   // largest real part
SR=3   // smallest real part

} eigSort_t;

#include "phist_operator.h"
#ifdef __cplusplus
extern "C" {
#endif
#include "phist_gen_s.h"
#include "phist_jada_decl.h"
#include "phist_gen_d.h"
#include "phist_jada_decl.h"
#include "phist_gen_c.h"
#include "phist_jada_decl.h"
#include "phist_gen_z.h"
#include "phist_jada_decl.h"
#ifdef __cplusplus
}
#endif

#endif
