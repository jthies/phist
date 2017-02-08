#ifndef __LANCZOS_H__
#define __LANCZOS_H__

#include <mpi.h>
extern "C"{
extern void imtql1_(int *, matdt_t *, matdt_t *, int *);
extern void imtql1f_(int *, matdt_t *, matdt_t *, int *);
}

#endif
