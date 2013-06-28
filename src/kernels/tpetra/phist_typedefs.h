#ifndef PHIST_TYPEDEFS_H
#define PHIST_TYPEDEFS_H

#ifdef __cplusplus
#ifndef NO_INCLUDES_IN_HEADERS
#include <cstddef>
#endif

// we want ptrdiff_t (aka long long int on 64 bit systems) as local index,
// but a bug in Trilinos prevents us from using it right now. until then,
// we use int as local index type

//! type of node-local indices
//typedef std::ptrdiff_t lidx_t;
typedef int lidx_t;

//! type of global indices
typedef std::ptrdiff_t gidx_t;

#else
#ifndef NO_INCLUDES_IN_HEADERS
#include <stddef.h>
#endif

//! type of node-local indices
//typedef ptrdiff_t lidx_t;
typedef int lidx_t;

//! type of global indices
typedef ptrdiff_t gidx_t;

#endif

#include "phist_void_aliases.h"

#endif
