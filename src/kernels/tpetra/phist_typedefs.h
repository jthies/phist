#ifndef PHIST_TYPEDEFS_H
#define PHIST_TYPEDEFS_H

#ifdef __cplusplus
#include <cstddef>

//! type of node-local indices
typedef std::ptrdiff_t lidx_t;

//! type of global indices
typedef std::ptrdiff_t gidx_t;

#else
#include <stddef.h>

//! type of node-local indices
typedef ptrdiff_t lidx_t;

//! type of global indices
typedef ptrdiff_t gidx_t;

#endif

#endif
