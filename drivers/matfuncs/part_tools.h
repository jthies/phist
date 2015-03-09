#include "matfuncs.h"

//! very simple reindexing routine for 2D problems on a cartesian
//! grid (like graphene or matpde). To initialize, call it with
//! arguments n1 and n2 ("grid size") (not thread-safe). To
//! use it, pass in the linear-index 'row' and returns the new
//! global index. 
//!
//! Currently only works for the very simple situation where the
//! domain can be split into equally sized square blocks, e.g.
//! 100x200 on 2 procs, 200x200 on 4 procs etc. Will attempt to
//! find such a partitioning on the init call, otherwise turn off
//! partitioning alltogether with a warning.
ghost_gidx_t iperm2d( ghost_gidx_t row, gidx_t arg2 );

//! same as the above for 3D problems like Anderson model problem.
ghost_gidx_t iperm3d( ghost_gidx_t row, gidx_t arg2, gidx_t arg3 );

