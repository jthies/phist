#include "matfuncs.h"

//! very simple reindexing routine for 2D problems on a cartesian
//! grid (like graphene or matpde). To initialize, call it with
//! arguments n1 and n2 ("grid size") (not thread-safe). To
//! use it, pass in 
//! arg2==-1: new linear global index, returns original GID (inverse perm)
//! arg2==+1: original linear GID, returns new linear GID.
//!
//! example: 
//!           4 5 6 7
//!          --------                    2 3 | 6 7
//! original: 0 1 2 3   perm2d(orig,+1): 0 1 | 4 5
//!
//! so perm2d(4,-1) on MPI rank 1 gives 2
//!
//! Currently only works for the very simple situation where the
//! domain can be split into equally sized square blocks, e.g.
//! 100x200 on 2 procs, 200x200 on 4 procs etc. Will attempt to
//! find such a partitioning on the init call, otherwise turn off
//! partitioning alltogether with a warning.
ghost_gidx perm2d( ghost_gidx row, gidx_t arg2 );

//! same as the above for 3D problems like Anderson model problem.
ghost_gidx perm3d( ghost_gidx row, gidx_t arg2, gidx_t arg3 );

