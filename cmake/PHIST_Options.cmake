option(PHIST_USE_CCACHE "Enable/disable using ccache to speed up repeated builds if you use the GNU compiler and ccache is found." ON)

# performance-related settings
option(PHIST_TRY_TO_PIN_THREADS "Try to pin the threads to cores, the exact implementation depends on the kernel library, hence the 'TRY'. A kernel library may or may not bind the threads without further warnings." On)
option(PHIST_TRY_TO_RESPECT_NUMA "Try to first-touch multi-vector data after allocation, makes the code run better on NUMA systems but temporary vectors become more expensive." On)

# output behavior
if(CMAKE_BUILD_TYPE MATCHES "Release")
        set(PHIST_OUTLEV 2 CACHE STRING "PHIST verbosity (release 2, 'INFO')")
else()
        set(PHIST_OUTLEV 4 CACHE STRING "PHIST verbosity (debug 4, 'EXTREME')")
endif()

if(CMAKE_BUILD_TYPE MATCHES "Debug")
  option(PHIST_TESTING "Enable potentially expensive checks in the code" ON)
else()
  option(PHIST_TESTING "Enable potentially expensive checks in the code" OFF)
endif()


if (PHIST_KERNEL_LIB STREQUAL "eigen")
  option(PHIST_ENABLE_MPI "Enable MPI within PHIST" OFF)
else()
  option(PHIST_ENABLE_MPI "Enable MPI within PHIST" ON)
endif()
option(PHIST_ENABLE_OPENMP "Enable OpenMP within PHIST" ON)
option(PHIST_USE_LIKWID "Enable instrumentation for likwid-perfctr" OFF)
option(PHIST_BUILD_MIC "Build for Intel MIC" OFF)

# reproducible random number generator?
# TODO - this is not quite working yet, once it is, make it the default
option(PHIST_BUILTIN_RNG
        "use random number generator provided by PHIST, which may be solower than the one implemented by the kernel lib. On the other hand, it gives reproducible results across runs, independent of the exact runtime setup."
        ON)

# use GHOST if found (must be ON if ghost is the kernel lib)
option(PHIST_USE_GHOST
        "try to find the GHOST library even if it is not the kernel lib" 
        OFF)

# try to find MPACK for high precision lapack routines
option(PHIST_USE_MPACK
        "try to find the MPACK library" 
        ON)

# try to find third-party solver packages (depending on kernel lib),
# for instance for epetra/tpetra and GHOST we have an interface to
# Anasazi and Belos.
option(PHIST_USE_SOLVER_TPLS
        "try to find supported third-party libraries that implement iterative methods for the kernel library used."
        ON)

# try to find third-party preconditioning packages (depending on kernel lib),
# for instance if the kernel lib is epetra, look for Ifpack and ML.
option(PHIST_USE_PRECON_TPLS
        "try to find supported third-party libraries that implement preconditioners for the kernel library used."
        ON)

# try to find third-party partitioning packages (depending on kernel lib),
# for instance if the kernel lib is epetra, look for Isorropia, builtin supports ParMETIS 
# and ColPack etc.
option(PHIST_USE_GRAPH_TPLS
        "try to find supported third-party libraries that implement graph partitioning and reordering for the kernel library used."
        ON)

# this option allows to disable all Trilinos libraries at once, overriding
# any of the choices PHIST_USE_*_TPLS
option(PHIST_USE_TRILINOS_TPLS
        "enable/disable using any optional Trilinos libraries like Zoltan, Belos, Kokkos etc."
        ON)

# option for using checkpoint/restart + automatic fault-tolerance features
option(PHIST_USE_CRAFT
                "enables/desables the use of checkpoint/restart and auto FT features using CRAFT"
                OFF)

#####################################################
# CODE INSTRUMENTATION
#####################################################

option(PHIST_TIMEMONITOR "Gather detailed function-level timings" On)
if( PHIST_TIMEMONITOR )
  message(STATUS "activated TimeMonitor timing monitoring")
endif()
option(PHIST_TIMEMONITOR_PERLINE "Gather timing information for each line with a PHIST_CHK_* macro!" Off)
option(PHIST_PERFCHECK "Check performance of operations (which have a performance model specification)" Off)
if( PHIST_PERFCHECK )
  message(STATUS "activated PerfCheck timing monitoring")
  if( PHIST_PERFCHECK AND PHIST_TIMEMONITOR )
    # these currently cannot coexist
    set(PHIST_TIMEMONITOR Off)
  endif()
  option(PHIST_PERFCHECK_SEPARATE_OUTPUT "write perfcheck output to separate files for each MPI rank (useful for analysing heterogenous runs)" OFF)
  option(PHIST_PERFCHECK_REALISTIC "Switch from ideal performance predictions to more realistic ones (e.g. for views)!" OFF)

  # these can be used to make the roofline model more accurate for compute-bound operations (not really applicable in 
  # phist right now, if we get such operations for e.g. wider mvecs we should maybe measure the flop rate instead)
  set(PHIST_PEAK_FLOPS_CPU -1 CACHE STRING "Peak GFlop/s rate of the CPU component (for the roofline model in PERFCHECK 
 only. If you don't know this value, just leave it at -1 and it will be ignored)")
  set(PHIST_PEAK_FLOPS_GPU -1 CACHE STRING "Peak GFlop/s rate of the GPU component (for the roofline model in PERFCHECK 
 only. If you don't know this value, just leave it at -1 and it will be ignored)")

endif()

