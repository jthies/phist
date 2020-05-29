#!/usr/bin/env bash
set -e
set -o pipefail

# we need a recent bash version
if [ "${BASH_VERSINFO}" -lt 4 ]; then
  echo "This script requires at least bash-4.0 to run."
  exit 1
fi

## default options and declarations
# kernel lib
KERNELS="builtin" # ghost epetra tpetra
PRGENV="gcc-5.1.0-openmpi" # intel-13.0.1-mpich gcc-4.9.2-openmpi
FLAGS="default" # optional-libs
#optional CMake flags 
# BENCH_LARGE_N=-1 disables benchmarks to speed up build jobs!
# XSDK_eNABLE_Fortran makes sure the Fortran interfaces are generated and tested
ADD_CMAKE_FLAGS="-DPHIST_BENCH_LARGE_N=-1 -DXSDK_ENABLE_Fortran=ON" 
WORKSPACE="$PWD/.."
VECT_EXT="native"
TRILINOS_VERSION="git"
CUDA_VERSION=8.0.6
# list of modules to load
# GCC_SANITIZE flag for debug mode, disabled for CUDA
# note: with any GCC>=7 I get false positives and non-zero exit codes
# despite setting environment variables such as ASAN_OPTIONS="exitcode=0", 
# so I disable the sanitizers completely for now
#SANITIZER="address"
SANITIZER=""


## parse command line arguments
usage() { echo "Usage: $0 [-k <builtin|ghost|epetra|tpetra|petsc|eigen>] [-e PrgEnv/<module-string>] [-f <optional-libs>]"; \
          echo "       [-c <cmake flags to be added>] [-v <SSE|AVX|AVX2|CUDA>] [-t <trilinos version>] [-w <workspace-dir>]" 1>&2; exit 1; }

function update_error { 
if [[ "${error}" = "0" ]]; then
  error=$1
fi
}

while getopts "k:e:f:c:v:w:t:h" o; do
    case "${o}" in
        k)
            KERNELS=${OPTARG}
            ((KERNELS == "builtin" || KERNELS == "ghost" || KERNELS == "tpetra" || KERNELS == "epetra" || KERNELS == "petsc" || KERNELS == "eigen")) || usage
            ;;
        e)
            PRGENV=${OPTARG}
            ;;
        f)
            FLAGS=${OPTARG}
            ;;
        c)
            
            ADD_CMAKE_FLAGS+=" ${OPTARG}"
            ;;
        v)
            VECT_EXT=${OPTARG}
            ;;
        w)
            WORKSPACE=${OPTARG}
            ;;
        t)
            TRILINOS_VERSION=${OPTARG}
            ;;
        h)
            usage
            ;;
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

echo "USER: ${USER}"
echo "HOST: ${HOSTNAME}"
echo "SOURCE DIR: ${PWD}"

## prepare system for compilation
# configure modulesystem
module() { eval `/usr/bin/modulecmd bash $*`; }
source /tools/modulesystem/spack_KP/share/spack/setup-env.sh

# load modules
module load PrgEnv/${PRGENV}||exit -1

if [[ "$FLAGS" = *optional-libs* ]]; then

  ADD_CMAKE_FLAGS+=" -DPHIST_USE_GRAPH_TPLS:BOOL=ON"
  ADD_CMAKE_FLAGS+=" -DPHIST_USE_PRECON_TPLS:BOOL=ON"
  ADD_CMAKE_FLAGS+=" -DPHIST_USE_SOLVER_TPLS:BOOL=ON"
else
  ADD_CMAKE_FLAGS+=" -DPHIST_USE_GRAPH_TPLS:BOOL=OFF"
  ADD_CMAKE_FLAGS+=" -DPHIST_USE_PRECON_TPLS:BOOL=OFF"
# note: presently epetra and tpetra require Belos for the interface to TSQR.
# PHIST can do without mvec_QR, so without optional-libs this is tested.
  ADD_CMAKE_FLAGS+=" -DPHIST_USE_SOLVER_TPLS:BOOL=OFF"
fi
if [ "${VECT_EXT}" = "CUDA" ]; then
  SANITIZER=""
  nvidia-smi
  export CUDA_VISIBLE_DEVICES=0
fi
module list

# be verbose from here on
set -x

if [[ $PRGENV =~ gcc* ]]; then
  if [[ $KERNELS =~ tpetra ]] && [[ $VECT_EXT =~ "CUDA" ]]; then
      export CXX=mpicxx
      export CC=mpicc
      # note: the trilinos module should set OMPI_CXX=nvcc_wrapper for us, phist/cmake will check that
      export OMPI_MPICXX_CXXFLAGS="-expt-extended-lambda"
      # we need a little more verbose output in the release build because
      # otherwise the word "Tesla" will not appear in the output and a check after running the tests will fail.
      ADD_CMAKE_FLAGS+=" -DPHIST_OUTLEV=3"
  else
    export FC=gfortran CC=gcc CXX=g++
  fi
  if [[ "${VECT_EXT}" != "CUDA" ]] && [[ "${PRGENV}" != "gcc-7.2.0-openmpi" ]]; then
    ADD_CMAKE_FLAGS+=" -DPHIST_USE_CCACHE=OFF"
#    ADD_CMAKE_FLAGS+=" -DPHIST_USE_CCACHE=ON"
#    export CCACHE_DIR=/localdata1/f_jkessx/.ccache/
  else
    ADD_CMAKE_FLAGS+=" -DPHIST_USE_CCACHE=OFF"
  fi
  if [[ "$PRGENV" =~ gcc-7* ]]; then
    # there's a problem in jada-tests, test STestSchurDecomp* with gcc 7.2.0 and -fsanitize=address, see #218
    SANITIZER=""
  fi
elif [[ "$PRGENV" =~ intel* ]]; then
  export FC=ifort CC=icc CXX=icpc
  
  # make CMake find and use MKL:
fi

# "gcc -fsanitize=address" requires this
ulimit -v unlimited

echo "Options: KERNEL_LIB=${KERNELS}, PRGENV=${PRGENV}, FLAGS=${FLAGS}, ADD_CMAKE_FLAGS='${ADD_CMAKE_FLAGS}', VECT_EXT=${VECT_EXT}"

## actually build and run tests
error=0

# release build including doc
if [ "$KERNELS" = "ghost" ]; then
  if [ "$FLAGS" = "optional-libs" ]; then
    POSTFIX=_optional-libs
  fi
  # this is the easiest way to make phist find ghost+dependencies
  export CMAKE_PREFIX_PATH=${WORKSPACE}/install-${PRGENV}-Release-${VECT_EXT}${POSTFIX}/lib/ghost:$CMAKE_PREFIX_PATH
  export PKG_CONFIG_PATH=${WORKSPACE}/install-${PRGENV}-Release-${VECT_EXT}${POSTFIX}/lib/pkgconfig:$PKG_CONFIG_PATH
  # also set the LD_LIBRARY_PATH appropriately
  export LD_LIBRARY_PATH=${WORKSPACE}/install-${PRGENV}-Release-${VECT_EXT}${POSTFIX}/lib/ghost:${WORKSPACE}/install-${PRGENV}-Release-${VECT_EXT}/lib/essex-physics:$LD_LIBRARY_PATH
fi

# let ctest print all output if there was an error!
export CTEST_OUTPUT_ON_FAILURE=1
INSTALL_PREFIX=../install_${KERNELS}_${PRGENV}_Release_${FLAGS// /_}
mkdir build_${KERNELS}_${PRGENV}_Release_${FLAGS// /_}; cd $_
cmake -DCMAKE_BUILD_TYPE=Release  \
      -DPHIST_KERNEL_LIB=$KERNELS \
      -DINTEGRATION_BUILD=On      \
      -DCMAKE_INSTALL_PREFIX=$INSTALL_PREFIX \
      ${ADD_CMAKE_FLAGS} \
      ..                                || update_error ${LINENO}
echo "make doc"
make doc &> doxygen.log                 || update_error ${LINENO}

echo "Make libs"
make -j 24 libs || make -j 1 libs       || update_error ${LINENO}

echo "Make drivers and test executables"
make -j 24 || make -j 1 || update_error ${LINENO}

echo "Running tests. Output is compressed and written to test.log.gz"
make check &> test.log                  || update_error ${LINENO}
if [ "${VECT_EXT}" = "CUDA" ]; then
  echo "Check if it actually ran on our Tesla card"
  fgrep "Tesla" test.log             || update_error ${LINENO}
fi
gzip test.log                           || update_error ${LINENO}
echo "Install..."
make install &> install.log             || update_error ${LINENO}
echo "Check installation"
make test_install || update_error ${LINENO}

cd ..

# debug build
if [ "$KERNELS" = "ghost" ]; then
  # this is the easiest way to make phist find ghost+dependencies
  export CMAKE_PREFIX_PATH=${WORKSPACE}/install-${PRGENV}-Debug-${VECT_EXT}/lib/ghost:$CMAKE_PREFIX_PATH
  export PKG_CONFIG_PATH=${WORKSPACE}/install-${PRGENV}-Debug-${VECT_EXT}/lib/pkgconfig:$PKG_CONFIG_PATH
  # also set the LD_LIBRARY_PATH appropriately
  export LD_LIBRARY_PATH=${WORKSPACE}/install-${PRGENV}-Debug-${VECT_EXT}/lib/ghost:${WORKSPACE}/install-${PRGENV}-Debug-${VECT_EXT}/lib/essex-physics:$LD_LIBRARY_PATH
fi
mkdir build_${KERNELS}_${PRGENV}_Debug_${FLAGS// /_}; cd $_
cmake -DCMAKE_BUILD_TYPE=Debug    \
      -DPHIST_KERNEL_LIB=$KERNELS \
      -DINTEGRATION_BUILD=On      \
      -DGCC_SANITIZE=$SANITIZER   \
      ${ADD_CMAKE_FLAGS} \
      ..                                || update_error ${LINENO}
make -j 24 || make                      || update_error ${LINENO}
echo "Running tests. Output is compressed and written to test.log.gz"
make check &> test.log                  || update_error ${LINENO}
if [ "${VECT_EXT}" = "CUDA" ]; then
  echo "Check if it actually ran on our Tesla card"
  fgrep "Tesla" test.log             || update_error ${LINENO}
fi
gzip test.log                           || update_error ${LINENO}
make audit                              || update_error ${LINENO}
cd ..


# return error code
exit $error
