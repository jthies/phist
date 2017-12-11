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
ADD_CMAKE_FLAGS="-DPHIST_BENCH_LARGE_N=-1" #optional CMake flags # -1 disables benchmarks to speed up build jobs!
WORKSPACE="$PWD/.."
VECT_EXT="native"
TRILINOS_VERSION="11.12.1"
# list of modules to load
MODULES_BASIC="cmake cppcheck gcovr doxygen"
# GCC_SANITIZE flag for debug mode, disabled for CUDA
SANITIZER="address"


## parse command line arguments
usage() { echo "Usage: $0 [-k <builtin|ghost|epetra|tpetra|petsc|Eigen>] [-e <PrgEnv/module-string>] [-f <optional-libs>]"; \
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
            ((KERNELS == "builtin" || KERNELS == "ghost" || KERNELS == "tpetra" || KERNELS == "epetra" || KERNELS == "petsc" || KERNELS == "Eigen")) || usage
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

echo "Options: KERNEL_LIB=${KERNELS}, PRGENV=${PRGENV}, FLAGS=${FLAGS}, ADD_CMAKE_FLAGS='${ADD_CMAKE_FLAGS}', VECT_EXT=${VECT_EXT}"

declare -A MODULES_KERNELS
MODULES_KERNELS=( 
  ["builtin"]=""
  ["ghost"]="gsl"
  ["epetra"]="trilinos/trilinos-${TRILINOS_VERSION}"
  ["tpetra"]="trilinos/trilinos-${TRILINOS_VERSION}" 
  ["petsc"]="petsc" 
  ["Eigen"]="Eigen" )

declare -A MODULES_KERNELS_OPTIONAL
MODULES_KERNELS_OPTIONAL=(
  ["builtin"]="ColPack parmetis trilinos"
  ["ghost"]="ColPack trilinos/trilinos-11.12.1"
  ["epetra"]=""
  ["tpetra"]=""
  ["petsc"]="" 
  ["Eigen"]="" )


## prepare system for compilation
# configure modulesystem
export MODULEPATH=/tools/modulesystem/modulefiles
module() { eval `/usr/bin/modulecmd bash $*`; }

# load modules
module load "PrgEnv/$PRGENV"
for m in $MODULES_BASIC; do module load $m; done
for m in ${MODULES_KERNELS["$KERNELS"]}; do module load $m; done
if [[ "$FLAGS" = *optional-libs* ]]; then
  for m in ${MODULES_KERNELS_OPTIONAL["$KERNELS"]}; do module load $m; done
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
  module load cuda
  SANITIZER=""
  nvidia-smi
  export CUDA_VISIBLE_DEVICES=0
fi
module list

# be verbose from here on
set -x

if [[ "$PRGENV" =~ gcc* ]]; then
  export FC=gfortran CC=gcc CXX=g++
  module load lapack
  if [ "${VECT_EXT}" != "CUDA" && "${PRGENV}" != "gcc-7.2.0-openmpi"]; then
    module load ccache
    ADD_CMAKE_FLAGS+="-DPHIST_USE_CCACHE=ON"
  else
    ADD_CMAKE_FLAGS+="-DPHIST_USE_CCACHE=OFF"
  fi
elif [[ "$PRGENV" =~ intel* ]]; then
  export FC=ifort CC=icc CXX=icpc
  
  module load mkl
  # make CMake find and use MKL:
  ADD_CMAKE_FLAGS+="-DBLA_VENDOR=Intel10_64lp"
fi

# Make Intel OpenMP and MKL deterministic (e.g. calculate identical results on different MPI procs)
export KMP_DETERMINISTIC_REDUCTION=1 MKL_CBWR="COMPATIBLE"

# "gcc -fsanitize=address" requires this
ulimit -v unlimited

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
make doc &> doxygen.log                 || update_error ${LINENO}
make -j 24 || make                      || update_error ${LINENO}
echo "Running tests. Output is compressed and written to test.log.gz"
make check &> test.log                  || update_error ${LINENO}
if [ "${VECT_EXT}" = "CUDA" ]; then
  echo "Check if it actually ran on our Tesla card"
  fgrep "1x Tesla" test.log             || update_error ${LINENO}
fi
gzip test.log                           || update_error ${LINENO}
echo "Install..."
make install &> install.log             || update_error ${LINENO}
echo "Check installation with pkg-config project"
mkdir jdqr_pkg_config; cd $_
PKG_CONFIG_PATH=../$INSTALL_PREFIX/lib/pkgconfig:$PKG_CONFIG_PATH cmake ../../exampleProjects/jdqr_pkg_config || update_error ${LINENO}
make || update_error ${LINENO}
cd ..
echo "Check installation with CMake project"
mkdir jdqr_cmake; cd $_
CMAKE_PREFIX_PATH=../$INSTALL_PREFIX:$CMAKE_PREFIX_PATH cmake ../../exampleProjects/jdqr_cmake || update_error ${LINENO}
make || update_error ${LINENO}
cd ..

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
  fgrep "1x Tesla" test.log             || update_error ${LINENO}
fi
gzip test.log                           || update_error ${LINENO}
make audit                              || update_error ${LINENO}
cd ..


# return error code
exit $error
