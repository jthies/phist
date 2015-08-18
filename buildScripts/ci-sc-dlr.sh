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
FLAGS="" # optional-libs
ADD_CMAKE_FLAGS="" #optional CMake flags
# list of modules to load
MODULES_BASIC="cmake ccache cppcheck lapack"

declare -A MODULES_KERNELS
MODULES_KERNELS=( 
  ["builtin"]=""
  ["ghost"]="gsl"
  ["epetra"]="trilinos"
  ["tpetra"]="trilinos" )

declare -A MODULES_KERNELS_OPTIONAL
MODULES_KERNELS_OPTIONAL=(
  ["builtin"]="parmetis trilinos"
  ["ghost"]="trilinos"
  ["epetra"]=""
  ["tpetra"]="" )


## parse command line arguments
usage() { echo "Usage: $0 [-k <builtin|ghost|epetra|tpetra>] [-e <PrgEnv/module-string>] [-f <optional-libs>]"; \
          echo "       [-c <cmake flags to be added>" 1>&2; exit 1; }

while getopts "k:e:f:c:h" o; do
    case "${o}" in
        k)
            KERNELS=${OPTARG}
            ((KERNELS == "builtin" || KERNELS == "ghost" || KERNELS == "tpetra" || KERNELS == "epetra")) || usage
            ;;
        e)
            PRGENV=${OPTARG}
            ;;
        f)
            FLAGS=${OPTARG}
            ;;
        c)
            ADD_CMAKE_FLAGS=${OPTARG}
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

echo "Options: KERNEL_LIB=${KERNELS}, PRGENV=${PRGENV}, FLAGS=${FLAGS}, ADD_CMAKE_FLAGS='${ADD_CMAKE_FLAGS}'"

## prepare system for compilation
# configure modulesystem
export MODULEPATH=/tools/modulesystem/modulefiles
module() { eval `/usr/bin/modulecmd bash $*`; }

# load modules
module load "PrgEnv/$PRGENV"
for m in $MODULES_BASIC; do module load $m; done
for m in ${MODULES_KERNELS["$KERNELS"]}; do module load $m; done
if [[ "$FLAGS" = *"optional-libs"* ]]; then
  for m in ${MODULES_KERNELS_OPTIONAL["$KERNELS"]}; do module load $m; done
fi
module list

# be verbose from here on
set -x

# use ccache to speed up build
if [[ "$PRGENV" = "gcc"* ]]; then
  export FC="ccache gfortran" CC="ccache gcc" CXX="ccache g++"
elif [[ "$PRGENV" = "intel"* ]]; then
  export FC=ifort CC=icc CXX=icpc
fi

# Make Intel OpenMP and MKL deterministic (e.g. calculate identical results on different MPI procs)
export KMP_DETERMINISTIC_REDUCTION=1 MKL_CBWR="COMPATIBLE"

# "gcc -fsanitize=address" requires this
ulimit -v unlimited

## actually build and run tests
error=0

# release build including doc
if [[ "$KERNELS" = "ghost" ]]; then
  # this is the easiest way to make phist find ghost+dependencies
  export CMAKE_PREFIX_PATH=$PWD/../install-${PRGENV}-Release/lib/ghost:$CMAKE_PREFIX_PATH
  export PKG_CONFIG_PATH=$PWD/../install-${PRGENV}-Release/lib/pkgconfig:$PKG_CONFIG_PATH
  # also set the LD_LIBRARY_PATH appropriately
  export LD_LIBRARY_PATH=$PWD/../install-${PRGENV}-Release/lib/ghost:$PWD/../install-${PRGENV}-Release/lib/essex-physics:$LD_LIBRARY_PATH
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
      ..                                || error=1
make doc &> doxygen.log                 || error=1
make -j 24 || make                      || error=1
echo "Running tests. Output is compressed and written to test.log.gz"
make check 2>&1 | gzip -c > test.log.gz || error=1

make install &> install.log             || error=1
echo "Check installation with pkg-config project"
mkdir jdqr_pkg_config; cd $_
PKG_CONFIG_PATH=../$INSTALL_PREFIX/lib/pkgconfig:$PKG_CONFIG_PATH cmake ../../exampleProjects/jdqr_pkg_config || error=1
make || error=1
cd ..
echo "Check installation with CMake project"
mkdir jdqr_cmake; cd $_
CMAKE_PREFIX_PATH=../$INSTALL_PREFIX:$CMAKE_PREFIX_PATH cmake ../../exampleProjects/jdqr_cmake || error=1
make || error=1
cd ..

cd ..

# debug build
if [[ "$KERNELS" = "ghost" ]]; then
  # this is the easiest way to make phist find ghost+dependencies
  export CMAKE_PREFIX_PATH=$PWD/../install-${PRGENV}-Debug/lib/ghost:$CMAKE_PREFIX_PATH
  export PKG_CONFIG_PATH=$PWD/../install-${PRGENV}-Debug/lib/pkgconfig:$PKG_CONFIG_PATH
  # also set the LD_LIBRARY_PATH appropriately
  export LD_LIBRARY_PATH=$PWD/../install-${PRGENV}-Debug/lib/ghost:$PWD/../install-${PRGENV}-Debug/lib/essex-physics:$LD_LIBRARY_PATH
fi
mkdir build_${KERNELS}_${PRGENV}_Debug_${FLAGS// /_}; cd $_
cmake -DCMAKE_BUILD_TYPE=Debug    \
      -DPHIST_KERNEL_LIB=$KERNELS \
      -DINTEGRATION_BUILD=On      \
      -DGCC_SANITIZE=address      \
      ${ADD_CMAKE_FLAGS} \
      ..                                || error=1
make -j 24 || make                      || error=1
echo "Running tests. Output is compressed and written to test.log.gz"
make check 2>&1 | gzip -c > test.log.gz || error=1
make audit                              || error=1
cd ..


# return error code
exit $error
