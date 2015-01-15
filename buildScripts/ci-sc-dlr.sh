#!/usr/bin/env bash
set -e

## default options and declarations
# kernel lib
KERNELS="builtin" # ghost epetra tpetra
PRGENV="gcc-4.9.2-openmpi" # intel-13.0.1-mpich gcc-4.8.2-openmpi
FLAGS="" # optional-libs

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
usage() { echo "Usage: $0 [-k <builtin|ghost|epetra|tpetra>] [-e <PrgEnv/module-string>] [-f <optional-libs>]" 1>&2; exit 1; }

while getopts "k:e:f:h" o; do
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
        h)
            usage
            ;;
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

echo "Options: KERNEL_LIB=${KERNELS}, PRGENV=${PRGENV}, FLAGS=${FLAGS}"

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

# use ccache to speed up build
if [[ "$PRGENV" = "gcc"* ]]; then
  export FC="ccache gfortran" CC="ccache gcc" CXX="ccache g++"
elif [[ "$PRGENV" = "intel"* ]]; then
  export FC=ifort CC=icc CXX=icpc
fi

# "gcc -fsanitize=address" requires this
ulimit -v unlimited

# do not execute meaningless complex tests if the kernel lib doesn't support complex arithmetic
if [[ "$KERNELS" = "tpetra" ]]; then
  export CMPLX_TESTS=ON
else
  export CMPLX_TESTS=OFF
fi

## actually build and run tests
error=0

# release build
mkdir build_${KERNELS}_${PRGENV}_release_${FLAGS// /_}; cd $_
cmake -DCMAKE_BUILD_TYPE=Release  \
      -DPHIST_KERNEL_LIB=$KERNELS \
      -DPHIST_ENABLE_COMPLEX_TESTS=${CMPLX_TESTS} \
      ..                            || error=1
make doc                            || error=1
make -j 6 || make                   || error=1
make test                           || error=1
cd ..

# debug build
mkdir build_${KERNELS}_${PRGENV}_debug_${FLAGS// /_}; cd $_
cmake -DCMAKE_BUILD_TYPE=Debug    \
      -DPHIST_KERNEL_LIB=$KERNELS \
      -DPHIST_ENABLE_COMPLEX_TESTS=${CMPLX_TESTS} \
      -DINTEGRATION_BUILD=On      \
      -DGCC_SANITIZE=address      \
      ..                            || error=1
make -j 6 || make                   || error=1
make test                           || error=1
make audit                          || error=1
cd ..


# return error code
exit $error
