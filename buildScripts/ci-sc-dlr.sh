#!/usr/bin/env bash
set -e

## default options and declarations
# kernel lib
KERNELS="builtin" # ghost epetra tpetra
FLAGS="" # " optional-libs"

# list of modules to load
MODULES_PRGENV="PrgEnv/gcc-4.9.2-openmpi"
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
  ["tpatra"]="" )


## parse command line arguments
usage() { echo "Usage: $0 [-k <builtin|ghost|epetra|tpetra>] [-e <modules-string>] [-f <optional-libs>]" 1>&2; exit 1; }

while getopts "k:e:f:h" o; do
    case "${o}" in
        k)
            KERNELS=${OPTARG}
            ((KERNELS == "builtin" || KERNELS == "ghost" || KERNELS == "tpetra" || KERNELS == "epetra")) || usage
            ;;
        e)
            MODULES_PRGENV=${OPTARG}
            ;;
        f)
            FLAGS=" ${OPTARG}"
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

echo "Options: KERNEL_LIB=${KERNELS}, MODULES_PRGENV=${MODULES_PRGENV}, FLAGS=${FLAGS}"
exit 0
## prepare system for compilation
# configure modulesystem
export MODULEPATH=/tools/modulesystem/modulefiles
module() { eval `/usr/bin/modulecmd bash $*`; }

# load modules
for m in $MODULES_PRGENV; do module load $m; done
for m in $MODULES_BASIC; do module load $m; done
for m in ${MODULES_KERNELS["$KERNELS"]}; do module load $m; done
if [[ "$FLAGS" = *"optional-libs"* ]]; then
  for m in ${MODULES_KERNELS_OPTIONAL["$KERNELS"]}; do module load $m; done
fi
module list

# detect compiler
if [[ "$MODULES_PRGENV" = *"gcc"* ]]; then
  COMPILER=gcc-$(gcc -dumpversion)
elif [[ "$MODULES_PRGENV" = *"intel"* ]]; then
  COMPILER=intel-$(icc -dumpversion)
else
  echo "Unknown compiler!"
  exit 1
fi

# use ccache to speed up build
if [[ "$COMPILER" = "gcc-"* ]]; then
  export FC="ccache gfortran" CC="ccache gcc" CXX="ccache g++"
elif [[ "$COMPILER" = "intel-"* ]]; then
  export FC="ccache ifort" CC="ccache icc" CXX="ccache icpc"
fi

# "gcc -fsanitize=address" requires this
ulimit -v unlimited


## actually build and run tests
error=0

# release build
mkdir build_${KERNELS}_${COMPILER}_release${FLAGS// /_}; cd $_
cmake -DCMAKE_BUILD_TYPE=Release  \
      -DPHIST_KERNEL_LIB=$KERNELS \
      ..                            || error=1
make doc                            || error=1
make -j 6 || make                   || error=1
make test                           || error=1
cd ..

# debug build
mkdir build_${KERNELS}_${COMPILER}_debug${FLAGS// /_}; cd $_
cmake -DCMAKE_BUILD_TYPE=Debug    \
      -DPHIST_KERNEL_LIB=$KERNELS \
      -DINTEGRATION_BUILD=On      \
      -DGCC_SANITIZE=address      \
      ..                            || error=1
make -j 6 || make                   || error=1
make test                           || error=1
make audit                          || error=1
cd ..


# return error code
exit $error
