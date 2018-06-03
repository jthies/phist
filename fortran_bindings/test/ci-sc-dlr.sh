#!/bin/bash
# configure modulesystem to load newer bash!
export MODULEPATH=/tools/modulesystem/modulefiles
module() { eval `/usr/bin/modulecmd bash $*`; }
export FC="mpif90"

BUILD_TYPE=Release
FLAGS=default
KERNELS=builtin
PRGENV=gcc-5.3.0-openmpi
ADD_CMAKE_FLAGS=""
TRILINOS_VERSION=git

while getopts "k:e:f:c:t:b:h" o; do
    case "${o}" in
        k)
            KERNELS=${OPTARG}
            ((KERNELS == "builtin" || KERNELS == "ghost" || KERNELS == "tpetra" || KERNELS == "epetra" || KERNELS == "petsc")) || usage
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
        t)
            TRILINOS_VERSION=${OPTARG}
            ;;
        b)
            BUILD_TYPE=${OPTARG}
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

echo "Options: KERNEL_LIB=${KERNELS}, PRGENV=${PRGENV}, FLAGS=${FLAGS}, ADD_CMAKE_FLAGS='${ADD_CMAKE_FLAGS}', BUILD_TYPE=${BUILD_TYPE}"

declare -A MODULES_KERNELS_OPTIONAL
MODULES_KERNELS_OPTIONAL=(
  ["builtin"]="ColPack parmetis trilinos/trilinos-${TRILINOS_VERSION}"
  ["ghost"]="ColPack trilinos/trilinos-11.12.1"
  ["epetra"]=""
  ["tpetra"]=""
  ["petsc"]="" 
  ["Eigen"]="" )

IDSTRING=${KERNELS}_${PRGENV}_${BUILD_TYPE}_${FLAGS}

module add PrgEnv/${PRGENV} cmake lapack

if [[ "$KERNELS" = "epetra" || "$KERNELS" = "tpetra" ]]; then
        # epetra can get very slow when compiled with OpenMP
        export OMP_NUM_THREADS=1
        module add trilinos/trilinos-${TRILINOS_VERSION}
fi

if [[ "$FLAGS" = *optional-libs* ]]; then
  for m in ${MODULES_KERNELS_OPTIONAL["$KERNELS"]}; do module load $m; done
fi
# make sure the phist installation is found

export CMAKE_PREFIX_PATH=/localdata/f_buildn/ESSEX_workspace/phist/FLAGS/${FLAGS}/KERNEL_LIB/${KERNELS}/PRGENV/${PRGENV}/label/SC-030083L/install_${IDSTRING}:${CMAKE_PREFIX_PATH}
export LD_LIBRARY_PATH=/localdata/f_buildn/ESSEX_workspace/phist/FLAGS/${FLAGS}/KERNEL_LIB/${KERNELS}/PRGENV/${PRGENV}/label/SC-030083L/install_${IDSTRING}/lib:${LD_LIBRARY_PATH}

echo "CMAKE_PREFIX_PATH='$CMAKE_PREFIX_PATH'"

export FC="mpif90"
mkdir build_${IDSTRING}
cd build_${IDSTRING}
cmake ..  || exit 1
make cf   || exit 2
make -j     || exit 3
make test || exit 4
exit 0
