#!/bin/bash

INSTALL_PREFIX=$1
if [[ "$#" = "2" ]]; then
  GHOST_PREFIX=$2
fi
error=0

function update_error { 
if [[ "${error}" = "0" ]]; then
  error=$1
fi
}

echo "Check installation in $INSTALL_PREFIX"

echo "... with pkg-config project"

rm -rf jdqr_pkg_config; mkdir jdqr_pkg_config; cd $_
PKG_CONFIG_PATH_OLD=${PKG_CONFIG_PATH}
PKG_CONFIG_PATH=$INSTALL_PREFIX/lib/pkgconfig:$PKG_CONFIG_PATH
if [[ "$#" = "2" ]]; then
  PKG_CONFIG_PATH=$GHOST_PREFIX/lib/pkgconfig:$PKG_CONFIG_PATH
fi
echo "PKG_CONFIG_PATH='$PKG_CONFIG_PATH'"
CC=`pkg-config --variable=cc phist`
CXX=`pkg-config --variable=cxx phist`
FC=`pkg-config --variable=fc phist`

cmake -DCMAKE_CXX_COMPILER=$CXX \
      -DCMAKE_C_COMPILER=$CC \
      -DCMAKE_Fortran_COMPILER=$FC \
      ../../exampleProjects/jdqr_pkg_config || update_error $LINENO
make || update_error $LINENO
./Djdqr || update_error $LINENO
cd ..
PKG_CONFIG_PATH=${PKG_CONFIG_PATH_OLD}

echo "... with CMake project"
rm -rf jdqr_cmake; mkdir jdqr_cmake; cd $_
# cmake in this project should automatically use the compilers specified in phistConfig.cmake
unset CC
unset CXX
unset FC
CMAKE_PREFIX_PATH=$INSTALL_PREFIX:$CMAKE_PREFIX_PATH
cmake ../../exampleProjects/jdqr_cmake || update_error $LINENO
make || update_error $LINENO
./Djdqr || update_error $LINENO
cd ..

exit ${error}
