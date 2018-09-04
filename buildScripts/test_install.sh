#!/bin/bash

INSTALL_PREFIX=$1
CC=$2
CXX=$3
FC=$4

echo "Check installation in $INSTALL_PREFIX"

echo "... with pkg-config project"
mkdir jdqr_pkg_config; cd $_
PKG_CONFIG_PATH=$INSTALL_PREFIX/lib/pkgconfig:$PKG_CONFIG_PATH
cmake ../../exampleProjects/jdqr_pkg_config || exit $LINENO
make || exit $LINENO
./Djdqr || exit $LINENO
cd ..
echo "... with CMake project"
mkdir jdqr_cmake; cd $_
CMAKE_PREFIX_PATH=$INSTALL_PREFIX:$CMAKE_PREFIX_PATH
cmake ../../exampleProjects/jdqr_cmake || exit $LINENO
make || exit $LINENO
./Djdqr || exit $LINENO
cd ..

exit 0
