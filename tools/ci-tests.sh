#!/bin/sh
#
# Run Travis tests. This has to be run in the main src folder.

# Test CMake build
mkdir cmake-build
cd cmake-build/
cmake .. || exit 1
make VERBOSE=1 || exit 1
make test || exit 1
cd ..

# Test autotools build
./autogen.sh
mkdir build
cd build/
../configure MEX=mkoctfile MEXTOPS='--mex' || exit 1
make V=1 || exit 1
make check V=1 || exit 1

