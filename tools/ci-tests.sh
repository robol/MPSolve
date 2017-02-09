#!/bin/sh
#
# Run Travis tests. This has to be run in the main src folder.

# Test CMake build
mkdir cmake-build
cd cmake-build/
cmake ..
make VERBOSE=1
make test
cd ..

# Test autotools build
./autogen.sh
mkdir build
cd build/
../configure MEX=mkoctfile MEXTOPS='--mex'
make V=1
make check V=1

