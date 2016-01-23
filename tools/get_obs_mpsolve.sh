#!/bin/bash
#

if [ -d "windows-mpsolve" ]; then
  rm -r windows-mpsolve
fi

mkdir windows-mpsolve
cd windows-mpsolve

cp ../download-mingw-rpm.py .

./download-mingw-rpm.py --no-clean --deps mingw32-libqt4 mingw32-libgmp10 mingw32-winpthreads mingw32-libpng
./download-mingw-rpm.py --no-deps --no-clean -p home:lrobol mpsolve

mv  usr/i686-w64-mingw32/sys-root/mingw/bin MPSolve
# rm -r usr/
# rm -r cache/

cd MPSolve

cp ../../installer.nsis .
makensis installer.nsis


