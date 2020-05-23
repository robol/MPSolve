#!/bin/sh

#
# Simple script to automate building on Windows
#

set -e

function log {
    echo -n "[MSYS2-MAKE] "
    echo "$*"
}


if [ "$1" = "" ]; then
  echo "Please specify the path to the MPSolve tarball";
  exit 1
fi

echo -n "Update packages? [yn]: "
read ans

if [ "$ans" = "y" ] || [ "$ans" = "Y" ]; then
	log "Updating dependencies"
	pacman -Syu
	pacman -S make gmp-devel zip \
       		mingw-w64-x86_64-{gcc,gcc-fortran,make,nsis,pkg-config,gtk3}
	pacman -S --disable-download-timeout git autoconf libtool mingw-w64-x86_64-{gmp,qt5}
fi

mpsolve_tarball=$(readlink -f $1)
mpsolve_version=$(echo ${mpsolve_tarball} | sed "s/.*mpsolve-//" | sed "s/\.tar\.gz//")

log "Building MPSolve, version ${mpsolve_version}; tarball at ${mpsolve_tarball}"

# Move to a temporary directory
tempdir=$(mktemp -d)

oldpath=$(pwd)

log "Moving to ${tempdir}"
cd ${tempdir}

log "Decompressing MPSolve"
tar xf ${mpsolve_tarball}
cd mpsolve-*

log "Starting to build"
./configure --disable-examples 
make -j4 
make install DESTDIR="$tempdir/install"

cd ../install/mingw64/bin
log "Copying dependencies here"
windeployqt.exe xmpsolve.exe

# Copy required DLLS
log "Copying required DLLs"
for dll in $(ldd xmpsolve.exe | grep -v  -i SYSTEM | grep -v Windows | cut -d '>' -f2 | cut -d '(' -f1); do
  if [ ! -f $(basename ${dll}) ]; then
    cp -v ${dll} .
  fi
done

for dll in $(ldd mpsolve.exe | grep -v  -i SYSTEM | grep -v Windows | cut -d '>' -f2 | cut -d '(' -f1); do
  if [ ! -f $(basename ${dll}) ]; then
    cp -v ${dll} .
  fi
done

log "Generating SetupMPSolve.exe"
cd ../
cp -v ../../mpsolve*/COPYING .
cp ../../mpsolve*/tools/installer.nsis .
makensis.exe installer.nsis
mv SetupMPSolve.exe ../../

cd ../
mv mingw64 MPSolve-${mpsolve_version}

zip -r MPSolve-${mpsolve_version}.zip MPSolve-${mpsolve_version}

mv MPSolve-${mpsolve_version}.zip ../
cd ../

cp -v MPSolve-${mpsolve_version}.zip ${oldpath}/
cp -v SetupMPSolve.exe ${oldpath}/SetupMPSolve-${mpsolve_version}.exe

log "You may want to remove ${tempdir}, now that everything is finished"
