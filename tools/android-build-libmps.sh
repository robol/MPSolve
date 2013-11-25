#!/bin/bash
#
# This script can be used to build the Android toolchain in
# $srcdir/android-ext, along with a copy of libgmp and libmps.
#
# This installation can be used to compile xmpsolve for Android
# by using the QMake project file provided in src/xmpsolve.
#
# Please run this script relative to the top srcdir folder, i.e.,:
# run:
# $ ./tools/android-build-libmps.sh
#
# Author: Leonardo Robol <leonardo.robol@sns.it>
# Date: 24 November 2013

GMP_VERSION="5.1.3"
ANDROID_ARCH="arm-linux-androideabi"

function step {
  echo -e "> $1"
}

function die {
  echo $1
  exit 1
}

# Preliminary steps
SCRIPT=$(readlink -f $0)
SCRIPTPATH=$(dirname $SCRIPT)
SRCDIR=$(echo $SCRIPTPATH | sed "s-/tools\$--" )

step "Checking if you already have an android-ext dir"

if [ -d "$SRCDIR/android-ext" ]; then
  echo    ""
  echo    "> You appear to have android-ext. Do you wish to continue anyway overwriting"
  echo -n "> your current installation? [y/n]: "
  read ANS
  if [ "$ANS" != "y" ]; then
    echo "Aborting."
    exit 0
  else
    rm -r $SRCDIR/android-ext
  fi
fi

step "Copying the Android toolchain in place"

$ANDROID_NDK_ROOT/build/tools/make-standalone-toolchain.sh \
	--install-dir=$SRCDIR/android-ext \
	|| die "Cannot build a standalone toolchain. Please check your NDK installation"

# Set up environment for Android
export PATH="$SRCDIR/android-ext/bin:$PATH"
export PKG_CONFIG_LIBDIR="$SRCDIR/android-ext/lib/pkgconfig"

step "Downloading a copy of libgmp"

cd $SRCDIR/android-ext
mkdir tarballs && cd tarballs

( wget -q http://ftp.gmplib.org/gmp/gmp-$GMP_VERSION.tar.bz2 && \
tar xf gmp-$GMP_VERSION.tar.bz2 && cd gmp-$GMP_VERSION )|| \
	die "Cannot download GMP. Check your Internet connectivity"

step "Building GMP for Android"

( cd $SRCDIR/android-ext/tarballs/gmp-$GMP_VERSION && \
	./configure --host=$ANDROID_ARCH --prefix=$SRCDIR/android-ext && make -j4 && make install ) ||
 	die "Cannot build and install GMP, aborting."

step "Building a copy of libmps"

( mkdir mpsolve-build && cd mpsolve-build && \
  $SRCDIR/configure --host=$ANDROID_ARCH --prefix=$SRCDIR/android-ext --disable-examples \
	--disable-ui CFLAGS="-DMPS_USE_BUILTIN_COMPLEX" GMP_CFLAGS="-I$SRCDIR/android-ext/include" \
	GMP_LIBS="-L$SRCDIR/android-ext/lib -lgmp" && make -j4 && make install ) || \
	die "Cannot build MPSolve for Android, aborting."

