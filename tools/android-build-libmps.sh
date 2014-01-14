#!/bin/bash
#
# This script can be used to build the Android toolchain in
# $srcdir/android-ext-$ANDROID_ARCH, along with a copy of libgmp and libmps.
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
ANDROID_ARCH="arm-linux-androideabi-4.8"
# ANDROID_ARCH="x86-4.8"

if [ "$1" = "" ]; then
  # Make the user choose architecture he/she wants: 
  echo "> Select the desired architecture: "
  echo " 1. x86-4.8"
  echo " 2. arm-linux-androideabi-4.8"
  echo " 3. mipsel-linux-android-4.8"
  echo ""
  echo -n "> Choice: "
  read ANS

  case $ANS in
   1)
     ANDROID_ARCH="x86-4.8"
     ;;
   2)
     ANDROID_ARCH="arm-linux-androideabi-4.8"
     ;;
   3)
     ANDROID_ARCH="mipsel-linux-android-4.8"
     ;;
   *)
     echo "Invalid choice, aborting"
     exit 1
     ;;
  esac

else
  ANDROID_ARCH="$1"
fi

ANDROID_BUILD_DIR="android-ext-$ANDROID_ARCH"
ANDROID_BUILD_ARCH=$ANDROID_ARCH

if [ "$ANDROID_ARCH" == "x86-4.8" ]; then
  ANDROID_BUILD_ARCH="i686-linux-android"
fi

if [ "$ANDROID_ARCH" == "arm-linux-androideabi-4.8" ]; then
  ANDROID_BUILD_ARCH="arm-linux-androideabi"
fi

if [ "$ANDROID_ARCH" == "mipsel-linux-android-4.8" ]; then
  ANDROID_BUILD_ARCH="mipsel-linux-android"
fi

echo " *** Android build system for libmps *** "
echo " "
echo " > Building for architecture $ANDROID_ARCH in $ANDROID_BUILD_DIR"
echo ""

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

step "Checking if you already have an $ANDROID_BUILD_DIR dir"

if [ -d "$SRCDIR/$ANDROID_BUILD_DIR" ]; then
  echo    ""
  echo    "> You appear to have $ANDROID_BUILD_DIR. Do you wish to continue anyway overwriting"
  echo -n "> your current installation? [y/n]: "
  read ANS
  if [ "$ANS" != "y" ]; then
    echo "Aborting."
    exit 0
  else
    rm -r $SRCDIR/$ANDROID_BUILD_DIR
  fi
fi

step "Copying the Android toolchain in place"

$ANDROID_NDK_ROOT/build/tools/make-standalone-toolchain.sh --toolchain=$ANDROID_ARCH \
	--install-dir=$SRCDIR/$ANDROID_BUILD_DIR \
	|| die "Cannot build a standalone toolchain. Please check your NDK installation"

# Set up environment for Android
export PATH="$SRCDIR/$ANDROID_BUILD_DIR/bin:$PATH"
export PKG_CONFIG_LIBDIR="$SRCDIR/$ANDROID_BUILD_DIR/lib/pkgconfig"

step "Downloading a copy of libgmp"

cd $SRCDIR/$ANDROID_BUILD_DIR
mkdir tarballs && cd tarballs

( wget  ftp://ftp.gmplib.org/pub/gmp-$GMP_VERSION/gmp-$GMP_VERSION.tar.bz2 && \
tar xf gmp-$GMP_VERSION.tar.bz2 && cd gmp-$GMP_VERSION )|| \
	die "Cannot download GMP. Check your Internet connectivity"

step "Building GMP for Android"

( cd $SRCDIR/$ANDROID_BUILD_DIR/tarballs/gmp-$GMP_VERSION && \
	./configure --host=$ANDROID_BUILD_ARCH --prefix=$SRCDIR/$ANDROID_BUILD_DIR && make -j4 && make install ) ||
 	die "Cannot build and install GMP, aborting."

step "Building a copy of libmps"

( mkdir mpsolve-build && cd mpsolve-build && \
  $SRCDIR/configure --host=$ANDROID_BUILD_ARCH --prefix=$SRCDIR/$ANDROID_BUILD_DIR --disable-examples \
	--disable-ui CFLAGS="-DMPS_USE_BUILTIN_COMPLEX" GMP_CFLAGS="-I$SRCDIR/$ANDROID_BUILD_DIR/include" \
	GMP_LIBS="-L$SRCDIR/$ANDROID_BUILD_DIR/lib -lgmp" && make -j4 && make install ) || \
	die "Cannot build MPSolve for Android, aborting."

