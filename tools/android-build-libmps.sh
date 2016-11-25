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

GMP_VERSION="6.1.1"
ADDITIONAL_OPTIONS=""

if [ "$1" = "" ]; then
  # Make the user choose architecture he/she wants:
  echo "> Select the desired ABI: "
  echo " 1. x86"
  echo " 2. x86_64"
  echo " 3. armeabi"
  echo " 4. armeabi-v7a"
  echo " 5. arm64-v8a"
  echo " 6. mips"
  echo " 7. mips64"
  echo " 8. All the above"
  echo ""
  echo -n "> Choice: "
  read ANS
else
  ANS=$1
fi

  case $ANS in
   1)
     ANDROID_ARCH="x86"
     ANDROID_BUILD_DIR="android-ext-x86"
     ;;
   2)
     ANDROID_ARCH="x86_64"
     ANDROID_BUILD_DIR="android-ext-x86_64"
     ;;
   3)
     ANDROID_ARCH="arm"
     ANDROID_BUILD_DIR="android-ext-armeabi"
     ;;
   4)
     ANDROID_ARCH="arm"
     ADDITIONAL_OPTIONS="--enable-vfp3"
     ANDROID_BUILD_DIR="android-ext-armeabi-v7a"
     ;;

   5)
     ANDROID_ARCH="arm64"
     ANDROID_BUILD_DIR="android-ext-arm64-v8a"
     ;;
   6)
     ANDROID_ARCH="mips"
     ANDROID_BUILD_DIR="android-ext-mips"
     ;;
   7)
     ANDROID_ARCH="mips64"
     ADDITIONAL_GMP_OPTIONS="--disable-assembly"
     ANDROID_BUILD_DIR="android-ext-mips64"
     ;;
   8)
     # Call the script with all the possible values
     $0 1 || exit 2
     $0 2 || exit 3
     $0 3 || exit 4
     $0 4 || exit 5
     $0 5 || exit 6
     $0 6 || exit 7
     $0 7 || exit 7
     ;;
   *)
     echo "Invalid choice, aborting"
     exit 1
     ;;
  esac

ANDROID_BUILD_ARCH=$ANDROID_ARCH

if [ "$ANDROID_ARCH" == "x86" ]; then
  ANDROID_BUILD_ARCH="i686-linux-android"
fi

if [ "$ANDROID_ARCH" == "x86_64" ]; then
  ANDROID_BUILD_ARCH="x86_64-linux-android"
fi

if [ "$ANDROID_ARCH" == "arm64" ]; then
  ANDROID_BUILD_ARCH="aarch64-linux-android"
fi

if [ "$ANDROID_ARCH" == "arm" ]; then
  ANDROID_BUILD_ARCH="arm-linux-androideabi"
fi

if [ "$ANDROID_ARCH" == "mips" ]; then
  ANDROID_BUILD_ARCH="mipsel-linux-android"
fi

if [ "$ANDROID_ARCH" == "mips64" ]; then
  ANDROID_BUILD_ARCH="mips64el-linux-android"
fi

echo "> Android build system for libmps starting"
echo " "
echo "> Building for architecture $ANDROID_ARCH in $ANDROID_BUILD_DIR"
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

if [ -d "$SRCDIR/$ANDROID_BUILD_DIR/tarballs" ]; then
    # Always remove MPSolve
    rm -rf $SRCDIR/$ANDROID_BUILD_DIR/tarballs/

    step "Assuming your toolchain is ok. Remove $ANDROID_BUILD_DIR to restart from scratch"
else

step "Copying the Android toolchain in place"

    $ANDROID_NDK_ROOT/build/tools/make_standalone_toolchain.py --arch "$ANDROID_ARCH" \
	--install-dir "$SRCDIR/$ANDROID_BUILD_DIR" --force \
	|| die "Cannot build a standalone toolchain. Please check your NDK installation"

fi

# Set up environment for Android
export PATH="$SRCDIR/$ANDROID_BUILD_DIR/bin:$PATH"
export PKG_CONFIG_LIBDIR="$SRCDIR/$ANDROID_BUILD_DIR/lib/pkgconfig"

cd $SRCDIR/$ANDROID_BUILD_DIR
mkdir tarballs && cd tarballs

if [ ! -r "$SRCDIR/$ANDROID_BUILD_DIR/lib/libgmp.a" ]; then

  step "Downloading a copy of libgmp"

  ( wget  ftp://ftp.gmplib.org/pub/gmp-$GMP_VERSION/gmp-$GMP_VERSION.tar.bz2 && \
  tar xf gmp-$GMP_VERSION.tar.bz2 && cd gmp-$GMP_VERSION )|| \
      	die "Cannot download GMP. Check your Internet connectivity"

  step "Building GMP for Android"

  ( cd $SRCDIR/$ANDROID_BUILD_DIR/tarballs/gmp-$GMP_VERSION && \
        	./configure --enable-cxx --host=$ANDROID_BUILD_ARCH \
                            --prefix=$SRCDIR/$ANDROID_BUILD_DIR \
                              $ADDITIONAL_GMP_OPTIONS \
                              && make -j4 && make install ) ||
 	    die "Cannot build and install GMP, aborting."

else

  step "Skipping GMP build, already present. Remove $SRCDIR/$ANDROID_BUILD_DIR/lib/libgmp.a to perform this step"

fi

step "Building a copy of libmps"
step "Additional options for this build: $ADDITIONAL_OPTIONS"

( mkdir mpsolve-build && cd mpsolve-build && \
  $SRCDIR/configure --host=$ANDROID_BUILD_ARCH --prefix=$SRCDIR/$ANDROID_BUILD_DIR --disable-examples \
	--disable-ui --disable-documentation --disable-examples \
        CFLAGS="-DMPS_USE_BUILTIN_COMPLEX" GMP_CFLAGS="-I$SRCDIR/$ANDROID_BUILD_DIR/include" \
	$ADDITIONAL_OPTIONS \
	GMP_LIBS="-L$SRCDIR/$ANDROID_BUILD_DIR/lib -lgmp" && make -j4 && make install ) || \
	die "Cannot build MPSolve for Android, aborting."

