#!/bin/sh
#
# Autoconfiguration script for MPSolve.
#
# Author: Leonardo Robol <leo@robol.it>

if [ ! -d "m4" ]; then
  mkdir m4
fi

AUTORECONF=$(which autoreconf)
if [ $? -ne 0 ]; then
  echo "It appears that Autotools is not correctly installed on this system."
  echo "Please install it and re-run this script."
  exit 1
fi

which libtool > /dev/null
if [ $? -ne 0 ]; then
  echo "It appears that Libtool is not correctly installed on this system."
  echo "Please install it and re-run this script."
  exit 2
fi

$AUTORECONF -i --force || exit 3
