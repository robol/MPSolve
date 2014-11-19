#!/bin/sh
#
# Autoconfiguration script for MPSolve.
#
# Author: Leonardo Robol <leo@robol.it>

echo "
 ****************************************
 *** MPSolve is initializing autoconf ***
 ***                                  ***
 *** Please use this command only if  ***
 *** you checked out the source from  ***
 *** the git repository!              ***
 ****************************************
"

if [ ! -d "m4" ]; then
  mkdir m4
fi

AUTORECONF=$(which autoreconf)
if [ $? -ne 0 ]; then
  echo "It appears that Autotools is not correctly installed on this system."
  echo "Please install it and re-run this script."
  exit 1
fi

$AUTORECONF -i --force || exit 3
