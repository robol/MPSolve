#!/bin/sh
#
# Autoconfiguratio script for MPSolve.
#
# Author: Leonardo Robol <leo@robol.it>

if [ ! -d "m4" ]; then
  mkdir m4
fi

autoreconf -i --force || exit 1
