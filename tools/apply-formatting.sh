#!/bin/sh
#

find src/{libmps,tests} include \( -name \*.c -or -name \*.h \) \
    -exec uncrustify --replace -c tools/mpsolve-style.cfg {} \;

# Cleanup backup files.
find . -name \*~ -delete
