#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.67])
AC_INIT([MPSolve], [3.0], [robol@poisson.phc.unipi.it])
AC_CONFIG_SRCDIR([Makefile.am])
AC_CONFIG_HEADERS([config.h])

LT_INIT

AM_INIT_AUTOMAKE([dist-bzip2])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AM_PROG_CC_C_O

MPS_CFLAGS="-I../../include -DNICE_DEBUG -std=c99"
UNISOLVE_CFLAGS="-I../../include -std=c99"
SECSOLVE_CFLAGS="-I../../include -DNICE_DEBUG -std=c99"

AC_SUBST(MPS_CFLAGS)
AC_SUBST(UNISOLVE_CFLAGS)
AC_SUBST(SECSOLVE_CFLAGS)

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lgmp':
AC_CHECK_LIB([gmp], [mpc_set])
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [pow])
# FIXME: Replace `main' with a function in `-lmps':
AC_CHECK_LIB([pthread], [pthread_create])

# Checks for header files.
AC_CHECK_HEADERS([float.h limits.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([memset modf pow sqrt strchr])

AC_CONFIG_FILES([Makefile
                 examples/Makefile
                 src/Makefile
		 src/unisolve/Makefile
		 src/secsolve/Makefile
		 src/libmps/Makefile])
AC_OUTPUT
