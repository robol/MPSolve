#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.67])
AC_INIT([MPSolve], [3.0], [robol@poisson.phc.unipi.it])
AC_CONFIG_SRCDIR([Makefile.am])
AC_CONFIG_HEADERS([config.h])

AC_CONFIG_MACRO_DIR([m4])

LT_INIT

AM_INIT_AUTOMAKE([dist-bzip2])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AM_PROG_CC_C_O

GENERAL_CFLAGS="-I../../include -ffast-math -std=c99"

# Determine CFLAGS based on the options that the user has specified.
# The supported build type are
#
# --without-debug Do not code debug prints (that is not the default)
# 

AC_ARG_ENABLE([debug],
	[AS_HELP_STRING([--disable-debug], [Disable debug prints in the code])])

AS_IF([test "x$enable_debug" != "xno"], [
	GENERAL_CFLAGS="$GENERAL_CFLAGS -DNICE_DEBUG"
])

MPS_CFLAGS="$GENERAL_CFLAGS"
UNISOLVE_CFLAGS="$GENERAL_CFLAGS"
SECSOLVE_CFLAGS="$GENERAL_CFLAGS"

AC_SUBST(MPS_CFLAGS)
AC_SUBST(UNISOLVE_CFLAGS)
AC_SUBST(SECSOLVE_CFLAGS)

# Checks for libraries.
AC_CHECK_LIB([gmp], [mpc_set])
AC_CHECK_HEADERS([gmp.h], [mpsolve_found_gmp_headers=yes; break;])

AC_CHECK_LIB([m], [pow])

AC_CHECK_LIB([pthread], [pthread_create])
AC_CHECK_HEADERS([pthread.h], [mpsolve_found_pthread_headers=yes; break;])

# Errors if libraries are not found, in this case GMP and Pthreads.
AS_IF([test "x$mpsolve_found_gmp_headers" != "xyes"], 
	[AC_MSG_ERROR([Unable to find gmp headers, pleas install libgmp-dev])])

AS_IF([test "x$mpsolve_found_pthread_headers" != "xyes"],
	[AC_MSG_ERROR([Unable to find pthreads headers, please install libpthreads-dev])])


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
