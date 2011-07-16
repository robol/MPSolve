#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.67])
AC_INIT([MPSolve], [3.0], [robol@poisson.phc.unipi.it])
AC_CONFIG_SRCDIR([Makefile.am])
AC_CONFIG_HEADERS([config.h])

# Just for now, won't be in final release, probabily,
# but it's just too useful for debugging. 
# AC_DISABLE_SHARED
AC_CONFIG_MACRO_DIR([m4])

LT_INIT

AM_INIT_AUTOMAKE([dist-bzip2])

# m4_ifdef([AM_SILENT_RULES],[AM_SILENT_RULES([yes])])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AM_PROG_CC_C_O
AC_PROG_CC_C99

GENERAL_CFLAGS=" -DNOMPTEMP"

# Determine CFLAGS based on the options that the user has specified.
# The supported build type are
#
# --without-debug Do not code debug prints (that is not the default)
# 

AC_ARG_ENABLE([debug],
	AS_HELP_STRING([--disable-debug], [Disable debug prints in the code]),
	[enable_debug=$enableval have_debug=$enableval],
	[enable_debug=yes have_debug=yes])
AC_MSG_RESULT([$enable_debug])

AS_IF([test x$enable_debug != xno], [
	GENERAL_CFLAGS="$GENERAL_CFLAGS -DNICE_DEBUG"
])

MPS_CFLAGS="$GENERAL_CFLAGS"
UNISOLVE_CFLAGS="$GENERAL_CFLAGS"
SECSOLVE_CFLAGS="$GENERAL_CFLAGS"
UNISOLVE_LIBS="-lm -lgmp"
SECSOLVE_LIBS="-lm -lgmp"

# Check the presence of the check modules used to perform
# the checking
PKG_CHECK_MODULES([CHECKING], [check >= 0.9.4], 
  [enable_check=yes], [enable_check=no])
AS_IF([test "x$enable_check" = "xyes"], [
	CHECK_PROGRAMS="check_secsolve check_unisolve"
])
AC_SUBST(CHECK_PROGRAMS)
AC_SUBST(CHECK_CFLAGS)
AC_SUBST(CHECK_LIBS)

# Check the existance of mkoctfile to compile octave 
# module
AC_CHECK_PROGS([MKOCTFILE], [mkoctfile])
enable_octave=no
AS_IF([test "$MKOCTFILE" != ""], [enable_octave=yes OCTAVE_MODULES_PROGRAMS=mps_roots.oct])
AC_SUBST(OCTAVE_MODULES_PROGRAMS)

# Check for the presence of "mex", the Matlab compiler
# and if it's present, define MATLAB_MODULES_PROGRAMS
# to be mps_roots.mexglx if we are on a 32 bit pc
# and mps_roots.mexglxa64 if we are on a 64 bit machine.
AC_CHECK_PROGS([MEX], [mex])

# We need to test if MEX is the matlab compiler
# or the polish version of the pdftex.
AC_MSG_CHECKING(if mex is the matlab compiler)
mex_is_matlab_output=`mex -help 2>&1 | grep MATLAB`
mex_is_matlab_retcode=$?
if [[ "$mex_is_matlab_retcode" != "0" ]]; then
  MEX=""
  AC_MSG_RESULT(no)
else
  AC_MSG_RESULT(yes)
fi

# If we have find a valid Matlab compiler then add
# MATLAB to the target of compilation
enable_matlab=no
AS_IF([test "$MEX" != ""], [enable_matlab=yes
  MATLAB_MODULES_PROGRAMS=mps_roots])
AC_SUBST(MATLAB_MODULES_PROGRAMS)

AC_SUBST(MPS_CFLAGS)
AC_SUBST(UNISOLVE_CFLAGS)
AC_SUBST(SECSOLVE_CFLAGS)
AC_SUBST(UNISOLVE_LIBS)
AC_SUBST(SECSOLVE_LIBS)

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
		 examples/octave/Makefile
                 examples/matlab/Makefile
                 src/Makefile
                 src/unisolve/Makefile
                 src/secsolve/Makefile
                 src/libmps/Makefile
                 src/tests/Makefile])


AC_OUTPUT

echo "

MPSolve configuration:

        Source location:	${srcdir}
        C compiler:		${CC}
        Debug enabled:		$enable_debug
        Check enabled:		$enable_check"

# Check Octave module
if [test x$enable_octave = xyes]; then
        echo "	Octave module:		$enable_octave (generates examples/octave/$OCTAVE_MODULES_PROGRAMS)"
else
        echo "	Octave module:		no"
fi

# Check MATLAB (TM) module
if [test x$enable_matlab = xyes]; then
        echo "	MATLAB (tm) module:	$enable_octave (generates examples/matlab/$MATLAB_MODULES_PROGRAMS.mex)"
else
        echo "	MATLAB (tm) module:	no"
fi

echo "
Type 'make' to compile MPSolve, and then make install
to install the binaries, the library and the headers
system-wide.

"
