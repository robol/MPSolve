#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.67])
AC_INIT([MPSolve], [3.0.0], [robol@poisson.phc.unipi.it])
AC_CONFIG_SRCDIR([Makefile.am])
AC_CONFIG_HEADERS([config.h])

AC_LIBTOOL_WIN32_DLL
AC_CONFIG_MACRO_DIR([m4])

LT_INIT

AM_INIT_AUTOMAKE([dist-bzip2])

m4_ifdef([AM_SILENT_RULES],[AM_SILENT_RULES([yes])])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AM_PROG_CC_C_O
AC_PROG_CC_C99

# Check if we are building on Windows, so we can add additional
# hacks to make MPSolve work there. 
AC_MSG_CHECKING([for Win32 platform])
case "$host" in
  *-*-mingw*|*-*-cygwin*)
    platform_win32=yes
    ;;
  *)
    platform_win32=no
    ;;
esac
AC_MSG_RESULT([$platform_win32])

# That's make detecting Windows real easy...
if test x$platform_win32 = xyes; then
  CFLAGS="$CFLAGS -D__WINDOWS"
fi

# Add GMP_CFLAGS and GMP_LIBS to CFLAGS and LDFLAGS in case they are set.
CFLAGS="$GMP_CFLAGS $CFLAGS"
LDFLAGS="$GMP_LIBS $LDFLAGS"

GENERAL_CFLAGS="-D_REENTRANT -Wstrict-prototypes -Wall $CFLAGS"
# AS_IF([test x$platform_win32 != "yes"], [
#  GENERAL_LIBS="$GENERAL_LIBS -no-undefined -Wl,--export-all-symbols
#])

# Determine CFLAGS based on the options that the user has specified.
# The supported build type are
#
# --without-debug Do not code debug prints (that is not the default)
# 
AC_ARG_ENABLE([debug],
	AS_HELP_STRING([--disable-debug], [Disable debug prints in the code]),
	[enable_debug=$enableval have_debug=$enableval],
	[enable_debug=yes have_debug=yes])
# AC_MSG_RESULT([$enable_debug])

# Select CFLAGS based on debug enabled or not
AS_IF([test x$enable_debug != xno], [
	GENERAL_CFLAGS="$GENERAL_CFLAGS -DNICE_DEBUG"
])

# Determine if examples are desired. If that's the case check for the library
# that are installed on the system and see which examples can be added. 
AC_ARG_ENABLE([examples],
	AS_HELP_STRING([--disable-examples], [Disable additional examples and modules bundled with MPSolve]),
	[build_examples=$enableval have_examples=$enableval],
	[build_examples=yes have_examples=yes])

# Checks for libraries.
AC_CHECK_LIB([gmp], [__gmpz_init2])
AC_CHECK_HEADERS([gmp.h], [mpsolve_found_gmp_headers=yes; break;])

AC_CHECK_LIB([tcmalloc_minimal], [malloc])

AC_CHECK_LIB([m], [pow])

AX_PTHREAD([], [
  AC_MSG_ERROR([Pthread library not found. Please set PTHREAD_CFLAGS and PTHREAD_LIBS correctly for your setup])
])

AC_SUBST(PTHREAD_CFLAGS)
AC_SUBST(PTHREAD_LIBS)
GENERAL_CFLAGS="$GENERAL_CFLAGS $PTHREAD_CFLAGS"
GENERAL_LIBS="$PTHREAD_LIBS $GENERAL_LIBS $LDFLAGS"
CC="$PTHREAD_CC"

# Errors if libraries are not found, in this case GMP and Pthreads.
AS_IF([test "x$mpsolve_found_gmp_headers" != "xyes"],
      [AC_MSG_ERROR([Unable to find gmp headers, pleas install libgmp-dev])])

# Check for check B-)
PKG_CHECK_MODULES([CHECK], [check >= 0.9.4], [have_check=yes], [have_check=no])
CHECK_CFLAGS="$GENERAL_CFLAGS"
CHECK_LIBS="$GENERAL_LIBS"
CHECK_PROGRAMS="test_mpsolve check_dpe check_secsolve check_unisolve check_multithread check_cluster"

AC_SUBST(CHECK_PROGRAMS)
AC_SUBST(CHECK_CFLAGS)
AC_SUBST(CHECK_LIBS)

# Check the existance of mkoctfile to compile octave 
# module
AS_IF([test x$build_examples == xyes], [
	AC_CHECK_PROGS([MKOCTFILE], [mkoctfile])
	enable_octave=no
	AS_IF([test "$MKOCTFILE" != ""], [enable_octave=yes OCTAVE_MODULES_PROGRAMS=mps_roots.oct])
	AC_SUBST(OCTAVE_MODULES_PROGRAMS)
])

# Check for the presence of "mex", the Matlab compiler
# and if it's present, define MATLAB_MODULES_PROGRAMS
# to be mps_roots.mexglx if we are on a 32 bit pc
# and mps_roots.mexglxa64 if we are on a 64 bit machine.
AS_IF([test x$build_examples == xyes], [
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
])

# Check for GTK libraries to enable the example that is given
# in examples/gtk
AS_IF([test x$build_examples == xyes], [
	ROOTS_OF_UNITY_PROGRAMS=roots_of_unity
	PKG_CHECK_MODULES([GTK], [gtk+-3.0], [], [
		PKG_CHECK_MODULES([GTK], [gtk+-2.0], [], [
			ROOTS_OF_UNITY_PROGRAMS=""
		])
	])

	AC_SUBST([GTK_CFLAGS])
	AC_SUBST([GTK_LIBS])
	AC_SUBST([GTK_FOUND])
	AC_SUBST([ROOTS_OF_UNITY_PROGRAMS])
])
			

# Check for doxygen to generate documentation
AC_CHECK_PROGS([DOXYGEN], [doxygen])

if [test -z "$DOXYGEN"]; then
  AC_MSG_WARN([Doxygen not found, documentation will not be generated])
fi
AM_CONDITIONAL([HAVE_DOXYGEN], [test -n "$DOXYGEN"])

# AS_IF([test "x$mpsolve_found_pthread_headers" != "xyes"],
#	[AC_MSG_ERROR([Unable to find pthreads headers, please install libpthreads-dev])])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_TYPE_SIZE_T

# Check for getline() function
AC_CHECK_FUNCS(getline)

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([memset modf pow sqrt strchr])

MPS_CFLAGS="$GENERAL_CFLAGS -D_MPS_PRIVATE"
MPS_LDFLAGS="$GENERAL_LIBS"
MPSOLVE_CFLAGS="$GENERAL_CFLAGS"
MPSOLVE_LIBS="$GENERAL_LIBS"

if test "$platform_win32" = yes; then
  MPS_LDFLAGS="$MPS_LDFLAGS -no-undefined -Wl,--export-all-symbols -Wl,--output-def,.libs/libmps-3.dll.def"
fi

AC_SUBST(MPS_LDFLAGS)
AC_SUBST(MPS_CFLAGS)
AC_SUBST(MPSOLVE_CFLAGS)
AC_SUBST(MPSOLVE_LIBS)

AC_CONFIG_FILES([Makefile
                 examples/Makefile
		 examples/octave/Makefile
                 examples/matlab/Makefile
		 examples/gtk/Makefile
                 src/Makefile
		 src/mpsolve/Makefile
                 src/libmps/Makefile
                 src/tests/Makefile
		 doc/Makefile
		 doc/Doxyfile])

AC_OUTPUT

echo "

MPSolve configuration:

        Source location:	${srcdir}
        C compiler:		${CC}
	CFLAGS:			${MPS_CFLAGS}
        LDFLAGS:                ${MPS_LDFLAGS}
        Debug enabled:		$enable_debug
        Check enabled:		$have_check"

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
