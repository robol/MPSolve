#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.67])
AC_INIT([MPSolve], [3.1.2], [robol@poisson.phc.unipi.it])
AC_CONFIG_SRCDIR([Makefile.am])
AC_CONFIG_HEADERS([config.h])

AC_CONFIG_MACRO_DIR([m4])

LT_INIT

AM_INIT_AUTOMAKE([dist-bzip2])

m4_ifdef([AM_SILENT_RULES],[AM_SILENT_RULES([yes])])

#
# ===== GENERAL SECTION =====
#
# This section contains checks for compilers, basic headers files, ...
#

# Checks for needed compilers
AC_PROG_CXX
AC_PROG_CC
AM_PROG_CC_C_O
AC_PROG_CC_C99

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

AC_SUBST(GMP_CFLAGS)
AC_SUBST(GMP_LIBS)

GENERAL_CFLAGS="-D_REENTRANT -Wstrict-prototypes -Wall $CFLAGS"

#
# ===== OPTION PARSING =====
#
# Add options to the configure script so that the user can select the features
# he's interested in. 
#

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

# Add another option, enable/disable UI. 
AC_ARG_ENABLE([ui],
	AS_HELP_STRING([--disable-ui], [Don't build the graphical user interface for MPSolve]),
	[build_ui=$enableval have_ui=$enableval],
	[build_ui=yes have_ui=yes])
	
# Determine if the Graphical debugger must be integrated in MPSolve
AC_ARG_ENABLE([graphical-debugger],
    AS_HELP_STRING([--disable-graphical-debugger], [Don't bundle the GTK graphical debugger in MPSolve]),
    [build_graphical_debugger=$enableval have_graphical_debugger=$enableval],
    [build_graphical_debugger=yes have_graphical_debugger=yes])

# 
# ===== REQUIRED LIBRARIES =====
#
# This section checks for the libraries that are needed to build MPSolve. These are
# required libraries and if these are not found on the system there's no chance that
# we can continue ignoring the error. 
#

# We need the multiprecicsion GMP library. 
AC_CHECK_LIB([gmp], [__gmpz_init2])
AC_CHECK_HEADERS([gmp.h], [mpsolve_found_gmp_headers=yes; break;])

# Raise an error if GMP is not correctly installed on this system.
AS_IF([test "x$mpsolve_found_gmp_headers" != "xyes"],
      [AC_MSG_ERROR([Unable to find gmp headers, pleas install libgmp-dev])])

# The general C math library. 
AC_CHECK_LIB([m], [pow])

# At the moment MPSolve cannot be built without the pthread library for POSIX threads.
# In principle it should be possible to use an abstraction layer to support also Windows
# cricital section, but this have not be done yet. Note, however, that is still possible
# to run MPSolve on Windows by using pthread compatilibity layer that are freely available
# on the Internet. 
AX_PTHREAD([], [
  AC_MSG_ERROR([Pthread library not found. Please set PTHREAD_CFLAGS and PTHREAD_LIBS correctly for your setup])
])

# Set the correct PTHREAD flags and, if needed, change the compiler to one that is pthread-enabled.
AC_SUBST(PTHREAD_CFLAGS)
AC_SUBST(PTHREAD_LIBS)
GENERAL_CFLAGS="$GENERAL_CFLAGS $PTHREAD_CFLAGS"
GENERAL_LIBS="$PTHREAD_LIBS $GENERAL_LIBS $LDFLAGS"
CC="$PTHREAD_CC"

AC_SUBST(CC)

#
# ===== OPTIONAL LIBRARIES =====
#
# This category contain libraries that can be used by MPSolve (to add features or be
# more efficient) but are not strictly required to build the software. 

# Check for the MPFR multiprecision floating point library. This could be, in principle,
# an alternative to GMP mpf_t type. This is still not implemented in MPSolve, though. 
AC_CHECK_LIB([mpfr], [mpfr_neg])
AC_CHECK_HEADERS([mpfr.h])

# This is a custom malloc implementation with thread-caching capabilities. Using this library
# will make MPSolve much faster. 
AC_CHECK_LIB([tcmalloc_minimal], [malloc])

# Check for the check module, so MPSolve tests can be run. 
PKG_CHECK_MODULES([CHECK], [check >= 0.9.4], [have_check=yes], [have_check=no])

# If check was fuond enable the tests that need to be built
AS_IF([test "x$have_check" = "xyes"], [ 
  CHECK_PROGRAMS="check_dpe check_secsolve check_unisolve check_multithread check_cluster"
])
  
AC_SUBST(CHECK_PROGRAMS)
AC_SUBST(CHECK_CFLAGS)
AC_SUBST(CHECK_LIBS)

# Check for the GTK libraries on this system. We can use both GTK 2.0 and GTK 3.0
have_gtk="1"
PKG_CHECK_MODULES([GTK], [gtk+-3.0], [], [
    PKG_CHECK_MODULES([GTK], [gtk+-2.0], [], [
        have_gtk="0"
    ])
])

AC_DEFINE([HAVE_GTK], $have_gtk, [Defined if GTK is available])
AC_SUBST(GTK_CFLAGS)
AC_SUBST(GTK_LIBS)

# Check that we have also the cairo graphic library. This will be already bundled
# in GTK3 but is not to be taken for granted in GTK2. 
PKG_CHECK_MODULES([CAIRO], [cairo], [have_cairo="1"], [have_cairo="0"])

AC_DEFINE([HAVE_CAIRO], $have_cairo, [Defined if CAIRO is available])
AC_SUBST(CAIRO_CFLAGS)
AC_SUBST(CAIRO_LIBS)

AS_IF([test x$have_gtk == x1 && test x$have_graphical_debugger == xyes && test x$have_cairo == x1], [
  AC_DEFINE([HAVE_GRAPHICAL_DEBUGGER], 1, [Defined if the graphical debugger should be bundled in MPSolve])
])

#
# ===== EXAMPLES ======
#
# This section contains the checks that are needed in case the user wants to build 
# the examples contained in the examples/ directory
#

# Check the existance of mkoctfile to compile octave modules.
AS_IF([test x$build_examples == xyes], [
	AC_CHECK_PROGS([MKOCTFILE], [mkoctfile])
	enable_octave=no
	
	# Detect Octave version. We are not able to build the Octave module
	# if Octave is older than 3.0
	AC_CHECK_PROGS([OCTAVE], [octave])
	AC_MSG_CHECKING([if Octave is newer than 3.0])
	OCTAVE_VERSION=`echo "disp(OCTAVE_VERSION)" | $OCTAVE -qf`	
	AX_COMPARE_VERSION([$OCTAVE_VERSION], [ge], [3.0.0], [
			AC_MSG_RESULT([yes])
			
			# If we have Octave that is sufficiently recent, we can enable the Octave module
        	AS_IF([test "$MKOCTFILE" != ""], [enable_octave=yes OCTAVE_MODULES_PROGRAMS=mps_roots.oct])			
		], [
			AC_MSG_RESULT([no])
		])

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
AS_IF([test x$build_examples == xyes && test x$have_gtk == x1 && test x$have_cairo == x1], [
	ROOTS_OF_UNITY_PROGRAMS=roots_of_unity

	AC_SUBST([GTK_CFLAGS])
	AC_SUBST([GTK_LIBS])
	AC_SUBST([GTK_FOUND])
	AC_SUBST([ROOTS_OF_UNITY_PROGRAMS])
])

#
# ===== MPSOLVE GRAPHICAL INTERFACE =====
#
# Check for Qt to build the MPSolve GUI. 
#

# Search qmake and Qt related things, so we can build
# MPSolve UI. Don't build the UI if the user has manually
# disabled it. 
QT4_PROGRAMS=
AS_IF([test x$build_ui == xyes], [
	AC_CHECK_PROGS([QMAKE_QT4], [qmake-qt4 qmake])
	if test "x$QMAKE_QT4" != "x"; then
		# Check that this is really the qmake-qt4 command and
		# not the Qt3 one. 
		AC_MSG_CHECKING([if qmake is at least for Qt 4.0.0])
		QMAKE_VERSION=`$QMAKE_QT4 -version | grep "Qt version" | cut -d ' ' -f4`
		AC_SUBST([QMAKE_VERSION])
		AX_COMPARE_VERSION([$QMAKE_VERSION], [ge], [4.0.0], [
			AC_MSG_RESULT([yes])
			QT4_PROGRAMS=xmpsolve
		], [
			AC_MSG_RESULT([no])
		])
	fi

	AC_CHECK_PROGS([UIC_QT4], [uic-qt4])
	if test "x$UIC_QT4" == "x"; then
		QT4_PROGRAMS=""
	fi

	if test "x$QT4_PROGRAMS" == "x"; then
		build_ui=no
		have_ui=no
	fi
])
AC_SUBST(QT4_PROGRAMS)

#
# ===== MPSOLVE DOCUMENTATION ======
#
# Check that we have doxygen on this sytem. If that's the case, build the documentation for
# MPSolve. 

AC_CHECK_PROGS([DOXYGEN], [doxygen])

if [test -z "$DOXYGEN"]; then
  AC_MSG_WARN([Doxygen not found, documentation will not be generated])
fi
AM_CONDITIONAL([HAVE_DOXYGEN], [test -n "$DOXYGEN"])

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
    src/xmpsolve/Makefile
    src/xmpsolve/xmpsolve.pro
    src/xmpsolve/xmpsolve.desktop
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
        echo "	MATLAB (tm) module:	$enable_matlab (generates examples/matlab/$MATLAB_MODULES_PROGRAMS.mex)"
else
        echo "	MATLAB (tm) module:	no"
fi

if [test x$build_ui = xyes]; then
	echo "	Graphical UI:		$build_ui"
else
	echo "	Graphical UI:		no"
fi

echo "
Type 'make' to compile MPSolve, and then make install
to install the binaries, the library and the headers
system-wide.

"
