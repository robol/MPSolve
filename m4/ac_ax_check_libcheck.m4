# -*- Autoconf -*-
# @stelib-file@
# serial 2 ax_check_libcheck.m4

#
# Copyright (C) 2010 Stefano Lattarini.
# Written by Stefano Lattarini <stefano.lattarini@gmail.com>.
# This file is part of C Ratfor.
#
# Define macro `AX_CHECK_LIBCHECK' which looks for the `check' library
# (unit testing library for C) and related header files.
#
# Copying and distribution of this file, with or without modification, are
# permitted in any medium without royalty provided the copyright notice and
# this notice are preserved.
#


#######################
##  Internal macros  ##
#######################

AC_DEFUN([_AX_LIBCHECK_HELP_STRING],
    [AS_HELP_STRING(
        [--with-libcheck@<:@=DIR@:>@ @<:@default=auto@:>@],
        [Look for the `check' library and headers in the directory DIR,
         and issue a fatal error if they cannot be found or used.
         DIR can have the special values "no" (no search for libcheck is
         performed), "yes" (the library is searched in the default
         locations, a fatal error is issued if it isn't found) and "auto"
         (like "yes", but issue a warning instead of a fatal error if the
         library is not found).])])

AC_DEFUN([_AX_LIBCHECK_NOT_FOUND_ERROR],
    [AC_MSG_ERROR([m4_do([library `check' not found or not usable, ],
                         [but the user explicitly required it])])])

AC_DEFUN([_AX_LIBCHECK_NOT_FOUND_WARN],
    [AC_MSG_WARN([m4_do([library `check' not found or not usable; ],
                        [unit tests will be skipped])])])

AC_DEFUN([_AX_HEADER_CHECK_H],
    [AC_CHECK_HEADER([check.h],
                     [ax_have_check_hdr=yes],
                     [ax_have_check_hdr=no])])

AC_DEFUN([_AX_CHECK_LIBCHECK],
    [m4_pushdef([$0_testfunc], [suite_create])
    # check for libcheck (no tongue-in-cheek)
    AC_SEARCH_LIBS($0_testfunc, [check libcheck])
    AS_CASE([$ac_cv_search_]$0_testfunc,
        [none*required],
            [ax_have_check_lib=yes],
        [no],
            [ax_have_check_lib=no],
        dnl* default
            [ax_have_check_lib=yes]
            [AS_VAR_APPEND([ax_libcheck_ldflags],
                           [" ${ac_cv_search_]$0_testfunc[}"])])
    m4_popdef([$0_testfunc])])


#######################
##  Public macro(s)  ##
#######################

#
# AX_CHECK_LIBCHECK
# -----------------
#
# Look for the `check' library (unit testing library for C) and related
# header files.
#
# If the library is accessible:
#  - define (with AC_DEFINE) the preprocessor symbol `HAVE_CHECK' to `1';
#  - define (with AC_SUBST) `LIBCHECK_CPPFLAGS' to the preprocessor's
#    flag(s) necessary to include the library header file(s);
#  - define (with AC_SUBST) `LIBCHECK_LDFLAGS' to the compiler's flag(s)
#    necessary to link to the library.
#
# If the library is *not* accessible:
#  - define (with AC_DEFINE) the preprocessor symbol `HAVE_CHECK' to `0';
#  - define (with AC_SUBST) `LIBCHECK_CPPFLAGS' and `LIBCHECK_LDFLAGS' to
#    the empty string.
#
# In particular, note that if the preprocessor symbol `HAVE_LIBCHECK' is
# defined to `1', the user code can `#include <check.h>' without having to
# look at the `HAVE_CHECK_H' symbol.
#
# The search for theh `check' library will honour the `--with-libcheck'
# configure option (full help for this option string is set by the internal
# macro `_AX_LIBCHECK_HELP_STRING' above).
#
AC_DEFUN_ONCE([AX_CHECK_LIBCHECK],
    [dnl* add `--with-libcheck' option to configure
    AC_ARG_WITH([libcheck], [_AX_LIBCHECK_HELP_STRING()],
                [:], [with_libcheck=auto])
    ax_libcheck_cppflags=''
    ax_libcheck_ldflags=''
    AC_MSG_CHECKING([whether/how to use library `check'])
    AS_CASE([$with_libcheck],
        [no],   [AC_MSG_RESULT([disabled])],
        [yes],  [AC_MSG_RESULT([mandatory, search in default paths])],
        [auto], [AC_MSG_RESULT([optional, search in default paths])],
                [AC_MSG_RESULT([mandatory, look in $with_libcheck])
                 ax_libcheck_cppflags="-I$with_libcheck/include"
                 ax_libcheck_ldflags="-L$with_libcheck/lib"])
    AS_IF([test x"$with_libcheck" != x"no"],
          [dnl* Save original values of CPPFLAGS and LDFLAGS.
           _ax_check_libcheck_save_CPPFLAGS=$CPPFLAGS
           _ax_check_libcheck_save_LDFLAGS=$LDFLAGS
           dnl* Temporarly redefine CPPFLAGS and LDFLAGS.
           AS_VAR_APPEND([CPPFLAGS], [" $ax_libcheck_cppflags"])
           AS_VAR_APPEND([LDFLAGS], [" $ax_libcheck_ldflags"])
           dnl* Check for <check.h> header, then check for the `check'
           dnl* library itself.
           _AX_HEADER_CHECK_H
           _AX_CHECK_LIBCHECK
           dnl* Reset CPPFLAGS and LDFLAGS.
           CPPFLAGS=$_ax_check_libcheck_save_CPPFLAGS
           LDFLAGS=$_ax_check_libcheck_save_LDFLAGS
           AS_UNSET([_ax_check_libcheck_save_CPPFLAGS])
           AS_UNSET([_ax_check_libcheck_save_LDFLAGS])])
    AS_IF([test ":$ax_have_check_hdr:$ax_have_check_lib:" = :yes:yes:],
          [ax_have_check=yes],
          [ax_have_check=no])
    AS_IF([test x"$ax_have_check" = x"yes"],
          [AC_DEFINE([HAVE_CHECK], [1],
                     [Define to 1 if the `check' library is available])],
          [ax_libcheck_cppflags='' ax_libcheck_ldflags=''
          AS_CASE([$with_libcheck],
                  [auto], [_AX_LIBCHECK_NOT_FOUND_WARN()],
                    [no], [:],
                          [_AX_LIBCHECK_NOT_FOUND_ERROR()])])
    AC_SUBST([LIBCHECK_LDFLAGS], [$ax_libcheck_ldflags])
    AC_SUBST([LIBCHECK_CPPLAGS], [$ax_libcheck_cppflags])])

# vim: ft=m4 ts=4 sw=4 et
