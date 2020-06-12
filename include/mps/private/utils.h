/*
 * This file is part of MPSolve 3.2.1
 *
 * Copyright (C) 2001-2020, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Leonardo Robol <leonardo.robol@unipi.it>
 */

#ifndef _MPS_UTILS_H
#define _MPS_UTILS_H

MPS_BEGIN_DECLS

#include <mps/mps.h>

/* This function is currently implemented in parser.c for historical reasons. */
char * build_equivalent_rational_string (mps_context * ctx, const char * line,
                                         long int * exponent, int * sign);

char * mps_utils_strip_string (mps_context * ctx, const char * input);
char * mps_utils_build_equivalent_rational_string (mps_context * ctx,
                                                   const char * input);

/* functions in newton.c */
int mps_intlog2 (int n);

/* function in strndup.c */
#ifndef HAVE_STRNDUP
char * mps_strndup (const char * source, size_t n);
#else
#define mps_strndup strndup
#endif

MPS_END_DECLS

#endif
