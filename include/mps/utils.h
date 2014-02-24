/*
 * This file is part of MPSolve 3.1.5
 *
 * Copyright (C) 2001-2014, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Leonardo Robol <robol@mail.dm.unipi.it>
 */

#ifndef _MPS_UTILS_H
#define _MPS_UTILS_H

#ifdef __cplusplus
extern "C" {
#endif

#include <mps/mps.h>


char * mps_utils_strip_string (mps_context * ctx, const char * input);
char * mps_utils_build_equivalent_rational_string (mps_context * ctx,
                                                   const char * input);


#ifdef __cplusplus
}
#endif

#endif
