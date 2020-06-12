/*
 * This file is part of MPSolve 3.2.1
 *
 * Copyright (C) 2001-2020, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Leonardo Robol <leonardo.robol@unipi.it>
 */

/**
 * @file
 * @brief Routines that check the emptyness of the intersection of several sets.
 */

#ifndef MPS_TOUCH_H_
#define MPS_TOUCH_H_

MPS_BEGIN_DECLS

/* functions in touch.c */
mps_boolean mps_ftouchnwt (mps_context * s, double * frad, int n, int i, int j);
mps_boolean mps_dtouchnwt (mps_context * s, rdpe_t * drad, int n, int i, int j);
mps_boolean mps_mtouchnwt (mps_context * s, rdpe_t * drad, int n, int i, int j);
mps_boolean mps_ftouchreal (mps_context * s, int n, int i);
mps_boolean mps_dtouchreal (mps_context * s, int n, int i);
mps_boolean mps_mtouchreal (mps_context * s, int n, int i);
mps_boolean mps_ftouchimag (mps_context * s, int n, int i);
mps_boolean mps_dtouchimag (mps_context * s, int n, int i);
mps_boolean mps_mtouchimag (mps_context * s, int n, int i);
mps_boolean mps_ftouchunit (mps_context * s, int n, int i);
mps_boolean mps_dtouchunit (mps_context * s, int n, int i);
mps_boolean mps_mtouchunit (mps_context * s, int n, int i);

/* functions in validation.c */
void mps_validate_inclusions (mps_context * ctx);

MPS_END_DECLS

#endif /* MPS_TOUCH_H_ */

