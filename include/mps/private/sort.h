/*
 * This file is part of MPSolve 3.1.7
 *
 * Copyright (C) 2001-2019, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Leonardo Robol <leonardo.robol@sns.it>
 */

/**
 * @file
 * @brief Implementation of sorting routines for MPSolve.
 */

#ifndef MPS_SORT_H_
#define MPS_SORT_H_

MPS_BEGIN_DECLS

void mps_fsort (mps_context * s);
void mps_dsort (mps_context * s);
void mps_msort (mps_context * s);

MPS_END_DECLS

#endif /* MPS_SORT_H_ */


