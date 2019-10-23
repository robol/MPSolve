/*
 * This file is part of MPSolve 3.1.8
 *
 * Copyright (C) 2001-2019, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Leonardo Robol <leonardo.robol@sns.it>
 */

/**
 * @file
 *
 * @brief This file contains the entry points of the various algorithms implemented in
 * MPSolve.
 */

#ifndef MPS_ALGORITHMS_H_
#define MPS_ALGORITHMS_H_

MPS_BEGIN_DECLS

/* This is the standard MPSolve algorithm used also in MPSolve 2.2
 * The version implemented here is modified to use the new framework. */
void mps_standard_mpsolve (mps_context * s);

/* This is the new algorithm inserted in MPSolve 3.0, that uses secular
 * equations to solve polynomial ones. */
void mps_secular_ga_mpsolve (mps_context * s);

MPS_END_DECLS

#endif /* endif MPS_ALGORITHMS_H_ */
