/*
 * This file is part of MPSolve 3.2.0
 *
 * Copyright (C) 2001-2020, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Leonardo Robol <leonardo.robol@unipi.it>
 */

/**
 * @file
 * @brief Implementation of Newton correction computation.
 */

#ifndef MPS_NEWTON_H_
#define MPS_NEWTON_H_

MPS_BEGIN_DECLS

void mps_fnewton (mps_context * st, mps_polynomial * p,
                  mps_approximation * root, cplx_t corr);
void mps_dnewton (mps_context * st, mps_polynomial * p,
                  mps_approximation * root, cdpe_t corr);
void mps_mnewton (mps_context * st, mps_polynomial * p,
                  mps_approximation * root, mpc_t corr, long int wp);

MPS_END_DECLS

#endif /* MPS_NEWTON_H_ */
