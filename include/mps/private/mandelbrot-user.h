/*
 * This file is part of MPSolve 3.1.5
 *
 * Copyright (C) 2001-2015, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Leonardo Robol <leonardo.robol@sns.it>
 */

/**
 * @file
 * @brief
 */

#ifndef MPS_MANDELBROT_USER_H_
#define MPS_MANDELBROT_USER_H_

MPS_BEGIN_DECLS

void mps_fnewton_usr (mps_context * st, mps_polynomial * poly, mps_approximation * root, cplx_t corr);
void mps_dnewton_usr (mps_context * st, mps_polynomial * poly, mps_approximation * root, cdpe_t corr);
void mps_mnewton_usr (mps_context * st, mps_polynomial * poly, mps_approximation * root, mpc_t corr, long int wp);
mps_boolean mps_feval_usr (mps_context * ctx, mps_polynomial * p, cplx_t x, cplx_t value, double * error);
mps_boolean mps_deval_usr (mps_context * ctx, mps_polynomial * p, cdpe_t x, cdpe_t value, rdpe_t error);
mps_boolean mps_meval_usr (mps_context * ctx, mps_polynomial * p, mpc_t x, mpc_t value, rdpe_t error);

MPS_END_DECLS

#endif /* MPS_MANDELBROT_USER_H_ */

