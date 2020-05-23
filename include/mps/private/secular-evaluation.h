/*
 * This file is part of MPSolve 3.1.9
 *
 * Copyright (C) 2001-2020, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Leonardo Robol <leonardo.robol@unipi.it>
 */

/**
 * @file
 * @brief Evaluation of secular equations.
 */

#include <mps/mps.h>

#ifndef MPS_SECULAR_EVALUATION_H_
#define MPS_SECULAR_EVALUATION_H_

MPS_BEGIN_DECLS

/* Functions in secular-evaluation.c */
mps_boolean mps_secular_feval (mps_context * s, mps_polynomial * p, cplx_t x, cplx_t value);
mps_boolean mps_secular_feval_with_error (mps_context * s, mps_polynomial * p,
                                          cplx_t x, cplx_t value, double * error);
mps_boolean mps_secular_deval (mps_context * s, mps_polynomial * p, cdpe_t x, cdpe_t value);
mps_boolean mps_secular_deval_derivative (mps_context * s, mps_polynomial * p, cdpe_t x, cdpe_t value);
mps_boolean mps_secular_deval_with_error (mps_context * s, mps_polynomial * p,
                                          cdpe_t x, cdpe_t value, rdpe_t error);
mps_boolean mps_secular_meval (mps_context * s, mps_polynomial * p, mpc_t x, mpc_t value);
mps_boolean mps_secular_meval_with_error (mps_context * s, mps_polynomial * p,
                                          mpc_t x, mpc_t value, rdpe_t error);
mps_boolean mps_secular_feval_derivative (mps_context * s, mps_polynomial * p, cplx_t x, cplx_t value);

MPS_END_DECLS

#endif /* MPS_SECULAR_EVALUATION_H_ */

