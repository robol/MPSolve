/*
 * This file is part of MPSolve 3.1.5
 *
 * Copyright (C) 2001-2014, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Leonardo Robol <robol@mail.dm.unipi.it>
 */

/**
 * @file
 * @brief Horner implementation for Monomial polynomials. 
 */

#ifndef MPS_HORNER_H_
#define MPS_HORNER_H_

MPS_BEGIN_DECLS

/* These two routines are implemented in newton.c */
void mps_parhorner (mps_context * st, int n, mpc_t x, mpc_t p[],
                    mps_boolean b[], mpc_t s, int n_thread);
void mps_aparhorner (mps_context * st, int n, rdpe_t x, rdpe_t p[],
                     mps_boolean b[], rdpe_t s, int n_thread);

/* The following routines are implemented in newton.c */
void mps_fhorner (mps_context * s, mps_monomial_poly * p, cplx_t x, cplx_t value);
void mps_fhorner_with_error (mps_context * s, mps_monomial_poly * p, cplx_t x, 
			     cplx_t value, double * relative_error);
void mps_dhorner (mps_context * s, mps_monomial_poly * p, cdpe_t x, cdpe_t value);
void mps_dhorner_with_error (mps_context * s, mps_monomial_poly * p, cdpe_t x, cdpe_t value, rdpe_t relative_error);
void mps_mhorner (mps_context * s, mps_monomial_poly * p, mpc_t x, mpc_t value);
void mps_mhorner_with_error (mps_context * s, mps_monomial_poly * p, 
			     mpc_t x, mpc_t value, rdpe_t relative_error, long int wp);
void mps_mhorner_with_error2 (mps_context * s, mps_monomial_poly * p, mpc_t x, 
			      mpc_t value, rdpe_t relative_error, long int wp);

MPS_END_DECLS

#endif /* endif MPS_HORNER_H_ */
