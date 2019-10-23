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
 * @brief Implementation of radius computation.
 */

#ifndef MPS_RADII_H_
#define MPS_RADII_H_

MPS_BEGIN_DECLS

/* Functions in general-radius.c */
void mps_fradii (mps_context * s, mps_polynomial * p, double * fradii);
void mps_dradii (mps_context * s, mps_polynomial * p, rdpe_t * dradii);
void mps_mradii (mps_context * s, mps_polynomial * p, rdpe_t * dradii);

/* Functions in monomial-radius.c */
void mps_monomial_fradii (mps_context * s, double * fradii);
void mps_monomial_dradii (mps_context * s, rdpe_t * dradii);
void mps_monomial_mradii (mps_context * s, rdpe_t * dradii);

MPS_END_DECLS

#endif /* MPS_RADII_H_ */
