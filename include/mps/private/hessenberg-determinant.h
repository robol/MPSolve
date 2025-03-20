/*
 * This file is part of MPSolve 3.2.2
 *
 * Copyright (C) 2001-2020, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Leonardo Robol <leonardo.robol@unipi.it>
 */

/**
 * @file
 *
 * @brief Implementation of determinant computation for Hessenberg matrices.
 */

#ifndef MPS_HESSENBERG_DETERMINANT_H_
#define MPS_HESSENBERG_DETERMINANT_H_

#include <mps/mps.h>

MPS_BEGIN_DECLS

void mps_fhessenberg_determinant (mps_context * ctx, cplx_t * hessenberg_matrix, size_t n, cplx_t output,
				  long int * exponent);
void mps_fhessenberg_shifted_determinant (mps_context * ctx, cplx_t * hessenberg_matrix, 
					  const cplx_t shift, size_t n, cplx_t output, 
					  long int * exponent);

void mps_dhessenberg_determinant (mps_context * ctx, cdpe_t * hessenberg_matrix, size_t n, cdpe_t output);
void mps_dhessenberg_shifted_determinant (mps_context * ctx, cdpe_t * hessenberg_matrix, const cdpe_t shift,
                                          size_t n, cdpe_t output);

void mps_mhessenberg_determinant (mps_context * ctx, mpc_t * hessenberg_matrix, size_t n,
                                  mpc_t output, rdpe_t error);
void mps_mhessenberg_shifted_determinant (mps_context * ctx, mpc_t * hessenberg_matrix, mpc_t shift,
                                          size_t n, mpc_t output, rdpe_t error);

MPS_END_DECLS

#endif /* endif MPS_HESSENBERG_DETERMINANT */

