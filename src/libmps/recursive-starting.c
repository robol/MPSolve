/*
 * This file is part of MPSolve 3.1.5
 *
 * Copyright (C) 2001-2014, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Leonardo Robol <leonardo.robol@sns.it>
 */

#include <mps/mps.h>

/**
 * @brief Select appropriate starting point for the approximation of the roots
 * of the given polynomial by applying a divide-and-conquer strategy described
 * in {TODO: Reference missing}. 
 *
 * @param ctx The current mps_context. 
 * @param poly The polynomial whose roots should be approximated.
 */
void
mps_recursive_fstart (mps_context * ctx, mps_polynomial * poly)
{
  mps_monomial_poly * mp = MPS_MONOMIAL_POLY (poly);

  /* Compute the starting radii by the Newton polygon. This function will store the result in
   * ctx->partitioning and ctx->n_radii. We need some refactoring to make sure that the results
   * will be stored locally and we will be able to recursively apply this strategy. */
  mps_fcompute_starting_radii (ctx, poly->degree, NULL, 0.0, 0.0, ctx->eps_out, mp->fap);
}

/**
 * @brief Select appropriate starting point for the approximation of the roots
 * of the given polynomial by applying a divide-and-conquer strategy described
 * in {TODO: Reference missing}. 
 *
 * @param ctx The current mps_context. 
 * @param poly The polynomial whose roots should be approximated.
 */
void
mps_recursive_dstart (mps_context * ctx, mps_polynomial * poly)
{
  
}

/**
 * @brief Select appropriate starting point for the approximation of the roots
 * of the given polynomial by applying a divide-and-conquer strategy described
 * in {TODO: Reference missing}. 
 *
 * @param ctx The current mps_context. 
 * @param poly The polynomial whose roots should be approximated.
 */
void
mps_recursive_mstart (mps_context * ctx, mps_polynomial * poly)
{
  
}
