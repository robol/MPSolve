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
  MPS_DEBUG_THIS_CALL (ctx);

  int i;
  mpq_t tmp_r, tmp_i;

  if (! MPS_IS_MONOMIAL_POLY (poly))
    {
      MPS_DEBUG_WITH_INFO (ctx, "Falling back to starting_strategy: MPS_STARTING_STRATEGY_DEFAULT, since the"
			   "input is not a monomial poly");

      mps_context_select_starting_strategy (ctx, MPS_STARTING_STRATEGY_DEFAULT);
      mps_polynomial_fstart (ctx, poly);
      return;
    }

  if (poly->degree < 200)
    {
      (*poly->fstart)(ctx, poly);
      return;
    }

  mps_monomial_poly * mp = MPS_MONOMIAL_POLY (poly);

  /* Compute the starting radii by the Newton polygon. This function will store the result in
   * ctx->partitioning and ctx->n_radii. We need some refactoring to make sure that the results
   * will be stored locally and we will be able to recursively apply this strategy. */
  mps_starting_configuration c = 
    mps_fcompute_starting_radii (ctx, poly->degree, NULL, 0.0, 0.0, ctx->eps_out, mp->fap);

  /* Find out a middle point in the Newton polygon. */
  int middle = 0; 
  int cursor = 0;
  while (middle < ctx->n / 2 && cursor < c.n_radii)
    {
      middle = c.partitioning[++cursor];
    }

  if (middle == 0 || middle == ctx->n)
    middle = ctx->n / 2; 

  /* Split in two polynomials */
  mps_context * rctx = mps_context_new ();
  
  mps_monomial_poly *left = mps_monomial_poly_new (rctx, middle);
  mps_monomial_poly *right = mps_monomial_poly_new (rctx, poly->degree - middle);

  /* Copy some configuration from the originating context */
  mps_context_add_debug_domain (rctx, ctx->debug_level);
  mps_context_select_algorithm (rctx, ctx->algorithm);
  mps_context_select_starting_strategy (rctx, MPS_STARTING_STRATEGY_RECURSIVE);

  MPS_DEBUG_WITH_INFO (ctx, "Divided the polynomial into two polynomials of degree %d and %d", 
		       middle, poly->degree - middle);

  mpq_init (tmp_r);
  mpq_init (tmp_i);

  /* Fill in the first polynomial */
  for (i = 0; i <= middle; i++)
    {
      mps_monomial_poly_set_coefficient_q (rctx, left, i, 
					   mp->initial_mqp_r[i],
					   mp->initial_mqp_i[i]);
    }

  mpq_neg (tmp_r, mp->initial_mqp_r[middle]);
  mpq_neg (tmp_i, mp->initial_mqp_i[middle]);
  mps_monomial_poly_set_coefficient_q (rctx, right, 0, tmp_r, tmp_i);

  for (i = middle + 1; i <= poly->degree; i++)
    {
      mps_monomial_poly_set_coefficient_q (rctx, right, i - middle, 
					   mp->initial_mqp_r[i],
					   mp->initial_mqp_i[i]);
    }
  
  mpq_clear (tmp_r);
  mpq_clear (tmp_i);

  /* Recursively solve both polynomials. */
  mps_context_set_input_poly (rctx, MPS_POLYNOMIAL (left));
  mps_mpsolve (rctx);
  mps_approximation ** left_approximations = mps_context_get_approximations (rctx);

  mps_context_set_input_poly (rctx, MPS_POLYNOMIAL (right));
  mps_mpsolve (rctx);
  mps_approximation ** right_approximations = mps_context_get_approximations (rctx);

  /* Use these approximations as starting points for the original polynomial. */
  mps_context_set_input_poly (ctx, poly);
  for (i = 0; i < poly->degree; i++)
    {
      mps_approximation * appr = (i < middle) ? left_approximations[i] : 
	right_approximations[i - middle];
      cplx_set (ctx->root[i]->fvalue, appr->fvalue);      
      free (appr);
    }

  free (left_approximations);
  free (right_approximations);

  mps_monomial_poly_free (rctx, MPS_POLYNOMIAL (left));
  mps_monomial_poly_free (rctx, MPS_POLYNOMIAL (right));
  mps_context_free (rctx);

  mps_starting_configuration_clear (ctx, &c);

  MPS_DEBUG_WITH_INFO (ctx, "Completed recursive repositioning of polynomials");
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
