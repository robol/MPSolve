/*
 * This file is part of MPSolve 3.1.9
 *
 * Copyright (C) 2001-2020, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Leonardo Robol <leonardo.robol@unipi.it>
 */

#include <mps/mps.h>

/**
 * @brief Create the monomial poly obtained by mp given by the coefficients
 * between degree \f$i\f$ and \f$j\f$. The resulting polynomial will have
 * degree \f$j - i\f$. 
 *
 * @param ctx The current mps_context. 
 * @param mp  The original monomial poly that should be sliced. 
 * @param i   The lower degree from which the coefficients should be copied. 
 * @param j   The higher degree from which the coefficients should be copied. 
 *
 * @return A newly allocated mps_monomial_poly.
 */
static mps_monomial_poly *
mps_slice_polynomial (mps_context * ctx, mps_monomial_poly * mp, 
		      int i, int j)
{
  /* Perform basic sanity checks on the indices. */
  if (j < i)
    return NULL;

  mps_monomial_poly * sliced_poly = mps_monomial_poly_new (ctx, j - i);
  int k;
  
  if (MPS_STRUCTURE_IS_RATIONAL (MPS_POLYNOMIAL (mp)->structure))
    {
      for (k = i; k <= j; k++)
	{
	  mps_monomial_poly_set_coefficient_q (ctx, sliced_poly, k - i, 
					       mp->initial_mqp_r[k], 
					       mp->initial_mqp_i[k]);
	}
    }
  else
    {
      for (k = i; k <= j; k++)
	{
	  mps_monomial_poly_set_coefficient_f (ctx, sliced_poly, k - i, mp->mfpc[k]);
	}      
    }

  return sliced_poly;
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
mps_recursive_fstart (mps_context * ctx, mps_polynomial * poly, mps_approximation ** approximations)
{
  MPS_DEBUG_THIS_CALL (ctx);

#ifndef DISABLE_DEBUG
  clock_t * recursive_start_timer = mps_start_timer();
#endif

  int i;

  if (! MPS_IS_MONOMIAL_POLY (poly))
    {
      MPS_DEBUG_WITH_INFO (ctx, "Falling back to starting_strategy: MPS_STARTING_STRATEGY_DEFAULT, since the"
			   "input is not a monomial poly");

      mps_context_select_starting_strategy (ctx, MPS_STARTING_STRATEGY_DEFAULT);
      mps_polynomial_fstart (ctx, poly, approximations);
      return;
    }

  if (poly->degree < 50)
    {
      (*poly->fstart)(ctx, poly, approximations);
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
  
  mps_monomial_poly * left = mps_slice_polynomial (rctx, mp, 0, middle);
  mps_monomial_poly * right = mps_slice_polynomial (rctx, mp, middle, poly->degree);

  /* Copy some configuration from the originating context */
  mps_context_add_debug_domain (rctx, ctx->debug_level);
  mps_context_select_algorithm (rctx, ctx->algorithm);
  // mps_context_select_starting_strategy (rctx, MPS_STARTING_STRATEGY_RECURSIVE);
  
  /* In every case we don't really need a high precision to find
   * good approximations. Here we're going for 16 bits of precision, 
   * but more reasoning could be put in this choice in the future. */
  mps_context_set_output_prec (rctx, 16);

  MPS_DEBUG_WITH_INFO (ctx, "Divided the polynomial into two polynomials of degree %d and %d", 
		       middle, poly->degree - middle);

  /* Recursively solve both polynomials. */
  mps_context_set_input_poly (rctx, MPS_POLYNOMIAL (left));
  mps_mpsolve (rctx);
  mps_approximation ** left_approximations = mps_context_get_approximations (rctx);

  mps_context_set_input_poly (rctx, MPS_POLYNOMIAL (right));
  mps_mpsolve (rctx);
  mps_approximation ** right_approximations = mps_context_get_approximations (rctx);

  /* Use these approximations as starting points for the original polynomial. */
  for (i = 0; i < poly->degree; i++)
    {
      mps_approximation * appr = (i < middle) ? left_approximations[i] : 
	right_approximations[i - middle];
      cplx_set (approximations[i]->fvalue, appr->fvalue);      
      free (appr);
    }

  free (left_approximations);
  free (right_approximations);

  mps_monomial_poly_free (rctx, MPS_POLYNOMIAL (left));
  mps_monomial_poly_free (rctx, MPS_POLYNOMIAL (right));
  mps_context_free (rctx);

  mps_starting_configuration_clear (ctx, &c);

#ifndef DISABLE_DEBUG  
  long int total_time = mps_stop_timer (recursive_start_timer);
  if (ctx->debug_level & MPS_DEBUG_TIMINGS)
    MPS_DEBUG (ctx, "Used %ld ms for the recursive starting strategy", 
	       total_time);
#endif
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
mps_recursive_dstart (mps_context * ctx, mps_polynomial * poly, mps_approximation ** approximations)
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
mps_recursive_mstart (mps_context * ctx, mps_polynomial * poly, mps_approximation ** approximations)
{
  
}
