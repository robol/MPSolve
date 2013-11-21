/*
 * This file is part of MPSolve 3.0
 *
 * Copyright (C) 2001-2013, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors: 
 *   Leonardo Robol <leonardo.robol@sns.it>
 */

#include <mps/mps.h>

static inline void 
switch_to_mp (mps_context * ctx)
{
  switch (ctx->algorithm)
    {

    case MPS_ALGORITHM_SECULAR_GA:
      mps_secular_switch_phase (ctx, mp_phase);
      break;

    case MPS_ALGORITHM_STANDARD_MPSOLVE:
      ctx->lastphase = mp_phase;
      mps_mp_set_prec (ctx, 2 * DBL_MANT_DIG);
      mps_prepare_data (ctx, ctx->mpwp);
      break;
    }
}

/**
 * @brief This function can be called to validate the inclusion radii and
 * cluster analysis for a limited precision polynomial. 
 *
 * In the current implementation MPSolve treat a limited precision polynomial
 * as an infinite one, and then gives a poteriori bounds to the approximations
 * using this function.
 *
 * @param ctx The current mps_context
 */
void 
mps_validate_inclusions (mps_context * ctx)
{
  MPS_DEBUG_THIS_CALL (ctx);
  
  mps_polynomial * current_poly = ctx->active_poly;
  int i = 0;
  mpc_t value;

  if (ctx->lastphase != mp_phase)
    switch_to_mp (ctx);

  /* Quick hackish way of making sure that we are using the higher precision available on
   * this sytem that is lower of the input precision. */
  long int prec = current_poly->prec;

  mpc_init2 (value, prec);
  mps_polynomial_raise_data (ctx, current_poly, prec);

  for (i = 0; i < ctx->n; i++)
    {
      /* Discard current information on the inclusion radius */
      ctx->root[i]->frad = DBL_MAX;
      rdpe_set (ctx->root[i]->drad, RDPE_MAX);

      mpc_set_prec (ctx->root[i]->mvalue, prec);
      mps_polynomial_mnewton (ctx, current_poly, ctx->root[i], 
			      value, prec);
    }

  mpc_clear (value);

  /* Perform a step of cluster analysis */
  mps_cluster_analysis (ctx, current_poly);
}
