/*
 * This file is part of MPSolve 3.0
 *
 * Copyright (C) 2001-2013, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors: 
 *   Leonardo Robol <robol@mail.dm.unipi.it>
 */

#include <mps/mps.h>

/**
 * @brief Perform a step of Aberth method in Jacobi-style.
 *
 * @param ctx The context in which this instance of MPSolve is running.
 * @param p The polynomial on which Aberth method should be applied.
 * @param approximations The starting approximations for the method. They will be
 *        updated with the new values when the function returns. 
 */
mps_boolean
mps_fjacobi_aberth_step (mps_context * ctx, mps_polynomial * p)
{
  cplx_t * faberth_corrections = NULL;
  mps_boolean again = false;
  int i = 0;

  faberth_corrections = cplx_valloc (ctx->n);

  for (i = 0; i < ctx->n; i++)
    {
      cplx_t abcorr;
      mps_polynomial_fnewton (ctx, p, ctx->root[i], faberth_corrections[i]);
      mps_faberth (ctx, i, abcorr);

      cplx_mul_eq (abcorr, faberth_corrections[i]);
      cplx_sub (abcorr, cplx_one, abcorr);
      cplx_div_eq (faberth_corrections[i], abcorr);
    }

  /* Update approximations with the new corrections */
  for (i = 0; i < ctx->n; i++)
    {
      cplx_sub_eq (ctx->root[i]->fvalue, faberth_corrections[i]);
      again |= ctx->root[i]->again;
    }

  cplx_vfree (faberth_corrections);
  return again;
}

/**
 * @brief Perform a packet of Aberth iterations on the approximation
 * to the roots of p. 
 *
 * @params ctx Current MPSolve context.
 * @params p The polynomial whose roots should be approximated.
 */
void
mps_faberth_packet (mps_context * ctx, mps_polynomial * p)
{
  double * radii = double_valloc (ctx->n);

  do 
    {
      if (ctx->debug_level & MPS_DEBUG_APPROXIMATIONS)
        MPS_DEBUG (ctx, "Performing a step of Aberth");
    } while (mps_fjacobi_aberth_step (ctx, p));

  mps_fradii (ctx, radii);
  mps_fcluster (ctx, radii, 2 * ctx->n);

  mps_fmodify (ctx, false);

  double_vfree (radii);
}