/*
 * This file is part of MPSolve 3.2.2
 *
 * Copyright (C) 2001-2020, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Dario Andrea Bini <bini@dm.unipi.it>
 *   Giuseppe Fiorentino <fiorent@dm.unipi.it>
 *   Leonardo Robol <leonardo.robol@unipi.it>
 */

#include <mps/mps.h>
#include <math.h>

#define MPS_STARTING_SIGMA (0.66 * (PI / ctx->n))
#define pi2 6.283184

void
mps_general_fstart (mps_context * ctx, mps_polynomial * p, mps_approximation ** approximations)
{
  int i;
  double sigma, ang;

  if (ctx->random_seed)
    sigma = drand ();
  else
    {
      sigma = ctx->last_sigma = MPS_STARTING_SIGMA;
    }

  /* In the case of user-defined polynomial choose as starting
   * approximations equally spaced points in the unit circle.  */
  ang = pi2 / ctx->n;
  for (i = 0; i < ctx->n; i++)
    {
      cplx_set_d (approximations[i]->fvalue, cos (ang * i + sigma),
                  sin (ang * i + sigma));
    }
}

void
mps_general_dstart (mps_context * ctx, mps_polynomial * p, mps_approximation ** approximations)
{
  int i;
  double sigma, ang;

  if (ctx->random_seed)
    sigma = drand ();
  else
    {
      sigma = ctx->last_sigma = MPS_STARTING_SIGMA;
    }

  /* In the case of user-defined polynomial choose as starting
   * approximations equally spaced points in the unit circle.  */
  ang = pi2 / ctx->n;
  for (i = 0; i < ctx->n; i++)
    {
      cdpe_set_d (approximations[i]->dvalue, cos (ang * i + sigma),
                  sin (ang * i + sigma));
    }
}

void
mps_general_mstart (mps_context * ctx, mps_polynomial * p, mps_approximation ** approximations)
{
  int i;
  double sigma, ang;

  if (ctx->random_seed)
    sigma = drand ();
  else
    {
      sigma = ctx->last_sigma = MPS_STARTING_SIGMA;
    }

  /* In the case of user-defined polynomial choose as starting
   * approximations equally spaced points in the unit circle.  */
  ang = pi2 / ctx->n;
  for (i = 0; i < ctx->n; i++)
    {
      cplx_t tmp;
      cplx_set_d (tmp, cos (ang * i + sigma),
                  sin (ang * i + sigma));
      mpc_set_cplx (approximations[i]->mvalue, tmp);
    }
}
