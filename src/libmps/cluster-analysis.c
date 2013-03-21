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

void
mps_cluster_analysis (mps_context * ctx, mps_polynomial * p)
{
  switch (ctx->lastphase)
  {
    case float_phase:
      {
        double * radii = double_valloc (ctx->n);

        mps_fradii (ctx, p, radii);
        mps_fcluster (ctx, radii, 2 * ctx->n);
        mps_fmodify (ctx, false);

        free (radii);
        break;
      }

    case dpe_phase:
      {
        rdpe_t * radii = rdpe_valloc (ctx->n);

        mps_dradii (ctx, p, radii);
        mps_dcluster (ctx, radii, 2 * ctx->n);
        mps_dmodify (ctx, false);

        free (radii);
        break;
      }

    case mp_phase:
      {
        rdpe_t * radii = rdpe_valloc (ctx->n);

        mps_mradii (ctx, p, radii);
        mps_mcluster (ctx, radii, 2 * ctx->n);
        mps_mmodify (ctx, false);

        free (radii);
        break;
      }

    case no_phase:
      break;
  }
}

