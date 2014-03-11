/*
 * This file is part of MPSolve 3.1.5
 *
 * Copyright (C) 2001-2014, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Dario Andrea Bini <bini@dm.unipi.it>
 *   Giuseppe Fiorentino <fiorent@dm.unipi.it>
 *   Leonardo Robol <robol@mail.dm.unipi.it>
 */


#include <mps/mps.h>
#include <math.h>

#define pi2 6.283184

/**
 * @brief Compute some sensible starting points for the given secular equation (floating point version).
 *
 * @param s The current mps_context.
 * @param sec The secular equation for which the starting points shall be computed.
 */
void
mps_secular_fstart (mps_context * s, mps_secular_equation * sec)
{
  MPS_DEBUG_THIS_CALL (s);

  int i;
  int n = MPS_POLYNOMIAL (sec)->degree;

  for (i = 0; i < n; i++)
    {
      if (!MPS_ROOT_STATUS_IS_COMPUTED (s->root[i]->status))
        {
          cplx_set_d (s->root[i]->fvalue, cos (i * n) * DBL_EPSILON * 4.0 * cplx_mod (sec->bfpc[i]),
                      sin (i * n) * DBL_EPSILON * 4.0 * cplx_mod (sec->bfpc[i]));

          s->root[i]->frad += cplx_mod (s->root[i]->fvalue);
          cplx_add_eq (s->root[i]->fvalue, sec->bfpc[i]);

          if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)
            MPS_DEBUG_CPLX (s, s->root[i]->fvalue, "s->froot[%d]", i);
        }
    }
}

/**
 * @brief Compute some sensible starting points for the given secular equation (DPE version).
 *
 * @param s The current mps_context.
 * @param sec The secular equation for which the starting points shall be computed.
 */
void
mps_secular_dstart (mps_context * s, mps_secular_equation * sec)
{
  MPS_DEBUG_THIS_CALL (s);

  int l;
  int n = MPS_POLYNOMIAL (sec)->degree;

  for (l = 0; l < MPS_POLYNOMIAL (sec)->degree; l++)
    {
      if (!MPS_ROOT_STATUS_IS_COMPUTED (s->root[l]->status))
        {
          rdpe_t bmod;

          cdpe_mod (bmod, sec->bdpc[l]);
          rdpe_mul_eq_d (bmod, 4.0 * DBL_EPSILON);

          rdpe_mul_d (cdpe_Re (s->root[l]->dvalue), bmod, cos (l * n));
          rdpe_mul_d (cdpe_Im (s->root[l]->dvalue), bmod, sin (l * n));

          rdpe_add_eq (s->root[l]->drad, bmod);
          cdpe_add_eq (s->root[l]->dvalue, sec->bdpc[l]);

          if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)
            {
              MPS_DEBUG_CDPE (s, s->root[l]->dvalue, "s->droot[%d]", l);
            }
        }
    }
}

/**
 * @brief Compute some sensible starting points for the given secular equation (MP version).
 *
 * @param s The current mps_context.
 * @param sec The secular equation for which the starting points shall be computed.
 */
void
mps_secular_mstart (mps_context * s, mps_secular_equation * sec)
{
  MPS_DEBUG_THIS_CALL (s);

  int l;
  int n = MPS_POLYNOMIAL (sec)->degree;

  for (l = 0; l < MPS_POLYNOMIAL (sec)->degree; l++)
    {
      if (!MPS_ROOT_STATUS_IS_COMPUTED (s->root[l]->status))
        {
          rdpe_t bmod;
          cdpe_t ctmp;

          mpc_rmod (bmod, sec->bmpc[l]);
          rdpe_mul_eq (bmod, s->mp_epsilon);
          rdpe_mul_eq_d (bmod, 4.0);

          rdpe_mul_d (cdpe_Re (ctmp), bmod, cos (l * n));
          rdpe_mul_d (cdpe_Im (ctmp), bmod, sin (l * n));

          mpc_set_cdpe (s->root[l]->mvalue, ctmp);
          rdpe_add_eq (s->root[l]->drad, bmod);

          mpc_add_eq (s->root[l]->mvalue, sec->bmpc[l]);
        }
    }
}
