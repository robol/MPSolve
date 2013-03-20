/*
 * This file is part of MPSolve 3.0
 *
 * Copyright (C) 2001-2013, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors: 
 *   Dario Andrea Bini <bini@dm.unipi.it>
 *   Giuseppe Fiorentino <fiorent@dm.unipi.it>
 *   Leonardo Robol <robol@mail.dm.unipi.it>
 */


#include <mps/mps.h>

mps_approximation *
mps_approximation_new (mps_context * s)
{
  mps_approximation * appr = mps_new (mps_approximation);
  mpc_init2 (appr->mvalue, s->mpwp);
  appr->again = true;
  appr->approximated = false;

  appr->status = MPS_ROOT_STATUS_CLUSTERED;
  appr->attrs  = MPS_ROOT_ATTRS_NONE;
  appr->inclusion = MPS_ROOT_INCLUSION_UNKNOWN;

  return appr;
}

void
mps_approximation_free (mps_context * s, mps_approximation * appr)
{
  mpc_clear (appr->mvalue);
  free (appr);
}

mps_approximation *
mps_approximation_copy (mps_context * ctx, mps_approximation * original)
{
  mps_approximation *new = mps_approximation_new (ctx);
  mpc_set_prec (new->mvalue, mpc_get_prec (original->mvalue));
  mpc_set (new->mvalue, original->mvalue);
  rdpe_set (new->drad, original->drad);
  cdpe_set (new->dvalue, original->dvalue);
  cplx_set (new->fvalue, original->fvalue);
  new->frad = original->frad;
  new->wp = original->wp;
  new->status = original->status;
  new->attrs = original->attrs;
  new->inclusion = original->inclusion;
  return new;
}
