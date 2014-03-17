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

mps_approximation *
mps_approximation_new (mps_context * s)
{
  mps_approximation * appr = mps_new (mps_approximation);

  mpc_init2 (appr->mvalue, s->mpwp);
  appr->again = true;
  appr->approximated = false;

  appr->status = MPS_ROOT_STATUS_CLUSTERED;
  appr->attrs = MPS_ROOT_ATTRS_NONE;
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

/* Public accessor functions */
void
mps_approximation_get_fvalue (mps_context * ctx, mps_approximation * approximation, cplx_t output)
{
  cplx_set (output, approximation->fvalue);
}

void
mps_approximation_get_dvalue (mps_context * ctx, mps_approximation * approximation, cdpe_t output)
{
  cdpe_set (output, approximation->dvalue);
}

void
mps_approximation_get_mvalue (mps_context * ctx, mps_approximation * approximation, mpc_t output)
{
  mpc_set_prec (output, mpc_get_prec (approximation->mvalue));
  mpc_set (output, approximation->mvalue);
}

double
mps_approximation_get_frad (mps_context * ctx, mps_approximation * approximation)
{
  return approximation->frad;
}

void
mps_approximation_get_drad (mps_context * ctx, mps_approximation * approximation, rdpe_t output)
{
  rdpe_set (output, approximation->drad);
}

mps_root_status
mps_approximation_get_status (mps_context * ctx, mps_approximation * approximation)
{
  return approximation->status;
}

mps_root_attrs
mps_approximation_get_attrs (mps_context * ctx, mps_approximation * approximation)
{
  return approximation->attrs;
}

mps_root_inclusion
mps_approximaiton_get_inclusion (mps_context * ctx, mps_approximation * approximation)
{
  return approximation->inclusion;
}

mps_boolean
mps_approximation_get_again (mps_context * ctx, mps_approximation * approximation)
{
  return approximation->again;
}

/* Public setters functions */
void
mps_approximation_set_fvalue (mps_context * ctx, mps_approximation * approximation, const cplx_t value)
{
  cplx_set (approximation->fvalue, value);
}

void
mps_approximation_set_dvalue (mps_context * ctx, mps_approximation * approximation, const cdpe_t value)
{
  cdpe_set (approximation->dvalue, value);
}

void
mps_approximation_set_mvalue (mps_context * ctx, mps_approximation * approximation, const mpc_t value)
{
  /* Ensure that we have a sufficient precision to store value correctly */
  if (mpc_get_prec (value) > approximation->wp)
    {
      mpc_set_prec (approximation->mvalue, mpc_get_prec (value));
      approximation->wp = mpc_get_prec (approximation->mvalue);
    }

  mpc_set (approximation->mvalue, value);
}

void
mps_approximation_set_frad (mps_context * ctx, mps_approximation * approximation, const double frad)
{
  approximation->frad = frad;
}

void
mps_approximation_set_drad (mps_context * ctx, mps_approximation * approximation, const rdpe_t drad)
{
  rdpe_set (approximation->drad, drad);
}

void
mps_approximation_set_status (mps_context * ctx, mps_approximation * approximation, const mps_root_status status)
{
  approximation->status = status;
}

void
mps_approximation_set_attrs (mps_context * ctx, mps_approximation * approximation, const mps_root_attrs attrs)
{
  approximation->attrs = attrs;
}

void
mps_approximation_set_inclusion (mps_context * ctx, mps_approximation * approximation, const mps_root_inclusion inclusion)
{
  approximation->inclusion = inclusion;
}

void
mps_approximation_set_again (mps_context * ctx, mps_approximation * approximation, const mps_boolean again)
{
  approximation->again = again;
}
