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
#include <string.h>

void
_mps_polynomial_free (mps_context * ctx, mps_polynomial * p)
{
  free (p);
}

long int
_mps_polynomial_raise_data (mps_context * ctx, mps_polynomial * p, long int wp)
{
  return wp;
}

void
_mps_polynomial_get_leading_coefficient (mps_context * ctx, mps_polynomial * p, mpc_t lc)
{
  mpc_set_ui (lc, 1U, 0U);
}

void 
mps_polynomial_init (mps_context * ctx, mps_polynomial * p)
{
  p->type_name = NULL;
  p->thread_safe = true;
  p->prec = 0;
  p->feval = NULL;
  p->deval = NULL;
  p->meval = NULL;
  p->fstart = mps_general_fstart;
  p->dstart = mps_general_dstart;
  p->mstart = mps_general_mstart;
  p->free = _mps_polynomial_free;
  p->raise_data = _mps_polynomial_raise_data;
  p->fnewton = NULL;
  p->dnewton = NULL;
  p->mnewton = NULL;
  p->get_leading_coefficient = _mps_polynomial_get_leading_coefficient;
}

mps_polynomial * 
mps_polynomial_new (mps_context * ctx)
{
  mps_polynomial *p = mps_new (mps_polynomial);
  mps_polynomial_init (ctx, p);
  return p;
}

mps_boolean
mps_polynomial_check_type (mps_polynomial * p, const char * type_name)
{
  return (p->type_name && (strcmp (p->type_name, type_name) == 0));
}

mps_boolean
mps_polynomial_feval (mps_context * ctx, mps_polynomial * p, cplx_t x, cplx_t value, double * error)
{
  return (*p->feval)(ctx, p, x, value, error);
}

mps_boolean
mps_polynomial_deval (mps_context * ctx, mps_polynomial * p, cdpe_t x, cdpe_t value, rdpe_t error)
{
  return (*p->deval)(ctx, p, x, value, error);
}

mps_boolean 
mps_polynomial_meval (mps_context * ctx, mps_polynomial * p, mpc_t x, mpc_t value, rdpe_t error)
{
  return (*p->meval)(ctx, p, x, value, error);
}

void
mps_polynomial_fstart (mps_context * ctx, mps_polynomial * p)
{
  ctx->operation = MPS_OPERATION_STARTING_POINTS_FP;
  (*p->fstart) (ctx, p);
}

void 
mps_polynomial_dstart (mps_context * ctx, mps_polynomial * p)
{
  ctx->operation = MPS_OPERATION_STARTING_POINTS_DPE;
  (*p->dstart) (ctx, p);
}

void 
mps_polynomial_mstart (mps_context * ctx, mps_polynomial * p)
{  
  ctx->operation = MPS_OPERATION_STARTING_POINTS_MP;
  (*p->mstart) (ctx, p);
}

void 
mps_polynomial_free (mps_context * ctx, mps_polynomial * p)
{
  (*p->free)(ctx, p);
}

long int
mps_polynomial_raise_data (mps_context * ctx, mps_polynomial * p, long int wp)
{
  return (*p->raise_data)(ctx, p, wp);
}

void 
mps_polynomial_fnewton (mps_context * ctx, mps_polynomial *p, 
                        mps_approximation * root, cplx_t corr)
{
  (*p->fnewton)(ctx, p, root, corr);
}

void 
mps_polynomial_dnewton (mps_context * ctx, mps_polynomial *p, 
                        mps_approximation * root, cdpe_t corr)
{
  (*p->dnewton) (ctx, p, root, corr);
}

void 
mps_polynomial_mnewton (mps_context * ctx, mps_polynomial *p, 
                        mps_approximation * root, mpc_t corr)
{
  (*p->mnewton) (ctx, p, root, corr);
}

void
mps_polynomial_get_leading_coefficient (mps_context * ctx, mps_polynomial * p,
                                        mpc_t leading_coefficient)
{
  (*p->get_leading_coefficient) (ctx, p, leading_coefficient);
}
