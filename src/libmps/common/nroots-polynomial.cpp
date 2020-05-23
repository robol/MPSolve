/*
 * This file is part of MPSolve 3.1.9
 *
 * Copyright (C) 2001-2014, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Leonardo Robol <robol@mail.dm.unipi.it>
 */

#include <mps/mps.h>

using namespace mps;

NRootsPolynomial::NRootsPolynomial (mps_context * ctx, int n) : 
  Polynomial (ctx, "NRootsPolynomial")
{
  /* Set the degree inside the mps_polynomial struct */
  degree = n;
}

mps_boolean
NRootsPolynomial::eval (mps_context * ctx, cplx_t x, cplx_t value, double * error)
{
  cplx_pow_si (value, x, degree);
  cplx_sub_eq (value, cplx_one);
  return true;
}

mps_boolean
NRootsPolynomial::eval (mps_context * ctx, cdpe_t x, cdpe_t value, rdpe_t error)
{
  cdpe_pow_si (value, x, degree);
  cdpe_sub_eq (value, cdpe_one);
  return true;
}

mps_boolean
NRootsPolynomial::eval (mps_context * ctx, mpc_t x, mpc_t value, rdpe_t error)
{
  mpc_pow_si (value, x, degree);
  mpc_sub_eq_ui (value, 1U, 0U);

  mpc_rmod (error, value);
  rdpe_add_eq (error, rdpe_one);

  rdpe_mul_eq (error, ctx->mp_epsilon);

  return true;
}

void 
NRootsPolynomial::newton (mps_context * ctx, mps_approximation * a, cplx_t x)
{
}

void
NRootsPolynomial::newton (mps_context * ctx, mps_approximation * a, cdpe_t x)
{
  rdpe_t ax; 
  cdpe_mod (ax, a->dvalue);

  cdpe_pow_si (x, a->dvalue, 1 - degree);
  cdpe_neg_eq (x);
  cdpe_add_eq (x, a->dvalue);

  cdpe_div_eq_d (x, degree);

  cdpe_mod (a->drad, x);
  rdpe_mul_eq_d (a->drad, degree);

  rdpe_mul_eq_d (ax, DBL_EPSILON * 4.0);
  a->again = rdpe_gt (a->drad, ax);
}

void
NRootsPolynomial::newton (mps_context * ctx, mps_approximation * a, mpc_t x, long int wp)
{
}

