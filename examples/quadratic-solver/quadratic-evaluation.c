/*
 * This file is part of MPSolve 3.0
 *
 * Copyright (C) 2001-2013, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Leonardo Robol <robol@mail.dm.unipi.it>
 */

#include "quadratic-poly.h"

#include <float.h>
#define LOG2     0.69314718055994530941

mps_boolean
mps_quadratic_poly_feval (mps_context * ctx, mps_polynomial * p, cplx_t x, cplx_t value, double * error)
{
  int i;
  int m = (int) (log (p->degree + 1.0) / LOG2);
  double ax = cplx_mod (x);
  cplx_t tmp;

  if ((1 << m) <= p->degree)
    m++;

  cplx_set (value, cplx_one);

  if (error)
    *error = 0.0;

  for (i = 1; i <= m; i++)
    {
      cplx_sqr (tmp, value);
      cplx_mul (value, x, tmp);
      cplx_add_eq (value, cplx_one);

      if (error)
        *error = *error * ax + cplx_mod (value);
    }

  if (error)
    *error *= DBL_EPSILON;

  return true;
}

mps_boolean
mps_quadratic_poly_deval (mps_context * ctx, mps_polynomial * p, cdpe_t x, cdpe_t value, rdpe_t error)
{
  int i;
  int m = (int) (log (p->degree + 1.0) / LOG2);
  rdpe_t ax, rtmp;
  cdpe_t tmp;

  if ((1 << m) <= p->degree)
    m++;

  cdpe_mod (ax, x);
  cdpe_set (value, cdpe_one);
  cdpe_mod (error, value);

  for (i = 1; i <= m; i++)
    {
      cdpe_sqr (tmp, value);
      cdpe_mul (value, x, tmp);
      cdpe_add_eq (value, cdpe_one);

      rdpe_mul_eq (error, ax);
      cdpe_mod (rtmp, value);
      rdpe_add_eq (error, rtmp);
    }

  rdpe_mul_eq_d (error, DBL_EPSILON);
  return true;
}

mps_boolean
mps_quadratic_poly_meval (mps_context * ctx, mps_polynomial * p, mpcf_t x, mpcf_t value, rdpe_t error)
{
  int i;
  int m = (int) (log (p->degree + 1.0) / LOG2);
  rdpe_t ax, rtmp;
  mpcf_t tmp;
  long int wp = mpcf_get_prec (x);

  /* Correct the working precision in case of a limited 
     precision polynomial (quite unlikely
   * in the quadratic case, but still. */
  if (p->prec > 0 && p->prec < wp)
    wp = p->prec;

  if ((1 << m) <= p->degree)
    m++;

  mpcf_init2 (tmp, wp);

  mpcf_rmod (ax, x);
  mpcf_set_ui (value, 1U, 0U);
  mpcf_rmod (error, value);

  for (i = 1; i <= m; i++)
    {
      mpcf_sqr (tmp, value);
      mpcf_mul (value, x, tmp);
      mpcf_add_eq_ui (value, 1U, 0U);

      rdpe_mul_eq (error, ax);
      mpcf_rmod (rtmp, value);
      rdpe_add_eq (error, rtmp);
    }

  rdpe_set_2dl (rtmp, 1.0, -wp);
  rdpe_mul_eq (error, rtmp);

  mpcf_clear (tmp);

  return true;
}
