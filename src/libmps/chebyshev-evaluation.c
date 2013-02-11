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

mps_boolean
mps_chebyshev_poly_meval (mps_context * ctx, mps_polynomial * poly, mpc_t x, mpc_t value, rdpe_t error)
{
  long int wp = mpc_get_prec (x);
  mps_chebyshev_poly * cpoly = MPS_CHEBYSHEV_POLY (poly);
  int i;

  mpc_t t0, t1, ctmp, ctmp2;
  rdpe_t ax, rtmp, rtmp2;

  mpc_rmod (ax, x);
  rdpe_set (error, rdpe_zero);

  /* Make sure that we have sufficient precision to perform the computation */
  mps_polynomial_raise_data (ctx, poly, wp);

  mpc_init2 (t0, wp);
  mpc_init2 (t1, wp);
  mpc_init2 (ctmp, wp);
  mpc_init2 (ctmp2, wp);

  mpc_set (value, cpoly->mfpc[0]);
  mpc_set_ui (t0, 1U, 0U);
  if (poly->degree == 0) 
    {
      return true;
    }

  mpc_set (t1, x);
  mpc_mul (ctmp, cpoly->mfpc[1], x);
  mpc_add_eq (value, ctmp);

  mpc_rmod (rtmp, ctmp);
  rdpe_add_eq (error, rtmp);

  for (i = 2; i <= poly->degree; i++) 
    {
      mpc_mul (ctmp, x, t1);
      mpc_mul_eq_ui (ctmp, 2U);
      mpc_rmod (rtmp, ctmp);
      mpc_sub_eq (ctmp, t0);

      mpc_rmod (rtmp2, t0);
      rdpe_add_eq (rtmp, rtmp2);

      mpc_mul (ctmp2, ctmp, cpoly->mfpc[i]);
      mpc_add_eq (value, ctmp2);

      rdpe_mul_eq (rtmp, ax);
      rdpe_add_eq (error, rtmp);

      mpc_set (t0, t1);
      mpc_set (t1, ctmp);
    }

  mpc_clear (t0);
  mpc_clear (t1);
  mpc_clear (ctmp);
  mpc_clear (ctmp2);

  rdpe_set_2dl (rtmp, 2.0, -wp);
  rdpe_mul_eq (error, rtmp);

  return true;
}
