/*
 * This file is part of MPSolve 3.1.8
 *
 * Copyright (C) 2001-2019, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Leonardo Robol <leonardo.robol@sns.it>
 */

#include <mps/mps.h>

mps_boolean
mps_chebyshev_poly_meval (mps_context * ctx, mps_polynomial * poly, mpcf_t x, mpcf_t value, rdpe_t error)
{
  long int wp = mpcf_get_prec (x);

  /* Lower the working precision in case of limited precision coefficients
   * in the input polynomial. */
  if (poly->prec > 0 && poly->prec < wp)
    wp = poly->prec;

  mps_chebyshev_poly * cpoly = MPS_CHEBYSHEV_POLY (poly);
  int i;

  mpcf_t t0, t1, ctmp, ctmp2;
  rdpe_t ax, rtmp, rtmp2;

  mpcf_rmod (ax, x);
  rdpe_set (error, rdpe_zero);

  /* Make sure that we have sufficient precision to perform the computation */
  mps_polynomial_raise_data (ctx, poly, wp);

  mpcf_init2 (t0, wp);
  mpcf_init2 (t1, wp);
  mpcf_init2 (ctmp, wp);
  mpcf_init2 (ctmp2, wp);

  mpcf_set (value, cpoly->mfpc[0]);
  mpcf_set_ui (t0, 1U, 0U);
  if (poly->degree == 0)
    {
      return true;
    }

  mpcf_set (t1, x);
  mpcf_mul (ctmp, cpoly->mfpc[1], x);
  mpcf_add_eq (value, ctmp);

  mpcf_rmod (rtmp, ctmp);
  rdpe_add_eq (error, rtmp);

  for (i = 2; i <= poly->degree; i++)
    {
      mpcf_mul (ctmp, x, t1);
      mpcf_mul_eq_ui (ctmp, 2U);
      mpcf_rmod (rtmp, ctmp);
      mpcf_sub_eq (ctmp, t0);

      mpcf_rmod (rtmp2, t0);
      rdpe_add_eq (rtmp, rtmp2);

      mpcf_mul (ctmp2, ctmp, cpoly->mfpc[i]);
      mpcf_add_eq (value, ctmp2);

      rdpe_mul_eq (rtmp, ax);
      rdpe_add_eq (error, rtmp);

      mpcf_set (t0, t1);
      mpcf_set (t1, ctmp);
    }

  mpcf_clear (t0);
  mpcf_clear (t1);
  mpcf_clear (ctmp);
  mpcf_clear (ctmp2);

  rdpe_set_2dl (rtmp, 2.0, -wp);
  rdpe_mul_eq (error, rtmp);

  return true;
}
