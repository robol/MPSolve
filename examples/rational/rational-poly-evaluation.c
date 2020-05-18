#include <mps/mps.h>
#include "rational-poly.h"
#include <float.h>

mps_boolean mps_rational_poly_feval (mps_context * ctx, mps_polynomial * p, 
				     cplx_t x, cplx_t value, double * error)
{
  mps_rational_poly * rp = MPS_RATIONAL_POLY (p);
  int i;
  cplx_t ftmp1;

  cplx_set (value, cplx_zero);
  *error = 0;

  /* Sum all the terms of the form c_i / (a_i + x * b_i)^p */
  for (i = 0; i < rp->nn; i++)
    {
      cplx_mul (ftmp1, rp->b[i], x);
      cplx_add_eq (ftmp1, rp->a[i]);
      cplx_pow_si (ftmp1, ftmp1, rp->pp);
      cplx_div (ftmp1, rp->c[i], ftmp1);

      *error += cplx_mod (ftmp1);

      cplx_add_eq (value, ftmp1);
    }

  *error *= DBL_EPSILON * 6;
}

mps_boolean mps_rational_poly_deval (mps_context * ctx, mps_polynomial * p, 
				     cdpe_t x, cdpe_t value, rdpe_t error)
{
  mps_rational_poly * rp = MPS_RATIONAL_POLY (p);
  int i;
  cdpe_t dtmp1;
  rdpe_t rtmp1;

  cdpe_set (value, cdpe_zero);
  rdpe_set (error, rdpe_zero);

  /* Sum all the terms of the form c_i / (a_i + x * b_i)^p */
  for (i = 0; i < rp->nn; i++)
    {
      cdpe_t a, b, c;

      cdpe_set_x (a, rp->a[i]);
      cdpe_set_x (b, rp->b[i]);
      cdpe_set_x (c, rp->c[i]);

      cdpe_mul (dtmp1, b, x);
      cdpe_add_eq (dtmp1, a);
      cdpe_pow_si (dtmp1, dtmp1, rp->pp);
      cdpe_div (dtmp1, c, dtmp1);

      cdpe_mod (rtmp1, dtmp1);
      rdpe_add_eq (error, rtmp1);

      cdpe_add_eq (value, dtmp1);
    }

  rdpe_mul_eq_d (error, DBL_EPSILON * 6);
}

mps_boolean mps_rational_poly_meval (mps_context * ctx, mps_polynomial * p, 
				     mpcf_t x, mpcf_t value, rdpe_t error)
{
  mps_rational_poly * rp = MPS_RATIONAL_POLY (p);
  int i;
  mpcf_t mtmp1, prod, a, b, c;
  rdpe_t rtmp1;

  long int wp = mpcf_get_prec (value);

  mpcf_set_ui (value, 0U, 0U);
  rdpe_set (error, rdpe_zero);

  mpcf_init2 (a, wp);
  mpcf_init2 (b, wp);
  mpcf_init2 (c, wp);
  mpcf_init2 (mtmp1, wp);
  mpcf_init2 (prod, wp);

  mpcf_set_ui (prod, 1U, 0U);

  /* Sum all the terms of the form c_i / (a_i + x * b_i)^p */
  for (i = 0; i < rp->nn; i++)
    {
      mpcf_set_d (a, cplx_Re (rp->a[i]), cplx_Im (rp->a[i]));
      mpcf_set_d (b, cplx_Re (rp->b[i]), cplx_Im (rp->b[i]));
      mpcf_set_d (c, cplx_Re (rp->c[i]), cplx_Im (rp->c[i]));

      mpcf_mul (mtmp1, b, x);
      mpcf_add_eq (mtmp1, a);
      mpcf_pow_si (mtmp1, mtmp1, rp->pp);
      mpcf_mul_eq (prod, mtmp1);
      mpcf_div (mtmp1, c, mtmp1);

      mpcf_rmod (rtmp1, mtmp1);
      rdpe_add_eq (error, rtmp1);

      mpcf_add_eq (value, mtmp1);
    }

  mpcf_mul_eq (value, prod);

  mpcf_clear (a);
  mpcf_clear (b);
  mpcf_clear (c);
  mpcf_clear (mtmp1);

  rdpe_mul_eq_d (error, DBL_EPSILON * 6);  
}
