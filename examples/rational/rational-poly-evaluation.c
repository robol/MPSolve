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
				     mpc_t x, mpc_t value, rdpe_t error)
{
  mps_rational_poly * rp = MPS_RATIONAL_POLY (p);
  int i;
  mpc_t mtmp1, prod, a, b, c;
  rdpe_t rtmp1;

  long int wp = mpc_get_prec (value);

  mpc_set_ui (value, 0U, 0U);
  rdpe_set (error, rdpe_zero);

  mpc_init2 (a, wp);
  mpc_init2 (b, wp);
  mpc_init2 (c, wp);
  mpc_init2 (mtmp1, wp);
  mpc_init2 (prod, wp);

  mpc_set_ui (prod, 1U, 0U);

  /* Sum all the terms of the form c_i / (a_i + x * b_i)^p */
  for (i = 0; i < rp->nn; i++)
    {
      mpc_set_d (a, cplx_Re (rp->a[i]), cplx_Im (rp->a[i]));
      mpc_set_d (b, cplx_Re (rp->b[i]), cplx_Im (rp->b[i]));
      mpc_set_d (c, cplx_Re (rp->c[i]), cplx_Im (rp->c[i]));

      mpc_mul (mtmp1, b, x);
      mpc_add_eq (mtmp1, a);
      mpc_pow_si (mtmp1, mtmp1, rp->pp);
      mpc_mul_eq (prod, mtmp1);
      mpc_div (mtmp1, c, mtmp1);

      mpc_rmod (rtmp1, mtmp1);
      rdpe_add_eq (error, rtmp1);

      mpc_add_eq (value, mtmp1);
    }

  mpc_mul_eq (value, prod);

  mpc_clear (a);
  mpc_clear (b);
  mpc_clear (c);
  mpc_clear (mtmp1);

  rdpe_mul_eq_d (error, DBL_EPSILON * 6);  
}
