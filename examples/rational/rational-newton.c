#include <mps/mps.h>
#include "rational-poly.h"

void mps_rational_poly_fnewton (mps_context * ctx, mps_polynomial * p,
				mps_approximation * root, cplx_t corr)
{
  /* TODO: Implement this function */
}

void mps_rational_poly_dnewton (mps_context * ctx, mps_polynomial * p,
				mps_approximation * root, cdpe_t corr)
{
  /* TODO: Implement this function */
}

void mps_rational_poly_mnewton (mps_context * ctx, mps_polynomial * p,
				mps_approximation * root, mpc_t corr, 
				long int wp)
{
  int i, nn, pp;
  mps_rational_poly * rp = MPS_RATIONAL_POLY (p);

  nn = rp->nn;
  pp = rp->pp;

  /* These variables are used to store the sum of the absolute
   * values of the different pieces of the Newton correction. */
  rdpe_t asum, asumd1, asumd2, rtmp1;

  mpc_t a_plus_bx, a_plus_bxp, a, b, c;
  mpc_t f, fd1, fd2, x;

  mpc_init2 (x, wp);
  mpc_init2 (a_plus_bx, wp);
  mpc_init2 (a_plus_bxp, wp);
  mpc_init2 (a, wp);
  mpc_init2 (b, wp);
  mpc_init2 (c, wp);
  mpc_init2 (f, wp);
  mpc_init2 (fd1, wp);
  mpc_init2 (fd2, wp);

  rdpe_set (asum, rdpe_zero);
  rdpe_set (asumd1, rdpe_zero);
  rdpe_set (asumd2, rdpe_zero);

  mpc_set_ui (corr, 0U, 0U);
  mpc_set_ui (f, 0U, 0U);
  mpc_set_ui (fd1, 0U, 0U);
  mpc_set_ui (fd2, 0U, 0U);

  mps_approximation_get_mvalue (ctx, root, x);

  for (i = 0; i < nn; i++)
    {
      /* Grab a copy of the coefficients of the rational function. */
      /* Note that, in principle, if integers or rationals coefficients
       * are available we should use them here, but they are not
       * implemented yet. */
      mpc_set_d (a, cplx_Re (rp->a[i]), cplx_Im (rp->a[i]));
      mpc_set_d (b, cplx_Re (rp->b[i]), cplx_Im (rp->b[i]));
      mpc_set_d (c, cplx_Re (rp->c[i]), cplx_Im (rp->c[i]));

      /* Compute 1 / (a + bx) and its p-th power */
      mpc_mul (a_plus_bx, b, x);
      mpc_add_eq (a_plus_bx, a);
      mpc_inv_eq (a_plus_bx);
      mpc_pow_si (a_plus_bxp, a_plus_bx, pp);

      /* Compute c / (a + bx)^p, b / (a + bx) and 
       * cb / (a + bx)^{p+1}. */
      mpc_mul (a, c, a_plus_bxp);
      mpc_mul (b, b, a_plus_bx);
      mpc_mul (c, a, b);

      mpc_mul_eq_ui (b, pp);
      mpc_mul_eq_ui (c, pp);

      /* Store the moduli of these elements. */
      mpc_rmod (rtmp1, a);
      rdpe_add_eq (asum, rtmp1);
      mpc_rmod (rtmp1, b);
      rdpe_add_eq (asumd1, rtmp1);
      mpc_rmod (rtmp1, c);
      rdpe_add_eq (asumd2, rtmp1);

      /* Eventually compute the sums that we need. */
      mpc_add_eq (f, a);
      mpc_add_eq (fd1, b);
      mpc_add_eq (fd2, c);      
    }

  rdpe_mul_eq (asumd1, asum);
  rdpe_sub (asumd1, asumd1, asumd2);
  rdpe_abs (asumd1, asumd1);
  rdpe_div_eq (asum, asumd1);

  mpc_set (corr, f);
  mpc_mul_eq (fd1, f);
  mpc_sub_eq (fd1, fd2);
  mpc_div_eq (corr, fd1);

  /* Compute the radius for the new approximation */
  rdpe_set_2dl (rtmp1, 1.0, -wp);
  rdpe_mul_eq (asum, rtmp1);
  rdpe_mul_eq_d (asum, 4.0);
  mpc_rmod (rtmp1, corr);
  rdpe_add_eq (asum, rtmp1);
  rdpe_mul_eq_d (asum, nn);

  /* TODO: We should check the root-neighborhood condition here. */
  mps_approximation_set_drad (ctx, root, asum);

  mpc_clear (a_plus_bx);
  mpc_clear (a_plus_bxp);
  mpc_clear (a);
  mpc_clear (b);
  mpc_clear (c);
  mpc_clear (f);
  mpc_clear (fd1);
  mpc_clear (fd2);
  mpc_clear (x);
}
