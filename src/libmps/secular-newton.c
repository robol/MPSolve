/*
 * newton.c
 *
 *  Created on: 15/giu/2011
 *      Author: leonardo
 */

#include <mps/secular.h>
#include <mps/interface.h>
#include <mps/link.h>
#include <mps/debug.h>

void
mps_secular_fnewton(mps_status* s, cplx_t x, double *rad, cplx_t corr,
    mps_boolean * again)
{
  int i;
  cplx_t ctmp, ctmp2, pol, fp, sumb;
  double dtmp;
  *again = true;

  mps_secular_equation* sec = (mps_secular_equation*) s->user_data;

  cplx_set(pol, cplx_zero);
  cplx_set(fp, cplx_zero);
  cplx_set(sumb, cplx_zero);

  for (i = 0; i < sec->n; i++)
    {
      /* Compute z - b_i */
      cplx_sub(ctmp, x, sec->bfpc[i]);

      if (cplx_eq_zero(ctmp))
        {
          *again = false;
          cplx_set(corr, cplx_zero);
          return;
        }

      /* Compute (z-b_i)^{-1} */
      cplx_inv_eq(ctmp);

      /* Compute sum of (z-b_i)^{-1} */
      cplx_add_eq(sumb, ctmp);

      /* Compute a_i / (z - b_i) */
      cplx_mul(ctmp2, sec->afpc[i], ctmp);

      /* Add a_i / (z - b_i) to pol */
      cplx_add_eq(pol, ctmp2);

      /* Compute a_i / (z - b_i)^2 */
      cplx_mul_eq(ctmp2, ctmp);

      /* Add it to fp */
      cplx_sub_eq(fp, ctmp2);
    }

  /* Compute secular function */
  cplx_sub_eq(pol, cplx_one);

  /* If S(z) is the secular equation and
   * |S(z)| < eps => |z - z_0| < eps(1 + u) + (n+1)u
   * where z_0 is the real root and u the machine precision,
   * that we can assume that z_0 is a pseudo root, i.e. solution
   * to a problem with small perturbed coefficients. */
  dtmp = cplx_mod(pol) * (DBL_EPSILON + 1) + (s->n + 1);

  /* Compute newton correction */
  cplx_mul(ctmp, pol, sumb);
  cplx_add_eq(fp, ctmp);

  if (!cplx_eq(fp, cplx_zero))
    cplx_div(corr, pol, fp);
  else
    cplx_set(corr, pol);

  /* dtmp here is the guaranteed upper bound to the evaluation of the
   * secular equation   */
  if (dtmp < 2 * DBL_EPSILON)
    {
      *again = false;
    }

  /* Radius is n * newt_corr, if it's better that Gerschgorin's one */
  dtmp = cplx_mod(corr) * sec->n;
  *rad = dtmp;

  /* If the correction is not useful in the current precision do
   * not iterate more   */
  if (*again && (cplx_mod(corr) < cplx_mod(x) * DBL_EPSILON))
    {
      *again = false;
    }

}

void
mps_secular_dnewton(mps_status* s, cdpe_t x, rdpe_t rad, cdpe_t corr,
    mps_boolean * again)
{
  int i;
  *again = true;

  mps_secular_equation* sec = (mps_secular_equation*) s->user_data;

  cdpe_t pol, fp, sumb, ctmp, ctmp2;
  rdpe_t rtmp, rtmp2, apol;

  cdpe_set(pol, cdpe_zero);
  cdpe_set(fp, cdpe_zero);
  cdpe_set(sumb, cdpe_zero);
  rdpe_set(rad, rdpe_zero);
  rdpe_set(apol, rdpe_zero);

  for (i = 0; i < sec->n; i++)
    {
      /* Compute z - b_i */
      cdpe_sub(ctmp, x, sec->bdpc[i]);

      if (cdpe_eq_zero(ctmp))
        {
          cdpe_set(corr, cdpe_zero);
          *again = false;
          return;
        }

      /* Invert it, i.e. compute 1 / (z - b_i) */
      cdpe_inv_eq(ctmp);

      /* Compute sum of 1 / (z - b_i) */
      cdpe_add_eq(sumb, ctmp);

      /* Compute a / (z - b_i) and its modulus */
      cdpe_mul(ctmp2, sec->adpc[i], ctmp);
      cdpe_add_eq(pol, ctmp2);
      cdpe_mod(rtmp, ctmp2);
      rdpe_add_eq(apol, rtmp);

      /* Compute a / (z - b_i)^2 and add it to the first derivative */
      cdpe_mul_eq(ctmp2, ctmp);
      cdpe_sub_eq(fp, ctmp2);
    }

  /* Compute poly */
  cdpe_sub_eq(pol, cdpe_one);

  /* Compute correction */
  cdpe_mul(ctmp, pol, sumb);
  cdpe_add_eq(fp, ctmp);

  if (!cdpe_eq(fp, cdpe_zero))
    {
      cdpe_div(corr, pol, fp);
    }
  else
    {
      cdpe_set(corr, pol);
    }

  /* Compute radius as n * | corr | */
  cdpe_mod(rad, corr);
  rdpe_mul_eq_d(rad, s->n);

  /* Compute \sum_i | a_i / (z - b_i) | + 1
   * and check if the secular equation is smaller
   * than this multiplied for n * epsilon */
  rdpe_add_eq(apol, rdpe_one);
  rdpe_mul_eq_d(apol, (sec->n + 3) * DBL_EPSILON);

  cdpe_mod(rtmp, pol);
  //    if (rdpe_lt(rtmp, apol))
  //            *again = false;

  /* If newton correction is less than
   * the modules of |x| multiplied for
   * for epsilon stop */
  if (*again)
    {
      /* Computation of |x| and |corr| */
      cdpe_mod(rtmp, corr);
      cdpe_mod(rtmp2, x);
      rdpe_mul_eq_d(rtmp2, sec->n * DBL_EPSILON);

      /* If |corr| < |x| * DBL_EPSILON then stop */
      if (rdpe_lt(rtmp, rtmp2))
        *again = false;
    }

}

void
mps_secular_mnewton(mps_status* s, mpc_t x, rdpe_t rad, mpc_t corr,
    mps_boolean * again)
{
  int i, j;
  mps_boolean x_is_b = false;

  /* Set again to true. If the convergence will be proved
   * during the iteration it will be set to false */
  *again = true;

  /* Get a pointer to the secular equation */
  mps_secular_equation* sec = (mps_secular_equation*) s->user_data;

  /* Declare temporary variables */
  mpc_t sumb, pol, fp, ctmp, ctmp2;
  cdpe_t cdtmp;
  rdpe_t rtmp, rtmp2;

  /* Set working precision */
  mpc_init2(sumb, s->mpwp);
  mpc_init2(pol, s->mpwp);
  mpc_init2(fp, s->mpwp);
  mpc_init2(ctmp, s->mpwp);
  mpc_init2(ctmp2, s->mpwp);

  /* Adjust precision of coefficients */
  if (s->mpwp != mpc_get_prec(sec->ampc[0]))
    {
      for (i = 0; i < sec->n; i++)
        {
          mpc_set_prec(sec->ampc[i], s->mpwp);
          mpc_set_prec(sec->bmpc[i], s->mpwp);
        }
    }

  /* Set some starting values */
  mpc_set_d(sumb, 0, 0);
  mpc_set_d(pol, 0, 0);
  mpc_set_d(fp, 0, 0);

  for (i = 0; i < sec->n; i++)
    {
      /* Compute z - b_i */
      mpc_sub(ctmp, x, sec->bmpc[i]);

      /* Keep away the case where the difference is zero */
      if (mpc_eq_zero(ctmp))
        {
          /* We are in the case where x = b_i so set s = i*/
          j = i;
          mpc_set_ui(corr, 0U, 0U);
          for (i = 0; i < sec->n; i++)
            {
              if (i == j)
                continue;
              mpc_sub(ctmp, sec->bmpc[j], sec->bmpc[i]);
              mpc_add(sumb, sec->ampc[i], sec->ampc[j]);
              mpc_div_eq(sumb, ctmp);
              mpc_add_eq(corr, sumb);
            }

          mpc_sub_eq_ui(corr, 1U, 0U);
          mpc_inv_eq(corr);
          mpc_mul_eq(corr, sec->ampc[j]);

          MPS_DEBUG_MPC(s, 15, sec->ampc[j], "sec->ampc[%d]", j)
          MPS_DEBUG_MPC(s, 15, corr, "Correction")

          x_is_b = true;
          break;
        }

      /* Compute (z-b_i)^{-1} */
      mpc_inv_eq(ctmp);

      /* Multiply sum of (z-b_i)^{-1} */
      mpc_add_eq(sumb, ctmp);

      /* Compute a_i / (z - b_i)  */
      mpc_mul(ctmp2, sec->ampc[i], ctmp);

      /* Add a_i / (z - b_i) to pol */
      mpc_add_eq(pol, ctmp2);

      /* Compute a_i / (z - b_i)^2 */
      mpc_mul_eq(ctmp2, ctmp);

      /* Add it to fp */
      mpc_sub_eq(fp, ctmp2);
    }

  /* If x != b_i for every b_i finalize the computation */
  if (!x_is_b)
    {
      /* Subtract one from pol */
      mpc_sub_eq_ui(pol, 1, 0);

      /* Compute correction */
      mpc_mul_eq(sumb, pol);
      mpc_add_eq(fp, sumb);
      if (!mpc_eq_zero(fp))
        mpc_div(corr, pol, fp);
      else
        mpc_set(corr, pol);
    }

  /* Compute radius */
  mpc_get_cdpe(cdtmp, corr);
  cdpe_mod(rtmp2, cdtmp);
  rdpe_mul_eq_d(rtmp2, (double) sec->n);

  rdpe_set(rad, rtmp2);

  /* Compute guaranteed modulus of pol */
  mpc_get_cdpe(cdtmp, pol);
  cdpe_mod(rtmp, cdtmp);
  rdpe_add(rtmp2, s->mp_epsilon, rdpe_one);
  rdpe_mul_eq_d(rtmp2, 1 + s->n);
  rdpe_mul_eq(rtmp, rtmp2);

  /* Epsilon for us */
  rdpe_mul_d(rtmp2, s->mp_epsilon, 2);

  /* If |S(x)| < eps stop */
  if (rdpe_lt(rtmp, rtmp2))
    *again = false;

  /* Check if newton correction is less than
   * the modules of x for s->prec_out, and if
   * that's the case, stop. */
  if (*again)
    {
      mpc_get_cdpe(cdtmp, x);
      cdpe_mod(rtmp, cdtmp);
      rdpe_mul_eq(rtmp, s->mp_epsilon);

      mpc_get_cdpe(cdtmp, corr);
      cdpe_mod(rtmp2, cdtmp);

      if (rdpe_lt(rtmp2, rtmp))
        {
          *again = false;
        }
    }

  mpc_clear(fp);
  mpc_clear(pol);
  mpc_clear(sumb);
  mpc_clear(ctmp);
  mpc_clear(ctmp2);
}
