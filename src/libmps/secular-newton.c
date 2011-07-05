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
#include <math.h>

void
mps_secular_fnewton(mps_status* s, cplx_t x, double *rad, cplx_t corr,
    mps_boolean * again)
{
  int i;
  cplx_t ctmp, ctmp2, pol, fp, sumb;
  double dtmp, g_corr;
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
  cplx_div(corr, pol, fp);

  MPS_DEBUG_CPLX(s, corr, "pol/fp");
  cplx_mul(ctmp, corr, sumb);
  MPS_DEBUG_CPLX(s, ctmp, "pre-real sigma");
  cplx_add_eq(ctmp, cplx_one);
  MPS_DEBUG_CPLX(s, ctmp, "real sigma");

  if (!cplx_eq(ctmp, cplx_zero))
    cplx_div_eq(corr, ctmp);


  /* We compute the following values in order to give a guaranteed
   * Newton inclusion circle:
   *
   * 1) theta = |fp| * (n + 8 + 3*sqrt(2)) * u
   *  This is used so fp - theta is less than the real
   *  value of sum_i a_i / (x - b_i)^2
   *
   * 2) ssp = pol / (fp - theta);
   *  This is the guaranteed S/S' that we can compute
   *  and that will be used to compute
   *
   * 3) gamma = ssp*sec*(n + 2 + sqrt(2))*u +
   *     sec * (pol/fp - ssp) =
   *     = sec * (ssp * (n+2+sqrt(2))*u + (pol/fp - ssp))
   *  This is used in the next step to give a minoration
   *  sigma of 1 + S/s' * (\sum_i a_i / (x-b_i))
   *
   * 4) sigma = |1 + ssp*sumb| - gamma
   *  Guaranteed 1 + S/s' * (\sum_i a_i / (x-b_i))
   *
   * 5) g_corr = ssp / sigma
   *  That is, finally, the guaranteed newton correction.
   */
  {

    double theta, gamma, sigma;
    cplx_t ssp, pol_div_fp, gamma_tmp, sigma_tmp;

    theta = cplx_mod (fp) * (s->n + 8 + 3 * sqrt(2)) * DBL_EPSILON;

    /* Compute ssp */
    cplx_sub(ssp, fp, cplx_d(theta,0));
    cplx_inv_eq(ssp);
    cplx_mul_eq(ssp, pol);

    /* Compute gamma */
    cplx_mul_d(gamma_tmp, ssp, (s->n + 2 + sqrt(2)) * DBL_EPSILON);
    cplx_div(pol_div_fp, pol, fp);
    cplx_sub_eq(pol_div_fp, ssp);
    gamma = cplx_mod(pol) * (cplx_mod(gamma_tmp) + cplx_mod(pol_div_fp));

    /* Computation of ssp * sumb in order to compute sigma */
    cplx_mul(sigma_tmp, ssp, sumb);
    cplx_add_eq(sigma_tmp, cplx_one);

    sigma = fabs(cplx_mod(sigma_tmp)) - gamma;

    if (sigma > 0)
        g_corr = cplx_mod(ssp) / sigma;
    else
        g_corr = DBL_MAX;
  }

  MPS_DEBUG(s, "Computed guaranteed newton radius = %e", g_corr * s->n);
  MPS_DEBUG(s, "Non guaranteed newton radius: %e", cplx_mod(corr) * s->n);

  /* Radius is n * newt_corr, if it's better that Gerschgorin's one */
  dtmp = g_corr * sec->n + DBL_EPSILON;
  if (dtmp < *rad)
    *rad = dtmp;

  /* dtmp here is the guaranteed upper bound to the evaluation of the
   * secular equation   */
  if (dtmp < 2 * DBL_EPSILON)
      *again = false;

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

  if (s->computation_style == 'm')
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
  cdpe_t cdtmp, cdtmp2;
  rdpe_t rtmp, rtmp2, g_corr;

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
          /* We are in the case where x = b_i so set j = i*/
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

  /* Get in cdtmp the sum of a_i / (z - b_i) that will be useful later */
  mpc_get_cdpe(cdtmp, pol);
  cdpe_mod(rtmp, cdtmp);

  /* If x != b_i for every b_i finalize the computation */
  if (!x_is_b)
    {
      /* Subtract one from pol */
      mpc_sub_eq_ui(pol, 1, 0);

      /* Compute correction */
      mpc_mul(ctmp2, sumb, pol);
      mpc_add(ctmp, fp, ctmp2);
      if (!mpc_eq_zero(ctmp))
        mpc_div(corr, pol, ctmp);
      else
        mpc_set(corr, pol);
    }


  /* We compute the following values in order to give a guaranteed
   * Newton inclusion circle:
   *
   * 1) theta = |fp| * (n + 8 + 3*sqrt(2)) * u
   *  This is used so fp - theta is less than the real
   *  value of sum_i a_i / (x - b_i)^2
   *
   * 2) ssp = pol / (fp - theta);
   *  This is the guaranteed S/S' that we can compute
   *  and that will be used to compute
   *
   * 3) gamma = ssp*sec*(n + 2 + sqrt(2))*u +
   *     sec * (pol/fp - ssp) =
   *     = sec * (ssp * (n+2+sqrt(2))*u + (pol/fp - ssp))
   *  This is used in the next step to give a minoration
   *  sigma of 1 + S/s' * (\sum_i a_i / (x-b_i))
   *
   * 4) sigma = 1 + ssp*sumb - gamma
   *  Guaranteed 1 + S/s' * (\sum_i a_i / (x-b_i))
   *
   * 5) g_corr = ssp / sigma
   *  That is, finally, the guaranteed newton correction.
   */
  {
    rdpe_t theta, ssp, gamma, sigma;
    rdpe_t fp_mod, pol_mod, sumb_mod;
    cdpe_t cdpe_tmp;

    /* Get the modulus of fp */
    mpc_get_cdpe(cdpe_tmp, fp);
    cdpe_mod(fp_mod, cdpe_tmp);

    /* Compute theta */
    rdpe_mul_d(theta, s->mp_epsilon, s->n + 8 + 3*sqrt(2));
    rdpe_mul_eq(theta, fp_mod);

    /* Compute ssp */
    mpc_get_cdpe(cdpe_tmp, fp);
    cdpe_mod(pol_mod, cdpe_tmp);
    rdpe_sub(ssp, fp_mod, theta);
    rdpe_inv_eq(ssp);
    rdpe_mul_eq(ssp, pol_mod);

    /* Compute gamma */
    {
      rdpe_t tmp;
      /* Compute -pol/fp */
      rdpe_div(gamma, pol_mod, fp_mod);
      rdpe_mul_eq_d(gamma, -1);

      /* And ssp - pol/fp and store it in tmp */
      rdpe_add_eq(gamma, ssp);

      rdpe_mul_d(tmp, ssp, s->n + 2 + sqrt(2));
      rdpe_mul_eq(tmp, s->mp_epsilon);

      /* Finalize computation */
      rdpe_add_eq(gamma, tmp);
      rdpe_mul_eq(gamma, pol_mod);
    }

    /* Compute sigma */
    mpc_get_cdpe(cdpe_tmp, sumb);
    cdpe_mod(sumb_mod, cdpe_tmp);
    rdpe_mul(sigma, sumb_mod, ssp);
    rdpe_sub_eq(sigma, gamma);
    rdpe_add_eq(sigma, rdpe_one);

    /* Compute g_corr */
    rdpe_div(g_corr, ssp, sigma);
  }

  /* Compute non-guaranteed newton correction */
  {
      cdpe_t ctmp;
      mpc_get_cdpe(ctmp, corr);
      cdpe_mod(rtmp, ctmp);
      rdpe_mul_eq_d(rtmp, s->n);
  }

  /* Radius is s->n * g_corr */
  rdpe_mul_eq_d(g_corr, s->n);

  MPS_DEBUG_RDPE(s, rtmp, "Non-guaranteed newton correction");
  MPS_DEBUG_RDPE(s, g_corr, "Computed newton correction");

  /* Set the radius, if convenient. */
  if (rdpe_lt(g_corr, s->drad[i]))
    rdpe_set(rad, g_corr);

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

  /* Final cleanup */
  mpc_clear(fp);
  mpc_clear(pol);
  mpc_clear(sumb);
  mpc_clear(ctmp);
  mpc_clear(ctmp2);
}
