/*
 * newton.c
 *
 *  Created on: 15/giu/2011
 *      Author: leonardo
 */

#include <mps/types.h>
#include <mps/secular.h>
#include <mps/core.h>
#include <mps/link.h>
#include <mps/debug.h>
#include <limits.h>
#include <math.h>

#define MPS_2SQRT2 2.82842712474619009760f

void
mps_secular_fnewton (mps_status * s, cplx_t x, double *rad, cplx_t corr,
                     mps_boolean * again, void * user_data,
		     mps_boolean skip_radius_computation)
{
  int i;
  cplx_t ctmp, ctmp2, pol, fp, sumb;
  double apol, prod_b = 1.0, new_rad = 0.0f;
  double asum = 0.0f, asum_on_apol;
  mps_secular_iteration_data * data = user_data;
  mps_secular_equation *sec = (mps_secular_equation *) s->secular_equation;

  /* First set again to true */
  *again = true;

  cplx_set (pol, cplx_zero);
  cplx_set (fp, cplx_zero);
  cplx_set (sumb, cplx_zero);

  for (i = 0; i < sec->n; i++)
    {
      /* Compute z - b_i */
      cplx_sub (ctmp, x, sec->bfpc[i]);

      /* Check if we are in the case where z == b_i and return,
       * without doing any further iteration */
      if (cplx_eq_zero (ctmp))
	{
	  *again = false;
          return;
	}

      /* Computation of prod [ (z - b_i) / (z - z_j) ] */
      prod_b *= cplx_mod (ctmp);
      cplx_sub (ctmp2, x, s->froot[i]);
      if (!cplx_eq_zero (ctmp2))
        {
	  prod_b /= cplx_mod (ctmp2);
        }

      /* Compute (z-b_i)^{-1} */
      cplx_inv_eq (ctmp);

      /* Compute sum of (z-b_i)^{-1} */
      cplx_add_eq (sumb, ctmp);

      /* Compute a_i / (z - b_i) */
      cplx_mul (ctmp2, sec->afpc[i], ctmp);
      
      /* Compute the sum of module of (a_i/(z-b_i)) * (i + 2) */
      asum += cplx_mod (ctmp2) * (i + 2);

      /* Add a_i / (z - b_i) to pol */
      cplx_add_eq (pol, ctmp2);

      /* Compute a_i / (z - b_i)^2a */
      cplx_mul_eq (ctmp2, ctmp);

      /* Add it to fp */
      cplx_sub_eq (fp, ctmp2);
    }

  /* Compute secular function */
  cplx_sub_eq (pol, cplx_one);

  /* Compute the module of pol */
  apol = cplx_mod (pol);

  if (apol < 0)
    {
      if (data)
	{
	  MPS_DEBUG (s, "Setting again to false on root %ld for root neighbourhood", data->k);
	}
      *again = false;
      return;
    }


  /* Compute newton correction */
  cplx_mul (corr, pol, sumb);
  cplx_add_eq (corr, fp);
  if (cplx_eq_zero (corr))
      cplx_set (corr, pol);
  else
    cplx_div (corr, pol, corr);

  asum_on_apol = asum / apol;

  /* If the approximation falls in the root neighbourhood then we can stop */
  if ((asum_on_apol + 1) * MPS_2SQRT2 * DBL_EPSILON > 1 ||
      (asum_on_apol < 0))
    {
      if (data && s->debug_level & MPS_DEBUG_APPROXIMATIONS)
	{
	  MPS_DEBUG (s, "Setting again to false on root %ld for root neighbourhood", data->k);
	}
      *again = false;
    }

  /* If the correction is not useful in the current precision do
   * not iterate more */
  if (*again && (cplx_mod (corr) < 8.0 * cplx_mod (x) * DBL_EPSILON))
    {
      if (data && s->debug_level & MPS_DEBUG_APPROXIMATIONS)
	{
	  MPS_DEBUG (s, "Setting again to false on root %ld for small Newton correction", data->k);
	}
      *again = false;
    }

  if (!*again || skip_radius_computation)
    {
      return;
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
   * 4) sigma = |1 + ssp*sumb| - gamma
   *  Guaranteed 1 + S/s' * (\sum_i a_i / (x-b_i))
   *
   * 5) g_corr = ssp / sigma
   *  That is, finally, the guaranteed newton correction.
   */
  if (!skip_radius_computation)
    {
      double theta, gamma, sigma, g_corr;
      cplx_t ssp, pol_div_fp, gamma_tmp, sigma_tmp, c_theta;

      theta = cplx_mod (fp) * (s->n + 8 + 3 * sqrt (2)) * DBL_EPSILON;

      /* Compute ssp */
      cplx_set_d (c_theta, theta, 0);
      cplx_sub (ssp, fp, c_theta);
      cplx_inv_eq (ssp);
      cplx_mul_eq (ssp, pol);

      /* Compute gamma */
      cplx_mul_d (gamma_tmp, ssp, (s->n + 2 + sqrt (2)) * DBL_EPSILON);
      cplx_div (pol_div_fp, pol, fp);
      cplx_sub_eq (pol_div_fp, ssp);
      gamma = cplx_mod (pol) * (cplx_mod (gamma_tmp) + cplx_mod (pol_div_fp));

      /* Computation of ssp * sumb in order to compute sigma */
      cplx_mul (sigma_tmp, ssp, sumb);
      cplx_add_eq (sigma_tmp, cplx_one);

      sigma = cplx_mod (sigma_tmp) - gamma;

      if (sigma > 0)
        {
          g_corr = cplx_mod (ssp) / sigma;
	  new_rad = g_corr * s->n;
        }

      /* Computation of radius with Gerschgorin */
      new_rad += cplx_mod (x) * DBL_EPSILON;

      /* Correct the old radius with the move that we are doing
       * and check if the new proposed radius is preferable. */
      if (new_rad < *rad || (*rad == 0) || (!data))
	*rad = new_rad;
    }
}

void
mps_secular_dnewton (mps_status * s, cdpe_t x, rdpe_t rad, cdpe_t corr,
                     mps_boolean * again, void * user_data,
		     mps_boolean skip_radius_computation)
{
  int i;

  mps_secular_equation *sec = (mps_secular_equation *) s->secular_equation;
  mps_secular_iteration_data * data = user_data;

  cdpe_t pol, fp, sumb, ctmp, ctmp2, old_x;
  rdpe_t rtmp, rtmp2, apol, prod_b, asum, asum_on_apol;

  *again = true;

  cdpe_set (old_x, x);
  cdpe_set (pol, cdpe_zero);
  cdpe_set (fp, cdpe_zero);
  cdpe_set (sumb, cdpe_zero);
  rdpe_set (apol, rdpe_zero);
  rdpe_set (asum, rdpe_zero);
  rdpe_set (prod_b, rdpe_one);

  for (i = 0; i < sec->n; i++)
    {
      /* Compute z - b_i */
      cdpe_sub (ctmp, x, sec->bdpc[i]);

      /* Compute prod [ (z - b_i) / (z - z_j) ] */
      cdpe_mod (rtmp, ctmp);
      rdpe_mul_eq (prod_b, rtmp);
      cdpe_sub (ctmp2, x, s->droot[i]);
      cdpe_mod (rtmp2, ctmp2);
      if (!rdpe_eq_zero (rtmp2))
        {
	  rdpe_div_eq (prod_b, rtmp2);
        }

      /* Alternative computation if x is one of the b_i */
      if (cdpe_eq_zero (ctmp))
	{
	  *again = false;
	  return;
	}

      /* Invert it, i.e. compute 1 / (z - b_i) */
      cdpe_inv_eq (ctmp);

      /* Compute sum of 1 / (z - b_i) */
      cdpe_add_eq (sumb, ctmp);

      /* Compute a / (z - b_i) and its modulus */
      cdpe_mul (ctmp2, sec->adpc[i], ctmp);
      cdpe_add_eq (pol, ctmp2);
      cdpe_mod (rtmp, ctmp2);
      rdpe_mul_eq_d (rtmp, i + 2);
      rdpe_add_eq (asum, rtmp);

      /* Compute a / (z - b_i)^2 and add it to the first derivative */
      cdpe_mul_eq (ctmp2, ctmp);
      cdpe_sub_eq (fp, ctmp2);
    }

  /* Compute poly */
  cdpe_sub_eq (pol, cdpe_one);
  cdpe_mod (apol, pol);

  /* Compute correction */
  cdpe_mul (corr, pol, sumb);
  cdpe_add_eq (corr, fp);
  if (cdpe_eq_zero (corr))
    cdpe_set (corr, pol);
  else
    cdpe_div (corr, pol, corr);

  rdpe_div (asum_on_apol, asum, apol);

  /* If newton correction is less than
   * the modules of |x| multiplied for
   * for epsilon stop */
  /* Computation of |x| and |corr| */
  /* cdpe_mod (rtmp, corr); */
  /* cdpe_mod (rtmp2, x); */
  /* rdpe_mul_eq_d (rtmp2, 4.0 * DBL_EPSILON); */
  
  /* /\* If |corr| < |x| * DBL_EPSILON then stop *\/ */
  /* if (rdpe_lt (rtmp, rtmp2)) */
  /*   { */
  /*     if (data && (s->debug_level & MPS_DEBUG_APPROXIMATIONS)) */
  /* 	{ */
  /* 	  MPS_DEBUG (s, "Setting again on root %ld to false because the Newton correction is too small", data->k); */
  /* 	  MPS_DEBUG_CDPE (s, corr, "Newton correction"); */
  /* 	} */
  /*     *again = false; */
  /*   } */

  rdpe_add (rtmp, rdpe_one, asum_on_apol);
  rdpe_mul_eq_d (rtmp, MPS_2SQRT2 * DBL_EPSILON);
  if (rdpe_ge (rtmp, rdpe_one))
    {
      if (data && (s->debug_level & MPS_DEBUG_APPROXIMATIONS))
	{
	  MPS_DEBUG (s, "Setting again on root %ld to false because the approximation is in the root neighbourhood", data->k);
	}
      *again = false;
    }

  if (!*again || skip_radius_computation)
    return;

  /* Computation of radius with Gerschgorin */
  rdpe_t new_rad;

  /* Compute the guaranteed radius */
  /* rdpe_set (new_rad, apol); */
  /* rdpe_mul_eq_d (new_rad, s->n); */
  /* rdpe_mul_eq (new_rad, prod_b); */

  /* /\* rdpe_set (rtmp, asum); *\/ */
  /* /\* rdpe_div_eq (rtmp, apol); *\/ */
  /* rdpe_add (rtmp, rdpe_one, asum_on_apol); */
  /* rdpe_add_eq_d (rtmp, 3 * s->n); */
  /* rdpe_mul_eq_d (rtmp, DBL_EPSILON); */
  /* rdpe_add_eq (rtmp, rdpe_one); */
  /* rdpe_mul_eq (new_rad, rtmp); */

  /* /\* Correct the old radius with the move that we are doing */
  /*  * and check if the new proposed radius is preferable. *\/ */
  /* if (data) */
  /*   data->radius_set = true; */


  /* Compute radius as n * newt_corr */
  cdpe_mod (new_rad, corr);
  rdpe_mul_eq_d (new_rad, s->n);

  cdpe_mod (rtmp, x);
  rdpe_mul_eq_d (rtmp, 4.0 * DBL_EPSILON);
  rdpe_add_eq (new_rad, rtmp);

  if (rdpe_lt (new_rad, rad) || !data)
    rdpe_set (rad, new_rad);
  else 
    data->radius_set = false;
}

void
mps_secular_mnewton (mps_status * s, mpc_t x, rdpe_t rad, mpc_t corr,
                     mps_boolean * again, void * user_data,
		     mps_boolean skip_radius_computation)
{
  int i;
  mps_secular_iteration_data *data = (mps_secular_iteration_data*) user_data;

  /* Set again to true. If the convergence will be proved
   * during the iteration it will be set to false */
  *again = true;

  /* Get a pointer to the secular equation */
  mps_secular_equation *sec = (mps_secular_equation *) s->secular_equation;

  /* Declare temporary variables */
  mpc_t sumb, pol, fp, ctmp, ctmp2;
  cdpe_t cdtmp, cdtmp2;
  rdpe_t rtmp, rtmp2, apol, asum, prod_b, asum_on_apol;

  /* Set working precision */
  mpc_init2 (sumb, s->mpwp);
  mpc_init2 (pol, s->mpwp);
  mpc_init2 (fp, s->mpwp);
  mpc_init2 (ctmp, s->mpwp);
  mpc_init2 (ctmp2, s->mpwp);

  /* Set some starting values */
  mpc_set_d (sumb, 0, 0);
  mpc_set_d (pol, 0, 0);
  mpc_set_d (fp, 0, 0);
  rdpe_set (prod_b, rdpe_one);

  rdpe_set (asum, rdpe_zero);
  rdpe_set (apol, rdpe_zero);

  for (i = 0; i < sec->n; i++)
    {
      /* Compute z - b_i */
      mpc_sub (ctmp, x, sec->bmpc[i]);

      if (mpc_eq_zero (ctmp))
	{
	  *again = false;
	  goto mnewton_exit;
	}

      mpc_get_cdpe (cdtmp2, ctmp);

      cdpe_mod (rtmp2, cdtmp2);
      rdpe_mul_eq (prod_b, rtmp2);

      mpc_sub (ctmp2, x, s->mroot[i]);
      mpc_get_cdpe (cdtmp2, ctmp2);
      if ((data && (data->k != i)) ||
	  ((!data && (!cdpe_eq_zero (cdtmp2)))))
        {
	  cdpe_mod (rtmp2, cdtmp2);
          rdpe_div_eq (prod_b, rtmp2);
        }

      /* Compute (z-b_i)^{-1} */
      mpc_inv_eq (ctmp);

      /* Add to the sum of (z-b_i)^{-1} */
      mpc_add_eq (sumb, ctmp);

      /* Compute a_i / (z - b_i)  */
      mpc_mul (ctmp2, sec->ampc[i], ctmp);

      /* Compute the sum of |a_i/(z-b_i)| */
      mpc_get_cdpe (cdtmp, ctmp2);
      cdpe_mod (rtmp, cdtmp);
      rdpe_add_eq (asum, rtmp);

      /* Add a_i / (z - b_i) to pol */
      mpc_add_eq (pol, ctmp2);

      /* Compute a_i / (z - b_i)^2 */
      mpc_mul_eq (ctmp2, ctmp);

      /* Add it to fp */
      mpc_sub_eq (fp, ctmp2);
    }

   /* If x != b_i for every b_i finalize the computation */
  /* Subtract one from pol */
  MPS_DEBUG_MPC (s, 10, pol, "pol");
  mpc_sub_eq_ui (pol, 1, 0);
  MPS_DEBUG_MPC (s, 10, pol, "pol");
      
  /* Compute correction */
  mpc_mul (ctmp2, sumb, pol);
  mpc_add (ctmp, fp, ctmp2);
  if (!mpc_eq_zero (ctmp))
    mpc_div (corr, pol, ctmp);
  else
    {
      MPS_DEBUG (s, "The derivative is null!");
      mpc_set (corr, pol);
    }

  /* Get the module of pol */
  mpc_get_cdpe (cdtmp, pol);
  cdpe_mod (apol, cdtmp);

  /* Computation of radius with Gerschgorin */
  rdpe_t new_rad;
  cdpe_t x_cdpe;
  rdpe_t ax;

  mpc_get_cdpe (x_cdpe, x);
  cdpe_mod (ax, x_cdpe);

  /* new_rad = (apol * s->n * prod_b * (1 + (3 * s->n + (asum_on_apol + 1)) * 4 * DBL_EPSILON)) + (cplx_mod (x) * 4 * DBL_EPSILON); */

  /* Compute the guaranteed radius */
  rdpe_mul (rtmp, asum, s->mp_epsilon);
  rdpe_add (new_rad, apol, rtmp);
  rdpe_mul_eq_d (new_rad, s->n);
  rdpe_mul_eq (new_rad, prod_b);

  /* MPS_DEBUG_RDPE (s, prod_b, "prod_b"); */
  /* MPS_DEBUG_RDPE (s, new_rad, "apol * n * prod_b"); */

  rdpe_mul (rtmp2, asum, s->mp_epsilon);
  rdpe_sub (rtmp2, apol, rtmp2);
  if (rdpe_le (rtmp2, rdpe_zero))
    {
      *again = false;
      goto mnewton_exit;
    }

  /* This is asum / apol */
  rdpe_div (asum_on_apol, asum, rtmp2);

  /* rdpe_add (rtmp, asum_on_apol, rdpe_one); */
  /* rdpe_add_eq_d (rtmp, 3 * s->n); */
  /* rdpe_add_eq (rtmp, rdpe_one); */
  /* rdpe_mul_eq (rtmp, s->mp_epsilon); */
  /* rdpe_mul_eq_d (rtmp, 4); */
  /* rdpe_add_eq (rtmp, rdpe_one); */
  /* rdpe_mul_eq (new_rad, rtmp); */

  mpc_get_cdpe (cdtmp, corr);
  cdpe_mod (new_rad, cdtmp);
  rdpe_mul_eq_d (new_rad, s->n);

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
   * 5) g_c orr = ssp / sigma
   *  That is, finally, the guaranteed newton correction.
   */
  /* if (s->mpsolve_ptr == MPS_MPSOLVE_PTR (mps_standard_mpsolve)) */
  {
    rdpe_t theta, ssp, gamma, sigma;
    rdpe_t fp_mod, pol_mod, sumb_mod, g_corr;
    cdpe_t cdpe_tmp, pol_cdpe, fp_cdpe, ssp_cdpe, cdpe_tmp2;

    /* Get the modulus of fp and pol */
    mpc_get_cdpe (fp_cdpe, fp);
    cdpe_mod (fp_mod, fp_cdpe);
    mpc_get_cdpe (pol_cdpe, pol);
    cdpe_mod (pol_mod, pol_cdpe);

    /* Compute theta */
    rdpe_mul_d (theta, s->mp_epsilon, s->n + 2 + 2 * sqrt (2));
    rdpe_mul_eq (theta, fp_mod);

    /* Compute ssp */
    mpc_get_cdpe (cdpe_tmp, fp);
    cdpe_set (cdpe_tmp2, cdpe_zero);
    rdpe_set (cdpe_Re (cdpe_tmp2), theta);
    cdpe_sub_eq (cdpe_tmp, cdpe_tmp2);
    cdpe_div (ssp_cdpe, pol_cdpe, cdpe_tmp);
    cdpe_mod (ssp, ssp_cdpe);

    /* Compute gamma */
    {
      rdpe_t tmp;
      cdpe_t pol_div_fp;

      /* Compute pol/fp */
      cdpe_div (pol_div_fp, pol_cdpe, fp_cdpe);
      cdpe_sub_eq (pol_div_fp, ssp_cdpe);
      rdpe_mul_d (gamma, s->mp_epsilon, s->n + 2 + sqrt (2));
      rdpe_mul_eq (gamma, ssp);
      cdpe_mod (tmp, pol_div_fp);
      rdpe_add_eq (gamma, tmp);
      rdpe_mul_eq (gamma, pol_mod);
    }

    /* Compute sigma */
    mpc_get_cdpe (cdpe_tmp, sumb);
    cdpe_mod (sumb_mod, cdpe_tmp);
    rdpe_mul (sigma, sumb_mod, ssp);
    rdpe_sub_eq (sigma, gamma);
    rdpe_add_eq (sigma, rdpe_one);

    /* Compute g_corr */
    rdpe_div (g_corr, ssp, sigma);

    /* Compute non-guaranteed newton correction */
    mpc_get_cdpe (cdpe_tmp, corr);
    cdpe_mod (rtmp, cdpe_tmp);
    rdpe_mul_eq_d (rtmp, s->n);

    /* Radius is s->n * g_corr */
    rdpe_mul_eq_d (g_corr, s->n);
    if (rdpe_eq_zero (g_corr))
      {
        rdpe_set_2dl (g_corr, 1.0, LONG_MIN);
      }

    /* Set the radius, if convenient. */
    if (rdpe_gt (sigma, rdpe_zero))
      {
        /* MPS_DEBUG (s, "Setting newton correction");
           MPS_DEBUG_RDPE (s, sigma, "sigma"); */
        rdpe_set (new_rad, g_corr);
      }
  }

  rdpe_mul (rtmp, ax, s->mp_epsilon);
  rdpe_add_eq (new_rad, rtmp);

  if (rdpe_lt (new_rad, rad))
    {
      /* MPS_DEBUG_RDPE (s, asum, "asum"); */
      /* MPS_DEBUG_RDPE (s, apol, "apol"); */
      /* MPS_DEBUG_RDPE (s, new_rad, "Setting rad"); */
      rdpe_set (rad, new_rad);
      /* MPS_DEBUG_MPC (s, 40, s->mroot[data->k], "s->mroot[%ld]", data->k); */
      /* MPS_DEBUG_MPC (s, 40, corr, "corr"); */
    }
  
  /* Check if newton correction is less than
   * the modules of x for s->output_config->prec, and if
   * that's the case, stop. */
  if (*again)
     {
       /* Check if the newton correction is small enough */
       rdpe_mul (rtmp, ax, s->mp_epsilon);

       mpc_get_cdpe (cdtmp, corr);
       cdpe_mod (rtmp2, cdtmp);
       
       if (rdpe_lt (rtmp2, rtmp))
	 {
	   if (data) 
	     {
	       MPS_DEBUG (s, "Stopping the iterations on root %ld because newton correction is smaller than machine precision",
			  data->k);
	     }
	   *again = false;
	 }

       /* Or if the secular equation is badly conditioned here */
       rdpe_set (rtmp, asum);
       rdpe_div_eq (rtmp, apol);
       rdpe_add_eq (rtmp, rdpe_one);
       rdpe_mul_eq (rtmp, s->mp_epsilon);
       if (rdpe_ge (rtmp, rdpe_one))
	 {
	   MPS_DEBUG_RDPE (s, rtmp, "Relative error on secular equation");
	   if (data)
	     {
	       MPS_DEBUG (s, "Setting again on root %ld to false because the approximation is in the root neighbourhood", data->k);
	     }
	   *again = false;
	 }       
     }

  if (*again && !skip_radius_computation)
    {
      mpc_get_cdpe (cdtmp, corr);
      cdpe_mod (rad, cdtmp);
      rdpe_mul_eq_d (rad, s->n);
    }

  /* Final cleanup */
 mnewton_exit:
  mpc_clear (fp);
  mpc_clear (pol);
  mpc_clear (sumb);
  mpc_clear (ctmp);
  mpc_clear (ctmp2);
}
