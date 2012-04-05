/*
 * newton.c
 *
 *  Created on: 15/giu/2011
 *      Author: leonardo
 */

#include <mps/mps.h>
#include <limits.h>
#include <math.h>

#define MPS_2SQRT2 2.82842712474619009760

void
mps_secular_fnewton (mps_status * s, cplx_t x, double *rad, cplx_t corr,
                     mps_boolean * again, void * user_data,
		     mps_boolean skip_radius_computation)
{
  int i;
  cplx_t ctmp, ctmp2, pol, fp, sumb;
  double apol, new_rad = 0.0;
  double asum = 0.0, asum_on_apol, ax = cplx_mod (x);
  double asum2 = 0.0, asumb = 0.0;
  double local_error, local_error2;
  mps_secular_iteration_data * data = user_data;
  mps_secular_equation *sec = (mps_secular_equation *) s->secular_equation;

  cplx_t * afpc = data->local_afpc;
  cplx_t * bfpc = data->local_bfpc;

  /* First set again to true */
  *again = true;

  cplx_set (pol, cplx_zero);
  cplx_set (fp, cplx_zero);
  cplx_set (sumb, cplx_zero);

  for (i = 0; i < sec->n; i++)
    {
      /* Compute z - b_i */
      cplx_sub (ctmp, x, bfpc[i]);

      /* Check if we are in the case where z == b_i and return,
       * without doing any further iteration */
      if (cplx_eq_zero (ctmp))
	{
	  *again = false;
          return;
	}

      /* Compute (z-b_i)^{-1} */
      cplx_inv_eq (ctmp);
      if (isinf (cplx_Re (ctmp)) || 
	  isinf (cplx_Re (ctmp)))
	{
	  *again = false;
	  return;
	}

      /* Local error computation */
      local_error = (ax + cplx_mod (bfpc[i])) / (pow (cplx_mod (ctmp), 2));
      local_error2 = local_error * cplx_mod (ctmp);
      asumb += local_error;

      /* Compute sum of (z-b_i)^{-1} */
      cplx_add_eq (sumb, ctmp);

      /* Compute a_i / (z - b_i) */
      cplx_mul (ctmp2, afpc[i], ctmp);
      
      /* Compute the sum of module of (a_i/(z-b_i)) * (i + 2) */
      asum += (local_error + cplx_mod (ctmp2)) * cplx_mod (afpc[i]);

      /* Add a_i / (z - b_i) to pol */
      cplx_add_eq (pol, ctmp2);

      /* Compute a_i / (z - b_i)^2a */
      cplx_mul_eq (ctmp2, ctmp);
      asum2 += (local_error2 + cplx_mod (ctmp2)) * cplx_mod (afpc[i]);

      /* Add it to fp */
      cplx_sub_eq (fp, ctmp2);
    }

  /* Compute secular function */
  cplx_sub_eq (pol, cplx_one);
  asum += 1.0;

  if (isnan (cplx_Re (pol)))
    abort();

  /* Compute the module of pol */
  apol = cplx_mod (pol);

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
      if (data && s->debug_level & MPS_DEBUG_PACKETS)
	MPS_DEBUG (s, "Setting again to false on root %ld for root neighbourhood", data->k);
      *again = false;
    }

  /* If the correction is not useful in the current precision do
   * not iterate more */
  if (*again && (cplx_mod (corr) < sec->n * ax * DBL_EPSILON))
    {
      if (data && s->debug_level & MPS_DEBUG_PACKETS)
	MPS_DEBUG (s, "Setting again to false on root %ld for small Newton correction", data->k);
      *again = false;
    }

  /* We compute the following values in order to give a guaranteed
   * Newton inclusion circle. */
  if (*again && !skip_radius_computation)
    { 
      double g_pol = apol + (DBL_EPSILON * asum * MPS_2SQRT2);
      double g_fp = cplx_mod (fp) - (DBL_EPSILON * asum2 * MPS_2SQRT2);
      double g_sumb = cplx_mod (sumb) + (DBL_EPSILON * asumb * MPS_2SQRT2);

      /* pthread_mutex_lock (data->gs_mutex); */
      /* MPS_DEBUG (s, "pol_%ld = %e", data->k, cplx_mod (pol)); */
      /* MPS_DEBUG (s, "fp_%ld = %e", data->k, cplx_mod (fp)); */
      /* MPS_DEBUG (s, "sumb_%ld = %e", data->k, cplx_mod (sumb)); */
      /* MPS_DEBUG (s, "g_pol_%ld = %e", data->k, g_pol); */
      /* MPS_DEBUG (s, "g_fp_%ld = %e", data->k, g_fp); */
      /* MPS_DEBUG (s, "g_sumb_%ld = %e", data->k, g_sumb); */

      new_rad = g_pol / (g_fp - g_pol * g_sumb) * sec->n + ax * 4.0 * DBL_EPSILON;

      /* MPS_DEBUG (s, "new_rad_%ld = %e", data->k, new_rad); */
      /* pthread_mutex_unlock (data->gs_mutex); */

      if (*again && new_rad < *rad && !(g_fp < 0 || new_rad < 0))
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
  rdpe_t rtmp, rtmp2, apol, asum, asum_on_apol;
  rdpe_t asumb, asum2, ax;

  *again = true;

  cdpe_set (old_x, x);
  cdpe_set (pol, cdpe_zero);
  cdpe_set (fp, cdpe_zero);
  cdpe_set (sumb, cdpe_zero);
  rdpe_set (apol, rdpe_zero);
  rdpe_set (asum, rdpe_zero);
  rdpe_set (asum2, rdpe_zero);
  rdpe_set (asumb, rdpe_zero);

  cdpe_mod (ax, x);

  for (i = 0; i < sec->n; i++)
    {
      /* Compute z - b_i */
      cdpe_sub (ctmp, x, sec->bdpc[i]);

      /* Compute prod [ (z - b_i) / (z - z_j) ] */
      cdpe_mod (rtmp, ctmp);
      rdpe_mul_eq_d (rtmp, (i + 10));
      rdpe_add_eq (asumb, rtmp);

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
      rdpe_mul_eq_d (rtmp, i + 10);
      rdpe_add_eq (asum, rtmp);

      /* Compute a / (z - b_i)^2 and add it to the first derivative */
      cdpe_mul_eq (ctmp2, ctmp);
      cdpe_mod (rtmp, ctmp2);
      rdpe_mul_eq_d (rtmp, (i + 10));
      rdpe_add_eq (asum2, rtmp);
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
  cdpe_mod (rtmp, corr); 
  cdpe_mod (rtmp2, x); 
  rdpe_mul_eq_d (rtmp2, s->n * DBL_EPSILON); 
  
  /* If |corr| < |x| * DBL_EPSILON then stop */
  if (rdpe_lt (rtmp, rtmp2)) 
    { 
      if (data && (s->debug_level & MPS_DEBUG_PACKETS)) 
   	{ 
   	  MPS_DEBUG (s, "Setting again on root %ld to false because the Newton correction is too small", data->k); 
   	  MPS_DEBUG_CDPE (s, corr, "Newton correction"); 
   	} 
      *again = false; 
    }

  rdpe_add (rtmp, rdpe_one, asum_on_apol);
  rdpe_mul_eq_d (rtmp, MPS_2SQRT2 * DBL_EPSILON);
  if (rdpe_ge (rtmp, rdpe_one))
    {
      if (data && (s->debug_level & MPS_DEBUG_PACKETS))
	{
	  MPS_DEBUG (s, "Setting again on root %ld to false because the approximation is in the root neighbourhood", data->k);
	}
      *again = false;
    }

  /* We compute the following values in order to give a guaranteed
   * Newton inclusion circle. */
  if (*again && !skip_radius_computation)
    { 
      rdpe_t new_rad;
      rdpe_t g_pol, g_fp, g_sumb;
      rdpe_t rtmp;
      
      rdpe_mul_d (g_pol, asum, DBL_EPSILON * MPS_2SQRT2);
      rdpe_add_eq_d (g_pol, 1.0);
      rdpe_mul_eq (g_pol, apol);

      rdpe_mul_d (g_fp, asum2, -DBL_EPSILON * MPS_2SQRT2);
      rdpe_add_eq_d (g_fp, 1.0);
      cdpe_mod (rtmp, fp);
      rdpe_mul_eq (g_fp, rtmp);

      rdpe_mul_d (g_sumb, asumb, DBL_EPSILON * MPS_2SQRT2);
      rdpe_add_eq_d (g_sumb, 1.0);
      cdpe_mod (rtmp, sumb);
      rdpe_mul_eq (g_sumb, rtmp);

      rdpe_mul (new_rad, g_pol, g_sumb);
      rdpe_sub (new_rad, g_fp, new_rad);
      rdpe_div (new_rad, g_pol, new_rad);
      rdpe_mul_eq_d (new_rad, s->n);

      MPS_DEBUG_RDPE (s, g_pol, "g_pol");
      MPS_DEBUG_RDPE (s, g_fp, "g_fp");
      MPS_DEBUG_RDPE (s, g_sumb,  "g_sumb");

      rdpe_mul_d (rtmp, ax, DBL_EPSILON * 4.0);
      rdpe_add_eq (new_rad, rtmp);
      
      if (*again && rdpe_lt (new_rad, rad) && !(rdpe_le (g_fp, rdpe_zero) || rdpe_le (new_rad, rdpe_zero)))
	rdpe_set (rad, new_rad);

      MPS_DEBUG_RDPE (s, new_rad, "Computed rad for roots %ld", data->k);
    }
}

void
mps_secular_mnewton (mps_status * s, mpc_t x, rdpe_t rad, mpc_t corr,
		     mps_boolean * again, void * user_data,
		     mps_boolean skip_radius_computation)
{
  int i;
  mpc_t ctmp, ctmp2, pol, fp, sumb;
  cdpe_t cdtmp;
  rdpe_t apol, new_rad, rtmp, ax, afp, rtmp2;
  rdpe_t asum_on_apol, acorr;

  rdpe_t asum, asum2, asumb, diff;
  rdpe_t asum_eps, asum2_eps, asumb_eps;

  rdpe_t local_error, local_error2;

  mpc_t * ampc;
  mpc_t * bmpc;

  /* mps_secular_equation * sec = s->secular_equation; */
  mps_secular_iteration_data * data = user_data;

  ampc = data->local_ampc;
  bmpc = data->local_bmpc;

  /* printf ("x_%ld = ", data->k); mpc_outln_str (stdout, 10, 15, x); printf ("\n"); fflush(stdout); */

  mpc_get_cdpe (cdtmp, x);
  cdpe_mod (ax, cdtmp);

  *again = true;

  mpc_init2 (ctmp, s->mpwp);
  mpc_init2 (ctmp2, s->mpwp);
  mpc_init2 (pol, s->mpwp);
  mpc_init2 (fp, s->mpwp);
  mpc_init2 (sumb, s->mpwp);

  rdpe_set (asum, rdpe_zero);
  rdpe_set (asum2, rdpe_zero);
  rdpe_set (asumb, rdpe_zero);

  for (i = 0; i < s->n; i++)
    {
      /* Compute z - b_i */
      mpc_sub (ctmp, x, bmpc[i]);
      mpc_rmod (diff, ctmp);

      /* Check if we are in the case where x == b_i. If that's
       * the case return without doing anything more */
      if (mpc_eq_zero (ctmp))
	{
	  *again = false;
	  goto mnewton_cleanup;
	}

      /* Invert x - b_i */
      mpc_inv_eq (ctmp);

      mpc_rmod (local_error, bmpc[i]);
      rdpe_add_eq (local_error, ax);
      mpc_rmod (rtmp, ctmp);
      rdpe_inv (local_error2, rtmp);
      rdpe_mul_eq (rtmp, rtmp);
      rdpe_div_eq (local_error, rtmp);
      rdpe_add_eq (asumb, local_error);
      rdpe_mul_eq (local_error2, local_error);
      
      /* Computation of the sum of x - b_i */
      mpc_add_eq (sumb, ctmp);

      /* Sum the module of (1 / (x - b_i)) to asumb. Will be used later
       * for radius computation */
      /* mpc_rmod (rtmp, bmpc[i]); */
      /* rdpe_add_eq (rtmp, ax); */
      /* rdpe_div_eq (rtmp, diff); */
      /* rdpe_add_eq_d (rtmp, (i + 2) + 8); */
      /* rdpe_div_eq (rtmp, diff); */
      /* rdpe_add_eq (asumb, rtmp); */

      /* Compute a_i / (x - b_i) */
      mpc_mul (ctmp2, ctmp, ampc[i]);

      /* Add it to the evaluation of the secular equation, 
       * that we call pol */
      mpc_add_eq (pol, ctmp2);

      mpc_rmod (rtmp, ctmp2);
      rdpe_add_eq (rtmp, local_error);
      mpc_rmod (rtmp2, ampc[i]);
      rdpe_mul_eq (rtmp, rtmp2);
      rdpe_add_eq (asum, rtmp);

      /* Get its module and add it to asum */
      /* mpc_rmod (rtmp, bmpc[i]); */
      /* rdpe_add_eq (rtmp, ax); */
      /* rdpe_div_eq (rtmp, diff); */
      /* rdpe_add_eq_d (rtmp, (i + 2) + 8); */
      /* mpc_rmod (rtmp2, ctmp2); */
      /* rdpe_mul_eq (rtmp, rtmp2); */
      /* rdpe_add_eq (asum, rtmp); */

      /* Computing the derivative S'(x) */
      mpc_mul_eq (ctmp, ctmp2);
      mpc_sub_eq (fp, ctmp);

      /* Add its module at asum2, that will be used later
       * for radius computation */
      /* mpc_rmod (rtmp, bmpc[i]); */
      /* rdpe_add_eq (rtmp, ax); */
      /* rdpe_div_eq (rtmp, diff); */
      /* rdpe_add_eq_d (rtmp, (i  + 2) + 8); */
      /* mpc_rmod (rtmp2, ctmp); */
      /* rdpe_mul_eq (rtmp, rtmp2); */
      /* rdpe_add_eq (asum2, rtmp); */

      mpc_rmod (rtmp, ctmp);
      rdpe_add_eq (rtmp, local_error2);
      rdpe_mul_eq (rtmp, rtmp2);
      rdpe_add_eq (asum2, rtmp);
    }

  /* Finalize the computation of S(x) */
  mpc_sub_eq_ui (pol, 1U, 0U);
  rdpe_add_eq (asum, rdpe_one);

  /* Compute the local error */
  rdpe_mul (asum_eps, asum, s->mp_epsilon);
  rdpe_mul (asum2_eps, asum2, s->mp_epsilon);
  rdpe_mul (asumb_eps, asumb, s->mp_epsilon);

  /* Compute newton correction */
  mpc_mul (ctmp, sumb, pol);
  mpc_add (ctmp, fp, ctmp);
  if (!mpc_eq_zero (ctmp))
    mpc_div (corr, pol, ctmp);
  else
    {
      if (s->debug_level & MPS_DEBUG_PACKETS)
	MPS_DEBUG (s, "The derivative is null!");
      mpc_set (corr, pol);
    }

  /* Compute some modules */
  mpc_rmod (acorr, corr);
  mpc_rmod (apol, pol);
  mpc_rmod (afp, fp);

  /* Compute asum / apol */
  rdpe_div (asum_on_apol, asum, apol);

  /* Check if more iteration on this root are needed or not */
  rdpe_add (rtmp, rdpe_one, asum_on_apol);
  rdpe_mul_eq (rtmp, s->mp_epsilon);
  if (rdpe_gt (rtmp, rdpe_one))
    {
      *again = false;
      if (s->debug_level & MPS_DEBUG_PACKETS)
	MPS_DEBUG (s, "Stopping Aberth iterations due to root neighourhood");
      goto mnewton_cleanup;
    }

  /* Check if the newton correction is small with respect to the
   * current precision. */
  rdpe_mul (rtmp, ax, s->mp_epsilon);
  rdpe_mul_eq_d (rtmp, s->n);
  if (rdpe_lt (acorr, rtmp))
    {
      *again = false;
      if (s->debug_level & MPS_DEBUG_PACKETS)
	MPS_DEBUG (s, "Stopping Aberth iterations due to small Newton correction");
      goto mnewton_cleanup;
    }

  if (!skip_radius_computation)
    {
      rdpe_t g_pol, g_fp, g_sumb, g_den;

      /* We need a guaranteed value of the module of S(x) ... */
      rdpe_add (g_pol, asum_eps, apol);

      /* ... of the module of S'(x), but guaranteed with a lower
       * bound and not upper, ... */
      rdpe_sub (g_fp, afp, asum2_eps);

      /* ... and finally of \sum 1 / (x - b_i) */
      rdpe_add (g_sumb, asumb, asumb_eps);

      /* Compute the guaranteed newton correction. To do this we need a guaranteed
       * denominator for the fraction. Let's do this in MP, so we can have a hope to find
       * a usable value.*/
      {
	mpf_t den, fp_mod, sumb_mod, ftmp;
	mpf_init2 (den, s->mpwp);
	mpf_init2 (fp_mod, s->mpwp);
	mpf_init2 (sumb_mod, s->mpwp);
	mpf_init2 (ftmp, s->mpwp);

	/* Get the guaranteed module of S'(x) */
	mpc_mod (fp_mod, fp);
	mpf_set_rdpe (ftmp, asum2_eps);
	mpf_sub_eq (fp_mod, ftmp); 

	/* And the guaranteed \sum 1 / (x - b_i) */
	mpc_mod (sumb_mod, sumb);
	mpf_set_rdpe (ftmp, asumb_eps);
	mpf_add_eq (sumb_mod, ftmp); 

	/* And the module of S(x) */
	mpc_mod (den, pol);
	mpf_set_rdpe (ftmp, asum_eps); 
	mpf_add_eq (ftmp, den);
	mpf_mul_eq (ftmp, sumb_mod);

	/* Substract them */
	mpf_sub (den, fp_mod, ftmp);
	mpf_get_rdpe (g_den, den);

	/* Getting error on g_den */
	rdpe_mul (rtmp, asumb, apol);
	rdpe_add_eq (rtmp, afp);
	rdpe_mul_eq (rtmp, s->mp_epsilon);

	rdpe_sub_eq (g_den, rtmp);

	mpc_mul (ctmp, pol, sumb);

	mpf_clear (ftmp);
	mpf_clear (sumb_mod);
	mpf_clear (fp_mod);
	mpf_clear (den);
      }

      /* If this is lower than zero, we haven't been able to give a 
       * usable results, so let's quit. */
      if (rdpe_lt (g_den, rdpe_zero))
	{
	  if (s->debug_level & MPS_DEBUG_PACKETS)
	    MPS_DEBUG (s, "Cannot give a guaranteed correction");
	  goto mnewton_cleanup;
	}

      /* In the other case compute the radius */
      rdpe_div (new_rad, g_pol, g_den);
      rdpe_mul_eq_d (new_rad, s->n);

      if (rdpe_lt (new_rad, rad))
	rdpe_set (rad, new_rad);

      /* MPS_DEBUG_MPC (s, 15, x, "x"); */
      /* MPS_DEBUG_RDPE (s, new_rad, "Computed newton radius for root"); */

    }
  
 mnewton_cleanup:

  mpc_clear (pol);
  mpc_clear (fp);
  mpc_clear (sumb);
  mpc_clear (ctmp);
  mpc_clear (ctmp2);
}

