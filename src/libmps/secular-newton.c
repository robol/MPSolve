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
#define KAPPA (sec->n + 7 * 1.4151135)
#define MPS_SQRT2 1.4142135623

void
mps_secular_fnewton (mps_status * s, mps_approximation * root, cplx_t corr,
                     void * user_data,
		     mps_boolean skip_radius_computation)
{
  int i;
  cplx_t ctmp, ctmp2, pol, fp, sumb;
  double apol;
  double asum = 0.0, asum_on_apol, ax = cplx_mod (root->fvalue);
  double asumb = 0.0;
  mps_secular_iteration_data * data = user_data;
  mps_secular_equation *sec = (mps_secular_equation *) s->secular_equation;

  cplx_t * afpc = data->local_afpc;
  cplx_t * bfpc = data->local_bfpc;

  cplx_t x;
  cplx_set (x, root->fvalue);

  /* First set again to true */
  root->again = true;

  cplx_set (pol, cplx_zero);
  cplx_set (fp, cplx_zero);
  cplx_set (sumb, cplx_zero);
  cplx_set (corr, cplx_zero);

  for (i = 0; i < sec->n; i++)
    {
      /* Compute z - b_i */
      cplx_sub (ctmp, x, bfpc[i]);

      /* Check if we are in the case where z == b_i and 
       * if that's the case perform the fallback computation
       * for the Newton correction */
      if (cplx_eq_zero (ctmp))
	{
	  int k;
	  double acorr;
	  double asum = 0;

	  for (k = 0; k < sec->n; k++)
	    {
	      if (i != k)
		{
		  cplx_sub (ctmp, bfpc[i], bfpc[k]);
		  cplx_add (ctmp2, afpc[i], afpc[k]);
		  cplx_div_eq (ctmp2, ctmp);
		  cplx_add_eq (corr, ctmp2);

		  asum += cplx_mod (ctmp2);
		}
	    }

	  if (!cplx_eq_zero (corr))
	    {
	      /* mps_boolean in_root_neighborhood = (asum * DBL_EPSILON * KAPPA > cplx_mod (corr)); */
	      cplx_div (corr, afpc[i], corr);
	      
	      acorr = cplx_mod (corr);
	      if (acorr < ax * DBL_EPSILON)
		{
		  root->again = false;  
		  /* root->approximated = true;    */
		}

	      /* root->frad = acorr * (1 + asum * KAPPA * DBL_EPSILON) * sec->n; */
	    }


	  return;
	}

      /* Compute (z-b_i)^{-1} */
      cplx_inv_eq (ctmp);

      /* Local error computation */
      asumb += cplx_mod (ctmp);

      /* Compute sum of (z-b_i)^{-1} */
      cplx_add_eq (sumb, ctmp);

      /* Compute a_i / (z - b_i) */
      cplx_mul (ctmp2, afpc[i], ctmp);
      
      /* Compute the sum of module of (a_i/(z-b_i)) * (i + 2) */
      asum += cplx_mod (ctmp2);

      /* Add a_i / (z - b_i) to pol */
      cplx_add_eq (pol, ctmp2);

      /* Compute a_i / (z - b_i)^2a */
      cplx_mul_eq (ctmp2, ctmp);

      /* Add it to fp */
      cplx_sub_eq (fp, ctmp2);
    }

  /* Compute secular function */
  cplx_sub_eq (pol, cplx_one);
  asum += 1.0;

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
  if ((asum_on_apol + 1) * KAPPA * DBL_EPSILON > 1)
    {
      if (data && s->debug_level & MPS_DEBUG_PACKETS)
	MPS_DEBUG (s, "Setting again to false on root %ld for root neighbourhood", data->k);
      root->again = false;

      if ((KAPPA * asum / cplx_mod (fp) < 1))
	{
	  if (data && s->debug_level & MPS_DEBUG_PACKETS)
	    {
	      MPS_DEBUG (s, "Setting approximated = true on root %ld for small conditioning number", data->k);
	    }
	  root->approximated = true;
	}
    }
  else if ((cplx_mod (corr) < MPS_SQRT2 * ax * DBL_EPSILON))
    {
      if (data && s->debug_level & MPS_DEBUG_PACKETS)
	MPS_DEBUG (s, "Setting approximated to true on root %ld for small Newton correction", data->k);
      root->again = false;  
      root->approximated = true;
    }

  if (!cplx_eq_zero (corr) && root->again)
    {
      double new_rad = cplx_mod (corr) * (1 + DBL_EPSILON * KAPPA * (asum + 1) / apol) * sec->n;

      if (new_rad < root->frad)
	root->frad = new_rad;
    }
}

void
mps_secular_dnewton (mps_status * s, mps_approximation * root, cdpe_t corr,
                     void * user_data,
		     mps_boolean skip_radius_computation)
{
  int i;

  mps_secular_equation *sec = (mps_secular_equation *) s->secular_equation;
  mps_secular_iteration_data * data = user_data;

  cdpe_t pol, fp, sumb, ctmp, ctmp2, x;
  rdpe_t rtmp, rtmp2, apol, asum, asum_on_apol;
  rdpe_t asumb, ax;

  root->again = true;
  
  cdpe_set (x, root->dvalue);
  cdpe_set (pol, cdpe_zero);
  cdpe_set (fp, cdpe_zero);
  cdpe_set (sumb, cdpe_zero);
  rdpe_set (apol, rdpe_zero);
  rdpe_set (asum, rdpe_zero);
  rdpe_set (asumb, rdpe_zero);

  cdpe_mod (ax, x);

  for (i = 0; i < sec->n; i++)
    {
      /* Compute z - b_i */
      cdpe_sub (ctmp, x, sec->bdpc[i]);

      /* Compute prod [ (z - b_i) / (z - z_j) ] */
      cdpe_mod (rtmp, ctmp);
      rdpe_add_eq (asumb, rtmp);

      /* Check if we are in the case where z == b_i and 
       * if that's the case perform the fallback computation
       * for the Newton correction */
      if (cdpe_eq_zero (ctmp))
	{
	  root->again = true;
	  int k;
	  rdpe_t acorr;
	  rdpe_t sigma;

	  cdpe_set (corr, cdpe_zero);
	  rdpe_set (sigma, rdpe_zero);

	  for (k = 0; k < sec->n; k++)
	    {
	      if (i != k)
		{
		  cdpe_sub (ctmp, sec->bdpc[i], sec->bdpc[k]);
		  cdpe_add (ctmp2, sec->adpc[i], sec->adpc[k]);
		  cdpe_div_eq (ctmp2, ctmp);
		  cdpe_add_eq (corr, ctmp2);

		  cdpe_mod (rtmp, ctmp2);
		  rdpe_add_eq (sigma, rtmp);
		}
	    }

	  cdpe_div (corr, sec->adpc[i], corr);
	  cdpe_mod (acorr, corr); 
	  /* rdpe_mul_d (rtmp, acorr, sec->n * (1 + sec->n * DBL_EPSILON));  */
	  /* if (rdpe_lt (rtmp, s->root[i]->drad))  */
	  /*   rdpe_set (s->root[i]->drad, rtmp);  */

	  if (root->again)
	    {
	      rdpe_mul_d (rtmp, ax, DBL_EPSILON * MPS_SQRT2);
	      if (rdpe_lt (acorr, rtmp))
		{
		  /* root->again = false; */
		}
	    }

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

  rdpe_add (rtmp, rdpe_one, asum_on_apol);
  rdpe_mul_eq_d (rtmp, KAPPA * DBL_EPSILON);
  if (rdpe_ge (rtmp, rdpe_one))
    {
      if (data && (s->debug_level & MPS_DEBUG_PACKETS))
	{
	  MPS_DEBUG (s, "Setting again on root %ld to false because the approximation is in the root neighbourhood", data->k);
	}
      root->again = false;

      cdpe_mod (rtmp, fp);
      rdpe_div (rtmp, asumb, rtmp);
      rdpe_mul_eq_d (rtmp, KAPPA);

      if (rdpe_le (rtmp, rdpe_one))
	{
	  if (data && s->debug_level & MPS_DEBUG_PACKETS)
	    MPS_DEBUG (s, "Setting approximated = true on root %ld for small conditioning number", data->k);
	  root->approximated = true;
	}
    }
  else 
    {      
      cdpe_mod (rtmp, corr); 
      cdpe_mod (rtmp2, x); 
      rdpe_mul_eq_d (rtmp2, MPS_SQRT2 * DBL_EPSILON); 

      /* If |corr| < |x| * DBL_EPSILON then stop */
      if (rdpe_lt (rtmp, rtmp2)) 
	{ 
	  if (data && (s->debug_level & MPS_DEBUG_PACKETS)) 
	    { 
	      MPS_DEBUG (s, "Setting again on root %ld to false because the Newton correction is too small", data->k); 
	      MPS_DEBUG_CDPE (s, corr, "Newton correction"); 
	    } 
	  root->approximated = true;
	}
    }

  if (!cdpe_eq_zero (corr) && root->again)
    {
      rdpe_mul_d (rtmp, asum_on_apol, DBL_EPSILON * KAPPA);
      rdpe_add_eq (rtmp, rdpe_one);
      cdpe_mod (rtmp2, corr);
      rdpe_mul_eq (rtmp2, rtmp);

      if (rdpe_lt (rtmp2, root->drad))
	rdpe_set (root->drad, rtmp2);
    }
}

void
mps_secular_mnewton (mps_status * s, mps_approximation * root, mpc_t corr,
		     void * user_data,
		     mps_boolean skip_radius_computation)
{
  int i;
  mpc_t ctmp, ctmp2, pol, fp, sumb;
  cdpe_t cdtmp;
  rdpe_t apol, rtmp, ax, afp, rtmp2;
  rdpe_t asum_on_apol, acorr;

  rdpe_t asum, asumb, diff;
  rdpe_t asum_eps, asumb_eps;

  mpc_t * ampc;
  mpc_t * bmpc;

  mps_secular_equation * sec = s->secular_equation;
  mps_secular_iteration_data * data = user_data;

  ampc = data->local_ampc;
  bmpc = data->local_bmpc;

  /* printf ("x_%ld = ", data->k); mpc_outln_str (stdout, 10, 15, x); printf ("\n"); fflush(stdout); */

  mpc_get_cdpe (cdtmp, root->mvalue);
  cdpe_mod (ax, cdtmp);

  root->again = true;

  mpc_init2 (ctmp, s->mpwp);
  mpc_init2 (ctmp2, s->mpwp);
  mpc_init2 (pol, s->mpwp);
  mpc_init2 (fp, s->mpwp);
  mpc_init2 (sumb, s->mpwp);

  rdpe_set (asum, rdpe_zero);
  rdpe_set (asumb, rdpe_zero);

  for (i = 0; i < sec->n; i++)
    {
      /* Compute z - b_i */
      mpc_sub (ctmp, root->mvalue, bmpc[i]);
      mpc_rmod (diff, ctmp);

      /* Check if we are in the case where x == b_i. If that's
       * the case return without doing anything more */
      if (mpc_eq_zero (ctmp))
	{
	  /* root->again = true; */
	  int k;
	  /* rdpe_t acorr; */
	  /* rdpe_t sigma; */

	  mpc_set_ui (corr, 0U, 0U);
	  /* rdpe_set (sigma, rdpe_zero); */

	  for (k = 0; k < sec->n; k++)
	    {
	      if (i != k)
		{
		  mpc_sub (ctmp, sec->bmpc[i], sec->bmpc[k]);
		  mpc_add (ctmp2, sec->ampc[i], sec->ampc[k]);
		  mpc_div_eq (ctmp2, ctmp);
		  mpc_add_eq (corr, ctmp2);

		  /* mpc_rmod (rtmp, ctmp2); */
		  /* rdpe_add_eq (sigma, rtmp); */
		}
	    }

	  /* mpc_rmod (rtmp, corr); */
	  /* rdpe_mul_eq_d (sigma, KAPPA); */
	  /* rdpe_mul_eq (sigma, s->mp_epsilon); */
	  /* if (rdpe_lt (rtmp, sigma)) */
	  /*   root->again = false; */

	  mpc_div (corr, sec->ampc[i], corr);
	  /* mpc_rmod (acorr, corr); */
 
	  /* rdpe_mul_d (rtmp, acorr, sec->n * (1 + sec->n));  */
	  /* rdpe_mul_eq (rtmp, s->mp_epsilon); */
	  /* if (rdpe_lt (rtmp, s->root[i]->drad))  */
	  /*   rdpe_set (s->root[i]->drad, rtmp);  */

	  /* rdpe_mul (rtmp, ax, s->mp_epsilon); */
	  /* if (root->again && rdpe_lt (acorr, rtmp)) */
	  /*   root->again = false; */

	  /* rdpe_mul (rtmp, sigma, s->mp_epsilon); */
	  /* rdpe_mul_eq_d (rtmp, KAPPA); */
	  /* rdpe_add_eq (rtmp, rdpe_one); */
	  /* rdpe_mul (rad, rtmp, acorr); */

	  goto mnewton_cleanup;
	}

      /* Invert x - b_i */
      mpc_inv_eq (ctmp);
      mpc_rmod (rtmp, ctmp);
      rdpe_add_eq (asumb, rtmp);
      
      /* Computation of the sum of x - b_i */
      mpc_add_eq (sumb, ctmp);

      /* Compute a_i / (x - b_i) */
      mpc_mul (ctmp2, ctmp, ampc[i]);

      /* Add it to the evaluation of the secular equation, 
       * that we call pol */
      mpc_add_eq (pol, ctmp2);

      mpc_rmod (rtmp, ctmp2);
      rdpe_add_eq (asum, rtmp);

      /* Computing the derivative S'(x) */
      mpc_mul_eq (ctmp, ctmp2);
      mpc_sub_eq (fp, ctmp);
    }

  /* Finalize the computation of S(x) */
  mpc_sub_eq_ui (pol, 1U, 0U);
  rdpe_add_eq (asum, rdpe_one);

  /* Compute the local error */
  rdpe_mul (asum_eps, asum, s->mp_epsilon);
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
  rdpe_mul_eq_d (rtmp, MPS_2SQRT2 * sec->n);
  if (rdpe_gt (rtmp, rdpe_one))
    {
      root->again = false;
      if (s->debug_level & MPS_DEBUG_PACKETS)
	MPS_DEBUG (s, "Stopping Aberth iterations due to root neighborhood");

      mpc_rmod (rtmp2, fp);
      rdpe_div (rtmp, asumb, rtmp2);
      rdpe_mul_eq_d (rtmp, KAPPA);

      if (rdpe_le (rtmp, rdpe_one))
	{
	  if (data && s->debug_level & MPS_DEBUG_PACKETS)
	    MPS_DEBUG (s, "Setting approximated = true on root %ld for small conditioning number", data->k);
	  root->approximated = true;
	}

      goto mnewton_cleanup;
    }
  else 
    {
      /* Check if the newton correction is small with respect to the
       * current precision. */
      rdpe_mul (rtmp, ax, s->mp_epsilon);
      rdpe_mul_eq_d (rtmp, MPS_SQRT2);
      if (rdpe_lt (acorr, rtmp))
	{
	  root->approximated = true;
	  if (s->debug_level & MPS_DEBUG_PACKETS)
	    MPS_DEBUG (s, "Stopping Aberth iterations due to small Newton correction");
	  goto mnewton_cleanup;
	}
    }


  mpc_rmod (rtmp2, corr);
  if (!rdpe_eq_zero (rtmp2) && root->again)
    {
      rdpe_mul_d (rtmp, asum_on_apol, DBL_EPSILON * KAPPA);
      rdpe_add_eq (rtmp, rdpe_one);
      rdpe_mul_eq (rtmp2, rtmp);
      if (rdpe_lt (rtmp2, root->drad))
	rdpe_set (root->drad, rtmp2);
    }
  
 mnewton_cleanup:

  mpc_clear (pol);
  mpc_clear (fp);
  mpc_clear (sumb);
  mpc_clear (ctmp);
  mpc_clear (ctmp2);
}

