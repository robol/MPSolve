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
#define KAPPA (log2(sec->n) + 7 * 1.4151135 + 1)
#define MPS_SQRT2 1.4142135623

/* We need some special codes to identify the meaning of the exit
 * status of mps_secular_fparallel_sum() and its DPE and MP
 * versions. */
#define MPS_PARALLEL_SUM_SUCCESS -1
#define MPS_PARALLEL_SUM_FAILED  -2

/**
 * @brief Perform the evaluation of the Newton correction with the formula
 * obtained implicitly by the secular equation using the parallel
 * algorithm.
 * 
 * @param s The mps_status associated to the current computation
 * @param root The approximation that shall be used as evaluation point
 * @param n The length of the terms that should be summed by the function
 * @param afpc A pointer to the first floating point a_i coefficient
 * @param bfpc A pointer to the first floating point b_i coefficient
 * @param pol The complex value where the result of the evaluation of S(x) will
 *   be stored
 * @param fp The complex value where the result of the evaluation of S'(x) will be
 *   stored
 * @param sumb The complex value where the result of the evaluation of the sum of 
 *   the terms \f$\frac{1}{(x - b_i)}\f$ will be stored
 * @param asum A pointer to a double where the sum of the module of the complex terms
 *   \f$\frac{a_i}{(x - b_i)}\f$ will be saved
 * @param asum2 A pointer to a double where the sum of the module of the complex terms
 *   \f$\frac{a_i}{(x - b_i)^2}\f$ will be stored
 * @param asumb A pointer to a double where the sum of the moduli of the complex terms
 *   \f$\frac{1}{x - b_i}$ will be stored
 * @return If the returned value is a positive integer it is the index of the b_i term 
 *   that is equal to x, suggesting that the alternate evaluation algorithm must be used. 
 *   If it is MPS_PARALLEL_SUM_FAILED then a floting point was encountered in the computation, 
 *   while MPS_PARALLEL_SUM_SUCCESS indicates that the evaluation was successful. 
 */
int
mps_secular_fparallel_sum (mps_status * s, mps_approximation * root, int n, cplx_t * afpc, cplx_t * bfpc,
			   cplx_t pol, cplx_t fp, cplx_t sumb, double * asum, double * asum2, double * asumb)
{
  if (n <= 2)
    {
      int i;
      cplx_t ctmp, ctmp2;
      for (i = 0; i < n; i++)
	{
	  /* Compute z - b_i */
	  cplx_sub (ctmp, root->fvalue, bfpc[i]);
	  
	  /* Check if we are in the case where z == b_i and return,
	   * without doing any further iteration */
	  if (cplx_eq_zero (ctmp))
	    {
	      return i;
	    }

	  /* Compute (z-b_i)^{-1} */
	  cplx_inv_eq (ctmp);
	  if (isinf (cplx_Re (ctmp)) || 
	      isinf (cplx_Re (ctmp)))
	    {
	      root->again = false;	     
	      return MPS_PARALLEL_SUM_FAILED;
	    }

	  /* Compute sum of (z-b_i)^{-1} */
	  cplx_add_eq (sumb, ctmp);

	  /* Add the module of the term to asumb */
	  *asumb += cplx_mod (ctmp);

	  /* Compute a_i / (z - b_i) */
	  cplx_mul (ctmp2, afpc[i], ctmp);

	  /* Compute the sum of module of (a_i/(z-b_i)) * (i + 2) */
	  *asum += cplx_mod (ctmp2);

	  /* Add a_i / (z - b_i) to pol */
	  cplx_add_eq (pol, ctmp2);

	  /* Compute a_i / (z - b_i)^2a */
	  cplx_mul_eq (ctmp2, ctmp);

	  /* Add the module of the term of the derivative to the sum */
	  *asum2 += cplx_mod (ctmp2);

	  /* Add it to fp */
	  cplx_sub_eq (fp, ctmp2);
	}
      
      return MPS_PARALLEL_SUM_SUCCESS;
    }
  else 
    {
      int i = n/2, k;
      if ((k = mps_secular_fparallel_sum (s, root, i, afpc, bfpc, pol, fp, sumb, asum, asum2, asumb)) >= 0)
	{
	  return k;
	}
      if ((k = mps_secular_fparallel_sum (s, root, n-i, afpc + i, bfpc + i, pol, fp, sumb, asum, asum2, asumb)) >= 0)
	{
	  return i + k;
	}
      
      return MPS_PARALLEL_SUM_SUCCESS;
    }
}

void
mps_secular_fnewton (mps_status * s, mps_approximation * root, cplx_t corr,
                     void * user_data,
		     mps_boolean skip_radius_computation)
{
  int i;
  cplx_t ctmp, ctmp2, pol, fp, sumb;
  double apol, acorr, afp;
  double asum = 0.0, asum2 = 0.0, asumb = 0.0, asum_on_apol, ax = cplx_mod (root->fvalue);
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

  if ((i = mps_secular_fparallel_sum (s, root, sec->n, s->secular_equation->afpc, 
				      s->secular_equation->bfpc, pol, 
				      fp, sumb, &asum, &asum2, &asumb)) >= 0)
    {
      int k;
      asum = 0.0;

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

      cplx_sub_eq (corr, cplx_one);

      if (!cplx_eq_zero (corr))
	{
	  cplx_div (corr, afpc[i], corr);
	      
	  acorr = cplx_mod (corr);
	  if (acorr < ax * DBL_EPSILON)
	    {
	      root->again = false;
	      if (asum * KAPPA * DBL_EPSILON < acorr)
		root->approximated = true;
	    }
	}
      else
	root->again = false;
      
      return;
    }

  if (i == MPS_PARALLEL_SUM_FAILED)
    {
      s->root_status[data->k] = MPS_ROOT_STATUS_NOT_FLOAT;
      root->again = false;
      return;
    }

  /* Compute secular function */
  cplx_sub_eq (pol, cplx_one);
  asum += 1.0;

  /* Compute the module of pol */
  apol = cplx_mod (pol);
  afp = cplx_mod (fp);

  /* Compute newton correction */
  cplx_mul (corr, pol, sumb);
  cplx_add_eq (corr, fp);
  if (cplx_eq_zero (corr))
      cplx_set (corr, pol);
  else
    cplx_div (corr, pol, corr);

  asum_on_apol = asum / apol;
  acorr = cplx_mod (corr);

  /* If the approximation falls in the root neighbourhood then we can stop */
  if ((asum_on_apol + 1) * KAPPA * DBL_EPSILON > 1)
    {
      if (data && s->debug_level & MPS_DEBUG_PACKETS)
	MPS_DEBUG (s, "Setting again to false on root %ld for root neighbourhood", data->k);
      root->again = false;

      if ((KAPPA * asum / afp < 1))
	{
	  if (data && s->debug_level & MPS_DEBUG_PACKETS)
	    {
	      MPS_DEBUG (s, "Setting approximated = true on root %ld for small conditioning number", data->k);
	    }
	  root->approximated = true;
	}
    }
  else if (acorr < MPS_SQRT2 * ax * DBL_EPSILON)
    {
      if (data && s->debug_level & MPS_DEBUG_PACKETS)
	MPS_DEBUG (s, "Setting approximated to true on root %ld for small Newton correction", data->k);
      root->again = false;  
      root->approximated = true;
    }

  if (!cplx_eq_zero (corr) && root->again)
    {
      double new_rad = acorr * s->n * (1 + KAPPA * DBL_EPSILON * asum_on_apol);

      if ((new_rad > 0) && (new_rad < root->frad))
	root->frad = new_rad;
    }
}

/**
 * @brief Perform the evaluation of the DPE Newton correction with the formula
 * obtained implicitly by the secular equation using the parallel
 * algorithm.
 * 
 * @param s The mps_status associated to the current computation
 * @param root The approximation that shall be used as evaluation point
 * @param n The length of the terms that should be summed by the function
 * @param adpc A pointer to the first DPE a_i coefficient
 * @param bdpc A pointer to the first DPE b_i coefficient
 * @param pol The CDPE value where the result of the evaluation of S(x) will
 *   be stored
 * @param fp The CDPE value where the result of the evaluation of S'(x) will be
 *   stored
 * @param sumb The CDPE value where the result of the evaluation of the sum of 
 *   the terms \f$\frac{1}{(x - b_i)}\f$ will be stored
 * @param asum The RDPE value where the sum of the module of the complex terms
 *   \f$\frac{a_i}{(x - b_i)}\f$ will be saved
 * @param asum2 The RDPE value where the sum of the module of the complex terms
 *   \f$\frac{a_i}{(x - b_i)^2}\f$ will be stored
 * @param asumb The RDPE value where the sum of the moduli of the complex terms
 *   \f$\frac{1}{x - b_i}$ will be stored
 * @return If the returned value is a positive integer it is the index of the b_i term 
 *   that is equal to x, suggesting that the alternate evaluation algorithm must be used. 
 *   If it is MPS_PARALLEL_SUM_FAILED then a floting point was encountered in the computation, 
 *   while MPS_PARALLEL_SUM_SUCCESS indicates that the evaluation was successful. 
 */
int
mps_secular_dparallel_sum (mps_status * s, mps_approximation * root, int n, cdpe_t * adpc, cdpe_t * bdpc,
			   cdpe_t pol, cdpe_t fp, cdpe_t sumb, rdpe_t asum, rdpe_t asum2, rdpe_t asumb)
{
  if (n <= 2)
    {
      int i;
      cdpe_t ctmp, ctmp2;
      rdpe_t rtmp;
      for (i = 0; i < n; i++)
	{
	  /* Compute z - b_i */
	  cdpe_sub (ctmp, root->dvalue, bdpc[i]);
	  
	  /* Check if we are in the case where z == b_i and return,
	   * without doing any further iteration */
	  if (cdpe_eq_zero (ctmp))
	    {
	      return i;
	    }
	    
	  /* Add the module of 1 / (x - b_i) to asumb */
	  cdpe_mod (rtmp, ctmp);
	  rdpe_add_eq (asumb, rtmp);

	  /* Compute (z-b_i)^{-1} */
	  cdpe_inv_eq (ctmp);

	  /* Compute sum of (z-b_i)^{-1} */
	  cdpe_add_eq (sumb, ctmp);

	  /* Compute a_i / (z - b_i) */
	  cdpe_mul (ctmp2, adpc[i], ctmp);

	  /* Compute the sum of module of (a_i/(z-b_i)) * (i + 2) */
	  cdpe_mod (rtmp, ctmp2);
	  rdpe_add_eq (asum, rtmp);

	  /* Add a_i / (z - b_i) to pol */
	  cdpe_add_eq (pol, ctmp2);

	  /* Compute a_i / (z - b_i)^2a */
	  cdpe_mul_eq (ctmp2, ctmp);
	  
	  /* Add the module of a_i / (x - b_i)^2 to asum2 */
	  cdpe_mod (rtmp, ctmp2);
	  rdpe_add_eq (asum2, rtmp);

	  /* Add it to fp */
	  cdpe_sub_eq (fp, ctmp2);
	}
      
      return MPS_PARALLEL_SUM_SUCCESS;
    }
  else 
    {
      int i = n/2, k;
      if ((k = mps_secular_dparallel_sum (s, root, i, adpc, bdpc, pol, fp, sumb, asum, asum2, asumb)) >= 0)
	{
	  return k;
	}
      if ((k = mps_secular_dparallel_sum (s, root, n-i, adpc + i, bdpc + i, pol, fp, sumb, asum, asum2, asumb)) >= 0)
	{
	  return i + k;
	}
      
      return MPS_PARALLEL_SUM_SUCCESS;
    }
}

void
mps_secular_dnewton (mps_status * s, mps_approximation * root, cdpe_t corr,
                     void * user_data,
		     mps_boolean skip_radius_computation)
{
  int i;
  cdpe_t ctmp, ctmp2, pol, fp, sumb, x;
  rdpe_t apol, asum, asumb, asum2, asum_on_apol, ax, rtmp, rtmp2, acorr;
  mps_secular_iteration_data * data = user_data;
  mps_secular_equation *sec = (mps_secular_equation *) s->secular_equation;

  cdpe_set (x, root->dvalue);
  rdpe_set (asum, rdpe_zero);
  rdpe_set (asumb, rdpe_zero);
  rdpe_set (asum2, rdpe_zero);
  cdpe_mod (ax, x);

  /* First set again to true */
  root->again = true;

  cdpe_set (pol, cdpe_zero);
  cdpe_set (fp, cdpe_zero);
  cdpe_set (sumb, cdpe_zero);
  cdpe_set (corr, cdpe_zero);

  if ((i = mps_secular_dparallel_sum (s, root, sec->n, s->secular_equation->adpc, s->secular_equation->bdpc, 
				       pol, fp, sumb, asum, asum2, asumb)) != MPS_PARALLEL_SUM_SUCCESS)
    {
      int k;
      rdpe_set (asum, rdpe_zero);

      for (k = 0; k < sec->n; k++)
	{
	  if (i != k)
	    {
	      cdpe_sub (ctmp, sec->bdpc[i], sec->bdpc[k]);
	      cdpe_add (ctmp2, sec->adpc[i], sec->adpc[k]);
	      cdpe_div_eq (ctmp2, ctmp);
	      cdpe_add_eq (corr, ctmp2);

	      cdpe_mod (rtmp, ctmp2);
	      rdpe_add_eq (asum, rtmp);
	    }
	}

      cdpe_sub_eq (corr, cdpe_one);

      if (!cdpe_eq_zero (corr))
	{
	  cdpe_div (corr, sec->adpc[i], corr);
	      
	  cdpe_mod (rtmp, corr);
	  rdpe_mul_d (rtmp2, ax, DBL_EPSILON);
	  if (rdpe_lt (rtmp, rtmp2))
	    {
	      root->again = false;
	      
	      rdpe_mul_d (rtmp2, asum, KAPPA * DBL_EPSILON);
	      if (rdpe_lt (rtmp2, rtmp))
		  root->approximated = true;
	    }
	}

      return;
    }

  /* Compute secular function */
  cdpe_sub_eq (pol, cdpe_one);
  rdpe_add_eq (asum, rdpe_one);

  /* Compute the module of pol */
  cdpe_mod (apol, pol);

  /* Compute newton correction */
  cdpe_mul (corr, pol, sumb);
  cdpe_add_eq (corr, fp);
  if (cdpe_eq_zero (corr))
      cdpe_set (corr, pol);
  else
    cdpe_div (corr, pol, corr);

  rdpe_div (asum_on_apol, asum, apol);

  /* If the approximation falls in the root neighbourhood then we can stop */
  rdpe_add (rtmp, asum_on_apol, rdpe_one);
  rdpe_mul_eq_d (rtmp, KAPPA * DBL_EPSILON);
  if (rdpe_gt (rtmp, rdpe_one))
    {
      if (data && s->debug_level & MPS_DEBUG_PACKETS)
	MPS_DEBUG (s, "Setting again to false on root %ld for root neighbourhood", data->k);
      root->again = false;

      rdpe_mul_d (rtmp, asum, KAPPA);
      cdpe_mod (rtmp2, fp);
      if (rdpe_lt (rtmp, rtmp2))
	{
	  if (data && s->debug_level & MPS_DEBUG_PACKETS)
	    {
	      MPS_DEBUG (s, "Setting approximated = true on root %ld for small conditioning number", data->k);
	    }
	  root->approximated = true;
	}
    }
  else 
    {
      cdpe_mod (acorr, corr);
      rdpe_mul_d (rtmp2, ax, MPS_SQRT2 * DBL_EPSILON);
      if (rdpe_lt (acorr, rtmp2))
	{
	  if (data && s->debug_level & MPS_DEBUG_PACKETS)
	    MPS_DEBUG (s, "Setting approximated to true on root %ld for small Newton correction", data->k);
	  root->again = false;  
	  root->approximated = true;
	}
    }

  if (!cdpe_eq_zero (corr) && root->again)
    {
      rdpe_t new_rad;
      cdpe_mod (new_rad, corr);
      rdpe_mul_d (rtmp, asum_on_apol, KAPPA * DBL_EPSILON);
      rdpe_add_eq (rtmp, rdpe_one);
      rdpe_mul_eq_d (rtmp, sec->n);
      rdpe_mul_eq (new_rad, rtmp);

      if (rdpe_lt (new_rad, root->drad))
	rdpe_set (root->drad, new_rad);
    }
}

/**
 * @brief Perform the evaluation of the MP Newton correction with the formula
 * obtained implicitly by the secular equation using the parallel
 * algorithm.
 * 
 * @param s The mps_status associated to the current computation
 * @param root The approximation that shall be used as evaluation point
 * @param n The length of the terms that should be summed by the function
 * @param ampc A pointer to the first MP a_i coefficient
 * @param bmpc A pointer to the first MP b_i coefficient
 * @param pol The MP value where the result of the evaluation of S(x) will
 *   be stored
 * @param fp The MP value where the result of the evaluation of S'(x) will be
 *   stored
 * @param sumb The MP value where the result of the evaluation of the sum of 
 *   the terms \f$\frac{1}{(x - b_i)}\f$ will be stored
 * @param asum The RDPE value where the sum of the module of the complex terms
 *   \f$\frac{a_i}{(x - b_i)}\f$ will be saved
 * @param asum2 The RDPE value where the sum of the module of the complex terms
 *   \f$\frac{a_i}{(x - b_i)^2}\f$ will be stored
 * @param asumb The RDPE value where the sum of the moduli of the complex terms
 *   \f$\frac{1}{x - b_i}$ will be stored
 * @return If the returned value is a positive integer it is the index of the b_i term 
 *   that is equal to x, suggesting that the alternate evaluation algorithm must be used. 
 *   If it is MPS_PARALLEL_SUM_FAILED then a floting point was encountered in the computation, 
 *   while MPS_PARALLEL_SUM_SUCCESS indicates that the evaluation was successful. 
 */
int
mps_secular_mparallel_sum (mps_status * s, mps_approximation * root, int n, mpc_t * ampc, mpc_t * bmpc,
			   mpc_t pol, mpc_t fp, mpc_t sumb, rdpe_t asum, rdpe_t asum2, rdpe_t asumb)
{
  if (n <= 4)
    {
      int i;
      mpc_t ctmp, ctmp2;
      rdpe_t rtmp;

      mpc_init2 (ctmp, s->mpwp);
      mpc_init2 (ctmp2, s->mpwp);

      for (i = 0; i < n; i++)
	{
	  /* Compute z - b_i */
	  mpc_sub (ctmp, root->mvalue, bmpc[i]);
	  
	  /* Check if we are in the case where z == b_i and return,
	   * without doing any further iteration */
	  if (mpc_eq_zero (ctmp))
	    {
	      return i;
	    }

	  /* Compute (z-b_i)^{-1} */
	  mpc_inv_eq (ctmp);
	  
	  /* Add the module of 1 / (x - b_i) to asumb */
	  mpc_rmod (rtmp, ctmp);
	  rdpe_add_eq (asumb, rtmp);

	  /* Compute sum of (z-b_i)^{-1} */
	  mpc_add_eq (sumb, ctmp);

	  /* Compute a_i / (z - b_i) */
	  mpc_mul (ctmp2, ampc[i], ctmp);

	  /* Compute the sum of module of (a_i/(z-b_i)) * (i + 2) */
	  mpc_rmod (rtmp, ctmp2);
	  rdpe_add_eq (asum, rtmp);

	  /* Add a_i / (z - b_i) to pol */
	  mpc_add_eq (pol, ctmp2);

	  /* Compute a_i / (z - b_i)^2a */
	  mpc_mul_eq (ctmp2, ctmp);
	  
	  /* Add the module of a_i / (x - b_i)^2 to asum2 */
	  mpc_rmod (rtmp, ctmp2);
	  rdpe_add_eq (asum2, rtmp);

	  /* Add it to fp */
	  mpc_sub_eq (fp, ctmp2);
	}

      mpc_clear (ctmp);
      mpc_clear (ctmp2);
      
      return MPS_PARALLEL_SUM_SUCCESS;
    }
  else 
    {
      int i = n/2, k;

      if ((k = mps_secular_mparallel_sum (s, root, i, ampc, bmpc, pol, fp, sumb, asum, asum2, asumb)) >= 0)
	{
	  return k;
	}
      if ((k = mps_secular_mparallel_sum (s, root, n-i, ampc + i, bmpc + i, pol, fp, sumb, asum, asum2, asumb)) >= 0)
	{
	  return i + k;
	}
      
      return MPS_PARALLEL_SUM_SUCCESS;
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

  rdpe_t asum, asumb, asum2;
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
  rdpe_set (asum2, rdpe_zero);

  if ((i = mps_secular_mparallel_sum (s, root, sec->n, ampc, bmpc, pol, fp, sumb, asum, asum2, asumb)) != MPS_PARALLEL_SUM_SUCCESS)
    {
      int k;
      
      mpc_set_ui (corr, 0U, 0U);

      for (k = 0; k < sec->n; k++)
	{
	  if (i != k)
	    {
	      mpc_sub (ctmp, sec->bmpc[i], sec->bmpc[k]);
	      mpc_add (ctmp2, sec->ampc[i], sec->ampc[k]);
	      mpc_div_eq (ctmp2, ctmp);
	      mpc_add_eq (corr, ctmp2);
	    }
	}

      mpc_set_ui (ctmp, 1U, 0U);
      mpc_sub_eq (corr, ctmp);

      mpc_div (corr, sec->ampc[i], corr);
      mpc_rmod (acorr, corr);
 
      rdpe_mul (rtmp, ax, s->mp_epsilon); 
      if (root->again && rdpe_lt (acorr, rtmp)) 
        {
	  root->again = false;
	  
	  /* Mark the root as approximated only if the Newton correction
	   * computation was well conditioned. */
	  rdpe_mul_d (rtmp, asum, KAPPA);
	  rdpe_mul_eq (rtmp, s->mp_epsilon);
	  if (rdpe_lt (rtmp, acorr))
	    root->approximated = true;
	}

      goto mnewton_cleanup;
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

