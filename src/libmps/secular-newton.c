/*
 * This file is part of MPSolve 3.0
 *
 * Copyright (C) 2001-2012, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors: 
 *   Dario Andrea Bini <bini@dm.unipi.it>
 *   Giuseppe Fiorentino <fiorent@dm.unipi.it>
 *   Leonardo Robol <robol@mail.dm.unipi.it>
 */


#include <mps/mps.h>
#include <limits.h>
#include <math.h>

#define MPS_2SQRT2 2.82842712474619009760
#define KAPPA (MPS_POLYNOMIAL (sec)->degree * log2(MPS_POLYNOMIAL (sec)->degree) + 7 * 1.4151135 + 1)
#define KAPPA_LINEAR (MPS_POLYNOMIAL (sec)->degree + 7 * 1.4142135623)
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
 * @param s The mps_context associated to the current computation
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
 * @return If the returned value is a positive integer it is the index of the \f$b_i\f$ term 
 *   that is equal to \f$x\f$, suggesting that the alternate evaluation algorithm must be used. 
 *   If it is MPS_PARALLEL_SUM_FAILED then a floting point was encountered in the computation, 
 *   while MPS_PARALLEL_SUM_SUCCESS indicates that the evaluation was successful. 
 */
int
mps_secular_fparallel_sum (mps_context * s, mps_approximation * root, int n, cplx_t * afpc, cplx_t * bfpc,
			   cplx_t pol, cplx_t fp, cplx_t sumb, double * asum)
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

	  /* Compute a_i / (z - b_i) */
	  cplx_mul (ctmp2, afpc[i], ctmp);

	  /* Compute the sum of module of (a_i/(z-b_i)) */
	  *asum += fabs (cplx_Re (ctmp2)) + fabs (cplx_Im (ctmp2));

	  /* Add a_i / (z - b_i) to pol */
	  cplx_add_eq (pol, ctmp2);

	  /* Compute a_i / (z - b_i)^2a */
	  cplx_mul_eq (ctmp2, ctmp);

	  /* Add it to fp */
	  cplx_sub_eq (fp, ctmp2);
	}
      
      return MPS_PARALLEL_SUM_SUCCESS;
    }
  else 
    {
      int i = n/2, k;
      if ((k = mps_secular_fparallel_sum (s, root, i, afpc, bfpc, pol, fp, sumb, asum)) >= 0)
	{
	  return k;
	}
      if ((k = mps_secular_fparallel_sum (s, root, n-i, afpc + i, bfpc + i, pol, fp, sumb, asum)) >= 0)
	{
	  return i + k;
	}
      
      return MPS_PARALLEL_SUM_SUCCESS;
    }
}

void
mps_secular_fnewton (mps_context * s, mps_polynomial * p, mps_approximation * root, cplx_t corr)
{
  int i;
  cplx_t ctmp, ctmp2, pol, fp, sumb;
  double apol, acorr;
  double asum = 0.0, asum_on_apol, ax = cplx_mod (root->fvalue);
  mps_secular_equation *sec = MPS_SECULAR_EQUATION (p);

  cplx_t *afpc, *bfpc;

  afpc = sec->afpc;
  bfpc = sec->bfpc;

  cplx_t x;
  cplx_set (x, root->fvalue);

  /* First set again to true */
  root->again = true;

  cplx_set (pol, cplx_zero);
  cplx_set (fp, cplx_zero);
  cplx_set (sumb, cplx_zero);
  cplx_set (corr, cplx_zero);

  if ((i = mps_secular_fparallel_sum (s, root, MPS_POLYNOMIAL (sec)->degree, sec->afpc, 
				      sec->bfpc, pol, 
				      fp, sumb, &asum)) >= 0)
    {
      int k;
      asum = 0.0;

      cplx_set (corr, cplx_zero);

      for (k = 0; k < MPS_POLYNOMIAL (sec)->degree; k++)
  	    {
  	       if (i != k)
  	        {
  	          cplx_sub (ctmp, bfpc[i], bfpc[k]);
  	          cplx_add (ctmp2, afpc[i], afpc[k]);
  	          cplx_inv_eq (ctmp);
              cplx_mul_eq (ctmp2, ctmp);
  	          cplx_add_eq (corr, ctmp2);

              asum += fabs (cplx_Re (ctmp2)) + fabs (cplx_Im (ctmp2));
           }
  	    }

      if (i == MPS_PARALLEL_SUM_FAILED)
        {
          root->status = MPS_ROOT_STATUS_NOT_FLOAT;
          root->again = false;
          return;
        }

      cplx_sub_eq (corr, cplx_one);

      if (!cplx_eq_zero (corr))
      	{
      	  cplx_div (corr, afpc[i], corr);
      	      
      	  acorr = cplx_mod (corr);
      	  if (acorr < ax * DBL_EPSILON)
      	    {
      	      root->again = false;
      	    }
      	}
      else
	      root->again = false;
      
      return;
    }

  if (i == MPS_PARALLEL_SUM_FAILED)
    {
      root->status = MPS_ROOT_STATUS_NOT_FLOAT;
      root->again = false;
      return;
    }

  /* Compute secular function */
  cplx_sub_eq (pol, cplx_one);
  asum += 1.0;

  /* Compute the module of pol */
  apol = cplx_mod (pol);
  /* afp = cplx_mod (fp); */

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
      if (s->debug_level & MPS_DEBUG_PACKETS)
	      MPS_DEBUG (s, "Setting again to false on root for root neighbourhood");
      root->again = false;
    }
  else if (acorr < MPS_SQRT2 * ax * DBL_EPSILON)
    {
      if (s->debug_level & MPS_DEBUG_PACKETS)
	      MPS_DEBUG (s, "Setting approximated to true on root for small Newton correction");
      root->again = false;  
      root->approximated = true;
    }

  if (!cplx_eq_zero (corr) && root->again)
    {
      double new_rad = acorr * s->n * (1 + KAPPA * DBL_EPSILON * asum_on_apol) +
        ax * DBL_EPSILON * 4.0;

      if ((new_rad > 0) && (new_rad < root->frad))
	      root->frad = new_rad;
    }
}

/**
 * @brief Perform the evaluation of the DPE Newton correction with the formula
 * obtained implicitly by the secular equation using the parallel
 * algorithm.
 * 
 * @param s The mps_context associated to the current computation
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
 * @return If the returned value is a positive integer it is the index of the \f$b_i\f$ term 
 *   that is equal to \f$x\f$, suggesting that the alternate evaluation algorithm must be used. 
 *   If it is MPS_PARALLEL_SUM_FAILED then a floting point was encountered in the computation, 
 *   while MPS_PARALLEL_SUM_SUCCESS indicates that the evaluation was successful. 
 */
int
mps_secular_dparallel_sum (mps_context * s, mps_approximation * root, int n, cdpe_t * adpc, cdpe_t * bdpc,
			   cdpe_t pol, cdpe_t fp, cdpe_t sumb, rdpe_t asum)
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
	    
	  /* Compute (z-b_i)^{-1} */
	  cdpe_inv_eq (ctmp);

	  /* Compute sum of (z-b_i)^{-1} */
	  cdpe_add_eq (sumb, ctmp);

	  /* Compute a_i / (z - b_i) */
	  cdpe_mul (ctmp2, adpc[i], ctmp);

	  /* Compute the sum of module of (a_i/(z-b_i)) * (i + 2) */
          rdpe_abs (rtmp, cdpe_Re (ctmp2));
	  rdpe_add_eq (asum, rtmp);
          rdpe_abs (rtmp, cdpe_Im (ctmp2));
          rdpe_add_eq (asum, rtmp);

	  /* Add a_i / (z - b_i) to pol */
	  cdpe_add_eq (pol, ctmp2);

	  /* Compute a_i / (z - b_i)^2a */
	  cdpe_mul_eq (ctmp2, ctmp);
	  
	  /* Add it to fp */
	  cdpe_sub_eq (fp, ctmp2);
	}
      
      return MPS_PARALLEL_SUM_SUCCESS;
    }
  else 
    {
      int i = n/2, k;
      if ((k = mps_secular_dparallel_sum (s, root, i, adpc, bdpc, pol, fp, sumb, asum)) >= 0)
	{
	  return k;
	}
      if ((k = mps_secular_dparallel_sum (s, root, n-i, adpc + i, bdpc + i, pol, fp, sumb, asum)) >= 0)
	{
	  return i + k;
	}
      
      return MPS_PARALLEL_SUM_SUCCESS;
    }
}

void
mps_secular_dnewton (mps_context * s, mps_polynomial * p, mps_approximation * root, cdpe_t corr)
{
  int i;
  cdpe_t ctmp, ctmp2, pol, fp, sumb, x;
  rdpe_t apol, asum, asum_on_apol, ax, rtmp, rtmp2, acorr;
  mps_secular_equation *sec = MPS_SECULAR_EQUATION (p);

  cdpe_set (x, root->dvalue);
  rdpe_set (asum, rdpe_zero);
  cdpe_mod (ax, x);

  /* First set again to true */
  root->again = true;

  cdpe_set (pol, cdpe_zero);
  cdpe_set (fp, cdpe_zero);
  cdpe_set (sumb, cdpe_zero);
  cdpe_set (corr, cdpe_zero);

  if ((i = mps_secular_dparallel_sum (s, root, MPS_POLYNOMIAL (sec)->degree, sec->adpc, sec->bdpc, 
				       pol, fp, sumb, asum)) != MPS_PARALLEL_SUM_SUCCESS)
    {
      int k;

      cdpe_set (corr, cdpe_zero);

      for (k = 0; k < MPS_POLYNOMIAL (sec)->degree; k++)
      	{
      	  if (i != k)
      	    {
      	      cdpe_sub (ctmp, sec->bdpc[i], sec->bdpc[k]);
      	      cdpe_add (ctmp2, sec->adpc[i], sec->adpc[k]);
              cdpe_inv_eq (ctmp);
      	      cdpe_mul_eq (ctmp2, ctmp);
      	      cdpe_add_eq (corr, ctmp2);
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
	    }
	}
      else
	root->again = false;

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
      if (s->debug_level & MPS_DEBUG_PACKETS)
	      MPS_DEBUG (s, "Setting again to false on root for root neighbourhood");
      root->again = false;
    }
  else 
    {
      cdpe_mod (acorr, corr);
      rdpe_mul_d (rtmp2, ax, MPS_SQRT2 * DBL_EPSILON);
      if (rdpe_lt (acorr, rtmp2))
	{
	  if (s->debug_level & MPS_DEBUG_PACKETS)
	    MPS_DEBUG (s, "Setting approximated to true on root for small Newton correction");
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
      rdpe_mul_eq_d (rtmp, MPS_POLYNOMIAL (sec)->degree);
      rdpe_mul_eq (new_rad, rtmp);

      if (rdpe_lt (new_rad, root->drad))
	      rdpe_set (root->drad, new_rad);
    }
}

/**
 * @brief Perform the evaluation of the Newton correction with the formula
 * obtained implicitly by the secular equation using the parallel
 * algorithm.
 * 
 * @param s The mps_context associated to the current computation
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
 * @return If the returned value is a positive integer it is the index of the \f$b_i\f$ term 
 *   that is equal to \f$x\f$, suggesting that the alternate evaluation algorithm must be used. 
 *   If it is MPS_PARALLEL_SUM_FAILED then a floting point was encountered in the computation, 
 *   while MPS_PARALLEL_SUM_SUCCESS indicates that the evaluation was successful. 
 */
int
mps_secular_mparallel_sum (mps_context * s, mps_approximation * root, int n, mpc_t * ampc, mpc_t * bmpc,
			   mpc_t pol, mpc_t fp, mpc_t sumb, rdpe_t asum)
{
  long int wp = mpc_get_prec (pol);
  
  if (n <= 2)
    {
      int i;
      mpc_t ctmp, ctmp2;

      mpc_init2 (ctmp, wp);
      mpc_init2 (ctmp2, wp);

      for (i = 0; i < n; i++)
	{
	  rdpe_t rtmp;

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

	  /* Compute sum of (z-b_i)^{-1} */
	  mpc_add_eq (sumb, ctmp);

	  /* Compute a_i / (z - b_i) */
	  mpc_mul (ctmp2, ampc[i], ctmp);

	  /* Compute the sum of module of (a_i/(z-b_i)) */
	  mpc_rmod (rtmp, ctmp2);
	  rdpe_add_eq (asum, rtmp);

	  /* Add a_i / (z - b_i) to pol */
	  mpc_add_eq (pol, ctmp2);

	  /* Compute a_i / (z - b_i)^2a */
	  mpc_mul_eq (ctmp2, ctmp);

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
      if ((k = mps_secular_mparallel_sum (s, root, i, ampc, bmpc, pol, fp, sumb, asum)) >= 0)
	{
	  return k;
	}
      if ((k = mps_secular_mparallel_sum (s, root, n-i, ampc + i, bmpc + i, pol, fp, sumb, asum)) >= 0)
	{
	  return i + k;
	}
      
      return MPS_PARALLEL_SUM_SUCCESS;
    }
}

void
mps_secular_mnewton (mps_context * s, mps_polynomial * p, mps_approximation * root, mpc_t corr)
{
  int i;
  long int wp = mpc_get_prec (root->mvalue);
  mpc_t ctmp, ctmp2, pol, fp, sumb, x;
  rdpe_t apol, acorr, rtmp, epsilon;
  rdpe_t asum, asum_on_apol, ax, axeps;
  mps_secular_equation *sec = MPS_SECULAR_EQUATION (p);

  /* Init MP variables */
  mpc_init2 (x, wp);
  mpc_init2 (ctmp, wp);
  mpc_init2 (ctmp2, wp);
  mpc_init2 (pol, wp);
  mpc_init2 (fp, wp);
  mpc_init2 (sumb, wp);

  mpc_set (x, root->mvalue);

  rdpe_set (asum, rdpe_zero);
  mpc_rmod (ax, root->mvalue);

  /* Setup a reasonable epsilon to use for the checks such as |corr| < |x| * eps */
  rdpe_set_2dl (epsilon, 1.0, 1 - wp);
  rdpe_mul (axeps, ax, epsilon);
  rdpe_mul_eq_d (axeps, 4.0);

  mpc_t *ampc, *bmpc;

  ampc = sec->ampc;
  bmpc = sec->bmpc;

  /* First set again to true */
  root->again = true;

  mpc_set_ui (pol, 0U, 0U);
  mpc_set_ui (fp, 0U, 0U);
  mpc_set_ui (sumb, 0U, 0U);
  mpc_set_ui (corr, 0U, 0U);

  if ((i = mps_secular_mparallel_sum (s, root, MPS_POLYNOMIAL (sec)->degree, sec->ampc, 
				      sec->bmpc, pol, 
				      fp, sumb, asum)) >= 0)
    {
      int k;

      rdpe_set (asum, rdpe_zero);
      mpc_set_ui (corr, 0U, 0U);

      for (k = 0; k < MPS_POLYNOMIAL (sec)->degree; k++)
  	    {
  	       if (i != k)
  	        {
		  rdpe_t rtmp;

  	          mpc_sub (ctmp, bmpc[i], bmpc[k]);
  	          mpc_add (ctmp2, ampc[i], ampc[k]);
  	          mpc_inv_eq (ctmp);
		  mpc_mul_eq (ctmp2, ctmp);
  	          mpc_add_eq (corr, ctmp2);

		  mpc_rmod (rtmp, ctmp2);
		  rdpe_add_eq (asum, rtmp);
		}
  	    }

      mpc_sub_eq_ui (corr, 1U, 0U);

      if (!mpc_eq_zero (corr))
      	{
      	  mpc_div (corr, ampc[i], corr);
	  mpc_rmod (acorr, corr);

	  if (rdpe_lt (acorr, axeps))
      	    {
      	      root->again = false;
      	    }
      	}
      else
	root->again = false;
      
      goto mnewton_cleanup;
    }

  /* Compute secular function */
  mpc_sub_eq_ui (pol, 1U, 0U);
  rdpe_add_eq (asum, rdpe_one);

  /* Compute the module of pol */
  mpc_rmod (apol, pol);

  /* Compute newton correction */
  mpc_mul (corr, pol, sumb);
  mpc_add_eq (corr, fp);
  if (mpc_eq_zero (corr))
    {
      mpc_set (corr, pol);
      root->again = false;
      goto mnewton_cleanup;
    }
  else
    mpc_div (corr, pol, corr);

  rdpe_div (asum_on_apol, asum, apol);
  mpc_rmod (acorr, corr);

  /* If the approximation falls in the root neighbourhood then we can stop */
  /* In floating point, (asum_on_apol + 1) * KAPPA * DBL_EPSILON > 1 */
  rdpe_add (rtmp, asum_on_apol, rdpe_one);
  rdpe_mul_eq_d (rtmp, KAPPA);
  rdpe_mul_eq (rtmp, epsilon);
  if (rdpe_gt (rtmp, rdpe_one))
    {
      if (s->debug_level & MPS_DEBUG_PACKETS)
	      MPS_DEBUG (s, "Setting again to false on root for root neighbourhood");
      root->again = false;
    }
  else if (rdpe_lt (acorr, axeps))
    {
      if (s->debug_level & MPS_DEBUG_PACKETS)
	      MPS_DEBUG (s, "Setting approximated to true on root for small Newton correction");
      root->again = false;  
      root->approximated = true;
    }
  
  if (!mpc_eq_zero (corr) && root->again)
    {
      rdpe_t new_rad;
      rdpe_t rad_epsilon;

      rdpe_mul_d (new_rad, acorr, s->n);
      rdpe_mul (rad_epsilon, epsilon, asum_on_apol);
      rdpe_mul_eq_d (rad_epsilon, KAPPA);
      rdpe_add_eq (rad_epsilon, rdpe_one);
      rdpe_mul_eq (new_rad, rad_epsilon);

      rdpe_add_eq (new_rad, axeps);

      if (rdpe_lt (new_rad, root->drad))
	      rdpe_set (root->drad, new_rad);
    }

 mnewton_cleanup:
  mpc_clear (ctmp);
  mpc_clear (ctmp2);
  mpc_clear (pol);
  mpc_clear (fp);
  mpc_clear (sumb);
  mpc_clear (x);
}
