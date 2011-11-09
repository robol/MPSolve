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
#include <limits.h>
#include <math.h>

void
mps_secular_fnewton (mps_status * s, cplx_t x, double *rad, cplx_t corr,
                     mps_boolean * again, void * user_data)
{
  int i;
  cplx_t ctmp, ctmp2, pol, fp, sumb;
  double apol, prod_b = 1.0, new_rad;
  double asum = 0.0f;
  mps_secular_iteration_data * data = user_data;

  /* First set again to true */
  *again = true;

  mps_secular_equation *sec = (mps_secular_equation *) s->secular_equation;

  cplx_set (pol, cplx_zero);
  cplx_set (fp, cplx_zero);
  cplx_set (sumb, cplx_zero);

  for (i = 0; i < sec->n; i++)
    {
      /* Compute z - b_i */
      cplx_sub (ctmp, x, sec->bfpc[i]);

      if (cplx_eq_zero (ctmp))
          return;

      /* Computation of prod [ (z - b_i) / (z - z_j) ] */
      prod_b *= cplx_mod (ctmp);
      cplx_sub (ctmp2, x, sec->bfpc[i]);
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

  /* Compute newton correction */
  cplx_mul (corr, pol, sumb);
  cplx_add_eq (corr, fp);
  if (cplx_eq_zero (corr))
      cplx_set (corr, pol);
  else
    cplx_div (corr, pol, corr);

  /* If the approximation falls in the root neighbourhood then we can stop */
  if ((asum / apol + 1) * DBL_EPSILON > 1)
    {
      MPS_DEBUG (s, "Setting again to false on root %ld for root neighbourhood", data->k);
      *again = false;
    }

  /* If the correction is not useful in the current precision do
   * not iterate more   */
  if (*again && (cplx_mod (corr) < cplx_mod (x) * DBL_EPSILON))
    {
      MPS_DEBUG (s, "Setting again to false on root %ld for small Newton correction", data->k);
      *again = false;
    }

  /* Computation of radius with Gerschgorin */
  new_rad = (apol * s->n * prod_b * (1 + (3 * s->n + (asum/apol + 1)) * DBL_EPSILON)) + (cplx_mod (x) * DBL_EPSILON);

  /* Correct the old radius with the move that we are doing
   * and check if the new proposed radius is preferable. */
  if (new_rad < *rad || (*rad == 0) || (!data))
      *rad = new_rad;
}

void
mps_secular_dnewton (mps_status * s, cdpe_t x, rdpe_t rad, cdpe_t corr,
                     mps_boolean * again, void * user_data)
{
  int i;
  int* k = (int*) user_data;

  mps_secular_equation *sec = (mps_secular_equation *) s->secular_equation;
  mps_secular_iteration_data * data = user_data;

  cdpe_t pol, fp, sumb, ctmp, ctmp2, old_x;
  rdpe_t rtmp, rtmp2, apol, g_corr, prod_b, asum;
  mps_boolean x_is_b = false;

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
      rdpe_mul_eq_d (prod_b, 1 + 4.0 * DBL_EPSILON);
      cdpe_sub (ctmp2, x, s->droot[i]);
      cdpe_mod (rtmp2, ctmp2);
      if (!rdpe_eq_zero (rtmp2))
        {
	  rdpe_mul_eq_d (rtmp2, 1.0 - 4 * DBL_EPSILON);
	  rdpe_div_eq (prod_b, rtmp2);
        }

      /* Alternative computation if x is one of the b_i */
      if (cdpe_eq_zero (ctmp))
	  return;

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

  /* Computation of radius with Gerschgorin */
  rdpe_t new_rad;

  /* Compute the guaranteed radius */
  rdpe_set (new_rad, apol);
  rdpe_mul_eq_d (new_rad, s->n);
  rdpe_mul_eq (new_rad, prod_b);

  rdpe_set (rtmp, asum);
  rdpe_div_eq (rtmp, apol);
  rdpe_add_eq (rtmp, rdpe_one);
  rdpe_add_eq_d (rtmp, 3 * s->n);
  rdpe_mul_eq_d (rtmp, DBL_EPSILON);
  rdpe_mul_eq (new_rad, rtmp);

  /* Correct the old radius with the move that we are doing
   * and check if the new proposed radius is preferable. */
  if (data)
    data->radius_set = true;

  if (rdpe_lt (new_rad, rad) || !data)
    rdpe_set (rad, new_rad);
  else 
    data->radius_set = false;
  
  /* If newton correction is less than
   * the modules of |x| multiplied for
   * for epsilon stop */
  if (*again)
    {
      /* Computation of |x| and |corr| */
      cdpe_mod (rtmp, corr);
      cdpe_mod (rtmp2, x);
      rdpe_mul_eq_d (rtmp2, sec->n * DBL_EPSILON);

      /* If |corr| < |x| * DBL_EPSILON then stop */
      if (rdpe_lt (rtmp, rtmp2))
        {
	  MPS_DEBUG (s, "dit Setting again on root %d to false because the Newton correction is too small", data->k);
          *again = false;
        }

      rdpe_set (rtmp, asum);
      rdpe_div_eq (rtmp, apol);
      rdpe_add_eq (rtmp, rdpe_one);
      rdpe_mul_eq_d (rtmp, DBL_EPSILON);
      if (rdpe_ge (rtmp, rdpe_one))
	{
	  MPS_DEBUG (s, "dit Setting again on root %d to false because the approximation is in the root neighbourhood", data->k);
	  *again = false;
	}
    }

}

void
mps_secular_mnewton (mps_status * s, mpc_t x, rdpe_t rad, mpc_t corr,
                     mps_boolean * again, void * user_data)
{
  int i, j;
  mps_boolean x_is_b = false;
  // int * k = ((int*) user_data);
  mps_secular_iteration_data *data = (mps_secular_iteration_data*) user_data;

  /* Set again to true. If the convergence will be proved
   * during the iteration it will be set to false */
  *again = true;

  /* Get a pointer to the secular equation */
  mps_secular_equation *sec = (mps_secular_equation *) s->secular_equation;

  /* Declare temporary variables */
  mpc_t sumb, pol, fp, ctmp, ctmp2;
  cdpe_t cdtmp, cdtmp2, prod_b;
  rdpe_t rtmp, rtmp2, s_eps;

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
  cdpe_set (prod_b, cdpe_one);
  rdpe_set (s_eps, rdpe_zero);

  for (i = 0; i < sec->n; i++)
    {
      /* Compute z - b_i */
      mpc_sub (ctmp, x, sec->bmpc[i]);

      /* Keep away the case where the difference is zero */
      if (mpc_eq_zero (ctmp))
        {
          /* We are in the case where x = b_i so set j = i */
          j = i;
          mpc_set_ui (corr, 0U, 0U);

          for (i = 0; i < sec->n; i++)
            {
              if (i == j)
                continue;
              mpc_sub (ctmp, sec->bmpc[j], sec->bmpc[i]);
              mpc_add (sumb, sec->ampc[i], sec->ampc[j]);
              mpc_div_eq (sumb, ctmp);
              mpc_add_eq (corr, sumb);
            }

          mpc_sub_eq_ui (corr, 1U, 0U);
          mpc_inv_eq (corr);
          mpc_mul_eq (corr, sec->ampc[j]);

          x_is_b = true;
          break;
        }

      mpc_get_cdpe (cdtmp2, ctmp);

      cdpe_mul_eq (prod_b, cdtmp2);
      mpc_sub (ctmp2, x, s->mroot[i]);
      mpc_get_cdpe (cdtmp2, ctmp2);
      if (!cdpe_eq_zero (cdtmp2))
        {
          cdpe_div_eq (prod_b, cdtmp2);
        }

      /* Compute (z-b_i)^{-1} */
      mpc_inv_eq (ctmp);

      /* Add to the sum of (z-b_i)^{-1} */
      mpc_add_eq (sumb, ctmp);

      /* Compute a_i / (z - b_i)  */
      mpc_mul (ctmp2, sec->ampc[i], ctmp);

      /* Add a_i / (z - b_i) to pol */
      mpc_add_eq (pol, ctmp2);

      /* Compute the right epsilon for this element */
      mpc_get_cdpe (cdtmp, ctmp2);
      cdpe_mod (rtmp, cdtmp);
      if (user_data)
	{
	  rdpe_mul_eq (rtmp, s->secular_equation->dregeneration_epsilon[data->k]);
	}
      else
	{
	  rdpe_mul_eq (rtmp, s->mp_epsilon);
	  rdpe_mul_eq_d (rtmp, 4);
	}

      rdpe_add_eq (s_eps, rtmp);

      /* Compute a_i / (z - b_i)^2 */
      mpc_mul_eq (ctmp2, ctmp);

      /* Add it to fp */
      mpc_sub_eq (fp, ctmp2);
    }

  /* If x != b_i for every b_i finalize the computation */
  if (!x_is_b)
    {
      /* Subtract one from pol */
      mpc_sub_eq_ui (pol, 1, 0);
      rdpe_add_eq (s_eps, s->mp_epsilon);
      
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
    }

  /* Computation of radius with Gerschgorin */
  rdpe_t new_rad;
  cdpe_t x_cdpe;

  mpc_get_cdpe (x_cdpe, x);
  mpc_get_cdpe (cdtmp, pol);
  cdpe_mod (new_rad, cdtmp);

  /* Compute relative error on pol using new_rad, so we can
   * compute the error bound using 1 + s_eps * fl(radius) 
   */
  rdpe_div_eq (s_eps, new_rad);

  rdpe_mul_eq_d (new_rad, 1 + 4 * DBL_EPSILON);

  /* FIXME: Add the right epsilon here */
   rdpe_add (rtmp, s_eps, rdpe_one); 
   rdpe_mul_eq (new_rad, rtmp); 

  /* If the relative error is greater than one
   * follow an alternative approach, at least if
   * we have the index on which we are iterating.
   *
   * Skip this step if the integer k of the index we
   * are iterating on was not given in input, meaning
   * that we should not consider the conditioning. */
   cdpe_mod (rtmp, prod_b);
   rdpe_mul_eq_d (rtmp, 1 + 4 * DBL_EPSILON);
   rdpe_mul_eq (new_rad, rtmp);
   rdpe_mul_eq_d (new_rad, 1 + DBL_EPSILON);
   rdpe_mul_eq_d (new_rad, s->n);

   /* Add the representation error to the radius */
   mpc_get_cdpe (cdtmp, x);
   cdpe_mod (rtmp, cdtmp);
   rdpe_mul_eq (rtmp, s->mp_epsilon);
   rdpe_add_eq (new_rad, rtmp);

  /* Correct the old radius with the move that we are doing
   * and check if the new proposed radius is preferable. */
   if (data)
     data->radius_set = true;

   if (rdpe_lt (new_rad, rad) || (!data))
     rdpe_set (rad, new_rad);
   else
     data->radius_set = false;

   /* Compute guaranteed modulus of pol */
   mpc_get_cdpe (cdtmp, pol);
   cdpe_mod (rtmp, cdtmp);
   rdpe_add (rtmp2, s->mp_epsilon, rdpe_one);
   rdpe_mul_eq_d (rtmp2, 1 + s->n);
   rdpe_mul_eq (rtmp, rtmp2);

   /* Epsilon for us */
   rdpe_mul_d (rtmp2, s->mp_epsilon, 2);
   
   /* If |S(x)| < eps stop */
   if (rdpe_lt (rtmp, s->mp_epsilon))
     *again = false;

   /* Check if newton correction is less than
    * the modules of x for s->output_config->prec, and if
    * that's the case, stop. */
   if (*again)
     {
       mpc_get_cdpe (cdtmp, x);
       cdpe_mod (rtmp, cdtmp);
       rdpe_mul_eq (rtmp, s->mp_epsilon);

       mpc_get_cdpe (cdtmp, corr);
       cdpe_mod (rtmp2, cdtmp);
       
       if (rdpe_lt (rtmp2, rtmp))
	 {
	   *again = false;
	 }
     }

  /* Final cleanup */
   mpc_clear (fp);
   mpc_clear (pol);
   mpc_clear (sumb);
   mpc_clear (ctmp);
   mpc_clear (ctmp2);
}
