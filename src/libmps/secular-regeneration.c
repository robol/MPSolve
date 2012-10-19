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
#include <math.h>
#include <string.h>

/**
 * @brief Update the working precision of a root, i.e. the variable
 * <code>s->root[i]->wp</code> with the given precision rounded to the 
 * closer multiple of 64 (since GMP does not handle intermediate precisions).
 *
 * @param s The <code>mps_context</code> of the computation.
 * @param i The index of the root whose precision must be updated.
 * @param wp The precision to set.
 */
long int
mps_secular_ga_update_root_wp (mps_context * s, int i, long int wp, mpc_t * bmpc)
{
  mps_secular_equation * sec = s->secular_equation;  
  mps_polynomial * p = s->active_poly;
  int j;

  s->root[i]->wp = ((wp - 1) / 64 + 1) * 64;
  
  MPS_LOCK (s->data_prec_max);
  if (s->data_prec_max.value < s->root[i]->wp)
    s->data_prec_max.value = s->root[i]->wp;
  MPS_UNLOCK (s->data_prec_max);

  if (s->debug_level & MPS_DEBUG_MEMORY)  
    MPS_DEBUG (s, "Setting wp for root %d to %ld bits", i, s->root[i]->wp);  

  if (sec->bmpc == bmpc)
    {
      for (j = 0; j < s->n; j++)
	{
	  pthread_mutex_lock (&sec->bmpc_mutex[j]); 
	  if (mpc_get_prec (bmpc[j]) < s->root[i]->wp) 
	    mpc_set_prec (bmpc[j], s->root[i]->wp);
	  pthread_mutex_unlock (&sec->bmpc_mutex[j]); 
	}
    }
  else
    {
      for (j = 0; j < s->n; j++)
	{
	  if (mpc_get_prec (bmpc[j]) < s->root[i]->wp) 
	    mpc_set_prec (bmpc[j], s->root[i]->wp);
	}
    }
  
  pthread_mutex_lock (&sec->ampc_mutex[i]);
  if (mpc_get_prec (sec->ampc[i]) < s->root[i]->wp)  
    mpc_set_prec (sec->ampc[i], s->root[i]->wp);  
  pthread_mutex_unlock (&sec->ampc_mutex[i]);

  mps_polynomial_raise_data (s, p, s->root[i]->wp);
  
  return s->root[i]->wp;
}

/**
 * @brief This routines is used to check if a root has changed from the last regeneration, 
 * in floating point phases. 
 *
 * If a root is approximated or isolated and does not differ much (i.e. less than the machine
 * epsilon) from the approximation that was present a cycle ago, than it's not necessary to 
 * recompute the value of the polynomial in that point, so <code>root_changed[i]</code> is set
 * to false
 *
 * @param s The <code>mps_context</code> of the computation.
 * @param old_b A vector of the old \f$b_i\f$ coefficients.
 * @param old_mb The vector of the old \f$b_i\f$ in <code>mpc_t</code> version. Must be set
 * to NULL if we are not in a multiprecision phase. 
 */
mps_boolean *
mps_secular_ga_find_changed_roots (mps_context * s, cdpe_t * old_b, mpc_t * old_mb)
{
  MPS_DEBUG_THIS_CALL;

  cdpe_t diff;
  int i;
  int changed_roots = 0;

  mps_secular_equation * sec = s->secular_equation;
  mps_boolean * root_changed = mps_boolean_valloc (s->n);

  mpc_t mdiff;
  if (old_mb)
    mpc_init2 (mdiff, mpc_get_prec (old_mb[0]));
  
  for (i = 0; i < s->n; i++)
    {
      if (s->just_raised_precision && MPS_INPUT_CONFIG_IS_MONOMIAL (s->input_config))
	{
	  root_changed[i] = true;
	  continue;
	}

      /* Do multiprecision sub if we are in mp_phase, otherwise go with plain
       * CDPE. */
      if (old_mb)
	{
	  mpc_sub (mdiff, old_mb[i], sec->bmpc[i]);
	  mpc_get_cdpe (diff, mdiff);
	  mpc_get_cdpe (sec->bdpc[i], sec->bmpc[i]);
	}
      else
	cdpe_sub (diff, old_b[i], sec->bdpc[i]);

      root_changed[i] = ! cdpe_eq_zero (diff);
      if ((s->debug_level & MPS_DEBUG_REGENERATION) && !root_changed[i])
	MPS_DEBUG (s, "b_%d hasn't changed, so p(b_%d) will not be recomputed", i, i);

      if (root_changed[i])
	changed_roots++;
    }

  if (old_mb)
    mpc_clear (mdiff);

  if (changed_roots != 0)
    MPS_DEBUG (s, "%d of %d approximations are different from last regeneration", changed_roots, s->n);

  return root_changed;
}

struct __mps_secular_ga_regenerate_coefficients_monomial_data {
  mps_context * s;
  cdpe_t * old_b;
  mpc_t * old_mb;
  mpc_t * bmpc;
  mps_boolean * root_changed;
  rdpe_t * root_epsilon;
  mps_boolean * success;
  int i;
};

void *
__mps_secular_ga_regenerate_coefficients_monomial_worker (void * data_ptr)
{
  struct __mps_secular_ga_regenerate_coefficients_monomial_data * data = data_ptr;

  mps_context * s = data->s;
  mpc_t * old_mb = data->old_mb;
  mpc_t * bmpc = data->bmpc + (mps_thread_get_id (s, s->pool) * s->n);
  mps_boolean * root_changed = data->root_changed;
  
  /* Pointers to the secular equation and the monomial_poly */
  mps_secular_equation * sec = s->secular_equation;
  mps_polynomial * p = s->active_poly;
  int i = data->i, j;

  /* Multiprecision variables */
  mpc_t mprod_b, ctmp, mdiff, lc, my_b;

  /* Floating point temporary variables */
  rdpe_t root_epsilon;

  /* The precision of the temporary variables at the start of the computation. We can set
   * this to s->mpwp; */
  long int coeff_wp = s->root[i]->wp;

  /* This variable is true if the regeneration succeeded. */
  mps_boolean success = true;

  /* mps_secular_raise_coefficient_precision (s, coeff_wp); */

  switch (s->lastphase)
    {
      /* If we are in floating point then the roots are know to
       * DBL_EPSILON precision */
    case float_phase:
    case dpe_phase:
      rdpe_set_d (root_epsilon, DBL_EPSILON);
      break;
      /* But if we are in multiprecision we can check the epsilon by looking
       * at s->mp_epsilon */
    case mp_phase:
      rdpe_set (root_epsilon, s->mp_epsilon);
      break;
    default:
      break;
    }

  /* Init multiprecision values */
  mpc_init2 (mprod_b, coeff_wp);
  mpc_init2 (ctmp, coeff_wp);
  mpc_init2 (mdiff, coeff_wp);
  mpc_init2 (lc, coeff_wp);
  mpc_init2 (my_b, coeff_wp);

  mpc_set_si (lc, -1, 0);
  if (MPS_INPUT_CONFIG_IS_MONOMIAL (s->input_config))
    {
      mps_monomial_poly *mp = MPS_MONOMIAL_POLY (p);
      mps_with_lock (mp->mfpc_mutex[s->n],
		     mpc_div_eq (lc, mp->mfpc[s->n]);
		     );
    }
      
  /*
   * The new coefficients of the secular equation can be computed
   * starting from the evaluation of the polynomial changed of
   * sign and divided for the product of the difference of the
   * b_i, i.e.:
   * 
   *   a_i = -p(b_i) / \prod_{i \neq j} (b_i - b_j)
   *
   */
  if (root_changed[i])
    {
      rdpe_t relative_error, rtmp;
      cdpe_t cpol;
      mpc_t tx;

      /* Set up a temporary memory location to hold the value of b_i, since we need
       * to play with its precision. This is not doable directly because it will
       * disturb other threads at work. */
      mpc_init2 (tx, s->root[i]->wp);

      mps_secular_ga_update_root_wp (s, i, s->root[i]->wp, bmpc);

      /* pthread_mutex_lock (&sec->bmpc_mutex[i]); */
      /* mpc_set_prec (sec->bmpc[i], s->root[i]->wp); */
      mpc_set (tx, sec->bmpc[i]);
      mps_polynomial_meval (s, p, tx, sec->ampc[i], relative_error);
      /* pthread_mutex_unlock (&sec->bmpc_mutex[i]); */

      if (s->debug_level & MPS_DEBUG_REGENERATION)
	{
	  MPS_DEBUG_MPC  (s, mpc_get_prec (sec->ampc[i]), sec->ampc[i], "p(b_%d)", i);
	  MPS_DEBUG_RDPE (s, relative_error, "Absolute error on p(b_%d) evaluation", i);
	}

      mpc_get_cdpe (cpol, sec->ampc[i]);
      cdpe_mod (rtmp, cpol);
      rdpe_div_eq (relative_error, rtmp);

      if (s->debug_level & MPS_DEBUG_REGENERATION)
	MPS_DEBUG_RDPE (s, relative_error, "Relative_error on p(b_%d) evaluation", i);

      /* mpc_get_cdpe (cpol, sec->bmpc[i]);  */
      /* cdpe_mod (rtmp, cpol);  */
      /* rdpe_mul_eq (rtmp, root_epsilon);  */

      while (rdpe_gt (relative_error, root_epsilon) && (s->root[i]->wp < s->n * s->mpwp))
	{
	  /* Update the working precision of the selected root with a realistic estimate of the
	   * required precision to get a result exact to machine precision */
	  mps_secular_ga_update_root_wp (s, i, 1 + s->root[i]->wp + (rdpe_Esp (relative_error) - rdpe_Esp (root_epsilon)), bmpc);
	  
	  /* Try to recompute the polynomial with the augmented precision and see if now relative_error matches */
	  /* pthread_mutex_lock (&sec->bmpc_mutex[i]); */
	  mpc_set_prec (tx, s->root[i]->wp);
	  mps_polynomial_meval (s, p, tx, sec->ampc[i], relative_error);
	  /* pthread_mutex_unlock (&sec->bmpc_mutex[i]); */
	  
	  if (s->debug_level & MPS_DEBUG_REGENERATION)
	    {
	      MPS_DEBUG_MPC  (s, mpc_get_prec (sec->ampc[i]), sec->ampc[i], "p(b_%d)", i);
	      MPS_DEBUG_RDPE (s, relative_error, "Absolute error on p(b_%d) evaluation", i);
	      MPS_DEBUG_MPC (s, mpc_get_prec (sec->ampc[i]), sec->ampc[i], "p(b_%d)", i);
	    }

	  mpc_get_cdpe (cpol, sec->ampc[i]);
	  cdpe_mod (rtmp, cpol);
	  rdpe_div_eq (relative_error, rtmp);
	      
	  if (s->debug_level & MPS_DEBUG_REGENERATION)
	    {
	      MPS_DEBUG_RDPE (s, relative_error, "Relative_error on p(b_%d) evaluation", i);
	    }
	}

      if (mpc_get_prec (mdiff) < s->root[i]->wp)
	{
	  mpc_set_prec (mdiff, s->root[i]->wp);
	  mpc_set_prec (mprod_b, s->root[i]->wp);
	  mpc_set_prec (lc, s->root[i]->wp);
	  mpc_set_si (lc, -1, 0);
	  if (MPS_INPUT_CONFIG_IS_MONOMIAL (s->input_config))
	    {
	      mps_monomial_poly *mp = MPS_MONOMIAL_POLY (p);
	      mps_with_lock (mp->mfpc_mutex[s->n],
			     mpc_div_eq (lc, mp->mfpc[s->n]);
			     );
	    }
	  mpc_set_prec (my_b, s->root[i]->wp);
	}

      pthread_mutex_lock (&sec->bmpc_mutex[i]);
      mpc_set (my_b, sec->bmpc[i]);
      pthread_mutex_unlock (&sec->bmpc_mutex[i]);
      
      /* Compute the difference of the b_i */
      mpc_set_ui (mprod_b, 1U, 0U);

      for (j = 0; j < s->n; ++j)
	    {
	      if (i == j)
		continue;
	      
	      /* Lock i-th and -j-th mutex to gain control of the the
	       * b_i of the secular equation. */
	      /* pthread_mutex_lock (&sec->bmpc_mutex[j]); */
	      mpc_sub (mdiff, my_b, bmpc[j]);
	      /* pthread_mutex_unlock (&sec->bmpc_mutex[j]);  */
	      
	      /* If the difference is zero then regeneration cannot succeed, and means
	       * that we need more precision in the roots */
	      if (mpc_eq_zero (mdiff))
		{
		  MPS_DEBUG (s, "Regeneration of the coefficients failed because sec->bdpc[%d] == sec->bdpc[%d]", i, j);
		  success = false;
		  goto monomial_regenerate_exit;
		}
	      mpc_mul_eq (mprod_b, mdiff);
	    }
      
      /* Actually divide the result and store it in
       * a_i, as requested. */
      mpc_div_eq (sec->ampc[i], mprod_b);
      mpc_mul_eq (sec->ampc[i], lc);
      
      /* Debug computed coefficients */
      if (s->debug_level & MPS_DEBUG_REGENERATION)
	{
	  MPS_DEBUG_MPC (s, s->mpwp, sec->ampc[i], "a_%d", i);
	  MPS_DEBUG_MPC (s, s->mpwp, sec->bmpc[i], "b_%d", i);
	}

      mpc_clear (tx);
      
    } /* Close the case where the coefficient are not approximated or isolated */
      else
	{
	  mpc_set_ui (mprod_b, 1U, 0U);
	  for (j = 0; j < MPS_POLYNOMIAL (sec)->degree; j++) {
	    if (root_changed[j] && i != j) {
	      mpc_sub (mdiff, bmpc[i], old_mb[j]); 
	      mpc_mul_eq (mprod_b, mdiff); 
	      mpc_sub (mdiff, bmpc[i], bmpc[j]);
	      mpc_div_eq (mprod_b, mdiff);
	    }
	  }
	  
	  mpc_mul_eq (sec->ampc[i], mprod_b);
	}

 monomial_regenerate_exit:
  /* Clear requested storage */
  mpc_clear (mdiff);
  mpc_clear (mprod_b);
  mpc_clear (ctmp);
  mpc_clear (lc);
  mpc_clear (my_b);
  /* mps_boolean_vfree (root_changed); */

  if (!success)
    *data->success = false;
  return NULL;
}

/**
 * @brief Compute the new secular equation coefficients based on the monomial input
 * in <code>s->monomial_poly</code>.
 *
 * @param s The <code>mps_context</code> of the computation
 * @param old_b The old \f$b_i\f$ coefficients of the secular equation used to recompute
 * the value of \f$a_i\f$ in the case where, denoting with \f$b_i\f$ the new coefficients, 
 * with \f$\tilde b_i\f$ the old ones and with \f$u\f$ the machine precision:
 * \f[|b_i - \tilde b_i| < u\f]
 * so there is no need to recompute the value of \f$p(b_i)\f$.
 * @param old_mb The MP version of <code>old_b</code>, or NULL if we are not in MP. 
 * @param root_changed A vector of booleans that is <codefalse</code> on the components that
 * did not changed from the last regeneration. 
 */
mps_boolean
mps_secular_ga_regenerate_coefficients_monomial (mps_context * s, cdpe_t * old_b, mpc_t * old_mb, mps_boolean * root_changed)
{
  MPS_DEBUG_THIS_CALL;

  int i, j;
  mps_secular_equation * sec = s->secular_equation;
  long current_wp = mpc_get_prec (sec->bmpc[0]);
  mps_boolean success = true;

  struct __mps_secular_ga_regenerate_coefficients_monomial_data * data = 
    mps_newv (struct __mps_secular_ga_regenerate_coefficients_monomial_data, s->n);

  MPS_DEBUG (s, "Regenerating coefficients from monomial input");

  if (!s->bmpc)
    {
      if (s->debug_level & MPS_DEBUG_MEMORY)
	MPS_DEBUG (s, "Allocating space for thread-local bmpc coefficients");
      s->bmpc = mps_newv (mpc_t, s->n * s->pool->n);
      mpc_vinit2 (s->bmpc, s->n * s->pool->n, s->mpwp);
    }

  for (i = 0; i < s->pool->n; i++)
    {
      for (j = 0; j < s->n; j++)
	{
	  mpc_set_prec (s->bmpc[i * s->n + j], current_wp);
	  mpc_set (s->bmpc[i * s->n + j], sec->bmpc[j]);
	}
    }

  for (i = s->n - 1; i >= 0; i--)
    {
      data[i].i = i;
      data[i].old_b = old_b;
      data[i].old_mb = old_mb;
      data[i].root_changed = root_changed;
      data[i].s = s;
      data[i].success = &success;
      data[i].bmpc = s->bmpc;
      mps_thread_pool_assign (s, s->pool, __mps_secular_ga_regenerate_coefficients_monomial_worker,   
       			      data + i);
      /* __mps_secular_ga_regenerate_coefficients_monomial_worker (data + i);  */
    }

  mps_thread_pool_wait (s, s->pool);

  /* mpc_vclear (bmpc, s->n * s->pool->n); */
  /* free (bmpc); */

  free (data);

  return success;
}

/**
 * @brief Regenerate the coefficients \f$a_i\f$ based on
 * the \f$b_i\f$.
 *
 * The coefficients used are the one in multiprecision, so
 * if floating point coefficient regeneration is desired, one
 * must take care to copy the floating point coefficient
 * in the multiprecision before calling this function, and
 * to copy them back after the computation.
 *
 * A useful function that does all the things above, doing
 * the rigth thing based on <code>s->lastphase</code>, where
 * <code>s</code> is a pointer to the <code>mps_context</code>
 * is <code>mps_secular_ga_regenerate_coefficients()</code> and
 * is the one that should be used
 *
 * @param s The mps_context of the computation.
 * @param old_b Old \f$b_i\f$ coefficients from which we are regenerating from. These are used
 * to adjust the a_i for already approximated roots.
 * @param old_mb Old \f$b_i\f$ in the multiprecision version. These are used only in the case
 * where <code>s->lastphase == mp_phase</code> to determine if the roots have changed. Otherwise
 * it must be set to NULL.
 */
mps_boolean
mps_secular_ga_regenerate_coefficients_mp (mps_context * s, cdpe_t * old_b, mpc_t * old_mb)
{
  /* Declaration and initialization of the multprecision
   * variables that are used only in that case */
  mps_boolean success = true;

  /* Get the list of the changed roots so we can compute the value of the 
   * polynomial only in that approximations. */
  mps_boolean * root_changed = mps_secular_ga_find_changed_roots (s, old_b, old_mb);

  /* Regenerate the coefficients of the secular equation starting from the monomial input */
  success = mps_secular_ga_regenerate_coefficients_monomial (s, old_b, old_mb, root_changed);

  if (!success)
      MPS_DEBUG (s, "Regeneration of the coefficients failed");

  mps_boolean_vfree (root_changed);

  return success;
}

/**
 * @brief Regenerate \f$a_i\f$ and \f$b_i\f$ setting
 * \f$b_i = z_i\f$, i.e. the current root approximation
 * and recomputing \f$a_i\f$ accordingly.
 *
 * @param s The mps_context of the computation.
 */
mps_boolean
mps_secular_ga_regenerate_coefficients (mps_context * s)
{
  MPS_DEBUG_THIS_CALL;

  s->operation = MPS_OPERATION_REGENERATION;

  cplx_t *old_b, *old_a;
  cdpe_t *old_db, *old_da;
  mpc_t *old_ma, *old_mb;
  mps_secular_equation *sec;
  int i;
  mps_boolean successful_regeneration = true;

  sec = (mps_secular_equation *) s->secular_equation;

  MPS_DEBUG_WITH_INFO (s, "Regenerating coefficients");

  old_mb = mpc_valloc (s->n);
  for (i = 0; i < s->n; i++)
    mpc_init2 (old_mb[i], s->root[i]->wp);

  /* Start timer and add execution time to the total counter */
#ifndef DISABLE_DEBUG
  clock_t *my_clock = mps_start_timer ();
#endif

  switch (s->lastphase)
    {
      /* If we are in the float phase regenerate coefficients
       * starting from floating point */
    case float_phase:

      s->mpwp = MPS_SECULAR_EQUIVALENT_FP_PRECISION;

      /* Allocate old_a and old_b */
      old_a = cplx_valloc (s->n);
      old_b = cplx_valloc (s->n);
      old_db = cdpe_valloc (s->n);

      /* Copy the old coefficients, and set the new
       * b_i with the current roots approximations. */
      for (i = 0; i < s->n; i++)
        {
          cplx_set (old_a[i], sec->afpc[i]);
          cplx_set (old_b[i], sec->bfpc[i]);
	  cdpe_set_x (old_db[i], old_b[i]);
	  mpc_set_cplx (old_mb[i], old_b[i]);
	  mpc_set_cplx (sec->bmpc[i], s->root[i]->fvalue);
        }

      mps_secular_ga_update_coefficients (s);

      /* Regeneration */
      if (!(successful_regeneration = mps_secular_ga_regenerate_coefficients_mp (s, old_db, old_mb)))
        {
          for (i = 0; i < s->n; i++)
            {
              cplx_set (sec->afpc[i], old_a[i]);
              cplx_set (sec->bfpc[i], old_b[i]);
            }
        }
      else
        {
	  mps_secular_ga_update_coefficients (s);
          for (i = 0; i < s->n; i++)
            {
	      /* We may risk that NaN or inf have been introduced because of huge
	       * coefficients computed, so let's check it and in the case of failure 
	       * switch to DPE. */
	      if (cplx_check_fpe (sec->afpc[i]) || cplx_check_fpe (sec->bfpc[i]) ||
		  (cplx_mod (sec->afpc[i]) > 1.0e300) ||
		  (cplx_mod (sec->bfpc[i]) > 1.0e300))
		{
		  successful_regeneration = false;
		  if (s->debug_level & MPS_DEBUG_REGENERATION)
		      MPS_DEBUG (s, "Found floating point exception in regenerated coefficients, reusing old ones.");
		  
		  for (i = 0; i < s->n; i++)
		    {
		      cplx_set (sec->afpc[i], old_a[i]);
		      cplx_set (sec->bfpc[i], old_b[i]);
		    }
		  break;
		}

	      if (s->debug_level & MPS_DEBUG_REGENERATION)
		{
		  MPS_DEBUG_CPLX (s, sec->afpc[i], "sec->afpc[%d]", i);	      
		  MPS_DEBUG_CPLX (s, sec->bfpc[i], "sec->bfpc[%d]", i);
		}

	      mpc_set_cplx (s->root[i]->mvalue, s->root[i]->fvalue);
            }
        }

      cplx_vfree (old_a);
      cplx_vfree (old_b);
      cdpe_vfree (old_db);

      mps_secular_set_radii (s);

      break;

      /* If this is the DPE phase regenerate DPE coefficients */
    case dpe_phase:

      s->mpwp = MPS_SECULAR_EQUIVALENT_FP_PRECISION;

      /* Allocate old_a and old_b */
      old_da = cdpe_valloc (s->n);
      old_db = cdpe_valloc (s->n);

      /* Copy the old coefficients, and set the new
       * b_i with the current roots approximations. */
      for (i = 0; i < s->n; i++)
        {
          cdpe_set (old_da[i], sec->adpc[i]);
          cdpe_set (old_db[i], sec->bdpc[i]);
          cdpe_set (sec->bdpc[i], s->root[i]->dvalue);
	  mpc_set_cdpe (old_mb[i], old_db[i]);
          mpc_set_cdpe (sec->bmpc[i], sec->bdpc[i]);
        }

      mps_secular_ga_update_coefficients (s);

      /* Regeneration */
      if (!(successful_regeneration = mps_secular_ga_regenerate_coefficients_mp (s, old_db, old_mb)))
        {
	  MPS_DEBUG (s, "Regeneration failed");
	  for (i = 0; i < s->n; i++)
            {
              cdpe_set (sec->adpc[i], old_da[i]);
              cdpe_set (sec->bdpc[i], old_db[i]);
	      mpc_set_cdpe (old_mb[i], old_db[i]);
	      mpc_set_cdpe (sec->ampc[i], old_da[i]);
	      mpc_set_cdpe (sec->bmpc[i], old_db[i]);
            }

	  mps_secular_ga_update_coefficients (s);
	  break;
	}
      else
        {
	  mps_secular_ga_update_coefficients (s);
          for (i = 0; i < s->n; i++)
	    mpc_set_cdpe (s->root[i]->mvalue, s->root[i]->dvalue);
          mps_secular_set_radii (s);
        }
      
      if (s->debug_level & MPS_DEBUG_REGENERATION)
	{
	  for (i = 0; i < s->n; i++)
	    {
	      MPS_DEBUG_CDPE (s, sec->bdpc[i], "sec->bdpc[%d]", i);
	      MPS_DEBUG_CDPE (s, sec->adpc[i], "sec->adpc[%d]", i);
	    }
	}

      /* Free data */
      cdpe_vfree (old_da);
      cdpe_vfree (old_db);

      break;

    case mp_phase:

      /* Allocate old_a and old_b */
      old_ma = mpc_valloc (s->n);
      old_db = cdpe_valloc (s->n);

      mpc_vinit2 (old_ma, s->n, s->mpwp);

      /* Copy the old coefficients, and set the new
       * b_i with the current roots approximations. */
      for (i = 0; i < s->n; i++)
        {
          mpc_set (old_ma[i], sec->ampc[i]);
          mpc_set (old_mb[i], sec->bmpc[i]);
          mpc_set (sec->bmpc[i], s->root[i]->mvalue);
	  mpc_get_cdpe (old_db[i], old_mb[i]);
        }

      mps_secular_ga_update_coefficients (s);

      /* Regeneration */
      if (mps_secular_ga_regenerate_coefficients_mp (s, old_db, old_mb))
        {
	  mps_secular_ga_update_coefficients (s);
          /* Finally set radius according to new computed a_i coefficients,
           * if they are convenient   */
          mps_secular_set_radii (s);
        }

      if (s->debug_level & MPS_DEBUG_REGENERATION)
	{
	  MPS_DEBUG (s, "Dumping regenerated coefficients");
	  for (i = 0; i < s->n; i++)
	    {
	      MPS_DEBUG_MPC(s, s->mpwp, sec->ampc[i], "ampc[%d]", i);
	      MPS_DEBUG_MPC(s, s->mpwp, sec->bmpc[i], "bmpc[%d]", i);
	    }
	}
       
      mpc_vclear (old_ma, s->n);
      mpc_vfree (old_ma);
      rdpe_vfree (old_db);

      break;

    default:
      break;

    }                           /* End of switch (s->lastphase) */

  /* Sum execution time to the total counter */
#ifndef DISABLE_DEBUG
  s->regeneration_time += mps_stop_timer (my_clock);
#endif

  if (successful_regeneration)
    {
      MPS_DEBUG (s, "Setting again to true");
      for (i = 0; i < s->n; i++) 
	s->root[i]->again = true;
    }

  mpc_vclear (old_mb, s->n);
  mpc_vfree (old_mb);

  return successful_regeneration;
}
