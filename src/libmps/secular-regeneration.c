#include <mps/mps.h>
#include <math.h>
#include <string.h>

/**
 * @brief Update the working precision of a root, i.e. the variable
 * <code>s->rootwp[i]</code> with the given precision rounded to the 
 * closer multiple of 64 (since GMP does not handle intermediate precisions).
 *
 * @param s The <code>mps_status</code> of the computation.
 * @param i The index of the root whose precision must be updated.
 * @param wp The precision to set.
 */
long int
mps_secular_ga_update_root_wp (mps_status * s, int i, long int wp)
{
  pthread_mutex_lock (&s->precision_mutex);

  mps_secular_equation * sec = s->secular_equation;  
  mps_monomial_poly * p = s->monomial_poly;  

   s->rootwp[i] = ((wp - 1) / 64 + 1) * 64;

   if (s->data_prec_max < s->rootwp[i])   
     s->data_prec_max = s->rootwp[i];   

    if (s->debug_level & MPS_DEBUG_MEMORY)  
      MPS_DEBUG (s, "Setting wp for root %d to %ld bits", i, s->rootwp[i]);  

    pthread_mutex_lock (&sec->bmpc_mutex[i]);
    if (mpc_get_prec (sec->bmpc[i]) < s->rootwp[i])  
      mpc_set_prec (sec->bmpc[i], s->rootwp[i]);  
    pthread_mutex_unlock (&sec->bmpc_mutex[i]);

    pthread_mutex_lock (&sec->ampc_mutex[i]);
    if (mpc_get_prec (sec->ampc[i]) < s->rootwp[i])  
      mpc_set_prec (sec->ampc[i], s->rootwp[i]);  
    pthread_mutex_unlock (&sec->ampc_mutex[i]);
    
    if (p != NULL) 
      mps_monomial_poly_raise_precision (s, p, s->rootwp[i]); 
    
   pthread_mutex_unlock (&s->precision_mutex);
  return s->rootwp[i];
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
 * @param s The <code>mps_status</code> of the computation.
 * @param old_b A vector of the old \f$b_i\f$ coefficients.
 * @param old_mb The vector of the old \f$b_i\f$ in <code>mpc_t</code> version. Must be set
 * to NULL if we are not in a multiprecision phase. 
 */
mps_boolean *
mps_secular_ga_find_changed_roots (mps_status * s, cdpe_t * old_b, mpc_t * old_mb)
{
  cdpe_t diff;
  rdpe_t rtmp, rtmp2;
  rdpe_t root_epsilon;
  int i;

  mps_secular_equation * sec = s->secular_equation;
  mps_boolean * root_changed = mps_boolean_valloc (s->n);

  mpc_t mdiff;

  /* Distinguish the case where we are in multiprecision phase */
  if (old_mb == NULL)
    rdpe_set_d (root_epsilon, DBL_EPSILON);
  else
    {
      rdpe_set (root_epsilon, s->mp_epsilon);

      mpc_init2 (mdiff, s->mpwp);
    }
  
  for (i = 0; i < s->n; i++)
    {
      if ((s->status[i][0] != 'a' && s->status[i][0] != 'o' && s->status[i][0] != 'i') ||
	  s->just_raised_precision)
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
	}
      else
	cdpe_sub (diff, old_b[i], sec->bdpc[i]);

      cdpe_mod (rtmp, diff);
      cdpe_mod (rtmp2, sec->bdpc[i]);
      rdpe_div_eq (rtmp, rtmp2);
      if (rdpe_le (rtmp, root_epsilon))
	{
	  if (s->debug_level & MPS_DEBUG_REGENERATION)
	    {
	      MPS_DEBUG (s, "b_%d hasn't changed, so p(b_%d) will not be recomputed", i, i);
	    }
	  root_changed[i] = false;
	}
      else
	{
	  root_changed[i] = true;
	}
    }

  if (old_mb)
    mpc_clear (mdiff);

  return root_changed;
}

struct __mps_secular_ga_regenerate_coefficients_monomial_data {
  mps_status * s;
  cdpe_t * old_b;
  mps_boolean * root_changed;
  rdpe_t * root_epsilon;
  mps_boolean * success;
  int i;
};

void *
__mps_secular_ga_regenerate_coefficients_monomial_worker (void * data_ptr)
{
  struct __mps_secular_ga_regenerate_coefficients_monomial_data * data = data_ptr;

  mps_status * s = data->s;
  cdpe_t * old_b = data->old_b;
  mps_boolean * root_changed = data->root_changed;
  
  /* Pointers to the secular equation and the monomial_poly */
  mps_secular_equation * sec = s->secular_equation;
  mps_monomial_poly * p = s->monomial_poly;

  /* Multiprecision variables */
  mpc_t mprod_b, ctmp, mdiff, lc, my_b;

  /* Floating point temporary variables */
  cdpe_t prod_b, diff, cdtmp;
  rdpe_t root_epsilon;

  /* The precision of the temporary variables at the start of the computation. We can set
   * this to s->mpwp; */
  long int coeff_wp = s->data_prec_max;

  /* This variable is true if the regeneration succeeded. */
  mps_boolean success = true;

  int i = data->i, j;

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
  mps_with_lock (p->mfpc_mutex[s->n],
		 mpc_div_eq (lc, p->mfpc[s->n]);
		 );
      
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

      mps_secular_ga_update_root_wp (s, i, s->rootwp[i]); 

      pthread_mutex_lock (&sec->bmpc_mutex[i]);
      mps_mhorner_with_error2 (s, p, sec->bmpc[i], sec->ampc[i], relative_error, s->rootwp[i]);
      pthread_mutex_unlock (&sec->bmpc_mutex[i]);

      mpc_get_cdpe (cpol, sec->ampc[i]);
      cdpe_mod (rtmp, cpol);
      rdpe_div_eq (relative_error, rtmp);

      if (s->debug_level & MPS_DEBUG_REGENERATION)
	MPS_DEBUG_RDPE (s, relative_error, "Relative_error on p(b_%d) evaluation", i);
      
      while (rdpe_gt (relative_error, root_epsilon))
	{
	  /* Update the working precision of the selected root with a realistic estimate of the
	   * required precision to get a result exact to machine precision */
	  mps_secular_ga_update_root_wp (s, i, 1 + s->rootwp[i] + (rdpe_Esp (relative_error) - rdpe_Esp (root_epsilon)));
	  
	  /* Try to recompute the polynomial with the augmented precision and see if now relative_error matches */
	  pthread_mutex_lock (&sec->bmpc_mutex[i]);
	  mps_mhorner_with_error2 (s, p, sec->bmpc[i], sec->ampc[i], relative_error, s->rootwp[i]);
	  pthread_mutex_unlock (&sec->bmpc_mutex[i]);
	  
	  mpc_get_cdpe (cpol, sec->ampc[i]);
	  cdpe_mod (rtmp, cpol);
	  rdpe_div_eq (relative_error, rtmp);
	  
	  if (s->debug_level & MPS_DEBUG_REGENERATION)
	    MPS_DEBUG_RDPE (s, relative_error, "Relative_error on p(b_%d) evaluation", i);
	}

      if (mpc_get_prec (mdiff) < s->rootwp[i])
	{
	  mpc_set_prec (mdiff, s->rootwp[i]);
	  mpc_set_prec (mprod_b, s->rootwp[i]);
	  mpc_set_prec (lc, s->rootwp[i]);
	  mpc_set_prec (my_b, s->rootwp[i]);
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
	      pthread_mutex_lock (&sec->bmpc_mutex[j]);
	      if (mpc_get_prec (sec->bmpc[j]) < mpc_get_prec (sec->bmpc[i]))
		mpc_set_prec (sec->bmpc[j], mpc_get_prec (sec->bmpc[i]));
	      /* pthread_mutex_lock (&sec->bmpc_mutex[i]);  */
	      mpc_sub (mdiff, sec->bmpc[i], sec->bmpc[j]);
	      /* pthread_mutex_unlock (&sec->bmpc_mutex[i]);  */
	      pthread_mutex_unlock (&sec->bmpc_mutex[j]);
	      
	      /* If the difference is zero than regeneration cannot succeed, and means
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
	      MPS_DEBUG_MPC (s, s->rootwp[i], sec->ampc[i], "a_%d", i);
	      MPS_DEBUG_MPC (s, s->rootwp[i], sec->bmpc[i], "b_%d", i);
	    }
      
    } /* Close the case where the coefficient are not approximated or isolated */
      else
	{
	  /* We should not recompute the polynomial here since the approximation
	   * hasn't changed. */
	  cdpe_set (prod_b, cdpe_one); 
	  for (j = 1; j < s->n; j++) 
	    { 
	      if (i == j || !root_changed[i])
		continue; 
	      
	      cdpe_sub (diff, old_b[j], sec->bdpc[i]); 
	      cdpe_sub (cdtmp, sec->bdpc[j], sec->bdpc[i]); 
	      cdpe_div_eq (diff, cdtmp); 
	      cdpe_mul_eq (prod_b, diff);
	    }
	  
	  mpc_set_cdpe (mprod_b, prod_b); 
	  mpc_div_eq (sec->ampc[i], mprod_b);
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
 * @param s The <code>mps_status</code> of the computation
 * @param old_b The old \f$b_i\f$ coefficients of the secular equation used to recompute
 * the value of \f$a_i\f$ in the case where, denoting with \f$b_i\f$ the new coefficients, 
 * with \f$\tilde b_i\f$ the old ones and with \f$u\f$ the machine precision:
 * \f[|b_i - \tilde b_i| < u\f]
 * so there is no need to recompute the value of \f$p(b_i)\f$.
 * @param root_changed A vector of booleans that is <codefalse</code> on the components that
 * did not changed from the last regeneration. 
 */
mps_boolean
mps_secular_ga_regenerate_coefficients_monomial (mps_status * s, cdpe_t * old_b, mps_boolean * root_changed)
{
  MPS_DEBUG_THIS_CALL;

  int i;
  mps_boolean success = true;

  struct __mps_secular_ga_regenerate_coefficients_monomial_data * data = 
    mps_newv (struct __mps_secular_ga_regenerate_coefficients_monomial_data, s->n);

  for (i = s->n - 1; i >= 0; i--)
    {
      data[i].i = i;
      data[i].old_b = old_b;
      data[i].root_changed = root_changed;
      data[i].s = s;
      data[i].success = &success;
      mps_thread_pool_assign (s, s->pool, __mps_secular_ga_regenerate_coefficients_monomial_worker, 
       			      data + i);   
      /* __mps_secular_ga_regenerate_coefficients_monomial_worker (data + i); */
      
    }

  mps_thread_pool_wait (s, s->pool);

  free (data);
  mps_boolean_vfree (root_changed);

  return success;
}

/**
 * @brief This function recomputes the \f$a_i\f$ coefficients of the secular equation based
 * on the new \f$b_i\f$ that are set in <code>s->secular_equation->bmpc[i]</code>.
 *
 * @param s The <code>mps_status</code> of the computation.
 * @param old_b The old \f$b_i\f$ coefficients of the secular equation used to recompute
 * the value of \f$a_i\f$ in the case where, denoting with \f$b_i\f$ the new coefficients, 
 * with \f$\tilde b_i\f$ the old ones and with \f$u\f$ the machine precision:
 * \f[|b_i - \tilde b_i| < u\f]
 * so there is no need to recompute the value of \f$S(b_i) * \prod_{j \neq i} (b_i - b_j)\f$.
 * @param root_changed A vector of booleans that is <codefalse</code> on the components that
 * did not changed from the last regeneration. 
 */
mps_boolean 
mps_secular_ga_regenerate_coefficients_secular (mps_status * s, cdpe_t * old_b, mps_boolean * root_changed)
{
  MPS_DEBUG_THIS_CALL;
  
  /* The secular equation to update */
  mps_secular_equation * sec = s->secular_equation;

  /* Multiprecision local variables */
  mpc_t prod_b, sec_ev;
  mpc_t ctmp, btmp;

  int i, j;

  cdpe_t cdtmp, cdtmp2;
  rdpe_t sec_eps, rtmp, rtmp2, eps_tmp;
  rdpe_t root_epsilon;

  mps_boolean success = true;

  /* Create an old secular_equation that will be used to evaluate the old secular
  * equation in the new nodes.*/
  mps_secular_equation  * old_sec = mps_new (mps_secular_equation);
  rdpe_t error;

  /* Determine the epsilon of the roots, that is a bound to the precision that
   * we desire con a_i coefficients of the secular equation. */
  switch (s->lastphase)
    {
    case float_phase:
    case dpe_phase:
      rdpe_set_d (root_epsilon, DBL_EPSILON);
      break;
    case mp_phase:
      rdpe_set (root_epsilon, s->mp_epsilon);
      break;
    default:
      mps_error (s, 1, "Phase not supported found in mps_secular_ga_regenerate_coefficients_secular");
    }

  /* Init multiprecision variables */
  mpc_init2 (prod_b, s->mpwp);
  mpc_init2 (sec_ev, s->mpwp);
  mpc_init2 (ctmp, s->mpwp);
  mpc_init2 (btmp, s->mpwp);

  /* If the input was rational generate multiprecision floating
   * point coefficient from the ones that the user originally
   * provided */
  if (MPS_INPUT_CONFIG_IS_RATIONAL (s->input_config) ||
      MPS_INPUT_CONFIG_IS_INTEGER (s->input_config))
    {
      MPS_DEBUG_WITH_INFO (s, "Regenerating coefficients from the multiprecision input");
      for (i = 0; i < s->n; i++)
	{
	  mpf_set_q (mpc_Re (sec->initial_ampc[i]), sec->initial_ampqrc[i]);
	  mpf_set_q (mpc_Im (sec->initial_ampc[i]), sec->initial_ampqic[i]);
	  mpf_set_q (mpc_Re (sec->initial_bmpc[i]), sec->initial_bmpqrc[i]);
	  mpf_set_q (mpc_Im (sec->initial_bmpc[i]), sec->initial_bmpqic[i]);
	}
    }

  /* Compute the new a_i */
  for (i = 0; i < s->n; i++)
    {
      mpc_set_ui (prod_b, 1, 0);
      mpc_set_ui (sec_ev, 0, 0);

      /* Set the regeneration of epsilon of this root to zero and keep
       * the module of b[i] in here since will be used later */
      rdpe_set (sec->dregeneration_epsilon[i], rdpe_zero);
      rdpe_set (sec_eps, rdpe_zero);
      mpc_get_cdpe (cdtmp, sec->bmpc[i]);
      cdpe_mod (eps_tmp, cdtmp);

      for (j = 0; j < sec->n; j++)
	{
	  /* Compute 1 / (b_i - old_b_j) */
	  mpc_sub (btmp, sec->bmpc[i], sec->initial_bmpc[j]);

	  /* Compute the local relative error that is introduce here, as 
	   * eps = (|b_i| + |old_b_j|) * mp_eps / (b_i - old_b_j) + mp_eps*/
	  mpc_get_cdpe (cdtmp, sec->initial_bmpc[j]);
	  cdpe_mod (rtmp, cdtmp);
	  rdpe_add_eq (rtmp, eps_tmp);
	  mpc_get_cdpe (cdtmp2, btmp);
	  cdpe_mod (rtmp2, cdtmp2);
	  rdpe_div_eq (rtmp, rtmp2);
	  rdpe_mul_eq (sec->dregeneration_epsilon[i], rtmp);

	  /* If b - old_b is zero, abort the computation */
	  if (mpc_eq_zero (btmp))
	    {
	      success = false;
	      MPS_DEBUG_WITH_INFO (s,
				   "Cannot regenerate coefficients, reusing old ones and setting best_approx to true.");
	      s->secular_equation->best_approx = true;
	      goto regenerate_m_exit;
	    }

	  mpc_inv (ctmp, btmp);

	  /* Add a_j / (b_i - old_b_j) to sec_ev */
	  mpc_mul_eq (ctmp, sec->initial_ampc[j]);
	  mpc_add_eq (sec_ev, ctmp);

	  /* Save the module of a_j / (b_i - old_b_j) to compute the
	   * relative error in the secular equation evalutation later */
	  mpc_get_cdpe (cdtmp, ctmp);
	  cdpe_mod (rtmp, cdtmp);
	  rdpe_add_eq (sec_eps, rtmp);

	  /* Multiply prod_b for
	   * b_i - b_j if i \neq j and prod_old_b
	   * for b_i - old_b_i.  */
	  mpc_mul_eq (prod_b, btmp);
	  if (i != j)
	    {
	      mpc_sub (ctmp, sec->bmpc[i], sec->bmpc[j]);
	      mpc_div_eq (prod_b, ctmp);

	      /* Compute the error in here */
	      mpc_get_cdpe (cdtmp, sec->bmpc[j]);
	      cdpe_mod (rtmp, cdtmp);
	      rdpe_add_eq (rtmp, eps_tmp);
	      mpc_get_cdpe (cdtmp2, ctmp);
	      cdpe_mod (rtmp2, cdtmp2);
	      rdpe_div_eq (rtmp, rtmp2);
	      rdpe_add_eq (sec->dregeneration_epsilon[i], rtmp);
	    }
	}

      /* Compute the new a_i as sec_ev * prod_old_b / prod_b */
      mpc_sub_eq_ui (sec_ev, 1, 0);
      mpc_mul (sec->ampc[i], sec_ev, prod_b);

      /* Compute the error obtained */
      rdpe_add_eq (sec_eps, rdpe_one);
      rdpe_mul_eq (sec->dregeneration_epsilon[i], s->mp_epsilon);
      rdpe_mul_eq (sec_eps, s->mp_epsilon);

      /* Relative error */
      mpc_get_cdpe (cdtmp, sec_ev);
      cdpe_mod (rtmp, cdtmp);
      rdpe_div_eq (sec_eps, rtmp);

      /* Sum the two for the moltiplication */
      rdpe_add_eq (sec->dregeneration_epsilon[i], sec_eps);
    }

  /* Fill its multiprecision fields */
  old_sec->ampc = sec->initial_ampc;
  old_sec->bmpc = sec->initial_bmpc;
  old_sec->n = sec->n;

  for (i = 0; i < s->n; i++)
    {
      rdpe_t rtmp;
      cdpe_t cdtmp, cprod_b;
      int j;

      /* Compute the new a_i as S(b_i) * \prod (x - old_b_i) / \prod (x - b_i) */
      mps_secular_meval_with_error (s, old_sec, sec->bmpc[i], sec->ampc[i], error);

      mpc_get_cdpe (cdtmp, s->mroot[i]);
      cdpe_mod (rtmp, cdtmp);
      rdpe_div_eq (error, rtmp);

      while (rdpe_gt (error, root_epsilon))
	{
	  mps_secular_ga_update_root_wp (s, i, 1 + s->rootwp[i] + (rdpe_Esp (error) - rdpe_Esp (root_epsilon)));

	  mps_secular_meval_with_error (s, old_sec, sec->bmpc[i], sec->ampc[i], error);

	  mpc_get_cdpe (cdtmp, sec->bmpc[i]);
	  cdpe_mod (rtmp, cdtmp);
	  rdpe_div_eq (error, rtmp);

	  if (s->debug_level & MPS_DEBUG_REGENERATION)
	    MPS_DEBUG_RDPE (s, error, "Relative_error on p(b_%d) evaluation", i);
	}

      /* Now that we have evaluated correctly the secular equation we can
       * compute prod_b to finalize the computation */
      cdpe_set (cprod_b, cdpe_one);
      for (j = 0; j < s->n; j++)
	{
	  mpc_sub (ctmp, sec->bmpc[i], old_sec->bmpc[j]);
	  mpc_get_cdpe (cdtmp, ctmp);
	  cdpe_mul_eq (cprod_b, cdtmp);

	  if (i != j)
	    {
	      mpc_sub (ctmp, sec->bmpc[i], sec->bmpc[j]);
	      mpc_get_cdpe (cdtmp, ctmp);
	      cdpe_div_eq (cprod_b, cdtmp);
	    }
	}
      
      mpc_set_cdpe (prod_b, cprod_b);      
      mpc_mul_eq (sec->ampc[i], prod_b);
    }

 regenerate_m_exit:
  
  /* Free data */
  free (old_sec);
  mpc_clear (prod_b);
  mpc_clear (sec_ev);
  mpc_clear (ctmp);
  mpc_clear (btmp);

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
 * <code>s</code> is a pointer to the <code>mps_status</code>
 * is <code>mps_secular_ga_regenerate_coefficients()</code> and
 * is the one that should be used
 *
 * @param s The mps_status of the computation.
 * @param old_b Old \f$b_i\f$ coefficients from which we are regenerating from. These are used
 * to adjust the a_i for already approximated roots.
 * @param old_mb Old \f$b_i\f$ in the multiprecision version. These are used only in the case
 * where <code>s->lastphase == mp_phase</code> to determine if the roots have changed. Otherwise
 * it must be set to NULL.
 */
mps_boolean
mps_secular_ga_regenerate_coefficients_mp (mps_status * s, cdpe_t * old_b, mpc_t * old_mb)
{
  /* Declaration and initialization of the multprecision
   * variables that are used only in that case */
  mps_boolean success = true;

  /* Get the list of the changed roots so we can compute the value of the 
   * polynomial only in that approximations. */
  mps_boolean * root_changed = mps_secular_ga_find_changed_roots (s, old_b, old_mb);

  // TODO: Remove this
  /* int i;    */
  /* for (i = 0; i < s->n; i++)    */
  /*   root_changed[i] = true;    */
  
  if (MPS_INPUT_CONFIG_IS_MONOMIAL (s->input_config))
    {
      /* Regenerate the coefficients of the secular equation starting from the monomial input */
      success = mps_secular_ga_regenerate_coefficients_monomial (s, old_b, root_changed);
    }
  else if (MPS_INPUT_CONFIG_IS_SECULAR (s->input_config))
    {
      /* Regenerate coefficients starting from the old coefficients store in <code>sec->initial_bmpc[i]</code>
       * evaluating the old secular equation in the new points. */
      success = mps_secular_ga_regenerate_coefficients_secular (s, old_b, root_changed);
    }

  if (!success)
      MPS_DEBUG (s, "Regeneration of the coefficients failed");

  return success;
}

/**
 * @brief Regenerate \f$a_i\f$ and \f$b_i\f$ setting
 * \f$b_i = z_i\f$, i.e. the current root approximation
 * and recomputing \f$a_i\f$ accordingly.
 *
 * @param s The mps_status of the computation.
 */
mps_boolean
mps_secular_ga_regenerate_coefficients (mps_status * s)
{
  MPS_DEBUG_THIS_CALL;

  cplx_t *old_b, *old_a;
  cdpe_t *old_db, *old_da;
  mpc_t *old_ma, *old_mb;
  mps_secular_equation *sec;
  int i;
  mps_boolean successful_regeneration = true;

  sec = (mps_secular_equation *) s->secular_equation;

  /* Copy the old regeneration epsilon that may come handy
   * in the case the regeneration does not work */
  rdpe_t *old_dregeneration_epsilon = rdpe_valloc (s->n);

  for (i = 0; i < s->n; i++)
    rdpe_set (old_dregeneration_epsilon[i], sec->dregeneration_epsilon[i]);

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
	  mpc_set_cplx (sec->bmpc[i], s->froot[i]);
        }

      mps_secular_ga_update_coefficients (s);

      /* Regeneration */
      if (!(successful_regeneration = mps_secular_ga_regenerate_coefficients_mp (s, old_db, NULL)))
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
	      sec->fregeneration_epsilon[i] = rdpe_get_d (sec->dregeneration_epsilon[i]) + DBL_EPSILON;

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
            }

          mps_secular_set_radii (s);
        }

      cplx_vfree (old_a);
      cplx_vfree (old_b);
      cdpe_vfree (old_db);

      mps_secular_fstart (s, s->n, NULL, 0, 0, s->eps_out);

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
          cdpe_set (sec->bdpc[i], s->droot[i]);
          mpc_set_cdpe (sec->bmpc[i], sec->bdpc[i]);
        }

      mps_secular_ga_update_coefficients (s);

      /* Regeneration */
      if (!(successful_regeneration = mps_secular_ga_regenerate_coefficients_mp (s, old_db, NULL)))
        {
	  MPS_DEBUG (s, "Regeneration failed");
	  for (i = 0; i < s->n; i++)
            {
              cdpe_set (sec->adpc[i], old_da[i]);
              cdpe_set (sec->bdpc[i], old_db[i]);
	      
	      mpc_set_cdpe (sec->ampc[i], old_da[i]);
	      mpc_set_cdpe (sec->bmpc[i], old_db[i]);
            }

	  mps_secular_ga_update_coefficients (s);
	  MPS_DEBUG (s ,"Ahi ahi ahi");
	  break;
	}
      else
        {
	  mps_secular_ga_update_coefficients (s);
          for (i = 0; i < s->n; i++)
	      rdpe_add_eq_d (sec->dregeneration_epsilon[i], DBL_EPSILON);
          mps_secular_set_radii (s);
        }
      
      MPS_DEBUG_CDPE (s, sec->bdpc[0], "sec->bdpc[%d]", 0);
      MPS_DEBUG_CDPE (s, sec->adpc[0], "sec->adpc[%d]", 0);

      /* Free data */
      cdpe_vfree (old_da);
      cdpe_vfree (old_db);

      mps_secular_dstart (s, s->n, NULL, 
			  (__rdpe_struct *) rdpe_zero,
                          (__rdpe_struct *) rdpe_zero, s->eps_out);

      break;

    case mp_phase:

      /* Allocate old_a and old_b */
      old_ma = mpc_valloc (s->n);
      old_mb = mpc_valloc (s->n);
      old_db = cdpe_valloc (s->n);

      mpc_vinit2 (old_ma, s->n, 2 * s->mpwp);
      mpc_vinit2 (old_mb, s->n, 2 * s->mpwp);

      /* Copy the old coefficients, and set the new
       * b_i with the current roots approximations. */
      for (i = 0; i < s->n; i++)
        {
          mpc_set (old_ma[i], sec->ampc[i]);
          mpc_set (old_mb[i], sec->bmpc[i]);
          mpc_set (sec->bmpc[i], s->mroot[i]);
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

      if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)
	{
	  MPS_DEBUG (s, "Dumping regenerated coefficients");
	  for (i = 0; i < s->n; i++)
	    {
	      MPS_DEBUG_MPC(s, 15, sec->ampc[i], "sec->ampc[%d]", i);
	      MPS_DEBUG_MPC(s, 15, sec->bmpc[i], "sec->bmpc[%d]", i);
	    }
	}
       
      mpc_vclear (old_ma, s->n);
      mpc_vclear (old_mb, s->n);
      mpc_vfree (old_ma);
      mpc_vfree (old_mb);
      rdpe_vfree (old_db);

      /* if (MPS_INPUT_CONFIG_IS_SECULAR (s->input_config)) */
      mps_secular_mstart (s, s->n, NULL,
			  (__rdpe_struct *) rdpe_zero, 
			  (__rdpe_struct *) rdpe_zero, s->eps_out);
      /* else */
      /* 	mps_mrestart (s); */

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
      for (i = 0; i < s->n; i++)
	  s->again[i] = true;
      mps_secular_set_radii (s);
    }
  else
    {
      for (i = 0; i < s->n; i++)
	rdpe_set (sec->dregeneration_epsilon[i], old_dregeneration_epsilon[i]);
    }

  rdpe_vfree (old_dregeneration_epsilon);

  return successful_regeneration;
}
