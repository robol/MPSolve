#include <mps/debug.h>
#include <mps/core.h>
#include <mps/link.h>
#include <mps/secular.h>
#include <mps/debug.h>
#include <math.h>
#include <string.h>

/**
 * @brief Determine the precision required to bound the relative
 * error on the regenerated coefficients.
 */
int
mps_secular_ga_required_regenerations_bits (mps_status * s)
{
  rdpe_t root_epsilon;
  rdpe_t pol_eps;
  rdpe_t regeneration_epsilon;
  rdpe_t total_eps;
  int required_bits;
  int i, j, wp;

  /* Workaround to make setting the multiplier easy */
  char * multiplier_env = getenv("MULTIPLIER");
  int multiplier = 2;
  if (multiplier_env)
    sscanf (multiplier_env, "%d", &multiplier);

  rdpe_t required_eps;
  rdpe_set_2dl (required_eps, 1.0, -s->mpwp);

  wp = s->mpwp;

  if (MPS_INPUT_CONFIG_IS_MONOMIAL (s->input_config))
    {
      mps_monomial_poly * p = s->monomial_poly;
      mps_secular_equation * sec = s->secular_equation;

      do {

	rdpe_set_2dl (root_epsilon, 1.0, -wp);
	rdpe_set (total_eps, rdpe_zero);

	if (s->debug_level & MPS_DEBUG_REGENERATION)
	  {
	    MPS_DEBUG (s, "Trying to compute regeneration error with wp = %d", wp);
	    MPS_DEBUG_RDPE (s, root_epsilon, "Epsilon in this wp is set to eps");
	  }

	for (i = 0; i < s->n; i++)
	  {
	    rdpe_t rtmp, rtmp2, apol;
	    cdpe_t ss, pol, ctmp;
	    
	    rdpe_set (regeneration_epsilon, root_epsilon);

	    /* Perform horner with checking of the error */
	    cdpe_set (pol, p->dpc[s->n]);
	    for (j = s->n - 1; j >= 0; j--)
	      {
		cdpe_mul (ss, sec->bdpc[i], pol);
		cdpe_add_eq (ss, p->dpc[j]);
		cdpe_div (ctmp, pol, ss);
		cdpe_mod (rtmp, ctmp);

		rdpe_set (rtmp2, root_epsilon);
		rdpe_add_eq (rtmp2, regeneration_epsilon);
		rdpe_mul_eq (rtmp, rtmp2);
		rdpe_add_eq (regeneration_epsilon, rtmp);

		cdpe_div (ctmp, p->dpc[j], ss);
		cdpe_mod (rtmp, ctmp);
		rdpe_mul_eq (rtmp, regeneration_epsilon);
		rdpe_add_eq (regeneration_epsilon, rtmp);

		cdpe_set (pol, ss);
	      }

	    /* cdpe_set (pol, p->dpc[s->n]); */
	    /* for (j = s->n - 1; j > 0; j--) */
	    /*   { */
	    /* 	cdpe_mul_eq (pol, sec->bdpc[i]); */
	    /* 	cdpe_add_eq (pol, p->dpc[j]); */
	    /*   } */
	    /* cdpe_add_eq (pol, p->dpc[0]); */

	    /* rdpe_set (apol, p->dap[s->n]); */
	    /* cdpe_mod (rtmp2, sec->bdpc[i]); */
	    /* for (j = s->n - 1; j > 0; j--) */
	    /*   { */
	    /* 	rdpe_mul_eq (apol, rtmp2); */
	    /* 	rdpe_add_eq (apol, p->dap[j]); */
	    /*   } */
	    /* rdpe_add_eq (apol, p->dap[0]); */

	    /* cdpe_mod (rtmp, pol); */
	    /* rdpe_div_eq (apol, rtmp); */
	    /* rdpe_mul (regeneration_epsilon, apol, root_epsilon); */
	    
	    /* Check if the new relative error is bigger than the 
	     * previous one. */
	    if (rdpe_gt (regeneration_epsilon, total_eps))
	      {
		rdpe_set (total_eps, regeneration_epsilon);
	      }
	  }

	if (s->debug_level & MPS_DEBUG_REGENERATION)
	  {
	    MPS_DEBUG_RDPE (s, regeneration_epsilon, "regeneration_epsilon");
	    MPS_DEBUG_RDPE (s, required_eps, "required epsilon");
	  }

      } while (rdpe_gt (regeneration_epsilon, required_eps) && (wp *= 2));

      return 2 * wp;
    }
  else if (MPS_INPUT_CONFIG_IS_SECULAR (s->input_config))
    return multiplier * s->mpwp;
  
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
 * @param bits The bits of precision used in regeneration
 */
int
mps_secular_ga_regenerate_coefficients_mp (mps_status * s, int bits)
{
  /* Declaration and initialization of the multprecision
   * variables that are used only in that case */
  int i, j;
  mps_boolean success = true;
  int coeff_wp = bits;
  int old_wp = s->mpwp;
  mps_secular_equation *sec = s->secular_equation;
  mps_monomial_poly *p = s->monomial_poly;
  rdpe_t eps_tmp, rtmp, sec_eps, rtmp2;
  cdpe_t cdtmp, cdtmp2;
  rdpe_t root_epsilon;

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
    }

  mps_secular_raise_coefficient_precision (s, coeff_wp);

  if (MPS_INPUT_CONFIG_IS_MONOMIAL (s->input_config))
    {
      mpc_t diff, prod_b, m_one, ctmp;      
      MPS_DEBUG_WITH_INFO (s, "Regenerating coefficients from polynomial input");

      /* Init multiprecision values */
      mpc_init2 (diff, coeff_wp);
      mpc_init2 (prod_b, coeff_wp);
      mpc_init2 (m_one, coeff_wp);
      mpc_init2 (ctmp, coeff_wp);

      s->mpwp = coeff_wp;

      mpc_set_ui (m_one, 0U, 0U);
      mpc_sub_eq (m_one, p->mfpc[s->n]);

      
      /*
       * The new coefficients of the secular equation can be computed
       * starting from the evaluation of the polynomial changed of
       * sign and divided for the product of the difference of the
       * b_i, i.e.:
       * 
       *   a_i = -p(b_i) / \prod_{i \neq j} (b_i - b_j)
       *
       */
      if (MPS_INPUT_CONFIG_IS_INTEGER (s->input_config) 
	  || MPS_INPUT_CONFIG_IS_RATIONAL (s->input_config))
	{
	  for (i = 0; i < s->n + 1; ++i)
	    {
	      mpf_set_q (mpc_Re (p->mfpc[i]), p->initial_mqp_r[i]);
	      mpf_set_q (mpc_Im (p->mfpc[i]), p->initial_mqp_i[i]);

	      /* MPS_DEBUG_MPC (s, 15, p->mfpc[i], "p->mfpc[%d]", i); */
	    }
	}

      for (i = 0; i < s->n; ++i)
	{
	  /* Set the epsilon on this a_i to zero */
	  rdpe_set (sec->dregeneration_epsilon[i], rdpe_zero);
	  
	  /* Evaluate the polynomial with absolute values to 
	   * compute the error */
	  mpc_get_cdpe (cdtmp, sec->bmpc[i]);
	  cdpe_mod (rtmp, cdtmp);
	  mps_aparhorner (s, s->n, rtmp, p->dap, p->spar, sec_eps, 0);

	  rdpe_set (sec_eps, p->dap[s->n]);
	  mpc_get_cdpe (cdtmp, sec->bmpc[i]);
	  cdpe_mod (rtmp2, cdtmp);
	  for (j = s->n - 1; j >= 0; j--)
	    {
	      rdpe_mul (rtmp, sec_eps, rtmp2);
	      rdpe_add (sec_eps, rtmp, p->dap[j]);
	    }
	  rdpe_mul_eq_d (sec_eps, s->n * 4);
	  rdpe_mul_eq (sec_eps, s->mp_epsilon);

	  /* Evaluate the polynomial with horner */
	  mps_parhorner (s, s->n, sec->bmpc[i], p->mfpc, p->spar, sec->ampc[i], 0);
	  mpc_set (sec->ampc[i], p->mfpc[s->n]); 
	  for (j = s->n - 1; j >= 0; j--) 
	    { 
	      mpc_div (ctmp, p->mfpc[j], sec->ampc[i]);
	      mpc_add_eq (ctmp, sec->bmpc[i]); 
	      mpc_mul_eq (sec->ampc[i], ctmp); 
	    }

	  /* MPS_DEBUG_MPC (s, coeff_wp, sec->bmpc[i], "sec->bmpc[%d]", i); */
	  /* MPS_DEBUG_MPC (s, coeff_wp, sec->ampc[i], "sec->ampc[%d]", i); */

	  /* Add to sec_eps the error induced by the approximation on the roots */
	   rdpe_set (rtmp, rdpe_zero); 
	   for (j = 0; j <= s->n; j++) 
	     { 
	       rdpe_mul_d (rtmp2, p->dap[j], j); 
	       rdpe_add_eq (rtmp, rtmp2); 
	     } 
	   rdpe_mul_eq (rtmp, root_epsilon); 
	   rdpe_add_eq (sec_eps, rtmp);

	  /* Compute relative error in polynomial horner */
	  mpc_get_cdpe (cdtmp, sec->ampc[i]);
	  cdpe_mod (rtmp, cdtmp);
	  rdpe_div_eq (sec_eps, rtmp);

	  /* Prepare some data for error computation, precisely |b_i| */
	  mpc_get_cdpe (cdtmp, sec->bmpc[i]);
	  cdpe_mod (rtmp, cdtmp);

	  /* Compute the difference of the b_i */
	  mpc_set_ui (prod_b, 1U, 0U);
	  for (j = 0; j < s->n; ++j)
	    {
	      if (i == j)
		continue;

	      mpc_sub (diff, sec->bmpc[i], sec->bmpc[j]);
	      mpc_mul_eq (prod_b, diff);

	      /* Compute the local error here */
	      mpc_get_cdpe (cdtmp2, diff);
	      cdpe_mod (rtmp2, cdtmp2);
	      mpc_get_cdpe (cdtmp2, sec->bmpc[j]);
	      cdpe_mod (eps_tmp, cdtmp2);
	      rdpe_add_eq (eps_tmp, rtmp);
	      rdpe_div_eq (eps_tmp, rtmp2);
	      rdpe_add_eq (sec->dregeneration_epsilon[i], eps_tmp);
	    }

	  /* MPS_DEBUG_MPC (s, 15, prod_b, "prod_b"); */
	  
	  /* Actually divide the result and store it in
	   * a_i, as requested. */
	  mpc_div_eq (sec->ampc[i], prod_b);
	  // mpc_mul_ui (sec->ampc[i], sec->ampc[i], -1U);
	  mpc_mul_eq (sec->ampc[i], m_one);

	  /* Debug computed coefficients */
	  if (s->debug_level & MPS_DEBUG_REGENERATION)
	    {
	      MPS_DEBUG_MPC (s, 10, sec->ampc[i], "a_%d", i);
	      MPS_DEBUG_MPC (s, 10, sec->bmpc[i], "b_%d", i);
	    }

	  /* MPS_DEBUG_RDPE (s, sec->dregeneration_epsilon[i], "error on b_%d differences", i); */
	  /* MPS_DEBUG_RDPE (s, sec_eps, "sec_eps"); */

	  /* Finalize error computation */
	  // rdpe_mul_eq (sec->dregeneration_epsilon[i], s->mp_epsilon);
	  rdpe_mul_eq (sec->dregeneration_epsilon[i], s->mp_epsilon);
	  rdpe_add_eq (sec->dregeneration_epsilon[i], sec_eps);

	  /* if (s->debug_level & MPS_DEBUG_REGENERATION) */
	  /*   { */
	  /*     MPS_DEBUG_RDPE (s, sec->dregeneration_epsilon[i], */
	  /* 		      "Relative error on a_%d", i); */
	  /*   } */
	}

      /* Clear requested storage */
      mpc_clear (diff);
      mpc_clear (prod_b);
      mpc_clear (m_one);
      mpc_clear (ctmp);
    }
  else if (MPS_INPUT_CONFIG_IS_SECULAR (s->input_config))
    {
      mpc_t prod_b, sec_ev;
      mpc_t ctmp, btmp;

      /* Init multiprecision variables */
      mpc_init2 (prod_b, coeff_wp);
      mpc_init2 (sec_ev, coeff_wp);
      mpc_init2 (ctmp, coeff_wp);
      mpc_init2 (btmp, coeff_wp);

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

	  /* if (s->debug_level & MPS_DEBUG_REGENERATION) */
	  /*   { */
	  /*     MPS_DEBUG_RDPE (s, sec_eps, "Relative error on sec_ev(b_%d)", i); */
	  /*     MPS_DEBUG_RDPE (s, sec->dregeneration_epsilon[i], */
	  /* 		      "Relative error on a_%d", i); */
	  /*   } */
	  
	  /* rdpe_set_2dl (rtmp, 1.0, 1.0);  */
	  /* if (rdpe_gt (sec->dregeneration_epsilon[i], rtmp))   */
	  /*   {   */
	  /*     success = false;   */
	  /*     goto regenerate_m_exit;   */
	  /*   } */
	}

    regenerate_m_exit:

      /* Free data */
      mpc_clear (prod_b);
      mpc_clear (sec_ev);
      mpc_clear (ctmp);
      mpc_clear (btmp);


    } /* End of if (MPS_INPUT_CONFIG_IS_SECULAR (s->input_config)) */

  s->mpwp = old_wp;
  
  mps_secular_raise_coefficient_precision (s, s->mpwp);
  rdpe_set_2dl (s->mp_epsilon, 1.0, -s->mpwp + 1);

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
  int i, j, bits;
  mps_boolean successful_regeneration = false;

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

      /* Copy the old coefficients, and set the new
       * b_i with the current roots approximations. */
      for (i = 0; i < s->n; i++)
        {
          cplx_set (old_a[i], sec->afpc[i]);
          cplx_set (old_b[i], sec->bfpc[i]);
	  mpc_set_cplx (sec->bmpc[i], s->froot[i]);
        }

      mps_secular_ga_update_coefficients (s);

      /* Regeneration */
      bits = mps_secular_ga_required_regenerations_bits (s);
      if (!(successful_regeneration = mps_secular_ga_regenerate_coefficients_mp (s, bits)))
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
		    {
		      MPS_DEBUG (s, "Found floating point exception in regenerated coefficients, reusing old ones.");
		    }
		  
		  for (i = 0; i < s->n; i++)
		    {
		      cplx_set (sec->afpc[i], old_a[i]);
		      cplx_set (sec->bfpc[i], old_b[i]);
		    }
		  break;
		}

	      MPS_DEBUG_CPLX (s, sec->afpc[i], "sec->afpc[%d]", i);	      
	      MPS_DEBUG_CPLX (s, sec->bfpc[i], "sec->bfpc[%d]", i);
            }

          mps_secular_set_radii (s);
        }

      cplx_vfree (old_a);
      cplx_vfree (old_b);

      mps_secular_fstart (s, s->n, 0, 0, 0, s->eps_out);

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
      bits = mps_secular_ga_required_regenerations_bits (s);
      if (!(successful_regeneration = mps_secular_ga_regenerate_coefficients_mp (s, bits)))
        {
            for (i = 0; i < s->n; i++)
            {
              cdpe_set (sec->adpc[i], old_da[i]);
              cdpe_set (sec->bdpc[i], old_db[i]);
            }
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

      mps_secular_dstart (s, s->n, 0, (__rdpe_struct *) rdpe_zero,
                          (__rdpe_struct *) rdpe_zero, s->eps_out);

      break;

    case mp_phase:

      /* Allocate old_a and old_b */
      old_ma = mpc_valloc (s->n);
      old_mb = mpc_valloc (s->n);

      mpc_vinit2 (old_ma, s->n, 2 * s->mpwp);
      mpc_vinit2 (old_mb, s->n, 2 * s->mpwp);

      /* Copy the old coefficients, and set the new
       * b_i with the current roots approximations. */
      for (i = 0; i < s->n; i++)
        {
          mpc_set (old_ma[i], sec->ampc[i]);
          mpc_set (old_mb[i], sec->bmpc[i]);
          mpc_set (sec->bmpc[i], s->mroot[i]);
        }

      mps_secular_ga_update_coefficients (s);

      /* Regeneration */
      bits = mps_secular_ga_required_regenerations_bits (s);
      if (mps_secular_ga_regenerate_coefficients_mp (s, bits))
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

      mps_secular_mstart (s, s->n, 0, (__rdpe_struct *) rdpe_zero,
                          (__rdpe_struct *) rdpe_zero, s->eps_out);
      break;

    default:
      break;

    }                           /* End of switch (s->lastphase) */

  /* If the coefficients have been regenerated, all the roots are candidates
   * for being iterated more */
  for (int i = 0; i < s->n; ++i)
      s->again[i] = true;

  /* Sum execution time to the total counter */
#ifndef DISABLE_DEBUG
  s->regeneration_time += mps_stop_timer (my_clock);
#endif

  if (successful_regeneration)
    {
      for (i = 0; i < s->n; i++)
	  s->again[i] = true;
    }
  else
    {
      for (i = 0; i < s->n; i++)
	rdpe_set (sec->dregeneration_epsilon[i], old_dregeneration_epsilon[i]);
    }


  rdpe_vfree (old_dregeneration_epsilon);
  mps_secular_set_radii (s);

  return successful_regeneration;
}
