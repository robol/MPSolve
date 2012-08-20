/*
 * secular-ga.c
 *
 *  Created on: 15/giu/2011
 *      Author: leonardo
 */

#include <mps/mps.h>
#include <math.h>
#include <string.h>

/**
 * @brief Update all the coefficients of the secular equation, and their
 * moduli, using the recomputed one stored in the multiprecision version.
 * 
 * @param s The mps_status of the computation.
 */
void
mps_secular_ga_update_coefficients (mps_status * s)
{
  int i;
  mps_secular_equation * sec = s->secular_equation;
  for (i = 0; i < s->n; ++i)
    {
      mpc_get_cplx (sec->afpc[i], sec->ampc[i]);
      mpc_get_cplx (sec->bfpc[i], sec->bmpc[i]);

      mpc_get_cdpe (sec->adpc[i], sec->ampc[i]);
      mpc_get_cdpe (sec->bdpc[i], sec->bmpc[i]);

      cdpe_mod (sec->aadpc[i], sec->adpc[i]);
      cdpe_mod (sec->abdpc[i], sec->bdpc[i]);

      sec->aafpc[i] = cplx_mod (sec->afpc[i]);
      sec->abfpc[i] = cplx_mod (sec->bfpc[i]);
    }
}



/**
 * @brief Check if iterations can terminate, i.e. if newton 
 * isolation has been reached, if the target was approximate. 
 * If the target was approximation then 
 * <code>mps_secular_improve ()</code> should be used to reach
 * the required precision.
 *
 * @param s The mps_status of the computation.
 */
mps_boolean
mps_secular_ga_check_stop (mps_status * s)
{
  MPS_DEBUG_THIS_CALL;

  int i;

    /* if  (!MPS_INPUT_CONFIG_IS_FP (s->secular_equation->input_structure) */
    /*   && s->lastphase != mp_phase) */
    /* return false; */

  for (i = 0; i < s->n; i++)
    {
      switch (s->lastphase)
        {
          /* Float case */
        case float_phase:
          if (!MPS_ROOT_STATUS_IS_COMPUTED (s, i))
            {
              MPS_DEBUG_WITH_INFO (s, "Root %d is not isolated, nor approximated, so we can't stop now.", i);
              return false;
            }
          break;

          /* Multiprecision and DPE case are the same, since the radii
           * are always RDPE. */
        case mp_phase:
          if (!MPS_ROOT_STATUS_IS_COMPUTED (s, i))
            {
              MPS_DEBUG_WITH_INFO (s, "Root %d is not isolated, nor approximated, so we can't stop now.", i);
	      MPS_DEBUG_WITH_INFO (s, "Status of root %d: %s", i, MPS_ROOT_STATUS_TO_STRING (s->root_status[i]));
              return false;
            }
          break;
        case dpe_phase:
	  MPS_DEBUG (s, "Status of root %d: %s", i, MPS_ROOT_STATUS_TO_STRING (s->root_status[i]));
          if (!MPS_ROOT_STATUS_IS_COMPUTED (s, i))
            {
              MPS_DEBUG_WITH_INFO (s, "Root %d is not isolated, nor approximated, so we can't stop now.", i);
              return false;
            }
          break;

        default:
          break;

        }
    }

  MPS_DEBUG_WITH_INFO (s, "Stop conditions were satisfied");
  return true;
}

/**
 * @brief Load the original coefficients of the secular equation.
 * 
 * This routine is used to load in the fields ampc and bmpc the original
 * coefficients of the secular equation, mainly in the improvement of the
 * isolated roots where the original coefficients (and not the rigenerated
 * one) are used. 
 */
void
mps_secular_ga_load_initial_coefficients (mps_status * s)
{
  MPS_DEBUG_THIS_CALL;

  int i;
  mps_secular_equation * sec = s->secular_equation;

 for (i = 0; i < s->n; i++)
    {
      if (MPS_INPUT_CONFIG_IS_SECULAR (s->input_config))
	{
	  if (MPS_INPUT_CONFIG_IS_FP (s->input_config))
	    {
	      mpc_set (sec->ampc[i],
		       sec->initial_ampc[i]);
	      mpc_set (sec->bmpc[i],
		       sec->initial_bmpc[i]);
	    }
	  else
	    {
	      mpc_set_q (sec->ampc[i],
			 sec->initial_ampqrc[i],
			 sec->initial_ampqic[i]);

	      mpc_set_q (sec->bmpc[i],
			 sec->initial_bmpqrc[i],
			 sec->initial_bmpqic[i]);
	    }
	}
    }
}

/**
 * @brief Improve the isolated roots. Should be called after isolation
 * is reached and confirmed by <code>mps_secular_ga_check_stop ()</code>.
 *
 * @param s The mps_status of the computation.
 */
void
mps_secular_ga_improve (mps_status * s)
{
  MPS_DEBUG_THIS_CALL;

  int i;
  mps_secular_equation  * sec = s->secular_equation;

  if (s->lastphase != mp_phase)
    mps_secular_switch_phase (s, mp_phase);

  /* Before improving with newton we shall reuse
   * the original coefficients */
  mps_secular_ga_load_initial_coefficients (s);
 
  mpc_t nwtcorr;
  cdpe_t ctmp;
  rdpe_t rtmp, old_rad, abroot;

  mpc_init2 (nwtcorr, s->mpwp);

  int starting_precision = s->mpwp;

  mps_secular_iteration_data it_data;
  it_data.local_ampc = sec->ampc;
  it_data.local_bmpc = sec->bmpc;

  for (i = 0; i < s->n; i++)
    {
      int j;

      /* Reset precision of coefficients and of the root we're 
       * interested in. */
      if (s->mpwp != starting_precision)
        {
          mps_secular_raise_coefficient_precision (s, starting_precision);
          mpc_set_prec (s->root[i]->mvalue, starting_precision);
          mpc_set_prec (nwtcorr, starting_precision);
          s->mpwp = starting_precision;
        }

      mpc_get_cdpe (ctmp, s->root[i]->mvalue);
      cdpe_mod (rtmp, ctmp);
      rdpe_div_eq (rtmp, s->root[i]->drad);

      /* Find correct digits and maximum number of iterations */
      int correct_digits = rdpe_log10 (rtmp) - 1;
      if (s->debug_level & MPS_DEBUG_IMPROVEMENT)
	MPS_DEBUG (s, "Root %d has %d correct digits", i, correct_digits);
      int iterations =
        log (1.0 * s->output_config->prec / correct_digits / LOG2_10) * LOG2_10 + 1;
      iterations = (iterations > 0) ? iterations : 0;

      if (s->debug_level & MPS_DEBUG_IMPROVEMENT)
        {
          if (iterations != 0)
	    {
              MPS_DEBUG (s, "Performing not more than %d iterations on root %d",
                         iterations, i);
	    }
          else
	    {
              MPS_DEBUG (s, "Not improving root %d, since it is already approximated", i);
	    }
        }

      for (j = 0; j < iterations; j++)
        { 
	  rdpe_set (old_rad, s->root[i]->drad);
          mps_secular_mnewton (s, s->root[i], nwtcorr,
                               &it_data, false);

	  /* Compute quadratic radius */
	  mpc_get_cdpe (ctmp, s->root[i]->mvalue); 
	  cdpe_mod (rtmp, ctmp); 
	  rdpe_div_eq (old_rad, rtmp); 
	  rdpe_mul_eq (old_rad, old_rad); 
	  rdpe_mul_eq (old_rad, rtmp);

	  /* Apply newton correction */
	  mpc_sub_eq (s->root[i]->mvalue, nwtcorr);

	  mpc_get_cdpe (ctmp, s->root[i]->mvalue);
	  cdpe_mod (abroot, ctmp);
	  
	  correct_digits *= 2;
	  /* rdpe_set_dl (s->root[i]->drad, 1.0, -correct_digits + 1); */
	  /* rdpe_add_eq (s->root[i]->drad, s->mp_epsilon); */
	  /* rdpe_mul_eq (s->root[i]->drad, abroot); */

          /* Debug iterations */
          if (s->debug_level & MPS_DEBUG_IMPROVEMENT)
            {
              MPS_DEBUG_MPC (s, 10, s->root[i]->mvalue, "s->mroot[%d]", i);
              MPS_DEBUG_RDPE (s, s->root[i]->drad, "s->drad[%d]", i);
            }

          /* Check if the approximation is already good. */
          mpc_get_cdpe (ctmp, s->root[i]->mvalue);
          cdpe_mod (rtmp, ctmp);
          rdpe_div (rtmp, s->root[i]->drad, rtmp);

          if (rdpe_le (rtmp, s->eps_out))
            {
              s->root_status[i] = MPS_ROOT_STATUS_APPROXIMATED;
              break;
            }
          else
            {
              s->mpwp *= 2;
              mps_secular_raise_coefficient_precision (s, s->mpwp);
              mpc_set_prec (nwtcorr, s->mpwp);
              mpc_set_prec (s->root[i]->mvalue, s->mpwp);
	      mps_secular_ga_load_initial_coefficients (s);
            }
        }

      /* Since we have passed the bound of the maximum allowed iterations
       * and quadratic convergence was guaranteed, the root is now 
       * approximated. */
      s->root_status[i] = MPS_ROOT_STATUS_APPROXIMATED;
    }
  mpc_clear (nwtcorr);

}

/**
 * @brief MPSolve main function for the secular equation solving
 * using Gemignani's approach.
 *
 * @param s The mps_status of the computation.
 */
void 
mps_secular_ga_mpsolve (mps_status * s)
{
  int roots_computed = 0;
  int packet;
  int iteration_per_packet = s->max_it;
  int i;
  mps_boolean skip_check_stop = false;
  mps_boolean just_regenerated = false;
  mps_secular_equation *sec = mps_secular_equation_from_status (s);

  if (s->output_config->search_set != MPS_SEARCH_SET_COMPLEX_PLANE)
    {
      mps_error (s, 1, "The restricted search set is not supported using the algorithm MPS_ALGORITHM_SECULAR_GA.");
      return;
    }

  /* Check if the secular equation is allocate or if only the
   * polynomial is present. In the last case, allocate an empty
   * secular equation to hold the data during the computation. */
  if (!s->secular_equation)
    {
      s->secular_equation = mps_secular_equation_new_raw (s, s->monomial_poly->n);
      sec = mps_secular_equation_from_status (s);
    }

  mps_allocate_data (s);

  s->just_raised_precision = true;

#ifndef DISABLE_DEBUG
  /* Reset all time counters */
  s->regeneration_time = 0;
  s->fp_iteration_time = 0;
  s->dpe_iteration_time = 0;
  s->mp_iteration_time = 0;
  clock_t *total_clock = mps_start_timer ();
#endif

  /* Set the output desired for the output */
  rdpe_set_2dl (s->eps_out, 1.0, -s->output_config->prec);

  /* Set degree and allocate polynomial-related variables
   * to allow initializitation to be performed. */
  s->deg = s->n = sec->n;

  /* Manually set FILE* pointer for streams.
   * More refined options will be added later. */
  s->outstr = s->rtstr = stdout;
  packet = 0;

  /* Set the maximum possible radius even in DPE, since
   * we may be starting directly from the DPE phase */
  for (i = 0; i < s->n; i++)
    {
      s->root[i]->frad = DBL_MAX;
      rdpe_set (s->root[i]->drad, RDPE_BIG); 
      /* rdpe_set_d (s->root[i]->drad, DBL_MAX) */
    }
  
  /* Set initial cluster structure as no cluster structure. */
  mps_cluster_reset (s);

  /* Set phase */
  s->lastphase = s->input_config->starting_phase;

  /* Set the number of roots in and out from the target set. Since we do not
   * support this yet, we can say that all the roots are in the set. */
  s->count[0] = s->n;
  s->count[1] = 0;
  s->count[2] = 0;

  /* If the input was polynomial we need to determined the secular
   * coefficients */
  if (MPS_INPUT_CONFIG_IS_MONOMIAL (s->input_config))
    {
      mps_monomial_poly *p = s->monomial_poly;

      for (i = 0; i < s->n; i++)
	cplx_set (sec->bfpc[i], cplx_zero);

      /* Check data first */
      if (s->input_config->starting_phase == no_phase)
	{
	  char which_case;
	  mps_check_data (s, &which_case);

	  if (mps_status_has_errors (s))
	    return;

	  MPS_DEBUG_WITH_INFO (s, "Check data suggests starting phase should be %s", (which_case == 'f') ? "floating point" : "DPE phase");

	  if (which_case == 'f')
	    s->lastphase = float_phase;
	  else
	    s->lastphase = dpe_phase;
	}
      else
	s->lastphase = s->input_config->starting_phase;

      MPS_DEBUG_WITH_INFO (s, "Computing starting points");
      if (s->lastphase == float_phase)
	  mps_fstart (s, s->n, NULL, 0.0, 0.0, s->eps_out, p->fap);
      else
	mps_dstart (s, s->n, NULL, (__rdpe_struct *) rdpe_zero,
		    (__rdpe_struct *) rdpe_zero, s->eps_out,
		    p->dap);

      /* Check if we can manage to perform the recomputatio of the
       * coefficients. If in floating point, switch do DPE if it fail.
       */
      if (!mps_secular_ga_regenerate_coefficients (s))
	{
	  if (s->lastphase == float_phase)
	    {
	      MPS_DEBUG(s, "Switching to DPE phase since initial regeneration of the coefficients did not succeed.");
	      s->lastphase = dpe_phase;
	      mps_dstart (s, s->n, 0, (__rdpe_struct *) rdpe_zero,
			  (__rdpe_struct *) rdpe_zero, s->eps_out,
			  p->dap);
	      if (!mps_secular_ga_regenerate_coefficients (s))
		{
		  MPS_DEBUG_WITH_INFO (s, "Initial generation of the secular equation coefficients did not succeed");
		  return;
		}
	    }
	  else
	    {
	      MPS_DEBUG_WITH_INFO (s, "Initial generation of the secular equation coefficients did not succeed");
	      return;
	    }
	  just_regenerated = true;
	}
    }
  else
    {
      MPS_DEBUG_WITH_INFO (s, "Generated initial coefficients for the secular equation");
      s->lastphase = float_phase;
      for (i = 0; i < s->n; i++)
	s->root[i]->wp = 53;
    }

  /* Select initial approximations using the custom secular
   * routine and based on the phase selected by the user. */
  if (!just_regenerated)
    {
      MPS_DEBUG_WITH_INFO (s, "Computing starting points");
      switch (s->lastphase)
	{
	case  float_phase:
	  mps_secular_fstart (s, s->n, NULL, 0.0, 0.0, s->eps_out);
	  break;

	case dpe_phase:
	  mps_secular_dstart (s, s->n, NULL, (__rdpe_struct *) rdpe_zero, 
			      (__rdpe_struct *) rdpe_zero, s->eps_out); 
	  break; 

	case mp_phase:
	  mps_secular_mstart (s, s->n, NULL, (__rdpe_struct *) rdpe_zero,
			      (__rdpe_struct *) rdpe_zero, s->eps_out);
	  break;

	default: 
	  break;
	}
    }

  /* Set initial radius */
  mps_secular_set_radii (s);

  for (i = 0; i < s->n; i++)
    {
      s->root[i]->again = true;
    }

  /* Cycle until approximated */
  do
    {
      skip_check_stop = false;
      s->secular_equation->best_approx = false;

      /* Perform an iteration of floating point Aberth method */
      switch (s->lastphase)
        {
        case float_phase:
          MPS_DEBUG_WITH_INFO (s, "Starting floating point iterations");
          roots_computed = mps_secular_ga_fiterate (s, iteration_per_packet, just_regenerated);
          /* If the computation fails we need to switch to DPE so do not
           * break here, but continue the cycle. */
          if (roots_computed != -1)
	    break;

        case dpe_phase:
          MPS_DEBUG_WITH_INFO (s, "Starting DPE iterations");
          roots_computed = mps_secular_ga_diterate (s, iteration_per_packet, just_regenerated);
          break;

        case mp_phase:
          MPS_DEBUG_WITH_INFO (s, "Starting MP iterations");
          roots_computed = mps_secular_ga_miterate (s, iteration_per_packet, just_regenerated);
          break;

        default:
          break;
        }

      /* Increase the packet counter */
      packet++;
	  
      /* Check thet we haven't passed the maximum number of allowed iterations */
      if (packet > s->max_pack)
	{
	  mps_error (s, 1, "Maximum number of iteration passed. Aborting.");
	  return;
	}
      
      /* Check if all roots were approximated with the
       * given input precision                      */      
      if (!just_regenerated)
	{
	  if (mps_secular_ga_check_stop (s)) 
	    break; 
	  else 
	    skip_check_stop = true; 
	}

       /* If the iterations has ended in less than 2 * not_computed_roots iterations
	* and we have just regenerated the coefficients, we should increase precision. */
       if (sec->best_approx)
	 {
	   skip_check_stop = false;
	   
	   /* Going to multiprecision if we're not there yet */
	   if (s->lastphase != mp_phase)
	     mps_secular_switch_phase (s, mp_phase);
          else
            {
              /* Raising precision otherwise */
              mps_secular_raise_precision (s, 2 * s->mpwp);
            }

	   if (MPS_INPUT_CONFIG_IS_MONOMIAL (s->input_config)) 
	     { 
	       MPS_DEBUG (s, "Performing restart phase");
	       mps_secular_restart (s);
	     }
	   
	   MPS_DEBUG (s, "Triggering regeneration");
	   if (!mps_secular_ga_regenerate_coefficients (s)) 
	     {
	       MPS_DEBUG (s, "Regeneration failed");
	     }

	   just_regenerated = true;
	   sec->best_approx = false;

	   /* Set the packet counter to zero, we are restarting */
	   packet = 0;
        }

      /* If we can't stop recompute coefficients in higher precision and
       * continue to iterate, unless the best approximation possible in
       * this precision has been reached. In that case increase the precision
       * of the computation. */
       s->just_raised_precision = false;

       /* Check if all the roots are approximated or, if we have done more than 4 packets
	* of iterations without finding all of them, if at least we are near to the result. */
       if (!just_regenerated)
	 {
	   if (mps_secular_ga_regenerate_coefficients (s))
	     {
	       just_regenerated = true;
	       skip_check_stop = false;
	     }
	   else
	     {
	       MPS_DEBUG (s, "Raising precision because regeneration failed");
	       
	       skip_check_stop = false;
	       
	       /* Going to multiprecision if we're not there yet */
	       if (s->lastphase != mp_phase)
		 {
		   mps_secular_switch_phase (s, mp_phase);
		 }
	       else
		 {
		   /* Raising precision otherwise */
		   mps_secular_raise_precision (s, 2 * s->mpwp);
		   mps_secular_ga_regenerate_coefficients (s);
		 }
	       
	       just_regenerated = true;
	       sec->best_approx = false;
	       
	       /* Set the packet counter to zero, we are restarting */
	       packet = 0;
	     }
	 }

       just_regenerated = false;
    }
  while (skip_check_stop || !mps_secular_ga_check_stop (s));

  mps_copy_roots (s);

  /* Finally improve the roots if approximation is required */
  if (s->output_config->goal == MPS_OUTPUT_GOAL_APPROXIMATE)
    {
      clock_t *my_timer = mps_start_timer ();
      mps_improve (s);
      unsigned int improve_time = mps_stop_timer (my_timer);
      if (s->debug_level & MPS_DEBUG_TIMINGS)
	{
	  MPS_DEBUG (s, "mps_improve took %u ms", improve_time);
	}
    }

  mps_restore_data (s);  

  /* Debug total time taken but only if debug is enabled */
#ifndef DISABLE_DEBUG
  long total_time = mps_stop_timer (total_clock);
  if (s->debug_level & MPS_DEBUG_TIMINGS)
    {
      MPS_DEBUG (s, "Time used for regeneration: %ld ms",
                 s->regeneration_time);
      MPS_DEBUG (s, "Time used in floating point iterations: %ld ms",
                 s->fp_iteration_time);
      MPS_DEBUG (s, "Time used in DPE iterations: %ld ms",
                 s->dpe_iteration_time);
      MPS_DEBUG (s, "Time used in multiprecision iterations: %ld ms",
                 s->mp_iteration_time);
      MPS_DEBUG (s, "Total time using MPSolve: %ld ms",
                 total_time);
    }
#endif

}
