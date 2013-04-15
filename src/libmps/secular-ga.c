/*
 * This file is part of MPSolve 3.0
 *
 * Copyright (C) 2001-2013, Dipartimento di Matematica "L. Tonelli", Pisa.
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
 * @brief Update all the coefficients of the secular equation, and their
 * moduli, using the recomputed one stored in the multiprecision version.
 * 
 * @param s The mps_context of the computation.
 */
void
mps_secular_ga_update_coefficients (mps_context * s)
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
 * @param s The mps_context of the computation.
 */
mps_boolean
mps_secular_ga_check_stop (mps_context * s)
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
          if (!MPS_ROOT_STATUS_IS_COMPUTED (s->root[i]->status))
            {
              MPS_DEBUG_WITH_INFO (s, "Root %d is not isolated, nor approximated, so we can't stop now.", i);
              return false;
            }
          break;

          /* Multiprecision and DPE case are the same, since the radii
           * are always RDPE. */
        case mp_phase:
          if (!MPS_ROOT_STATUS_IS_COMPUTED (s->root[i]->status))
            {
              MPS_DEBUG_WITH_INFO (s, "Root %d is not isolated, nor approximated, so we can't stop now.", i);
              MPS_DEBUG_WITH_INFO (s, "Status of root %d: %s", i, MPS_ROOT_STATUS_TO_STRING (s->root[i]->status));
              return false;
            }
          break;
        case dpe_phase:
          MPS_DEBUG (s, "Status of root %d: %s", i, MPS_ROOT_STATUS_TO_STRING (s->root[i]->status));
          if (!MPS_ROOT_STATUS_IS_COMPUTED (s->root[i]->status))
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
 * @brief MPSolve main function for the secular equation solving
 * using Gemignani's approach.
 *
 * @param s The mps_context of the computation.
 */
void 
mps_secular_ga_mpsolve (mps_context * s)
{
  int roots_computed = 0;
  int packet;
  int i;
  mps_boolean skip_check_stop = false;
  mps_boolean just_regenerated = false;
  mps_secular_equation *sec = mps_secular_equation_from_status (s);

  if (s->output_config->search_set != MPS_SEARCH_SET_COMPLEX_PLANE)
    {
      mps_error (s, "The restricted search set is not supported using the algorithm MPS_ALGORITHM_SECULAR_GA.");
      return;
    }

  if (s->output_config->root_properties != MPS_OUTPUT_PROPERTY_NONE)
    {
      mps_error (s, "The root properties detection is not supported using the algorithm MPS_ALGORITHM_SECULAR_GA.");
      return;
    }

  /* Check if the secular equation is allocated or if only the
   * polynomial is present. In the last case, allocate an empty
   * secular equation to hold the data during the computation. */
  if (!sec)
    {
      s->secular_equation = mps_secular_equation_new_raw (s, s->n);
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
  s->deg = s->n = s->active_poly->degree;

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

  /* If the input was polynomial we need to determine the secular
   * coefficients */
  if (!MPS_IS_SECULAR_EQUATION (s->active_poly))
    {
      mps_polynomial *p = s->active_poly;

      for (i = 0; i < s->n; i++)
        {
          cplx_set (sec->bfpc[i], cplx_zero);
          cdpe_set (sec->bdpc[i], cdpe_zero);
        }

      /* Check data first */
      if (s->input_config->starting_phase == no_phase)
        {
          char which_case;
          mps_check_data (s, &which_case);

          if (mps_context_has_errors (s))
            {
              mps_stop_timer (total_clock);
              return;
            }

          MPS_DEBUG_WITH_INFO (s, "Check data suggests starting phase should be %s", (which_case == 'f') ? "floating point" : "DPE phase");

          if (which_case == 'f')
            s->lastphase = float_phase;
          else
            s->lastphase = dpe_phase;
        }
      else
        s->lastphase = s->input_config->starting_phase;

      MPS_DEBUG_WITH_INFO (s, "Computing starting points and performing first Aberth packet");

      /* Perform a packet of Aberth iterations */
      switch (s->lastphase)
      {
        case float_phase:
          mps_polynomial_fstart (s, p);

          if (p->fnewton)
            mps_faberth_packet (s, p, false);

          break;

        case dpe_phase:
          mps_polynomial_dstart (s, p);

          if (p->dnewton)
            mps_daberth_packet (s, p, false);
          break;

        default:
          mps_error (s, "Unrecognized starting phase");
          mps_stop_timer (total_clock);
          return;
      }

      mps_cluster_analysis (s, p);

      if (mps_secular_ga_check_stop (s))
        goto cleanup;

      /* In the case where we started in DPE but the initial approximation are
       * representable as standard floating point numbers, go back to float_phase. */
      if (s->lastphase == dpe_phase)
        {
          mps_boolean really_need_dpe = false;
          rdpe_t module;

          for (i = 0; i < s->n; i++)
            {
              cdpe_mod (module, s->root[i]->dvalue);
              if (rdpe_gt (module, rdpe_maxd) || rdpe_lt (module, rdpe_mind))
                really_need_dpe = true;
            }

          if (! really_need_dpe) 
          {
            MPS_DEBUG_WITH_INFO (s, "Going back to float_phase because all the approximations "
              "are representable as standard floating point numbers.");
            s->lastphase = float_phase;
            for (i = 0; i < s->n; i++)
              {
                cdpe_get_x (s->root[i]->fvalue, s->root[i]->dvalue);
                s->root[i]->frad = DBL_MAX;
              }
          }
      }

      /* Check if we can manage to perform the recomputation of the
       * coefficients. If in floating point, switch do DPE if it fail.
       */
      if (!mps_secular_ga_regenerate_coefficients (s))
        {
          if (s->lastphase == float_phase)
            {
              MPS_DEBUG(s, "Switching to DPE phase since initial regeneration of the coefficients did not succeed.");
              s->lastphase = dpe_phase;
              mps_polynomial_dstart (s, p);

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
      
      for (i = 0; i < s->n; i++)
        {
          cplx_set (s->secular_equation->afpc[i], MPS_SECULAR_EQUATION (s->active_poly)->afpc[i]);
          cplx_set (s->secular_equation->bfpc[i], MPS_SECULAR_EQUATION (s->active_poly)->bfpc[i]);
          cdpe_set (s->secular_equation->adpc[i], MPS_SECULAR_EQUATION (s->active_poly)->adpc[i]);
          cdpe_set (s->secular_equation->bdpc[i], MPS_SECULAR_EQUATION (s->active_poly)->bdpc[i]);
          mpc_set (s->secular_equation->ampc[i], MPS_SECULAR_EQUATION (s->active_poly)->ampc[i]);
          mpc_set (s->secular_equation->bmpc[i], MPS_SECULAR_EQUATION (s->active_poly)->bmpc[i]);

          MPS_POLYNOMIAL (s->secular_equation)->degree = s->active_poly->degree;
          MPS_POLYNOMIAL (s->secular_equation)->structure = MPS_STRUCTURE_COMPLEX_FP;
        }

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
          mps_secular_fstart (s, sec);
          break;

        case dpe_phase:
          mps_secular_dstart (s, sec);
          break; 

        case mp_phase:
          mps_secular_mstart (s, sec);
          break;

        default: 
          break;
        }
    }

  for (i = 0; i < s->n; i++)
    {
      s->root[i]->again = true;
      s->root[i]->approximated = false;
    }

  /* Cycle until approximated */
  do
    {
      skip_check_stop = false;
      s->best_approx = false;

      /* Perform an iteration of floating point Aberth method */
      switch (s->lastphase)
        {
        case float_phase:
          MPS_DEBUG_WITH_INFO (s, "Starting floating point iterations");

        if (s->jacobi_iterations)
          roots_computed = mps_faberth_packet (s, MPS_POLYNOMIAL (sec), just_regenerated);
        else
          roots_computed = mps_secular_ga_fiterate (s, s->max_it, just_regenerated);

          /* If the computation fails we need to switch to DPE so do not
           * break here, but continue the cycle. */
          if (roots_computed != -1)
            break;

        case dpe_phase:
          MPS_DEBUG_WITH_INFO (s, "Starting DPE iterations");

        if (s->jacobi_iterations)
          roots_computed = mps_daberth_packet (s, MPS_POLYNOMIAL (sec), just_regenerated);
        else   
          roots_computed = mps_secular_ga_diterate (s, s->max_it, just_regenerated);

          break;

        case mp_phase:
          MPS_DEBUG_WITH_INFO (s, "Starting MP iterations");

        if (s->jacobi_iterations)
          roots_computed = mps_maberth_packet (s, MPS_POLYNOMIAL (sec), just_regenerated);
        else 
          roots_computed = mps_secular_ga_miterate (s, s->max_it, just_regenerated);

          break;

        default:
          break;
        }

      /* Increase the packet counter */
      packet++;
          
      /* Check thet we haven't passed the maximum number of allowed iterations */
      if (packet > s->max_pack)
        {
          mps_error (s, "Maximum number of iteration passed. Aborting.");
          mps_stop_timer (total_clock);
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
       if (s->best_approx)
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

           if (MPS_IS_MONOMIAL_POLY (s->active_poly))
             {   
               MPS_DEBUG (s, "Performing restart phase");
               mps_secular_restart (s);
             }

           if (!mps_secular_ga_regenerate_coefficients (s)) 
             {
               MPS_DEBUG (s, "Regeneration failed");
             }
           else
             just_regenerated = true;

           /* just_regenerated = true; */
           s->best_approx = false;

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
         {
           if (mps_secular_ga_regenerate_coefficients (s))
             {
               skip_check_stop = false;
               just_regenerated = true;
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
               
               /* just_regenerated = true; */
               s->best_approx = false;
               
               /* Set the packet counter to zero, we are restarting */
               packet = 0;
             }
         }
    }
  while (skip_check_stop || !mps_secular_ga_check_stop (s));

  cleanup:

  mps_copy_roots (s);
  mps_dump (s);

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
