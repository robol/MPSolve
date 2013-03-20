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


/**
 * @file
 * @brief File with the implementation of the driver routines
 * for MPSolve. 
 */

#include <math.h>
#include <float.h>
#include <mps/mps.h>

/**
 * @brief Main routine of the program that implements the algorithm
 * in the standard polynomial version.
 *
 * The program is divided into many parts
 * - Check the correctness of data, scale coefficients if
 *   needed, and select cases: the variable <code>which_case</code> is
 *   <code>'f'</code> or <code>'d'</code>  according to float or dpe case.
 * - Call msolve or dsolve according to the value of which_case.
 * - Allocate MP variables mfpc, mroot, drad (if needed).
 * - Start MPsolve loop
 *   - prepare data according to the current precision
 *      and to the data_type (density/sparsity/user)  
 *   - Call msolve with the current precision
 * - check for termination
 */
void
mps_standard_mpsolve (mps_context * s)
{
  int i, nzc;
  char which_case;
  mps_boolean d_after_f, computed;
  clock_t *my_timer = mps_start_timer ();

  mps_allocate_data (s);

  if (s->DOLOG)
    s->debug_level |= MPS_DEBUG_TRACE;

  /* == 1 ==  Setup variables, i.e. copy coefficients
     into dpr, dpc and similar. */
  mps_setup (s);

  s->lastphase = no_phase;
  computed = false;
  s->over_max = false;

  /* == 2 ==  Resume from pre-computed roots */
  if (s->resume) 
    {
      mps_error (s, 1, "Resume not supported yet");
      return;
    }

  /* == 3 ==  Check data and get starting phase */
  if (s->skip_float)
    which_case = 'd';
  else
    which_case = 'f';

  /* This variable is true if we need a dpe phase after the
   * float phase */
  d_after_f = false;

  /* Check if a dpe phase is needed and deflate polynomial */
  mps_check_data (s, &which_case);

  /* Check for errors in check data */
  if (mps_context_has_errors (s))
    return;

  rdpe_set_2dl (s->eps_out, 1.0, - s->output_config->prec);

  if (s->DOLOG)
    fprintf (s->logstr, "Which_case = %c, skip_float= %d\n", which_case,
             s->skip_float);

  /* == 4 ==  Float phase */
  if (which_case == 'f')
    {
      if (s->DOLOG)
        fprintf (s->logstr, "Float phase ...\n");
      mps_fsolve (s, &d_after_f);
      s->lastphase = float_phase;

      if (s->DOLOG)
        mps_dump (s);

      computed = mps_check_stop (s);
      if (computed && s->output_config->goal != MPS_OUTPUT_GOAL_APPROXIMATE)
        goto exit_sub;
      /* stop for COUNT and ISOLATE goals */
    }

  /* == 5 ==  DPE phase */
  if (which_case == 'd' || d_after_f)
    {                           /* DPE phase */
      if (s->DOLOG)
        fprintf (s->logstr, "DPE phase ...\n");
      /* If we are arriving from a float phase copy the floating points
       * roots approximations in the DPE root approximations. */
      if (d_after_f)
        for (i = 0; i < s->n; i++)
          {
            rdpe_set_d (s->root[i]->drad, s->root[i]->frad);
            cdpe_set_x (s->root[i]->dvalue, s->root[i]->fvalue);
          }
      s->lastphase = dpe_phase;
      mps_dsolve (s, d_after_f);

      if (s->DOLOG)
        mps_dump (s);

      computed = mps_check_stop (s);
      if (computed && s->output_config->goal != MPS_OUTPUT_GOAL_APPROXIMATE)
        goto exit_sub;
    }

  /* == 6 ==   Allocate MP variables mfpc, mroot, drad, mfppc, mfppc1
   * (the real input case is not implemented yet ) */
  MPS_DEBUG (s, "Starting MP phase");

  s->lastphase = mp_phase;
  
  /* ==== 6.1 initialize mp variables */
  mps_mp_set_prec (s, 2 * DBL_MANT_DIG);

  /* Prepare data according to the current working precision */
  mps_prepare_data (s, s->mpwp);

  /* ==== 6.2 set initial values for mp variables */
  for (i = 0; i < s->n; i++)
    {
      if (which_case == 'd' || d_after_f)
        mpc_set_cdpe (s->root[i]->mvalue, s->root[i]->dvalue);
      else
        {
          mpc_set_cplx (s->root[i]->mvalue, s->root[i]->fvalue);
          rdpe_set_d (s->root[i]->drad, s->root[i]->frad);
        }
    }
  if (computed && s->output_config->goal == MPS_OUTPUT_GOAL_APPROXIMATE)
    {
      MPS_DEBUG (s, "Exiting since the approximation are computed and the goal is MPS_OUTPUT_GOAL_APPROXIMATE");
      goto exit_sub;
    }

  MPS_DEBUG (s, "s->mpwp = %ld, s->mpwp_max = %ld", s->mpwp, s->mpwp_max);
  MPS_DEBUG (s, "s->input_config->prec = %ld", s->active_poly->prec);

  /* == 7 ==  Start MPsolve loop */
  s->mpwp = DBL_MANT_DIG;
  while (!computed && s->mpwp < s->mpwp_max && (s->active_poly->prec == 0 || s->mpwp
                                                < s->active_poly->prec))
    {
      s->mpwp *= 2;

      if (s->active_poly->prec != 0 && s->mpwp > s->active_poly->prec)
        s->mpwp = s->active_poly->prec + (int) (log (4.0 * s->n) / LOG2);

      if (s->mpwp > s->mpwp_max)
        {
          s->mpwp = s->mpwp_max;
          s->over_max = true;
        }

      if (s->DOLOG)
        fprintf (s->logstr, "MAIN: mp_loop: mpwp=%ld\n", s->mpwp);

      /* == 7.1 ==   prepare data according to the current precision */
      mps_mp_set_prec (s, s->mpwp);
      mps_prepare_data (s, s->mpwp);

      /* == 7.2 ==   Call msolve with the current precision */
      if (s->DOLOG)
        fprintf (s->logstr, "MAIN: now call msolve nclust=%ld\n", s->clusterization->n);
      mps_msolve (s);
      s->lastphase = mp_phase;

      /* if (s->DOLOG) dump(logstr); */

      if (s->DOLOG)
        {                       /* count isolated zeros */
          nzc = 0;
          for (i = 0; i < s->n; i++)
            {
              if (s->root[i]->status == MPS_ROOT_STATUS_ISOLATED || 
                  s->root[i]->status == MPS_ROOT_STATUS_APPROXIMATED)
                nzc++;
            }
          fprintf (s->logstr, "MAIN: isolated %d roots\n", nzc);
          fprintf (s->logstr, "MAIN: after msolve check stop\n");
        }

      /* == 7.3 ==  Check the stop condition */
      computed = mps_check_stop (s);
      mps_mmodify (s, true);

      /* == 7.4 ==  reset the status vector */
      for (i = 0; i < s->n; i++)
        if (s->root[i]->status == MPS_ROOT_STATUS_NEW_CLUSTERED)
          s->root[i]->status = MPS_ROOT_STATUS_CLUSTERED;
    }

  /* == 8 ==  Check for termination */
  if (!computed)
    {
      if (s->over_max)
        {
          s->over_max = true;
          /* mps_error (s, 1, "Reached the maximum working precision"); */
          MPS_DEBUG (s, "Reached the maximum working precision");
          goto exit_sub;
        }
      else
        {
          /* mps_warn (s, "Reached the input precision"); */
          MPS_DEBUG (s, "Reached the input precision");
          goto exit_sub;
        }
          
    }

exit_sub:

  /* == 9 ==  Check inclusion disks */
  if (computed && s->clusterization->n < s->n)
    if (!mps_inclusion (s))
      {
        mps_error (s, 1, "Unable to compute inclusion disks");
        return;
      }

  /* == 10 ==  Refine roots */
  if (computed && !s->over_max && s->output_config->goal == MPS_OUTPUT_GOAL_APPROXIMATE)
    {
      s->lastphase = mp_phase;
      mps_improve (s);
    }

  /* == 11 ==  Restore to highest used precision */
  if (s->lastphase == mp_phase)
    mps_restore_data (s);

  long total_time = mps_stop_timer (my_timer);
  MPS_DEBUG (s, "Total time using MPSolve: %lu ms", total_time);

  /* Finally copy the roots ready for output */
  mps_copy_roots (s);
}

/**
  * @brief Setup vectors and variables
  */
void
mps_setup (mps_context * s)
{
  int i;
  mps_monomial_poly *p = MPS_MONOMIAL_POLY(s->active_poly);
  mpf_t mptemp;
  mpc_t mptempc;

  if (s->DOLOG)
    {
      /* fprintf (s->logstr, "Goal      = %5s\n", s->goal); */
      /* fprintf (s->logstr, "Data type = %3s\n", s->data_type); */
      fprintf (s->logstr, "Degree    = %d\n", s->n);
      fprintf (s->logstr, "Input prec.  = %ld digits\n", (long) (s->active_poly->prec
                                                                 * LOG10_2));
      fprintf (s->logstr, "Output prec. = %ld digits\n", (long) (s->output_config->prec
                                                                 * LOG10_2));
    }

  /* setup temporary vectors */
   if (MPS_DENSITY_IS_SPARSE (s->active_poly->density)) 
     for (i = 0; i <= MPS_POLYNOMIAL (p)->degree; i++) 
       { 
         p->fap[i] = 0.0;
         p->fpr[i] = 0.0;
         rdpe_set (p->dap[i], rdpe_zero); 
         cplx_set (p->fpc[i], cplx_zero); 
         rdpe_set (p->dpr[i], rdpe_zero); 
         cdpe_set (p->dpc[i], cdpe_zero); 
       }

  /* Indexes of the first (and only) cluster start from
   * 0 and reach n */
  mps_cluster_reset (s);

  /* set input and output epsilon */
  rdpe_set_2dl (s->eps_in, 1.0, 1 - s->active_poly->prec);
  rdpe_set_2dl (s->eps_out, 1.0, 1 - s->output_config->prec);

  /* precision of each root */
  for (i = 0; i < s->n; i++)
    s->root[i]->wp = 53;

  /* output order info */
  for (i = 0; i < s->n; i++)
    s->order[i] = i;

  /* reset root counts */
  s->count[0] = s->count[1] = s->count[2] = 0;

  /* compute DPE approximations */
  if (!MPS_IS_MONOMIAL_POLY (s->active_poly))
    return;                     /* nothing to do */

  /* init temporary mp variables */
  mpf_init2 (mptemp, DBL_MANT_DIG);
  mpc_init2 (mptempc, DBL_MANT_DIG);

  /* main loop */
  s->skip_float = false;
  for (i = 0; i <= s->n; i++)
    {

      if (MPS_DENSITY_IS_SPARSE (s->active_poly->density) && !p->spar[i])
        continue;

      if (MPS_STRUCTURE_IS_REAL (s->active_poly->structure))
        {
          if (MPS_STRUCTURE_IS_RATIONAL (s->active_poly->structure) ||
              MPS_STRUCTURE_IS_INTEGER (s->active_poly->structure))
            {
              mpf_set_q (mptemp, p->initial_mqp_r[i]);
              mpf_get_rdpe (p->dpr[i], mptemp);
              /*#G GMP 2.0.2 bug begin */
              if (rdpe_sgn (p->dpr[i]) != mpq_sgn (p->initial_mqp_r[i]))
                rdpe_neg_eq (p->dpr[i]);
              /*#G GMP bug end */
            }

          if (MPS_STRUCTURE_IS_FP (s->active_poly->structure))
            mpf_get_rdpe (p->dpr[i], mpc_Re (p->mfpc[i]));
          
          cdpe_set_e (p->dpc[i], p->dpr[i], rdpe_zero);

          /* compute dap[i] and check for float phase */
          rdpe_abs (p->dap[i], p->dpr[i]);
          rdpe_abs (p->dap[i], p->dpr[i]);
          if (rdpe_gt (p->dap[i], rdpe_maxd)
              || rdpe_lt (p->dap[i], rdpe_mind))
            s->skip_float = true;

        }
      else if (MPS_STRUCTURE_IS_COMPLEX (s->active_poly->structure))
        {
          if (MPS_STRUCTURE_IS_RATIONAL (s->active_poly->structure) ||
              MPS_STRUCTURE_IS_INTEGER (s->active_poly->structure))
            {
              mpc_set_q (mptempc, p->initial_mqp_r[i], p->initial_mqp_i[i]);
              mpc_get_cdpe (p->dpc[i], mptempc);
              /*#G GMP 2.0.2 bug begin */
              if (rdpe_sgn (cdpe_Re (p->dpc[i])) != mpq_sgn (p->initial_mqp_r[i]))
                rdpe_neg_eq (cdpe_Re (p->dpc[i]));
              if (rdpe_sgn (cdpe_Im (p->dpc[i])) != mpq_sgn (p->initial_mqp_i[i]))
                rdpe_neg_eq (cdpe_Im (p->dpc[i]));
              /*#G GMP bug end */
            }
          else if (MPS_STRUCTURE_IS_FP (s->active_poly->structure))
              mpc_get_cdpe (p->dpc[i], p->mfpc[i]);
          
          /* compute dap[i] */
          cdpe_mod (p->dap[i], p->dpc[i]);
          if (rdpe_gt (p->dap[i], rdpe_maxd)
              || rdpe_lt (p->dap[i], rdpe_mind))
            s->skip_float = true;
          
        }
    }

  /* free temporary mp variables */
  mpf_clear (mptemp);
  mpc_clear (mptempc);

  /* adjust input data type */
  if (MPS_STRUCTURE_IS_FP (s->active_poly->structure) && s->skip_float)
    s->active_poly->structure = MPS_STRUCTURE_IS_REAL (s->active_poly->structure) ? 
      MPS_STRUCTURE_REAL_BIGFLOAT : MPS_STRUCTURE_COMPLEX_BIGFLOAT;

  /* prepare floating point vectors */
  if (!s->skip_float)
    for (i = 0; i <= MPS_POLYNOMIAL (p)->degree; i++)
      {
        if (MPS_DENSITY_IS_SPARSE (s->active_poly->density) || !p->spar[i])
          continue;
        if (MPS_STRUCTURE_IS_REAL (s->active_poly->structure))
          {
            p->fpr[i] = rdpe_get_d (p->dpr[i]);
            p->fap[i] = fabs (p->fpr[i]);
            cplx_set_d (p->fpc[i], p->fpr[i], 0.0);
          }
        else
          {
            cdpe_get_x (p->fpc[i], p->dpc[i]);
            p->fap[i] = cplx_mod (p->fpc[i]);
          }
      }
}

/**
 * @brief Check consistency of data and makes some basic adjustments.
 *
 * This routine check, for example, if there are zero roots in the polynomial
 * (i.e. no costant term) and deflates the polynomial if necessary (shifting
 * the coefficients). 
 *
 * It sets the value of the parameter <code>which_case</code> to <code>'f'</code>
 * if a floating point phase is enough, or to <code>'d'</code> if
 * a <code>dpe</code> phase is needed.
 * 
 * @param s The <code>mps_context</code> associated with the current computation.
 * @param which_case the address of the variable which_case;
 */
void
mps_check_data (mps_context * s, char *which_case)
{
  rdpe_t min_coeff, max_coeff, tmp;
  mps_monomial_poly *p = MPS_MONOMIAL_POLY (s->active_poly);
  cdpe_t ctmp;
  int i;

  /* case of user-defined polynomial */
  if (!MPS_IS_MONOMIAL_POLY (s->active_poly))
    {
      if (s->output_config->multiplicity)
        mps_error (s, 1,
                   "Multiplicity detection not yet implemented for user polynomial");
      if (s->output_config->root_properties)
        mps_error (s, 1,
                   "Real/imaginary detection not yet implemented for user polynomial");
      *which_case = 'd';
      return;
    }

  /* Check consistency of input */
  if (rdpe_eq (p->dap[s->n], rdpe_zero))
    {
      mps_warn (s, "The leading coefficient is zero");
      do
        (s->n)--;
      while (rdpe_eq (p->dap[s->n], rdpe_zero));
    }

  /* Compute min_coeff */
  if (rdpe_lt (p->dap[0], p->dap[s->n]))
    rdpe_set (min_coeff, p->dap[0]);
  else
    rdpe_set (min_coeff, p->dap[s->n]);

  /* Compute max_coeff and its logarithm */
  rdpe_set (max_coeff, p->dap[0]);
  for (i = 1; i <= s->n; i++)
    if (rdpe_lt (max_coeff, p->dap[i]))
      rdpe_set (max_coeff, p->dap[i]);
  s->lmax_coeff = rdpe_log (max_coeff);

  /*  Multiplicity and sep */
  if (s->output_config->multiplicity)
    {
    if (MPS_STRUCTURE_IS_INTEGER (s->active_poly->structure))
      {
        mps_compute_sep (s);
      }
    else if (MPS_STRUCTURE_IS_RATIONAL (s->active_poly->structure))
      {
        mps_warn (s, "The multiplicity option has not been yet implemented");
        s->sep = 0.0;
      }
    else 
      {
        mps_warn (s, "The input polynomial has neither integer nor rational");
        mps_warn (s, " coefficients: unable to compute multiplicities");
        s->sep = 0.0;
      }
    }

  /* Real/Imaginary detection */
  if (s->output_config->root_properties || 
      s->output_config->search_set == MPS_SEARCH_SET_REAL || 
      s->output_config->search_set == MPS_SEARCH_SET_IMAG)
    {
      if (MPS_STRUCTURE_IS_INTEGER (s->active_poly->structure))
        {
          mps_compute_sep (s);
        }
      else if (MPS_STRUCTURE_IS_RATIONAL (s->active_poly->structure))
        {
          mps_error (s, 1,
                    "The real/imaginary option has not been yet implemented for rational input");
          return;
          s->sep = 0.0;
        }
      else
        {
          mps_error (s, 1, "The input polynomial has neither integer nor rational "
                           "coefficients: unable to perform real/imaginary options");
          return;
          s->sep = 0.0;
        }
    }

  /* Select cases (dpe or floating point)
   * First normalize the polynomial (only the float version) */
  rdpe_div (tmp, max_coeff, min_coeff);
  rdpe_mul_eq_d (tmp, (double) (s->n + 1));
  rdpe_mul_eq (tmp, rdpe_mind);
  rdpe_div_eq (tmp, rdpe_maxd);

  if (rdpe_lt (tmp, rdpe_one))
    {
      mpc_t m_min_coeff;
      cdpe_t c_min_coeff;

      /* if  (n+1)*max_coeff/min_coeff < dhuge/dtiny -  float case */
      *which_case = 'f';
      rdpe_mul_eq (min_coeff, max_coeff);
      rdpe_mul (tmp, rdpe_mind, rdpe_maxd);
      rdpe_div (min_coeff, tmp, min_coeff);
      rdpe_sqrt_eq (min_coeff);

      rdpe_set (cdpe_Re (c_min_coeff), min_coeff);
      rdpe_set (cdpe_Im (c_min_coeff), rdpe_zero);

      mpc_init2 (m_min_coeff, mpc_get_prec (p->mfpc[0]));
      mpc_set_cdpe (m_min_coeff, c_min_coeff);

      /* min_coeff = sqrt(dhuge*dtiny/(min_coeff*max_coeff)) */
      for (i = 0; i <= s->n; i++)
        {
          /* Multiply the MP leading coefficient */
          mpc_mul_eq (p->mfpc[i], m_min_coeff);

          rdpe_mul (tmp, p->dap[i], min_coeff);
          p->fap[i] = rdpe_get_d (tmp);
          if (MPS_STRUCTURE_IS_COMPLEX (s->active_poly->structure))
            {
              cdpe_mul_e (ctmp, p->dpc[i], min_coeff);
              cdpe_get_x (p->fpc[i], ctmp);
            }
          else
            {
              /* fpr(i)=dpr(i)*min_coeff !! DARIO riattiva dopo impl. caso reale */
              cdpe_mul_e (ctmp, p->dpc[i], min_coeff);
              cdpe_get_x (p->fpc[i], ctmp);
            }
        }

        mpc_clear (m_min_coeff);
    }
  else
    *which_case = 'd';
}

/*********************************************************
 *      SUBROUTINE COMPUTE_SEP                            *
 *********************************************************/
void
mps_compute_sep (mps_context * s)
{
  s->sep = s->n * s->lmax_coeff;
  s->sep = -(s->sep) - s->n * (1 + log ((double) s->n)) / LOG2;
  if (s->DOLOG)
    fprintf (s->logstr, "Sep = %f\n", s->sep);
}
