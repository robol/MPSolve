/***********************************************************
 **       Multiprecision Polynomial Solver (MPSolve)       **
 **                 Version 2.2, May 2001                  **
 **                                                        **
 **                      Written by                        **
 **       Dario Andrea Bini and Giuseppe Fiorentino        **
 **       (bini@dm.unipi.it)  (fiorent@dm.unipi.it)        **
 **                                                        **
 ** (C) 2001, Dipartimento di Matematica, FRISCO LTR 21024 **
 ***********************************************************/

/**
 * @file
 * @brief File with the implementation of the driver routines
 * for MPSolve. 
 */

#include <math.h>
#include <float.h>
#include <mps/core.h>
#include <mps/interface.h>



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
mps_standard_mpsolve (mps_status * s)
{
  int i, nzc;
  char which_case;
  mps_boolean d_after_f, computed, over_max;
  clock_t *my_timer = mps_start_timer ();

  mps_allocate_data (s);

  if (s->DOLOG)
    s->debug_level |= MPS_DEBUG_TRACE;

  /* == 1 ==  Setup variables, i.e. copy coefficients
     into dpr, dpc and similar. */
  mps_setup (s);

  s->lastphase = no_phase;
  computed = false;
  over_max = false;

  /* == 2 ==  Resume from pre-computed roots */
  if (s->resume)                /* to complete */
    mps_error (s, 1, "Resume not supported yet");

  /* == 3 ==  Check data and get starting phase */
  if (s->skip_float)
    which_case = 'd';
  else
    which_case = 'f';

  /* This variable is true if we need a dpe phase after the
   * float phase */
  d_after_f = false;

  /* Check if a dpe phase is needed and deflate polynomial */
  if (s->data_type[0] != 'u')
    {
      mps_check_data (s, &which_case);
    }
  else
    {
      /* Check data but only if the user has provided a suitable
       * check_data routine. */
      if (s->check_data_usr != NULL)
        (*s->check_data_usr) (s, &which_case);
      else
        /* Otherwise fallback on the standard check_data routine */
        mps_check_data (s, &which_case);
    }

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
        mps_dump (s, s->logstr);

      computed = mps_check_stop (s);
      if (computed && s->goal[0] != 'a')
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
            rdpe_set_d (s->drad[i], s->frad[i]);
            cdpe_set_x (s->droot[i], s->froot[i]);
          }
      s->lastphase = dpe_phase;
      mps_dsolve (s, d_after_f);

      if (s->DOLOG)
        mps_dump (s, s->logstr);

      computed = mps_check_stop (s);
      if (computed && s->goal[0] != 'a')
        goto exit_sub;
    }

  /* == 6 ==   Allocate MP variables mfpc, mroot, drad, mfppc, mfppc1
   * (the real input case is not implemented yet ) */
  if (s->DOLOG)
    fprintf (s->logstr, "MP phase ...\n");

  /* ==== 6.1 initialize mp variables */
  mps_mp_set_prec (s, 2 * DBL_MANT_DIG);

  /* Prepare data according to the current working precision */
  mps_prepare_data (s, s->mpwp);

  /* ==== 6.2 set initial values for mp variables */
  for (i = 0; i < s->n; i++)
    {
      if (which_case == 'd' || d_after_f)
        mpc_set_cdpe (s->mroot[i], s->droot[i]);
      else
        {
          mpc_set_cplx (s->mroot[i], s->froot[i]);
          rdpe_set_d (s->drad[i], s->frad[i]);
        }
    }
  if (computed && s->goal[0] == 'a')
    goto exit_sub;

  /* == 7 ==  Start MPsolve loop */
  s->mpwp = DBL_MANT_DIG;
  while (!computed && s->mpwp < s->mpwp_max && (s->input_config->prec == 0 || s->mpwp
                                                < s->input_config->prec))
    {

      s->mpwp *= 2;

      if (s->input_config->prec != 0 && s->mpwp > s->input_config->prec)
        s->mpwp = s->input_config->prec + (int) (log (4.0 * s->n) / LOG2);

      if (s->mpwp > s->mpwp_max)
        {
          s->mpwp = s->mpwp_max;
          over_max = true;
        }

      if (s->DOLOG)
        fprintf (s->logstr, "MAIN: mp_loop: mpwp=%ld\n", s->mpwp);

      /* == 7.1 ==   prepare data according to the current precision */
      mps_mp_set_prec (s, s->mpwp);
      mps_prepare_data (s, s->mpwp);

      /* == 7.2 ==   Call msolve with the current precision */
      if (s->DOLOG)
        fprintf (s->logstr, "MAIN: now call msolve nclust=%d\n", s->nclust);
      mps_msolve (s);
      s->lastphase = mp_phase;

      /* if (s->DOLOG) dump(logstr); */

      if (s->DOLOG)
        {                       /* count isolated zeros */
          nzc = 0;
          for (i = 0; i < s->n; i++)
            {
              if (s->status[i][0] == 'i' || s->status[i][0] == 'a')
                nzc++;
            }
          fprintf (s->logstr, "MAIN: isolated %d roots\n", nzc);
          fprintf (s->logstr, "MAIN: after msolve check stop\n");
        }

      /* == 7.3 ==  Check the stop condition */
      computed = mps_check_stop (s);
      mps_mmodify (s);

      /* == 7.4 ==  reset the status vector */
      for (i = 0; i < s->n; i++)
        if (s->status[i][0] == 'C')
          s->status[i][0] = 'c';
    }

  /* == 8 ==  Check for termination */
  if (!computed)
    {
      if (over_max)
        mps_error (s, 1, "Reached the maximum working precision");
      else
        mps_warn (s, "Reached the input precision");
    }

exit_sub:

  /* == 9 ==  Check inclusion disks */
  if (computed && s->nclust < s->n)
    if (!mps_inclusion (s))
      mps_error (s, 1, "Unable to compute inclusion disks");

  /* == 10 ==  Refine roots */
  if (computed && !over_max && s->goal[0] == 'a')
    {
      s->lastphase = mp_phase;
      mps_improve (s);
    }

  /* == 11 ==  Restore to highest used precision */
  if (s->lastphase == mp_phase)
    mps_restore_data (s);

  MPS_DEBUG (s, "Total time using MPSolve: %lu ms", mps_stop_timer (my_timer));

  /* Finally copy the roots ready for output */
  mps_copy_roots (s);
}

/**
  * @brief Setup vectors and variables
  */
void
mps_setup (mps_status * s)
{
  int i;
  mps_monomial_poly *p = s->monomial_poly;
  tmpf_t mptemp;
  tmpc_t mptempc;

  if (s->DOLOG)
    {
      fprintf (s->logstr, "Goal      = %5s\n", s->goal);
      fprintf (s->logstr, "Data type = %3s\n", s->data_type);
      fprintf (s->logstr, "Degree    = %d\n", s->n);
      fprintf (s->logstr, "Input prec.  = %ld digits\n", (long) (s->input_config->prec
                                                                 * LOG10_2));
      fprintf (s->logstr, "Output prec. = %ld digits\n", (long) (s->output_config->prec
                                                                 * LOG10_2));
    }

  /* setup temporary vectors */
   if (s->data_type[0] == 's') 
     for (i = 0; i <= s->n; i++) 
       { 
         p->fap[i] = 0.0; 
         p->fpr[i] = 0.0; 
         rdpe_set (p->dap[i], rdpe_zero); 
         cplx_set (p->fpc[i], cplx_zero); 
         rdpe_set (p->dpr[i], rdpe_zero); 
         cdpe_set (p->dpc[i], cdpe_zero); 
       }  

  /* setup status and clusters so that there is only one cluster
   *  containing all the roots */
  for (i = 0; i < s->n; i++)
    {
      /* Set the i-th root as a clustered root */
      s->status[i][0] = 'c';
      /* Set the i-th root as 'uncertain', regarding the
       * inclusion in R or in iR. */
      s->status[i][1] = 'w';
      /* Set the i-th root as 'uncertain' (regarding the inclusion
       * in the search set). */
      s->status[i][2] = 'u';
      s->clust[i] = i;
    }

  /* Indexes of the first (and only) cluster start from
   * 0 and reach n */
  mps_cluster_reset (s);

  /* Check if the numer of thread is greater of the number of roots,
     and in that case decrease it */
  if (s->n_threads > s->n)
    s->n_threads = s->n;

  /* set input and output epsilon */
  rdpe_set_2dl (s->eps_in, 1.0, 1 - s->input_config->prec);
  rdpe_set_2dl (s->eps_out, 1.0, 1 - s->output_config->prec);

  /* precision of each root */
  for (i = 0; i < s->n; i++)
    s->rootwp[i] = 53;

  /* output order info */
  for (i = 0; i < s->n; i++)
    s->order[i] = i;

  /* reset root counts */
  s->count[0] = s->count[1] = s->count[2] = 0;

  /* compute DPE approximations */
  if (s->data_type[0] == 'u')
    return;                     /* nothing to do */

  /* init temporary mp variables */
  tmpf_init2 (mptemp, DBL_MANT_DIG);
  tmpc_init2 (mptempc, DBL_MANT_DIG);

  /* main loop */
  s->skip_float = false;
  for (i = 0; i <= s->n; i++)
    {

      if (s->data_type[0] == 's' && !p->spar[i])
        continue;

      switch (s->data_type[1])
        {                       /* switch 1 */

        case 'r':              /* Real */
          switch (s->data_type[2])
            {                   /* switch 2 */

            case 'i':          /* Real - Integer Coefs */
            case 'q':          /* Real - Rational Coefs */
              mpf_set_q (mptemp, p->initial_mqp_r[i]);
              mpf_get_rdpe (p->dpr[i], mptemp);
              /*#G GMP 2.0.2 bug begin */
              if (rdpe_sgn (p->dpr[i]) != mpq_sgn (p->initial_mqp_r[i]))
                rdpe_neg_eq (p->dpr[i]);
              /*#G GMP bug end */
              break;

            case 'f':          /* Real - Big/Float Coefs */
              mpf_get_rdpe (p->dpr[i], p->mfpr[i]);
              break;
            }                   /* switch 2 */

          cdpe_set_e (p->dpc[i], p->dpr[i], rdpe_zero);

          /* compute dap[i] and check for float phase */
          rdpe_abs (p->dap[i], p->dpr[i]);
          rdpe_abs (p->dap[i], p->dpr[i]);
          if (rdpe_gt (p->dap[i], rdpe_maxd)
              || rdpe_lt (p->dap[i], rdpe_mind))
            s->skip_float = true;

          break;                /* switch 1 */

        case 'c':              /* Complex */

          switch (s->data_type[2])
            {                   /* switch 3 */

            case 'i':          /* Complex - Integer Coefs */
            case 'q':          /* Complex - Rational Coefs */
              mpc_set_q (mptempc, p->initial_mqp_r[i], p->initial_mqp_i[i]);
              mpc_get_cdpe (p->dpc[i], mptempc);
              /*#G GMP 2.0.2 bug begin */
              if (rdpe_sgn (cdpe_Re (p->dpc[i])) != mpq_sgn (p->initial_mqp_r[i]))
                rdpe_neg_eq (cdpe_Re (p->dpc[i]));
              if (rdpe_sgn (cdpe_Im (p->dpc[i])) != mpq_sgn (p->initial_mqp_i[i]))
                rdpe_neg_eq (cdpe_Im (p->dpc[i]));
              /*#G GMP bug end */
              break;

            case 'f':          /* Complex - Big/Float Coefs */
              mpc_get_cdpe (p->dpc[i], p->mfpc[i]);
              break;

            }                   /* switch 3 */

          /* compute dap[i] */
          cdpe_mod (p->dap[i], p->dpc[i]);
          if (rdpe_gt (p->dap[i], rdpe_maxd)
              || rdpe_lt (p->dap[i], rdpe_mind))
            s->skip_float = true;

          break;

        }                       /* switch 1 */
    }                           /* for */

  /* free temporary mp variables */
  tmpf_clear (mptemp);
  tmpc_clear (mptempc);

  /* adjust input data type */
  if (s->data_type[2] == 'f' && s->skip_float)
    s->data_type[2] = 'b';

  /* prepare floating point vectors */
  if (!s->skip_float)
    for (i = 0; i <= p->n; i++)
      {
        if (s->data_type[0] == 's' || !p->spar[i])
          continue;
        if (s->data_type[1] == 'r')
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
 * @param s The <code>mps_status</code> associated with the current computation.
 * @param which_case the address of the variable which_case;
 */
void
mps_check_data (mps_status * s, char *which_case)
{
  rdpe_t min_coeff, max_coeff, tmp;
  mps_monomial_poly *p = s->monomial_poly;
  cdpe_t ctmp;
  int i;

  /* case of user-defined polynomial */
  if (s->data_type[0] == 'u')
    {
      if (s->goal[2] == 'm')
        mps_error (s, 1,
                   "Multiplicity detection not yet implemented for user polynomial");
      if (s->goal[3] != 'n')
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

  /* count number of zero roots (undeflated input polynomial) */
  s->zero_roots = 0;
  while (rdpe_eq (p->dap[s->zero_roots], rdpe_zero))
    (s->zero_roots)++;
  /* shift down input vectors */
  if (s->zero_roots)
    {
      for (i = 0; i <= s->n - s->zero_roots; i++)
        {
          rdpe_set (p->dap[i], p->dap[i + s->zero_roots]);
          p->fap[i] = p->fap[i + s->zero_roots];
          p->fpr[i] = p->fpr[i + s->zero_roots];
          cplx_set (p->fpc[i], p->fpc[i + s->zero_roots]);
          rdpe_set (p->dpr[i], p->dpr[i + s->zero_roots]);
          cdpe_set (p->dpc[i], p->dpc[i + s->zero_roots]);
          mpf_set (p->mfpr[i], p->mfpr[i + s->zero_roots]);
          mpc_set (p->mfpc[i], p->mfpc[i + s->zero_roots]);
          if (i < s->n - s->zero_roots)
            mpc_set (p->mfppc[i], p->mfppc[i + s->zero_roots]);
          mpq_set (p->initial_mqp_r[i], p->initial_mqp_r[i + s->zero_roots]);
          mpq_set (p->initial_mqp_i[i], p->initial_mqp_i[i + s->zero_roots]);
          p->spar[i] = p->spar[i + s->zero_roots];
        }
      s->n = s->n - s->zero_roots;
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
  if (s->goal[2] == 'm')
    switch (s->data_type[2])
      {
      case 'i':
        mps_compute_sep (s);
        break;
      case 'q':
        mps_warn (s, "The multiplicity option has not been yet implemented");
        s->sep = 0.0;
        break;
      default:
        mps_warn (s, "The input polynomial has neither integer nor rational");
        mps_warn (s, " coefficients: unable to compute multiplicities");
        s->sep = 0.0;
        break;
      }

  /* Real/Imaginary detection */
  if (s->goal[3] != 'n' || s->goal[1] == 'R' || s->goal[1] == 'I')
    switch (s->data_type[2])
      {
      case 'i':
        mps_compute_sep (s);
        break;
      case 'q':
        mps_warn (s,
                  "The  real/imaginary option has not been yet implemented");
        s->sep = 0.0;
        break;
      default:
        mps_warn (s, "The input polynomial has neither integer nor rational");
        mps_warn (s,
                  " coefficients: unable to perform real/imaginary options");
        s->sep = 0.0;
        break;
      }

  /* Select cases (dpe or floating point)
   * First normalize the polynomial (only the float version) */
  rdpe_div (tmp, max_coeff, min_coeff);
  rdpe_mul_eq_d (tmp, (double) (s->n + 1));
  rdpe_mul_eq (tmp, rdpe_mind);
  rdpe_div_eq (tmp, rdpe_maxd);
  if (rdpe_lt (tmp, rdpe_one))
    {
      /* if  (n+1)*max_coeff/min_coeff < dhuge/dtiny -  float case */
      *which_case = 'f';
      rdpe_mul_eq (min_coeff, max_coeff);
      rdpe_mul (tmp, rdpe_mind, rdpe_maxd);
      rdpe_div (min_coeff, tmp, min_coeff);
      rdpe_sqrt_eq (min_coeff);
      /* min_coeff = sqrt(dhuge*dtiny/(min_coeff*max_coeff)) */
      for (i = 0; i <= s->n; i++)
        {
          rdpe_mul (tmp, p->dap[i], min_coeff);
          p->fap[i] = rdpe_get_d (tmp);
          if (s->data_type[1] == 'c')
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
    }
  else
    *which_case = 'd';
}

/*********************************************************
 *      SUBROUTINE COMPUTE_SEP                            *
 *********************************************************/
void
mps_compute_sep (mps_status * s)
{
  s->sep = s->n * s->lmax_coeff;
  s->sep = -(s->sep) - s->n * (1 + log ((double) s->n)) / LOG2;
  if (s->DOLOG)
    fprintf (s->logstr, "Sep = %f\n", s->sep);
}
