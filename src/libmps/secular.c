/*
 * mps_secular.c
 *
 *  Created on: 10/apr/2011
 *      Author: leonardo
 */

#include <mps/mps.h>

/**
 * @brief Utility routine that dumps the DPE coefficients of the secular equation
 * passed as second argument. It is only used for debugging purpose.
 *
 * @param s The mps_status of the computation
 * @param sec The secular equation whose DPE coefficients will be dumped
 */
void
mps_secular_dump (mps_status * s, mps_secular_equation * sec)
{
  int i;
  MPS_DEBUG (s, "Dumping secular equation:");

  switch (s->lastphase)
    {
    case float_phase:
      for (i = 0; i < s->n; ++i)
	{
	  MPS_DEBUG_CPLX (s, sec->afpc[i], "sec->afpc[%d]", i);
	  MPS_DEBUG_CPLX (s, sec->bfpc[i], "sec->bfpc[%d]", i);
	}
      break;
    case dpe_phase:
      for (i = 0; i < sec->n; i++)
	{
	  MPS_DEBUG_CDPE (s, sec->adpc[i], "sec->adpc[%d]", i);
	  MPS_DEBUG_CDPE (s, sec->bdpc[i], "sec->bdpc[%d]", i);
	}
      break;
    case mp_phase:
      for (i = 0; i < s->n; ++i)
	{
	  MPS_DEBUG_MPC (s, 20, sec->ampc[i], "sec->ampc[%d]", i);
	  MPS_DEBUG_MPC (s, 20, sec->bmpc[i], "sec->bmpc[%d]", i);
	}
      break;
    default:
      break;
    }
}

void
mps_secular_restart (mps_status * s)
{
  MPS_DEBUG_THIS_CALL;

  int i;

  switch (s->lastphase)
    {
    case float_phase:
      for (i = 0; i < s->n; i++)
	mpc_set_cplx (s->mroot[i], s->froot[i]);
      break;
    case dpe_phase:
      for (i = 0; i < s->n; i++)
	mpc_set_cdpe (s->mroot[i], s->droot[i]);
      break;
    default:
      break;
    }

  mps_mrestart (s);

  for (i = 0; i < s->n; i++)
    {
      mpc_get_cplx (s->froot[i], s->mroot[i]);
      mpc_get_cdpe (s->droot[i], s->mroot[i]);
    }
}

/**
 * @brief Deflate a secular equation lowering the degree of the
 * polynomial that represent it, if that is possible.
 *
 * Please note the <code>s->n</code> and <code>s->deg</code>
 * will not be touched by this routine, so you should check
 * that they are set according to <code>sec->n</code> if deflation
 * take place.
 * 
 * @see <code>mps_status_set_degree ()</code>
 *
 * @param s The mps_status of the computation
 * @param sec The secular equation that will be deflated.
 */
void
mps_secular_deflate (mps_status * s, mps_secular_equation * sec)
{
  int i, j, k;

  for (i = 0; i < sec->n; i++)
    {
      for (j = i + 1; j < sec->n; j++)
        {
          /* If the input is floating point check on the
           * DPE input */
          if (MPS_INPUT_CONFIG_IS_FP (s->input_config))
            {
              /* Do not deflate in floating point, since it is not working
               * correctly right now */
              MPS_DEBUG_WITH_INFO (s, "Floating point deflation still has some rough edges, so it's disabled");
              return;
              if (mpc_eq
                  (sec->initial_bmpc[i], sec->initial_bmpc[j],
                   s->mpwp_max) == 0)
                {
                  MPS_DEBUG_MPC (s, 10, sec->initial_bmpc[i],
                                 "sec->initial_bmpc[%d]", i);
                  MPS_DEBUG_MPC (s, 10, sec->initial_bmpc[j],
                                 "sec->initial_bmpc[%d]", j);
                  mpc_add_eq (sec->initial_ampc[i], sec->initial_ampc[j]);

                  for (k = j; k < sec->n - 1; k++)
                    {
                      mpc_set (sec->initial_ampc[k],
                               sec->initial_ampc[k + 1]);
                      mpc_set (sec->initial_bmpc[k],
                               sec->initial_bmpc[k + 1]);
                    }

                  sec->n--;
                  j--;
                  i--;
                }
            }
          /* Otherwise, in the case of rational or integer input
           * (that are handled in the same way) use initial_*mqpc
           * values */
          else if (MPS_INPUT_CONFIG_IS_INTEGER (s->input_config) ||
                   MPS_INPUT_CONFIG_IS_RATIONAL (s->input_config))
            {
              if (mpq_equal (sec->initial_bmpqrc[i], sec->initial_bmpqrc[j])
                  && mpq_equal (sec->initial_bmpqic[i],
                                sec->initial_bmpqic[j]))
                {
                  MPS_DEBUG_WITH_INFO (s,
                                       "Coefficients b[%d] and b[%d] are equal, deflating",
                                       i, j);
                  mpq_add (sec->initial_ampqrc[i], sec->initial_ampqrc[i],
                           sec->initial_ampqrc[j]);
                  mpq_add (sec->initial_ampqic[i], sec->initial_ampqic[i],
                           sec->initial_ampqic[j]);


                  /* Copy other coefficients back of one position */
                  for (k = j; k < sec->n - 1; k++)
                    {
                      mpq_set (sec->initial_ampqrc[k],
                               sec->initial_ampqrc[k + 1]);
                      mpq_set (sec->initial_ampqic[k],
                               sec->initial_ampqic[k + 1]);
                    }

                  sec->n--;
                  j--;
                }
            }
        }
    }

  /* If the input was rational or integer, we need to reset the dpe coefficients
   * according to it */
  if (MPS_INPUT_CONFIG_IS_INTEGER (s->input_config) ||
      MPS_INPUT_CONFIG_IS_RATIONAL (s->input_config))
    {
      mpf_t ftmp;
      mpf_init (ftmp);

      /* Set DPE coefficients */
      for (i = 0; i < sec->n; i++)
        {
          mpf_set_q (ftmp, sec->initial_ampqrc[i]);
          mpf_get_rdpe (cdpe_Re (sec->adpc[i]), ftmp);

          mpf_set_q (ftmp, sec->initial_ampqic[i]);
          mpf_get_rdpe (cdpe_Im (sec->adpc[i]), ftmp);

          mpf_set_q (ftmp, sec->initial_bmpqrc[i]);
          mpf_get_rdpe (cdpe_Re (sec->bdpc[i]), ftmp);

          mpf_set_q (ftmp, sec->initial_bmpqic[i]);
          mpf_get_rdpe (cdpe_Im (sec->bdpc[i]), ftmp);
        }

      mpf_clear (ftmp);
    }

  /* If the input was floating point update the coefficients using initial_*mpc
   * values */
  if (MPS_INPUT_CONFIG_IS_FP (s->input_config))
    {
      for (i = 0; i < sec->n; i++)
        {
          mpc_get_cdpe (sec->adpc[i], sec->ampc[i]);
          mpc_get_cdpe (sec->bdpc[i], sec->bmpc[i]);
        }

    }

  MPS_DEBUG (s, "Secular equation deflated to degree %lu", sec->n);
}

/**
 * @brief Raw version of mps_secular_equation_new that only
 * allocate space for the coefficients but relies on the user
 * to fill their values.
 *
 * @param s The mps_status of the computation.
 * @param n The degree of the new secular equation to be created.
 */
mps_secular_equation *
mps_secular_equation_new_raw (mps_status * s, unsigned long int n)
{
  int i;
  mps_secular_equation *sec =
    (mps_secular_equation *) mps_malloc (sizeof (mps_secular_equation));

  /* Allocate floating point coefficients */
  sec->afpc = cplx_valloc (n);
  sec->bfpc = cplx_valloc (n);

  /* Allocate complex dpe coefficients of the secular equation */
  sec->adpc = cdpe_valloc (n);
  sec->bdpc = cdpe_valloc (n);

  /* Allocate multiprecision complex coefficients of the secular equation */
  sec->ampc = mpc_valloc (n);
  sec->bmpc = mpc_valloc (n);
  sec->initial_ampc = mpc_valloc (n);
  sec->initial_bmpc = mpc_valloc (n);
  sec->initial_ampqrc = mpq_valloc (n);
  sec->initial_bmpqrc = mpq_valloc (n);
  sec->initial_ampqic = mpq_valloc (n);
  sec->initial_bmpqic = mpq_valloc (n);

  /* Allocate space for the moduli of the coefficients */
  sec->aadpc = rdpe_valloc (n);
  sec->abdpc = rdpe_valloc (n);
  sec->aafpc = double_valloc (n);
  sec->abfpc = double_valloc (n);

  /* Init multiprecision arrays */
  mpc_vinit2 (sec->ampc, n, s->mpwp);
  mpc_vinit2 (sec->bmpc, n, s->mpwp);
  mpc_vinit2 (sec->initial_ampc, n, s->mpwp);
  mpc_vinit2 (sec->initial_bmpc, n, s->mpwp);
  mpq_vinit (sec->initial_ampqrc, n);
  mpq_vinit (sec->initial_bmpqrc, n);
  mpq_vinit (sec->initial_ampqic, n);
  mpq_vinit (sec->initial_bmpqic, n);

  /* Epsilon arrays */
  sec->dregeneration_epsilon = rdpe_valloc (n);
  sec->fregeneration_epsilon = double_valloc (n);

  /* Set the epsilon array to zero */
  for (i = 0; i < n; i++)
    {
      sec->fregeneration_epsilon[i] = 0.0f;
      rdpe_set (sec->dregeneration_epsilon[i], rdpe_zero);
    }

  sec->n = n;

  /* Set up the mutexes for thread safety */
  sec->ampc_mutex = mps_newv (pthread_mutex_t, sec->n);
  sec->bmpc_mutex = mps_newv (pthread_mutex_t, sec->n);

  for (i = 0; i < n; i++)
    {
      pthread_mutex_init (&sec->ampc_mutex[i], NULL);
      pthread_mutex_init (&sec->bmpc_mutex[i], NULL);
    }
  
  return sec;
}

/**
 * @brief Create a new secular equation struct
 * 
 * @param s The mps_status of the computation.
 * @param afpc The floating point complex numerator coefficients.
 * @param bfpc The floating point complex denominator coefficients.
 * @param n The degree of the secular equation.
 */
mps_secular_equation *
mps_secular_equation_new (mps_status * s, cplx_t * afpc, cplx_t * bfpc,
                          unsigned long int n)
{
  int i;

  /* Allocate the space for the new struct */
  mps_secular_equation *sec = mps_secular_equation_new_raw (s, n);

  /* Copy the complex coefficients passed as argument */
  for (i = 0; i < n; i++)
    {
      /* a_i coefficients */
      cplx_set (sec->afpc[i], afpc[i]);

      /* b_i coefficients */
      cplx_set (sec->bfpc[i], bfpc[i]);
    }

  sec->n = n;
  mps_secular_deflate (s, sec);

  for (i = 0; i < sec->n; i++)
    {
      cdpe_init (sec->adpc[i]);
      cdpe_set_x (sec->adpc[i], sec->afpc[i]);

      mpc_set_cplx (sec->ampc[i], sec->afpc[i]);

      cdpe_init (sec->bdpc[i]);
      cdpe_set_x (sec->bdpc[i], sec->bfpc[i]);

      mpc_set_cplx (sec->bmpc[i], sec->bfpc[i]);
    }

  return sec;
}

/**
 * @brief Free a secular equation and the data in it.
 *
 * @param s The secular equation to be freed.
 */
void
mps_secular_equation_free (mps_secular_equation * s)
{
  /* Free internal data */
  cplx_vfree (s->afpc);
  cplx_vfree (s->bfpc);

  cdpe_vfree (s->adpc);
  cdpe_vfree (s->bdpc);

  mpc_vclear (s->ampc, s->n);
  mpc_vclear (s->bmpc, s->n);

  mpc_vfree (s->ampc);
  mpc_vfree (s->bmpc);

  rdpe_vfree (s->aadpc);
  rdpe_vfree (s->abdpc);
  double_vfree (s->aafpc);
  double_vfree (s->abfpc);

  /* And old coefficients */
  mpc_vclear (s->initial_ampc, s->n);
  mpc_vclear (s->initial_bmpc, s->n);
  mpq_vclear (s->initial_ampqrc, s->n);
  mpq_vclear (s->initial_bmpqrc, s->n);
  mpq_vclear (s->initial_ampqic, s->n);
  mpq_vclear (s->initial_bmpqic, s->n);
  mpc_vfree (s->initial_ampc);
  mpc_vfree (s->initial_bmpc);
  mpq_vfree (s->initial_ampqrc);
  mpq_vfree (s->initial_bmpqrc);
  mpq_vfree (s->initial_ampqic);
  mpq_vfree (s->initial_bmpqic);

  /* Epsilon arrays */
  rdpe_vfree (s->dregeneration_epsilon);
  free (s->fregeneration_epsilon);

  /* Mutexes */
  free (s->ampc_mutex);
  free (s->bmpc_mutex);

  /* ...and then release it */
  free (s);
}


/**
 * @brief Evaluate secular equation in the point x.
 *
 * The evalutation will be done in floating point and this routine
 * is used only for debugging purpose.
 * 
 * @param s The mps_status of the computation.
 * @param x The point in which the secular equation will be evaluated.
 * @param sec_ev The result of the evalutation (output variable). 
 */
void
mps_secular_evaluate (mps_status * s, cplx_t x, cplx_t sec_ev)
{
  cplx_t ctmp;
  int i;
  mps_secular_equation *sec = (mps_secular_equation *) s->secular_equation;
  cplx_set (sec_ev, cplx_zero);

  for (i = 0; i < s->n; i++)
    {
      /* Compute 1 / (x - b_i) */
      cplx_sub (ctmp, x, sec->bfpc[i]);
      cplx_inv_eq (ctmp);

      /* Compute a_i / (x - b_i) */
      cplx_mul_eq (ctmp, sec->afpc[i]);

      /* Sum to the secular eqation */
      cplx_add_eq (sec_ev, ctmp);
    }

  cplx_sub_eq (sec_ev, cplx_one);
}

/**
 * @brief Secular version of <code>mps_check_data ()</code> that
 * does nothing except to set the starting case according to
 * <code>sec->starting_case</code>.
 */
void
mps_secular_check_data (mps_status * s, char *which_case)
{
  /* While we can't found a good criterion to check
   * the possibility to start in pure floating point we
   * use the DPE version. */
  *which_case = (s->input_config->starting_phase == float_phase) ? 'f' : 'd';
}

/**
 * @brief Raise precision of the coefficient of the secular equation
 * (not the roots and neither the precison of the system) to <code>wp</code>.
 *
 * @param s The mps_status of the computation.
 * @param wp The bits of precision to which the coefficients will be set.
 *
 * @see <code>mps_secular_raise_root_precision ()</code>
 * @see <code>mps_secular_raise_precision ()</code>
 */
void
mps_secular_raise_coefficient_precision (mps_status * s, int wp)
{
  MPS_DEBUG_THIS_CALL;

  int i;
  mps_secular_equation *sec = (mps_secular_equation *) s->secular_equation;

  for (i = 0; i < s->n; i++)
    {
      mpc_set_prec (sec->ampc[i], wp);
      mpc_set_prec (sec->bmpc[i], wp);
    }

  if (MPS_INPUT_CONFIG_IS_MONOMIAL (s->input_config))
    mps_monomial_poly_raise_precision (s, s->monomial_poly, wp);

  if (MPS_INPUT_CONFIG_IS_SECULAR (s->input_config) && (!MPS_INPUT_CONFIG_IS_FP (s->input_config)))
    {
      for (i = 0; i < s->n; i++)
	{
	  mpf_set_q (mpc_Re (sec->ampc[i]), sec->initial_ampqrc[i]);
	  mpf_set_q (mpc_Im (sec->ampc[i]), sec->initial_ampqic[i]);
	  mpf_set_q (mpc_Re (sec->bmpc[i]), sec->initial_bmpqrc[i]);
	  mpf_set_q (mpc_Im (sec->bmpc[i]), sec->initial_bmpqic[i]);
	}
    }

  rdpe_set_2dl (s->mp_epsilon, 1.0, -wp);
  MPS_DEBUG_WITH_INFO (s, "Precision of the coefficients is now at %d bits", wp);
  if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)
    MPS_DEBUG_RDPE (s, s->mp_epsilon, "Machine epsilon is s->mp_epsilon");
}

/**
 * @brief Raise precision of the roots (not the coefficients nor the
 * system) to <code>wp</code> bits.
 * 
 * @param s The mps_status of the computation.
 * @param wp The bits of precision to which the roots will be set.
 *
 * @see <code>mps_secular_raise_coefficient_precision ()</code>
 * @see <code>mps_secular_raise_precision ()</code>
 */
void
mps_secular_raise_root_precision (mps_status * s, int wp)
{
  MPS_DEBUG_THIS_CALL;
  int i;

  for (i = 0; i < s->n; i++)
    {
      mpc_set_prec (s->mroot[i], wp);
    }
}

/**
 * @brief Raise (or lower) the precision of the coefficients to
 * <code>wp</code> bits. This will change precision of the coefficients
 * via <code>mps_secular_raise_coefficient_precision ()</code>,
 * of the roots via <code>mps_secular_raise_root_precision ()</code>
 * and set <code>s->mpwp</code>.
 *
 * @param s The mps_status of the computation.
 * @param wp The bits of precision to which all the computation will
 * be brought.
 */
void
mps_secular_raise_precision (mps_status * s, int wp)
{
  MPS_DEBUG_THIS_CALL;

  mps_secular_raise_coefficient_precision (s, wp);
  mps_secular_raise_root_precision (s, wp);
  s->mpwp = wp;

  s->just_raised_precision = true;
}

/**
 * @brief Prepare data for the iteration in the new phase specified
 * in the second parameter.
 *
 * Note that for now this function is only able to handle switch
 * from floating point phases (i.e. float_phase or dpe_phase) to
 * multiprecision, and not coming back.
 *
 * @param s The mps_status of the computation.
 * @param phase The phase to switch the computation to.
 */
void
mps_secular_switch_phase (mps_status * s, mps_phase phase)
{
  MPS_DEBUG_THIS_CALL;

  s->just_raised_precision = true;

  int i = 0;
  mps_secular_equation *sec = (mps_secular_equation *) s->secular_equation;
  if (phase == mp_phase)
    {
      /* Debug the approximations that we have now before going
       * to the multiprecision phase */
      if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)
	{
	  MPS_DEBUG (s, "Dumping current approximations before starting MP");
	  mps_dump (s);
	}

      mps_secular_raise_precision (s, MPS_SECULAR_STARTING_MP_PRECISION);
      switch (s->lastphase)
        {
        case float_phase:
          /* Copy the approximated roots and the
           * secular equation coefficients */
          for (i = 0; i < s->n; i++)
            {
              mpc_set_cplx (s->mroot[i], s->froot[i]);
              mpc_set_cplx (sec->ampc[i], sec->afpc[i]);
              mpc_set_cplx (sec->bmpc[i], sec->bfpc[i]);
              rdpe_set_d (s->drad[i], s->frad[i]);
            }
          break;

        case dpe_phase:
          /* Copy the coefficients and the approximated
           * roots into the multiprecision values    */
          for (i = 0; i < s->n; i++)
            {
              mpc_set_cdpe (s->mroot[i], s->droot[i]);
              mpc_set_cdpe (sec->ampc[i], sec->adpc[i]);
              mpc_set_cdpe (sec->bmpc[i], sec->bdpc[i]);
            }

        default:
          break;

        }

      /* Set lastphase to mp_phase */
      s->lastphase = mp_phase;

      /* Set epsilon */
      rdpe_set_2dl (s->mp_epsilon, 1.0, -s->mpwp + 1);
    }
  else
    {
      fprintf (stderr, "mps_secular_switch_phase is only able to manage\n"
               "switches from float_phase or dpe_phase to mp_phase. Aborting.");
      exit (EXIT_FAILURE);
    }
}

/**
 * @brief Update radii of the roots according to the coefficients
 * of the secular equation in this moment, if they are better of
 * the radii present now.
 *
 * @param s The mps_status of the computation.
 */
void
mps_secular_set_radii (mps_status * s)
{
  MPS_DEBUG_THIS_CALL;

  int i;
  mps_secular_equation *sec = (mps_secular_equation *) s->secular_equation;

  /* DPE and multiprecision implementation */
  rdpe_t rad, rad_eps;
  cdpe_t ctmp;
  rdpe_t * drad = rdpe_valloc (s->n);

  mpc_t mtmp;
  mpc_init2 (mtmp, mps_status_get_data_prec_max (s));

  /* Check if the Gerschgorin's radii are more convenient */
  for (i = 0; i < s->n; i++)
    {
      if (s->lastphase == mp_phase)
	rdpe_set (rad_eps, s->mp_epsilon);
      else 
	rdpe_set_d (rad_eps, DBL_EPSILON);

      rdpe_mul_eq_d (rad_eps, s->n * 4);
      rdpe_add_eq (rad_eps, rdpe_one);
      
      mpc_get_cdpe (ctmp, sec->ampc[i]);
      cdpe_mod (rad, ctmp);
      
      rdpe_mul_eq (rad, rad_eps);
      
      rdpe_mul_eq_d (rad, (double) s->n);
      
      rdpe_set (drad[i], rad);
    }

  switch (s->lastphase)
    {
    case float_phase:
       { 
	 for (i = 0; i < s->n; i++) 
	   { 
	     rdpe_set_d (s->drad[i], s->frad[i]); 
	     mpc_set_d  (s->mroot[i], cplx_Re (s->froot[i]),  
			 cplx_Im (s->froot[i])); 
	   } 
	 
	 mps_mcluster (s, drad, 2.0 * s->n);  
	 mps_fmodify (s, false);  

	 for (i = 0; i < s->n; i++)
	   s->frad[i] = rdpe_get_d (s->drad[i]); 
       }
       break;

    case dpe_phase:
      mps_mcluster (s, drad, 2.0 * s->n);
      mps_dmodify (s, false);
      break;
      
    case mp_phase:
      mps_mcluster (s, drad, 2.0 * s->n);
      mps_mmodify (s, false);
      break;

    default:
      break;
    }
      
  rdpe_vfree (drad);

  mpc_clear (mtmp);
}
