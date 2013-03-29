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

/**
 * @brief Utility routine that dumps the DPE coefficients of the secular equation
 * passed as second argument. It is only used for debugging purpose.
 *
 * @param s The mps_context of the computation
 * @param sec The secular equation whose DPE coefficients will be dumped
 */
void
mps_secular_dump (mps_context * s, mps_secular_equation * sec)
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
      for (i = 0; i < MPS_POLYNOMIAL (sec)->degree; i++)
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
mps_secular_restart (mps_context * s)
{
  MPS_DEBUG_THIS_CALL;

  int i;

  switch (s->lastphase)
    {
    case float_phase:
      for (i = 0; i < s->n; i++)
        mpc_set_cplx (s->root[i]->mvalue, s->root[i]->fvalue);
      break;
    case dpe_phase:
      for (i = 0; i < s->n; i++)
        mpc_set_cdpe (s->root[i]->mvalue, s->root[i]->dvalue);
      break;
    default:
      break;
    }

  mps_mrestart (s);

  for (i = 0; i < s->n; i++)
    {
      mpc_get_cplx (s->root[i]->fvalue, s->root[i]->mvalue);
      mpc_get_cdpe (s->root[i]->dvalue, s->root[i]->mvalue);
    }
}

/**
 * @brief Deflate a secular equation lowering the degree of the
 * polynomial that represent it, if that is possible.
 *
 * Please note the <code>s->n</code> and <code>s->deg</code>
 * will not be touched by this routine, so you should check
 * that they are set according to <code>MPS_POLYNOMIAL (sec)->degree</code> if deflation
 * takes place.
 * 
 * @see <code>mps_context_set_degree ()</code>
 *
 * @param s The mps_context of the computation
 * @param sec The secular equation that will be deflated.
 */
void
mps_secular_deflate (mps_context * s, mps_secular_equation * sec)
{
  int i, j, k;

  /* If the input is floating point check on the
   * DPE input */
  if (MPS_STRUCTURE_IS_FP (MPS_POLYNOMIAL (sec)->structure))
    {
      /* Do not deflate in floating point, since it is not working
       * correctly right now */
      MPS_DEBUG_WITH_INFO (s, "Floating point deflation still has some rough edges, so it's disabled");
      return;
    }

  for (i = 0; i < MPS_POLYNOMIAL (sec)->degree; i++)
    {
      for (j = i + 1; j < MPS_POLYNOMIAL (sec)->degree; j++)
        {
          /* Otherwise, in the case of rational or integer input
           * (that are handled in the same way) use initial_*mqpc
           * values */
          if (MPS_STRUCTURE_IS_INTEGER (MPS_POLYNOMIAL (sec)->structure) ||
              MPS_STRUCTURE_IS_RATIONAL (MPS_POLYNOMIAL (sec)->structure))
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
                  for (k = j; k < MPS_POLYNOMIAL (sec)->degree - 1; k++)
                    {
                      mpq_set (sec->initial_ampqrc[k],
                               sec->initial_ampqrc[k + 1]);
                      mpq_set (sec->initial_ampqic[k],
                               sec->initial_ampqic[k + 1]);
                    }

                  MPS_POLYNOMIAL (sec)->degree--;
                  j--;
                }
            }
        }
    }

  /* If the input was rational or integer, we need to reset the dpe coefficients
   * according to it */
  if (MPS_STRUCTURE_IS_INTEGER (MPS_POLYNOMIAL (sec)->structure) ||
      MPS_STRUCTURE_IS_RATIONAL (MPS_POLYNOMIAL (sec)->structure))
    {
      mpf_t ftmp;
      mpf_init (ftmp);

      /* Set DPE coefficients */
      for (i = 0; i < MPS_POLYNOMIAL (sec)->degree; i++)
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
  if (MPS_STRUCTURE_IS_FP (MPS_POLYNOMIAL (sec)->structure))
    {
      for (i = 0; i < MPS_POLYNOMIAL (sec)->degree; i++)
        {
          mpc_get_cdpe (sec->adpc[i], sec->ampc[i]);
          mpc_get_cdpe (sec->bdpc[i], sec->bmpc[i]);
        }

    }

  MPS_DEBUG (s, "Secular equation deflated to degree %d", MPS_POLYNOMIAL (sec)->degree);
}

/**
 * @brief Raw version of mps_secular_equation_new that only
 * allocate space for the coefficients but relies on the user
 * to fill their values.
 *
 * @param s The mps_context of the computation.
 * @param n The degree of the new secular equation to be created.
 */
mps_secular_equation *
mps_secular_equation_new_raw (mps_context * s, unsigned long int n)
{
  int i;
  mps_secular_equation *sec =
    (mps_secular_equation *) mps_malloc (sizeof (mps_secular_equation));

  mps_polynomial_init (s, MPS_POLYNOMIAL (sec));

  MPS_POLYNOMIAL (sec)->type_name = "mps_secular_equation";

  /* Hook up the overloaded methods for secular equations */
  mps_polynomial * p = MPS_POLYNOMIAL (sec);
  p->feval = mps_secular_poly_feval_with_error;
  p->deval = mps_secular_poly_deval_with_error;
  p->meval = mps_secular_poly_meval_with_error;
  p->fstart = mps_secular_poly_fstart;
  p->dstart = mps_secular_poly_dstart;
  p->mstart = mps_secular_poly_mstart;
  p->free = mps_secular_equation_free;
  p->raise_data = mps_secular_raise_coefficient_precision;
  p->fnewton = mps_secular_fnewton;
  p->dnewton = mps_secular_dnewton;
  p->mnewton = mps_secular_mnewton;

  /* Allocate floating point coefficients */
  sec->afpc = cplx_valloc (n);
  sec->bfpc = cplx_valloc (n);

  /* Allocate complex dpe coefficients of the secular equation */
  sec->adpc = cdpe_valloc (n);
  sec->bdpc = cdpe_valloc (n);

  /* Allocate multiprecision complex coefficients of the secular equation */
  /* Prepare the double buffer */
  struct mps_secular_equation_double_buffer *db = &sec->db;
  db->ampc1 = mpc_valloc (n);
  db->ampc2 = mpc_valloc (n);
  db->bmpc1 = mpc_valloc (n);
  db->bmpc2 = mpc_valloc (n);

  sec->ampc = db->ampc1;
  sec->bmpc = db->bmpc1;

  db->active = 1;

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
  mpc_vinit2 (sec->db.ampc1, n, s->mpwp);
  mpc_vinit2 (sec->db.bmpc1, n, s->mpwp);
  mpc_vinit2 (sec->db.ampc2, n, s->mpwp);
  mpc_vinit2 (sec->db.bmpc2, n, s->mpwp);
  mpc_vinit2 (sec->initial_ampc, n, s->mpwp);
  mpc_vinit2 (sec->initial_bmpc, n, s->mpwp);
  mpq_vinit (sec->initial_ampqrc, n);
  mpq_vinit (sec->initial_bmpqrc, n);
  mpq_vinit (sec->initial_ampqic, n);
  mpq_vinit (sec->initial_bmpqic, n);

  MPS_POLYNOMIAL (sec)->degree = n;

  /* Set up the mutexes for thread safety */
  sec->ampc_mutex = mps_newv (pthread_mutex_t, MPS_POLYNOMIAL (sec)->degree);
  sec->bmpc_mutex = mps_newv (pthread_mutex_t, MPS_POLYNOMIAL (sec)->degree);

  for (i = 0; i < n; i++)
    {
      pthread_mutex_init (&sec->ampc_mutex[i], NULL);
      pthread_mutex_init (&sec->bmpc_mutex[i], NULL);
    }

  pthread_mutex_init (&sec->precision_mutex, NULL);
  
  return sec;
}

/**
 * @brief Create a new secular equation struct
 * 
 * @param s The mps_context of the computation.
 * @param afpc The floating point complex numerator coefficients.
 * @param bfpc The floating point complex denominator coefficients.
 * @param n The degree of the secular equation.
 */
mps_secular_equation *
mps_secular_equation_new (mps_context * s, cplx_t * afpc, cplx_t * bfpc,
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

  MPS_POLYNOMIAL (sec)->degree = n;
  mps_secular_deflate (s, sec);

  for (i = 0; i < MPS_POLYNOMIAL (sec)->degree; i++)
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
 * @param ctx The secular equation to be freed.
 * @param p The secular equation casted to a mps_polynomial
 */
void
mps_secular_equation_free (mps_context *ctx, mps_polynomial * p)
{
  mps_secular_equation *s = MPS_SECULAR_EQUATION (p);

  /* Free internal data */
  cplx_vfree (s->afpc);
  cplx_vfree (s->bfpc);

  cdpe_vfree (s->adpc);
  cdpe_vfree (s->bdpc);

  mpc_vclear (s->db.ampc1, MPS_POLYNOMIAL (s)->degree);
  mpc_vclear (s->db.bmpc1, MPS_POLYNOMIAL (s)->degree);
  mpc_vclear (s->db.ampc2, MPS_POLYNOMIAL (s)->degree);
  mpc_vclear (s->db.bmpc2, MPS_POLYNOMIAL (s)->degree);

  mpc_vfree (s->db.ampc1);
  mpc_vfree (s->db.ampc2);
  mpc_vfree (s->db.bmpc1);
  mpc_vfree (s->db.bmpc2);

  rdpe_vfree (s->aadpc);
  rdpe_vfree (s->abdpc);
  double_vfree (s->aafpc);
  double_vfree (s->abfpc);

  /* And old coefficients */
  mpc_vclear (s->initial_ampc, MPS_POLYNOMIAL (s)->degree);
  mpc_vclear (s->initial_bmpc, MPS_POLYNOMIAL (s)->degree);
  mpq_vclear (s->initial_ampqrc, MPS_POLYNOMIAL (s)->degree);
  mpq_vclear (s->initial_bmpqrc, MPS_POLYNOMIAL (s)->degree);
  mpq_vclear (s->initial_ampqic, MPS_POLYNOMIAL (s)->degree);
  mpq_vclear (s->initial_bmpqic, MPS_POLYNOMIAL (s)->degree);
  mpc_vfree (s->initial_ampc);
  mpc_vfree (s->initial_bmpc);
  mpq_vfree (s->initial_ampqrc);
  mpq_vfree (s->initial_bmpqrc);
  mpq_vfree (s->initial_ampqic);
  mpq_vfree (s->initial_bmpqic);

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
 * @param s The mps_context of the computation.
 * @param x The point in which the secular equation will be evaluated.
 * @param sec_ev The result of the evalutation (output variable). 
 */
void
mps_secular_evaluate (mps_context * s, cplx_t x, cplx_t sec_ev)
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
 *
 * @param s The current mps_context.
 * @param which_case Pointer to a char that will be set to 'f' or 'd' depending
*  on the chosen start phase. 
 */
void
mps_secular_check_data (mps_context * s, char *which_case)
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
 * @param s The mps_context of the computation.
 * @param p The secular equation casted to a mps_polynomial.
 * @param wp The bits of precision to which the coefficients will be set.
 *
 * @see <code>mps_secular_raise_root_precision ()</code>
 * @see <code>mps_secular_raise_precision ()</code>
 */
long int
mps_secular_raise_coefficient_precision (mps_context * s, mps_polynomial * p, long int wp)
{
  MPS_DEBUG_THIS_CALL;

  int i;
  mps_secular_equation *sec = MPS_SECULAR_EQUATION (p);

  mpc_t * raising_ampc;
  mpc_t * raising_bmpc;

  pthread_mutex_lock (&sec->precision_mutex);

  if (wp < mpc_get_prec (sec->ampc[0]))
    {
      pthread_mutex_unlock (&sec->precision_mutex);
      return mpc_get_prec (sec->ampc[0]);
    }

  if (sec->db.active == 1)
    {
      raising_ampc = sec->db.ampc2;
      raising_bmpc = sec->db.bmpc2;
    }
  else
    {
      raising_ampc = sec->db.ampc1;
      raising_bmpc = sec->db.bmpc1;
    }

  for (i = 0; i < s->n; i++)
    {
      mpc_set_prec (raising_ampc[i], wp);
      if (!MPS_STRUCTURE_IS_FP (s->active_poly->structure))
        {
          mpf_set_q (mpc_Re (raising_ampc[i]), sec->initial_ampqrc[i]);
          mpf_set_q (mpc_Im (raising_ampc[i]), sec->initial_ampqic[i]);
        }
      else
        mpc_set (raising_ampc[i], sec->ampc[i]);

      mpc_set_prec (raising_bmpc[i], wp);
      if (!MPS_STRUCTURE_IS_FP (s->active_poly->structure))
        {
          mpf_set_q (mpc_Re (raising_bmpc[i]), sec->initial_bmpqrc[i]);
          mpf_set_q (mpc_Im (raising_bmpc[i]), sec->initial_bmpqic[i]);
        }
      else
        mpc_set (raising_bmpc[i], sec->bmpc[i]);
    }

  sec->ampc = raising_ampc;
  sec->bmpc = raising_bmpc;

  sec->db.active = (sec->db.active % 2) + 1;

  pthread_mutex_unlock (&sec->precision_mutex);

  if (s->debug_level & MPS_DEBUG_MEMORY)
    MPS_DEBUG_WITH_INFO (s, "Precision of the coefficients is now at %ld bits", wp);

  return mpc_get_prec (sec->ampc[0]);
}

/**
 * @brief Raise precision of the roots (not the coefficients nor the
 * system) to <code>wp</code> bits.
 * 
 * @param s The mps_context of the computation.
 * @param wp The bits of precision to which the roots will be set.
 *
 * @see <code>mps_secular_raise_coefficient_precision ()</code>
 * @see <code>mps_secular_raise_precision ()</code>
 */
void
mps_secular_raise_root_precision (mps_context * s, int wp)
{
  MPS_DEBUG_THIS_CALL;
  int i;

  for (i = 0; i < s->n; i++)
    {
      mpc_set_prec (s->root[i]->mvalue, wp);
    }
}

/**
 * @brief Raise (or lower) the precision of the coefficients to
 * <code>wp</code> bits. This will change precision of the coefficients
 * via <code>mps_secular_raise_coefficient_precision ()</code>,
 * of the roots via <code>mps_secular_raise_root_precision ()</code>
 * and set <code>s->mpwp</code>.
 *
 * @param s The mps_context of the computation.
 * @param wp The bits of precision to which all the computation will
 * be brought.
 */
void
mps_secular_raise_precision (mps_context * s, int wp)
{
  MPS_DEBUG_THIS_CALL;

  int i;

  mps_secular_raise_coefficient_precision (s, MPS_POLYNOMIAL (s->secular_equation), wp);
  mps_secular_raise_root_precision (s, wp);

  s->mpwp = wp;
  rdpe_set_2dl (s->mp_epsilon, 1.0, -wp);

  s->just_raised_precision = true;


  for (i = 0; i < s->n; i++)
    {
      s->root[i]->approximated = false;
      s->root[i]->again = true;
    }
}

/**
 * @brief Prepare data for the iteration in the new phase specified
 * in the second parameter.
 *
 * Note that for now this function is only able to handle switch
 * from floating point phases (i.e. float_phase or dpe_phase) to
 * multiprecision, and not coming back.
 *
 * @param s The mps_context of the computation.
 * @param phase The phase to switch the computation to.
 */
void
mps_secular_switch_phase (mps_context * s, mps_phase phase)
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
              mpc_set_cplx (s->root[i]->mvalue, s->root[i]->fvalue);
              mpc_set_cplx (sec->ampc[i], sec->afpc[i]);
              mpc_set_cplx (sec->bmpc[i], sec->bfpc[i]);
              rdpe_set_d (s->root[i]->drad, s->root[i]->frad);
            }
          break;

        case dpe_phase:
          /* Copy the coefficients and the approximated
           * roots into the multiprecision values    */
          for (i = 0; i < s->n; i++)
            {
              mpc_set_cdpe (s->root[i]->mvalue, s->root[i]->dvalue);
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
 * @param s The mps_context of the computation.
 */
void
mps_secular_set_radii (mps_context * s)
{
  MPS_DEBUG_THIS_CALL;

  int i;
  mps_secular_equation *sec = (mps_secular_equation *) s->secular_equation;

  /* DPE and multiprecision implementation */
  rdpe_t rad, rad_eps, rtmp;
  cdpe_t ctmp;
  rdpe_t * drad = rdpe_valloc (s->n);

  mpc_t mtmp;
  mpc_init2 (mtmp, mps_context_get_data_prec_max (s));

  if (s->lastphase == mp_phase)
    rdpe_set (rad_eps, s->mp_epsilon);
  else 
    rdpe_set_d (rad_eps, DBL_EPSILON);

  rdpe_mul_eq_d (rad_eps, s->n * 4);
  rdpe_add_eq (rad_eps, rdpe_one);

  /* Check if the Gerschgorin's radii are more convenient */
  for (i = 0; i < s->n; i++)
    {
      mpc_get_cdpe (ctmp, sec->ampc[i]);
      cdpe_mod (rad, ctmp);

      rdpe_mul_eq (rad, rad_eps);

      rdpe_mul_eq_d (rad, (double) s->n);

      rdpe_set (drad[i], rad);

      mpc_rmod (rtmp, s->root[i]->mvalue);
      if (s->lastphase == mp_phase)
        rdpe_mul_eq (rtmp, s->mp_epsilon);
      else
        rdpe_mul_eq_d (rtmp, DBL_EPSILON);
      rdpe_mul_eq_d (rtmp, 4.0);
      rdpe_add_eq (drad[i], rtmp);
    }

  switch (s->lastphase)
    {
    case float_phase:
       { 
         for (i = 0; i < s->n; i++) 
           { 
             rdpe_set_d (s->root[i]->drad, s->root[i]->frad); 
             mpc_set_d  (s->root[i]->mvalue, cplx_Re (s->root[i]->fvalue),  
                         cplx_Im (s->root[i]->fvalue)); 
           } 
         
         mps_mcluster (s, drad, 2.0 * s->n);  
         mps_fmodify (s, false);  

         for (i = 0; i < s->n; i++)
           {
             s->root[i]->frad = rdpe_get_d (s->root[i]->drad); 
             if (s->root[i]->frad == 0.0)
               s->root[i]->frad += cplx_mod (s->root[i]->fvalue) * DBL_EPSILON;
           }
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

mps_boolean 
mps_secular_poly_feval_with_error (mps_context * ctx, mps_polynomial * p, cplx_t x, cplx_t value, double * error)
{
  cplx_t ctmp;
  int i;
  mps_secular_equation *sec = MPS_SECULAR_EQUATION (p);
  if (!mps_secular_feval_with_error (ctx, p, x, value, error))
    return false;

  *error /= cplx_mod (value);

  for (i = 0; i < p->degree; i++)
    {
      cplx_sub (ctmp, x, sec->bfpc[i]);
      cplx_mul_eq (value, ctmp);
    }

  cplx_mul_eq_d (value, -1.0);
  *error *= cplx_mod (value);
  return true;
}

mps_boolean 
mps_secular_poly_deval_with_error (mps_context * ctx, mps_polynomial * p, cdpe_t x, cdpe_t value, rdpe_t error)
{
  cdpe_t ctmp;
  rdpe_t rtmp;
  int i;
  mps_secular_equation *sec = MPS_SECULAR_EQUATION (p);
  if (!mps_secular_deval_with_error (ctx, p, x, value, error))
    return false;

  cdpe_mod (rtmp, value);
  rdpe_div_eq (error, rtmp);

  for (i = 0; i < p->degree; i++)
    {
      cdpe_sub (ctmp, x, sec->bdpc[i]);
      cdpe_mul_eq (value, ctmp);
    }

  cdpe_mul_eq_d (value, -1.0);

  cdpe_mod (rtmp, value);
  rdpe_mul_eq (error, rtmp);

  return true;
}

mps_boolean 
mps_secular_poly_meval_with_error (mps_context * ctx, mps_polynomial * p, mpc_t x, mpc_t value, rdpe_t error)
{
  mpc_t ctmp;
  rdpe_t rtmp;
  int i;
  mps_secular_equation *sec = MPS_SECULAR_EQUATION (p);

  if (!mps_secular_meval_with_error (ctx, p, x, value, error))
    {
      return false;
    }

  mpc_rmod (rtmp, value);
  rdpe_div_eq (error, rtmp);

  mpc_init2 (ctmp, mpc_get_prec (x));

  for (i = 0; i < p->degree; i++)
    {
      mpc_sub (ctmp, x, sec->bmpc[i]);
      mpc_mul_eq (value, ctmp);
    }

  mpc_set_si (ctmp, -1, 0);
  mpc_mul_eq (value, ctmp);

  mpc_clear (ctmp);

  mpc_rmod (rtmp, value);
  rdpe_mul_eq (error, rtmp);

  return true;
}

void 
mps_secular_poly_fstart (mps_context * ctx, mps_polynomial * p)
{
  mps_secular_fstart (ctx, MPS_SECULAR_EQUATION (p));
}

void 
mps_secular_poly_dstart (mps_context * ctx, mps_polynomial * p)
{
  mps_secular_dstart (ctx, MPS_SECULAR_EQUATION (p));
}

void 
mps_secular_poly_mstart (mps_context * ctx, mps_polynomial * p)
{
  mps_secular_mstart (ctx, MPS_SECULAR_EQUATION (p));
}
