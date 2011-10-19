/*
 * secular-ga.c
 *
 *  Created on: 15/giu/2011
 *      Author: leonardo
 */

#include <mps/debug.h>
#include <mps/core.h>
#include <mps/poly.h>
#include <mps/link.h>
#include <mps/secular.h>
#include <mps/debug.h>
#include <math.h>
#include <string.h>

/**
 * @file
 * @brief Implementation of secular routines for the GA algorithm
 */

/**
 * @brief Routine that performs a block of iteration
 * in floating point on the secular equation.
 *
 * @param s the pointer to the mps_status struct.
 * @param maxit Maximum number of iteration to perform.
 * @return The number of approximated roots after the iteration.
 */
int
mps_secular_ga_fiterate (mps_status * s, int maxit)
{
  MPS_DEBUG_THIS_CALL;

  int computed_roots = 0;
  int iterations = 0;
  int i;
  int nit = 0;
  int it_threshold;
  int old_cr;

#ifndef DISABLE_DEBUG
  clock_t *my_clock = mps_start_timer ();
#endif

  mps_secular_equation *sec = s->secular_equation;

  double old_rad;
  cplx_t old_root;

  /* Check the roots that are already isolated or approximated
   * and mark them as computed */
  for (i = 0; i < s->n; i++)
    if (!s->again[i])
      computed_roots++;

  /* Set the iterations threshold to 2 iterations
   * for every non approximated root. */
  it_threshold = 2 * (s->n - computed_roots);

  /* Save the number of already computed roots se we can check if this
   * packet has been an improvement on that */
  old_cr = computed_roots;
  if (s->debug_level & MPS_DEBUG_PACKETS)
    {
      MPS_DEBUG (s, "There are %d roots already approximated", old_cr);
    }

  while (computed_roots < s->n && iterations < maxit - 1)
    {
      cplx_t corr, abcorr;
      double modcorr;

      /* Increase iterations counter */
      iterations++;

      for (i = 0; i < s->n; i++)
        {
          if (s->again[i])
            {
              nit++;
              cplx_set (old_root, s->froot[i]);
              old_rad = s->frad[i];
              mps_secular_fnewton (s, s->froot[i], &s->frad[i], corr,
                                   &s->again[i], &i);

              /* Apply Aberth correction */
              mps_faberth (s, i, abcorr);
              cplx_mul_eq (abcorr, corr);
              cplx_sub (abcorr, cplx_one, abcorr);
              cplx_div (abcorr, corr, abcorr);
              cplx_sub_eq (s->froot[i], abcorr);

	      // MPS_DEBUG_CPLX (s, abcorr, "aberth_correction");
	      /* MPS_DEBUG (s, "Iteration %d", nit); */
	      /* MPS_DEBUG_CPLX (s, corr, "corr"); */
	      /* MPS_DEBUG_CPLX (s, abcorr, "abcorr"); */
	      /* mps_dump (s, s->logstr); */

              /* Check if we need to switch to DPE */
              if (isnan (cplx_Re (s->froot[i]))
                  || isinf (cplx_Re (s->froot[i]))
                  || isnan (cplx_Im (s->froot[i]))
                  || isinf (cplx_Im (s->froot[i])) || isnan (s->frad[i])
                  || isinf (s->frad[i])
		  || s->status[i][0] == 'x')
                {
		  if (s->status[i][0] != 'x')
		    {
		      MPS_DEBUG_WITH_INFO (s,
					   "Switching to DPE phase because NAN or INF was introduced in computation");
		    }
		  else
		    {
		      s->status[i][0] = 'c';
		      MPS_DEBUG_WITH_INFO (s, "Switching to DPE phase because there is an approximation not representable in double");
		    }
                  cplx_set (s->froot[i], old_root);
                  s->frad[i] = old_rad;
                  s->lastphase = dpe_phase;

                  /* Copy roots, radius and coefficients */
                  for (i = 0; i < s->n; i++)
                    {
                      cdpe_set_x (s->droot[i], s->froot[i]);
                      rdpe_set_d (s->drad[i], s->frad[i]);

		      cdpe_set_x (sec->adpc[i], sec->afpc[i]);
		      cdpe_set_x (sec->bdpc[i], sec->bfpc[i]);

		      MPS_DEBUG_CDPE (s, sec->adpc[i], "sec->adpc[%d]", i);
		      MPS_DEBUG_CDPE (s, sec->bdpc[i], "sec->bdpc[%d]", i);
                    }

#ifndef DISABLE_DEBUG
                  s->fp_iteration_time += mps_stop_timer (my_clock);
#endif
                  return -1;
                }

              /* Correct the radius */
              modcorr = cplx_mod (abcorr);
              s->frad[i] += modcorr;

              if (!s->again[i])
                computed_roots++;
            }
        }
    }

  /* Check if the roots are improvable in floating point */
  MPS_DEBUG_WITH_INFO (s, "Performed %d iterations with floating point arithmetic",
                       nit);

  if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)
      mps_dump (s, s->logstr);

  if (nit <= it_threshold || old_cr == computed_roots)
    {
      if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)
	{
	  MPS_DEBUG (s, "Setting approximation as best_approx");
	}
      s->secular_equation->best_approx = true;
    }
   mps_fcluster (s, 2.0 * s->n); 
   mps_fmodify (s); 

  for (i = 0; i < s->n; i++)
    {
      if (s->status[i][0] == 'C')
        s->status[i][0] = 'c';
    }

  /* Count time taken  */
#ifndef DISABLE_DEBUG
  s->fp_iteration_time += mps_stop_timer (my_clock);
#endif

  /* Return the number of approximated roots */
  return computed_roots;
}

/**
 * @brief Routine that performs a block of iteration
 * in DPE on the secular equation.
 *
 * @param s the pointer to the mps_status struct.
 * @param maxit Maximum number of iteration to perform.
 * @return The number of approximated roots after the iteration.
 */
int
mps_secular_ga_diterate (mps_status * s, int maxit)
{
  MPS_DEBUG_THIS_CALL;

  int computed_roots = 0;
  int iterations = 0;
  int i;
  int nit = 0;
  int it_threshold;
  int old_cr;

#ifndef DISABLE_DEBUG
  clock_t *my_clock = mps_start_timer ();
#endif

  /* Iterate with newton until we have good approximations
   * of the roots */

  for (i = 0; i < s->n; i++)
    {
      if (!s->again[i])
        computed_roots++;
    }

  /* Set the iterations threshold to 2 iterations
   * for every non approximated root. */
  it_threshold = 2 * (s->n - computed_roots);
  
  /* Save the old number of computed roots */
  old_cr = computed_roots;

  if (s->debug_level & MPS_DEBUG_PACKETS)
    MPS_DEBUG (s, "There are %d roots already approximated", computed_roots);

  /* Use this dump only for debugging purpose */
  /* mps_dump (s, s->logstr); */

  while (computed_roots < s->n && iterations < maxit - 1)
    {
      cdpe_t corr, abcorr;
      rdpe_t modcorr;

      /* Increase iterations counter */
      iterations++;

      for (i = 0; i < s->n; i++)
        {
          if (s->again[i])
            {
              nit++;
              mps_secular_dnewton (s, s->droot[i], s->drad[i], corr,
                                   &s->again[i], &i);
	      
              /* Apply Aberth correction */
              mps_daberth (s, i, abcorr);
              cdpe_mul_eq (abcorr, corr);
              cdpe_sub (abcorr, cdpe_one, abcorr);

              if (cdpe_ne (abcorr, cdpe_zero))
                {
                  cdpe_div (abcorr, corr, abcorr);
                  cdpe_sub_eq (s->droot[i], abcorr);

                  /* Correct the radius */
                  cdpe_mod (modcorr, abcorr);
                  rdpe_add_eq (s->drad[i], modcorr);
                }
              else
                {
                  s->again[i] = true;
                }

              if (!s->again[i])
                computed_roots++;
            }
        }
    }

  /* Check if no more than 2 iterations per root
   * were computed, and in that case state that
   * a coefficient regeneration won't be of much help */
  MPS_DEBUG_WITH_INFO (s, "Performed %d iterations", nit);
  if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)
    {
      mps_dump (s, s->logstr);
    }
  if (nit <= it_threshold || (old_cr == computed_roots))
    {
      if (s->debug_level & MPS_DEBUG_PACKETS)
	MPS_DEBUG (s, "Setting best_approx to true");
      s->secular_equation->best_approx = true;
    }

  mps_dcluster (s, 2.0 * s->n);
  mps_dmodify (s);

  for (i = 0; i < s->n; i++)
    {
      if (s->status[i][0] == 'C')
        s->status[i][0] = 'c';
      /* if (s->status[i][0] == 'a' || s->status[i][0] == 'i') */
      /*   s->again[i] = false; */
      /* else */
      /*   s->again[i] = true; */
    }

  /* These lines are used to debug the again vector, but are not useful
   * at the moment being */
  if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)
    {
      __MPS_DEBUG (s, "Again vector = ");
      for(i = 0; i < s->n; i++)
	{
	  fprintf (s->logstr, "%d ", s->again[i]);
	}
      fprintf (s->logstr, "\n");
    }

  /* Clock the routine */
#ifndef DISABLE_DEBUG
  s->dpe_iteration_time += mps_stop_timer (my_clock);
#endif

  /* Return the number of approximated roots */
  return computed_roots;
}

/**
 * @brief Routine that performs a block of iteration
 * in Multiprecision on the secular equation.
 *
 * @param s the pointer to the mps_status struct.
 * @param maxit Maximum number of iteration to perform.
 * @return The number of approximated roots after the iteration.
 */
int
mps_secular_ga_miterate (mps_status * s, int maxit)
{
  MPS_DEBUG_THIS_CALL;

  MPS_DEBUG_WITH_INFO (s, "Precision is at %ld bits", s->mpwp);
  MPS_DEBUG_RDPE (s, s->mp_epsilon, "Machine epsilon is s->mp_epsilon");

  int computed_roots = 0;
  int iterations = 0;
  int i, j, k;
  int nit = 0;
  int it_threshold;
  int old_cr;

  mpc_t corr, abcorr;
  cdpe_t ctmp;
  rdpe_t modcorr, rtmp;

#ifndef DISABLE_DEBUG
  clock_t *my_clock = mps_start_timer ();
#endif

  /* Init data with the right precision */
  mpc_init2 (corr, s->mpwp);
  mpc_init2 (abcorr, s->mpwp);

  /* Iterate with newton until we have good approximations
   * of the roots */
  
  for (i = 0; i < s->n; i++)
    {
      /* Set again to false if the root is already approximated */
      if (s->status[i][0] == 'a' || s->status[i][0] == 'i'
	  || s->status[i][0] == 'o')
	{
	  MPS_DEBUG_WITH_INFO (s, "Setting again[%d] to false", i);
	  s->again[i] = false;
	}

      if (s->again[i])
        s->rootwp[i] = s->mpwp;
      else
        computed_roots++;
    }

  if (s->debug_level & MPS_DEBUG_APPROXIMATIONS) 
    {
      MPS_DEBUG_WITH_INFO (s, "%d roots are already approximated on the start of miterate", computed_roots)
    }

  /* Set the iteration threshold to two times the remaining roots
   * to compute. */
  old_cr = computed_roots;
  it_threshold = 2 * (s->n - computed_roots);

  while (computed_roots < s->n && iterations < maxit)
    {
      /* Increase iterations counter */
      iterations++;

      for (i = 0; i < s->nclust; i++)
        {
          for (j = s->punt[i]; j < s->punt[i + 1]; j++)
            {
              k = s->clust[j];
              if (s->again[k])
                {
                  nit++;
                  mps_secular_mnewton (s, s->mroot[k], s->drad[k], corr,
                                       &s->again[k], &k);

                  /* Apply Aberth correction */
                  mps_maberth_s (s, k, i, abcorr);
                  mpc_mul_eq (abcorr, corr);
                  mpc_ui_sub (abcorr, 1, 0, abcorr);
                  mpc_div (abcorr, corr, abcorr);
                  mpc_sub_eq (s->mroot[k], abcorr);

                  /* Correct the radius */
                  mpc_get_cdpe (ctmp, abcorr);
                  cdpe_mod (modcorr, ctmp);
                  rdpe_add_eq (s->drad[k], modcorr);

                  mpc_get_cdpe (ctmp, s->mroot[k]);
                  cdpe_mod (rtmp, ctmp);
                  rdpe_div_eq (modcorr, rtmp);

                  if (!s->again[k])
                    computed_roots++;
                }
            }
        }
    }

  /* Deallocate multiprecision local variables */
  mpc_clear (abcorr);
  mpc_clear (corr);

  MPS_DEBUG_WITH_INFO (s, "Performed %d iterations", nit);
  if (nit <= it_threshold || computed_roots == old_cr)
    s->secular_equation->best_approx = true;

  /* Perform cluster analysis */
  mps_mcluster (s, 2.0 * s->n);
  mps_mmodify (s);

  for (i = 0; i < s->n; i++) 
    {
      if (s->status[i][0] == 'C') 
	s->status[i][0] = 'c';    
    }

  /* for (i = 0; i < s->n; i++) */
  /*   { */
  /*     if (s->status[i][0] == 'a' || s->status[i][0] == 'i') */
  /*       { */
  /*         s->again[i] = false; */
  /*       } */
  /*     else */
  /*       s->again[i] = true; */
  /*   } */


  /* These lines are used to debug the again vector, but are not useful
   * at the moment being */
  if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)
    {
      __MPS_DEBUG (s, "Again vector = ");
      for (i = 0; i < s->n; i++)
	{
	  fprintf (s->logstr, "%d ", s->again[i]);
	}
      fprintf (s->logstr, "\n");
      __MPS_DEBUG (s, "Status = ");
      for (i = 0; i < s->n; i++)
	{
	  fprintf (s->logstr, "%c ", s->status[i][0]);
	}
      // mps_dump (s, s->logstr);
    }
  
  /* Clock the routine */
#ifndef DISABLE_DEBUG
  s->mp_iteration_time += mps_stop_timer (my_clock);
#endif

  /* Return the number of approximated roots */
  return computed_roots;
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
 */
int
mps_secular_ga_regenerate_coefficients_mp (mps_status * s)
{
  /* Declaration and initialization of the multprecision
   * variables that are used only in that case */
  int i, j;
  mps_boolean success = true;
  int precision_increase_ratio = 2;
  int coeff_wp = precision_increase_ratio * s->mpwp;
  mps_secular_equation *sec = s->secular_equation;
  rdpe_t eps_tmp, rtmp, sec_eps, rtmp2;
  cdpe_t cdtmp, cdtmp2;

  mps_secular_raise_coefficient_precision (s, coeff_wp);

  if (MPS_REPRESENTATION_IS_MONOMIAL (s->input_config))
    {
      mpc_t diff, prod_b, m_one, ctmp;      
      MPS_DEBUG_WITH_INFO (s, "Regenerating coefficients from polynomial input");

      /* Init multiprecision values */
      mpc_init2 (diff, coeff_wp);
      mpc_init2 (prod_b, coeff_wp);
      mpc_init2 (m_one, coeff_wp);
      mpc_init2 (ctmp, coeff_wp);

      s->mpwp = coeff_wp;

      for (i = 0; i < s->n; ++i)
	{
	  mpc_set_prec (s->mfpc[i], coeff_wp);
	}

      mpc_set_ui (m_one, 0U, 0U);
      mpc_sub_eq (m_one, s->mfpc[s->n]);

      
      /*
       * The new coefficients of the secular equation can be computed
       * starting from the evaluation of the polynomial changed of
       * sign and divided for the product of the difference of the
       * b_i, i.e.:
       * 
       *   a_i = -p(b_i) / \prod_{i \neq j} (b_i - b_j)
       *
       */
      if (MPS_STRUCTURE_IS_INTEGER (s->input_config) 
	  || MPS_STRUCTURE_IS_RATIONAL (s->input_config))
	{
	  for (i = 0; i < s->n + 1; ++i)
	    {
	      mpf_set_q (mpc_Re (s->mfpc[i]), sec->initial_bmpqrc[i]);
	      mpf_set_q (mpc_Im (s->mfpc[i]), sec->initial_bmpqic[i]);

	      MPS_DEBUG_MPC (s, 15, s->mfpc[i], "s->mfpc[%d]", i);
	    }
	}

      for (i = 0; i < s->n; ++i)
	{
	  /* Evaluate the polynomial with absolute values to 
	   * compute the error */
	  mpc_get_cdpe (cdtmp, sec->bmpc[i]);
	  cdpe_mod (rtmp, cdtmp);
	  mps_aparhorner (s, s->n, rtmp, s->dap, s->spar, sec_eps, 0);

	  rdpe_set (sec_eps, s->dap[s->n]);
	  mpc_get_cdpe (cdtmp, sec->bmpc[i]);
	  cdpe_mod (rtmp2, cdtmp);
	  for (j = s->n - 1; j >= 0; j--)
	    {
	      rdpe_mul (rtmp, sec_eps, rtmp2);
	      rdpe_add (sec_eps, rtmp, s->dap[j]);
	    }
	  rdpe_mul_eq_d (sec_eps, s->n * 4);

	  /* Evaluate the polynomial with horner */
	  mps_parhorner (s, s->n, sec->bmpc[i], s->mfpc, s->spar, sec->ampc[i], 0);
	  mpc_set (sec->ampc[i], s->mfpc[s->n]); 
	  for (j = s->n - 1; j >= 0; j--) 
	    { 
	      mpc_div (ctmp, s->mfpc[j], sec->ampc[i]);
	      mpc_add_eq (ctmp, sec->bmpc[i]); 
	      mpc_mul_eq (sec->ampc[i], ctmp); 
	    } 
	  /* MPS_DEBUG_MPC (s, coeff_wp, sec->bmpc[i], "sec->bmpc[%d]", i); */
	  /* MPS_DEBUG_MPC (s, coeff_wp, sec->ampc[i], "sec->ampc[%d]", i); */

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

	  MPS_DEBUG_MPC (s, 15, prod_b, "prod_b");
	  
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

	  MPS_DEBUG_RDPE (s, sec->dregeneration_epsilon[i], "error on b_%d differences", i);
	  MPS_DEBUG_RDPE (s, sec_eps, "sec_eps");

	  /* Finalize error computation */
	  rdpe_add_eq (sec->dregeneration_epsilon[i], sec_eps);
	  rdpe_mul_eq (sec->dregeneration_epsilon[i], s->mp_epsilon);

	  if (s->debug_level & MPS_DEBUG_REGENERATION)
	    {
	      MPS_DEBUG_RDPE (s, sec->dregeneration_epsilon[i],
			      "Relative error on a_%d", i);
	    }
	}

      s->mpwp = coeff_wp / precision_increase_ratio;

      /* Clear requested storage */
      mpc_clear (diff);
      mpc_clear (prod_b);
      mpc_clear (m_one);
      mpc_clear (ctmp);
    }
  else if (MPS_REPRESENTATION_IS_SECULAR (s->input_config))
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
      if (MPS_STRUCTURE_IS_RATIONAL (s->input_config) ||
	  MPS_STRUCTURE_IS_INTEGER (s->input_config))
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

	      /* Compute the local relative error that ios introduce here, as 
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

	  if (s->debug_level & MPS_DEBUG_REGENERATION)
	    {
	      MPS_DEBUG_RDPE (s, sec_eps, "Relative error on sec_ev(b_%d)", i);
	      MPS_DEBUG_RDPE (s, sec->dregeneration_epsilon[i],
			      "Relative error on a_%d", i);
	    }
	}

    regenerate_m_exit:

      /* Free data */
      mpc_clear (prod_b);
      mpc_clear (sec_ev);
      mpc_clear (ctmp);
      mpc_clear (btmp);


    } /* End of if (MPS_REPRESENTATION_IS_SECULAR (s->input_config)) */
  
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
void
mps_secular_ga_regenerate_coefficients (mps_status * s)
{
  MPS_DEBUG_THIS_CALL;

  cplx_t *old_b, *old_a;
  cdpe_t *old_db, *old_da;
  mpc_t *old_ma, *old_mb;
  mps_secular_equation *sec;
  int i, j;

  /* Start timer and add execution time to the total counter */
#ifndef DISABLE_DEBUG
  clock_t *my_clock = mps_start_timer ();
#endif

  sec = (mps_secular_equation *) s->secular_equation;

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
          cplx_set (sec->bfpc[i], s->froot[i]);
          mpc_set_cplx (sec->bmpc[i], sec->bfpc[i]);
          mpc_set_cplx (sec->ampc[i], sec->afpc[i]);
        }

      /* Regeneration */
      if (!mps_secular_ga_regenerate_coefficients_mp (s))
        {
          MPS_DEBUG_WITH_INFO (s,
                               "Regeneration of coefficients failed, reusing old ones");
          for (i = 0; i < s->n; i++)
            {
              cplx_set (sec->afpc[i], old_a[i]);
              cplx_set (sec->bfpc[i], old_b[i]);
            }
        }
      else
        {
          for (i = 0; i < s->n; i++)
            {
              mpc_get_cplx (sec->bfpc[i], sec->bmpc[i]);
              mpc_get_cplx (sec->afpc[i], sec->ampc[i]);
	      sec->fregeneration_epsilon[i] = rdpe_get_d (sec->dregeneration_epsilon[i]) + DBL_EPSILON;

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

      /* Regeneration */
      if (!mps_secular_ga_regenerate_coefficients_mp (s))
        {
          MPS_DEBUG_WITH_INFO (s,
                               "Regeneration of the coefficients failed, reusing old ones.")
            for (i = 0; i < s->n; i++)
            {
              cdpe_set (sec->adpc[i], old_da[i]);
              cdpe_set (sec->bdpc[i], old_db[i]);
            }
        }
      else
        {
          for (i = 0; i < s->n; i++)
	    {
	      mpc_get_cdpe (sec->adpc[i], sec->ampc[i]);
	      rdpe_add_eq_d (sec->dregeneration_epsilon[i], DBL_EPSILON);
	    }
          mps_secular_set_radii (s);
        }

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

      /* Regeneration */
      if (mps_secular_ga_regenerate_coefficients_mp (s))
        {
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
}


/**
 * @brief Check if iterations can terminate, i.e. if newton 
 * isolation has been reached. If the target was approximation then 
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

    /* if  (!MPS_STRUCTURE_IS_FP (s->secular_equation->input_structure) */
    /*   && s->lastphase != mp_phase) */
    /* return false; */

  for (i = 0; i < s->n; i++)
    {
      switch (s->lastphase)
        {
          /* Float case */
        case float_phase:
          if (s->status[i][0] != 'i' && s->status[i][0] != 'a'
              && s->status[i][0] != 'o')
            {
              MPS_DEBUG_WITH_INFO (s,
                                   "Root %d is not isolated, nor approximated, so we can't stop now.",
                                   i);
              return false;
            }
          break;

          /* Multiprecision and DPE case are the same, since the radii
           * are always RDPE. */
        case mp_phase:
          if (s->status[i][0] != 'i' && s->status[i][0] != 'a'
              && s->status[i][0] != 'o')
            {
              MPS_DEBUG_WITH_INFO (s,
                                   "Root %d is not isolated, nor approximated, so we can't stop now.",
                                   i);
	      MPS_DEBUG_WITH_INFO (s,
				   "Status of root %d: %c", i, s->status[i][0]);
              return false;
            }
          break;
        case dpe_phase:
          if (s->status[i][0] != 'i' && s->status[i][0] != 'a'
              && s->status[i][0] != 'o')
            {
              MPS_DEBUG_WITH_INFO (s,
                                   "Root %d is not isolated, nor approximated, so we can't stop now.",
                                   i);
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
  for (i = 0; i < s->n; i++)
    {
      if (MPS_REPRESENTATION_IS_SECULAR (s->input_config))
	{
	  if (MPS_STRUCTURE_IS_FP (s->input_config))
	    {
	      mpc_set (s->secular_equation->ampc[i],
		       s->secular_equation->initial_ampc[i]);
	      mpc_set (s->secular_equation->bmpc[i],
		       s->secular_equation->initial_bmpc[i]);
	    }
	  else
	    {
	      mpc_set_q (s->secular_equation->ampc[i],
			 s->secular_equation->initial_ampqrc[i],
			 s->secular_equation->initial_ampqic[i]);

	      mpc_set_q (s->secular_equation->bmpc[i],
			 s->secular_equation->initial_bmpqrc[i],
			 s->secular_equation->initial_bmpqic[i]);
	    }
	}
      else if (MPS_REPRESENTATION_IS_MONOMIAL (s->input_config))
	{
	  
	}
      rdpe_set (sec->dregeneration_epsilon[i],
		s->mp_epsilon);
    }

  mpc_t nwtcorr;
  cdpe_t ctmp;
  rdpe_t rtmp, old_rad;

  mpc_init2 (nwtcorr, s->mpwp);

  int starting_precision = s->mpwp;

  for (i = 0; i < s->n; i++)
    {
      int j;

      /* Reset precision of coefficients and of the root we're 
       * interested in. */
      if (s->mpwp != starting_precision)
        {
          mps_secular_raise_coefficient_precision (s, starting_precision);
          mpc_set_prec (s->mroot[i], starting_precision);
          mpc_set_prec (nwtcorr, starting_precision);
          s->mpwp = starting_precision;
        }

      mpc_get_cdpe (ctmp, s->mroot[i]);
      cdpe_mod (rtmp, ctmp);
      rdpe_div_eq (rtmp, s->drad[i]);

      /* Find correct digits and maximum number of iterations */
      int correct_digits = rdpe_log10 (rtmp);
      if (s->debug_level & MPS_DEBUG_IMPROVEMENT)
	MPS_DEBUG (s, "Root %d has %d correct digits", i, correct_digits);
      int iterations =
        log (s->output_config->prec / correct_digits / LOG2_10) * LOG2_10 + 1;
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
	  rdpe_set (old_rad, s->drad[i]);
          mps_secular_mnewton (s, s->mroot[i], s->drad[i], nwtcorr,
                               &s->again[i], &i);
          mpc_sub_eq (s->mroot[i], nwtcorr);

	  /* Compute quadratic radius */
	  mpc_get_cdpe (ctmp, s->mroot[i]);
	  cdpe_mod (rtmp, ctmp);
	  rdpe_div_eq (old_rad, rtmp);
	  rdpe_mul_eq (old_rad, old_rad);
	  rdpe_mul_eq (old_rad, rtmp);

	  if (rdpe_lt (old_rad, s->drad[i]))  
	    rdpe_set (s->drad[i], old_rad);  

          /* Debug iterations */
          if (s->debug_level & MPS_DEBUG_IMPROVEMENT)
            {
              MPS_DEBUG_MPC (s, 10, s->mroot[i], "s->mroot[%d]", i);
              MPS_DEBUG_RDPE (s, s->drad[i], "s->drad[%d]", i);
            }

          /* Check if the approximation is already good. */
          mpc_get_cdpe (ctmp, s->mroot[i]);
          cdpe_mod (rtmp, ctmp);
          rdpe_div (rtmp, s->drad[i], rtmp);

          if (rdpe_le (rtmp, s->eps_out))
            {
              s->status[i][0] = 'a';
              break;
            }
          else
            {
              s->mpwp *= 2;
              mps_secular_raise_coefficient_precision (s, s->mpwp);
              mpc_set_prec (nwtcorr, s->mpwp);
              mpc_set_prec (s->mroot[i], s->mpwp);
            }
        }

      /* Since we have passed the bound of the maximum allowed iterations
       * and quadratic convergence was guaranteed, the root is now 
       * approximated. */
      s->status[i][0] = 'a';
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
  int iteration_per_packet = 10;
  int i;
  mps_boolean skip_check_stop = false;
  rdpe_t r_eps;
  mps_secular_equation *sec = mps_secular_equation_from_status (s);
  mps_phase phase = sec->starting_case;

  rdpe_set_d (r_eps, DBL_EPSILON);

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

  if (MPS_REPRESENTATION_IS_SECULAR (s->input_config))
    s->data_type = "uri";
  else if (MPS_REPRESENTATION_IS_MONOMIAL (s->input_config))
    s->data_type = "dri";

  /* We set the selected phase */
  s->lastphase = phase;

  /* Manually set FILE* pointer for streams.
   * More refined options will be added later. */
  s->outstr = s->rtstr = stdout;
  packet = 0;

  /* Set the maximum possible radius even in DPE, since
   * we may be starting directly from the DPE phase */
  for (i = 0; i < s->n; i++)
    {
      s->frad[i] = DBL_MAX;
      rdpe_set (s->drad[i], RDPE_MAX);
    }
  
  /* Set initial cluster structure as no cluster structure. */
  mps_cluster_reset (s);

  /* Set phase */
  s->lastphase = phase;

  /* If the input was polynomial we need to determined the secular
   * coefficients */
  if (MPS_REPRESENTATION_IS_MONOMIAL (s->input_config))
    {
      int nit;
      mps_boolean excep;

      if (s->lastphase == float_phase)
	mps_fstart (s, s->n, 0, 0.0, 0.0, s->eps_out, s->fap);
      else
	mps_dstart (s, s->n, 0, (__rdpe_struct *) rdpe_zero,
		    (__rdpe_struct *) rdpe_zero, s->eps_out,
		    s->dap);

      mps_secular_ga_regenerate_coefficients (s);
    }

  /* Select initial approximations using the custom secular
   * routine and based on the phase selected by the user. */
  switch (s->lastphase)
    {
    case  float_phase: 
       mps_secular_fstart (s, s->n, 0, 0.0, 0.0, s->eps_out); 
       break;

     case dpe_phase: 
       mps_secular_dstart (s, s->n, 0, (__rdpe_struct *) rdpe_zero, 
                           (__rdpe_struct *) rdpe_zero, s->eps_out); 
       break; 

     case mp_phase: 
       mps_secular_mstart (s, s->n, 0, (__rdpe_struct *) rdpe_zero, 
                           (__rdpe_struct *) rdpe_zero, s->eps_out); 
       break;

     default: 
       break;
    }

  /* Set initial radius */
  mps_secular_set_radii (s);

  for (i = 0; i < s->n; i++)
    {
      s->again[i] = true;
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
          roots_computed = mps_secular_ga_fiterate (s, iteration_per_packet);
          /* If the computation fails we need to switch to DPE so do not
           * break here, but continue the cycle. */
          if (roots_computed != -1)
            {
              MPS_DEBUG_WITH_INFO (s,
                                   "%d roots were computed to the best precision available",
                                   roots_computed);
              break;
            }

        case dpe_phase:
          MPS_DEBUG_WITH_INFO (s, "Starting DPE iterations");
          roots_computed = mps_secular_ga_diterate (s, iteration_per_packet);
          MPS_DEBUG_WITH_INFO (s,
                               "%d roots were computed to the best precision available",
                               roots_computed);
          break;

        case mp_phase:
          MPS_DEBUG_WITH_INFO (s, "Starting MP iterations");
          roots_computed = mps_secular_ga_miterate (s, iteration_per_packet);
          MPS_DEBUG_WITH_INFO (s,
                               "%d roots were computed to the best precision available",
                               roots_computed);
          break;

        default:
          break;
        }

      packet++;

      /* Check if all roots were approximated with the
       * given input precision                      */      
      if (mps_secular_ga_check_stop (s))
        {
          ; // break;
        }

      /* If we can't stop recompute coefficients in higher precision and
       * continue to iterate, unless the best approximation possible in
       * this precision has been reached. In that case increase the precision
       * of the computation. */
      if (s->lastphase != mp_phase)
	{
	  if (sec->best_approx)
	    mps_secular_ga_regenerate_coefficients (s);
	  else
	    skip_check_stop = true;
	}
      
      if (packet > s->max_pack)
	{
	  mps_error (s, 1, "Maximum number of iteration passed. Aborting.");
	  return;
	}

      if (s->secular_equation->best_approx && packet > 4)
        {
          /* Going to multiprecision if we're not there yet */
          if (s->lastphase != mp_phase)
            {
              mps_secular_switch_phase (s, mp_phase);
	      mps_secular_ga_regenerate_coefficients (s);
	    }
          else
            {
              /* Raising precision otherwise */
              mps_secular_raise_precision (s, 2 * s->mpwp);
              mps_secular_ga_regenerate_coefficients (s);
              skip_check_stop = true;

	      for (i = 0; i < s->n; ++i)
		{
		  if (s->status[i][0] != 'a' ||
		      s->status[i][0] != 'o' ||
		      s->status[i][0] != 'i')
		    {
		      s->again[i] = true;
		    }
		}
            }

          packet = 0;
        }

      if (s->lastphase != mp_phase)
	skip_check_stop = true;
    }
  while (skip_check_stop || !mps_secular_ga_check_stop (s));

  /* Finally improve the roots if approximation is required */
  if (s->goal[0] == 'a')
    {
      if (MPS_REPRESENTATION_IS_SECULAR (s->input_config))
	{
	  clock_t *my_timer = mps_start_timer ();
	  mps_secular_ga_improve (s);
	  unsigned int improve_time = mps_stop_timer (my_timer);
	  if (s->debug_level & MPS_DEBUG_TIMINGS)
	    {
	      MPS_DEBUG (s, "mps_secular_ga_improve took %u ms", improve_time);
	    }
	}
      else if (MPS_REPRESENTATION_IS_MONOMIAL (s->input_config))
	{
	  clock_t *my_timer = mps_start_timer ();
	  mps_secular_ga_improve (s);
	  unsigned int improve_time = mps_stop_timer (my_timer);
	  if (s->debug_level & MPS_DEBUG_TIMINGS)
	    {
	      MPS_DEBUG (s, "mps_improve took %u ms", improve_time);
	    }
	}
    }

  /* Debug total time taken but only if debug is enabled */
#ifndef DISABLE_DEBUG
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
                 mps_stop_timer (total_clock));
    }
#endif

}
