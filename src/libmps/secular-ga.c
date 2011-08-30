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
  MPS_DEBUG_THIS_CALL int computed_roots = 0;
  int iterations = 0;
  int i;
  int nit = 0;

  mps_secular_equation *sec = s->secular_equation;

  double old_rad;
  cplx_t old_root;

  /* Iterate with newton until we have good approximations
   * of the roots */
  mps_update (s);

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
				   &s->again[i]);

	      /* Apply Aberth correction */
	      mps_faberth (s, i, abcorr);
	      cplx_mul_eq (abcorr, corr);
	      cplx_sub (abcorr, cplx_one, abcorr);
	      cplx_div (abcorr, corr, abcorr);
	      cplx_sub_eq (s->froot[i], abcorr);

	      /* Check if we need to switch to DPE */
	      if (isnan (cplx_Re (s->froot[i]))
		  || isinf (cplx_Re (s->froot[i]))
		  || isnan (cplx_Im (s->froot[i]))
		  || isinf (cplx_Im (s->froot[i])) || isnan (s->frad[i])
		  || isinf (s->frad[i]))
		{
		  MPS_DEBUG (s,
			     "Switching to DPE phase because NAN or INF was introduced in computation");
		  cplx_set (s->froot[i], old_root);
		  s->frad[i] = old_rad;
		  s->lastphase = dpe_phase;

		  /* Copy roots and radius */
		  for (i = 0; i < s->n; i++)
		    {
		      cdpe_set_x (s->droot[i], s->froot[i]);
		      rdpe_set_d (s->drad[i], s->frad[i]);
		    }
		  return computed_roots;
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
  MPS_DEBUG (s, "Performed %d iterations with floating point arithmetic",
	     nit);
  if (nit <= 2 * s->n)
    s->secular_equation->best_approx = true;

  mps_fcluster (s, 2.0 * s->n);
  mps_fmodify (s);

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
  MPS_DEBUG_THIS_CALL int computed_roots = 0;
  int iterations = 0;
  int i;
  int nit = 0;

  /* Iterate with newton until we have good approximations
   * of the roots */
  mps_update (s);

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
				   &s->again[i]);

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
  if (nit <= 2 * s->n)
    {
      s->secular_equation->best_approx = true;
    }

  mps_dcluster (s, 2.0 * s->n);
  mps_dmodify (s);

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
  MPS_DEBUG_THIS_CALL int computed_roots = 0;
  int iterations = 0;
  int i;
  int nit = 0;

  mpc_t corr, abcorr;
  cdpe_t ctmp;
  rdpe_t modcorr;

  /* Init data with the right precision */
  mpc_init2 (corr, s->mpwp);
  mpc_init2 (abcorr, s->mpwp);

  /* Iterate with newton until we have good approximations
   * of the roots */
  mps_update (s);

  for (i = 0; i < s->n; i++)
    {
      if (s->again[i])
	s->rootwp[i] = s->mpwp;
    }

  while (computed_roots < s->n && iterations < maxit)
    {
      /* Increase iterations counter */
      iterations++;

      for (i = 0; i < s->n; i++)
	{
	  if (s->again[i])
	    {
	      nit++;
	      mps_secular_mnewton (s, s->mroot[i], s->drad[i], corr,
				   &s->again[i]);

	      /* Apply Aberth correction */
	      mps_maberth (s, i, abcorr);
	      mpc_mul_eq (abcorr, corr);
	      mpc_ui_sub (abcorr, 1, 0, abcorr);
	      mpc_div (abcorr, corr, abcorr);
	      mpc_sub_eq (s->mroot[i], abcorr);

	      /* Correct the radius */
	      mpc_get_cdpe (ctmp, abcorr);
	      cdpe_mod (modcorr, ctmp);
	      rdpe_add_eq (s->drad[i], modcorr);

	      if (!s->again[i])
		computed_roots++;
	    }
	}
    }

  /* Deallocate multiprecision local variables */
  mpc_clear (abcorr);
  mpc_clear (corr);

  MPS_DEBUG (s, "Performed %d iterations", nit);
  if (nit <= 2 * s->n)
    s->secular_equation->best_approx = true;

  /* Perform cluster analysis */
  mps_mcluster (s, 2.0 * s->n);
  mps_mmodify (s);

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
  int success = 1;
  int coeff_wp = 2 * s->mpwp;
  mpc_t prod_b, sec_ev;
  mpc_t ctmp, btmp;
  mps_secular_equation *sec = s->secular_equation;

regenerate_m_start:

  mpc_init2 (prod_b, coeff_wp);
  mpc_init2 (sec_ev, coeff_wp);
  mpc_init2 (ctmp, coeff_wp);
  mpc_init2 (btmp, coeff_wp);

  mps_secular_raise_coefficient_precision (s, coeff_wp);

  /* If the input was rational generate multiprecision floating
   * point coefficient from the ones that the user originally
   * provided */
  if (MPS_STRUCTURE_IS_RATIONAL (sec->input_structure) ||
      MPS_STRUCTURE_IS_INTEGER (sec->input_structure))
    {
      MPS_DEBUG (s,
		 "Regenerating coefficients from the multiprecision input");
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

      for (j = 0; j < sec->n; j++)
	{
	  /* Compute 1 / (b_i - old_b_j) */
	  mpc_sub (btmp, sec->bmpc[i], sec->initial_bmpc[j]);

	  /* If b - old_b is zero, abort the computation */
	  if (mpc_eq_zero (btmp))
	    {
	      success = -1;

	      MPS_DEBUG (s,
			 "Cannot regenerate coefficients, reusing old ones and setting best_approx to true.");
	      s->secular_equation->best_approx = true;

	      goto regenerate_m_exit;
	    }

	  mpc_inv (ctmp, btmp);

	  /* Add a_j / (b_i - old_b_j) to sec_ev */
	  mpc_mul_eq (ctmp, sec->initial_ampc[j]);
	  mpc_add_eq (sec_ev, ctmp);

	  /* Multiply prod_b for
	   * b_i - b_j if i \neq j and prod_old_b
	   * for b_i - old_b_i.  */
	  mpc_mul_eq (prod_b, btmp);
	  if (i != j)
	    {
	      mpc_sub (ctmp, sec->bmpc[i], sec->bmpc[j]);
	      mpc_div_eq (prod_b, ctmp);
	    }
	}

      /* Compute the new a_i as sec_ev * prod_old_b / prod_b */
      mpc_sub_eq_ui (sec_ev, 1, 0);
      mpc_mul (sec->ampc[i], sec_ev, prod_b);
    }


  mps_secular_raise_coefficient_precision (s, s->mpwp);


regenerate_m_exit:

  /* Free data */
  mpc_clear (prod_b);
  mpc_clear (sec_ev);
  mpc_clear (ctmp);
  mpc_clear (btmp);

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
  MPS_DEBUG_THIS_CALL cplx_t *old_b, *old_a;
  cdpe_t *old_db, *old_da;
  mpc_t *old_ma, *old_mb;
  mps_secular_equation *sec;
  int i, j;

  sec = (mps_secular_equation *) s->secular_equation;

  switch (s->lastphase)
    {

      /* If we are in the float phase regenerate coefficients
       * starting from floating point */
    case float_phase:

      s->mpwp = 56;

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
      if (mps_secular_ga_regenerate_coefficients_mp (s) != 0)
	{
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
	    }
	}

      cplx_vfree (old_a);
      cplx_vfree (old_b);

      /* mps_secular_fstart(s, s->n, 0, 0, 0, s->eps_out); */

      break;

      /* If this is the DPE phase regenerate DPE coefficients */
    case dpe_phase:

      s->mpwp = 56;

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
      if (mps_secular_ga_regenerate_coefficients_mp (s) != 0)
	{
	  for (i = 0; i < s->n; i++)
	    {
	      cdpe_set (sec->adpc[i], old_da[i]);
	      cdpe_set (sec->bdpc[i], old_db[i]);
	    }
	}
      else
	{
	  for (i = 0; i < s->n; i++)
	    mpc_get_cdpe (sec->adpc[i], sec->ampc[i]);
	}

      /* Free data */
      cdpe_vfree (old_da);
      cdpe_vfree (old_db);

      /* mps_secular_dstart(s, s->n, 0, (__rdpe_struct *) rdpe_zero,
         (__rdpe_struct *) rdpe_zero, s->eps_out); */

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


      /* This is too much verbose to be enabled even in the debugged version,
       * so do not display it (unless we are trying to catch some errors on
       * coefficient regeneration). */

      /*
         for (i = 0; i < s->n; i++)
         {
         MPS_DEBUG_MPC(s, 15, sec->ampc[i], "sec->ampc[%d]", i);
         MPS_DEBUG_MPC(s, 15, sec->bmpc[i], "sec->bmpc[%d]", i);
         }
       */

      mpc_vclear (old_ma, s->n);
      mpc_vclear (old_mb, s->n);
      mpc_vfree (old_ma);
      mpc_vfree (old_mb);

      mps_secular_mstart (s, s->n, 0, (__rdpe_struct *) rdpe_zero,
			  (__rdpe_struct *) rdpe_zero, s->eps_out);

      break;

    default:
      break;

    }				/* End of switch (s->lastphase) */
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
  MPS_DEBUG_THIS_CALL double frad = pow (10, -s->prec_out);
  rdpe_t drad, root_mod;
  cdpe_t root;
  int i;

  /* Set the maximum rad allowed */
  rdpe_set_2dl (drad, 1.0, -s->prec_out);

  /* If we are in floating point and output precision is higher of
   * what is reachable, standing in floating point, do not exit */
  if (frad == 0 && s->lastphase == float_phase)
    return false;

  for (i = 0; i < s->n; i++)
    {
      switch (s->lastphase)
	{

	  /* Float case */
	case float_phase:
	  if (s->status[i][0] != 'i' && s->status[i][0] != 'a')
	    return false;

	  /* Multiprecision and DPE case are the same, since the radii
	   * are always RDPE. */
	case mp_phase:
	  if (s->status[i][0] != 'i' && s->status[i][0] != 'a')
	    return false;
	  break;
	case dpe_phase:
	  if (s->status[i][0] != 'i' && s->status[i][0] != 'a')
	    return false;
	  break;

	default:
	  break;

	}
    }

  MPS_DEBUG (s, "Stop conditions were satisfied");
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

  /* Before improving with newton we shall reuse
   * the original coefficients */
  for (i = 0; i < s->n; i++)
    {
      if (MPS_STRUCTURE_IS_FP (s->secular_equation->input_structure))
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

  /* Improve the roots with newton */
  /* mps_improve (s);
     return; */
  
  // We should check the condition number here, or use the old mps_improve ()
  // functions that already considers it, adapting it to the new structure.
#define MAX_ITERATIONS 150

  mpc_t nwtcorr;
  cdpe_t ctmp;
  rdpe_t rtmp;
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

      for (j = 0; j < MAX_ITERATIONS; j++)
	{
	  mps_secular_mnewton (s, s->mroot[i], s->drad[i], nwtcorr,
			       &s->again[i]);
	  mpc_sub_eq (s->mroot[i], nwtcorr);

	  /* Debug iterations */
	  MPS_DEBUG_MPC (s, 10, s->mroot[i], "s->mroot[%d]", i);
	  MPS_DEBUG_RDPE (s, s->drad[i], "s->drad[%d]", i);

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
    }
  mpc_clear (nwtcorr);

#undef MAX_ITERATIONS
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
  mps_secular_equation *sec = mps_secular_equation_from_status (s);
  mps_phase phase = sec->starting_case;

  /* Set the output desired for the output */
  rdpe_set_2dl (s->eps_out, 1.0, -s->prec_out);

  /* Set degree and allocate polynomial-related variables
   * to allow initializitation to be performed. */
  s->deg = s->n = sec->n;
  s->data_type = "uri";

  /* We set the selected phase */
  s->lastphase = phase;

  /* Allocate other data */
  mps_allocate_data (s);

  for (i = 0; i < s->n; i++)
    {
      s->status[i][1] = 'w';
      s->status[i][2] = 'u';
    }

  /* Manually set FILE* pointer for streams.
   * More refined options will be added later. */
  s->outstr = s->rtstr = stdout;
  packet = 0;

  for (i = 0; i < s->n; i++)
    s->frad[i] = DBL_MAX;

  /* Set initial cluster structure as no cluster structure. */
  mps_cluster_reset (s);

  /* Set phase */
  s->lastphase = phase;

  /* Select initial approximations using the custom secular
   * routine and based on the phase selected by the user. */
  switch (s->lastphase)
    {
    case float_phase:
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

  /* Cycle until approximated */
  do
    {
      s->secular_equation->best_approx = false;
      /* Perform an iteration of floating point Aberth method */
      switch (s->lastphase)
	{
	case float_phase:
	  roots_computed = mps_secular_ga_fiterate (s, iteration_per_packet);
	  MPS_DEBUG (s, "%d roots were computed", roots_computed);
	  break;

	case dpe_phase:
	  roots_computed = mps_secular_ga_diterate (s, iteration_per_packet);
	  MPS_DEBUG (s, "%d roots were computed", roots_computed);
	  break;

	case mp_phase:
	  roots_computed = mps_secular_ga_miterate (s, iteration_per_packet);
	  MPS_DEBUG (s, "%d roots were computed", roots_computed);
	  break;

	default:
	  break;
	}

      packet++;

      /* Check if all roots were approximated with the
       * given input precision                      */
      /* if (mps_secular_ga_check_stop (s))
         {
         mps_improve (s);
         return;
         } */

      /* If we can't stop recompute coefficients in higher precision and
       * continue to iterate, unless the best approximation possible in
       * this precision has been reached. In that case increase the precision
       * of the computation. */
      if (!sec->best_approx)
	{
	  mps_secular_ga_regenerate_coefficients (s);
	}

      /* Instead of using else we recheck best approx because it could
       * have been set by the coefficient regeneration */
      if (sec->best_approx || packet > 5)
	{
	  /* Going to multiprecision if we're not there yet */
	  if (s->lastphase != mp_phase)
	    mps_secular_switch_phase (s, mp_phase);
	  else
	    {
	      /* Raising precision otherwise */
	      mps_secular_raise_precision (s, 2 * s->mpwp);
	      mps_secular_ga_regenerate_coefficients (s);
	    }

	  packet = 0;
	}
    }
  while (!mps_secular_ga_check_stop (s));

  /* Finally improve the roots if approximation is required */
  mps_secular_ga_improve (s);
}
