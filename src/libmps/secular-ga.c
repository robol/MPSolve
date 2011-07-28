/*
 * secular-ga.c
 *
 *  Created on: 15/giu/2011
 *      Author: leonardo
 */

#include <mps/core.h>
#include <mps/poly.h>
#include <mps/link.h>
#include <mps/secular.h>

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
mps_secular_ga_fiterate(mps_status* s, int maxit)
{
  MPS_DEBUG_THIS_CALL

  int computed_roots = 0;
  int iterations = 0;
  int i;
  int nit = 0;

  mps_secular_equation* sec = s->secular_equation;

  double old_rad;
  cplx_t old_root;

  /* Iterate with newton until we have good approximations
   * of the roots */
  /* Set again to true */
  for (i = 0; i < s->n; i++)
    s->again[i] = true;

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
              cplx_set(old_root, s->froot[i]);
              old_rad = s->frad[i];
              mps_secular_fnewton(s, s->froot[i], &s->frad[i], corr,
                  &s->again[i]);

              /* Apply Aberth correction */
              mps_faberth(s, i, abcorr);
              cplx_mul_eq(abcorr, corr);
              cplx_sub(abcorr, cplx_one, abcorr);
              cplx_div(abcorr, corr, abcorr);
              cplx_sub_eq(s->froot[i], abcorr);

              /* Check if we need to switch to DPE */
              if (isnan(cplx_Re(s->froot[i])) || isinf(cplx_Re(s->froot[i])) ||
                  isnan(cplx_Im(s->froot[i])) || isinf(cplx_Im(s->froot[i])) ||
                  isnan(s->frad[i]) || isinf(s->frad[i]))
                {
                  MPS_DEBUG(s, "Switching to DPE phase because NAN or INF was introduced in computation")
                  cplx_set(s->froot[i], old_root);
                  s->frad[i] = old_rad;
                  s->lastphase = dpe_phase;

                  /* Copy roots and radius */
                  for(i = 0; i < s->n; i++)
                  {
                      cdpe_set_x(s->droot[i], s->froot[i]);
                      rdpe_set_d(s->drad[i], s->frad[i]);
                  }
                  return computed_roots;
                }

              /* Correct the radius */
              modcorr = cplx_mod(abcorr);
              s->frad[i] += modcorr;

              if (!s->again[i])
                computed_roots++;
            }
        }
    }

  /* Check if the roots are improvable in floating point */
  MPS_DEBUG(s, "Performed %d iterations", nit)
  if (nit <= 2 * s->n)
      s->secular_equation->best_approx = true;

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
mps_secular_ga_diterate(mps_status* s, int maxit)
{
  MPS_DEBUG_THIS_CALL

  int computed_roots = 0;
  int iterations = 0;
  int i;
  int nit = 0;

  /* Iterate with newton until we have good approximations
   * of the roots */
  /* Set again to true */
  for (i = 0; i < s->n; i++)
    s->again[i] = true;

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
              mps_secular_dnewton(s, s->droot[i], s->drad[i], corr,
                  &s->again[i]);

              /* Apply Aberth correction */
              mps_daberth(s, i, abcorr);
              cdpe_mul_eq(abcorr, corr);
              cdpe_sub(abcorr, cdpe_one, abcorr);

              if (cdpe_ne(abcorr, cdpe_zero))
              {
                cdpe_div(abcorr, corr, abcorr);
                cdpe_sub_eq(s->droot[i], abcorr);

                /* Correct the radius */
                cdpe_mod(modcorr, abcorr);
                rdpe_add_eq(s->drad[i], modcorr);
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
mps_secular_ga_miterate(mps_status* s, int maxit)
{
  MPS_DEBUG_THIS_CALL

  int computed_roots = 0;
  int iterations = 0;
  int i;
  int nit = 0;

  mpc_t corr, abcorr;
  cdpe_t ctmp;
  rdpe_t modcorr, drad;

  rdpe_set_2dl(drad, 1.0, -s->prec_out);

  /* Init data with the right precision */
  mpc_init2(corr, s->mpwp);
  mpc_init2(abcorr, s->mpwp);

  /* Iterate with newton until we have good approximations
   * of the roots */
  /* Allocate again and set all to true */
  for (i = 0; i < s->n; i++)
    {
      if (rdpe_gt(s->drad[i], drad))
        s->again[i] = true;
      else
        {
          s->again[i] = false;
          computed_roots++;
        }
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
              mps_secular_mnewton(s, s->mroot[i], s->drad[i], corr,
                  &s->again[i]);

              /* Apply Aberth correction */
              mps_maberth(s, i, abcorr);
              mpc_mul_eq(abcorr, corr);
              mpc_ui_sub(abcorr, 1, 0, abcorr);
              mpc_div(abcorr, corr, abcorr);
              mpc_sub_eq(s->mroot[i], abcorr);

              /* Correct the radius */
              mpc_get_cdpe(ctmp, abcorr);
              cdpe_mod(modcorr, ctmp);
              rdpe_add_eq(s->drad[i], modcorr);

              if (!s->again[i])
                computed_roots++;
            }
        }
    }

  /* Deallocate multiprecision local variables */
  mpc_clear(abcorr);
  mpc_clear(corr);

  MPS_DEBUG(s, "Performed %d iterations", nit)
  if (nit <= 2 * s->n)
          s->secular_equation->best_approx = true;

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
 */
int
mps_secular_ga_regenerate_coefficients_mp(mps_status* s)
{
    /* Declaration and initialization of the multprecision
     * variables that are used only in that case */
    int i, j;
    int success = 0;
    int coeff_wp = 2 * s->mpwp;
    mpc_t prod_b, sec_ev;
    mpc_t ctmp, btmp;
    mps_secular_equation* sec = s->secular_equation;

    regenerate_m_start:

    mpc_init2(prod_b, coeff_wp);
    mpc_init2(sec_ev, coeff_wp);
    mpc_init2(ctmp, coeff_wp);
    mpc_init2(btmp, coeff_wp);

    mps_secular_raise_coefficient_precision(s, coeff_wp);

    /* Compute the new a_i */
    for (i = 0; i < s->n; i++)
      {
        mpc_set_ui(prod_b, 1, 0);
        mpc_set_ui(sec_ev, 0, 0);

        for (j = 0; j < sec->n; j++)
          {
            /* Compute 1 / (b_i - old_b_j) */
            mpc_sub(btmp, sec->bmpc[i], sec->old_bmpc[j]);

            /* If b - old_b is zero, abort the computation */
            if (mpc_eq_zero(btmp))
              {
                success = -1;

                MPS_DEBUG(s, "Cannot regenerate coefficients, reusing old ones and setting best_approx to true.")
                s->secular_equation->best_approx = true;

                goto regenerate_m_exit;
              }

            mpc_inv(ctmp, btmp);

            /* Add a_j / (b_i - old_b_j) to sec_ev */
            mpc_mul_eq(ctmp, sec->old_ampc[j]);
            mpc_add_eq(sec_ev, ctmp);

            /* Multiply prod_b for
             * b_i - b_j if i \neq j and prod_old_b
             * for b_i - old_b_i.  */
            mpc_mul_eq(prod_b, btmp);
            if (i != j)
              {
                mpc_sub(ctmp, sec->bmpc[i], sec->bmpc[j]);
                mpc_div_eq(prod_b, ctmp);
              }
          }

        /* Compute the new a_i as sec_ev * prod_old_b / prod_b */
        mpc_sub_eq_ui(sec_ev, 1, 0);
        mpc_mul(sec->ampc[i], sec_ev, prod_b);
      }


    mps_secular_raise_coefficient_precision(s, s->mpwp);


    regenerate_m_exit:

    /* Free data */
    mpc_clear(prod_b);
    mpc_clear(sec_ev);
    mpc_clear(ctmp);
    mpc_clear(btmp);

    return success;
}

/**
 * @brief Regenerate \f$a_i\f$ and \f$b_i\f$ setting
 * \f$b_i = z_i\f$, i.e. the current root approximation
 * and recomputing \f$a_i\f$ accordingly.
 */
void
mps_secular_ga_regenerate_coefficients(mps_status* s)
{
  MPS_DEBUG_THIS_CALL

  cplx_t *old_b, *old_a;
  cdpe_t *old_db, *old_da;
  mpc_t *old_ma, *old_mb;
  mps_secular_equation *sec;
  int i, j;

  sec = (mps_secular_equation*) s->secular_equation;

  switch (s->lastphase)
    {

  /* If we are in the float phase regenerate coefficients
   * starting from floating point */
  case float_phase:

    s->mpwp = 53;

    /* Allocate old_a and old_b */
    old_a = cplx_valloc(s->n);
    old_b = cplx_valloc(s->n);

    /* Copy the old coefficients, and set the new
     * b_i with the current roots approximations. */
    for (i = 0; i < s->n; i++)
      {
        cplx_set(old_a[i], sec->afpc[i]);
        cplx_set(old_b[i], sec->bfpc[i]);
        cplx_set(sec->bfpc[i], s->froot[i]);
        mpc_set_cplx(sec->bmpc[i], sec->bfpc[i]);
        mpc_set_cplx(sec->ampc[i], sec->afpc[i]);
      }

    // Regeneration
    if (mps_secular_ga_regenerate_coefficients_mp(s) != 0)
    {
        for(i = 0; i < s->n; i++)
        {
            cplx_set(sec->afpc[i], old_a[i]);
            cplx_set(sec->bfpc[i], old_b[i]);
        }
    }
    else
    {
        for(i = 0; i < s->n; i++)
            mpc_get_cplx(sec->afpc[i], sec->ampc[i]);
    }

    cplx_vfree(old_a);
    cplx_vfree(old_b);

    mps_secular_fstart(s, s->n, 0, 0, 0, s->eps_out);

    break;

    /* If this is the DPE phase regenerate DPE coefficients */
  case dpe_phase:

    s->mpwp = 53;

    /* Allocate old_a and old_b */
    old_da = cdpe_valloc(s->n);
    old_db = cdpe_valloc(s->n);

    /* Copy the old coefficients, and set the new
     * b_i with the current roots approximations. */
    for (i = 0; i < s->n; i++)
      {
        cdpe_set(old_da[i], sec->adpc[i]);
        cdpe_set(old_db[i], sec->bdpc[i]);
        cdpe_set(sec->bdpc[i], s->droot[i]);
        mpc_set_cdpe(sec->bmpc[i], sec->bdpc[i]);
      }

    // Regeneration
    if (mps_secular_ga_regenerate_coefficients_mp(s) != 0)
    {
        for(i = 0; i < s->n; i++)
        {
            cdpe_set(sec->adpc[i], old_da[i]);
            cdpe_set(sec->bdpc[i], old_db[i]);
        }
    }
    else
    {
        for(i = 0; i < s->n; i++)
            mpc_get_cdpe(sec->adpc[i], sec->ampc[i]);
    }

    /* Free data */
    cdpe_vfree(old_da);
    cdpe_vfree(old_db);

//    /* Debug new coefficients found */
//    for (i = 0; i < s->n; i++)
//      {
//        MPS_DEBUG_CDPE(s, sec->adpc[i], "sec->adpc[%d]", i);
//      }

    mps_secular_dstart(s, s->n, 0, (__rdpe_struct *) rdpe_zero,
        (__rdpe_struct *) rdpe_zero, s->eps_out);

    break;

  case mp_phase:

        /* Allocate old_a and old_b */
        old_ma = mpc_valloc(s->n);
        old_mb = mpc_valloc(s->n);

        mpc_vinit2(old_ma, s->n, 2 * s->mpwp);
        mpc_vinit2(old_mb, s->n, 2 * s->mpwp);

        /* Copy the old coefficients, and set the new
         * b_i with the current roots approximations. */
        for (i = 0; i < s->n; i++)
          {
            mpc_set(old_ma[i], sec->ampc[i]);
            mpc_set(old_mb[i], sec->bmpc[i]);
            mpc_set(sec->bmpc[i], s->mroot[i]);
          }

        // Regeneration
        if (!mps_secular_ga_regenerate_coefficients_mp(s))
        {
            MPS_DEBUG(s, "OOps");
        }


    /* This is too much verbose to be enabled even in the debugged version,
     * so do not display it (unless we are trying to catch some errors on
     * coefficient regeneration). */

//    for (i = 0; i < s->n; i++)
//      {
//        MPS_DEBUG_MPC(s, 15, sec->ampc[i], "sec->ampc[%d]", i);
//        MPS_DEBUG_MPC(s, 15, sec->bmpc[i], "sec->bmpc[%d]", i);
//      }

     mpc_vclear(old_ma, s->n);
     mpc_vclear(old_mb, s->n);
     mpc_vfree(old_ma);
     mpc_vfree(old_mb);

        mps_secular_mstart(s, s->n, 0, (__rdpe_struct *) rdpe_zero,
            (__rdpe_struct *) rdpe_zero, s->eps_out);

    break;

  default:
    break;

    } /* End of switch (s->lastphase)*/

  /* Finally set radius according to new computed a_i coefficients,
   * if they are convenient   */
  mps_secular_set_radii(s);
}


/*
 * @brief Check if iterations can terminate
 */
mps_boolean
mps_secular_ga_check_stop(mps_status* s)
{
  MPS_DEBUG_THIS_CALL

  double frad = pow(10, -s->prec_out);
  rdpe_t drad;
  int i;

  /* Set the maximum rad allowed */
  rdpe_set_2dl(drad, 1.0, -s->prec_out);

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
        if (s->frad[i] > frad)
          return false;

      /* Multiprecision and DPE case are the same, since the radii
       * are always RDPE. */
      case mp_phase:
      case dpe_phase:
        if (rdpe_log10(s->drad[i]) * LOG2_10 > -s->prec_out)
            return false;
        break;

      default:
        break;

        }
    }

  MPS_DEBUG(s, "Stop conditions were satisfied");
  return true;
}

/**
 * @brief MPSolve main function for the secular equation solving
 * using Gemignani's approach.
 */
void
mps_secular_ga_mpsolve(mps_status* s)
{
  int roots_computed = 0;
  int iteration_per_packet = 10;
  int i;
  mps_secular_equation* sec = mps_secular_equation_from_status(s);
  mps_phase phase = sec->starting_case;

  /* Set degree and allocate polynomial-related variables
   * to allow initializitation to be performed. */
  s->deg = s->n = sec->n;
  mps_status_allocate_poly_inplace(s, sec->n);
  s->data_type = "uri";

  /* We set the selected phase */
  s->lastphase = phase;

  /* Allocate other data */
  mps_allocate_data(s);

  /* Manually set FILE* pointer for streams.
   * More refined options will be added later. */
  s->outstr = s->rtstr = stdout;

  for (i = 0; i < s->n; i++)
    s->frad[i] = DBL_MAX;

  /* Set initial cluster structure as no cluster structure. */
  mps_cluster_reset(s);

  /* Copy initial coefficients so they can be used to compute */
  mps_secular_save_coefficients(s, sec);

  /* Set phase */
  s->lastphase = phase;

  /* Select initial approximations using the custom secular
   * routine and based on the phase selected by the user. */
  switch (s->lastphase)
    {
  case float_phase:
    mps_secular_fstart(s, s->n, 0, 0.0, 0.0, s->eps_out);
    break;

  case dpe_phase:
    mps_secular_dstart(s, s->n, 0, (__rdpe_struct *) rdpe_zero,
        (__rdpe_struct *) rdpe_zero, s->eps_out);
    break;

  case mp_phase:
    mps_secular_mstart(s, s->n, 0, (__rdpe_struct *) rdpe_zero,
        (__rdpe_struct *) rdpe_zero, s->eps_out);
    break;

  default:
    break;
    }

  /* Set initial radius */
  mps_secular_set_radii(s);

  /* Cycle until approximated */
  do
    {
      s->secular_equation->best_approx = false;
      /* Perform an iteration of floating point Aberth method */
      switch (s->lastphase)
        {
      case float_phase:
        roots_computed = mps_secular_ga_fiterate(s, iteration_per_packet);
        MPS_DEBUG(s, "%d roots were computed", roots_computed)
        break;

      case dpe_phase:
        roots_computed = mps_secular_ga_diterate(s, iteration_per_packet);
        MPS_DEBUG(s, "%d roots were computed", roots_computed)
        break;

      case mp_phase:
        roots_computed = mps_secular_ga_miterate(s, iteration_per_packet);
        MPS_DEBUG(s, "%d roots were computed", roots_computed)
        break;

      default:
        break;
        }

      /* Check if all roots were approximated with the
       * given input precision                      */
      if (mps_secular_ga_check_stop(s))
        return;

      /* If we can't stop recompute coefficients in higher precision and
       * continue to iterate, unless the best approximation possible in
       * this precision has been reached. In that case increase the precision
       * of the computation. */
      if (!sec->best_approx)
      {
        mps_secular_ga_regenerate_coefficients(s);
      }

      /* Instead of using else we recheck best approx because it could
       * have been set by the coefficient regeneration */
      if (sec->best_approx)
      {
          /* Going to multiprecision if we're not there yet */
          if (s->lastphase != mp_phase)
              mps_secular_switch_phase(s, mp_phase);
          else
          {
              /* Raising precision otherwise */
              mps_secular_raise_precision(s, 2 * s->mpwp);
              mps_secular_ga_regenerate_coefficients(s);
          }
      }
    }
  while (!mps_secular_ga_check_stop(s));
}
