/*
 * This file is part of MPSolve 3.2.2
 *
 * Copyright (C) 2001-2020, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Dario Andrea Bini <bini@dm.unipi.it>
 *   Giuseppe Fiorentino <fiorent@dm.unipi.it>
 *   Leonardo Robol <leonardo.robol@unipi.it>
 */

#include <mps/mps.h>
#include <math.h>
#include <string.h>

#ifndef log2
#define log2(x) (log (x) / LOG2)
#endif

/**
 * @brief Update the working precision of a root, i.e. the variable
 * <code>s->root[i]->wp</code> with the given precision rounded to the
 * closer multiple of the result of mps_context_get_minimum_precision()
 * (since GMP does not handle intermediate precisions).
 *
 * @param s The <code>mps_context</code> of the computation.
 * @param i The index of the root whose precision must be updated.
 * @param wp The precision to set.
 */
MPS_PRIVATE long int
mps_secular_ga_update_root_wp (mps_context * s, int i, long int wp, mpc_t * bmpc)
{
  mps_secular_equation * sec = s->secular_equation;
  mps_polynomial * p = s->active_poly;
  long int min_prec = mps_context_get_minimum_precision (s);

  s->root[i]->wp = ((wp - 1) / min_prec + 1) * min_prec;

  MPS_LOCK (s->data_prec_max);
  if (s->data_prec_max.value < s->root[i]->wp)
    s->data_prec_max.value = s->root[i]->wp;
  MPS_UNLOCK (s->data_prec_max);

  if (s->debug_level & MPS_DEBUG_MEMORY)
    MPS_DEBUG (s, "Setting wp for root %d to %ld bits", i, s->root[i]->wp);

  if (mpc_get_prec (sec->ampc[i]) < s->root[i]->wp)
    mpc_set_prec (sec->ampc[i], s->root[i]->wp);

  mps_polynomial_raise_data (s, p, s->root[i]->wp);

  return s->root[i]->wp;
}

/**
 * @brief This routines is used to check if a root has changed from the last regeneration,
 * in floating point phases.
 *
 * If a root is approximated or isolated and does not differ much (i.e. less than the machine
 * epsilon) from the approximation that was present a cycle ago, than it's not necessary to
 * recompute the value of the polynomial in that point, so <code>root_changed[i]</code> is set
 * to false
 *
 * @param s The <code>mps_context</code> of the computation.
 * @param old_b A vector of the old \f$b_i\f$ coefficients.
 * @param old_mb The vector of the old \f$b_i\f$ in <code>mpc_t</code> version. Must be set
 * to NULL if we are not in a multiprecision phase.
 */
static mps_boolean *
mps_secular_ga_find_changed_roots (mps_context * s, cdpe_t * old_b, mpc_t * old_mb)
{
  MPS_DEBUG_THIS_CALL (s);

  cdpe_t diff;
  int i;
  int changed_roots = 0;

  mps_secular_equation * sec = s->secular_equation;
  mps_boolean * root_changed = mps_boolean_valloc (s->n);

  mpc_t mdiff;
  if (old_mb)
    mpc_init2 (mdiff, mpc_get_prec (old_mb[0]));

  for (i = 0; i < s->n; i++)
    {
      /* FIXME: At the moment we do not have a complete analysis for the partial
       * regeneration strategy, so we use the general strategy all the time. */
      if (s->just_raised_precision || true)
        {
          root_changed[i] = true;
          continue;
        }

      /* Do multiprecision sub if we are in mp_phase, otherwise go with plain
       * CDPE. */
      if (old_mb)
        {
          mpc_sub (mdiff, old_mb[i], sec->bmpc[i]);
          mpc_get_cdpe (diff, mdiff);
          mpc_get_cdpe (sec->bdpc[i], sec->bmpc[i]);
        }
      else
        cdpe_sub (diff, old_b[i], sec->bdpc[i]);

      root_changed[i] = !cdpe_eq_zero (diff);
      if ((s->debug_level & MPS_DEBUG_REGENERATION) && !root_changed[i])
        MPS_DEBUG (s, "b_%d hasn't changed, so p(b_%d) will not be recomputed", i, i);

      if (root_changed[i])
        changed_roots++;
    }

  if (old_mb)
    mpc_clear (mdiff);

  if (changed_roots != 0)
    MPS_DEBUG (s, "%d of %d approximations are different from last regeneration", changed_roots, s->n);

  return root_changed;
}

/*! @cond PRIVATE */
struct __mps_secular_ga_regenerate_coefficients_monomial_data {
  mps_context * s;
  cdpe_t * old_b;
  mpc_t * old_mb;
  mpc_t * bmpc;
  mps_boolean * root_changed;
  rdpe_t * root_epsilon;
  mps_boolean * success;
  int i;
};
/*! @endcond */

static void *
__mps_secular_ga_regenerate_coefficients_monomial_worker (void * data_ptr)
{
  struct __mps_secular_ga_regenerate_coefficients_monomial_data * data = data_ptr;

  mps_context * s = data->s;
  mpc_t * restrict old_mb = data->old_mb;
  mpc_t * restrict bmpc = data->bmpc;
  mps_boolean * root_changed = data->root_changed;

  /* Pointers to the secular equation and the monomial_poly */
  mps_secular_equation * restrict sec = s->secular_equation;
  mps_polynomial * restrict p = s->active_poly;
  int i = data->i, j;

  /* Multiprecision variables */
  mpc_t mprod_b, ctmp, mdiff, lc, my_b;

  /* Floating point temporary variables */
  rdpe_t root_epsilon;

  /* The precision of the temporary variables at the start of the computation. We can set
   * this to s->mpwp; */
  long int coeff_wp = s->mpwp;

  /* This variable is true if the regeneration succeeded. */
  mps_boolean success = true;

  if (s->exit_required)
    return NULL;

  /* mps_secular_raise_coefficient_precision (s, coeff_wp); */

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

    default:
      break;
    }

  /* Init multiprecision values */
  mpc_init2 (mprod_b, coeff_wp);
  mpc_init2 (ctmp, coeff_wp);
  mpc_init2 (mdiff, coeff_wp);
  mpc_init2 (lc, coeff_wp);
  mpc_init2 (my_b, coeff_wp);

  mpc_set_si (lc, -1, 0);
  mps_polynomial_get_leading_coefficient (s, p, ctmp);
  mpc_div_eq (lc, ctmp);

  /*
   * The new coefficients of the secular equation can be computed
   * starting from the evaluation of the polynomial changed of
   * sign and divided for the product of the difference of the
   * b_i, i.e.:
   *
   *   a_i = -p(b_i) / \prod_{i \neq j} (b_i - b_j)
   *
   */
  if (root_changed[i])
    {
      rdpe_t relative_error, rtmp;
      cdpe_t cpol, cdiff, cprod_b;
      mpc_t tx;

      /* Set up a temporary memory location to hold the value of b_i, since we need
       * to play with its precision. This is not doable directly because it will
       * disturb other threads at work. */
      mpc_init2 (tx, s->root[i]->wp);

      /* Give a sensible minimum bound to the necessary precision */
      s->root[i]->wp = MAX (s->mpwp + log2 (s->n), s->root[i]->wp);
      mps_secular_ga_update_root_wp (s, i, s->root[i]->wp, bmpc);

      mpc_set (tx, bmpc[i]);
      mps_polynomial_meval (s, p, tx, sec->ampc[i], relative_error);

      mpc_get_cdpe (cpol, sec->ampc[i]);
      cdpe_mod (rtmp, cpol);
      rdpe_div_eq (relative_error, rtmp);

      mpc_set_si (lc, -1, 0);
      mps_polynomial_get_leading_coefficient (s, p, ctmp);
      mpc_div_eq (lc, ctmp);

      mpc_mul_eq (sec->ampc[i], lc);

      if (s->debug_level & MPS_DEBUG_REGENERATION)
        MPS_DEBUG_MPC (s, mpc_get_prec (sec->ampc[i]), sec->ampc[i], "p(b_%d)", i);

      if (s->debug_level & MPS_DEBUG_REGENERATION)
        MPS_DEBUG_RDPE (s, relative_error, "Relative_error on p(b_%d) evaluation", i);

      if (rdpe_gt (relative_error, root_epsilon))
        {
          long int required_precision = 1 + s->root[i]->wp + (rdpe_Esp (relative_error) - rdpe_Esp (root_epsilon));
	  if (s->debug_level & MPS_DEBUG_MEMORY)
	    MPS_DEBUG (s, "Setting precision to %ld", required_precision);

          /* Update the working precision of the selected root with a realistic estimate of the
           * required precision to get a result exact to machine precision */
          mps_secular_ga_update_root_wp (s, i, required_precision, bmpc);

          /* Try to recompute the polynomial with the augmented precision and see if now relative_error matches */
          mpc_set_prec (tx, s->root[i]->wp);
          mps_polynomial_meval (s, p, tx, sec->ampc[i], relative_error);

          mpc_get_cdpe (cpol, sec->ampc[i]);
          cdpe_mod (rtmp, cpol);
          rdpe_div_eq (relative_error, rtmp);

          mpc_set_si (lc, -1, 0);
          mps_polynomial_get_leading_coefficient (s, p, ctmp);
          mpc_div_eq (lc, ctmp);

          mpc_mul_eq (sec->ampc[i], lc);

          if (s->debug_level & MPS_DEBUG_REGENERATION)
            MPS_DEBUG_MPC (s, mpc_get_prec (sec->ampc[i]), sec->ampc[i], "p(b_%d)", i);

          if (s->debug_level & MPS_DEBUG_REGENERATION)
            {
              MPS_DEBUG_RDPE (s, relative_error, "Relative_error on p(b_%d) evaluation", i);
            }
        }

      if (mpc_get_prec (mprod_b) < s->root[i]->wp)
        {
          mpc_set_prec (mprod_b, s->root[i]->wp);
          mpc_set_prec (lc, s->root[i]->wp);
          mpc_set_si (lc, -1, 0);
          mps_polynomial_get_leading_coefficient (s, p, ctmp);
          mpc_div_eq (lc, ctmp);
        }

      pthread_mutex_lock (&sec->bmpc_mutex[i]);
      mpc_set (my_b, bmpc[i]);
      pthread_mutex_unlock (&sec->bmpc_mutex[i]);

      /* Compute the difference of the b_i */
      mpc_set_ui (mprod_b, 1U, 0U);

      cdpe_set (cprod_b, cdpe_one);

      for (j = 0; j < s->n; ++j)
        {
          if (i == j)
            continue;

          mpc_sub (mdiff, my_b, bmpc[j]);

          /* We'll make floating point multiplication to have
           * a faster implementation. */
          if (s->lastphase != mp_phase)
            mpc_get_cdpe (cdiff, mdiff);

          /* If the difference is zero then regeneration cannot succeed, and means
           * that we need more precision in the roots */
          if (mpc_eq_zero (mdiff))
            {
              MPS_DEBUG (s, "Regeneration of the coefficients failed because sec->bdpc[%d] == sec->bdpc[%d]", i, j);
              MPS_DEBUG_MPC (s, s->mpwp / LOG2_10 + 3, my_b, "b_%d", i);
              MPS_DEBUG_MPC (s, s->mpwp / LOG2_10 + 3, bmpc[j], "b_%d", j);
              success = false;
              goto monomial_regenerate_exit;
            }

          if (s->lastphase != mp_phase)
            cdpe_mul_eq (cprod_b, cdiff);
          else
            mpc_mul_eq (mprod_b, mdiff);
        }

      if (s->lastphase != mp_phase)
        mpc_set_cdpe (mprod_b, cprod_b);

      /* Actually divide the result and store it in
       * a_i, as requested. */
      mpc_div_eq (sec->ampc[i], mprod_b);

      /* Debug computed coefficients */
      if (s->debug_level & MPS_DEBUG_REGENERATION)
        {
          MPS_DEBUG_MPC (s, s->mpwp, sec->ampc[i], "a_%d", i);
          MPS_DEBUG_MPC (s, s->mpwp, sec->bmpc[i], "b_%d", i);
        }

      mpc_clear (tx);
    } /* Close the case where the coefficient are not approximated or isolated */
  else
    {
      mpc_set_ui (mprod_b, 1U, 0U);

      for (j = 0; j < MPS_POLYNOMIAL (sec)->degree; j++)
        {
          if (root_changed[j] && i != j)
            {
              mpc_sub (mdiff, bmpc[i], old_mb[j]);

              mpc_mul_eq (mprod_b, mdiff);

              mpc_sub (mdiff, bmpc[i], bmpc[j]);

              mpc_div_eq (mprod_b, mdiff);
            }
        }

      mpc_mul_eq (sec->ampc[i], mprod_b);
    }

monomial_regenerate_exit:
  /* Clear requested storage */
  mpc_clear (mdiff);
  mpc_clear (mprod_b);
  mpc_clear (ctmp);
  mpc_clear (lc);
  mpc_clear (my_b);
  /* mps_boolean_vfree (root_changed); */

  if (!success)
    *data->success = false;
  return NULL;
}

/**
 * @brief Compute the new secular equation coefficients based on the monomial input
 * in <code>s->monomial_poly</code>.
 *
 * @param s The <code>mps_context</code> of the computation
 * @param old_b The old \f$b_i\f$ coefficients of the secular equation used to recompute
 * the value of \f$a_i\f$ in the case where, denoting with \f$b_i\f$ the new coefficients,
 * with \f$\tilde b_i\f$ the old ones and with \f$u\f$ the machine precision:
 * \f[|b_i - \tilde b_i| < u\f]
 * so there is no need to recompute the value of \f$p(b_i)\f$.
 * @param old_mb The MP version of <code>old_b</code>, or NULL if we are not in MP.
 * @param root_changed A vector of booleans that is <code>false</code> on the components that
 * did not changed from the last regeneration.
 */
static mps_boolean
mps_secular_ga_regenerate_coefficients_monomial (mps_context * s, cdpe_t * old_b, mpc_t * old_mb, mps_boolean * root_changed)
{
  MPS_DEBUG_THIS_CALL (s);

  int i;
  mps_secular_equation * sec = s->secular_equation;
  mps_boolean success = true;

  struct __mps_secular_ga_regenerate_coefficients_monomial_data * data =
    mps_newv (struct __mps_secular_ga_regenerate_coefficients_monomial_data, s->n);

  MPS_DEBUG (s, "Regenerating coefficients from monomial input");

  for (i = s->n - 1; i >= 0; i--)
    {
      data[i].i = i;
      data[i].old_b = old_b;
      data[i].old_mb = old_mb;
      data[i].root_changed = root_changed;
      data[i].s = s;
      data[i].success = &success;
      data[i].bmpc = sec->bmpc;
      mps_thread_pool_assign (s, s->pool, __mps_secular_ga_regenerate_coefficients_monomial_worker,
                              data + i);
    }

  mps_thread_pool_wait (s, s->pool);

  free (data);

  return success;
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
 * <code>s</code> is a pointer to the <code>mps_context</code>
 * is <code>mps_secular_ga_regenerate_coefficients()</code> and
 * is the one that should be used
 *
 * @param s The mps_context of the computation.
 * @param old_b Old \f$b_i\f$ coefficients from which we are regenerating from. These are used
 * to adjust the a_i for already approximated roots.
 * @param old_mb Old \f$b_i\f$ in the multiprecision version. These are used only in the case
 * where <code>s->lastphase == mp_phase</code> to determine if the roots have changed. Otherwise
 * it must be set to NULL.
 */
mps_boolean
mps_secular_ga_regenerate_coefficients_mp (mps_context * s, cdpe_t * old_b, mpc_t * old_mb)
{
  /* Declaration and initialization of the multprecision
   * variables that are used only in that case */
  mps_boolean success = true;

  /* Get the list of the changed roots so we can compute the value of the
   * polynomial only in that approximations. */
  mps_boolean * root_changed = mps_secular_ga_find_changed_roots (s, old_b, old_mb);

  /* Regenerate the coefficients of the secular equation starting from the monomial input */
  success = mps_secular_ga_regenerate_coefficients_monomial (s, old_b, old_mb, root_changed);

  if (!success)
    MPS_DEBUG (s, "Regeneration of the coefficients failed");

  mps_boolean_vfree (root_changed);

  return success;
}

static int
__mps_compare_approximations (const void * approximation1, const void * approximation2)
{
  mps_approximation * a1 = *((mps_approximation**)approximation1);
  mps_approximation * a2 = *((mps_approximation**)approximation2);

  long int wp = mpc_get_prec (a1->mvalue);
  rdpe_t epsilon, rtmp;

  rdpe_set_2dl (epsilon, 1.0, -wp);

  int return_value = 0;
  mpc_t cmp;
  cdpe_t ccmp;

  mpc_init2 (cmp, wp);

  mpc_sub (cmp, a1->mvalue, a2->mvalue);
  mpc_get_cdpe (ccmp, cmp);

  mpc_add (cmp, a1->mvalue, a2->mvalue);
  mpc_rmod (rtmp, cmp);
  rdpe_mul_eq (epsilon, rtmp);

  rdpe_abs (rtmp, cdpe_Re (ccmp));
  if (rdpe_lt (rtmp, epsilon))
    {
      rdpe_abs (rtmp, cdpe_Im (ccmp));
      if (rdpe_lt (rtmp, epsilon))
        return_value = 0;
      else
        return_value = rdpe_lt (cdpe_Im (ccmp), rdpe_zero) ? -1 : 1;
    }
  else
    {
      return_value = rdpe_lt (cdpe_Re (ccmp), rdpe_zero) ? -2 : 2;
    }

  mpc_clear (cmp);

  return return_value;
}

/**
 * @brief Separate approximations that are too near with respect to the current precision.
 *
 * @param ctx The current mps_context.
 */
static void
mps_secular_ga_separate_approximations (mps_context * ctx)
{
  int i, cluster_base = -1;
  mpc_t perturbation;
  rdpe_t epsilon;

  /* We need to permutate the approximations in a way that
   * allow us to recognize the ones that are identical to the
   * current working precision. */
  mps_approximation ** current_approximations = mps_newv (mps_approximation *, ctx->n);

  mpc_init2 (perturbation, ctx->mpwp);

  if (ctx->lastphase == mp_phase)
    rdpe_set (epsilon, ctx->mp_epsilon);
  else
    rdpe_set_d (epsilon, DBL_EPSILON);

  for (i = 0; i < ctx->n; i++)
    {
      current_approximations[i] = ctx->root[i];
      mpc_set_prec (current_approximations[i]->mvalue, ctx->mpwp);
    }

  qsort (current_approximations, ctx->n, sizeof(mps_approximation*),
         __mps_compare_approximations);

  /* Check if we find a cluster of k equal approximations */
  for (i = 0; i < ctx->n - 1; i++)
    {
      if (cluster_base >= 0)
        {
          /* If the next approximation is different start tracking the cluster.
           * Stop event if it is equal but it's the last one. In that case increase
           * the value of i. */
          if ((__mps_compare_approximations (&current_approximations[i], &current_approximations[i + 1]) != 0)
              || (i == ctx->n - 2))
            {
              /* This is a workaround to handle the special case where the cluster ends
               * with the last approximation. */
              if (i == ctx->n - 2)
                i++;

              int k = i - cluster_base + 1;
              int j;

              if (ctx->debug_level & MPS_DEBUG_REGENERATION)
                {
                  MPS_DEBUG (ctx, "Found cluster of %d numerically equal roots: ", k);
                  for (j = 0; j < k; j++)
                    {
                      MPS_DEBUG_MPC (ctx, ctx->mpwp / LOG2_10 + 3,
                                     current_approximations[cluster_base + j]->mvalue,
                                     "Root %d in the cluster", j);
                    }
                }

              for (j = 0; j < k; j++)
                {
                  cdpe_t mod, e;
                  mpc_rmod (cdpe_Re (mod), current_approximations[cluster_base + j]->mvalue);

                  rdpe_set (cdpe_Im (mod), rdpe_zero);
                  rdpe_mul_eq (cdpe_Re (mod), epsilon);
                  rdpe_mul_eq_d (cdpe_Re (mod), 4.0);

                  /* We introduce a little rotation in the perturbation so they will
                   * be all different and disposed on a small circle around the cluster. */
                  rdpe_set_d (cdpe_Re (e), cos (1.0 * j / k * 2 * PI));
                  rdpe_set_d (cdpe_Im (e), sin (1.0 * j / k * 2 * PI));

                  cdpe_mul_eq (e, mod);
                  mpc_set_cdpe (perturbation, e);

                  if (ctx->debug_level & MPS_DEBUG_REGENERATION)
                    {
                      MPS_DEBUG (ctx, "Perturbing approximation since it is numerically equal to another one");
                      MPS_DEBUG_MPC (ctx, ctx->mpwp / LOG2_10,
                                     current_approximations[cluster_base + j]->mvalue, "Root before perturbation");
                    }

                  /* Add a small epsilon to the approximation so it will be different from
                   * the other ones in the cluster. */
                  mpc_add_eq (current_approximations[cluster_base + j]->mvalue, perturbation);

                  /* Correct the approximation radius to account the change that we have made
                   * to the approximation. */
                  rdpe_add_eq (current_approximations[cluster_base + j]->drad, cdpe_Re (mod));
                  current_approximations[cluster_base + j]->frad += rdpe_get_d (cdpe_Re (mod));

                  if (ctx->debug_level & MPS_DEBUG_REGENERATION)
                    {
                      MPS_DEBUG_MPC (ctx, ctx->mpwp / LOG2_10,
                                     current_approximations[cluster_base + j]->mvalue, "Root after perturbation");
                    }
                }

              cluster_base = -1;
            }
        }
      else
        {
          if (__mps_compare_approximations (&current_approximations[i], &current_approximations[i + 1]) == 0)
            {
              cluster_base = i;
            }
        }
    }

  free (current_approximations);
  mpc_clear (perturbation);
}

/**
 * @brief Regenerate \f$a_i\f$ and \f$b_i\f$ setting
 * \f$b_i = z_i\f$, i.e. the current root approximation
 * and recomputing \f$a_i\f$ accordingly.
 *
 * @param s The mps_context of the computation.
 */
MPS_PRIVATE mps_boolean
mps_secular_ga_regenerate_coefficients (mps_context * s)
{
  MPS_DEBUG_THIS_CALL (s);

  s->operation = MPS_OPERATION_REGENERATION;

  mps_secular_equation *sec;
  int i;
  mps_boolean successful_regeneration = true;

  sec = (mps_secular_equation*) s->secular_equation;

  MPS_DEBUG_WITH_INFO (s, "Regenerating coefficients");

  switch (s->lastphase)
    {
    case float_phase:
      for (i = 0; i < s->n; i++)
        {
          mpc_set_prec (s->root[i]->mvalue, DBL_MANT_DIG);
          mpc_set_cplx (s->root[i]->mvalue, s->root[i]->fvalue);
        }
      break;

    case dpe_phase:
      for (i = 0; i < s->n; i++)
        {
          mpc_set_prec (s->root[i]->mvalue, DBL_MANT_DIG);
          mpc_set_cdpe (s->root[i]->mvalue, s->root[i]->dvalue);
        }
      break;

    default:
      break;
    }

  mps_secular_ga_separate_approximations (s);

  /* Start timer and add execution time to the total counter */
#ifndef DISABLE_DEBUG
  clock_t *my_clock = mps_start_timer ();
#endif

  mps_regeneration_driver * driver = s->regeneration_driver;

  switch (s->lastphase)
    {
    /* If we are in the float phase regenerate coefficients
     * starting from floating point */
    case float_phase:
      driver->update_fsecular_equation (s, s->active_poly, s->root, sec);
      break;

    /* If this is the DPE phase regenerate DPE coefficients */
    case dpe_phase:
      driver->update_dsecular_equation (s, s->active_poly, s->root, sec);
      break;

    case mp_phase:
      driver->update_msecular_equation (s, s->active_poly, s->root, sec);
      break;

    default:
      break;
    }                           /* End of switch (s->lastphase) */

  mps_secular_set_radii (s);

  /* Sum execution time to the total counter */
#ifndef DISABLE_DEBUG
  s->regeneration_time += mps_stop_timer (my_clock);
#endif

  if (successful_regeneration)
    {
      MPS_DEBUG (s, "Setting again to true");
      for (i = 0; i < s->n; i++)
        s->root[i]->again = true;
    }

  return successful_regeneration;
}
