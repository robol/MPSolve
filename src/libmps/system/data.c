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
#include <float.h>

/**
 * @brief Globally set the current precision of mp variables
 *
 * @param s The <code>mps_context</code> of the computation.
 * @param prec The precision that is desired for the next MP computations.
 */
MPS_PRIVATE void
mps_mp_set_prec (mps_context * s, long int prec)
{
  long int min_prec = mps_context_get_minimum_precision (s);

  s->mpwp = (prec / min_prec + 1) * min_prec;
  rdpe_set_2dl (s->mp_epsilon, 1.0, -s->mpwp);

  if (s->debug_level & MPS_DEBUG_MEMORY)
    MPS_DEBUG_RDPE (s, s->mp_epsilon, "Increased precision to %ld bits. Machine epsilon set to eps", s->mpwp);
}

/**
 * @brief Allocate all the data needed by MPSolve. Must be called after setting
 * the degree of the polynomial (or, more generally, the number of root of the
 * equation) in <code>s->deg</code>.
 *
 * @param s The <code>mps_context</code> of the computation.
 */
MPS_PRIVATE void
mps_allocate_data (mps_context * s)
{
  MPS_DEBUG_THIS_CALL (s);
  int i;

  if (s->initialized)
    return;

  s->root = mps_newv (mps_approximation *, s->n);
  for (i = 0; i < s->n; i++)
    s->root[i] = mps_approximation_new (s);

  /* Reset the cluster structure, so we can start without assumption
   * on the location of the roots. */
  mps_cluster_reset (s);

  s->order = int_valloc (s->deg);

  s->fppc1 = cplx_valloc (s->deg + 1);

  s->mfpc1 = mpc_valloc (s->deg + 1);
  for (i = 0; i <= s->deg; i++)
    mpc_init2 (s->mfpc1[i], 0);

  s->mfppc1 = mpc_valloc (s->deg + 1);
  for (i = 0; i <= s->deg; i++)
    mpc_init2 (s->mfppc1[i], 0);

  /* temporary vectors */
  s->spar1 = mps_boolean_valloc (s->deg + 2);
  s->again_old = mps_boolean_valloc (s->deg);

  s->fap1 = double_valloc (s->deg + 1);
  s->fap2 = double_valloc (s->deg + 1);

  s->dap1 = rdpe_valloc (s->deg + 1);
  s->dpc1 = cdpe_valloc (s->deg + 1);
  s->dpc2 = cdpe_valloc (s->deg + 1);

  /* Setting some default here, that were not settable because we didn't know
   * the degree of the polynomial */
  for (i = 0; i < s->n; i++)
    s->root[i]->wp = DBL_DIG * LOG2_10;

  /* Init the mutex that need it */
  pthread_mutex_init (&s->precision_mutex, NULL);

  /* Other mt variables in status */
  MPS_INIT_LOCK (s->data_prec_max);

  s->initialized = true;
}

/**
 * @brief Raise precision performing a real computation of the data.
 *
 * @param s The <code>mps_context</code> of the computation.
 * @param prec The desired precision.
 * @return The precision set (that may be different from the one requested
 * since GMP works only with precision divisible by some given integer.
 */
MPS_PRIVATE long int
mps_raise_data (mps_context * s, long int prec)
{
  int k;
  mps_polynomial *p = s->active_poly;

  /* raise the precision of  mroot */
  for (k = 0; k < s->n; k++)
    mpc_set_prec (s->root[k]->mvalue, prec);

  /* raise the precision of auxiliary variables */
  for (k = 0; k < s->n + 1; k++)
    {
      mpc_set_prec (s->mfpc1[k], prec);
      mpc_set_prec (s->mfppc1[k], prec);
    }

  mps_polynomial_raise_data (s, p, prec);

  return mpc_get_prec(s->root[0]->mvalue);
}

/**
 * @brief The same of <code>mps_raise_data()</code> but using
 * raw routines of GMP, that will not change allocations.
 *
 * @param s The <code>mps_context</code> of the computation.
 * @param prec The desired precision.
 */
MPS_PRIVATE void
mps_raise_data_raw (mps_context * s, long int prec)
{
  int k;

  if (!MPS_IS_MONOMIAL_POLY (s->active_poly))
    return;

  mps_monomial_poly *p = MPS_MONOMIAL_POLY (s->active_poly);

  /* raise the precision of  mroot */
  for (k = 0; k < s->n; k++)
    mpc_set_prec_raw (s->root[k]->mvalue, prec);

  /* raise the precision of  mfpc */
  if (MPS_IS_MONOMIAL_POLY (s->active_poly))
    for (k = 0; k < s->n + 1; k++)
      mpc_set_prec_raw (p->mfpc[k], prec);

  /* Raise the precision of sparse vectors */
  if (MPS_DENSITY_IS_SPARSE (s->active_poly->density))
    for (k = 0; k < s->n; k++)
      if (p->spar[k + 1])
        mpc_set_prec_raw (p->mfppc[k], prec);

  /* raise the precision of auxiliary variables */
  for (k = 0; k < s->n + 1; k++)
    {
      mpc_set_prec_raw (s->mfpc1[k], prec);
      mpc_set_prec_raw (s->mfppc1[k], prec);
    }
}

/**
 * @brief Compute the mp_complex values of the coefficients of p(x)
 * with the  current precision of mpwds words, given the
 * rational or integer coefficients.
 *
 * @param s The <code>mps_context</code> of the computation.
 * @param prec The precision that should be set and to which the data should
 * be adjusted.
 */
MPS_PRIVATE void
mps_prepare_data (mps_context * s, long int prec)
{
  MPS_DEBUG_THIS_CALL (s);

  pthread_mutex_lock (&s->precision_mutex);

  if (s->debug_level & MPS_DEBUG_MEMORY)
    MPS_DEBUG (s, "Increasing working precision to %ld bits", prec);

  MPS_LOCK (s->data_prec_max);

  if (prec > s->data_prec_max.value)
    {
      s->data_prec_max.value = mps_raise_data (s, prec);
    }
  else
    {
      mps_polynomial_raise_data (s, s->active_poly, prec);
    }

  MPS_UNLOCK (s->data_prec_max);

  pthread_mutex_unlock (&s->precision_mutex);
}

/**
 * @brief Resets the data to the highest used precision
 *
 * @param s The <code>mps_context</code> of the computation.
 */
MPS_PRIVATE void
mps_restore_data (mps_context * s)
{
  MPS_LOCK (s->data_prec_max);
  if (s->debug_level & MPS_DEBUG_MEMORY)
    MPS_DEBUG (s, "Restore data to %ld bits", s->data_prec_max.value);

  if (s->data_prec_max.value)
    {
      MPS_UNLOCK (s->data_prec_max);
      mps_raise_data_raw (s, s->data_prec_max.value);
    }
  else
    MPS_UNLOCK (s->data_prec_max);
}

/**
 * @brief Free all the data allocated with <code>mps_allocate_data()</code>
 *
 * @param s The <code>mps_context</code> of the computation.
 */
MPS_PRIVATE void
mps_free_data (mps_context * s)
{
  int i;

  if (s->debug_level & MPS_DEBUG_MEMORY)
    {
      MPS_DEBUG (s, "Deallocating data");
    }

  if (s->bmpc)
    {
      mpc_vclear (s->bmpc, s->n * s->pool->n);
      free (s->bmpc);
      s->bmpc = NULL;
    }

  /* Release our reference to any active polynomial without
   * freeing it, since that is responsability of the user. 
   *
   * Note: This is really important since this context may be
   * recycled at a later time, and having a pointer to a possibly
   * not-anymore-valid polynomial will cause a lot of issues. */
  if (s->active_poly)
     s->active_poly = NULL;

  mps_clusterization_free (s, s->clusterization);
  s->clusterization = NULL;

  free (s->order);

  for (i = 0; i < s->n; i++)
    mps_approximation_free (s, s->root[i]);
  free (s->root);

  for (i = 0; i <= s->deg; i++)
    mpc_clear (s->mfpc1[i]);
  mpc_vfree (s->mfpc1);

  cplx_vfree (s->fppc1);
  for (i = 0; i <= s->deg; i++)
    {
      mpc_clear (s->mfppc1[i]);
    }

  pthread_mutex_destroy(&s->precision_mutex);

  free (s->mfppc1);

  /* free temporary vectors */
  free (s->spar1);
  free (s->again_old);

  free (s->fap1);
  free (s->fap2);

  rdpe_vfree (s->dap1);
  cdpe_vfree (s->dpc1);
  cdpe_vfree (s->dpc2);
}
