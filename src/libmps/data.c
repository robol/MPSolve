/************************************************************
 **                                                        **
 **             __  __ ___  ___      _                     **
 **            |  \/  | _ \/ __| ___| |_ _____             **
 **            | |\/| |  _/\__ \/ _ \ \ V / -_)            **
 **            |_|  |_|_|  |___/\___/_|\_/\___|            **
 **                                                        **
 **       Multiprecision Polynomial Solver (MPSolve)       **
 **               Version 2.9, September 2012              **
 **                                                        **
 **                      Written by                        **
 **                                                        **
 **     Dario Andrea Bini       <bini@dm.unipi.it>         **
 **     Giuseppe Fiorentino     <fiorent@dm.unipi.it>      **
 **     Leonardo Robol          <robol@mail.dm.unipi.it>   **
 **                                                        **
 **           (C) 2012, Dipartimento di Matematica         **
 ***********************************************************/

#include <mps/mps.h>
#include <float.h>

/**
 * @brief Globally set the current precision of mp variables
 *
 * @param s The <code>mps_context</code> of the computation.
 * @param prec The precision that is desired for the next MP computations. 
 */
void
mps_mp_set_prec (mps_context * s, long int prec)
{
  s->mpwp = (prec / 64 + 1) * 64;
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
void
mps_allocate_data (mps_context * s)
{
  MPS_DEBUG_THIS_CALL;
  int i;

  /* s->clust = int_valloc (s->deg); */
  /* s->punt = int_valloc (s->deg + 1); */
  /* s->clust_detached = int_valloc (s->deg); */

  s->root = mps_newv (mps_approximation*, s->n);
  for (i = 0; i < s->n; i++)
    s->root[i] = mps_approximation_new (s);

  /* s->status = (char (*)[3]) char_valloc (3 * s->deg); */
  
  s->root_status = mps_newv (mps_root_status, s->deg);
  s->root_attrs  = mps_newv (mps_root_attrs,  s->deg);
  s->root_inclusion = mps_newv (mps_root_inclusion, s->deg);

  mps_cluster_reset (s);

  s->order = int_valloc (s->deg);

  /* s->fppc = cplx_valloc (s->deg + 1); */
  s->fppc1 = cplx_valloc (s->deg + 1);

  s->mfpc1 = mpc_valloc (s->deg + 1);
  for (i = 0; i <= s->deg; i++)
    mpc_init2 (s->mfpc1[i], 0);

  /* s->mfppc = mpc_valloc (s->deg + 1); */
  /* for (i = 0; i <= s->deg; i++) */
  /*   mpc_init2 (s->mfppc[i], 0); */

  s->mfppc1 = mpc_valloc (s->deg + 1);
  for (i = 0; i <= s->deg; i++)
    mpc_init2 (s->mfppc1[i], 0);

  s->mfpc2 = mpc_valloc ((s->deg + 1) * s->n_threads);
  for (i = 0; i < (s->deg + 1) * s->n_threads; i++)
    mpc_init2 (s->mfpc2[i], 0);

  /* temporary vectors */
  s->spar1 = mps_boolean_valloc (s->deg + 2);
  s->spar2 = mps_boolean_valloc ((s->deg + 2) * s->n_threads);
  s->h = mps_boolean_valloc (s->deg + 2);
  s->again_old = mps_boolean_valloc (s->deg);

  s->oldpunt = int_valloc (s->deg + 1);
  s->clust_aux = int_valloc (s->deg + 1);
  s->punt_aux = int_valloc (s->deg + 1);
  s->punt_out = int_valloc (s->deg + 1);
  s->clust_out = int_valloc (s->deg + 1);

  s->fap1 = double_valloc (s->deg + 1);
  s->fap2 = double_valloc (s->deg + 1);

  s->dap1 = rdpe_valloc (s->deg + 1);
  s->dpc1 = cdpe_valloc (s->deg + 1);
  s->dpc2 = cdpe_valloc (s->deg + 1);

  s->fradii = double_valloc (s->deg + 1);
  s->partitioning = int_valloc (s->deg + 2);
  s->dradii = rdpe_valloc (s->deg + 1);

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
 * since GMP works only with precision divisible by 64bits.
 */
long int
mps_raise_data (mps_context * s, long int prec)
{
  int i, k;
  mps_monomial_poly *p = s->monomial_poly;

  /* raise the precision of  mroot */
  for (k = 0; k < s->n; k++)
    mpc_set_prec (s->root[k]->mvalue, prec);

  if (!MPS_INPUT_CONFIG_IS_USER (s->input_config))
    {
      /* raise the precision of  mfpc */
      for (k = 0; k < s->n + 1; k++)
        if (!MPS_INPUT_CONFIG_IS_SPARSE (s->input_config) || p->spar[k])
          mpc_set_prec (p->mfpc[k], prec);

      for (i = 0; i <= s->n; i++)
        if (!MPS_INPUT_CONFIG_IS_SPARSE (s->input_config) || p->spar[i])
          {
	    if (MPS_INPUT_CONFIG_IS_REAL (s->input_config))
	      {
		if (MPS_INPUT_CONFIG_IS_RATIONAL (s->input_config) ||
		    MPS_INPUT_CONFIG_IS_INTEGER (s->input_config))
		  {
		    mpf_set_q (mpc_Re (p->mfpc[i]), p->initial_mqp_r[i]);
                    mpf_set_ui (mpc_Im (p->mfpc[i]), 0);
                    /* GMP 2.0.2 bug begin */
                    if (mpf_sgn (mpc_Re (p->mfpc[i])) !=
                        mpq_sgn (p->initial_mqp_r[i]))
                      mpf_neg (mpc_Re (p->mfpc[i]), mpc_Re (p->mfpc[i]));
                    /* GMP bug end */
		  }
	      }
	    else if (MPS_INPUT_CONFIG_IS_COMPLEX (s->input_config))
	      {
		if (MPS_INPUT_CONFIG_IS_INTEGER (s->input_config) ||
		    MPS_INPUT_CONFIG_IS_RATIONAL (s->input_config))
		  {
                    mpc_set_q (p->mfpc[i], p->initial_mqp_r[i], p->initial_mqp_i[i]);
                    /* GMP 2.0.2 bug begin */
                    if (mpf_sgn (mpc_Re (p->mfpc[i])) !=
                        mpq_sgn (p->initial_mqp_r[i]))
                      mpf_neg (mpc_Re (p->mfpc[i]), mpc_Re (p->mfpc[i]));
                    if (mpf_sgn (mpc_Im (p->mfpc[i])) !=
                        mpq_sgn (p->initial_mqp_i[i]))
                      mpf_neg (mpc_Im (p->mfpc[i]), mpc_Im (p->mfpc[i]));
                    /* GMP bug end */
		  }
              }
          }
    }

  /* Raise the precision of p' */
  if (MPS_INPUT_CONFIG_IS_SPARSE (s->input_config))
    for (k = 0; k < s->n; k++)
      if (p->spar[k + 1])
        {
          mpc_set_prec (p->mfppc[k], prec);
          mpc_mul_ui (p->mfppc[k], p->mfpc[k + 1], k + 1);
        }

  /* raise the precision of auxiliary variables */
  for (k = 0; k < s->n + 1; k++)
    {
      mpc_set_prec (s->mfpc1[k], prec);
      mpc_set_prec (s->mfppc1[k], prec);
    }

  if (MPS_INPUT_CONFIG_IS_SPARSE (s->input_config))
    for (k = 0; k < (s->n + 1) * s->n_threads; k++)
      mpc_set_prec (s->mfpc2[k], prec);

  return mpc_get_prec (s->root[0]->mvalue);
}

/**
 * @brief The same of <code>mps_raise_data()</code> but using
 * raw routines of GMP, that will not change allocations. 
 *
 * @param s The <code>mps_context</code> of the computation.
 * @param prec The desired precision.
 */
void
mps_raise_data_raw (mps_context * s, long int prec)
{
  int k;
  mps_monomial_poly *p = s->monomial_poly;

  /* raise the precision of  mroot */
  for (k = 0; k < s->n; k++)
    mpc_set_prec_raw (s->root[k]->mvalue, prec);

  /* raise the precision of  mfpc */
  if (!MPS_INPUT_CONFIG_IS_USER (s->input_config))
    for (k = 0; k < s->n + 1; k++)
      mpc_set_prec_raw (p->mfpc[k], prec);

  /* Raise the precision of sparse vectors */
  if (MPS_INPUT_CONFIG_IS_SPARSE (s->input_config))
    for (k = 0; k < s->n; k++)
      if (p->spar[k + 1])
        mpc_set_prec_raw (p->mfppc[k], prec);

  /* raise the precision of auxiliary variables */
  for (k = 0; k < s->n + 1; k++)
    {
      mpc_set_prec_raw (s->mfpc1[k], prec);
      mpc_set_prec_raw (s->mfppc1[k], prec);
    }

  if (MPS_INPUT_CONFIG_IS_SPARSE (s->input_config))
    for (k = 0; k < (s->n + 1) * s->n_threads; k++)
      mpc_set_prec_raw (s->mfpc2[k], prec);
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
void 
mps_prepare_data (mps_context * s, long int prec)
{
  MPS_DEBUG_THIS_CALL;

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
      if (MPS_INPUT_CONFIG_IS_MONOMIAL (s->input_config))
	{
	  mps_monomial_poly_raise_precision (s, s->monomial_poly, prec);
	}
      else if (MPS_INPUT_CONFIG_IS_SECULAR (s->input_config))
	mps_secular_raise_coefficient_precision (s, prec);
    }

  MPS_UNLOCK (s->data_prec_max);

  pthread_mutex_unlock (&s->precision_mutex);
}

/**
 * @brief Resets the data to the highest used precision
 *
 * @param s The <code>mps_context</code> of the computation.
 */
void
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
void
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
    }

  mps_clusterization_free (s, s->clusterization);
  free (s->root_status);
  free (s->root_attrs);
  free (s->root_inclusion);
  free (s->order);

  /* free (s->fap); */
  /* rdpe_vfree (s->dap); */

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

  for (i = 0; i < (s->deg + 1) * s->n_threads; i++)
      mpc_clear (s->mfpc2[i]);

  free (s->mfppc1);
  free (s->mfpc2);

  /* free temporary vectors */
  free (s->spar1);
  free (s->spar2);
  free (s->h);
  free (s->again_old);

  free (s->oldpunt);
  free (s->clust_aux);
  free (s->punt_aux);
  free (s->punt_out);
  free (s->clust_out);

  free (s->fap1);
  free (s->fap2);

  rdpe_vfree (s->dap1);
  cdpe_vfree (s->dpc1);
  cdpe_vfree (s->dpc2);

  free (s->partitioning);
  free (s->fradii);
  rdpe_vfree (s->dradii);

  mps_thread_pool_free (s, s->pool);
}
