#include <mps/mps.h>

/**
 * @brief Worker for the fpolzer routine.
 */
void *
mps_thread_fpolzer_worker (void *data_ptr)
{
  mps_thread_worker_data *data = (mps_thread_worker_data *) data_ptr;
  mps_status *s = data->s;
  mps_monomial_poly *p = s->monomial_poly;
  int i, iter;
  cplx_t corr, abcorr, froot;
  double rad1, modcorr;
  mps_thread_job job;


  while (!(*data->excep) && (*data->nzeros) < s->n)
    {
      job = mps_thread_job_queue_next (s, data->queue);
      i = job.i;
      iter = job.iter;

      /* Check if we got over the maximum number of iterations */
      if (job.iter == MPS_THREAD_JOB_EXCEP)
        {
          (*data->excep) = true;
          return 0;
        }

      /* Lock this roots to make sure that we are the only one working on it */
      pthread_mutex_lock (&data->roots_mutex[i]);

      if (s->again[i])
        {

          /* Check if, while we were waiting, excep condition has been reached */
          if (*data->excep || !s->again[i] || (*data->nzeros > s->n))
            {
              pthread_mutex_unlock (&data->roots_mutex[i]);
              return 0;
            }

          (*data->it)++;

          rad1 = s->frad[i];

          /* Make a local copy of the root */
          pthread_mutex_lock (&data->aberth_mutex[i]);
          cplx_set (froot, s->froot[i]);
          pthread_mutex_unlock (&data->aberth_mutex[i]);

          if (!MPS_INPUT_CONFIG_IS_USER (s->input_config))
            {
              mps_fnewton (s, s->n, froot, &s->frad[i], corr, p->fpc, p->fap,
                           &s->again[i], false);
              if (iter == 0 && !s->again[i] && s->frad[i] > rad1 && rad1 != 0)
                s->frad[i] = rad1;
              /***************************************
               The above condition is needed to cope with the case
               where at the first iteration the starting point
               is already in the root neighbourhood and the actually
               computed radius is too big since the value of the first
               derivative is too small.
               In this case the previous radius bound, obtained by
               means of Rouche' is more reliable and strict
               **************************************/
            }
          else if (s->fnewton_usr != NULL)
            {
              (*s->fnewton_usr) (s, froot, &s->frad[i], corr, &s->again[i], NULL, false);
            }
          else
            {
              mps_fnewton_usr (s, froot, &s->frad[i], corr, &s->again[i]);
            }

          if (s->again[i] ||
              /* the correction is performed only if iter!=1 or rad(i)!=rad1 */
              MPS_INPUT_CONFIG_IS_USER (s->input_config) || iter != 0 || s->frad[i] != rad1)
            {
              mps_faberth (s, i, abcorr);

              cplx_mul_eq (abcorr, corr);
              cplx_sub (abcorr, cplx_one, abcorr);

              if (cplx_eq_zero (abcorr))
                {
                  MPS_DEBUG (s, "Aberth correction is zero");
                  cplx_set_d (abcorr, DBL_EPSILON, 0);
                }

              cplx_div (abcorr, corr, abcorr);
              cplx_sub_eq (froot, abcorr);
              modcorr = cplx_mod (abcorr);
              s->frad[i] += modcorr;

              pthread_mutex_lock (&data->aberth_mutex[i]);
              cplx_set (s->froot[i], froot);
              pthread_mutex_unlock (&data->aberth_mutex[i]);
            }

          /* check for new approximated roots */
          if (!s->again[i])
            {
              (*data->nzeros)++;
              if (*data->nzeros == s->n)
                {
                  pthread_mutex_unlock (&data->roots_mutex[i]);
                  return 0;
                }
            }
        }

      pthread_mutex_unlock (&data->roots_mutex[i]);

    }

  return NULL;
}

/**
 * @brief Drop-in replacement for the stock fpolzer routine.
 * This version adds multithread support.
 */
void
mps_thread_fpolzer (mps_status * s, int *it, mps_boolean * excep)
{
  int i, nzeros = 0, n_threads = s->n_threads;

  mps_thread_worker_data *data;
  pthread_mutex_t *aberth_mutex =
    (pthread_mutex_t *) mps_malloc (sizeof (pthread_mutex_t) * s->n);
  pthread_mutex_t *roots_mutex =
    (pthread_mutex_t *) mps_malloc (sizeof (pthread_mutex_t) * s->n);

  for (i = 0; i < s->n; i++)
    {
      pthread_mutex_init (roots_mutex + i, NULL);
      pthread_mutex_init (aberth_mutex + i, NULL);
    }

  /* Create a new job queue */
  mps_thread_job_queue *queue = mps_thread_job_queue_new (s);

  *it = 0;
  *excep = false;

  /* count the number of approximations in the root neighbourhood */
  for (i = 0; i < s->n; i++)
    if (!s->again[i])
      nzeros++;
  if (nzeros == s->n)
    {
      return;
    }

  data = (mps_thread_worker_data *) mps_malloc (sizeof (mps_thread_worker_data)
						* n_threads);

  for (i = 0; i < n_threads; i++)
    {
      data[i].it = it;
      data[i].nzeros = &nzeros;
      data[i].s = s;
      data[i].excep = excep;
      data[i].thread = i;
      data[i].n_threads = n_threads;
      data[i].aberth_mutex = aberth_mutex;
      data[i].roots_mutex = roots_mutex;
      data[i].queue = queue;
      /* pthread_create (&threads[i], NULL, &mps_thread_fpolzer_worker, */
                      /* data + i); */
      mps_thread_pool_assign (s, s->pool, mps_thread_fpolzer_worker, data + i);
    }

  mps_thread_pool_wait (s, s->pool);

  free (data);
  free (roots_mutex);
  free (aberth_mutex);
  mps_thread_job_queue_free (queue);
}

/**
 * @brief Multithread worker for mps_thread_dpolzer ()
 */
void *
mps_thread_dpolzer_worker (void *data_ptr)
{
  int iter, i;
  rdpe_t rad1, rtmp;
  cdpe_t corr, abcorr;

  /* Parse input data */
  mps_thread_worker_data *data = (mps_thread_worker_data *) data_ptr;
  mps_status *s = data->s;
  mps_monomial_poly *p = s->monomial_poly;
  mps_thread_job job;

  while (!(*data->excep) && (*data->nzeros < s->n))
    {
      job = mps_thread_job_queue_next (s, data->queue);
      i = job.i;
      iter = job.iter;

      /* Check if we got over the maximum number of iterations */
      if (job.iter == MPS_THREAD_JOB_EXCEP)
        {
          (*data->excep) = true;
          return 0;
        }

      /* Make sure that we are the only one iterating on this root */
      pthread_mutex_lock (&data->roots_mutex[i]);

      if (s->again[i])
        {
          /* Check if, while we were waiting, excep condition has been reached */
          if (*data->excep || !s->again[i] || (*data->nzeros > s->n))
            {
              pthread_mutex_unlock (&data->roots_mutex[i]);
              return 0;
            }

          (*data->it)++;
          rdpe_set (rad1, s->drad[i]);

          if (!MPS_INPUT_CONFIG_IS_USER (s->input_config))
            {
              mps_dnewton (s, s->n, s->droot[i], s->drad[i], corr, p->dpc,
                           p->dap, &s->again[i], false);
              if (iter == 0 && !s->again[i] && rdpe_gt (s->drad[i], rad1)
                  && rdpe_ne (rad1, rdpe_zero))
                rdpe_set (s->drad[i], rad1);
            }
          else if (s->dnewton_usr != NULL)
            {
              (*s->dnewton_usr) (s, s->droot[i], s->drad[i], corr,
                                 &s->again[i], NULL, false);
            }
          else
            {
              mps_dnewton_usr (s, s->droot[i], s->drad[i], corr,
                               &s->again[i]);
            }

          /************************************************
           The above condition is needed to manage with the case where
           at the first iteration the starting point is already in the
           root neighbourhood and the actually computed radius is too
           big since the value of the first derivative is too small.
           In this case the previous radius bound, obtained by means of
           Rouche' is more reliable and strict
           **********************************************/

          if (s->again[i] ||
              /* the correction is performed only if iter!=1 or rad(i)!=rad1 */
              MPS_INPUT_CONFIG_IS_USER (s->input_config) || iter != 0
              || rdpe_ne (s->drad[i], rad1))
            {
              mps_daberth (s, i, abcorr);
              cdpe_mul_eq (abcorr, corr);
              cdpe_sub (abcorr, cdpe_one, abcorr);
              if (cdpe_eq_zero (abcorr))
                {
                  MPS_DEBUG (s, "Aberth correction is zero.");
                  s->lastphase = dpe_phase;
                  cdpe_set_d (abcorr, DBL_EPSILON, 0);
                }

              cdpe_div (abcorr, corr, abcorr);
              cdpe_sub_eq (s->droot[i], abcorr);
              cdpe_mod (rtmp, abcorr);
              rdpe_add_eq (s->drad[i], rtmp);
            }

          /* check for new approximated roots */
          if (!s->again[i])
            {
              (*data->nzeros)++;
              if ((*data->nzeros) == s->n)
                {
                  pthread_mutex_unlock (&data->roots_mutex[i]);
                  return 0;
                }
            }
        }
      pthread_mutex_unlock (&data->roots_mutex[i]);
    }

  return NULL;
}

/**
 * @brief Multithread version of mps_dpolzer ().
 */
void
mps_thread_dpolzer (mps_status * s, int *it, mps_boolean * excep)
{
  mps_thread_worker_data *data;
  pthread_mutex_t *aberth_mutex, *roots_mutex;
  int i, nzeros = 0;

  /* initialize the iteration counter */
  *it = 0;
  *excep = false;

  /* count the number of approximations in the root neighbourhood */
  for (i = 0; i < s->n; i++)
    if (!s->again[i])
      nzeros++;
  if (nzeros == s->n)
    return;

  /* Prepare queue */
  mps_thread_job_queue *queue = mps_thread_job_queue_new (s);

  /* Allocate space for thread data */
  data = (mps_thread_worker_data *) mps_malloc (sizeof (mps_thread_worker_data)
                                            * s->n_threads);

  /* Allocate mutexes and init them */
  aberth_mutex = (pthread_mutex_t *) mps_malloc (sizeof (pthread_mutex_t) * s->n);
  roots_mutex = (pthread_mutex_t *) mps_malloc (sizeof (pthread_mutex_t) * s->n);
  for (i = 0; i < s->n; i++)
    {
      pthread_mutex_init (&aberth_mutex[i], NULL);
      pthread_mutex_init (&roots_mutex[i], NULL);
    }

  /* Start spawning thread */
  for (i = 0; i < s->n_threads; i++)
    {
      data[i].aberth_mutex = aberth_mutex;
      data[i].excep = excep;
      data[i].it = it;
      data[i].n_threads = s->n_threads;
      data[i].nzeros = &nzeros;
      data[i].queue = queue;
      data[i].roots_mutex = roots_mutex;
      data[i].s = s;
      data[i].thread = i;
      mps_thread_pool_assign (s, s->pool, mps_thread_dpolzer_worker, data + i);
    }

  /* Wait for the thread to complete */
  mps_thread_pool_wait (s, s->pool);

  free (aberth_mutex);
  free (roots_mutex);
  free (data);
  mps_thread_job_queue_free (queue);
}

/**
 * @brief Worker for the mpolzer routine.
 */
void *
mps_thread_mpolzer_worker (void *data_ptr)
{
  mps_thread_worker_data *data = (mps_thread_worker_data *) data_ptr;
  mps_status *s = data->s;
  mps_monomial_poly *p = s->monomial_poly;
  mps_thread_job job;
  int iter, l;
  mpc_t corr, abcorr, mroot, diff;
  rdpe_t eps, rad1, rtmp;
  cdpe_t ctmp;

  mpc_init2 (abcorr, s->mpwp);
  mpc_init2 (corr, s->mpwp);
  mpc_init2 (mroot, s->mpwp);
  mpc_init2 (diff, s->mpwp);

  rdpe_mul_d (eps, s->mp_epsilon, (double) 4 * s->n);

  /* Continue to iterate while exception condition has not
   * been reached and there more roots to approximate   */
  while ((*data->nzeros) < s->n)
    {
      /* Get next job for this thread */
      job = mps_thread_job_queue_next (s, data->queue);

      /* Set variables to be used in the rest of the code */
      iter = job.iter;

      /* Check if we exceeded the maximum number of iterations */
      if (job.iter == MPS_THREAD_JOB_EXCEP)
        {
          (*data->excep) = true;
	  goto endfun;
        }

      l = job.i;

      /* Lock roots_mutex to assure that we are the only thread
       * working on this root. Parallel computation on the same
       * root is not useful, since we would be performing the
       * same computations.                                  */
      pthread_mutex_lock (&data->roots_mutex[l]);

      /* MPS_DEBUG (s, "Iterating on root %d, iter %d", l, job.iter); */

      if (s->again[l])
        {
          /* Check if, while we were waiting, excep condition has been reached,
           * or all the zeros has been approximated.                         */
          if (*data->excep || (*data->nzeros) >= s->n)
            {
              pthread_mutex_unlock (&data->roots_mutex[l]);
              goto endfun;
            }

          /* Increment total iteration counter */
          (*data->it)++;

          /* Copy locally the root to work on */
          pthread_mutex_lock (&data->aberth_mutex[l]);
          mpc_set (mroot, s->mroot[l]);
          pthread_mutex_unlock (&data->aberth_mutex[l]);

          if (!MPS_INPUT_CONFIG_IS_USER (s->input_config))
            {
              /* sparse/dense polynomial */
              rdpe_set (rad1, s->drad[l]);
              mps_mnewton (s, s->n, mroot, s->drad[l], corr, p->mfpc,
                           p->mfppc, p->dap, p->spar, &s->again[l],
                           data->thread, false);
              if (iter == 0 && !s->again[l] && rdpe_gt (s->drad[l], rad1)
                  && rdpe_ne (rad1, rdpe_zero))
                rdpe_set (s->drad[l], rad1);

              /************************************************
               The above condition is needed to cope with the case
               where at the first iteration the starting point is
               already in the root neighbourhood and the actually
               computed radius is too big since the value of the
               first derivative is too small.
               In this case the previous radius bound, obtained by
               means of Rouche' is more reliable and strict
               ***********************************************/
            }
          else /* user's polynomial */ if (s->mnewton_usr != NULL)
            {
              (*s->mnewton_usr) (s, mroot, s->drad[l], corr, &s->again[l], NULL, false);
            }
          else
            {
              mps_mnewton_usr (s, mroot, s->drad[l], corr, &s->again[l]);
            }

          if (s->again[l] ||
              /* the correction is performed only if iter!=1 or rad[l]!=rad1 */
              MPS_INPUT_CONFIG_IS_USER (s->input_config) || iter != 0
              || rdpe_ne (s->drad[l], rad1))
            {
              /* Global lock to aberth step to reach a real Gauss-Seidel iteration */
              pthread_mutex_lock (data->global_aberth_mutex);

              /* Compute Aberth correction with locks so we can lock the
               * roots while reading them.                          */
              mps_maberth_s_wl (s, l, job.cluster_item->cluster, abcorr,
                                data->aberth_mutex);

              /* Apply aberth correction that has been computed */
              mpc_mul_eq (abcorr, corr);
              mpc_neg_eq (abcorr);
              mpc_add_eq_ui (abcorr, 1, 0);
              mpc_div (abcorr, corr, abcorr);
              mpc_sub_eq (mroot, abcorr);
              mpc_get_cdpe (ctmp, abcorr);
              cdpe_mod (rtmp, ctmp);
              rdpe_add_eq (s->drad[l], rtmp);

              /* Lock aberth_mutex and copy the computed root back
               * to its place                                   */
              pthread_mutex_lock (&data->aberth_mutex[l]);
              mpc_set (s->mroot[l], mroot);
              pthread_mutex_unlock (&data->aberth_mutex[l]);

              /* Go with others aberth iterations */
              pthread_mutex_unlock (data->global_aberth_mutex);
              sched_yield ();
            }

          /* check for new approximated roots */
          if (!s->again[l])
            {
              (*data->nzeros)++;
              if ((*data->nzeros) == s->n)
                {
                  pthread_mutex_unlock (&data->roots_mutex[l]);
                  goto endfun;
                }
            }

        }
      
      pthread_mutex_unlock (&data->roots_mutex[l]);

      /* MPS_DEBUG_MPC (s, 15, s->mroot[l], "s->mroot[%d]", l); */
      /* MPS_DEBUG_RDPE (s, s->drad[l], "s->drad[%d]", l); */

      if ((*data->nzeros) == s->n)
        {
          goto endfun;
        }
    }

endfun:                        /* free local MP variables */
  mpc_clear (corr);
  mpc_clear (abcorr);
  mpc_clear (mroot);
  mpc_clear (diff);

  return NULL;
}

/**
 * @brief Drop-in threaded replacement for the stock mpolzer.
 */
void
mps_thread_mpolzer (mps_status * s, int *it, mps_boolean * excep)
{
  int i, nzeros = 0, n_threads = s->n_threads;

  *it = 0;
  *excep = false;

  /* Check if we have already approxmiated roots */
  for (i = 0; i < s->n; i++)
    if (!s->again[i])
      nzeros++;
  if (nzeros == s->n)
    {
      return;
    }

  /* Lower the number of threads if there are a lot of approximated roots */
  if (s->n_threads > (s->n - nzeros))    
    n_threads = s->n - nzeros;    
  else    
    n_threads = s->n_threads;  

  MPS_DEBUG_WITH_INFO (s, "Spawning %d worker", n_threads);

  mps_thread_worker_data *data;

  /* Allocate and the init mutexes needed by the routine */
  pthread_mutex_t *roots_mutex =
    (pthread_mutex_t *) mps_malloc (sizeof (pthread_mutex_t) * s->n);
  pthread_mutex_t *aberth_mutex =
    (pthread_mutex_t *) mps_malloc (sizeof (pthread_mutex_t) * s->n);
  pthread_mutex_t global_aberth_mutex = PTHREAD_MUTEX_INITIALIZER;

  for (i = 0; i < s->n; i++)
    {
      pthread_mutex_init (&aberth_mutex[i], NULL);
      pthread_mutex_init (&roots_mutex[i], NULL);
    }

  /* Create a new work queue */
  mps_thread_job_queue *queue = mps_thread_job_queue_new (s);

  data = (mps_thread_worker_data *) mps_malloc (sizeof (mps_thread_worker_data)
						* n_threads);

  /* Set data to be passed to every thread and actually spawn the threads. */
  for (i = 0; i < n_threads; i++)
    {
      data[i].it = it;
      data[i].nzeros = &nzeros;
      data[i].s = s;
      data[i].excep = excep;
      data[i].thread = i;
      data[i].n_threads = n_threads;
      data[i].aberth_mutex = aberth_mutex;
      data[i].global_aberth_mutex = &global_aberth_mutex;
      data[i].queue = queue;
      data[i].roots_mutex = roots_mutex;
      mps_thread_pool_assign (s, s->pool, mps_thread_mpolzer_worker, data + i);
    }

  /* Wait for the threads to complete */
  mps_thread_pool_wait (s, s->pool);

  /* Free data and exit */
  free (data);
  for (i = 0; i < s->n; i++)
    {
      pthread_mutex_destroy (&roots_mutex[i]);
      pthread_mutex_destroy (&aberth_mutex[i]);
    }
  free (roots_mutex);
  free (aberth_mutex);
  mps_thread_job_queue_free (queue);
}
