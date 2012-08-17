/************************************************************
 **                                                        **
 **             __  __ ___  ___      _                     **
 **            |  \/  | _ \/ __| ___| |_ _____             **
 **            | |\/| |  _/\__ \/ _ \ \ V / -_)            **
 **            |_|  |_|_|  |___/\___/_|\_/\___|            **
 **                                                        **
 **       Multiprecision Polynomial Solver (MPSolve)       **
 **                 Version 2.9, April 2011                **
 **                                                        **
 **                      Written by                        **
 **                                                        **
 **     Dario Andrea Bini       <bini@dm.unipi.it>         **
 **     Giuseppe Fiorentino     <fiorent@dm.unipi.it>      **
 **     Leonardo Robol          <robol@mail.dm.unipi.it>   **
 **                                                        **
 **           (C) 2011, Dipartimento di Matematica         **
 ***********************************************************/

#include <mps/mps.h>
#include <math.h>
#include <string.h>

/* Routine called to perform the floating point newton iterations. */
void *
__mps_secular_ga_fiterate_worker (void* data_ptr)
{
  mps_thread_worker_data *data = (mps_thread_worker_data *) data_ptr;
  mps_status *s = data->s;
  mps_secular_equation *sec = s->secular_equation;
  int i;
  cplx_t corr, abcorr;
  double modcorr;
  mps_thread_job job;

  mps_secular_iteration_data it_data;

  it_data.local_afpc = cplx_valloc (s->n);
  it_data.local_bfpc = cplx_valloc (s->n);
  for (i = 0; i < s->n; i++)
    {
      cplx_set (it_data.local_afpc[i], sec->afpc[i]);
      cplx_set (it_data.local_bfpc[i], sec->bfpc[i]);
    }

  while ((*data->nzeros < s->n) && !*data->excep)
    {
      job = mps_thread_job_queue_next (s, data->queue);
      i = job.i;

      if (job.iter == MPS_THREAD_JOB_EXCEP || *data->nzeros >= s->n)
	goto cleanup;

      pthread_mutex_lock (&data->roots_mutex[i]);

      if (job.iter == MPS_THREAD_JOB_EXCEP || *data->nzeros >= s->n)
	{
	  pthread_mutex_unlock (&data->roots_mutex[i]);
	  goto cleanup;
	}

      if (s->again[i])
	{
	  /* Increment the number of performed iterations */
#if defined(__GCC__)
	  __sync_add_and_fetch (data->it, 1);
#else
	  pthread_mutex_lock (data->gs_mutex);
	  (*data->it)++;
	  pthread_mutex_unlock (data->gs_mutex);
#endif

	  it_data.k = i;
	  it_data.gs_mutex = data->gs_mutex;
	  // MPS_DEBUG_CPLX (s, s->froot[i], "s->froot[%d]", i);
	  cdpe_set_x (s->droot[i], s->froot[i]);
	  mps_secular_fnewton (s, s->froot[i], &s->frad[i], corr,
			       &s->again[i], &it_data, false);

	  if (s->root_status[i] == MPS_ROOT_STATUS_NOT_FLOAT)
	    {
	      *data->excep = true;
	      break;
	    }

	  /* Apply Aberth correction */
	  mps_faberth_wl (s, i, abcorr, data->aberth_mutex);
	  cplx_mul_eq (abcorr, corr);
	  cplx_sub (abcorr, cplx_one, abcorr);
	  cplx_div (abcorr, corr, abcorr);

	  if (isnan (cplx_Re (s->froot[i])) || isnan (cplx_Im (s->froot[i])))
	    {
	      *data->excep = true;
	      cdpe_get_x (s->froot[i], s->droot[i]);
	      break;
	    }

	  pthread_mutex_lock (&data->aberth_mutex[i]);
	  cplx_sub_eq (s->froot[i], abcorr);
	  pthread_mutex_unlock (&data->aberth_mutex[i]);

	  /* Correct the radius */
	  modcorr = cplx_mod (abcorr);
	  s->frad[i] += modcorr;

	  if (!s->again[i])
	    {
	      if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)
		MPS_DEBUG (s, "Root %d again was set to false on iteration %d by thread %d", i, *data->it, data->thread);

#if defined( __GCC__)
	      __sync_add_and_fetch (data->nzeros, 1);
#else
	      pthread_mutex_lock (data->gs_mutex);
	      (*data->nzeros)++;
	      pthread_mutex_unlock (data->gs_mutex);
#endif
	    }
	}

      pthread_mutex_unlock (&data->roots_mutex[i]);
    }

 cleanup:

  cplx_vfree (it_data.local_afpc);
  cplx_vfree (it_data.local_bfpc);

  return NULL;
}

/**
 * @brief Routine that performs a block of iteration
 * in floating point on the secular equation.
 *
 * @param s the pointer to the mps_status struct.
 * @param maxit Maximum number of iteration to perform.
 * @param just_regenerated true if this is the first iteration after a coefficient
 * regeneration. If just_regenerated is true and the iteration packet is completed
 * in less than 2 * (n - computed_roots) iterations that best_approx is set to true
 * in s->secular_equation so a raise in the precision will be triggered.
 * @return The number of approximated roots after the iteration.
 */
int
mps_secular_ga_fiterate (mps_status * s, int maxit, mps_boolean just_regenerated)
{
  int computed_roots = 0;
  int i;
  int nit = 0;
  int it_threshold;

  mps_boolean excep = false;

#ifndef DISABLE_DEBUG
  clock_t *my_clock = mps_start_timer ();
#endif

  mps_secular_equation *sec = s->secular_equation;

  s->operation = MPS_OPERATION_ABERTH_FP_ITERATIONS;

  mps_thread_worker_data *data;
  pthread_mutex_t *aberth_mutex =
    (pthread_mutex_t *) mps_malloc (sizeof (pthread_mutex_t) * s->n);
  pthread_mutex_t *roots_mutex =
    (pthread_mutex_t *) mps_malloc (sizeof (pthread_mutex_t) * s->n);

  pthread_mutex_t gs_mutex = PTHREAD_MUTEX_INITIALIZER;

  for (i = 0; i < s->n; i++)
    {
      pthread_mutex_init (roots_mutex + i, NULL);
      pthread_mutex_init (aberth_mutex + i, NULL);
    }

  data = mps_newv (mps_thread_worker_data, s->n_threads);

  MPS_DEBUG_THIS_CALL;

  sec->best_approx = false;

  /* Mark the approximated roots as ready for output */
  for (i = 0; i < s->n; i++)
    {
      /* Set again to false if the root is already approximated */
      if (MPS_ROOT_STATUS_IS_COMPUTED (s, i))
	{
	  if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)
	    {
	      MPS_DEBUG_WITH_INFO (s, "Setting again[%d] to false since the root is ready for output (or isolated)", i);
	    }
	  s->again[i] = false;
	}

      if (!s->again[i])
        computed_roots++;
    }

  /* Set the iterations threshold to 2 iterations
   * for every non approximated root. */
  it_threshold = 2 * (s->n - computed_roots);

  if (s->debug_level & MPS_DEBUG_PACKETS)
    {
      MPS_DEBUG (s, "There are %d roots with again set to false", computed_roots);
      MPS_DEBUG (s, "Iteration theshold set to %d iterations", it_threshold);
    }

  mps_thread_job_queue *queue = mps_thread_job_queue_new (s);

  for (i = 0; i < s->n_threads; i++)
    {
      data[i].it = &nit;
      data[i].nzeros = &computed_roots;
      data[i].s = s;
      data[i].thread = i;
      data[i].n_threads = s->n_threads;
      data[i].aberth_mutex = aberth_mutex;
      data[i].roots_mutex = roots_mutex;
      data[i].queue = queue;
      data[i].gs_mutex = &gs_mutex;
      data[i].excep = &excep;

       mps_thread_pool_assign (s, s->pool, __mps_secular_ga_fiterate_worker, data + i); 
      /* __mps_secular_ga_fiterate_worker (data + i); */
    }

  mps_thread_pool_wait (s, s->pool);

  /* Check if the roots are improvable in floating point */
  MPS_DEBUG_WITH_INFO (s, "Performed %d iterations with floating point arithmetic",
                       nit);

  if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)
      mps_dump (s);

  if (nit <= it_threshold && just_regenerated)
    {
      MPS_DEBUG_WITH_INFO (s, "Asking for regeneration since we stopped after a few iterations");
      s->secular_equation->best_approx = true;
    }

  if (excep)
    {
      MPS_DEBUG_WITH_INFO (s, "Switching to DPE arithmetic since there are roots not representable in standard floating point");
      for (i = 0; i < s->n; i++)
	{
	  cdpe_set_x (s->droot[i], s->froot[i]);
	  rdpe_set_d (s->drad[i], s->frad[i]);
	  s->root_status[i] = MPS_ROOT_STATUS_CLUSTERED;
	}
      s->lastphase = dpe_phase;
    }

  /* Compute the inclusion radii with Gerschgorin so we can compute
   * clusterizations for the roots. */
  /* mps_fradii (s, fradii); */
  /* mps_fcluster (s, fradii, 2.0 * s->n);  */
  /* mps_fmodify (s, false);  */

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

  if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)
    mps_dump (s);

  /* Count time taken  */
#ifndef DISABLE_DEBUG
  s->fp_iteration_time += mps_stop_timer (my_clock);
#endif

  mps_thread_job_queue_free (queue);
  free (data);
  free (roots_mutex);
  free (aberth_mutex);

  /* Return the number of approximated roots */
  return computed_roots;
}

/* Routine called to perform the floating point newton iterations with DPE */
void *
__mps_secular_ga_diterate_worker (void* data_ptr)
{
  mps_thread_worker_data *data = (mps_thread_worker_data *) data_ptr;
  mps_status *s = data->s;
  /* mps_secular_equation *sec = s->secular_equation; */
  int i;
  cdpe_t corr, abcorr, droot;
  rdpe_t modcorr;
  mps_thread_job job;

  mps_secular_iteration_data it_data;

  while ((*data->nzeros < s->n))
    {
      job = mps_thread_job_queue_next (s, data->queue);
      i = job.i;

      if (job.iter == MPS_THREAD_JOB_EXCEP)
        {
          return NULL;
        }

      pthread_mutex_lock (&data->roots_mutex[i]);

      if (s->again[i])
	{
          /* Lock this roots to make sure that we are the only one working on it */
	  cdpe_set (droot, s->droot[i]);

	  (*data->it)++;

	  it_data.k = i;
	  mps_secular_dnewton (s, droot, s->drad[i], corr,
			       &s->again[i], &it_data, false);

	  /* Apply Aberth correction */
	  mps_daberth_wl (s, i, abcorr, data->aberth_mutex);
	  cdpe_mul_eq (abcorr, corr);
	  cdpe_sub (abcorr, cdpe_one, abcorr);
	  cdpe_div (abcorr, corr, abcorr);
	  cdpe_sub_eq (droot, abcorr);

	  /* Correct the radius */
	  cdpe_mod (modcorr, abcorr);
	  rdpe_add_eq (s->drad[i], modcorr);

	  if (!s->again[i])
	    {
	      if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)
		MPS_DEBUG (s, "Root %d again was set to false on iteration %d by thread %d", i, *data->it, data->thread);
	      (*data->nzeros)++;
	    }

	  cdpe_set (s->droot[i], droot);
	}

      pthread_mutex_unlock (&data->roots_mutex[i]);
    }

  return NULL;
}

/**
 * @brief Routine that performs a block of iteration
 * in floating point on the secular equation using
 * CDPE
 *
 * @param s the pointer to the mps_status struct.
 * @param maxit Maximum number of iteration to perform.
 * @return The number of approximated roots after the iteration.
 * @param just_regenerated true if this is the first iteration after a coefficient
 * regeneration. If just_regenerated is true and the iteration packet is completed
 * in less than 2 * (n - computed_roots) iterations that best_approx is set to true
 * in s->secular_equation so a raise in the precision will be triggered.
 */
int
mps_secular_ga_diterate (mps_status * s, int maxit, mps_boolean just_regenerated)
{
  int computed_roots = 0;
  int i;
  int nit = 0;
  int it_threshold;

  s->operation = MPS_OPERATION_ABERTH_DPE_ITERATIONS;

#ifndef DISABLE_DEBUG
  clock_t *my_clock = mps_start_timer ();
#endif

  mps_secular_equation *sec = s->secular_equation;

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

  data = mps_newv (mps_thread_worker_data, s->n_threads);

  MPS_DEBUG_THIS_CALL;

  sec->best_approx = false;

  /* Mark the approximated roots as ready for output */
  for (i = 0; i < s->n; i++)
    {
      /* Set again to false if the root is already approximated */
      if (MPS_ROOT_STATUS_IS_COMPUTED (s, i))
	{
	  if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)
	    {
	      MPS_DEBUG_WITH_INFO (s, "Setting again[%d] to false since the root is ready for output (or isolated)", i);
	    }
	  s->again[i] = false;
	}

      if (!s->again[i])
        computed_roots++;
    }

  /* Set the iterations threshold to 2 iterations
   * for every non approximated root. */
  it_threshold = 2 * (s->n - computed_roots);

  if (s->debug_level & MPS_DEBUG_PACKETS)
    {
      MPS_DEBUG (s, "There are %d roots with again set to false", computed_roots);
      MPS_DEBUG (s, "Iteration theshold set to %d iterations", it_threshold);
    }

  mps_thread_job_queue *queue = mps_thread_job_queue_new (s);

  for (i = 0; i < s->n_threads; i++)
    {
      data[i].it = &nit;
      data[i].nzeros = &computed_roots;
      data[i].s = s;
      data[i].thread = i;
      data[i].n_threads = s->n_threads;
      data[i].aberth_mutex = aberth_mutex;
      data[i].roots_mutex = roots_mutex;
      data[i].queue = queue;

      mps_thread_pool_assign (s, s->pool, __mps_secular_ga_diterate_worker, data + i); 
    }

  mps_thread_pool_wait (s, s->pool);

  /* Check if the roots are improvable in floating point */
  MPS_DEBUG_WITH_INFO (s, "Performed %d iterations with CDPE arithmetic",
                       nit);

  if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)
      mps_dump (s);

  if (nit <= it_threshold && just_regenerated)
    {
      MPS_DEBUG_WITH_INFO (s, "Asking for regeneration since we stopped after a few iterations");
      s->secular_equation->best_approx = true;
    }

  /* Compute the inclusion radii with Gerschgorin so we can compute
   * clusterizations for the roots. */
  /* mps_dradii (s, dradii); */
  /* mps_dcluster (s, dradii, 2.0 * s->n);  */
  /* mps_dmodify (s, false);  */

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

  mps_thread_job_queue_free (queue);
  free (aberth_mutex);
  free (roots_mutex);

  free (data);

  /* Return the number of approximated roots */
  return computed_roots;
}

/* Routine called to perform the floating point newton iterations with MP */
void *
__mps_secular_ga_miterate_worker (void* data_ptr)
{
  mps_thread_worker_data *data = (mps_thread_worker_data *) data_ptr;
  mps_status *s = data->s;
  /* mps_secular_equation *sec = s->secular_equation; */
  int i;
  mpc_t corr, abcorr;
  mpc_t mroot; 
  rdpe_t modcorr;
  mps_thread_job job;

  mps_secular_iteration_data it_data;
  mps_cluster * cluster = NULL;

  mpc_init2 (corr, s->mpwp);
  mpc_init2 (abcorr, s->mpwp);
  mpc_init2 (mroot, s->mpwp); 

  /* Get a copy of the MP coefficients that is local to this thread */
  it_data.local_ampc = mpc_valloc (s->n);
  it_data.local_bmpc = mpc_valloc (s->n);
  mpc_vinit2 (it_data.local_ampc, s->n, s->mpwp);
  mpc_vinit2 (it_data.local_bmpc, s->n, s->mpwp);
  for (i = 0; i < s->n; i++)
    {
      mpc_set (it_data.local_ampc[i], s->secular_equation->ampc[i]);
      mpc_set (it_data.local_bmpc[i], s->secular_equation->bmpc[i]);
    }

  while ((*data->nzeros < s->n))
    {
      job = mps_thread_job_queue_next (s, data->queue);
      i = job.i;

      if (job.iter == MPS_THREAD_JOB_EXCEP || *data->nzeros >= s->n)
	goto cleanup;

      pthread_mutex_lock (&data->roots_mutex[i]);

      if (job.iter == MPS_THREAD_JOB_EXCEP || *data->nzeros >= s->n)
	{
	  pthread_mutex_unlock (&data->roots_mutex[i]);
	  goto cleanup;
	}

      /* printf ("Thread %d iterating on root %d\n", data->thread, i); */

      cluster = job.cluster_item->cluster;

      if (s->again[i])
	{
          /* Lock this roots to make sure that we are the only one working on it */
	  pthread_mutex_lock (&data->aberth_mutex[i]); 
	  mpc_set (mroot, s->mroot[i]); 
	  pthread_mutex_unlock (&data->aberth_mutex[i]); 

          /* Check if, while we were waiting, excep condition has been reached,
           * or all the zeros has been approximated.                         */
          if ((*data->nzeros) >= s->n)
            {
              pthread_mutex_unlock (&data->roots_mutex[i]);
              goto cleanup;
            }

	  /* pthread_mutex_lock (data->gs_mutex); */
	  (*data->it)++;
	  /* pthread_mutex_unlock (data->gs_mutex); */

	  it_data.k = i;
	  it_data.gs_mutex = data->gs_mutex;
	  mps_secular_mnewton (s, mroot, s->drad[i], corr,
			       &s->again[i], &it_data, false);

	  /* Apply Aberth correction */
	  mps_maberth_s_wl (s, i, cluster, abcorr, data->aberth_mutex);
	  mpc_mul_eq (abcorr, corr);
	  mpc_ui_sub (abcorr, 1U, 0U, abcorr); 
	  
	  if (!mpc_eq_zero (abcorr)) 
	    {
	      mpc_div (abcorr, corr, abcorr); 

	      pthread_mutex_lock (&data->aberth_mutex[i]);
	      mpc_sub_eq (mroot, abcorr); 
	      pthread_mutex_unlock (&data->aberth_mutex[i]);
	    } 
	  else
	    s->again[i] = true;

	  pthread_mutex_lock (&data->aberth_mutex[i]); 
	  mpc_set (s->mroot[i], mroot); 
	  pthread_mutex_unlock (&data->aberth_mutex[i]); 
	   
	  /* Correct the radius */
	  mpc_rmod (modcorr, abcorr);
	  rdpe_add_eq (s->drad[i], modcorr);

	  if (!s->again[i])
	    {
	      if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)
		MPS_DEBUG (s, "Root %d again was set to false on iteration %d by thread %d", i, *data->it, data->thread);

	      (*data->nzeros)++;
	    }
	}

      pthread_mutex_unlock (&data->roots_mutex[i]);
    }

 cleanup:
  mpc_clear (mroot);
  mpc_clear (abcorr);
  mpc_clear (corr);

  mpc_vclear (it_data.local_bmpc, s->n);
  mpc_vclear (it_data.local_ampc, s->n);
  mpc_vfree (it_data.local_ampc);
  mpc_vfree (it_data.local_bmpc);

  return NULL;
}

/**
 * @brief Routine that performs a block of iteration
 * in floating point on the secular equation using
 * CDPE
 *
 * @param s the pointer to the mps_status struct.
 * @param maxit Maximum number of iteration to perform. 
 * @param just_regenerated true if this is the first iteration after a coefficient
 * regeneration. If just_regenerated is true and the iteration packet is completed
 * in less than 2 * (n - computed_roots) iterations that best_approx is set to true
 * in s->secular_equation so a raise in the precision will be triggered.
 * @return The number of approximated roots after the iteration.
 */
int
mps_secular_ga_miterate (mps_status * s, int maxit, mps_boolean just_regenerated)
{
  int computed_roots = 0;
  int i;
  int nit = 0;
  int it_threshold;

  s->operation = MPS_OPERATION_ABERTH_MP_ITERATIONS;

#ifndef DISABLE_DEBUG
  clock_t *my_clock = mps_start_timer ();
#endif

  mps_secular_equation *sec = s->secular_equation;

  mps_thread_worker_data *data;
  pthread_mutex_t *aberth_mutex =
    (pthread_mutex_t *) mps_malloc (sizeof (pthread_mutex_t) * s->n);
  pthread_mutex_t *roots_mutex =
    (pthread_mutex_t *) mps_malloc (sizeof (pthread_mutex_t) * s->n);

  pthread_mutex_t gs_mutex=PTHREAD_MUTEX_INITIALIZER;

  for (i = 0; i < s->n; i++)
    {
      pthread_mutex_init (roots_mutex + i, NULL);
      pthread_mutex_init (aberth_mutex + i, NULL);
    }

  data = mps_newv (mps_thread_worker_data, s->n_threads);

  MPS_DEBUG_THIS_CALL;

  sec->best_approx = false;

  /* Mark the approximated roots as ready for output */
  for (i = 0; i < s->n; i++)
    {
      /* Set again to false if the root is already approximated */
      if (MPS_ROOT_STATUS_IS_COMPUTED (s, i))
	{
	  if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)
	    {
	      MPS_DEBUG_WITH_INFO (s, "Setting again[%d] to false since the root is ready for output (or isolated)", i);
	    }
	  s->again[i] = false;
	}

      if (!s->again[i])
        computed_roots++;
    }

  /* Set the iterations threshold to 2 iterations
   * for every non approximated root. */
  it_threshold = 2 * (s->n - computed_roots);

  if (s->debug_level & MPS_DEBUG_PACKETS)
    {
      MPS_DEBUG (s, "There are %d roots with again set to false", computed_roots);
      MPS_DEBUG (s, "Iteration theshold set to %d iterations", it_threshold);
    }

  mps_thread_job_queue *queue = mps_thread_job_queue_new (s);

  for (i = 0; i < s->n_threads; i++)
    {
      data[i].it = &nit;
      data[i].nzeros = &computed_roots;
      data[i].s = s;
      data[i].thread = i;
      data[i].n_threads = s->n_threads;
      data[i].aberth_mutex = aberth_mutex;
      data[i].roots_mutex = roots_mutex;
      data[i].queue = queue;
      data[i].gs_mutex = &gs_mutex;

      mps_thread_pool_assign (s, s->pool, __mps_secular_ga_miterate_worker, data + i); 
    }

  mps_thread_pool_wait (s, s->pool);

  /* Check if the roots are improvable in floating point */
  MPS_DEBUG_WITH_INFO (s, "Performed %d iterations with MP arithmetic",
                       nit);

  if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)
      mps_dump (s);

  if (nit <= it_threshold && just_regenerated)
    {
      MPS_DEBUG_WITH_INFO (s, "Asking for regeneration since we stopped after a few iterations");
      s->secular_equation->best_approx = true;
    }

  /* Compute the inclusion radii with Gerschgorin so we can compute
   * clusterizations for the roots. */
  /* mps_mradii (s, dradii); */
  /* mps_mcluster (s, dradii, 2.0 * s->n);  */
  /* mps_mmodify (s, false);  */

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
  s->mp_iteration_time += mps_stop_timer (my_clock);
#endif

  mps_thread_job_queue_free (queue);
  free (aberth_mutex);
  free (roots_mutex);

  free (data);

  /* Return the number of approximated roots */
  return computed_roots;
}


/* /\** */
/*  * @brief Routine that performs a block of iteration */
/*  * in Multiprecision on the secular equation. */
/*  * */
/*  * @param s the pointer to the mps_status struct. */
/*  * @param maxit Maximum number of iteration to perform. */
/*  * @return The number of approximated roots after the iteration. */
/*  *\/ */
/* int */
/* mps_secular_ga_miterate (mps_status * s, int maxit, mps_boolean just_regenerated) */
/* { */
/*   int computed_roots = 0; */
/*   int iterations = 0; */
/*   int i, k; */
/*   int nit = 0; */
/*   int it_threshold; */

/*   mpc_t corr, abcorr; */
/*   cdpe_t ctmp; */
/*   rdpe_t modcorr, rtmp; */
/*   rdpe_t * drad = rdpe_valloc (s->n); */

/*   mps_secular_equation * sec = s->secular_equation; */

/* #ifndef DISABLE_DEBUG */
/*   clock_t *my_clock = mps_start_timer (); */
/* #endif */

/*   /\* The data used to determined if the radius has been */
/*    * set and to intercommunicate with the iterator *\/ */
/*   mps_secular_iteration_data user_data; */

/*   /\* Cluster data used in iterations *\/ */
/*   mps_cluster_item * c_item; */
/*   mps_cluster * cluster; */
/*   mps_root * root; */

/*   MPS_DEBUG_WITH_INFO (s, "Precision is at %ld bits", s->mpwp); */
/*   MPS_DEBUG_RDPE (s, s->mp_epsilon, "Machine epsilon is s->mp_epsilon"); */

/*   MPS_DEBUG_THIS_CALL; */

/*   /\* Init data with the right precision *\/ */
/*   mpc_init2 (corr, s->mpwp); */
/*   mpc_init2 (abcorr, s->mpwp); */

/*   sec->best_approx = false; */

/*   /\* Iterate with newton until we have good approximations */
/*    * of the roots *\/ */
/*   for (i = 0; i < s->n; i++) */
/*     { */
/*       /\* Set again to false if the root is already approximated *\/ */
/*       if (MPS_ROOT_STATUS_IS_COMPUTED (s, i)) */
/* 	{ */
/* 	  MPS_DEBUG_WITH_INFO (s, "Setting again[%d] to false since the root is ready for output (or isolated)", i); */
/* 	  s->again[i] = false; */
/* 	} */

/*       if (s->again[i]) */
/*         s->rootwp[i] = s->mpwp; */
/*       else */
/*         computed_roots++; */
/*     } */

/*   if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)  */
/*     { */
/*       MPS_DEBUG_WITH_INFO (s, "%d roots are already approximated on the start of miterate", computed_roots); */
/*     } */

/*   /\* Set the iteration threshold to two times the remaining roots */
/*    * to compute. *\/ */
/*   it_threshold = 2 * (s->n - computed_roots); */

/*   while (computed_roots < s->n && iterations < maxit) */
/*     { */
/*       /\* Increase iterations counter *\/ */
/*       iterations++; */

/*       for (c_item = s->clusterization->first; c_item != NULL; c_item = c_item->next) */
/*         { */
/* 	  cluster = c_item->cluster; */
/*           for (root = cluster->first; root != NULL; root = root->next) */
/*             { */
/*               k = root->k; */
/*               if (s->again[k]) */
/*                 { */
/*                   nit++; */

/* 		  /\* Set the correct index *\/ */
/* 		  user_data.k = k; */
/*                   mps_secular_mnewton (s, s->mroot[k], s->drad[k], corr, */
/*                                        &s->again[k], &user_data, false); */

/* 		  /\* Apply Aberth correction *\/ */
/* 		  mps_maberth_s (s, k, cluster, abcorr); */
/* 		  mpc_mul_eq (abcorr, corr); */
/* 		  mpc_ui_sub (abcorr, 1U, 0U, abcorr); */

/* 		  if (!mpc_eq_zero (abcorr)) */
/* 		    { */
/* 		      mpc_div (abcorr, corr, abcorr); */
/* 		      mpc_sub_eq (s->mroot[k], abcorr); */

/* 		      /\* Correct the radius *\/ */
/* 		      mpc_get_cdpe (ctmp, abcorr); */
/* 		      cdpe_mod (modcorr, ctmp); */
/* 		      rdpe_add_eq (s->drad[k], modcorr); */
/* 		    } */
/* 		  else */
/* 		    s->again[k] = true; */

/* 		  mpc_get_cdpe (ctmp, s->mroot[k]); */
/* 		  cdpe_mod (rtmp, ctmp); */
/* 		  rdpe_mul_eq_d (rtmp, 8.0); */
/* 		  rdpe_mul_eq (rtmp, s->mp_epsilon); */

/* 		  rdpe_add_eq (s->drad[k], rtmp); */

/*                   if (!s->again[k]) */
/*                     computed_roots++; */
/*                 } */
/*             } */
/*         } */
/*     } */

/*   /\* Deallocate multiprecision local variables *\/ */
/*   mpc_clear (abcorr); */
/*   mpc_clear (corr); */

/*   if (s->debug_level & MPS_DEBUG_APPROXIMATIONS) */
/*     { */
/*       MPS_DEBUG (s, "Performed %d iterations", nit); */
/*     } */

/*   if (s->debug_level & MPS_DEBUG_APPROXIMATIONS) */
/*     mps_dump (s); */

/*   if (nit <= it_threshold && just_regenerated) */
/*     s->secular_equation->best_approx = true; */

/*   /\* Perform cluster analysis *\/ */
/*   mps_mradii (s, drad); */
/*   mps_mcluster (s, drad, 2.0 * s->n); */
/*   mps_mmodify (s, false); */

/*   /\* These lines are used to debug the again vector, but are not useful */
/*    * at the moment being *\/ */
/*   if (s->debug_level & MPS_DEBUG_APPROXIMATIONS) */
/*     { */
/*       __MPS_DEBUG (s, "Again vector = "); */
/*       for (i = 0; i < s->n; i++) */
/* 	{ */
/* 	  fprintf (s->logstr, "%d ", s->again[i]); */
/* 	} */
/*       fprintf (s->logstr, "\n"); */
/*       MPS_DEBUG (s, "Status = "); */
/*       for (i = 0; i < s->n; i++) */
/* 	{ */
/* 	  fprintf (s->logstr, " %4d: %s ", i, */
/* 		   MPS_ROOT_STATUS_TO_STRING (s->root_status[i])); */
/* 	} */
/*       fprintf (s->logstr, "\n"); */
/*     } */
  
/*   /\* Clock the routine *\/ */
/* #ifndef DISABLE_DEBUG */
/*   s->mp_iteration_time += mps_stop_timer (my_clock); */
/* #endif */

/*   rdpe_vfree (drad); */

/*   /\* Return the number of approximated roots *\/ */
/*   return computed_roots; */
/* } */
