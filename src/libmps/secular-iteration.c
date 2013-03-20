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
#include <math.h>
#include <string.h>

/* Routine called to perform the floating point newton iterations. */
void *
__mps_secular_ga_fiterate_worker (void* data_ptr)
{
  mps_thread_worker_data *data = (mps_thread_worker_data *) data_ptr;
  mps_context *s = data->s;
  int i;
  cplx_t corr, abcorr;
  double modcorr;
  mps_thread_job job;

  while (true)
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

      if (s->root[i]->again && !s->root[i]->approximated)
        {
          /* Increment the number of performed iterations */
#if defined(__GCC__)
          __sync_add_and_fetch (data->it, 1);
#else
          pthread_mutex_lock (data->gs_mutex);
          (*data->it)++;
          pthread_mutex_unlock (data->gs_mutex);
#endif
          cdpe_set_x (s->root[i]->dvalue, s->root[i]->fvalue);

          mps_secular_fnewton (s, MPS_POLYNOMIAL (s->secular_equation), s->root[i], corr);

          if (s->root[i]->status == MPS_ROOT_STATUS_NOT_FLOAT)
            {
              *data->excep = true;
              break;
            }

          /* Apply Aberth correction */
          mps_faberth_wl (s, i, abcorr, data->aberth_mutex);

          if (isnan (cplx_Re (abcorr)) || isnan (cplx_Im (abcorr))) 
            {
              s->root[i]->again = false;
              pthread_mutex_unlock (&data->roots_mutex[i]);
              continue;
            }

          cplx_mul_eq (abcorr, corr);
          cplx_sub (abcorr, cplx_one, abcorr);
          cplx_div (abcorr, corr, abcorr);

          if (cplx_check_fpe (abcorr))
            {
              s->root[i]->again = false;
              pthread_mutex_unlock (&data->roots_mutex[i]);
              continue;
            }

          if (!s->root[i]->again || s->root[i]->approximated)
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
          else 
            {
              pthread_mutex_lock (&data->aberth_mutex[i]);
              cplx_sub_eq (s->root[i]->fvalue, abcorr);
              pthread_mutex_unlock (&data->aberth_mutex[i]);

              /* Correct the radius */
              modcorr = cplx_mod (abcorr);
              s->root[i]->frad += modcorr;
            }

        }

      pthread_mutex_unlock (&data->roots_mutex[i]);
    }

 cleanup:

  return NULL;
}

/**
 * @brief Routine that performs a block of iteration
 * in floating point on the secular equation.
 *
 * @param s the pointer to the mps_context struct.
 * @param maxit Maximum number of iteration to perform.
 * @param just_regenerated true if this is the first iteration after a coefficient
 * regeneration. If just_regenerated is true and the iteration packet is completed
 * in less than 2 * (n - computed_roots) iterations that best_approx is set to true
 * in s->secular_equation so a raise in the precision will be triggered.
 * @return The number of approximated roots after the iteration.
 */
int
mps_secular_ga_fiterate (mps_context * s, int maxit, mps_boolean just_regenerated)
{
  int computed_roots = 0;
  int approximated_roots = 0;
  int root_neighborhood_roots = 0;
  int i;
  int nit = 0;
  int it_threshold = 0;
  mps_boolean excep = false;

#ifndef DISABLE_DEBUG
  clock_t *my_clock = mps_start_timer ();
#endif

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

  /* Mark the approximated roots as ready for output */
  for (i = 0; i < s->n; i++)
    {
      /* Set again to false if the root is already approximated. If a root is approximated but
       * it has less digits than the current precision don't stop the iterations on that component. */
      if (s->root[i]->status == MPS_ROOT_STATUS_ISOLATED ||
          s->root[i]->status == MPS_ROOT_STATUS_APPROXIMATED)
        {
          if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)
            {
              MPS_DEBUG_WITH_INFO (s, "Setting again[%d] to false since the root is ready for output (or isolated)", i);
            }
          s->root[i]->again = false;

          if (s->root[i]->status == MPS_ROOT_STATUS_APPROXIMATED)
            s->root[i]->approximated = true;
        }

      if (!s->root[i]->again || s->root[i]->approximated)
        computed_roots++;
    }

  MPS_DEBUG_WITH_INFO (s, "%d roots are already approximated at the start of the packet", computed_roots);

  it_threshold = s->n - computed_roots;

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

      mps_thread_pool_assign (s, s->pool,  
                              __mps_secular_ga_fiterate_worker, data + i);
    }

  mps_thread_pool_wait (s, s->pool);

  /* Check if the roots are improvable in floating point */
  MPS_DEBUG_WITH_INFO (s, "Performed %d iterations with floating point arithmetic",
                       nit);

  if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)
      mps_dump (s);

  /* Check if we need to get higher precision for the roots */
   s->best_approx = true; 
   for (i = 0; i < s->n; i++) 
     { 
       if (!s->root[i]->approximated) 
        s->best_approx = false; 
       if (s->root[i]->approximated) 
        approximated_roots++; 
       if (!s->root[i]->again) 
        root_neighborhood_roots++; 
     } 

  if (just_regenerated && (nit <= it_threshold))
    s->best_approx = true;

  MPS_DEBUG_WITH_INFO(s, "%d roots are approximated with the current precision", approximated_roots);
  MPS_DEBUG_WITH_INFO (s,"%d roots are in the root neighborhood", root_neighborhood_roots);
  MPS_DEBUG_WITH_INFO (s, "%d roots have reached a stop condition", computed_roots);

  if (excep)
    {
      MPS_DEBUG_WITH_INFO (s, "Switching to DPE arithmetic since there are roots not representable in standard floating point");
      for (i = 0; i < s->n; i++)
        {
          cdpe_set_x (s->root[i]->dvalue, s->root[i]->fvalue);
          rdpe_set_d (s->root[i]->drad, s->root[i]->frad);
          s->root[i]->status = MPS_ROOT_STATUS_CLUSTERED;
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
          fprintf (s->logstr, "%d ", s->root[i]->again);
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
  mps_context *s = data->s;
  int i;
  cdpe_t corr, abcorr, droot;
  rdpe_t modcorr;
  mps_thread_job job;

  while (true)
    {
      job = mps_thread_job_queue_next (s, data->queue);
      i = job.i;

      if (job.iter == MPS_THREAD_JOB_EXCEP)
        {
          return NULL;
        }

      pthread_mutex_lock (&data->roots_mutex[i]);

      if (s->root[i]->again && !s->root[i]->approximated)
        {
          /* Lock this roots to make sure that we are the only one working on it */
          cdpe_set (droot, s->root[i]->dvalue);

          (*data->it)++;

          mps_secular_dnewton (s, MPS_POLYNOMIAL (s->secular_equation), s->root[i], corr);

          /* Apply Aberth correction */
          mps_daberth_wl (s, i, abcorr, data->aberth_mutex);
          cdpe_mul_eq (abcorr, corr);
          cdpe_sub (abcorr, cdpe_one, abcorr);
          cdpe_div (abcorr, corr, abcorr);

          cdpe_sub_eq (droot, abcorr);

	       /* Correct the radius */
         if (s->root[i]->again)
           {
	           cdpe_mod (modcorr, abcorr);
	           rdpe_add_eq (s->root[i]->drad, modcorr);
           }

          if (!s->root[i]->again || s->root[i]->approximated)
            {
              if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)
                MPS_DEBUG (s, "Root %d again was set to false on iteration %d by thread %d", i, *data->it, data->thread);
              (*data->nzeros)++;
            }
          else
            cdpe_set (s->root[i]->dvalue, droot);
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
 * @param s the pointer to the mps_context struct.
 * @param maxit Maximum number of iteration to perform.
 * @return The number of approximated roots after the iteration.
 * @param just_regenerated true if this is the first iteration after a coefficient
 * regeneration. If just_regenerated is true and the iteration packet is completed
 * in less than 2 * (n - computed_roots) iterations that best_approx is set to true
 * in s->secular_equation so a raise in the precision will be triggered.
 */
int
mps_secular_ga_diterate (mps_context * s, int maxit, mps_boolean just_regenerated)
{
  int computed_roots = 0;
  int root_neighborhood_roots = 0;
  int approximated_roots = 0;
  int i;
  int nit = 0;
  int it_threshold = 0;

  s->operation = MPS_OPERATION_ABERTH_DPE_ITERATIONS;

#ifndef DISABLE_DEBUG
  clock_t *my_clock = mps_start_timer ();
#endif

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

  s->best_approx = false;

  /* Mark the approximated roots as ready for output */
  for (i = 0; i < s->n; i++)
    {
      /* Set again to false if the root is already approximated. If a root is approximated but
       * it has less digits than the current precision don't stop the iterations on that component. */
      if (s->root[i]->status == MPS_ROOT_STATUS_ISOLATED ||
          s->root[i]->status == MPS_ROOT_STATUS_APPROXIMATED)
        {
          if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)
            {
              MPS_DEBUG_WITH_INFO (s, "Setting again[%d] to false since the root is ready for output (or isolated)", i);
            }
          s->root[i]->again = false;
          s->root[i]->approximated = true;
        }

      if (!s->root[i]->again || s->root[i]->approximated)
        computed_roots++;
    }

  it_threshold = (s->n - computed_roots);

  MPS_DEBUG_WITH_INFO (s, "%d roots are already approximated at the start of the packet", computed_roots);

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

  /* Check if we need to get higher precision for the roots */
  s->best_approx = true;
  for (i = 0; i < s->n; i++)
    {
      if (!s->root[i]->approximated)
        s->best_approx = false;
      if (s->root[i]->approximated)
        approximated_roots++;
      if (!s->root[i]->again)
        root_neighborhood_roots++;
    }

  if (just_regenerated && (nit <= it_threshold))
    s->best_approx = true;

  MPS_DEBUG_WITH_INFO(s, "%d roots are approximated with the current precision", approximated_roots);
  MPS_DEBUG_WITH_INFO (s,"%d roots are in the root neighborhood", root_neighborhood_roots);
  MPS_DEBUG_WITH_INFO (s, "%d roots have reached a stop condition", computed_roots);

  /* These lines are used to debug the again vector, but are not useful
   * at the moment being */
  if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)
    {
      __MPS_DEBUG (s, "Again vector = ");
      for(i = 0; i < s->n; i++)
        {
          fprintf (s->logstr, "%d ", s->root[i]->again);
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
  mps_context *s = data->s;
  /* mps_secular_equation *sec = s->secular_equation; */
  int i;
  mpc_t corr, abcorr;
  mpc_t mroot; 
  rdpe_t modcorr;
  mps_thread_job job;

  mps_cluster * cluster = NULL;

  mpc_init2 (corr, s->mpwp);
  mpc_init2 (abcorr, s->mpwp);
  mpc_init2 (mroot, s->mpwp); 

  /* Get a copy of the MP coefficients that is local to this thread */
  while (true)
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

      if (s->root[i]->again && !s->root[i]->approximated)
        {
          /* Lock this roots to make sure that we are the only one working on it */
          pthread_mutex_lock (&data->aberth_mutex[i]); 
          mpc_set (mroot, s->root[i]->mvalue); 
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

          mps_secular_mnewton (s, MPS_POLYNOMIAL (s->secular_equation), s->root[i], corr);

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
            s->root[i]->again = true;


          if (!s->root[i]->again || s->root[i]->approximated)
            {
              if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)
                MPS_DEBUG (s, "Root %d again was set to false on iteration %d by thread %d", i, *data->it, data->thread);

              (*data->nzeros)++;
            }
          else 
            {
              pthread_mutex_lock (&data->aberth_mutex[i]); 
              mpc_set (s->root[i]->mvalue, mroot); 
              pthread_mutex_unlock (&data->aberth_mutex[i]); 
           
              /* Correct the radius */
              mpc_rmod (modcorr, abcorr);
              rdpe_add_eq (s->root[i]->drad, modcorr);

              mpc_rmod (modcorr, mroot);
              rdpe_mul_eq (modcorr, s->mp_epsilon);
              rdpe_add_eq (s->root[i]->drad, modcorr);
            }
        }

      pthread_mutex_unlock (&data->roots_mutex[i]);
    }

 cleanup:
  mpc_clear (mroot);
  mpc_clear (abcorr);
  mpc_clear (corr);

  return NULL;
}

/**
 * @brief Routine that performs a block of iteration
 * in floating point on the secular equation using
 * CDPE
 *
 * @param s the pointer to the mps_context struct.
 * @param maxit Maximum number of iteration to perform. 
 * @param just_regenerated true if this is the first iteration after a coefficient
 * regeneration. If just_regenerated is true and the iteration packet is completed
 * in less than 2 * (n - computed_roots) iterations that best_approx is set to true
 * in s->secular_equation so a raise in the precision will be triggered.
 * @return The number of approximated roots after the iteration.
 */
int
mps_secular_ga_miterate (mps_context * s, int maxit, mps_boolean just_regenerated)
{
  int computed_roots = 0;
  int approximated_roots = 0;
  int root_neighborhood_roots = 0;
  int i;
  int nit = 0;
  int it_threshold = 0;

  s->operation = MPS_OPERATION_ABERTH_MP_ITERATIONS;

#ifndef DISABLE_DEBUG
  clock_t *my_clock = mps_start_timer ();
#endif

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

  s->best_approx = false;

  it_threshold = s->n - computed_roots;

  /* Mark the approximated roots as ready for output */
  for (i = 0; i < s->n; i++)
    {
      /* Set again to false if the root is already approximated. If a root is approximated but
       * it has less digits than the current precision don't stop the iterations on that component. */
      if (s->root[i]->status == MPS_ROOT_STATUS_ISOLATED ||
          s->root[i]->status == MPS_ROOT_STATUS_APPROXIMATED)
        {
          if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)
            {
              MPS_DEBUG_WITH_INFO (s, "Setting again[%d] to false since the root is ready for output (or isolated)", i);
            }
          s->root[i]->again = false;
          s->root[i]->approximated = true;
        }

      if (!s->root[i]->again || s->root[i]->approximated)
        computed_roots++;
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

  /* Check if we need to get higher precision for the roots */
  s->best_approx = true;
  for (i = 0; i < s->n; i++)
    {
      if (!s->root[i]->approximated)
        s->best_approx = false;
      if (s->root[i]->approximated)
        approximated_roots++;
      if (!s->root[i]->again)
        root_neighborhood_roots++;
    }

  if (just_regenerated && (nit <= it_threshold))
    s->best_approx = true;

  MPS_DEBUG_WITH_INFO(s, "%d roots are approximated with the current precision", approximated_roots);
  MPS_DEBUG_WITH_INFO (s,"%d roots are in the root neighborhood", root_neighborhood_roots);
  MPS_DEBUG_WITH_INFO (s, "%d roots have reached a stop condition", computed_roots);

  /* These lines are used to debug the again vector, but are not useful
   * at the moment being */
  if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)
    {
      __MPS_DEBUG (s, "Again vector = ");
      for(i = 0; i < s->n; i++)
        {
          fprintf (s->logstr, "%d ", s->root[i]->again);
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

