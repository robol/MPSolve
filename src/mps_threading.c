/*
 * mps_threading.c
 *
 *  Created on: 19/mag/2011
 *      Author: leonardo
 */

/**
 * @file
 * @brief Implementation of coordination routines for multithreading.
 */

#include <mps/core.h>
#include <mps/threading.h>
#include <pthread.h>
#include <semaphore.h>
#include <sched.h>

/**
 * @brief Create a new mps_thread_job_queue that can
 * handle at most max_iter iterations for n_roots roots.
 */
mps_thread_job_queue*
mps_thread_job_queue_new(mps_status* s)
{
  /* Space allocation and related jobs */
  mps_thread_job_queue* q;
  q = (mps_thread_job_queue*) malloc(sizeof(mps_thread_job_queue));
  pthread_mutex_init(&q->mutex, NULL);

  /* Set initial data */
  q->i = 0;
  q->iter = 0;
  q->n_roots = s->n;
  q->max_iter = s->max_it;
  q->i_clust = 0;
  return q;
}

/*
 * @brief Free a mps_thread_job_queue previously allocated
 * with mps_thread_job_queue_new ().
 */
void
mps_thread_job_queue_free(mps_thread_job_queue* q)
{
  pthread_mutex_destroy(&q->mutex);
  free(q);
}

/**
 * @brief Obtain iter and i for the next available job.
 */
mps_thread_job
mps_thread_job_queue_next(mps_status* s, mps_thread_job_queue* q)
{
  mps_thread_job j;
  pthread_mutex_lock(&q->mutex);

  if (q->iter == MPS_THREAD_JOB_EXCEP)
    {
      j.iter = MPS_THREAD_JOB_EXCEP;
      pthread_mutex_unlock(&q->mutex);
      return j;
    }

  /* Get the next element of the cluster */
  j.i = (q->i)++;

  /* Check if the previous one was the last element in the
   * cluster, and if that's the case pass to the next one. */
  if (j.i == s->n)
    {
      j.iter = ++q->iter;
      q->i = 1;
      j.i = 0;

      /* Check if maximum number of iteration was reached and
       * if that was the case set j->iter to MPS_THREAD_JOB_EXCEP.  */
      if (j.iter == q->max_iter)
        {
          j.iter = MPS_THREAD_JOB_EXCEP;
          q->iter = MPS_THREAD_JOB_EXCEP;
          pthread_mutex_unlock(&q->mutex);
          return j;
        }
    }
  else
    j.iter = q->iter;

  /* Find out in which cluster we are, if we get over
   * the limit of the actual cluster, step one over. */
  if (s->punt[q->i_clust] == j.i)
    q->i_clust++;

  /* If this is the first root we are not stepping to
   * the next cluster, but to the first one */
  if (j.i == 0)
    q->i_clust = 0;

  j.i_clust = q->i_clust;

  pthread_mutex_unlock(&q->mutex);
  return j;
}

/**
 * @brief Worker for the fpolzer routine.
 */
void*
mps_thread_fpolzer_worker(void* data_ptr)
{
  mps_thread_worker_data* data = (mps_thread_worker_data*) data_ptr;
  mps_status* s = data->s;
  int i, iter;
  cplx_t corr, abcorr;
  double rad1, modcorr;
  mps_thread_job job;

  while (!(*data->excep))
    {
      job = mps_thread_job_queue_next(s, data->queue);
      i = job.i;
      iter = job.iter;

      /* Check if we got over the maximum number of iterations */
      if (job.iter == MPS_THREAD_JOB_EXCEP)
        {
          (*data->excep) = true;
          return 0;
        }


      if (s->again[i])
        {
          /* Lock this roots to make sure that we are the only one working on it */
          pthread_mutex_lock(&data->roots_mutex[i]);

          /* Check if, while we were waiting, excep condition has been reached */
          if (*data->excep)
            {
              pthread_mutex_unlock(&data->roots_mutex[i]);
              return 0;
            }

          (*data->it)++;

          rad1 = s->frad[i];
          if (s->data_type[0] != 'u')
            {
              mps_fnewton(s, s->n, s->froot[i], &s->frad[i], corr, s->fpc,
                  s->fap, &s->again[i]);
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
              (*s->fnewton_usr)(s, s->froot[i], &s->frad[i], corr, &s->again[i]);
            }
          else
            {
              mps_fnewton_usr(s, s->froot[i], &s->frad[i], corr, &s->again[i]);
            }

          if (s->again[i] ||
          /* the correction is performed only if iter!=1 or rad(i)!=rad1 */
          s->data_type[0] == 'u' || iter != 0 || s->frad[i] != rad1)
            {
              mps_faberth(s, i, abcorr);
              cplx_mul_eq(abcorr, corr);
              cplx_sub(abcorr, cplx_one, abcorr);
              cplx_div(abcorr, corr, abcorr);
              cplx_sub_eq(s->froot[i], abcorr);
              modcorr = cplx_mod(abcorr);
              s->frad[i] += modcorr;
            }

          /* check for new approximated roots */
          if (!s->again[i])
            {
              (*data->nzeros)++;
              if (*data->nzeros == s->n)
                {
                  pthread_mutex_unlock(&data->roots_mutex[i]);
                  return 0;
                }
            }

          pthread_mutex_unlock(&data->roots_mutex[i]);
        }

    }

  return 0;
}

/**
 * @brief Drop-in replacement for the stock fpolzer routine.
 * This version adds multithread support.
 */
void
mps_thread_fpolzer(mps_status* s, int* it, mps_boolean* excep)
{
  int i, nzeros = 0, n_threads = s->n_threads;
  pthread_t* threads = (pthread_t*) malloc(sizeof(pthread_t) * n_threads);
  mps_thread_worker_data* data;
  pthread_mutex_t aberth_mutex =
  PTHREAD_MUTEX_INITIALIZER;
  pthread_mutex_t* roots_mutex = (pthread_mutex_t*) malloc(sizeof(pthread_mutex_t) * s->n);

  for(i = 0; i < s->n; i++)
    pthread_mutex_init(roots_mutex + i, NULL);

  /* Create a new job queue */
  mps_thread_job_queue* queue = mps_thread_job_queue_new(s);

  *it = 0;
  *excep = false;

  /* count the number of approximations in the root neighbourhood */
  for (i = 0; i < s->n; i++)
    if (!s->again[i])
      nzeros++;
  if (nzeros == s->n)
    {
      free(threads);
      return;
    }

  data = (mps_thread_worker_data*) malloc(sizeof(mps_thread_worker_data)
      * n_threads);

  for (i = 0; i < n_threads; i++)
    {
      data[i].it = it;
      data[i].nzeros = &nzeros;
      data[i].s = s;
      data[i].excep = excep;
      data[i].thread = i;
      data[i].n_threads = n_threads;
      data[i].aberth_mutex = &aberth_mutex;
      data[i].roots_mutex = roots_mutex;
      data[i].queue = queue;
      pthread_create(&threads[i], NULL, &mps_thread_fpolzer_worker, data + i);
    }

  for (i = 0; i < n_threads; i++)
    {
      pthread_join(threads[i], NULL);
    }
  free(data);
  free(threads);
  free(roots_mutex);
  mps_thread_job_queue_free(queue);
}

/**
 * @brief Worker for the mpolzer routine.
 */
void*
mps_thread_mpolzer_worker(void* data_ptr)
{
  mps_thread_worker_data *data = (mps_thread_worker_data*) data_ptr;
  mps_status* s = data->s;
  mps_thread_job job;
  int iter, l, k;
  tmpc_t corr, abcorr, mroot, diff;
  rdpe_t eps, rad1, rtmp;
  cdpe_t ctmp, abcorr_cdpe;

  tmpc_init2(abcorr, s->mpwp);
  tmpc_init2(corr, s->mpwp);
  tmpc_init2(mroot, s->mpwp);
  tmpc_init2(diff, s->mpwp);

  rdpe_mul_d(eps, s->mp_epsilon, (double) 4 * s->n);

  /* initialize the iteration counter */
  (*data->excep) = false;

  /* Continue to iterate while exception condition has not
   * been reached and there more roots to approximate   */
  while (!(*data->excep) && (*data->nzeros) < s->n)
    {
      /* Get next job for this thread */
      job = mps_thread_job_queue_next(s, data->queue);

      /* Set variables to be used in the rest of the code */
      iter = job.iter;

      /* Check if we exceeded the maximum number of iterations */
      if (job.iter == MPS_THREAD_JOB_EXCEP)
        {
          (*data->excep) = true;
          return 0;
        }

      l = s->clust[job.i];
      if (s->again[l])
        {
          /* Lock roots_mutex to assure that we are the only thread
           * working on this root. Parallel computation on the same
           * root is not useful, since we would be performing the
           * same computations.                                  */
          pthread_mutex_lock(&data->roots_mutex[l]);

          /* Check if, while we were waiting, excep condition has been reached,
           * or all the zeros has been approximated.                         */
          if (*data->excep || (*data->nzeros) >= s->n)
            {
              pthread_mutex_unlock(&data->roots_mutex[l]);
              return 0;
            }

          /* Increment total iteration counter */
          (*data->it)++;

          /* Copy locally the root to work on */
          pthread_mutex_lock(&data->aberth_mutex[l]);
          mpc_set(mroot, s->mroot[l]);
          pthread_mutex_unlock(&data->aberth_mutex[l]);

          if (s->data_type[0] != 'u')
            {
              /* sparse/dense polynomial */
              rdpe_set(rad1, s->drad[l]);
              mps_mnewton(s, s->n, mroot, s->drad[l], corr, s->mfpc,
                  s->mfppc, s->dap, s->spar, &s->again[l], data->thread);
              if (iter == 0 && !s->again[l] && rdpe_gt(s->drad[l], rad1)
                  && rdpe_ne(rad1, rdpe_zero))
                rdpe_set(s->drad[l], rad1);

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
          else /* user's polynomial */
          if (s->mnewton_usr != NULL)
            {
              (*s->mnewton_usr)(s, mroot, s->drad[l], corr, &s->again[l]);
            }
          else
            {
              mps_mnewton_usr(s, mroot, s->drad[l], corr, &s->again[l]);
            }

          if (s->again[l] ||
          /* the correction is performed only if iter!=1 or rad[l]!=rad1 */
          s->data_type[0] == 'u' || iter != 0 || rdpe_ne(s->drad[l], rad1))
            {
//              mps_maberth_s(s, l, job.i_clust, abcorr);
              /* Compute Aberth correction manually so we can lock */
              cdpe_set(abcorr_cdpe, cdpe_zero);
              for(k = s->punt[job.i_clust]; k < s->punt[job.i_clust + 1]; k++)
                {
                  if (k == l)
                    continue;
                  pthread_mutex_lock(&data->aberth_mutex[k]);
                  mpc_sub(diff, mroot, s->mroot[k]);
                  pthread_mutex_unlock(&data->aberth_mutex[k]);
                  mpc_get_cdpe(ctmp, diff);
                  cdpe_inv_eq(ctmp);
                  cdpe_add_eq(abcorr_cdpe, ctmp);
                }
              mpc_set_cdpe(abcorr, abcorr_cdpe);

              mpc_mul_eq(abcorr, corr);
              mpc_neg_eq(abcorr);
              mpc_add_eq_ui(abcorr, 1, 0);
              mpc_div(abcorr, corr, abcorr);
              mpc_sub_eq(mroot, abcorr);
              mpc_get_cdpe(ctmp, abcorr);
              cdpe_mod(rtmp, ctmp);
              rdpe_add_eq(s->drad[l], rtmp);

              pthread_mutex_lock(&data->aberth_mutex[l]);
              mpc_set(s->mroot[l], mroot);
              pthread_mutex_unlock(&data->aberth_mutex[l]);
            }

          /* check for new approximated roots */
          if (!s->again[l])
            {
              (*data->nzeros)++;
              if ((*data->nzeros) == s->n)
                {
                  pthread_mutex_unlock(&data->roots_mutex[l]);
                  goto endfun;
                }
            }

          pthread_mutex_unlock(&data->roots_mutex[l]);
        }

      if ((*data->nzeros) == s->n)
        {
          goto endfun;
        }
    }

  endfun: /* free local MP variables */
  tmpc_clear(corr);
  tmpc_clear(abcorr);
  tmpc_clear(mroot);
  tmpc_clear(diff);

  return 0;
}

/**
 * @brief Drop-in threaded replacement for the stock mpolzer.
 */
void
mps_thread_mpolzer(mps_status* s, int *it, mps_boolean *excep)
{
  int i, nzeros = 0, n_threads = s->n_threads;
  pthread_t* threads = (pthread_t*) malloc(sizeof(pthread_t) * n_threads);
  mps_thread_worker_data* data;

  /* Allocate and the init mutexes needed by the routine */
  pthread_mutex_t* roots_mutex = (pthread_mutex_t*) malloc(sizeof(pthread_mutex_t) * s->n);
  pthread_mutex_t* aberth_mutex = (pthread_mutex_t*) malloc(sizeof(pthread_mutex_t) * s->n);
  for(i = 0; i < s->n; i++)
    {
      pthread_mutex_init(&aberth_mutex[i], NULL);
      pthread_mutex_init(&roots_mutex[i], NULL);
    }

  /* Create a new work queue */
  mps_thread_job_queue* queue = mps_thread_job_queue_new(s);

  *it = 0;
  *excep = false;

  /* Count the number of approximations in the root neighbourhood */
  for (i = 0; i < s->n; i++)
    if (!s->again[i])
      nzeros++;
  if (nzeros == s->n)
    {
      free(threads);
      return;
    }

  data = (mps_thread_worker_data*) malloc(sizeof(mps_thread_worker_data)
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
      data[i].queue = queue;
      data[i].roots_mutex = roots_mutex;
      pthread_create(&threads[i], NULL, &mps_thread_mpolzer_worker, data + i);
    }

  /* Wait for the threads to complete */
  for (i = 0; i < n_threads; i++)
    {
      pthread_join(threads[i], NULL);
    }

  /* Free data and exit */
  free(data);
  free(threads);
  for(i = 0; i < s->n; i++)
    {
      pthread_mutex_destroy(&roots_mutex[i]);
      pthread_mutex_destroy(&aberth_mutex[i]);
    }
  free(roots_mutex);
  free(aberth_mutex);
  mps_thread_job_queue_free(queue);
}

