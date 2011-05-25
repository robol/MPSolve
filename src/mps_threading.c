/*
 * mps_threading.c
 *
 *  Created on: 19/mag/2011
 *      Author: leonardo
 */

#include <mps/core.h>
#include <mps/threading.h>
#include <pthread.h>

/**
 * @brief Create a new mps_thread_pool.
 */
mps_thread_pool*
mps_thread_pool_new(int n)
{
  int i;
  mps_thread_pool *p = (mps_thread_pool*) malloc(sizeof(mps_thread_pool));
  pthread_mutex_init(&p->full_mutex, NULL);
  pthread_mutex_init(&p->ready_mutex, NULL);
  pthread_mutex_init(&p->aberth_mutex, NULL);
  pthread_cond_init(&p->full, NULL);
  p->ready = mps_boolean_valloc(n);
  p->threads = (pthread_t*) malloc(sizeof(pthread_t) * n);
  p->n = n;
  for (i = 0; i < n; i++)
    p->ready[i] = true;
  return p;
}

/**
 * @brief Free a mps_thread_pool previosuly allocated with
 * mps_thread_pool_new ()
 */
void
mps_thread_pool_free(mps_thread_pool* pool)
{
  mps_boolean_vfree(pool->ready);
  free(pool->threads);
  free(pool);
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

  for (iter = 0; iter < s->max_it; iter++)
    {
      for (i = data->thread; i < s->n; i += data->n_threads)
        {
          if (s->again[i])
            {
              (*data->it)++;

              rad1 = s->frad[i];
              if (s->data_type[0] != 'u')
                {
                  mps_fnewton(s, s->n, s->froot[i], &s->frad[i], corr, s->fpc,
                      s->fap, &s->again[i]);
                  if (iter == 0 && !s->again[i] && s->frad[i] > rad1 && rad1
                      != 0)
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
                  (*s->fnewton_usr)(s, s->froot[i], &s->frad[i], corr,
                      &s->again[i]);
                }
              else
                {
                  mps_fnewton_usr(s, s->froot[i], &s->frad[i], corr,
                      &s->again[i]);
                }

              if (s->again[i] ||
              /* the correction is performed only if iter!=1 or rad(i)!=rad1 */
              s->data_type[0] == 'u' || iter != 0 || s->frad[i] != rad1)
                {
                  pthread_mutex_lock(data->aberth_mutex);
                  mps_faberth(s, i, abcorr);
                  pthread_mutex_unlock(data->aberth_mutex);
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
                      return 0;
                    }
                }

            }
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
  pthread_mutex_t aberth_mutex = PTHREAD_MUTEX_INITIALIZER;

  *it = 0;
  *excep = false;

  /* count the number of approximations in the root neighbourhood */
  for (i = 0; i < s->n; i++)
    if (!s->again[i])
      nzeros++;
  if (nzeros == s->n)
    {
      free (threads);
      return;
    }

  *it = 0;
  *excep = false;

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
      pthread_create(&threads[i], NULL, &mps_thread_fpolzer_worker, data + i);
    }

  for (i = 0; i < n_threads; i++)
    {
      pthread_join(threads[i], NULL);
    }
  free (data);
  free(threads);

  *excep = true;
}

/**
 * @brief Worker for the mpolzer routine.
 */
void*
mps_thread_mpolzer_worker(void* data_ptr)
{
  mps_thread_worker_data *data = (mps_thread_worker_data*) data_ptr;
  mps_status* s = data->s;
  int i, j, iter, l;
  tmpc_t corr, abcorr;
  rdpe_t eps, rad1, rtmp;
  cdpe_t ctmp;

  tmpc_init2(abcorr, s->mpwp);
  tmpc_init2(corr, s->mpwp);

  rdpe_mul_d(eps, s->mp_epsilon, (double) 4 * s->n);

  /* initialize the iteration counter */
  (*data->excep) = false;

  /* Start Aberth's iterations */
  for (iter = 0; iter < s->max_it; iter++)
    { /* do_iter: */
      for (j = 0; j < s->nclust; j++)
        {
          for (i = data->thread; i < s->punt[j + 1] - s->punt[j]; i
              += data->n_threads)
            { /* do_indice: */
              l = s->clust[s->punt[j] + i];
              if (s->again[l])
                {
                  (*data->it)++;
                  if (s->data_type[0] != 'u')
                    {
                      /* sparse/dense polynomial */
                      rdpe_set(rad1, s->drad[l]);
                      mps_mnewton(s, s->n, s->mroot[l], s->drad[l], corr,
                          s->mfpc, s->mfppc, s->dap, s->spar, &s->again[l]);
                      if (iter == 0 && !s->again[l]
                          && rdpe_gt(s->drad[l], rad1) && rdpe_ne(rad1,
                          rdpe_zero))
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
                      (*s->mnewton_usr)(s, s->mroot[l], s->drad[l], corr,
                          &s->again[l]);
                    }
                  else
                    {
                      mps_mnewton_usr(s, s->mroot[l], s->drad[l], corr,
                          &s->again[l]);
                    }

                  if (s->again[l] ||
                  /* the correction is performed only if iter!=1 or rad[l]!=rad1 */
                  s->data_type[0] == 'u' || iter != 0 || rdpe_ne(s->drad[l],
                      rad1))
                    {
                      pthread_mutex_lock(data->aberth_mutex);
                      mps_maberth_s(s, l, j, abcorr);
                      pthread_mutex_unlock(data->aberth_mutex);
                      mpc_mul_eq(abcorr, corr);
                      mpc_neg_eq(abcorr);
                      mpc_add_eq_ui(abcorr, 1, 0);
                      mpc_div(abcorr, corr, abcorr);
                      mpc_sub_eq(s->mroot[l], abcorr);
                      mpc_get_cdpe(ctmp, abcorr);
                      cdpe_mod(rtmp, ctmp);
                      rdpe_add_eq(s->drad[l], rtmp);
                    }

                  /* check for new approximated roots */
                  if (!s->again[l])
                    {
                      (*data->nzeros)++;
                      if ((*data->nzeros) == s->n)
                        goto endfun;
                    }

                }
            }
        }
      if ((*data->nzeros) == s->n)
        goto endfun;
    }
  (*data->excep) = true;

  endfun: /* free local MP variables */
          tmpc_clear(corr);
          tmpc_clear(abcorr);

  return 0;
}

/**
 * @brief Drop-in replacement for the stock mpolzer.
 */
void
mps_thread_mpolzer(mps_status* s, int *it, mps_boolean *excep)
{
  int i, nzeros = 0, n_threads = s->n_threads;
  pthread_t* threads = (pthread_t*) malloc(sizeof(pthread_t) * n_threads);
  pthread_mutex_t aberth_mutex  = PTHREAD_MUTEX_INITIALIZER;
  mps_thread_worker_data* data;

  *it = 0;
  *excep = false;

  /* count the number of approximations in the root neighbourhood */
  for (i = 0; i < s->n; i++)
    if (!s->again[i])
      nzeros++;
  if (nzeros == s->n)
    {
      free (threads);
      return;
    }

  *it = 0;
  *excep = false;

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
      pthread_create(&threads[i], NULL, &mps_thread_mpolzer_worker, data + i);
    }

  for (i = 0; i < n_threads; i++)
    {
      pthread_join(threads[i], NULL);
    }
  free(data);
  free(threads);
}

