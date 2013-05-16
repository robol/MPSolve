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


#include <float.h>
#include <mps/mps.h>
#include <pthread.h>
#include <stdio.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef __WINDOWS
#include <windows.h>
#endif

/**
 * @brief Get number of logic cores on the local machine, or
 * 0 if that information is not available with the method
 * known to this implementations.
 */
int
mps_thread_get_core_number (mps_context * s)
{
  int cores = 0;
  char * cores_env = NULL;

#ifdef __WINDOWS
  SYSTEM_INFO windows_sys_info;
#endif

  if ((cores_env = getenv ("MPS_JOBS")) != NULL)
    {
      /* Give reasonable bounds to the possible values of MPS_JOBS */
      cores = MAX (1, MIN (MPS_MAX_CORES, atoi (cores_env)));

      return cores;
    }

  /* Test for POSIX platforms */
#ifdef HAVE_SYSCONF
  cores = sysconf (_SC_NPROCESSORS_ONLN);
#endif

#ifdef __WINDOWS
  GetSystemInfo (&windows_sys_info);
  cores = windows_sys_info.dwNumberOfProcessors;
#endif

  if (cores != 0)
    MPS_DEBUG_WITH_INFO (s, "Found %d cores on this system", cores);

  /* In case no runtime method of finding the available cores
   * worked out, select a fixed value. */
  if (cores <= 0)
  {
    cores = 8;
    if (s->debug_level & MPS_DEBUG_INFO)
    {
      MPS_DEBUG(s, "No runtime information about available cores found");
      MPS_DEBUG(s, "Selecting a fixed number of %d threads", cores);
      MPS_DEBUG(s, "Use the MPS_JOBS environment variable to override this value");
    }
  }

  return cores;
}


/**
 * @brief Create a new mps_thread_job_queue that can
 * handle at most max_iter iterations for n_roots roots.
 */
mps_thread_job_queue *
mps_thread_job_queue_new (mps_context * s)
{
  /* Space allocation and related jobs */
  mps_thread_job_queue *q;
  q = (mps_thread_job_queue *) mps_malloc (sizeof (mps_thread_job_queue));

  pthread_mutex_init (&q->mutex, NULL);

  /* Set initial data */
  q->iter = 0;
  q->n_roots = s->n;
  q->max_iter = s->max_it;
  q->cluster_item = s->clusterization->first;
  q->root = q->cluster_item->cluster->first;
  return q;
}

/*
 * @brief Free a mps_thread_job_queue previously allocated
 * with mps_thread_job_queue_new ().
 */
void
mps_thread_job_queue_free (mps_thread_job_queue * q)
{
  pthread_mutex_destroy (&q->mutex);
  free (q);
}

/**
 * @brief Obtain iter and i for the next available job.
 */
mps_thread_job
mps_thread_job_queue_next (mps_context * s, mps_thread_job_queue * q)
{
  mps_thread_job j;
  pthread_mutex_lock (&q->mutex);

  j.i = 0;
  j.cluster_item = NULL;

  if (q->iter == MPS_THREAD_JOB_EXCEP)
    {
      j.iter = MPS_THREAD_JOB_EXCEP;
    }
  else
    {
      /* Assigning the root */
      j.i = q->root->k;
      j.cluster_item = q->cluster_item;
      j.iter = q->iter;

      /* Get the next element of the cluster, incrementing the queue */
      q->root = q->root->next ;

      /* Check if the previous one was the last element in the
       * cluster, and if that's the case pass to the next one. */
      if (q->root == NULL)
        {
          q->cluster_item = q->cluster_item->next;

          /* If we got to the end of the clusterization restart from
           * the first cluster and dump the iteration counter. */
          if (q->cluster_item == NULL)
            {
              q->cluster_item = s->clusterization->first;
              q->iter++;
            }

          q->root = q->cluster_item->cluster->first;
      

          /* Check if maximum number of iteration was reached and
           * if that was the case set j->iter to MPS_THREAD_JOB_EXCEP.  */
          if (j.iter == q->max_iter)
            {
              j.iter = MPS_THREAD_JOB_EXCEP;
              q->iter = MPS_THREAD_JOB_EXCEP;
            }
        }
    }

  pthread_mutex_unlock (&q->mutex);
  return j;
}

void *
mps_thread_mainloop (void * thread_ptr)
{
  mps_thread * thread = (mps_thread *) thread_ptr;
  mps_thread_pool * pool = thread->pool;

  while (thread->alive)
    {
      /* Try to pop a work item from the queue, if available. */
      pthread_mutex_lock (&pool->work_completed_mutex);
      pthread_mutex_lock (&pool->queue_changed_mutex);

      if (pool->queue->first != NULL)
        {
          mps_thread_pool_queue_item * item = pool->queue->first;

          if (!thread->busy)
            {
              pool->busy_counter++;
              thread->busy = true;
            }

          /* Pop item from the queue and release the lock on it */
          pool->queue->first = item->next;
          if (item->next == NULL)
            pool->queue->last = item;

          pthread_mutex_unlock (&pool->queue_changed_mutex);
          pthread_mutex_unlock (&pool->work_completed_mutex);

          item->work (item->args);
          free (item);
        }
      else
      {
        /* Check if other threads are sleeping */
        if (thread->busy)
          {
            pool->busy_counter--;
            thread->busy = false;
          }
        pthread_cond_signal (&pool->work_completed_cond);
        pthread_mutex_unlock (&pool->work_completed_mutex);

        if (!thread->alive)
        {
          pthread_mutex_unlock (&pool->queue_changed_mutex);
          pthread_exit (NULL);
        }

        pthread_cond_wait (&pool->queue_changed, &pool->queue_changed_mutex);
        pthread_mutex_unlock (&pool->queue_changed_mutex);
      }
    }

  pthread_exit (NULL);
  return NULL;
}

/**
 * @brief Start the thread mainloop.
 */
void
mps_thread_start_mainloop (mps_context * s, mps_thread * thread)
{
  pthread_create (thread->thread, NULL, &mps_thread_mainloop, thread);
}

/**
 * @brief Limit the maximum number of threads that can be used in the thread pool.
 */
void mps_thread_pool_set_concurrency_limit (mps_context * s, mps_thread_pool * pool, 
                                            unsigned int concurrency_limit)
{
  if (pool == NULL)
    pool = s->pool;

  if (concurrency_limit == 0)
    concurrency_limit = mps_thread_get_core_number (s);

  if (concurrency_limit < pool->concurrency_limit)
  {
    mps_thread * old_first = pool->first;
    mps_thread * thread;
    int i = 0;

    for (thread = pool->first; i < (pool->concurrency_limit - concurrency_limit); thread = thread->next, i++);

    pool->first = thread;
    pool->n = concurrency_limit;

    i = 0;
    for (thread = old_first; i < (pool->concurrency_limit - concurrency_limit); i++)
    {
      mps_thread * next = thread->next;
      mps_thread_free (s, thread);
      thread = next;
    }
  }
  else
  {
    int i = 0;
    for (i = 0; i < concurrency_limit - pool->concurrency_limit; i++)
      mps_thread_pool_insert_new_thread (s, s->pool);
  }

  pool->concurrency_limit = concurrency_limit;
}

void
mps_thread_pool_assign (mps_context * s, mps_thread_pool * pool, 
                        mps_thread_work work, void * args)
{
  if (!pool)
    pool = s->pool;

  /* Insert the job in the queue */
  pthread_mutex_lock (&pool->queue_changed_mutex);

  mps_thread_pool_queue_item * item = mps_new (mps_thread_pool_queue_item);

  item->work = work;
  item->args = args;

  if (pool->queue->first == NULL)
    {
      pool->queue->first = pool->queue->last = item;
      item->next = NULL;
    }
  else
    {
      pool->queue->last->next = item;
      pool->queue->last = item;
      item->next = NULL;
    }

  pthread_cond_signal (&pool->queue_changed);
  pthread_mutex_unlock (&pool->queue_changed_mutex);
}

/**
 * @brief Wait for a thread pool to complete its jobs.
 */
void
mps_thread_pool_wait (mps_context * s, mps_thread_pool * pool)
{
  pthread_mutex_lock (&pool->work_completed_mutex);

  while (true)
    {
      if (pool->busy_counter == 0 && pool->queue->first == NULL)
        {
          pthread_mutex_unlock (&pool->work_completed_mutex);
          return;
        }
      else
        {
          pthread_cond_wait (&pool->work_completed_cond, &pool->work_completed_mutex);
        }
    }
}

/**
 * @brief Allocate a new <code>mps_thread</code> and start its mainloop.
 */
mps_thread * 
mps_thread_new (mps_context * s, mps_thread_pool * pool)
{  
  if (!pool)
    pool = s->pool;
  mps_thread * thread = mps_new (mps_thread);

  /* Set the initial values in the thread */
  thread->data = NULL;
  pthread_mutex_init (&thread->busy_mutex, NULL); 
  pthread_cond_init  (&thread->start_condition, NULL); 
  thread->thread = mps_new (pthread_t);
  thread->work = NULL;
  thread->args = NULL;
  thread->alive = true;
  thread->pool = pool;
  thread->busy = false;

  /* Start the thread mainloop */
  mps_thread_start_mainloop (s, thread);

  return thread;
}

/**
 * @brief Free a thread asking it to stop.
 */
void
mps_thread_free (mps_context * s, mps_thread * thread)
{
  /* Wait for the thread to finish its work, if it is doing something */
  /* pthread_mutex_lock (&thread->busy_mutex); */
  /* pthread_mutex_unlock (&thread->busy_mutex); */
  pthread_mutex_lock (&thread->pool->queue_changed_mutex);
  thread->alive = false;

  /* Start the thread, if it is not running */
  pthread_cond_broadcast (&thread->pool->queue_changed);
  pthread_mutex_unlock (&thread->pool->queue_changed_mutex);

  pthread_join (*thread->thread, NULL);

  free (thread->thread);
  free (thread);
}

/**
 * @brief Create a new thread and add it to the specified thread pool.
 */
void
mps_thread_pool_insert_new_thread (mps_context * s, mps_thread_pool * pool)
{
  if (!pool)
    pool = s->pool;

  mps_thread * thread = mps_thread_new (s, pool);
  
  thread->next = pool->first; 
  pool->first = thread;
  pool->n++;
}

/**
 * @brief Allocate a new thread pool and return a pointer to it, 
 * with a number of threads suitable for this system.
 */
mps_thread_pool *
mps_thread_pool_new (mps_context * s, int n_threads)
{
  mps_thread_pool * pool = mps_new (mps_thread_pool);
  int threads = mps_thread_get_core_number (s); 
  int i;

  if (n_threads != 0)
    threads = n_threads;

  pool->n = 0;
  pool->first = NULL;

  pool->queue = mps_new (mps_thread_pool_queue);
  pool->queue->first = pool->queue->last = NULL;

  pthread_mutex_init (&pool->queue_changed_mutex, NULL);
  pthread_cond_init (&pool->queue_changed, NULL);

  pthread_mutex_init (&pool->work_completed_mutex, NULL);
  pthread_cond_init (&pool->work_completed_cond, NULL);

  pool->busy_counter = 0;
  
  for (i = 0; i < threads; i++) 
    mps_thread_pool_insert_new_thread (s, pool); 

  pool->concurrency_limit = threads;

  mps_thread_pool_wait (s, pool);

  return pool;
}

/**
 * @brief Free a thread pool and all its threads, waiting for them
 * to terminate.
 */
void 
mps_thread_pool_free (mps_context * s, mps_thread_pool * pool)
{
  if (!pool)
    pool = s->pool;

  mps_thread * thread = pool->first;
  mps_thread * next_thread;

  mps_thread_pool_wait (s, pool);

  while (thread)
    {
      next_thread = thread->next;
      mps_thread_free (s, thread);
      thread = next_thread;
    }

  free (pool->queue);
  free (pool);
}

int mps_thread_get_id (mps_context * s, mps_thread_pool * pool)
{
  pthread_t self = pthread_self ();
  int i = 0;

  mps_thread * thread = pool->first;
  while (thread)
    {
      if (pthread_equal (*thread->thread, self))
        {
          return i;
        }
      i++;
      thread = thread->next;
    }
  return -1;
}
