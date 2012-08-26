/*
 * mps_threading.c
 *
 *  Created on: 19/mag/2011
 *      Author: leonardo
 */

#include <float.h>
#include <mps/mps.h>
#include <pthread.h>
#include <stdio.h>

/**
 * @brief Get number of logic cores on the local machine, or
 * 0 if that information is not available with the method
 * known to this implementations.
 */
int
mps_thread_get_core_number (mps_status * s)
{
  FILE *cpuinfo = fopen ("/proc/cpuinfo", "r");
  char buf;
  int cores = 0;

  /* If the metafile /proc/cpuinfo is not available
   * return 0                                    */
  if (!cpuinfo)
    {
      if (s->debug_level & MPS_DEBUG_MEMORY)
	MPS_DEBUG (s, "Found %d cores on this system", cores);
      return cores;
    }

  /* Check for newlines in /proc/cpuinfo, that should correspond
   * to logical cores.                                        */
  while ((buf = fgetc (cpuinfo)) != EOF)
    {
      if (buf == '\n')
        if (fgetc (cpuinfo) == '\n')
          cores++;
    }

  fclose (cpuinfo);

  if (s->debug_level & MPS_DEBUG_MEMORY)
    MPS_DEBUG (s, "Found %d cores on this system", cores);
  return cores;
}


/**
 * @brief Create a new mps_thread_job_queue that can
 * handle at most max_iter iterations for n_roots roots.
 */
mps_thread_job_queue *
mps_thread_job_queue_new (mps_status * s)
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
mps_thread_job_queue_next (mps_status * s, mps_thread_job_queue * q)
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
  int free_count;

  while (thread->alive)
    {
      /* Start by locking the busy mutex and wait for start condition. This
       * will unlock the mutex and wait for the thread to be signaled. */
      /* printf ("(thread %p) Finished...\n", thread); fflush(stdout); */

      pthread_mutex_lock (&thread->busy_mutex);

      sem_post (&thread->pool->free_count);
      sem_getvalue (&thread->pool->free_count, &free_count);
      thread->busy = false;

      /* printf("(thread %p) Now semaphore value is %d\n", thread, free_count);  */

      pthread_mutex_lock (&thread->pool->free_count_changed_mutex);
      pthread_cond_signal (&thread->pool->free_count_changed_cond);
      pthread_mutex_unlock (&thread->pool->free_count_changed_mutex);

      pthread_cond_wait  (&thread->start_condition, &thread->busy_mutex);
      pthread_mutex_unlock (&thread->busy_mutex);

      if (thread->alive)
	{
	  thread->work (thread->args);
	}

      /* printf("(thread %p) Realmente finito\n", thread); fflush(stdout); */
    }

  /* printf("Mi hanno buttato fuori\n"); fflush(stdout); */

  pthread_mutex_lock (&thread->pool->free_count_changed_mutex);
  pthread_mutex_lock (&thread->busy_mutex);

  sem_post (&thread->pool->free_count);
  thread->busy = false;

  pthread_mutex_unlock (&thread->busy_mutex);
  pthread_cond_signal (&thread->pool->free_count_changed_cond);
  pthread_mutex_unlock (&thread->pool->free_count_changed_mutex);
  pthread_exit (NULL);

  return NULL;
}

/**
 * @brief Start the thread mainloop.
 */
void
mps_thread_start_mainloop (mps_status * s, mps_thread * thread)
{
  pthread_create (thread->thread, NULL, &mps_thread_mainloop, thread);
}

/**
 * @brief Limit the maximum number of threads that can be used in the thread pool.
 */
void mps_thread_pool_set_concurrency_limit (mps_status * s, mps_thread_pool * pool, 
					    unsigned int concurrency_limit)
{
  int i;
  long int l_cl;

  if (!pool)
    pool = s->pool;

  if (pool->n < concurrency_limit)
    concurrency_limit = pool->n;

  /* We need to keep some threads occupied with nothing to do */  
  if (pool->concurrency_limit == 0 && concurrency_limit == 0)
    return;

  /* Update concurrency magic values */
  pool->concurrency_limit = (pool->concurrency_limit == 0) ? pool->n : pool->concurrency_limit;
  l_cl = (concurrency_limit == 0) ? pool->n : concurrency_limit;

  for (i = 0; i < pool->concurrency_limit - l_cl; i++)
    sem_wait (&pool->free_count);

  for (i = 0; i < l_cl - (long int) pool->concurrency_limit; i++)
    sem_post (&pool->free_count);

  pool->concurrency_limit = concurrency_limit;
}

void
mps_thread_pool_assign (mps_status * s, mps_thread_pool * pool, 
			mps_thread_work work, void * args)
{
  if (!pool)
    pool = s->pool;

  /* Assign work to the first free thread in the pool */
  mps_thread * thread = pool->first;

  /* Wait for free threads in the pool */
  /* int valp; */
  /* int i = 0; */
  /* sem_getvalue (&pool->free_count, &valp); */
  /* printf("Sem = %d\n", valp); */
  sem_wait (&pool->free_count);

  while (thread != NULL)
    {
      pthread_mutex_lock (&thread->busy_mutex);
      if (thread->busy == false)
	{
	  /* printf("Assigning to thread %p\n", thread); fflush(stdout);  */
	  thread->work = work;
	  thread->args = args;

	  thread->busy = true;
	  pthread_cond_signal (&thread->start_condition);
	  pthread_mutex_unlock (&thread->busy_mutex);

	  return;
	}
      else 
	{
	  pthread_mutex_unlock (&thread->busy_mutex);
	  thread = thread->next;
	}
    }

  /* printf ("Non ho trovato il thread libero\n");  */
}

/**
 * @brief Wait for a thread pool to complete its jobs.
 */
void
mps_thread_pool_wait (mps_status * s, mps_thread_pool * pool)
{
  int value;
  
  if (!pool)
    pool = s->pool;
  
  long int threads_to_wait = (pool->concurrency_limit == 0) ? pool->n : pool->concurrency_limit;

  do
    {
      pthread_mutex_lock (&pool->free_count_changed_mutex);
      sem_getvalue (&pool->free_count, &value);

      if (value != threads_to_wait)
	pthread_cond_wait (&pool->free_count_changed_cond, &pool->free_count_changed_mutex);

      pthread_mutex_unlock (&pool->free_count_changed_mutex);

    } while (value != threads_to_wait);

}

/**
 * @brief Allocate a new <code>mps_thread</code> and start its mainloop.
 */
mps_thread * 
mps_thread_new (mps_status * s, mps_thread_pool * pool)
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
  thread->busy = true;

  /* Start the thread mainloop */
  mps_thread_start_mainloop (s, thread);

  return thread;
}

/**
 * @brief Free a thread asking it to stop.
 */
void
mps_thread_free (mps_status * s, mps_thread * thread)
{
  /* Wait for the thread to finish its work, if it is doing something */
  /* pthread_mutex_lock (&thread->busy_mutex); */
  /* pthread_mutex_unlock (&thread->busy_mutex); */

  thread->alive = false;

  /* Start the thread, if it is not running */
  pthread_mutex_lock (&thread->busy_mutex);
  pthread_cond_signal (&thread->start_condition);
  pthread_mutex_unlock (&thread->busy_mutex);

  pthread_join (*thread->thread, NULL);

  free (thread->thread);
  free (thread);
}

/**
 * @brief Create a new thread and add it to the specified thread pool.
 */
void
mps_thread_pool_insert_new_thread (mps_status * s, mps_thread_pool * pool)
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
mps_thread_pool_new (mps_status * s, int n_threads)
{
  mps_thread_pool * pool = mps_new (mps_thread_pool);
  int threads = mps_thread_get_core_number (s); 
  int i;

  if (n_threads != 0)
    threads = n_threads;

  pool->n = 0;
  pool->first = NULL;

  pthread_mutex_init (&pool->free_count_changed_mutex, NULL);
  pthread_cond_init (&pool->free_count_changed_cond, NULL);
  
  sem_init (&pool->free_count, 0, 0);
  
  for (i = 0; i < threads; i++) 
    mps_thread_pool_insert_new_thread (s, pool); 

  pool->concurrency_limit = 0;

  mps_thread_pool_wait (s, pool);

  return pool;
}

/**
 * @brief Free a thread pool and all its threads, waiting for them
 * to terminate.
 */
void 
mps_thread_pool_free (mps_status * s, mps_thread_pool * pool)
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

  free (pool);
}

int mps_thread_get_id (mps_status * s, mps_thread_pool * pool)
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
