/*
 * threading.h
 *
 *  Created on: 19/mag/2011
 *      Author: leonardo
 */

/**
 * @file
 * @brief Multithreading iterations for MPSolve.
 */

#ifdef __cplusplus
extern "C"
{
#endif

#ifndef MPS_THREADING_H_
#define MPS_THREADING_H_

#include <pthread.h>
#include <semaphore.h>
#include <mps/mps.h>

#define MPS_THREAD_JOB_EXCEP -1

  /**
   * @brief A generic routine that can be performed by a <code>mps_thread</code>.
   */
  typedef void * (*mps_thread_work) (void *);

  /**
   * @brief A new job for <code>mps_thread_fsolve()</code>,
   * <code>mps_thread_dsolve()</code> or <code>mps_thread_msolve()</code>.
   *
   */
  struct mps_thread_job
  {
    /**
     * @brief The index if the root to iterate on.
     */
    int i;

    /**
     * @brief The iteration that will be performed on this root.
     */
    int iter;

    /**
     * @brief cluster_item The cluster element of <code>s->clusterization
     * that we are iterating on.
     */
    mps_cluster_item * cluster_item;
  };


  /**
   * @brief Struct holding a job queue.
   *
   * This structure can be used to coordinate the work in the different
   * thread during multithread computation in MPSolve.
   *
   * It must be allocated using <code>mps_thread_job_queue_new()</code>
   * and freed with <code>mps_thread_job_queue_free()</code>.
   * A new job can be requested with the routine
   * <code>mps_thread_job_queue_next()</code>.
   *
   * @see mps_thread_job_queue_next()
   */
  struct mps_thread_job_queue
  {

    /**
     * @brief Maximum number of iteration to perform before
     *  raising an exeption.
     */
    unsigned int max_iter;

    /**
     * @brief Number of the roots of this problem (i.e. degree of
     * the polynomial).
     */
    unsigned int n_roots;

    /**
     * @brief Iterations that is being performed right now.
     */
    int iter;

    /**
     * @brief Next root to iterate on.
     */
    mps_root * root;

    /**
     * @brief Element of <code>s->clusterization</code> that
     * we are iterating on.
     */
    mps_cluster_item * cluster_item;

    /**
     * @brief Internal mutex of the queue used to guarantee
     * exclusive access.
     */
    pthread_mutex_t mutex;
  };

  /**
   * @brief Data packed to be passed to a new thread that will
   * perform floating point, dpe or multiprecision iterations.
   */
  struct mps_thread_worker_data
  {

    /**
     * @brief Pointer to the integer that holds the number of zeros
     *  computed until now.
     */
    int *nzeros;

    /**
     * @brief Pointer to the integer that holds the number of iterations
     *  performed until now.
     */
    int *it;

    /**
     * @brief  The pointer to the <code>mps_status</code> struct.
     */
    mps_status *s;

    /**
     * @brief The index of this thread.
     */
    int thread;

    /**
     * @brief The total number of threads.
     */
    int n_threads;

    /**
     * @brief Pointer to the boolean excep value. Setting this to true
     *  cause the iteration to enter exception state.
     *
     *  If this state is reached all threads returns because no more
     *  iteration are needed / useful.
     */
    mps_boolean *excep;

    /**
     * @brief Array of <code>n</code> mutexes where <code>n = s->n</code>, i.e.
     * is the total number of roots of the polynomial.
     *
     * The mutex in position <code>i</code> gets locked when a
     * thread needs to read and/or write from/to
     * the i-th root.
     */
    pthread_mutex_t *aberth_mutex;

    /**
     * @brief Global aberth mutex used to coordinate all aberth
     * computations. 
     */
    pthread_mutex_t *global_aberth_mutex;

    /**
     * @brief Array of <code>n</code> mutexes that gets locked when a thread
     * start to iterate over a root. This is done to ensure that only a thread
     * at a time is iterating over a root.
     */
    pthread_mutex_t *roots_mutex;

    /**
     * @brief Pointer to the <code>mps_thread_job_queue</code> that the thread
     * may query for other work.
     */
    mps_thread_job_queue *queue;
  };

  /**
   * @brief A thread that is part of a thread pool.
   */
  struct mps_thread {

    /**
     * @brief Pool of which this thread is part.
     */
    mps_thread_pool * pool;

    /**
     * @brief The pthread_t assigned to the worked.
     */
    pthread_t * thread;

    /**
     * @brief The next thread in the pool, or NULL if this
     * is the last thread contained in it.
     */
    mps_thread * next;

    /**
     * @brief The data assigned to this thread, that sets 
     * the worker that he has to do.
     */
    mps_thread_worker_data * data;

    /**
     * @brief True if the thread is busy.
     */
    mps_boolean busy;

    /**
     * @brief Busy mutex of the thread. This is locked when the thread
     * is doing something, se we can emulate a join on it by
     * try to lock and unlock this mutex.
     */
    pthread_mutex_t busy_mutex;

    /**
     * @brief Condition that allow the thread to run. Before the thread
     * finish the busy state (unlocking the busy mutex) or when it is
     * created, it waits for the start condition to be true before
     * doing anything.
     */
    pthread_cond_t start_condition;

    /**
     * @brief A boolean value that is true if the thread must continue to
     * poll, or false if it is required to exit. Since the thread
     * may be waiting for work a call to pthread_cond_signal on
     * start condition may be required to make it exit after setting
     * this variable.
     */
    mps_boolean alive;

    /**
     * @brief The routine that must be called when the thread starts. 
     */
    mps_thread_work work;

    /**
     * @brief The argument to be passed to the thread.
     */
    void * args;
  };

  /**
   * @brief A thread pool that contains a set of <code>mps_thread</code>
   * and allow to manage them as a set of worker. 
   */
  struct mps_thread_pool {
    
    /**
     * @brief The numer of thread in the thread pool.
     */
    unsigned int n;

    /**
     * @brief A pointer to the first thread in the thread pool.
     */
    mps_thread * first;

    /**
     * @brief A semaphore counting the numer of free threads.
     */
    sem_t free_count;

    /** 
     * @brief Mutex associated with the semaphore changed condition. 
     */ 
    pthread_mutex_t free_count_changed_mutex; 
    
    /** 
     * @brief Condition that is signaled when the semaphore changes 
     * its value. 
     */ 
    pthread_cond_t free_count_changed_cond; 
  };

  /* EXPORTED ROUTINES */

  void * mps_thread_mainloop (void * thread_ptr);

  void mps_thread_start_mainloop (mps_status * s, mps_thread * thread);

  mps_thread * mps_thread_new (mps_status * s, mps_thread_pool * pool);

  void mps_thread_free (mps_status * s, mps_thread * thread);

  void mps_thread_pool_assign (mps_status * s, mps_thread_pool * pool, mps_thread_work work, void * args);

  void mps_thread_pool_insert_new_thread (mps_status * s, mps_thread_pool * pool);

  void mps_thread_pool_wait (mps_status * s, mps_thread_pool * pool);

  mps_thread_pool * mps_thread_pool_new (mps_status * s);

  void mps_thread_pool_free (mps_status * s, mps_thread_pool * pool);

  mps_thread_job_queue * mps_thread_job_queue_new (mps_status * s);

  void mps_thread_job_queue_free (mps_thread_job_queue * q);

  mps_thread_job mps_thread_job_queue_next (mps_status * s, mps_thread_job_queue * q);

  void mps_thread_fpolzer (mps_status * s, int *nit, mps_boolean * excep);

  void mps_thread_mpolzer (mps_status * s, int *nit, mps_boolean * excep);

  void mps_thread_dpolzer (mps_status * s, int *nit, mps_boolean * excep);

  int mps_thread_get_core_number (mps_status * s);

  /* MACROS */

  /**
   * @brief Get a pointer to an array of n+2 booleans
   * that is local to the thread.
   */
#define mps_thread_get_spar2(s, n_thread) (s->spar2 + (s->deg + 2) * (n_thread))

  /**
   * @brief Get a pointer to an array of n+1 multiprecision
   * that is local to the thread.
   */
#define mps_thread_get_mfpc2(s, n_thread) (s->mfpc2 + (s->deg + 1) * (n_thread))

  /**
   * @brief Get a pointer to an array of n+2 DPE
   * that is local to the thread.
   */
#define mps_thread_get_dap2(s, n_thread) (s->dap2 + (s->deg + 2) * (n_thread))


#ifdef __cplusplus
}
#endif

#endif                          /* MPS_THREADING_H_ */
