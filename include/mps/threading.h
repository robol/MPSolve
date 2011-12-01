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

  /* FORWARD DECLARATIONS */
  typedef struct  mps_thread_job mps_thread_job;
  typedef struct  mps_thread_job_queue mps_thread_job_queue;
  typedef struct  mps_thread_worker_data mps_thread_worker_data;

#include <pthread.h>
#include <mps/core.h>

#define MPS_THREAD_JOB_EXCEP -1

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
     * @brief Pointo the <code>mps_thread_job_queue</code> that the thread
     * may query for other work.
     */
    mps_thread_job_queue *queue;
  };

  /* EXPORTED ROUTINES */

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
