/*
 * threading.h
 *
 *  Created on: 19/mag/2011
 *      Author: leonardo
 */

#ifdef __cplusplus
extern "C" {
#endif

#ifndef THREADING_H_
#define THREADING_H_

#include <pthread.h>
#include <semaphore.h>
#include <mps/core.h>


#define MPS_THREAD_FULL -1
#define MPS_THREAD_AVAIL 0

typedef enum {
  MPS_FNEWTON_THREAD,
  MPS_DNEWTON_THREAD,
  MPS_MNWETON_THREAD,
} mps_thread_job_type;

/**
 * @brief A mps_thread_pool threadpool.
 */
typedef struct {
  pthread_t* threads;
  mps_boolean* ready;
  pthread_cond_t full;
  pthread_mutex_t full_mutex;
  pthread_mutex_t ready_mutex;
  pthread_mutex_t aberth_mutex;
  int n;
  int nzeros;
} mps_thread_pool;

typedef struct {
  int* nzeros;
  int* it;
  mps_status* s;
  int* iter;
  int thread;
  int n_threads;
  mps_boolean* excep;
} mps_thread_worker_data;

void
mps_thread_fpolzer(mps_status* s, int* nit, mps_boolean* excep);

void
mps_thread_mpolzer(mps_status* s, int *nit, mps_boolean *excep);


#ifdef __cplusplus
}
#endif

#endif /* THREADING_H_ */
