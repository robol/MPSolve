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
#include <mps/core.h>

#define MPS_THREAD_JOB_EXCEP -1

typedef struct {
  int i;
  int iter;
  int i_clust;
} mps_thread_job;

typedef struct {
  unsigned int max_iter;
  unsigned int n_roots;
  int iter;
  int i;
  int i_clust;
  pthread_mutex_t mutex;
} mps_thread_job_queue;

typedef struct {
  int* nzeros;
  int* it;
  mps_status* s;
  int thread;
  int n_threads;
  int* thread_completed;
  mps_boolean* excep;
  pthread_mutex_t* aberth_mutex;
  pthread_mutex_t* roots_mutex;
  mps_thread_job_queue* queue;
} mps_thread_worker_data;

void
mps_thread_fpolzer(mps_status* s, int* nit, mps_boolean* excep);

void
mps_thread_mpolzer(mps_status* s, int *nit, mps_boolean *excep);


#ifdef __cplusplus
}
#endif

#endif /* THREADING_H_ */
