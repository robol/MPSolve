#include <pthread.h>
#include <mps/mps.h>

/**
 * @file
 * @brief Implementation of some thread-safe types that can be easily used
 * with the macro MPS_LOCK() and MPS_UNLOCK().
 */

#ifndef __MPS_MT_TYPES
#define __MPS_MT_TYPES

#define MPS_LOCK(x) (pthread_mutex_lock (&(x).mutex));
    
#define MPS_UNLOCK(x) (pthread_mutex_unlock (&(x).mutex))

#define MPS_INIT_LOCK(x) (pthread_mutex_init (&(x).mutex, NULL))

  /**
   * @brief A thread safe version of mps_boolean. 
   *
   * Must be accessed using the macro MPS_LOCK (x) and
   * MPS_UNLOCK (x).
   */
  struct mps_boolean_mt {
    mps_boolean value;
    pthread_mutex_t mutex;
  };

  /**
   * @brief A thread safe version of mps_boolean. 
   *
   * Must be accessed using the macro MPS_LOCK (x) and
   * MPS_UNLOCK (x).
   */
  struct mps_long_int_mt {
    long int value;
    pthread_mutex_t mutex;
  };

#ifndef __cplusplus
    typedef struct mps_boolean_mt mps_boolean_mt;
    typedef struct mps_long_int_mt mps_long_int_mt;

#endif
#endif
