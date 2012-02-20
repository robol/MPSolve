/*
 * mps_interface.c
 *
 *  Created on: 05/apr/2011
 *      Author: Leonardo Robol <robol@poisson.phc.unipi.it>
 */


#include <mps/mps.h>
#include <mps/link.h>
#include <gmp.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#ifdef MPS_CATCH_FPE
#include <fenv.h>
int feenableexcept(int excepts);
#endif

/**
 * @brief Call the real polynomial (or secular equation, or whatever) solver
 * and do the computation.
 *
 * The algorithm used must be selected before this call with <code>mps_select_algorithm</code>
 * and the data (the coefficients, or whatever the algorithm may require) should be provided
 * after that.
 *
 * Roots can then be obtained with the functions <code>mps_status_get_roots_*</code>
 *
 */
void
mps_mpsolve (mps_status * s)
{
#ifdef MPS_CATCH_FPE
  feenableexcept (FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW); 
#endif

  (*s->mpsolve_ptr) (s);
}

void*
mps_caller (mps_status * s)
{
  s->mpsolve_ptr (s);
  return NULL;
}

mps_async_handle *
mps_mpsolve_async (mps_status * s)
{
  #ifdef MPS_CATCH_FPE
  feenableexcept (FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW); 
#endif

  mps_thread_pool * private_pool = mps_thread_pool_new (s, 1);
  mps_thread_pool_assign (s, private_pool, (mps_thread_work) mps_caller, s);
  
  return (mps_async_handle*) private_pool;
}

void
mps_mpsolve_wait (mps_status * s, mps_async_handle * handle)
{
  mps_thread_pool_wait (s, (mps_thread_pool*) handle);
  mps_thread_pool_free (s, (mps_thread_pool*) handle);
}

/**
 * @brief Allocator for memory to be used in mpsolve.
 */
void *
mps_malloc (size_t size)
{
  /* fprintf (stderr, "Allocating %lu bytes of memory\n", size); */
  register void *value = malloc (size);
  if (value == 0)
    {
      fprintf (stderr, "virtual memory exhausted");
      exit (1);
    }
  return value;
}

/**
 * @brief Allocate size bytes on the stack
 */
void *
mps_alloca (size_t size)
{
  register void *value = alloca (size);
  return value;
}
