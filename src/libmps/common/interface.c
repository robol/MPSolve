/*
 * This file is part of MPSolve 3.2.2
 *
 * Copyright (C) 2001-2020, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Dario Andrea Bini <bini@dm.unipi.it>
 *   Giuseppe Fiorentino <fiorent@dm.unipi.it>
 *   Leonardo Robol <leonardo.robol@unipi.it>
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
int feenableexcept (int excepts);
#endif

#if HAVE_CONFIG_H
# include <config.h>
#endif
#undef malloc
#undef realloc

void *malloc (size_t n);
void *realloc (void * ptr, size_t n);

/**
 * @brief Perform some preliminary checks and setup before starting the real
 * MPSolve loop.
 */
static void
mps_preliminary_setup (mps_context * ctx)
{
  /* Make sure that non thread safe polynomial implementations are handled
   * in a safe way. */
  if (!ctx->active_poly->thread_safe)
    {
      mps_thread_pool_set_concurrency_limit (ctx, NULL, 1);
    }
}

/**
 * @brief Call the real polynomial (or secular equation, or whatever) solver
 * and do the computation.
 *
 * The algorithm used must be selected before this call with <code>mps_select_algorithm</code>
 * and the data (the coefficients, or whatever the algorithm may require) should be provided
 * after that.
 *
 * Roots can then be obtained with the functions <code>mps_context_get_roots_*</code>
 *
 */
void
mps_mpsolve (mps_context * s)
{
#ifdef MPS_CATCH_FPE
  feenableexcept (FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif

  if (mps_context_has_errors (s))
    return;

  mps_preliminary_setup (s);

  (*s->mpsolve_ptr)(s);
}

static void*
mps_caller (mps_context * s)
{
  if (!mps_context_has_errors (s))
    {
      mps_preliminary_setup (s);
      s->mpsolve_ptr (s);
    }

  /* Call user defined callback if available */
  if (s->callback == NULL)
    return NULL;
  else
    return (*s->callback)(s, s->user_data);
}

void
mps_mpsolve_async (mps_context * s, mps_callback callback, void * user_data)
{
#ifdef MPS_CATCH_FPE
  feenableexcept (FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif

  /* Set up callbacks */
  s->callback = callback;
  s->user_data = user_data;

  mps_thread_pool * private_pool = mps_thread_pool_new (s, 1);
  mps_thread_pool_set_strict_async (private_pool, true);
  s->self_thread_pool = private_pool;
  mps_thread_pool_assign (s, private_pool, (mps_thread_work) mps_caller, s);
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
      abort ();
    }
  return value;
}

/**
 * @brief Reallocator for memory used in MPSolve.
 */
void *
mps_realloc (void * pointer, size_t size)
{
  register void *value = realloc (pointer, size);

  if (value == 0 && size > 0)
    {
      fprintf (stderr, "virtual memory exhausted");
      abort ();
    }
  return value;
}
