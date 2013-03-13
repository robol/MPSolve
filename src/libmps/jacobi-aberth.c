/*
 * This file is part of MPSolve 3.0
 *
 * Copyright (C) 2001-2013, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors: 
 *   Leonardo Robol <robol@mail.dm.unipi.it>
 */

#include <mps/mps.h>

struct __mps_fjacobi_aberth_step_data {
  mps_context * ctx;
  mps_polynomial * p;
  int i;
  pthread_mutex_t * mutexes;
};

static void * 
__mps_fjacobi_aberth_step_worker (void * data_ptr)
{
  struct __mps_fjacobi_aberth_step_data *data = (struct __mps_fjacobi_aberth_step_data *) data_ptr;
  cplx_t abcorr, corr;
  int i = data->i;
  mps_approximation * root = data->ctx->root[i];

  mps_polynomial_fnewton (data->ctx, data->p, root, corr);
  mps_faberth_wl (data->ctx, data->i, abcorr, data->mutexes);

  cplx_mul_eq (abcorr, corr);
  cplx_sub (abcorr, cplx_one, abcorr);

  if (cplx_eq_zero (abcorr) || cplx_check_fpe (abcorr))
    root->again = false;
  else 
    cplx_div (corr, corr, abcorr);

  if (root->again)
  {
    pthread_mutex_lock (&data->mutexes[i]);
    cplx_sub_eq (root->fvalue, corr);  
    pthread_mutex_unlock (&data->mutexes[i]);

    root->frad += cplx_mod (corr);
  }

  free (data);
  return NULL;
}

/**
 * @brief Perform a step of Aberth method in Jacobi-style.
 *
 * @param ctx The context in which this instance of MPSolve is running.
 * @param p The polynomial on which Aberth method should be applied.
 * @param approximations The starting approximations for the method. They will be
 *        updated with the new values when the function returns. 
 * @param nit Number of iterations performed in the packet. 
 */
mps_boolean
mps_fjacobi_aberth_step (mps_context * ctx, mps_polynomial * p, int * nit)
{
  mps_boolean again = false;
  int i = 0;
  pthread_mutex_t * mutexes = mps_newv (pthread_mutex_t, ctx->n);

  for (i = 0; i < ctx->n; i++)
    pthread_mutex_init (&mutexes[i], NULL);

  for (i = 0; i < ctx->n; i++)
    {
      if (ctx->root[i]->again)
      {
        struct __mps_fjacobi_aberth_step_data * data = mps_new (struct __mps_fjacobi_aberth_step_data);

        data->ctx = ctx;
        data->p = p;
        data->i = i;
        data->mutexes = mutexes;

        /* The worker is in charge of freeing the data that we have allocated
         * here, so we can't ignore this issue in the main thread. */
        mps_thread_pool_assign (ctx, ctx->pool, 
          __mps_fjacobi_aberth_step_worker, data);

        if (nit)
          (*nit)++;
      }
    }

  mps_thread_pool_wait (ctx, ctx->pool);

  /* Update again */
  for (i = 0; i < ctx->n; i++)
    {
      if (ctx->root[i]->again)
        {
          again = true;
          break;
        }
    }

  cplx_vfree (mutexes);
  return again;
}

/**
 * @brief Perform a packet of Aberth iterations on the approximation
 * to the roots of p. 
 *
 * @params ctx Current MPSolve context.
 * @params p The polynomial whose roots should be approximated.
 *
 * @return The number of approximated roots. 
 */
int
mps_faberth_packet (mps_context * ctx, mps_polynomial * p)
{
  int iterations = 0, i = 0, approximated_roots = 0, packet = 0, root_neighborhood_roots = 0;

  for (i = 0; i < ctx->n; i++)
    {
      if (MPS_ROOT_STATUS_IS_COMPUTED (ctx->root[i]->status))
        ctx->root[i]->again = false;

      if (MPS_ROOT_STATUS_IS_APPROXIMATED (ctx->root[i]->status))
        ctx->root[i]->approximated = true;
    }

  do 
    {
      packet++;

      if (ctx->debug_level & MPS_DEBUG_APPROXIMATIONS)
        MPS_DEBUG (ctx, "Carrying out a packet of Aberth iterations (packet = %d)", packet);

    } while (mps_fjacobi_aberth_step (ctx, p, &iterations) && packet <= ctx->max_it);

  MPS_DEBUG (ctx, "Performed %d iterations in floating point", iterations);

  /* Check if we need to get higher precision for the roots */
   ctx->best_approx = true; 
   for (i = 0; i < ctx->n; i++) 
     { 
       if (!ctx->root[i]->approximated) 
        ctx->best_approx = false; 
       if (ctx->root[i]->approximated) 
        approximated_roots++; 
       if (!ctx->root[i]->again) 
        root_neighborhood_roots++; 
     }

  return root_neighborhood_roots;
}

struct __mps_djacobi_aberth_step_data {
  mps_context * ctx;
  mps_polynomial * p;
  int i;
  cdpe_t * aberth_correction;
};

static void * 
__mps_djacobi_aberth_step_worker (void * data_ptr)
{
  struct __mps_djacobi_aberth_step_data *data = (struct __mps_djacobi_aberth_step_data *) data_ptr;
  cdpe_t abcorr;
  mps_approximation * root = data->ctx->root[data->i];

  mps_polynomial_dnewton (data->ctx, data->p, data->ctx->root[data->i], *data->aberth_correction);
  mps_daberth (data->ctx, data->i, abcorr);

  cdpe_mul_eq (abcorr, *data->aberth_correction);
  cdpe_sub (abcorr, cdpe_one, abcorr);
  cdpe_div (*data->aberth_correction, *data->aberth_correction, abcorr);

  if (root->again)
  {
    rdpe_t correction_module;
    cdpe_sub_eq (root->dvalue, *data->aberth_correction);
    cdpe_mod (correction_module, *data->aberth_correction);
    rdpe_add_eq (root->drad, correction_module);
  }

  free (data);
  return NULL;
}

/**
 * @brief Perform a step of Aberth method in Jacobi-style.
 *
 * @param ctx The context in which this instance of MPSolve is running.
 * @param p The polynomial on which Aberth method should be applied.
 * @param approximations The starting approximations for the method. They will be
 *        updated with the new values when the function returns. 
 */
mps_boolean
mps_djacobi_aberth_step (mps_context * ctx, mps_polynomial * p, int * nit)
{
  cdpe_t * daberth_corrections = NULL;
  mps_boolean again = false;
  int i = 0;

  daberth_corrections = cdpe_valloc (ctx->n);

  for (i = 0; i < ctx->n; i++)
    {
      if (ctx->root[i]->again)
      {
        struct __mps_djacobi_aberth_step_data * data = mps_new (struct __mps_djacobi_aberth_step_data);

        data->ctx = ctx;
        data->p = p;
        data->i = i;
        data->aberth_correction = &daberth_corrections[i];

        /* The worker is in charge of freeing the data that we have allocated
         * here, so we can't ignore this issue in the main thread. */
        mps_thread_pool_assign (ctx, ctx->pool, __mps_djacobi_aberth_step_worker, 
          data);

        if (nit)
          (*nit)++;
      }
    }

  mps_thread_pool_wait (ctx, ctx->pool);

  /* Update again */
  for (i = 0; i < ctx->n; i++)
    {
      if (ctx->root[i]->again)
        {
          again = true;
          break;
        }
    }

  cdpe_vfree (daberth_corrections);
  return again;
}

/**
 * @brief Perform a packet of Aberth iterations on the approximation
 * to the roots of p. 
 *
 * @params ctx Current MPSolve context.
 * @params p The polynomial whose roots should be approximated.
 */
int
mps_daberth_packet (mps_context * ctx, mps_polynomial * p)
{
  int iterations = 0, packet = 0, i, approximated_roots = 0,
    root_neighborhood_roots = 0;

  do 
    {
      packet++;

      if (ctx->debug_level & MPS_DEBUG_APPROXIMATIONS)
        MPS_DEBUG (ctx, "Carrying out a packet of Aberth iterations (packet = %d)", packet);

    } while (mps_djacobi_aberth_step (ctx, p, &iterations));

  /* Check if we need to get higher precision for the roots */
   ctx->best_approx = true; 
   for (i = 0; i < ctx->n; i++) 
     { 
       if (!ctx->root[i]->approximated) 
        ctx->best_approx = false; 
       if (ctx->root[i]->approximated) 
        approximated_roots++; 
       if (!ctx->root[i]->again) 
        root_neighborhood_roots++; 
     }

  return root_neighborhood_roots;
}