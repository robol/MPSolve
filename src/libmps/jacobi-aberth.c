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
  mps_approximation * root;
  cplx_t * correction;
};

static void * 
__mps_fjacobi_aberth_step_worker (void * data_ptr)
{
  struct __mps_fjacobi_aberth_step_data *data = (struct __mps_fjacobi_aberth_step_data *) data_ptr;

  cplx_t abcorr, corr;

  mps_context * ctx = data->ctx;
  mps_approximation * root = data->root;
  mps_polynomial * p = data->p;

  mps_polynomial_fnewton (ctx, p, root, corr);

  if (root->approximated)
    root->again = false;

  if (root->again)
  {
    mps_faberth (ctx, root, abcorr);
    cplx_mul_eq (abcorr, corr);
    cplx_sub (abcorr, cplx_one, abcorr);

    if (cplx_check_fpe (abcorr))
    {
      root->again = false;
      root->status = MPS_ROOT_STATUS_NOT_FLOAT;
    }

    if (cplx_eq_zero (abcorr))
      root->again = false;
    else 
      cplx_div (corr, corr, abcorr);

    cplx_set (*data->correction, corr);
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

  cplx_t * corrections = mps_newv (cplx_t, ctx->n);

  for (i = 0; i < ctx->n; i++)
    {
      if (ctx->root[i]->again)
      {
        struct __mps_fjacobi_aberth_step_data * data = mps_new (struct __mps_fjacobi_aberth_step_data);

        data->ctx = ctx;
        data->p = p;
        data->root = ctx->root[i];
        data->correction = &corrections[i];

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
          cplx_sub_eq (ctx->root[i]->fvalue, corrections[i]);
          ctx->root[i]->frad += cplx_mod (corrections[i]);
          again = true;
        }
    }

  cplx_vfree (corrections);

  return again;
}

/**
 * @brief Perform a packet of Aberth iterations on the approximation
 * to the roots of p. 
 *
 * @param ctx Current MPSolve context.
 * @param p The polynomial whose roots should be approximated.
 * @param just regenrated true if this packet is the first following a regeneration.
 *
 * @return The number of approximated roots. 
 */
int
mps_faberth_packet (mps_context * ctx, mps_polynomial * p, mps_boolean just_regenerated)
{
  int iterations = 0, i = 0, approximated_roots = 0, packet = 0, root_neighborhood_roots = 0;
  int it_threshold = ctx->n;

#ifndef DISABLE_DEBUG
  clock_t *my_clock = mps_start_timer ();
#endif

  for (i = 0; i < ctx->n; i++)
    {
      if (MPS_ROOT_STATUS_IS_APPROXIMATED (ctx->root[i]->status))
        ctx->root[i]->approximated = true;

      if (!ctx->root[i]->again)
        it_threshold--;
    }

  do 
    {
      packet++;

      if (ctx->debug_level & MPS_DEBUG_APPROXIMATIONS)
        MPS_DEBUG (ctx, "Carrying out a packet of floating point Aberth iterations (packet = %d)", packet);

    } while (mps_fjacobi_aberth_step (ctx, p, &iterations) && packet <= ctx->max_it);

  MPS_DEBUG_WITH_INFO (ctx, "Performed %d iterations in floating point", iterations);

  /* Check if we need to get higher precision for the roots */
  ctx->best_approx = true; 
  for (i = 0; i < ctx->n; i++) 
    { 
      if (! (ctx->root[i]->approximated || MPS_ROOT_STATUS_IS_COMPUTED (ctx->root[i]->status))) 
       ctx->best_approx = false; 
      if (ctx->root[i]->approximated) 
       approximated_roots++; 
      if (!ctx->root[i]->again) 
       root_neighborhood_roots++; 
    }

  MPS_DEBUG_WITH_INFO (ctx, "%d roots are approximated within the current precision", approximated_roots);
  MPS_DEBUG_WITH_INFO (ctx,"%d roots are in the root neighborhood", root_neighborhood_roots);

#ifndef DISABLE_DEBUG
  ctx->fp_iteration_time += mps_stop_timer (my_clock);
#endif

  return root_neighborhood_roots;
}

struct __mps_djacobi_aberth_step_data {
  mps_context * ctx;
  mps_polynomial * p;
  mps_approximation * root;
  cdpe_t * aberth_correction;
};

static void * 
__mps_djacobi_aberth_step_worker (void * data_ptr)
{
  struct __mps_djacobi_aberth_step_data *data = (struct __mps_djacobi_aberth_step_data *) data_ptr;
  cdpe_t abcorr;

  mps_context * ctx = data->ctx;
  mps_approximation * root = data->root;
  mps_polynomial * p = data->p;

  mps_polynomial_dnewton (ctx, p, root, *data->aberth_correction);

  if (root->approximated)
    root->again = false;  

  if (root->again)
  {
    mps_daberth (ctx, root, abcorr);
    cdpe_mul_eq (abcorr, *data->aberth_correction);
    cdpe_sub (abcorr, cdpe_one, abcorr);

    if (! cdpe_eq_zero (abcorr))
      cdpe_div (*data->aberth_correction, *data->aberth_correction, abcorr);
    else
      root->again = false;
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
        data->root = ctx->root[i];
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
          rdpe_t correction_module;
          again = true;

          cdpe_sub_eq (ctx->root[i]->dvalue, daberth_corrections[i]);
          cdpe_mod (correction_module, daberth_corrections[i]);
          rdpe_add_eq (ctx->root[i]->drad, correction_module);
        }
    }

  cdpe_vfree (daberth_corrections);
  return again;
}

/**
 * @brief Perform a packet of Aberth iterations on the approximation
 * to the roots of p. 
 *
 * @param ctx Current MPSolve context.
 * @param p The polynomial whose roots should be approximated.
 * @param just regenrated true if this packet is the first following a regeneration.
 */
int
mps_daberth_packet (mps_context * ctx, mps_polynomial * p, mps_boolean just_regenerated)
{
  int iterations = 0, i = 0, approximated_roots = 0, packet = 0, root_neighborhood_roots = 0;
  int it_threshold = ctx->n;

#ifndef DISABLE_DEBUG
  clock_t *my_clock = mps_start_timer ();
#endif  

  for (i = 0; i < ctx->n; i++)
    {
      if (MPS_ROOT_STATUS_IS_APPROXIMATED (ctx->root[i]->status))
        ctx->root[i]->approximated = true;

      if (!ctx->root[i]->again)
        it_threshold--;
    }

  do 
    {
      packet++;

      if (ctx->debug_level & MPS_DEBUG_APPROXIMATIONS)
        MPS_DEBUG (ctx, "Carrying out a packet of CDPE Aberth iterations (packet = %d)", packet);

    } while (mps_djacobi_aberth_step (ctx, p, &iterations));

  MPS_DEBUG_WITH_INFO (ctx, "Performed %d iterations in CDPE", iterations);

  /* Check if we need to get higher precision for the roots */
   ctx->best_approx = true; 
   for (i = 0; i < ctx->n; i++) 
     { 
       if (! (ctx->root[i]->approximated || MPS_ROOT_STATUS_IS_COMPUTED (ctx->root[i]->status))) 
        ctx->best_approx = false; 
       if (ctx->root[i]->approximated) 
        approximated_roots++; 
       if (!ctx->root[i]->again) 
        root_neighborhood_roots++; 
     }

  MPS_DEBUG_WITH_INFO (ctx, "%d roots are approximated within the current precision", approximated_roots);
  MPS_DEBUG_WITH_INFO (ctx,"%d roots are in the root neighborhood", root_neighborhood_roots);

#ifndef DISABLE_DEBUG
  ctx->fp_iteration_time += mps_stop_timer (my_clock);
#endif  

  return root_neighborhood_roots;
}



struct __mps_mjacobi_aberth_step_data {
  mps_context * ctx;
  mps_polynomial * p;
  mps_approximation * root;
  mpc_t * aberth_correction;
};

static void * 
__mps_mjacobi_aberth_step_worker (void * data_ptr)
{
  struct __mps_mjacobi_aberth_step_data *data = (struct __mps_mjacobi_aberth_step_data *) data_ptr;
  mpc_t corr, abcorr;

  mps_context * ctx = data->ctx;
  mps_approximation * root = data->root;
  mps_polynomial * p = data->p;

  mpc_init2 (corr, ctx->mpwp);
  mpc_init2 (abcorr, ctx->mpwp);

  mps_polynomial_mnewton (ctx, p, root, corr);

  if (root->approximated)
    root->again = false;

  if (root->again)
  {
    mps_maberth (ctx, root, abcorr);
    mpc_mul_eq (abcorr, corr);
    mpc_ui_sub (abcorr, 1U, 0U, abcorr);

    if (! mpc_eq_zero (abcorr)) 
      mpc_div (abcorr, corr, abcorr); 
    else
      root->again = false;

    mpc_set (*data->aberth_correction, abcorr);
  }

  mpc_clear (corr);
  mpc_clear (abcorr);

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
mps_mjacobi_aberth_step (mps_context * ctx, mps_polynomial * p, int * nit)
{
  mpc_t * maberth_corrections = NULL;
  mps_boolean again = false;
  int i = 0;

  maberth_corrections = mpc_valloc (ctx->n);
  mpc_vinit2 (maberth_corrections, ctx->n, ctx->mpwp);

  for (i = 0; i < ctx->n; i++)
    {
      if (ctx->root[i]->again)
      {
        struct __mps_mjacobi_aberth_step_data * data = mps_new (struct __mps_mjacobi_aberth_step_data);

        data->ctx = ctx;
        data->p = p;
        data->root = ctx->root[i];
        data->aberth_correction = &maberth_corrections[i];

        /* The worker is in charge of freeing the data that we have allocated
         * here, so we can't ignore this issue in the main thread. */
        mps_thread_pool_assign (ctx, ctx->pool, __mps_mjacobi_aberth_step_worker, 
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
          rdpe_t correction_module;
          again = true;

          mpc_sub_eq (ctx->root[i]->mvalue, maberth_corrections[i]);
          mpc_rmod (correction_module, maberth_corrections[i]);
          rdpe_add_eq (ctx->root[i]->drad, correction_module);
        }
    }

  mpc_vclear (maberth_corrections, ctx->n);
  mpc_vfree (maberth_corrections);

  return again;
}

/**
 * @brief Perform a packet of Aberth iterations on the approximation
 * to the roots of p. 
 *
 * @param ctx Current MPSolve context.
 * @param p The polynomial whose roots should be approximated.
 * @param just regenrated true if this packet is the first following a regeneration.
 */
int
mps_maberth_packet (mps_context * ctx, mps_polynomial * p, mps_boolean just_regenerated)
{
  int iterations = 0, i = 0, approximated_roots = 0, packet = 0, root_neighborhood_roots = 0;
  int it_threshold = ctx->n;

#ifndef DISABLE_DEBUG
  clock_t *my_clock = mps_start_timer ();
#endif  

  for (i = 0; i < ctx->n; i++)
    {
      if (MPS_ROOT_STATUS_IS_APPROXIMATED (ctx->root[i]->status))
        ctx->root[i]->approximated = true;

      if (!ctx->root[i]->again)
        it_threshold--;
    }

  do 
    {
      packet++;

      if (ctx->debug_level & MPS_DEBUG_APPROXIMATIONS)
        MPS_DEBUG (ctx, "Carrying out a packet of multiprecision Aberth iterations (packet = %d)", packet);

    } while (mps_mjacobi_aberth_step (ctx, p, &iterations));

  MPS_DEBUG_WITH_INFO (ctx, "Performed %d iterations in multiprecision", iterations);

  /* Check if we need to get higher precision for the roots */
   ctx->best_approx = true; 
   for (i = 0; i < ctx->n; i++) 
     { 
       if (! (ctx->root[i]->approximated || MPS_ROOT_STATUS_IS_COMPUTED (ctx->root[i]->status))) 
        ctx->best_approx = false;
       if (ctx->root[i]->approximated) 
        approximated_roots++; 
       if (!ctx->root[i]->again) 
        root_neighborhood_roots++; 
     }

  MPS_DEBUG_WITH_INFO (ctx, "%d roots are approximated within the current precision", approximated_roots);
  MPS_DEBUG_WITH_INFO (ctx,"%d roots are in the root neighborhood", root_neighborhood_roots);

#ifndef DISABLE_DEBUG
  ctx->fp_iteration_time += mps_stop_timer (my_clock);
#endif  

  return root_neighborhood_roots;
}

