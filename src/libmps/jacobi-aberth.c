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
  cplx_t * aberth_correction;
};

static void * 
__mps_fjacobi_aberth_step_worker (void * data_ptr)
{
  struct __mps_fjacobi_aberth_step_data *data = (struct __mps_fjacobi_aberth_step_data *) data_ptr;
  cplx_t abcorr;

  mps_polynomial_fnewton (data->ctx, data->p, data->ctx->root[data->i], *data->aberth_correction);
  mps_faberth (data->ctx, data->i, abcorr);

  cplx_mul_eq (abcorr, *data->aberth_correction);
  cplx_sub (abcorr, cplx_one, abcorr);
  cplx_div (*data->aberth_correction, *data->aberth_correction, abcorr);

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
mps_fjacobi_aberth_step (mps_context * ctx, mps_polynomial * p)
{
  cplx_t * faberth_corrections = NULL;
  mps_boolean again = false;
  int i = 0;

  faberth_corrections = cplx_valloc (ctx->n);

  for (i = 0; i < ctx->n; i++)
    {
      if (ctx->root[i]->again)
      {
        struct __mps_fjacobi_aberth_step_data * data = mps_new (struct __mps_fjacobi_aberth_step_data);

        data->ctx = ctx;
        data->p = p;
        data->i = i;
        data->aberth_correction = &faberth_corrections[i];

        /* The worker is in charge of freeing the data that we have allocated
         * here, so we can't ignore this issue in the main thread. */
        mps_thread_pool_assign (ctx, ctx->pool, __mps_fjacobi_aberth_step_worker, 
          data);
      }
    }

  mps_thread_pool_wait (ctx, ctx->pool);

  /* Update approximations with the new corrections */
  for (i = 0; i < ctx->n; i++)
    {
      cplx_sub_eq (ctx->root[i]->fvalue, faberth_corrections[i]);
      again |= ctx->root[i]->again;
    }

  cplx_vfree (faberth_corrections);
  return again;
}

/**
 * @brief Perform a packet of Aberth iterations on the approximation
 * to the roots of p. 
 *
 * @params ctx Current MPSolve context.
 * @params p The polynomial whose roots should be approximated.
 */
void
mps_faberth_packet (mps_context * ctx, mps_polynomial * p)
{
  double * radii = double_valloc (ctx->n);
  int iteration = 0;

  do 
    {
      iteration++;

      if (ctx->debug_level & MPS_DEBUG_APPROXIMATIONS)
        MPS_DEBUG (ctx, "Carrying out a packet of Aberth iterations (packet = %d)", iteration);

    } while (mps_fjacobi_aberth_step (ctx, p));

  mps_fradii (ctx, radii);
  mps_fcluster (ctx, radii, 2 * ctx->n);

  mps_fmodify (ctx, false);

  double_vfree (radii);
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

  mps_polynomial_dnewton (data->ctx, data->p, data->ctx->root[data->i], *data->aberth_correction);
  mps_daberth (data->ctx, data->i, abcorr);

  cdpe_mul_eq (abcorr, *data->aberth_correction);
  cdpe_sub (abcorr, cdpe_one, abcorr);
  cdpe_div (*data->aberth_correction, *data->aberth_correction, abcorr);

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
mps_djacobi_aberth_step (mps_context * ctx, mps_polynomial * p)
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
      }
    }

  mps_thread_pool_wait (ctx, ctx->pool);

  /* Update approximations with the new corrections */
  for (i = 0; i < ctx->n; i++)
    {
      cdpe_sub_eq (ctx->root[i]->dvalue, daberth_corrections[i]);
      again |= ctx->root[i]->again;
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
void
mps_daberth_packet (mps_context * ctx, mps_polynomial * p)
{
  rdpe_t * radii = rdpe_valloc (ctx->n);
  int iteration = 0;

  do 
    {
      iteration++;

      if (ctx->debug_level & MPS_DEBUG_APPROXIMATIONS)
        MPS_DEBUG (ctx, "Carrying out a packet of Aberth iterations (packet = %d)", iteration);

    } while (mps_djacobi_aberth_step (ctx, p));

  mps_dradii (ctx, radii);
  mps_dcluster (ctx, radii, 2 * ctx->n);

  mps_dmodify (ctx, false);

  rdpe_vfree (radii);
}