/*
 * This file is part of MPSolve 3.0
 *
 * Copyright (C) 2001-2012, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors: 
 *   Dario Andrea Bini <bini@dm.unipi.it>
 *   Giuseppe Fiorentino <fiorent@dm.unipi.it>
 *   Leonardo Robol <robol@mail.dm.unipi.it>
 */


#include <mps/mps.h>
#include <math.h>
#include <limits.h>

#ifndef log2
#define log2(x) (log(x) / LOG2)
#endif

static int 
get_approximated_bits (mps_approximation * appr)
{
  rdpe_t module;
  mpc_rmod (module, appr->mvalue);

  return (rdpe_log (module) - rdpe_log (appr->drad)) / LOG2 - 1;
}

static rdpe_t *
evaluate_root_conditioning (mps_context * ctx, mps_polynomial * p, mps_approximation ** appr, int n)
{
  int i;

  rdpe_t * root_conditioning = rdpe_valloc (n);

  for (i = 0; i < n; i++)
    {
      mpc_t value;
      rdpe_t error, module;
      mpc_init2 (value, appr[i]->wp);

      mps_polynomial_meval (ctx, p, appr[i]->mvalue, value, error);
      mpc_rmod (module, value);

      /* Get the relative error of this evaluation */
      if (! rdpe_eq_zero (module))
        rdpe_div_eq (error, module);
      else
        rdpe_set_d (error, DBL_EPSILON * p->degree);

      /* log2(error) + wp - log(n) is a good estimate of log(k) */
      rdpe_set_d (root_conditioning[i], rdpe_log (error) / LOG2 + appr[i]->wp - log2(n));
      rdpe_exp_eq (root_conditioning[i]);

      rdpe_div_eq (root_conditioning[i], appr[i]->drad);

      mpc_clear (value);
    }

  return root_conditioning;
}

static void
improve_root (mps_context * ctx, mps_polynomial * p, mps_approximation * root, long precision)
{
  mpc_t newton_correction;
  rdpe_t corr_mod, epsilon;

  mpc_set_prec (root->mvalue, precision);
  mpc_init2 (newton_correction, precision);

  mps_polynomial_mnewton (ctx, p, root, newton_correction);

  mpc_sub_eq (root->mvalue, newton_correction);
  mpc_rmod (corr_mod, newton_correction);
  rdpe_add_eq (root->drad, corr_mod);

  mpc_rmod (corr_mod, root->mvalue);
  rdpe_set_2dl (epsilon, 1.0, 2 - precision);

  rdpe_mul_eq (corr_mod, epsilon);
  rdpe_add_eq (root->drad, corr_mod);

  mpc_clear (newton_correction);
}

typedef struct {
  mps_context * ctx;
  mps_polynomial * p;
  mps_approximation * root;
  long int precision;
} __improve_root_data;

static void *
improve_root_wrapper (void * data_ptr)
{
  __improve_root_data *data = (__improve_root_data*) data_ptr;
  improve_root (data->ctx, data->p, data->root, data->precision);
  free (data);
  return NULL;
}


/**
 * @brief Improve all the approximations up to prec_out digits.
 *
 * For each approximation compute the value of sigma such that, given some
 * approximations \f$x_j\f$ of the roots, \f$r_j\f$ the values of the
 * inclusion radii and \f$d_i\f$ the number of correct digits:
 * \f[
 *   e_j < e_0 * \sigma^{2^j} \qquad \sigma=\frac{k}{k-1}=\frac{1}{1-t} \qquad k=\frac{1}{t}
 * \f]
 * and
 * \f[
 *   t = \min_j |z_i-z_j|-r_j
 * \f]
 * Then compute the number of digits needed for the j-th
 * iteration i.e., if \f$cond\f$ is the conditioning of the root: 
 * \f[
 *   d_j = \log(\frac{e_j}{|x|}) + cond
 * \f]
 * where 
 * \f[ 
 *   \log(\frac{e_j}{|x|}) = (f+g){2j} \qquad
 *   cond = \log(\frac{rad}{\epsilon})
 * \f] 
 * and
 * \f[ 
 *   cond \approx \lVert p \rVert (1+ \frac{|x_i|}{a_n \prod_{j \neq i} |x_i-x_j|}
 * \f]
 * and
 * \f[ 
 *   cond \approx \frac{r_i}{\epsilon |x_i|}
 * \f] 
 * for user-defined polynomials.
 *
 * <code>s->mpwp</code> denotes the number of bits of the current working
 * precision.
 *
 * @param ctx The mps_context associated with the computation.
 */
void
mps_improve (mps_context * ctx)
{
  int i;
  long int current_precision = 0L;
  int approximated_roots = 0;
  mps_polynomial * p = ctx->active_poly;
  rdpe_t * root_conditioning = NULL;
  ctx->operation = MPS_OPERATION_REFINEMENT;

  /* We need to be able to evaluate the Newton correction in a point
   * in order to perform the refinement. This is not necessary true
   * for custom polynomial types, so add a check in here */
  if (p->mnewton == NULL && p->density != MPS_DENSITY_USER)
    return;

  /* Set lastphase to mp */
  ctx->lastphase = mp_phase;

  /* Determine the conditioning of the roots */
  root_conditioning = evaluate_root_conditioning (ctx, p, ctx->root, ctx->n);

  /* We adopt the strategy of various iterations refinements on
   * the approximations by setting the precision of the input
   * polynomial in an increasing sequence. */

  /* We start by determining the minimum precision at which we can
   * extract some information. */
  current_precision = LONG_MAX;
  for (i = 0; i < ctx->n; i++)
    {
      if (ctx->root[i]->wp < current_precision)
        current_precision = ctx->root[i]->wp;

      if (MPS_ROOT_STATUS_IS_APPROXIMATED (ctx->root[i]->status))
        approximated_roots++;
    }

  /* Start by iterating on the roots that are not approximated, and
   * continue until we get all of them. */
  while (approximated_roots < ctx->n)
    {
      mps_polynomial_raise_data (ctx, p, current_precision);

      MPS_DEBUG (ctx, "Step of improvement");

      for (i = 0; i < ctx->n; i++)
        if (ctx->root[i]->status == MPS_ROOT_STATUS_ISOLATED)
          {
            /* Evaluate the necessary precision to iterate on this root. 
             * If the the current polynomial precision is enough, iterate on it. 
             * Otherwise, let it for the next round. */
            long int necessary_precision = get_approximated_bits (ctx->root[i]) + log2(ctx->n) + 
              rdpe_log(root_conditioning[i]) / LOG2;

            if (necessary_precision < current_precision)
            {
              __improve_root_data * data = mps_new (__improve_root_data);
              data->ctx = ctx;
              data->p = p;
              data->root = ctx->root[i];
              data->precision = current_precision;

              mps_thread_pool_assign (ctx, NULL, improve_root_wrapper, data);
            }
          }

      mps_thread_pool_wait (ctx, ctx->pool);

      for (i = 0; i < ctx->n; i++)
        if (get_approximated_bits (ctx->root[i]) >= ctx->output_config->prec)
          {
            ctx->root[i]->status = MPS_ROOT_STATUS_APPROXIMATED;
            approximated_roots++;

            if (ctx->debug_level & MPS_DEBUG_IMPROVEMENT)
              MPS_DEBUG (ctx, "Approximated roots = %d", approximated_roots);
          }

      /* Increase precision to reach the desired number of approximated roots */
      current_precision = 2 * current_precision;

      /* Check if we have gone too far with the precision, and we have gone over
       * the maximum precision allowed for this polynomial. */
      if (current_precision > p->prec && p->prec != 0)
        {
          ctx->over_max = true;
          goto cleanup;
        }

      /* Increase data prec max that will be useful to the end user to know
       * the precision needed to hold these approximations. */
      ctx->data_prec_max.value = current_precision;

      if (ctx->debug_level & MPS_DEBUG_IMPROVEMENT)
        MPS_DEBUG (ctx, "Increasing precision to %ld", current_precision);
    }

cleanup:    

  free (root_conditioning);
}
