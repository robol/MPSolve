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

pthread_mutex_t output_mutex = PTHREAD_MUTEX_INITIALIZER;

  typedef struct {
    int i;
    mps_context * s;
    long int base_wp;
  } __mps_improve_data;

void * mps_improve_root2 (void*);

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
 * @param s The mps_context associated with the computation.
 */
void
mps_improve (mps_context * s)
{
  int i;
  clock_t *my_timer = mps_start_timer ();

  s->operation = MPS_OPERATION_REFINEMENT;

  /* Set lastphase to mp */
  s->lastphase = mp_phase;

  __mps_improve_data *improve_data = mps_newv (__mps_improve_data, s->n);

  long int base_wp = s->mpwp;

  for (i = 0; i < s->n; i++)
    {
      improve_data[i].i = i;
      improve_data[i].s = s;
      improve_data[i].base_wp = base_wp;
      mps_thread_pool_assign (s, NULL, mps_improve_root2, improve_data + i);
      // mps_improve_root (s, i);
    }

  mps_thread_pool_wait (s, s->pool);
  free (improve_data);

  long improve_time = mps_stop_timer (my_timer);
  if (s->debug_level & MPS_DEBUG_TIMINGS)
    MPS_DEBUG (s, "Improvement of roots took %lu ms", improve_time);
}

void
mps_evaluate_root_conditioning (mps_context *ctx, mps_approximation *root, rdpe_t conditioning)
{
  rdpe_set (conditioning, rdpe_one);
}

void *
mps_improve_root2 (void * data_ptr)
{
  __mps_improve_data * data = (__mps_improve_data*) data_ptr;
  int i = data->i;

  int j;
  mps_context *ctx = data->s;

  mps_approximation * root = mps_approximation_copy (ctx, ctx->root[i]);
  rdpe_t aroot;
  long int wp = root->wp;

  /* Determine the number of steps necessary to have at least 
   * log10 (ctx->output_config->prec) correct digits. */
  mpc_rmod (aroot, root->mvalue);
  int correct_bits = rdpe_Esp (aroot) - rdpe_Esp (root->drad) - 1;
  int max_steps = mps_intlog2 (ctx->output_config->prec / correct_bits);

  mps_secular_iteration_data it_data;
  if (MPS_INPUT_CONFIG_IS_SECULAR (ctx->input_config))
    {
      mps_secular_equation *sec = MPS_SECULAR_EQUATION (ctx->active_poly);
      it_data.local_ampc = sec->ampc;
      it_data.local_bmpc = sec->bmpc;
    }

  mpc_t nwtcorr;
  mps_polynomial * p = ctx->active_poly;
  mpc_init2 (nwtcorr, wp);

  /* We need to decide the initial precision needed for the iterations based on 
   * the conditioning of the root. If the conditioning is 2^k then we need w + k
   * bits of precision to have the right digits in the output. */
  rdpe_t conditioning;
  mps_evaluate_root_conditioning (ctx, root, conditioning);
  wp += rdpe_Esp (conditioning);

  if (ctx->debug_level & MPS_DEBUG_IMPROVEMENT)
    MPS_DEBUG (ctx, "Starting to refine root %d", i);
  if (ctx->root_status[i] != MPS_ROOT_STATUS_ISOLATED || 
      ctx->root_status[i] == MPS_ROOT_STATUS_APPROXIMATED_IN_CLUSTER)
    {
      if (ctx->debug_level & MPS_DEBUG_IMPROVEMENT)
	MPS_DEBUG (ctx, "Not approximating root %d since it is already approximated", i);
      
      return NULL;
    }

  for (j = 0; j < max_steps; j++)
    {
      /* mps_prepare_data (ctx, wp); */

      mpc_set_prec (nwtcorr, wp);
      mpc_set_prec (root->mvalue, wp);
      root->wp = wp;

      if (ctx->mpwp < wp)
	{
	  mps_polynomial_raise_data (ctx, p, wp);

	  if (MPS_INPUT_CONFIG_IS_SECULAR (ctx->input_config))
	    {
	      mps_secular_equation *sec = MPS_SECULAR_EQUATION (p);
	      it_data.local_ampc = sec->ampc;
	      it_data.local_bmpc = sec->bmpc;
	    }

	  ctx->mpwp = wp;
	}
      
      if (MPS_INPUT_CONFIG_IS_MONOMIAL (ctx->input_config))
	{
	  mps_monomial_poly *mp = MPS_MONOMIAL_POLY (p);
	  mps_mnewton (ctx, ctx->n, root, nwtcorr, mp->mfpc, mp->mfppc, mp->dap, mp->spar,
		       mps_thread_get_id (ctx, ctx->pool), false);
	}
      else
	{
	  (*ctx->mnewton_usr) (ctx, root, nwtcorr, 
			       MPS_INPUT_CONFIG_IS_SECULAR (ctx->input_config) ? &it_data : NULL, false);
	}

      mpc_set_prec (root->mvalue, 2 * wp);

      mpc_sub_eq (root->mvalue, nwtcorr);   
       
      /* Double the number of correct bits */
      correct_bits = 2 * correct_bits - 1;

      /* Set a proper radius to the approximations */
      rdpe_set_2dl (root->drad, 2.0, - correct_bits); 
      rdpe_mul_eq (root->drad, aroot);

      /* Double the current precision */
      wp *= 2;
    }

  mps_approximation_free (ctx, ctx->root[i]);
  ctx->root[i] = root;

  mpc_clear (nwtcorr);

  return NULL;
}

