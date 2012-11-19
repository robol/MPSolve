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

pthread_mutex_t output_mutex = PTHREAD_MUTEX_INITIALIZER;

  typedef struct {
    int i;
    mps_context * s;
    long int base_wp;
  } __mps_improve_data;

void * mps_improve_root2 (void*);
void * mps_improve_root  (void*);

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
      // mps_thread_pool_assign (s, NULL, mps_improve_root2, improve_data + i);
      mps_improve_root (improve_data + i);
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
  rdpe_t aroot, rtmp;
  long int wp = root->wp;

  /* Determine the number of steps necessary to have at least 
   * log10 (ctx->output_config->prec) correct digits. */
  mpc_rmod (aroot, root->mvalue);
  int correct_bits = rdpe_Esp (aroot) - rdpe_Esp (root->drad) - 1;
  int conditioning_bits = root->wp - correct_bits; 
  int max_steps = mps_intlog2 (ctx->output_config->prec / correct_bits);

  if (max_steps <= 0)
    max_steps = INT_MAX;

  mps_secular_iteration_data it_data;
  if (ctx->secular_equation)
    {
      it_data.local_ampc = ctx->secular_equation->ampc;
      it_data.local_bmpc = ctx->secular_equation->bmpc;
    }

  mpc_t nwtcorr;
  mps_monomial_poly * p = ctx->monomial_poly;
  mpc_init2 (nwtcorr, wp);

  /* We need to decide the initial precision needed for the iterations based on 
   * the conditioning of the root. If the conditioning is 2^k then we need w + k
   * bits of precision to have the right digits in the output. */
  rdpe_t conditioning;
  mps_evaluate_root_conditioning (ctx, root, conditioning);
  wp += conditioning_bits; 

  if (ctx->debug_level & MPS_DEBUG_IMPROVEMENT)
    MPS_DEBUG (ctx, "Starting to refine root %d", i);
  if (ctx->root_status[i] != MPS_ROOT_STATUS_ISOLATED || 
      ctx->root_status[i] == MPS_ROOT_STATUS_APPROXIMATED_IN_CLUSTER)
    {
      if (ctx->debug_level & MPS_DEBUG_IMPROVEMENT)
	MPS_DEBUG (ctx, "Not approximating root %d since it is already approximated", i);
      
      return NULL;
    }

  for (j = 0; max_steps ; j++)
    {
      /* mps_prepare_data (ctx, wp); */

      mpc_set_prec (nwtcorr, wp);
      mpc_set_prec (root->mvalue, wp);
      root->wp = wp;

      /* Update the value of data_prec_max if needed */
      MPS_LOCK (ctx->data_prec_max);
      if (ctx->data_prec_max.value < root->wp)
	ctx->data_prec_max.value = root->wp;
      MPS_UNLOCK (ctx->data_prec_max);

      /* if (ctx->mpwp < wp) */
	{
	  if (MPS_INPUT_CONFIG_IS_MONOMIAL (ctx->input_config))
	    {
	      mps_monomial_poly_raise_precision (ctx, ctx->monomial_poly, wp);
	    }
	  else if (MPS_INPUT_CONFIG_IS_SECULAR (ctx->input_config))
	    {
	      mps_secular_raise_coefficient_precision (ctx, wp);
	      it_data.local_ampc = ctx->secular_equation->ampc;
	      it_data.local_bmpc = ctx->secular_equation->bmpc;
	    }

	  /* ctx->mpwp = wp; */
	}

      /* Save the old radius before the iteration */
      rdpe_set (rtmp, root->drad);
      rdpe_div_eq (rtmp, aroot);
      rdpe_sqr_eq (rtmp);
      rdpe_mul_eq (rtmp, aroot);
      
      if (MPS_INPUT_CONFIG_IS_MONOMIAL (ctx->input_config))
	{
	  mps_mnewton (ctx, ctx->n, root, nwtcorr, p->mfpc, p->mfppc, p->dap, p->spar,
		       mps_thread_get_id (ctx, ctx->pool), false);
	}
      else
	{
	  (*ctx->mnewton_usr) (ctx, root, nwtcorr, 
			       MPS_INPUT_CONFIG_IS_SECULAR (ctx->input_config) ? &it_data : NULL, false);
	}

      mpc_set_prec (root->mvalue, 2 * wp);
      mpc_sub_eq (root->mvalue, nwtcorr);

      MPS_DEBUG_MPC (ctx, 15, nwtcorr, "nwtcorr");

      if (rdpe_Esp (aroot) - rdpe_Esp (root->drad) - 1 > correct_bits)
	correct_bits = rdpe_Esp (aroot) - rdpe_Esp (root->drad) - 1;
       
      /* Double the number of correct bits */
      correct_bits = 2 * correct_bits - 1;
      MPS_DEBUG (ctx, "Correct bits for root %d = %d", i, correct_bits); 

      /* Set a proper radius to the approximations */
      rdpe_set_2dl (root->drad, 2.0, - correct_bits);   
      rdpe_mul_eq (root->drad, aroot);  
      
      if (rdpe_lt (rtmp, root->drad))  
	rdpe_set (root->drad, rtmp);  

      MPS_DEBUG_MPC (ctx, 45, root->mvalue, "Approximation");
      MPS_DEBUG_RDPE (ctx, root->drad, "Radius");

      rdpe_div (rtmp, root->drad, aroot);
      if (-rdpe_Esp (rtmp) > ctx->output_config->prec)
	break;

      /* if (correct_bits > ctx->output_config->prec) */
      /* 	break; */

      /* Double the current precision */
      wp *= 2;
    }

  mps_approximation_free (ctx, ctx->root[i]);
  ctx->root[i] = root;
  ctx->root_status[i] = MPS_ROOT_STATUS_APPROXIMATED;

  mpc_clear (nwtcorr);

  return NULL;
}

void *
mps_improve_root (void * data_ptr)
{
  __mps_improve_data * data = (__mps_improve_data*) data_ptr;
  int i = data->i;
  mps_context * s = data->s;
  int j, k, m;
  int mpwp = s->mpwp;
  long mpnb_in, mpnb_out;
  mpc_t mtmp;
  mpc_t nwtcorr;
  cdpe_t ctmp;
  rdpe_t tmp, t, st, sigma, newrad, oldrad, abroot, mp_epsilon;
  double f, g, cnd;
  mps_monomial_poly *p = s->monomial_poly;

  mps_secular_iteration_data it_data;
  if (s->secular_equation)
    {
      it_data.local_ampc = s->secular_equation->ampc;
      it_data.local_bmpc = s->secular_equation->bmpc;
    }

   if (s->debug_level & MPS_DEBUG_IMPROVEMENT) 
     { 
       MPS_DEBUG (s, "Refining the roots"); 
     }

  /* == 1 ==
   * compute the number mpnb_in of bits
   * corresponding to the given input precision.
   * Set mpnb_in=0 if the input precision is infinite (prec_in=0) */
  if (s->input_config->prec == 0)
    mpnb_in = 0;
  else
    mpnb_in = (long) (s->input_config->prec * LOG2_10 + log (4.0 * s->n) / LOG2);
  mpnb_out = (long) (s->output_config->prec * LOG2_10);

  /* == 2  ==
   * compute the coefficients of the polynomial as mpc_t with mpnb_in bits
   * only if the polynomial is not assigned as a straight line program and
   * the input precision is not infinite. */
  if (mpnb_in != 0)
    mps_mp_set_prec (s, mpnb_in);

  /* mpc_init2(mtmp, mpwp);  *//* puo' essere settato a precisione minima */
  mpc_init2 (mtmp, mpnb_out * 2);      /* puo' essere settato a precisione minima */
  mpc_init2 (nwtcorr, mpnb_out * 2);

  if (s->input_config->prec != 0 && !MPS_INPUT_CONFIG_IS_USER (s->input_config) && mpnb_in < s->mpwp)
    mps_prepare_data (s, mpnb_in);
  else
    {
      if (mpnb_out * 2 > s->mpwp)
	{
	  mps_mp_set_prec (s, mpnb_out * 2);
	  mps_prepare_data (s, mpnb_out * 2);
	}
    }


  /* == 3 ==
   * scan the approximations to apply Newton's iterations */
  /* for (i = 0; i < s->n; i++) */
    {
      if (s->debug_level & MPS_DEBUG_IMPROVEMENT)
        MPS_DEBUG (s, "Starting to refine root %d", i);
      if (s->root_status[i] != MPS_ROOT_STATUS_ISOLATED || 
	  s->root_status[i] == MPS_ROOT_STATUS_APPROXIMATED_IN_CLUSTER)
        {
	  if (s->debug_level & MPS_DEBUG_IMPROVEMENT)
	    MPS_DEBUG (s, "Not approximating root %d since it is already approximated", i);

          goto improve_clear;             /* Do not refine approximated roots */
        }

      /*  == 3.1 ==
       * for data_type[0]='d' compute  t=Min_j |root(i)-root(j)|-rad(j)-rad(i)
       * otherwise set t=5*n*rad[i] since the root is Newton-isolated.
       * This allows us to remove an O(n^2) complexity  */

      if (MPS_INPUT_CONFIG_IS_SPARSE (s->input_config))
        rdpe_mul_d (t, s->root[i]->drad, 5.0 * s->n);
      else
        {
          k = i + 1;
          if (i == s->n - 1)
            k = 0;
          mpc_sub (mtmp, s->root[k]->mvalue, s->root[i]->mvalue);
          mpc_get_cdpe (ctmp, mtmp);
          cdpe_mod (t, ctmp);
          rdpe_sub_eq (t, s->root[k]->drad);
          rdpe_sub_eq (t, s->root[i]->drad);
          for (j = 0; j < s->n; j++)
            if (j != i)
              {
                mpc_sub (mtmp, s->root[j]->mvalue, s->root[i]->mvalue);
                mpc_get_cdpe (ctmp, mtmp);
                cdpe_mod (tmp, ctmp);
                rdpe_sub_eq (tmp, s->root[i]->drad);
                rdpe_sub_eq (tmp, s->root[j]->drad);
                if (rdpe_gt (t, tmp))
                  rdpe_set (t, tmp);
              }
        }

      /*  == 3.2 ==
       * compute an  estimate of the condition number in terms of bits
       * as log_2(rad/(4*n*epsilon*|x|))       */

      rdpe_mul_d (tmp, s->root[i]->drad, 4.0 * s->n);
      mpc_get_cdpe (ctmp, s->root[i]->mvalue);
      cdpe_mod (abroot, ctmp);
      rdpe_div (tmp, tmp, abroot);

      cnd = s->root[i]->wp + rdpe_log (tmp) / LOG2 + 1;

      /* then evaluate the number of bits g,f */
      rdpe_div (t, s->root[i]->drad, t);
      rdpe_mul_eq_d (t, (double) s->n - 1);
      rdpe_sub (st, rdpe_one, t);
      rdpe_div (sigma, t, st);

      /* Workaround added by me to solve nan problems. Leo. */
      rdpe_set_2dl (mp_epsilon, 2.0, - mpwp);
      rdpe_add_eq (sigma, mp_epsilon); 

      g = -rdpe_log (sigma) / LOG2;
      rdpe_set (tmp, abroot);
      rdpe_mul_eq (tmp, sigma);
      rdpe_div (tmp, s->root[i]->drad, tmp);
      f = -rdpe_log (tmp) / LOG2;

      /* evaluate the upper bound m to the number of iterations
       * needed to reach the desired precision */
      m = (int) (log ((mpnb_out - f) / g) / LOG2) + 1;

      MPS_DEBUG (s, "A maximum of %d iterations will be performed to improve root %d", m, i);

      /*  == 4 ==      Start Newton */
      rdpe_set (oldrad, s->root[i]->drad);
      for (j = 1; j <= m; j++)
        {
          if (s->debug_level & MPS_DEBUG_IMPROVEMENT)
            MPS_DEBUG (s, "Iteration %d of the improvement of root %d", j, i);
          g *= 2;

	  /* { */
	  /*   rdpe_t rtmp; */
	  /*   mpc_rmod (rtmp, s->root[i]->mvalue); */
	  /*   int correct_digits = (-rdpe_log (s->root[i]->drad) - rdpe_log (rtmp)) / LOG2_10; */
	  /*   MPS_DEBUG_RDPE (s, s->root[i]->drad, "s->drad[%d]", i); */
	  /*   MPS_DEBUG_MPC (s, correct_digits, s->root[i]->mvalue, "mroot_%d", i); */
	  /* } */

	  /* Round it to 64 integers */
	  mpwp = (long) (f + g + cnd);

          if (mpwp > mpnb_in && mpnb_in != 0)
	    {
	      /* Lower the precision so it won't go over mpnb_in
	       * that would clearly get us to an error, for over estimating
	       * the precision of the input coefficients. */
	      mpwp = mpnb_in - 63;
	      break;
	    }

          mps_mp_set_prec (s, mpwp);

          mpc_clear (nwtcorr);
          mpc_init2 (nwtcorr, mpwp);

	  /* If using the standard MPSolve algorithm then use the old
	   * mps_prepare_data routine, otherwise use the one that
	   * raises the precision of the coefficients */
	  if (mpwp > s->mpwp)
	    mps_prepare_data (s, mpwp);

          if (MPS_INPUT_CONFIG_IS_MONOMIAL (s->input_config))
            {
              mps_mnewton (s, s->n, s->root[i], nwtcorr, p->mfpc, p->mfppc, p->dap, p->spar,
                           0, false);
            }
          else if (s->mnewton_usr != NULL)
            {
              (*s->mnewton_usr) (s, s->root[i], nwtcorr, 
				 MPS_INPUT_CONFIG_IS_SECULAR (s->input_config) ? &it_data : NULL, false);
            }
          else
            {
              mps_mnewton_usr (s, s->root[i], nwtcorr);
            }
          mpc_sub_eq (s->root[i]->mvalue, nwtcorr);

          /* correct radius, since the computed one is referred to the previous
           * approximation. Due to the quadratic convergence the new approximation
           * the radius is bounded by 2^(-g-f+1) */
          rdpe_set_2dl (newrad, 4.0, (long) (-g - f + 1));
          rdpe_set (tmp, abroot);
          rdpe_mul_eq (newrad, tmp);
          rdpe_mul_eq (tmp, s->eps_out);

	  if (rdpe_eq (s->root[i]->drad, rdpe_zero)) 
	    rdpe_set (s->root[i]->drad, newrad); 

	  if (rdpe_lt (newrad, s->root[i]->drad))    
	    rdpe_set (s->root[i]->drad, newrad);    
	   
	  mpc_rmod (tmp, s->root[i]->mvalue);
	  rdpe_mul_eq (tmp, s->eps_out);
	  rdpe_mul_eq_d (tmp, 4.0);
	  rdpe_add_eq (s->root[i]->drad, tmp);
	  
	  if (s->debug_level & MPS_DEBUG_IMPROVEMENT)
	    MPS_DEBUG_RDPE (s, s->root[i]->drad, "Radius of root %d at iteration %d", i, j);
	   
	  /* Check if the radius that we have obtained until now is good, and if
	   * we have passed the maximum allowed precision. */
          if (rdpe_lt (s->root[i]->drad, tmp) || 
	      (mpnb_in != 0 && mpwp >= mpnb_in))
	    {
	      if (mpwp >= mpnb_in && mpnb_in != 0)
		s->over_max = true;

	      if (s->debug_level & MPS_DEBUG_IMPROVEMENT)
		{
		  if (mpwp >= mpnb_in && mpnb_in != 0)
		    {
		      MPS_DEBUG (s, "Stopping newton iterations on root %d because we have reached input precision", i);
		    }
		  else
		    {
		      MPS_DEBUG (s, "Stopping newton iterations on root %d because radius is small enough", i);
		    }
		}
	      break;
	    }
        }
    }

    improve_clear:

      mpc_clear (nwtcorr);
      mpc_clear (mtmp);

      return NULL;
}
