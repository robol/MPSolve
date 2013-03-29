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
    mps_approximation * starting_approximation;
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

  MPS_DEBUG_WITH_INFO (s, "Lowering the number of threads to 1");
  mps_thread_pool_set_concurrency_limit (s, NULL, 1);

  for (i = 0; i < s->n; i++)
    {
      improve_data[i].s = s;
      improve_data[i].base_wp = base_wp;
      improve_data[i].i = i;
      improve_data[i].starting_approximation = mps_approximation_copy (s, s->root[i]);
    }

  for (i = 0; i < s->n; i++)
    {
      if (s->root[i]->status != MPS_ROOT_STATUS_ISOLATED || 
	       s->root[i]->status == MPS_ROOT_STATUS_APPROXIMATED_IN_CLUSTER)
      	{
      	  if (s->debug_level & MPS_DEBUG_IMPROVEMENT)
      	    MPS_DEBUG (s, "Not approximating root %i since it is already approximated", i);	 
      	}
      else
        mps_thread_pool_assign (s, NULL, mps_improve_root2, improve_data + i);
    }

  mps_thread_pool_wait (s, s->pool);

  for (i = 0; i < s->n; i++)
    {
      mps_approximation_free (s, s->root[i]);
      s->root[i] = improve_data[i].starting_approximation;
      s->root[i]->status = MPS_ROOT_STATUS_APPROXIMATED;
    }

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
  mps_context * ctx = data->s;
  mps_approximation * root = data->starting_approximation;
  mps_polynomial * p = ctx->active_poly;
  long int wp = root->wp;

  mpc_t nwtcorr;

  mpc_init2 (nwtcorr, 64);

  /* Get the number of correct digits that you have obtained until 
   * now. */
  rdpe_t root_mod;

  mpc_rmod (root_mod, root->mvalue);
  int correct_bits = MAX (0, rdpe_Esp (root_mod) - rdpe_Esp (root->drad) - 1);

  rdpe_t conditioning;
  rdpe_div (conditioning, root->drad, root_mod);
  rdpe_div_eq (conditioning, root->drad);

  int conditioning_bits = MAX (0, rdpe_log (conditioning)) / LOG2 + 
    MAX (mpc_get_prec (root->mvalue), root->wp);

  MPS_DEBUG_MPC (ctx, 15, root->mvalue, "Approximating root ");
  MPS_DEBUG (ctx, "Correct bits = %d", correct_bits);
  MPS_DEBUG (ctx, "Conditioning bits = %d", conditioning_bits);

  while (correct_bits <= ctx->output_config->prec || rdpe_eq_zero (root_mod))
  {
    rdpe_t nwtcorr_mod;

    wp = 2 * correct_bits + conditioning_bits;

    if (wp >= p->prec && p->prec != 0)
    {
      MPS_DEBUG (ctx, 
        "Reached maximum allowed precision due to limited input precision. Aborting improvement");
      MPS_DEBUG(ctx, "Precision = %ld", p->prec)
      ctx->over_max = true;
      break;
    }

    mpc_set_prec (nwtcorr, wp);
    mpc_set_prec (root->mvalue, wp);

    mps_polynomial_raise_data (ctx, p, wp);

    root->wp = wp;

    MPS_LOCK (ctx->data_prec_max);
    if (ctx->data_prec_max.value < wp)
      ctx->data_prec_max.value = wp;
    MPS_UNLOCK (ctx->data_prec_max);

    mps_polynomial_mnewton (ctx, p, root, nwtcorr);

    mpc_sub_eq (root->mvalue, nwtcorr);
    mpc_rmod (root_mod, root->mvalue);
    mpc_rmod (nwtcorr_mod, nwtcorr);

    rdpe_add_eq (root->drad, nwtcorr_mod);

    correct_bits = MAX (rdpe_Esp (root_mod) - rdpe_Esp (root->drad), 0);

    if (ctx->debug_level & MPS_DEBUG_IMPROVEMENT)
      MPS_DEBUG (ctx, "    Correct bits = %d", correct_bits);
  }

  if (correct_bits >= ctx->output_config->prec)
    root->status = MPS_ROOT_STATUS_APPROXIMATED;

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

   mps_polynomial * p = s->active_poly;
   mps_approximation * root = data->starting_approximation;

    if (s->debug_level & MPS_DEBUG_IMPROVEMENT)  
      {  
        MPS_DEBUG (s, "Refining the roots");  
      } 

   /* == 1 == 
    * compute the number mpnb_in of bits 
    * corresponding to the given input precision. 
    * Set mpnb_in=0 if the input precision is infinite (prec_in=0) */ 
   if (p->prec == 0) 
     mpnb_in = 0; 
   else 
     mpnb_in = (long) (p->prec * LOG2_10 + log (4.0 * s->n) / LOG2); 
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

   if (p->prec != 0 && mpnb_in < s->mpwp) 
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
       if (s->root[i]->status != MPS_ROOT_STATUS_ISOLATED ||  
 	  s->root[i]->status == MPS_ROOT_STATUS_APPROXIMATED_IN_CLUSTER) 
         { 
 	  if (s->debug_level & MPS_DEBUG_IMPROVEMENT) 
 	    MPS_DEBUG (s, "Not approximating root %d since it is already approximated", i); 

           goto improve_clear;             /* Do not refine approximated roots */ 
         } 

       /*  == 3.1 == 
        * for data_type[0]='d' compute  t=Min_j |root(i)-root(j)|-rad(j)-rad(i) 
        * otherwise set t=5*n*rad[i] since the root is Newton-isolated. 
        * This allows us to remove an O(n^2) complexity  */ 

       if (MPS_DENSITY_IS_SPARSE (p->density)) 
         rdpe_mul_d (t, root->drad, 5.0 * s->n); 
       else 
         { 
           k = i + 1; 
           if (i == s->n - 1) 
             k = 0; 
           mpc_sub (mtmp, s->root[k]->mvalue, root->mvalue); 
           mpc_get_cdpe (ctmp, mtmp); 
           cdpe_mod (t, ctmp); 
           rdpe_sub_eq (t, s->root[k]->drad); 
           rdpe_sub_eq (t, root->drad); 
           for (j = 0; j < s->n; j++) 
             if (j != i) 
               { 
                 mpc_sub (mtmp, s->root[j]->mvalue, root->mvalue); 
                 mpc_get_cdpe (ctmp, mtmp); 
                 cdpe_mod (tmp, ctmp); 
                 rdpe_sub_eq (tmp, root->drad); 
                 rdpe_sub_eq (tmp, s->root[j]->drad); 
                 if (rdpe_gt (t, tmp)) 
                   rdpe_set (t, tmp); 
               } 
         } 

       /*  == 3.2 == 
        * compute an  estimate of the condition number in terms of bits 
        * as log_2(rad/(4*n*epsilon*|x|))       */ 

       rdpe_mul_d (tmp, root->drad, 4.0 * s->n); 
       mpc_get_cdpe (ctmp, root->mvalue); 
       cdpe_mod (abroot, ctmp); 
       rdpe_div (tmp, tmp, abroot); 

       cnd = root->wp + rdpe_log (tmp) / LOG2 + 1; 

       /* then evaluate the number of bits g,f */ 
       rdpe_div (t, root->drad, t); 
       rdpe_mul_eq_d (t, (double) s->n - 1); 
       rdpe_sub (st, rdpe_one, t); 
       rdpe_div (sigma, t, st); 

       /* Workaround added by me to solve nan problems. Leo. */ 
       rdpe_set_2dl (mp_epsilon, 2.0, - mpwp); 
       rdpe_add_eq (sigma, mp_epsilon);  

       g = -rdpe_log (sigma) / LOG2; 
       rdpe_set (tmp, abroot); 
       rdpe_mul_eq (tmp, sigma); 
       rdpe_div (tmp, root->drad, tmp); 
       f = -rdpe_log (tmp) / LOG2; 

       /* evaluate the upper bound m to the number of iterations 
        * needed to reach the desired precision */ 
       m = (int) (log ((mpnb_out - f) / g) / LOG2) + 1; 

       MPS_DEBUG (s, "A maximum of %d iterations will be performed to improve root %d", m, i); 

       /*  == 4 ==      Start Newton */ 
       rdpe_set (oldrad, root->drad); 
       for (j = 1; j <= m; j++) 
         { 
           if (s->debug_level & MPS_DEBUG_IMPROVEMENT) 
             MPS_DEBUG (s, "Iteration %d of the improvement of root %d", j, i); 
           g *= 2; 

 	  /* { */ 
 	  /*   rdpe_t rtmp; */ 
 	  /*   mpc_rmod (rtmp, root->mvalue); */ 
 	  /*   int correct_digits = (-rdpe_log (root->drad) - rdpe_log (rtmp)) / LOG2_10; */ 
 	  /*   MPS_DEBUG_RDPE (s, root->drad, "s->drad[%d]", i); */ 
 	  /*   MPS_DEBUG_MPC (s, correct_digits, root->mvalue, "mroot_%d", i); */ 
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

      mpc_set_prec (nwtcorr, mpwp);
      mpc_set_prec (root->mvalue, mpwp);
      mps_polynomial_mnewton (s, p, root, nwtcorr);
      mpc_sub_eq (root->mvalue, nwtcorr); 

           /* correct radius, since the computed one is referred to the previous 
            * approximation. Due to the quadratic convergence the new approximation 
            * the radius is bounded by 2^(-g-f+1) */ 
           rdpe_set_2dl (newrad, 4.0, (long) (-g - f + 1)); 
           rdpe_set (tmp, abroot); 
           rdpe_mul_eq (newrad, tmp); 
           rdpe_mul_eq (tmp, s->eps_out); 

 	  if (rdpe_eq (root->drad, rdpe_zero))  
 	    rdpe_set (root->drad, newrad);  

 	  if (rdpe_lt (newrad, root->drad))     
 	    rdpe_set (root->drad, newrad);     
	   
 	  mpc_rmod (tmp, root->mvalue); 
 	  rdpe_mul_eq (tmp, s->eps_out); 
 	  rdpe_mul_eq_d (tmp, 4.0); 
 	  rdpe_add_eq (root->drad, tmp); 
	  
 	  if (s->debug_level & MPS_DEBUG_IMPROVEMENT) 
 	    MPS_DEBUG_RDPE (s, root->drad, "Radius of root %d at iteration %d", i, j); 
	   
 	  /* Check if the radius that we have obtained until now is good, and if 
 	   * we have passed the maximum allowed precision. */ 
           if (rdpe_lt (root->drad, tmp) ||  
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
