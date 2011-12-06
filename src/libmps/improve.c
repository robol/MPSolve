/***********************************************************
**       Multiprecision Polynomial Solver (MPSolve)       **
**                 Version 2.2, May 2001                  **
**                                                        **
**                      Written by                        **
**       Dario Andrea Bini and Giuseppe Fiorentino        **
**       (bini@dm.unipi.it)  (fiorent@dm.unipi.it)        **
**                                                        **
** (C) 2001, Dipartimento di Matematica, FRISCO LTR 21024 **
***********************************************************/

#include <mps/core.h>
#include <mps/secular.h>
#include <mps/debug.h>
#include <math.h>


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
 * @param s The mps_status associated with the computation.
 */
void
mps_improve (mps_status * s)
{
  int i, j, k, m;
  long mpnb_in, mpnb_out;
  mpc_t mtmp;
  mpc_t nwtcorr;
  cdpe_t ctmp;
  rdpe_t tmp, t, st, sigma, newrad, oldrad, abroot;
  double f, g, cnd;
  mps_boolean again;
  mps_monomial_poly *p = s->monomial_poly;
  clock_t *my_timer = mps_start_timer ();

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

  /* mpc_init2(mtmp, mpwp); *//* puo' essere settato a precisione minima */
  mpc_init2 (mtmp, mpnb_out * 2);      /* puo' essere settato a precisione minima */
  mpc_init2 (nwtcorr, mpnb_out * 2);

  if (s->input_config->prec != 0 && s->data_type[0] != 'u')
    mps_prepare_data (s, mpnb_in);
  else
    {
      mps_mp_set_prec (s, mpnb_out * 2);
      mps_prepare_data (s, mpnb_out * 2);
      if (MPS_INPUT_CONFIG_IS_SECULAR (s->input_config))
	mps_secular_raise_coefficient_precision (s, mpnb_out * 2);
    }


  /* == 3 ==
   * scan the approximations to apply Newton's iterations */
  for (i = 0; i < s->n; i++)
    {
      if (s->debug_level & MPS_DEBUG_IMPROVEMENT)
        MPS_DEBUG (s, "Starting to refine root %d", i);
      if (s->status[i][0] != 'i' || s->status[i][2] == 'o')
        {
	  if (s->debug_level & MPS_DEBUG_IMPROVEMENT)
	    MPS_DEBUG (s, "Not approximating root %d since it is already approximated", i);

          continue;             /* Do not refine approximated roots */
        }

      /*  == 3.1 ==
       * for data_type[0]='d' compute  t=Min_j |root(i)-root(j)|-rad(j)-rad(i)
       * otherwise set t=5*n*rad[i] since the root is Newton-isolated.
       * This allows us to remove an O(n^2) complexity  */

      if (s->data_type[0] == 's')
        rdpe_mul_d (t, s->drad[i], 5.0 * s->n);
      else
        {
          k = i + 1;
          if (i == s->n - 1)
            k = 0;
          mpc_sub (mtmp, s->mroot[k], s->mroot[i]);
          mpc_get_cdpe (ctmp, mtmp);
          cdpe_mod (t, ctmp);
          rdpe_sub_eq (t, s->drad[k]);
          rdpe_sub_eq (t, s->drad[i]);
          for (j = 0; j < s->n; j++)
            if (j != i)
              {
                mpc_sub (mtmp, s->mroot[j], s->mroot[i]);
                mpc_get_cdpe (ctmp, mtmp);
                cdpe_mod (tmp, ctmp);
                rdpe_sub_eq (tmp, s->drad[i]);
                rdpe_sub_eq (tmp, s->drad[j]);
                if (rdpe_gt (t, tmp))
                  rdpe_set (t, tmp);
              }
        }

      /*  == 3.2 ==
       * compute an  estimate of the condition number in terms of bits
       * as log_2(rad/(4*n*epsilon*|x|))       */

      rdpe_mul_d (tmp, s->drad[i], 4.0 * s->n);
      mpc_get_cdpe (ctmp, s->mroot[i]);
      cdpe_mod (abroot, ctmp);
      rdpe_div (tmp, tmp, abroot);
      cnd = s->rootwp[i] + rdpe_log (tmp) / LOG2 + 1;

      /* then evaluate the number of bits g,f */
      rdpe_div (t, s->drad[i], t);
      rdpe_mul_eq_d (t, (double) s->n - 1);
      rdpe_sub (st, rdpe_one, t);
      rdpe_div (sigma, t, st);
      g = -rdpe_log (sigma) / LOG2;
      rdpe_set (tmp, abroot);
      rdpe_mul_eq (tmp, sigma);
      rdpe_div (tmp, s->drad[i], tmp);
      f = -rdpe_log (tmp) / LOG2;

      /* evaluate the upper bound m to the number of iterations
       * needed to reach the desired precision */
      m = (int) (log ((mpnb_out - f) / g) / LOG2) + 1;

      /*  == 4 ==      Start Newton */
      rdpe_set (oldrad, s->drad[i]);
      for (j = 1; j <= m; j++)
        {
          if (s->debug_level & MPS_DEBUG_IMPROVEMENT)
            MPS_DEBUG (s, "Iteration %d of the improvement of root %d", j, i);
          g *= 2;
          s->mpwp = (long) (f + g + cnd);
          if (s->mpwp >= mpnb_in && mpnb_in != 0)
            s->mpwp = mpnb_in;

          mpc_clear (nwtcorr);
          mpc_init2 (nwtcorr, s->mpwp);

          mps_mp_set_prec (s, s->mpwp);

	  /* If using the standard MPSolve algorithm then use the old
	   * mps_prepare_data routine, otherwise use the one that
	   * raises the precision of the coefficients */
	  mps_prepare_data (s, s->mpwp);
          if (MPS_INPUT_CONFIG_IS_SECULAR (s->input_config))
	      mps_secular_raise_coefficient_precision (s, s->mpwp);

          if (s->data_type[0] != 'u')
            {
              mps_mnewton (s, s->n, s->mroot[i], s->drad[i],
                           nwtcorr, p->mfpc, p->mfppc, p->dap, p->spar,
                           &again, 0, false);
            }
          else if (s->mnewton_usr != NULL)
            {
              (*s->mnewton_usr) (s, s->mroot[i], s->drad[i], nwtcorr, &again, NULL, false);
            }
          else
            {
              mps_mnewton_usr (s, s->mroot[i], s->drad[i], nwtcorr, &again);
            }
          mpc_sub_eq (s->mroot[i], nwtcorr);

          /* correct radius, since the computed one is referred to the previous
           * approximation. Due to the quadratic convergence the new approximation
           * the radius is bounded by 2^(-g-f+1) */
          rdpe_set_2dl (newrad, 1.0, (long) (-g - f + 1));
          rdpe_set (tmp, abroot);
          rdpe_mul_eq (newrad, tmp);
          rdpe_mul_eq (tmp, s->eps_out);

          if (rdpe_eq (s->drad[i], rdpe_zero))
            rdpe_set (s->drad[i], newrad);
          if (rdpe_lt (newrad, s->drad[i]))
            rdpe_set (s->drad[i], newrad);

          if (rdpe_lt (s->drad[i], tmp) || s->mpwp == mpnb_in)
            break;              /* loop1 */
        }

      /* update the record working precision for root i */
      s->rootwp[i] = mpc_get_prec (s->mroot[i]);
    }

  mpc_clear (nwtcorr);
  mpc_clear (mtmp);

  long improve_time = mps_stop_timer (my_timer);
  if (s->debug_level & MPS_DEBUG_TIMINGS)
    MPS_DEBUG (s, "Improvement of roots took %lu ms", improve_time);
}
