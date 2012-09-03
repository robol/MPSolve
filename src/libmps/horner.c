/************************************************************
 **                                                        **
 **             __  __ ___  ___      _                     **
 **            |  \/  | _ \/ __| ___| |_ _____             **
 **            | |\/| |  _/\__ \/ _ \ \ V / -_)            **
 **            |_|  |_|_|  |___/\___/_|\_/\___|            **
 **                                                        **
 **       Multiprecision Polynomial Solver (MPSolve)       **
 **                 Version 2.9, April 2011                **
 **                                                        **
 **                      Written by                        **
 **                                                        **
 **     Dario Andrea Bini       <bini@dm.unipi.it>         **
 **     Giuseppe Fiorentino     <fiorent@dm.unipi.it>      **
 **     Leonardo Robol          <robol@mail.dm.unipi.it>   **
 **                                                        **
 **           (C) 2011, Dipartimento di Matematica         **
 ***********************************************************/

#include <mps/mps.h>

void
mps_mhorner_sparse (mps_status * s, mps_monomial_poly * p, mpc_t x, mpc_t value);

/**
 * @brief Compute the value of the polynomial <code>p</code> in the point <code>x</code>
 * and save it in <code>value</code>. If you need a bound to the relative error, try
 * mps_mhorner_with_error().
 *
 * @param s The <code>mps_status</code> of the computation.
 * @param p The <code>monomial_poly</code> to evaluate.
 * @param x The point where the polynomial will be evaluated.
 * @param value The multiprecision complex variable where the result will be stored.
 */
void 
mps_mhorner (mps_status * s, mps_monomial_poly * p, mpc_t x, mpc_t value)
{
  int j;

  if (MPS_INPUT_CONFIG_IS_SPARSE (s->input_config))
    {
      mps_mhorner_sparse (s, p, x, value);
    }
  else  
    { 
      mps_with_lock (p->mfpc_mutex[p->n],
		     mpc_set (value, p->mfpc[p->n]);
		     );

      for (j = p->n - 1; j >= 0; j--)
	{
	  mpc_mul_eq (value, x);
	  
	  pthread_mutex_lock (&p->mfpc_mutex[j]);
	  mpc_add_eq (value, p->mfpc[j]);
	  pthread_mutex_unlock (&p->mfpc_mutex[j]);
	}
     } 
}

/**
 * @brief Compute the value of the polynomial <code>p</code> in the point <code>x</code>
 * and save it in <code>value</code>. 
 *
 * A upper bound to the relative error of the evaluation will be stored in <code>relative_error</code>.
 * The error is computed using the formula
 * \f[
 *  n \frac{ap(|x|)}{|p(x)|} u
 * \f]
 * where \f$ap(x)\f$ is the polynomial with the coefficients equal to the moduli of the coefficients
 * of \f$p(x)\f$. 
 *
 * @param s The <code>mps_status</code> of the computation.
 * @param p The <code>monomial_poly</code> to evaluate.
 * @param x The point where the polynomial will be evaluated.
 * @param value The multiprecision complex variable where the result will be stored.
 * @param error The <code>RDPE</code> where the absolute error will be saved.
 * @param wp The working precision to use for the computation. If this value is <code>0</code> then <code>s->mpwp</code>
 * will be used.
 */
void
mps_mhorner_with_error2 (mps_status * s, mps_monomial_poly * p, mpc_t x, mpc_t value, rdpe_t error, long int wp)
{
  int i;
  rdpe_t apol, ax, u;
  cdpe_t cx;

  pthread_mutex_lock (&p->mfpc_mutex[0]);
  if (mpc_get_prec (p->mfpc[0]) < wp)
    {
      pthread_mutex_unlock (&p->mfpc_mutex[0]);
      mps_monomial_poly_raise_precision (s, p, wp);
    }
  else
    pthread_mutex_unlock (&p->mfpc_mutex[0]);

  if (mpc_get_prec (x) < wp)
    mpc_set_prec (x, wp);

  /* Set 4 * machine precision in u */
  rdpe_set_2dl (u, 1.0, 2 - wp);
  
  /* Compute the polynomial using horner */
  mps_mhorner (s, p, x, value);

  /* Compute ap(|x|) using horner */
  mpc_get_cdpe (cx, x);
  cdpe_mod (ax, cx);
      
  if (MPS_INPUT_CONFIG_IS_SPARSE (s->input_config))
    {
      mps_aparhorner (s, p->n + 1, ax, p->dap, p->spar, 
		      apol, mps_thread_get_id (s, s->pool));
    }
  else 
    {
      rdpe_set (apol, p->dap[p->n]);
      for (i = p->n - 1; i >= 0; i--)
	{
	  rdpe_mul_eq (apol, ax);
	  rdpe_add_eq (apol, p->dap[i]);
	}
    }

  /* Compute ap(|x|) / |p(x)| */
  mpc_get_cdpe (cx, value);
  cdpe_mod (ax, cx);

  rdpe_set (error, apol);
  rdpe_add_eq (error, ax);

  rdpe_mul_eq (error, u);
}

/**
 * @brief Compute the value of the polynomial <code>p</code> in the point <code>x</code>
 * and save it in <code>value</code>. 
 *
 * A upper bound to the relative error of the evaluation will be stored in <code>relative_error</code>.
 *
 * @param s The <code>mps_status</code> of the computation.
 * @param p The <code>monomial_poly</code> to evaluate.
 * @param x The point where the polynomial will be evaluated.
 * @param value The multiprecision complex variable where the result will be stored.
 * @param relative_error The <code>RDPE</code> where the relative error will be saved.
 * @param wp The working precision to use for the computation. If this value is <code>0</code> then <code>s->mpwp</code>
 * will be used.
 */
void 
mps_mhorner_with_error (mps_status * s, mps_monomial_poly * p, mpc_t x, mpc_t value, rdpe_t relative_error, long int wp)
{
  int j, my_wp;
  mpc_t ss;

  cdpe_t cdtmp;
  rdpe_t r_ss, r_value, pol_on_ss;
  rdpe_t my_eps, rtmp;

  if (wp == 0)
    my_wp = s->mpwp;
  else
    my_wp = wp;

  /* Set up precision related variables */
  rdpe_set_2dl (my_eps, 0.5, - my_wp);
  
  /* Init multiprecision temporary values */
  mpc_init2 (ss, my_wp);

  rdpe_set (relative_error, rdpe_zero);

  mpc_set (value, p->mfpc[p->n]);
  for (j = p->n - 1; j >= 0; j--)
    {
      /* Normal horner computation */
      mpc_mul (ss, value, x);
      mpc_add_eq (ss, p->mfpc[j]);

      /* Error estimate */
      mpc_get_cdpe (cdtmp, ss);
      cdpe_mod (r_ss, cdtmp);
      mpc_get_cdpe (cdtmp, value);
      cdpe_mod (r_value, cdtmp);

      rdpe_div (pol_on_ss, r_value, r_ss);
      rdpe_add (rtmp, relative_error, my_eps);
      rdpe_mul_eq (rtmp, pol_on_ss);
      rdpe_add_eq (relative_error, rtmp);

      rdpe_div (rtmp, p->dap[j], r_ss);
      rdpe_mul_eq (rtmp, my_eps);
      rdpe_add_eq (relative_error, rtmp);

      /* Update horner value */
      mpc_set (value, ss);
    }

  mpc_clear (ss);
}

/**
 * @brief Compute \f$p(z)\f$ by means of the Horner rule.
 *
 * @param st The pointer to the mps_status struct.
 * @param n The degree of the polynomial.
 * @param x The point in which the polynomial should be evaluated.
 * @param p Vector of the coefficients of the polynomial.
 * @param b Sparsity vector. If <code>b[i] == 0</code> then
 *  <code>p[i] == 0</code>.
 * @param s RDPE pointer in which the result will be stored
 * @param n_thread The index of the thread that is calling this
 *  routine. This is used to determine a safe memory are
 *  for temporary sparsity vectors.
 */
void
mps_mhorner_sparse (mps_status * s, mps_monomial_poly * p, mpc_t x,
		    mpc_t value)
{
  int m, j, i, i1, i2, q;
  mpc_t tmp, y;
  mps_boolean bi;

  /* Degree of the polynomial and sparsity vector */
  mps_boolean * b = p->spar;

  mps_boolean *spar2 = mps_boolean_valloc (p->n + 2);
  mpc_t *mfpc2 = mps_newv (mpc_t, p->n + 1);

  long int wp;

  MPS_DEBUG_THIS_CALL;

  /* pthread_mutex_lock (&p->mfpc_mutex[0]); */
  wp = mpc_get_prec (p->mfpc[0]);
  mpc_vinit2 (mfpc2, p->n + 1, wp);
  /* pthread_mutex_unlock (&p->mfpc_mutex[0]); */

  mpc_init2 (tmp, wp);
  mpc_init2 (y, wp);

  for (i = 0; i < p->n + 2; i++)
    spar2[i] = b[i];

  for (i = 0; i <= p->n; i++)
    if (b[i])
      {
	pthread_mutex_lock (&p->mfpc_mutex[i]);
	mpc_set (mfpc2[i], p->mfpc[i]);
	pthread_mutex_unlock (&p->mfpc_mutex[i]);
      }

  q = mps_intlog2 (p->n + 2);
  m = p->n + 1;
  mpc_set (y, x);
  for (j = 0; j < q; j++)
    {
      spar2[m] = false;
      m = (m + 1) >> 1;
      for (i = 0; i < m; i++)
	{
	  i2 = (i << 1) + 1;
	  i1 = i2 - 1;
	  bi = spar2[i1] || spar2[i2];
	  if (bi)
	    {
	      if (spar2[i1])
		if (spar2[i2])
		  {
		    mpc_mul (tmp, y, mfpc2[i2]);
		    mpc_add (mfpc2[i], mfpc2[i1], tmp);
		  }
		else
		  mpc_set (mfpc2[i], mfpc2[i1]);
	      else
		mpc_mul (mfpc2[i], y, mfpc2[i2]);
	    }
	  spar2[i] = bi;
	}
      spar2[m] = false;
      mpc_sqr_eq (y);
    }
  mpc_set (value, mfpc2[0]);

  mpc_clear (y);
  mpc_clear (tmp);

  mpc_vclear (mfpc2, p->n + 1);
  free (spar2);
  free (mfpc2);
}

/**
 * @brief Evaluate the polynomial p in the point x.
 *
 * @param s The <code>mps_status</code> of the computation.
 * @param p The <code>mps_monomial_poly</code> to evaluate.
 * @param x The point where the polynomial will be evaluated.
 * @param value The value computed by the function.
 */
void
mps_dhorner (mps_status * s, mps_monomial_poly * p, cdpe_t x, cdpe_t value)
{
  int j;

  cdpe_set (value, p->dpc[p->n]);
  for (j = p->n - 1; j >= 0; j--)
    {
      cdpe_mul_eq (value, x);
      cdpe_add_eq (value, p->dpc[j]);
    }
}



/**
 * @brief Evaluate the polynomial p in the point x, and give also a bound to the
 * relative error occured in the computation. 
 *
 * @param s The <code>mps_status</code> of the computation.
 * @param p The <code>mps_monomial_poly</code> to evaluate.
 * @param x The point where the polynomial will be evaluated.
 * @param value The value computed by the function.
 * @param error A bound to the absolute error of the computation.
 */
void
mps_dhorner_with_error (mps_status * s, mps_monomial_poly * p, cdpe_t x, cdpe_t value, rdpe_t error)
{
  rdpe_t ax;
  int j;

  mps_dhorner (s, p, x, value);
  
  cdpe_mod (ax, x);
  rdpe_set (error, p->dap[p->n]);
  for (j = p->n - 1; j >= 0; j--)
    {
      rdpe_mul_eq (error, ax);
      rdpe_add_eq (error, p->dap[j]);
    }

  rdpe_mul_eq_d (error, DBL_EPSILON);
}

/**
 * @brief Evaluate the polynomial p in the point x.
 *
 * @param s The <code>mps_status</code> of the computation.
 * @param p The <code>mps_monomial_poly</code> to evaluate.
 * @param x The point where the polynomial will be evaluated.
 * @param value The value computed by the function.
 */
void
mps_fhorner (mps_status * s, mps_monomial_poly * p, cplx_t x, cplx_t value)
{
  int j;

  cplx_set (value, p->fpc[p->n]);
  for (j = p->n - 1; j >= 0; j--)
    {
      cplx_mul_eq (value, x);
      cplx_add_eq (value, p->fpc[j]);
    }
}

/**
 * @brief Evaluate the polynomial p in the point x, and give also a bound to the
 * relative error occured in the computation. 
 *
 * @param s The <code>mps_status</code> of the computation.
 * @param p The <code>mps_monomial_poly</code> to evaluate.
 * @param x The point where the polynomial will be evaluated.
 * @param error A pointer to the location when an upper bound to the computation
 * error will be stored. 
 * @param value The value computed by the function.
 */
void
mps_fhorner_with_error (mps_status * s, mps_monomial_poly * p, cplx_t x, cplx_t value, double * error)
{
  int j;
  double ax = cplx_mod (x);

  mps_fhorner (s, p, x, value);

  *error = p->fap[p->n];
  for (j = p->n - 1; j >= 0; j--)
    {
      *error *= ax;
      *error += p->fap[j];
    }

  *error *= DBL_EPSILON;
}
