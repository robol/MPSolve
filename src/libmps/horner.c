#include <mps/core.h>

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
  mpc_set (value, p->mfpc[p->n]);
  for (j = s->n - 1; j >= 0; j--)
    {
      mpc_mul_eq (value, x);
      mpc_add_eq (value, p->mfpc[j]);
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

  if (mpc_get_prec (p->mfpc[0]) < wp)
    {
      for(i = 0; i < p->n; i++)
	mpc_set_prec (p->mfpc[i], wp);
    }

  /* Set 4 * machine precision in u */
  rdpe_set_2dl (u, 1.0, 2 - wp);
  
  /* Compute the polynomial using horner */
  mps_mhorner (s, p, x, value);
  
  /* Compute ap(|x|) using horner */
  mpc_get_cdpe (cx, x);
  cdpe_mod (ax, cx);

  rdpe_set (apol, p->dap[p->n]);
  for (i = s->n - 1; i >= 0; i--)
    {
      rdpe_mul_eq (apol, ax);
      rdpe_add_eq (apol, p->dap[i]);
    }

  /* Compute ap(|x|) / |p(x)| */
  mpc_get_cdpe (cx, value);
  cdpe_mod (ax, cx);

  rdpe_set (error, apol);
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
  for (j = s->n - 1; j >= 0; j--)
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
 * @brief Evaluate the polynomial p in the point x.
 *
 * @param s The <code>mps_status</code> of the computation.
 * @param p The <code>mps_monomial_poly</code> to evaluate.
 * @param z The point where the polynomial will be evaluated.
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
 * @param z The point where the polynomial will be evaluated.
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
 * @param z The point where the polynomial will be evaluated.
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
 * @param z The point where the polynomial will be evaluated.
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

/* Sparse version of horner... */
	  /* if (MPS_INPUT_CONFIG_IS_SPARSE (s->input_config)) */
	  /*   { */
	  /*     int skip = 0; */
	  /*     cdpe_set (pol, p->dpc[s->n]); */
	  /*     for (j = s->n - 1; j >= 0; j--) */
	  /* 	{ */
	  /* 	  if (!p->spar[j]) */
	  /* 	    { */
	  /* 	      int power = 0, t, remaining; */

	  /* 	      /\* Determine how many coefficients we shall skip *\/ */
	  /* 	      skip = j; */
	  /* 	      while (!p->spar[j] && (j >= 0)) { j--; } */
	  /* 	      skip = skip - j + 1; */

	  /* 	      /\* Smarter less complex way needed *\/ */
	  /* 	      cdpe_set (ss, cdpe_one); */
	  /* 	      while (power != skip) */
	  /* 		{ */
	  /* 		  cdpe_set (ctmp, sec->bdpc[i]); */
	  /* 		  t = 2; */
	  /* 		  remaining = skip - power; */
	  /* 		  while (t < remaining) */
	  /* 		    { */
	  /* 		      cdpe_mul_eq (ctmp, ctmp); */
	  /* 		      t *= 2; */
	  /* 		    } */
	  /* 		  t /= 2; */
	  /* 		  power += t; */
	  /* 		  cdpe_mul_eq (ss, ctmp); */
	  /* 		} */

	  /* 	      cdpe_add_eq (ss, p->dpc[j]); */
	  /* 	      cdpe_div (ctmp, pol, ss); */
	  /* 	      cdpe_mod (rtmp, ctmp); */
			
	  /* 	      rdpe_set (rtmp2, root_epsilon); */
	  /* 	      rdpe_add_eq (rtmp2, regeneration_epsilon); */
	  /* 	      rdpe_mul_eq (rtmp, rtmp2); */
	  /* 	      rdpe_add_eq (regeneration_epsilon, rtmp); */
		    
	  /* 	      cdpe_div (ctmp, p->dpc[j], ss); */
	  /* 	      cdpe_mod (rtmp, ctmp); */
	  /* 	      rdpe_mul_eq (rtmp, regeneration_epsilon); */
	  /* 	      rdpe_add_eq (regeneration_epsilon, rtmp); */
		    
	  /* 	      cdpe_set (pol, ss); */
	  /* 	    } */
	  /* 	} */
	  /*   } */
