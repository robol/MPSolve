#include <mps/core.h>

/**
 * @brief Return a newly allocated mps_monomial_poly of the given degree.
 */
mps_monomial_poly *
mps_monomial_poly_new (mps_status * s, long int degree)
{
  mps_monomial_poly  * mp = mps_malloc (sizeof (mps_monomial_poly));
  
  /* Set the degree of the polynomial */
  mp->n = degree;

  /* Allocate the space for the coefficients of the polynomial, all
  * the floating point versions. */
  mp->spar = mps_boolean_valloc (degree + 2);
  mp->fpc  = cplx_valloc (degree + 1);
  mp->fpr  = double_valloc (degree + 1);
  mp->dpr  = rdpe_valloc (degree + 1);
  mp->dpc  = cdpe_valloc (degree + 1);
  mp->mfpc = mpc_valloc (degree + 1);
  mp->mfpr = mpf_valloc (degree + 1);

  mpf_vinit (mp->mfpr, degree + 1);
  mpc_vinit (mp->mfpc, degree + 1);

  /* Allocate space for the moduli of the coefficients */
  mp->fap = double_valloc (degree + 1);
  mp->dap = rdpe_valloc (degree + 1);

  /* Allocate space for the coefficients of the derivative */
  mp->fppc = cplx_valloc (degree);
  mp->mfppc = mpc_valloc (degree + 1);
  mpc_vinit (mp->mfppc, degree + 1);

  /* Allocate space for the coefficients initially parsed as
   * exact */
  mp->initial_mqp_r = mpq_valloc (degree + 1);
  mp->initial_mqp_i = mpq_valloc (degree + 1);

  mpq_vinit (mp->initial_mqp_r, degree + 1);
  mpq_vinit (mp->initial_mqp_i, degree + 1);

  return mp;
}

/**
 * @brief Free a instance of <code>mps_monomial_poly</code> previously
 * allocated with <code>mps_monomial_poly_new()</code>.
 */
void
mps_monomial_poly_free (mps_status * s, mps_monomial_poly * mp)
{
  mps_boolean_vfree (mp->spar);
  double_vfree (mp->fpr);
  cplx_vfree (mp->fpc);
  rdpe_vfree (mp->dpr);
  cdpe_vfree (mp->dpc);

  double_vfree (mp->fap);
  rdpe_vfree (mp->dap);

  mpf_vclear (mp->mfpr, mp->n + 1);
  mpc_vclear (mp->mfpc, mp->n + 1);

  mpf_vfree (mp->mfpr);
  mpc_vfree (mp->mfpc);

  mpq_vclear (mp->initial_mqp_r, mp->n + 1);
  mpq_vclear (mp->initial_mqp_i, mp->n + 1);

  mpq_vfree (mp->initial_mqp_r);
  mpq_vfree (mp->initial_mqp_i);

  cplx_vfree (mp->fppc);
  mpc_vclear (mp->mfppc, mp->n + 1);
  mpc_vfree (mp->mfppc);

  free (mp);
}

/**
 * @brief Raise the precision bits of the multiprecision fields of the 
 * polynomial to selected value.
 *
 * @param s  The status of the computation.
 * @param mp The polynomial that need the precision raised.
 * @param prec The selected bits of precision.
 */
void
mps_monomial_poly_raise_precision (mps_status * s, mps_monomial_poly * mp, long int prec)
{
  int k;

  /* raise the precision of  mfpc */
  if (s->data_type[0] != 'u')
    for (k = 0; k < mp->n + 1; k++)
      {
	mpc_set_prec (mp->mfpc[k], prec);
      }

  /* Raise the precision of p' */
  if (s->data_type[0] == 's')
    for (k = 0; k < s->n; k++)
      if (mp->spar[k + 1])
        {
          mpc_set_prec (mp->mfppc[k], prec);
          mpc_mul_ui (mp->mfppc[k], mp->mfpc[k + 1], k + 1);
        }
}

/**
 * @brief This routine can be used to set the i-th coefficients of the 
 * polynomial with a multiprecision rational number. 
 *
 * @param s The mps_status the will be used for the computation. This shall be passed
 * along with the polynomial because it manages user interaction and can perform
 * error handling.
 *
 * @param mp The monomial_poly that will hold the coefficient.
 * @param i  The exponent of the monomial to which the coefficient is related.
 * @param real_part The real part of the coefficient.
 * @param imag_part The imaginary part of the coefficients.
 */
void
mps_monomial_poly_set_coefficient_q (mps_status * s, mps_monomial_poly * mp, long int i, 
				     mpq_t real_part, mpq_t imag_part)
{
  /* Set the coefficients first */
  mpq_set (mp->initial_mqp_r[i], real_part);
  mpq_set (mp->initial_mqp_i[i], imag_part);

  /* Then update the other coefficients */
  mpf_set_q (mpc_Re (mp->mfpc[i]), real_part);
  mpf_set_q (mpc_Im (mp->mfpc[i]), imag_part);

  mpc_get_cdpe (mp->dpc[i], mp->mfpc[i]);
  mpc_get_cplx (mp->fpc[i], mp->mfpc[i]);
}
