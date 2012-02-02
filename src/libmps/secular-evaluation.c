#include <mps/mps.h>

/**
 * @brief Evaluate a secular equation <code>sec</code> in the point <code>x</code>
 *
 * @param s The <code>mps_status</code> of the computation.
 * @param sec The secular equation to evaluate.
 * @param x The point in which the secular equation must be evaluated.
 * @param value The value of the secular equation in the pointer <code>x</code>.
 */
void
mps_secular_feval (mps_status * s, mps_secular_equation * sec, cplx_t x, cplx_t value)
{
  cplx_t ctmp;
  int i;
  
  cplx_set (value, cplx_zero);
  
  for (i = 0; i < s->n; i++)
    {
      cplx_sub (ctmp, x, sec->bfpc[i]);
      cplx_div (ctmp, sec->afpc[i], ctmp);
      cplx_add_eq (value, ctmp);
    }

  cplx_sub_eq (value, cplx_one);
}

/**
 * @brief Evaluate a secular equation <code>sec</code> in the point <code>x</code>, 
 * estimating the error on the evaluation.
 *
 * @param s The <code>mps_status</code> of the computation.
 * @param sec The secular equation to evaluate.
 * @param x The point in which the secular equation must be evaluated.
 * @param value The value of the secular equation in the pointer <code>x</code>.
 * @param error The absolute error on the evaluation.
 */
void
mps_secular_feval_with_error (mps_status * s, mps_secular_equation * sec, cplx_t x, cplx_t value,
			      double * error)
{
  cplx_t ctmp;
  int i;
  
  *error = 0.0f;
  cplx_set (value, cplx_zero);
  
  for (i = 0; i < s->n; i++)
    {
      cplx_sub (ctmp, x, sec->bfpc[i]);
      cplx_div (ctmp, sec->afpc[i], ctmp);
      cplx_add_eq (value, ctmp);
      *error += cplx_mod (ctmp) * (i + 2);
    }

  cplx_sub_eq (value, cplx_one);
  *error += 1.0f;

  *error *= 4.0f * DBL_EPSILON;
}

/**
 * @brief Evaluate a secular equation <code>sec</code> in the point <code>x</code>
 *
 * @param s The <code>mps_status</code> of the computation.
 * @param sec The secular equation to evaluate.
 * @param x The point in which the secular equation must be evaluated.
 * @param value The value of the secular equation in the point <code>x</code>.
 */
void
mps_secular_deval (mps_status * s, mps_secular_equation * sec, cdpe_t x, cdpe_t value)
{
  cdpe_t ctmp;
  int i;

  cdpe_set (value, cdpe_zero);

  for (i = 0; i < s->n; i++)
    {
      cdpe_sub (ctmp, x, sec->bdpc[i]);
      cdpe_div (ctmp, sec->adpc[i], ctmp);
      cdpe_add_eq (value, ctmp);
    }

  cdpe_sub_eq (value, cdpe_one);
}

/**
 * @brief Evaluate a secular equation <code>sec</code> in the point <code>x</code>
 *
 * @param s The <code>mps_status</code> of the computation.
 * @param sec The secular equation to evaluate.
 * @param x The point in which the secular equation must be evaluated.
 * @param value The value of the secular equation in the point <code>x</code>.
 * @param error A bound to the module of the relative error occurred in the computation.
 */
void
mps_secular_deval_with_error (mps_status * s, mps_secular_equation * sec, 
			      cdpe_t x, cdpe_t value, rdpe_t error)
{
  cdpe_t ctmp;
  rdpe_t rtmp;
  int i;

  cdpe_set (value, cdpe_zero);
  rdpe_set (error, rdpe_zero);

  for (i = 0; i < s->n; i++)
    {
      cdpe_sub (ctmp, x, sec->bdpc[i]);
      cdpe_div (ctmp, sec->adpc[i], ctmp);
      cdpe_mod (rtmp, ctmp);
      cdpe_add_eq (value, ctmp);
      rdpe_mul_eq_d (error, i + 2);
      rdpe_add_eq (error, rtmp);
    }

  cdpe_sub_eq (value, cdpe_one);
  rdpe_add_eq (error, rdpe_one);

  rdpe_mul_eq_d (error, 4.0f * DBL_EPSILON);
}


/**
 * @brief Evaluate a secular equation <code>sec</code> in the point <code>x</code>.
 *
 * @param s The <code>mps_status</code> of the computation.
 * @param sec The secular equation to evaluate.
 * @param x The point in which the sceular equation must be evaluated.
 * @param value The value of the secular equation in the point <code>x</code>.
 */
void
mps_secular_meval (mps_status * s, mps_secular_equation * sec, mpc_t x, mpc_t value)
{
  mpc_t ctmp;
  unsigned int wp = mpc_get_prec (x);
  int i;

  mpc_init2  (ctmp, wp);
  mpc_set_ui (value, 0U, 0U);

  for (i = 0; i < s->n; ++i)
    {
      mpc_sub (ctmp, x, sec->bmpc[i]);
      mpc_div (ctmp, sec->ampc[i], ctmp);
      mpc_add_eq (value, ctmp);
    }
  
  mpc_sub_eq_ui (value, 1U, 0U);
  
  mpc_clear (ctmp);
}

/**
 * @brief Evaluate a secular equation <code>sec</code> in the point <code>x</code>.
 *
 * @param s The <code>mps_status</code> of the computation.
 * @param sec The secular equation to evaluate.
 * @param x The point in which the sceular equation must be evaluated.
 * @param value The value of the secular equation in the point <code>x</code>.
 * @param error A bound to the absolute value of the error introduced in the computation.
 */
void
mps_secular_meval_with_error (mps_status * s, mps_secular_equation * sec, mpc_t x, mpc_t value, rdpe_t error)
{
  mpc_t ctmp;
  unsigned int wp = mpc_get_prec (x);
  int i;
  mpf_t ftmp;
  mpf_t merror;

  mpf_init2 (ftmp, wp);
  mpf_init2 (merror, wp);

  mpc_init2 (ctmp, wp);
  mpc_set_ui (value, 0U, 0U);

  mpf_set_ui (merror, 0U);

  for (i = 0; i < s->n; ++i)
    {
      mpc_sub (ctmp, x, sec->bmpc[i]);
      mpc_div (ctmp, sec->ampc[i], ctmp);
      mpc_add_eq (value, ctmp);

      mpc_mod (ftmp, ctmp);
      mpf_mul_eq_ui (ftmp, i + 2);
      mpf_add_eq (merror, ftmp);
    }
  
  mpc_sub_eq_ui (value, 1U, 0U);
  mpf_add_eq_ui (merror, 1U);

  mpf_set_2dl (ftmp, 1.0, 1 - wp);
  mpf_mul_eq (merror, ftmp);

  mpf_get_rdpe (error, merror);

  /* MPS_DEBUG_MPC (s, wp, value, "Evaluation of secular equation"); */
  /* MPS_DEBUG_MPF (s, wp, merror, "Error on the secular equation"); */

  /* mps_mhorner_with_error2 (s, s->monomial_poly, x, value, error, wp); */

  mpf_clear (merror);
  mpf_clear (ftmp);
  mpc_clear (ctmp);
}

