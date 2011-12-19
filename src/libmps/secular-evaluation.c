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
      *error += cplx_mod (ctmp);
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
  cdpe_t cdtmp;
  rdpe_t rtmp;
  unsigned int wp = mpc_get_prec (x);
  int i;

  mpc_init2  (ctmp, wp);
  mpc_set_ui (value, 0U, 0U);

  rdpe_set (error, rdpe_zero);

  for (i = 0; i < s->n; ++i)
    {
      mpc_sub (ctmp, x, sec->bmpc[i]);
      mpc_div (ctmp, sec->ampc[i], ctmp);
      mpc_add_eq (value, ctmp);

      mpc_get_cdpe (cdtmp, ctmp);
      cdpe_mod (rtmp, cdtmp);
      rdpe_add_eq (error, rtmp);
    }
  
  mpc_sub_eq_ui (value, 1U, 0U);
  rdpe_add_eq (error, rdpe_one);

  rdpe_set_2dl (rtmp, 1.0, 1 - (long int) wp);
  rdpe_mul_eq (error, rtmp);
  
  mpc_clear (ctmp);
}

