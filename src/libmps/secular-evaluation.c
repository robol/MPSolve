/*
 * This file is part of MPSolve 3.0
 *
 * Copyright (C) 2001-2013, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors: 
 *   Dario Andrea Bini <bini@dm.unipi.it>
 *   Giuseppe Fiorentino <fiorent@dm.unipi.it>
 *   Leonardo Robol <robol@mail.dm.unipi.it>
 */

#include <mps/mps.h>

/**
 * @brief Evaluate a secular equation <code>sec</code> in the point <code>x</code>
 *
 * @param s The <code>mps_context</code> of the computation.
 * @param sec The secular equation to evaluate.
 * @param x The point in which the secular equation must be evaluated.
 * @param value The value of the secular equation in the pointer <code>x</code>.
 */
mps_boolean
mps_secular_feval (mps_context * s, mps_polynomial * p, cplx_t x, cplx_t value)
{
  mps_secular_equation * sec = MPS_SECULAR_EQUATION (p);
  cplx_t ctmp;
  int i;
  
  cplx_set (value, cplx_zero);
  
  for (i = 0; i < s->n; i++)
    {
      cplx_sub (ctmp, x, sec->bfpc[i]);
      if (cplx_eq_zero (ctmp))
        return false;
      cplx_div (ctmp, sec->afpc[i], ctmp);
      cplx_add_eq (value, ctmp);
    }

  cplx_sub_eq (value, cplx_one);
  return true;
}

mps_boolean
mps_secular_feval_derivative (mps_context * s, mps_polynomial * p, cplx_t x, cplx_t value)
{
  mps_secular_equation * sec = MPS_SECULAR_EQUATION (p);
  cplx_t ctmp;
  int i;

  cplx_set (value, cplx_zero);

  for (i = 0; i < s->n; i++)
    {
      /* Compute 1 / (x - b_i) */
      cplx_sub (ctmp, x, sec->bfpc[i]);

      if (cplx_eq_zero (ctmp))
        return false;

      cplx_inv_eq (ctmp);
      cplx_mul_eq (ctmp, ctmp);

      /* Compute a_i / (x - b_i) */
      cplx_mul_eq (ctmp, sec->afpc[i]);

      /* Sum to the secular eqation */
      cplx_sub_eq (value, ctmp);
    }

  return true;
}

/**
 * @brief Evaluate a secular equation <code>sec</code> in the point <code>x</code>, 
 * estimating the error on the evaluation.
 *
 * @param s The <code>mps_context</code> of the computation.
 * @param sec The secular equation to evaluate.
 * @param x The point in which the secular equation must be evaluated.
 * @param value The value of the secular equation in the pointer <code>x</code>.
 * @param error The absolute error on the evaluation.
 */
mps_boolean
mps_secular_feval_with_error (mps_context * s, mps_polynomial * p, cplx_t x, cplx_t value,
                              double * error)
{
  mps_secular_equation * sec = MPS_SECULAR_EQUATION (p);
  cplx_t ctmp;
  int i;
  
  *error = 0.0f;
  cplx_set (value, cplx_zero);
  
  for (i = 0; i < s->n; i++)
    {
      cplx_sub (ctmp, x, sec->bfpc[i]);

      if (cplx_eq_zero (ctmp))
        return false;

      cplx_div (ctmp, sec->afpc[i], ctmp);
      cplx_add_eq (value, ctmp);
      *error += cplx_mod (ctmp) * (i + 2);
    }

  cplx_sub_eq (value, cplx_one);
  *error += 1.0f;

  *error *= 4.0f * DBL_EPSILON;

  return true;
}

/**
 * @brief Evaluate a secular equation <code>sec</code> in the point <code>x</code>
 *
 * @param s The <code>mps_context</code> of the computation.
 * @param sec The secular equation to evaluate.
 * @param x The point in which the secular equation must be evaluated.
 * @param value The value of the secular equation in the point <code>x</code>.
 */
mps_boolean
mps_secular_deval (mps_context * s, mps_polynomial * p, cdpe_t x, cdpe_t value)
{
  mps_secular_equation *sec = MPS_SECULAR_EQUATION (p);
  cdpe_t ctmp;
  int i;

  cdpe_set (value, cdpe_zero);

  for (i = 0; i < s->n; i++)
    {
      cdpe_sub (ctmp, x, sec->bdpc[i]);

      if (cdpe_eq_zero (ctmp))
        return false;

      cdpe_div (ctmp, sec->adpc[i], ctmp);
      cdpe_add_eq (value, ctmp);
    }

  cdpe_sub_eq (value, cdpe_one);

  return true;
}

mps_boolean
mps_secular_deval_derivative (mps_context * s, mps_polynomial * p, cdpe_t x, cdpe_t value)
{
  mps_secular_equation *sec= MPS_SECULAR_EQUATION (p);
  cdpe_t ctmp;
  int i;

  cdpe_set (value, cdpe_zero);

  for (i = 0; i < s->n; i++)
    {
      /* Compute 1 / (x - b_i) */
      cdpe_sub (ctmp, x, sec->bdpc[i]);

      if (cdpe_eq_zero (ctmp))
        return false;

      cdpe_inv_eq (ctmp);
      cdpe_mul_eq (ctmp, ctmp);

      /* Compute a_i / (x - b_i) */
      cdpe_mul_eq (ctmp, sec->adpc[i]);

      /* Sum to the secular eqation */
      cdpe_sub_eq (value, ctmp);
    }

  return true;
}


/**
 * @brief Evaluate a secular equation <code>sec</code> in the point <code>x</code>
 *
 * @param s The <code>mps_context</code> of the computation.
 * @param sec The secular equation to evaluate.
 * @param x The point in which the secular equation must be evaluated.
 * @param value The value of the secular equation in the point <code>x</code>.
 * @param error A bound to the module of the relative error occurred in the computation.
 */
mps_boolean
mps_secular_deval_with_error (mps_context * s, mps_polynomial * p,
                              cdpe_t x, cdpe_t value, rdpe_t error)
{
  mps_secular_equation * sec = MPS_SECULAR_EQUATION (p);
  cdpe_t ctmp;
  rdpe_t rtmp;
  int i;

  cdpe_set (value, cdpe_zero);
  rdpe_set (error, rdpe_zero);

  for (i = 0; i < s->n; i++)
    {
      cdpe_sub (ctmp, x, sec->bdpc[i]);
      if (cdpe_eq_zero (ctmp))
        return false;
      cdpe_div (ctmp, sec->adpc[i], ctmp);
      cdpe_mod (rtmp, ctmp);
      cdpe_add_eq (value, ctmp);
      rdpe_mul_eq_d (rtmp, i + 2);
      rdpe_add_eq (error, rtmp);
    }

  cdpe_sub_eq (value, cdpe_one);
  rdpe_add_eq (error, rdpe_one);

  rdpe_mul_eq_d (error, 4.0f * DBL_EPSILON);

  return true;
}


/**
 * @brief Evaluate a secular equation <code>sec</code> in the point <code>x</code>.
 *
 * @param s The <code>mps_context</code> of the computation.
 * @param sec The secular equation to evaluate.
 * @param x The point in which the sceular equation must be evaluated.
 * @param value The value of the secular equation in the point <code>x</code>.
 */
mps_boolean
mps_secular_meval (mps_context * s, mps_polynomial * p, mpc_t x, mpc_t value)
{
  mps_secular_equation * sec = MPS_SECULAR_EQUATION (p);
  mps_boolean success = true;
  mpc_t ctmp;
  unsigned int wp = mpc_get_prec (x);
  int i;

  mpc_init2  (ctmp, wp);
  mpc_set_ui (value, 0U, 0U);

  for (i = 0; i < s->n; ++i)
    {
      mpc_sub (ctmp, x, sec->bmpc[i]);
      if (mpc_eq_zero (ctmp))
        {
          success = false;
          goto cleanup;
        }

      mpc_div (ctmp, sec->ampc[i], ctmp);
      mpc_add_eq (value, ctmp);
    }
  
  mpc_sub_eq_ui (value, 1U, 0U);

 cleanup:
  mpc_clear (ctmp);
  return success;
}

/**
 * @brief Evaluate a secular equation <code>sec</code> in the point <code>x</code>.
 *
 * @param s The <code>mps_context</code> of the computation.
 * @param sec The secular equation to evaluate.
 * @param x The point in which the sceular equation must be evaluated.
 * @param value The value of the secular equation in the point <code>x</code>.
 * @param error A bound to the absolute value of the error introduced in the computation.
 */
mps_boolean
mps_secular_meval_with_error (mps_context * s, mps_polynomial * p, mpc_t x, mpc_t value, rdpe_t error)
{
  mps_secular_equation * sec = MPS_SECULAR_EQUATION (p);
  mpc_t ctmp;
  rdpe_t rtmp, ax;
  cdpe_t cdtmp;
  unsigned int wp = mpc_get_prec (x);
  int i;

  mps_boolean successful_evaluation = true;

  if (mpc_get_prec (sec->ampc[0]) < wp)
    mps_polynomial_raise_data (s, p, wp);

  mpc_init2  (ctmp, wp);
  mpc_set_ui (value, 0U, 0U);
  mpc_set_prec (value, wp);

  /* Get |x| */
  mpc_rmod (ax, x);

  rdpe_set (error, rdpe_zero);

  for (i = 0; i < s->n; i++)
    {
      mpc_sub (ctmp, x, sec->bmpc[i]);

      if (mpc_eq_zero (ctmp))
        {
          successful_evaluation = false;
          goto cleanup;
        }

      mpc_div (ctmp, sec->ampc[i], ctmp);
      mpc_add_eq (value, ctmp);

      mpc_get_cdpe (cdtmp, ctmp); 
      cdpe_mod (rtmp, cdtmp); 
      rdpe_mul_eq_d (rtmp, i + 2); 
      rdpe_add_eq (error, rtmp);
    }
  
  mpc_sub_eq_ui (value, 1U, 0U);
  rdpe_add_eq (error, rdpe_one);

  rdpe_set_2dl (rtmp, 4.0, 1 - (long int) wp);
  rdpe_mul_eq (error, rtmp);

 cleanup:
  
  mpc_clear (ctmp);

  return successful_evaluation;
}

