/*
 * This file is part of MPSolve 3.1.8
 *
 * Copyright (C) 2001-2019, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Dario Andrea Bini <bini@dm.unipi.it>
 *   Giuseppe Fiorentino <fiorent@dm.unipi.it>
 *   Leonardo Robol <leonardo.robol@sns.it>
 */


/********************************************************
   This file contains a sample of a user-defined polynomial.
   In this sample the (Mandelbrot) polynomial is defined by
   the relation
    p_{i+1}(x)=x*p_i(x)^2+1 where p_0(x)=1,        (1)
   and its derivative by
    p'_{i+1}(x)=p_i(x)^2+x*2*p_i(x)*p'_i(x).
   The degree is 1, 3, 7, 15, 31,... (2^i-1).
   The three programs below compute the values of p_i(x) and
   p'_i(x), by means of the above formulae in the float version,
   dpe version and mp version.
   The programs compute also the auxiliary variable apeps needed
   for testing the stop condition apeps>|p|,  and for computing
   the error bound rad (radius of the inclusion disk).
   The formulae for the computation of apeps and rad are
   obtained by means of a rounding error analysis of (1).
 **********************************************************/

#include <float.h>
#include <math.h>
#include <mps/mps.h>

/**
 * @brief User-defined program for the computation of \f$p\f$, \f$p'\f$.
 *
 * @param s The current mps_context
 * @param poly The mps_polynomial being solved.
 * @param root The approximation whose Newton correction shall be computed.
 * @param corr The output value where the newton correction will be stored.
 *
 * This sample computes the 'Mandelbrot polynomial by
 * means of the relation: p=1+x*p**2, starting with p=1
 */
void
mps_fnewton_usr (mps_context * s, mps_polynomial * poly, mps_approximation * root, cplx_t corr)
{
  cplx_t p, pp, pt, tmp;
  double ap, ax, eps;
  int i, m;

  cplx_t x;

  cplx_set (x, root->fvalue);
  double * rad = &root->frad;
  mps_boolean * again = &root->again;

  m = (int)(log (s->n + 1.0) / LOG2);
  if ((1 << m) <= s->n)
    m++;
  eps = (DBL_EPSILON * 4.0) * s->n;
  ax = cplx_mod (x);

  cplx_set (p, cplx_one);
  cplx_set (pp, cplx_zero);
  ap = 1.0;
  for (i = 1; i <= m; i++)
    {
      cplx_sqr (tmp, p);
      cplx_mul (pt, x, tmp);
      cplx_add_eq (pt, cplx_one);
      cplx_mul_eq (pp, x);
      cplx_mul_eq (pp, p);
      cplx_mul_eq_d (pp, 2.0);
      cplx_add_eq (pp, tmp);
      cplx_set (p, pt);
      ap = ap * ax + cplx_mod (p);
    }
  ap = ap * ax;
  cplx_div (corr, p, pp);

  *again = cplx_mod (p) > eps * ap * 3;

  *rad = s->n * (cplx_mod (p) + 3 * ap * eps) / cplx_mod (pp);
}

/**
 * @brief User-defined program for the computation of \f$p\f$, \f$p'\f$.
 *
 * @param s The current mps_context
 * @param poly The mps_polynomial being solved.
 * @param root The approximation whose Newton correction shall be computed.
 * @param corr The output value where the newton correction will be stored.
 *
 * This sample computes the 'Mandelbrot polynomial by
 * means of the relation: p=1+x*p**2, starting with p=1
 */
MPS_PRIVATE void
mps_dnewton_usr (mps_context * s, mps_polynomial * poly, mps_approximation * root, cdpe_t corr)
{
  cdpe_t p, pp, pt, tmp;
  rdpe_t ap, ax, eps, temp, apeps, atmp;
  int i, m;
  cdpe_t x;

  cdpe_set (x, root->dvalue);

  m = (int)(log (s->n + 1.0) / LOG2);
  if ((1 << m) <= s->n)
    m++;
  rdpe_set_d (eps, DBL_EPSILON);
  rdpe_mul_eq_d (eps, (double)4 * s->n);
  cdpe_mod (ax, x);

  cdpe_set (p, cdpe_one);
  cdpe_set (pp, cdpe_zero);
  rdpe_set (ap, rdpe_one);
  for (i = 1; i <= m; i++)
    {
      cdpe_sqr (tmp, p);
      cdpe_mul (pt, x, tmp);
      cdpe_add_eq (pt, cdpe_one);
      cdpe_mul_eq (pp, x);
      cdpe_mul_eq (pp, p);
      cdpe_mul_eq_d (pp, 2.0);
      cdpe_add_eq (pp, tmp);
      cdpe_set (p, pt);
      rdpe_mul_eq (ap, ax);
      cdpe_mod (atmp, p);
      rdpe_add_eq (ap, atmp);
    }
  rdpe_mul_eq (ap, ax);
  cdpe_div (corr, p, pp);

  cdpe_mod (temp, p);
  rdpe_mul (apeps, ap, eps);
  rdpe_mul_eq_d (apeps, 3.0);
  root->again = rdpe_gt (temp, apeps);

  rdpe_add (root->drad, temp, apeps);
  rdpe_mul_eq_d (root->drad, (double)s->n);
  cdpe_mod (temp, pp);
  rdpe_div_eq (root->drad, temp);
  if (rdpe_eq (root->drad, rdpe_zero))
    rdpe_mul (root->drad, ax, eps);
}

/**
 * @brief User-defined program for the computation of \f$p\f$, \f$p'\f$.
 *
 * @param s The current mps_context
 * @param poly The mps_polynomial being solved.
 * @param root The approximation whose Newton correction shall be computed.
 * @param corr The output value where the newton correction will be stored.
 *
 * This sample computes the 'Mandelbrot polynomial by
 * means of the relation: p=1+x*p**2, starting with p=1
 */
MPS_PRIVATE void
mps_mnewton_usr (mps_context * s, mps_polynomial * poly, mps_approximation * root, mpcf_t corr, long int wp)
{
  int i, m;
  rdpe_t ap, ax, eps, temp, apeps, atmp;
  cdpe_t ctmp;
  mpcf_t p, pp, pt, tmp;

  mpcf_init2 (p, s->mpwp);
  mpcf_init2 (pp, s->mpwp);
  mpcf_init2 (pt, s->mpwp);
  mpcf_init2 (tmp, s->mpwp);

  m = (int)(log (s->n + 1.0) / LOG2);
  if ((1 << m) <= s->n)
    m++;
  rdpe_set (eps, s->mp_epsilon);
  rdpe_mul_eq_d (eps, (double)4 * s->n);
  mpcf_get_cdpe (ctmp, root->mvalue);
  cdpe_mod (ax, ctmp);

  mpcf_set_ui (p, 1, 0);
  mpcf_set_ui (pp, 0, 0);
  rdpe_set (ap, rdpe_one);
  for (i = 1; i <= m; i++)
    {
      mpcf_sqr (tmp, p);
      mpcf_mul (pt, root->mvalue, tmp);
      mpcf_add_eq_ui (pt, 1, 0);
      mpcf_mul_eq (pp, root->mvalue);
      mpcf_mul_eq (pp, p);
      mpcf_mul_eq_ui (pp, 2);
      mpcf_add_eq (pp, tmp);
      mpcf_set (p, pt);
      rdpe_mul_eq (ap, ax);
      mpcf_get_cdpe (ctmp, p);
      cdpe_mod (atmp, ctmp);
      rdpe_add_eq (ap, atmp);
    }
  rdpe_mul_eq (ap, ax);
  mpcf_div (corr, p, pp);

  mpcf_get_cdpe (ctmp, p);
  cdpe_mod (temp, ctmp);
  rdpe_mul (apeps, ap, eps);
  rdpe_mul_eq_d (apeps, 3.0);
  root->again = rdpe_gt (temp, apeps);

  rdpe_add (root->drad, temp, apeps);
  rdpe_mul_eq_d (root->drad, (double)s->n);
  mpcf_get_cdpe (ctmp, pp);
  cdpe_mod (temp, ctmp);
  rdpe_div_eq (root->drad, temp);
  if (rdpe_eq (root->drad, rdpe_zero))
    rdpe_mul (root->drad, ax, eps);

  mpcf_clear (tmp);
  mpcf_clear (pt);
  mpcf_clear (pp);
  mpcf_clear (p);
}

mps_boolean
mps_feval_usr (mps_context * ctx, mps_polynomial * p, cplx_t x, cplx_t value, double * error)
{
  int i;
  int m = (int)(log (p->degree + 1.0) / LOG2);
  double ax = cplx_mod (x);
  cplx_t tmp;

  if ((1 << m) <= p->degree)
    m++;

  cplx_set (value, cplx_one);

  if (error)
    *error = 0.0;

  for (i = 1; i <= m; i++)
    {
      cplx_sqr (tmp, value);
      cplx_mul (value, x, tmp);
      cplx_add_eq (value, cplx_one);

      if (error)
        *error = *error * ax + cplx_mod (value);
    }

  if (error)
    *error *= DBL_EPSILON;

  return true;
}

mps_boolean
mps_deval_usr (mps_context * ctx, mps_polynomial * p, cdpe_t x, cdpe_t value, rdpe_t error)
{
  int i;
  int m = (int)(log (p->degree + 1.0) / LOG2);
  rdpe_t ax, rtmp;
  cdpe_t tmp;

  if ((1 << m) <= p->degree)
    m++;

  cdpe_mod (ax, x);
  cdpe_set (value, cdpe_one);
  cdpe_mod (error, value);

  for (i = 1; i <= m; i++)
    {
      cdpe_sqr (tmp, value);
      cdpe_mul (value, x, tmp);
      cdpe_add_eq (value, cdpe_one);

      rdpe_mul_eq (error, ax);
      cdpe_mod (rtmp, value);
      rdpe_add_eq (error, rtmp);
    }

  rdpe_mul_eq_d (error, DBL_EPSILON);
  return true;
}

mps_boolean
mps_meval_usr (mps_context * ctx, mps_polynomial * p, mpcf_t x, mpcf_t value, rdpe_t error)
{
  int i;
  int m = (int)(log (p->degree + 1.0) / LOG2);
  rdpe_t ax, rtmp;
  mpcf_t tmp;
  long int wp = mpcf_get_prec (x);

  /* Correct the working precision in case of a limited precision polynomial (quite unlikely
   * in the Mandelbrot case, but still. */
  if (p->prec > 0 && p->prec < wp)
    wp = p->prec;

  if ((1 << m) <= p->degree)
    m++;

  mpcf_init2 (tmp, wp);

  mpcf_rmod (ax, x);
  mpcf_set_ui (value, 1U, 0U);
  mpcf_rmod (error, value);

  for (i = 1; i <= m; i++)
    {
      mpcf_sqr (tmp, value);
      mpcf_mul (value, x, tmp);
      mpcf_add_eq_ui (value, 1U, 0U);

      rdpe_mul_eq (error, ax);
      mpcf_rmod (rtmp, value);
      rdpe_add_eq (error, rtmp);
    }

  rdpe_set_2dl (rtmp, 1.0, -wp);
  rdpe_mul_eq (error, rtmp);

  mpcf_clear (tmp);

  return true;
}
