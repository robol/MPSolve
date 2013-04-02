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

/******************************************************
*              SUBROUTINE FNEWTON_USR                 *   
******************************************************* 
* user-defined program for the computation of p, p',  *
* and ap corr = p/p1, rad = n*(ABS(p)+ap)/ABS(p'),    *
* again = ABS(p)>ap.                                  *
* This sample computes the 'Mandelbrot polynomial by  *
* means of the relation: p=1+x*p**2, starting with p=1*
******************************************************/
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

  m = (int) (log (s->n + 1.0) / LOG2);
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

/******************************************************
*              SUBROUTINE DNEWTON_USR                 *   
******************************************************* 
 DPE computation
******************************************************/
void
mps_dnewton_usr (mps_context * s, mps_polynomial * poly, mps_approximation * root, cdpe_t corr)
{
  cdpe_t p, pp, pt, tmp;
  rdpe_t ap, ax, eps, temp, apeps, atmp;
  int i, m;
  cdpe_t x;

  cdpe_set (x, root->dvalue);

  m = (int) (log (s->n + 1.0) / LOG2);
  if ((1 << m) <= s->n)
    m++;
  rdpe_set_d (eps, DBL_EPSILON);
  rdpe_mul_eq_d (eps, (double) 4 * s->n);
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
  rdpe_mul_eq_d (root->drad, (double) s->n);
  cdpe_mod (temp, pp);
  rdpe_div_eq (root->drad, temp);
  if (rdpe_eq (root->drad, rdpe_zero))
    rdpe_mul (root->drad, ax, eps);
}

/******************************************************
*              SUBROUTINE MNEWTON_USR                 *   
******************************************************* 
 multiprecision computation
******************************************************/
void
mps_mnewton_usr (mps_context * s, mps_polynomial * poly, mps_approximation * root, mpc_t corr)
{
  int i, m;
  rdpe_t ap, ax, eps, temp, apeps, atmp;
  cdpe_t ctmp;
  mpc_t p, pp, pt, tmp;

  mpc_init2 (p, s->mpwp);
  mpc_init2 (pp, s->mpwp);
  mpc_init2 (pt, s->mpwp);
  mpc_init2 (tmp, s->mpwp);

  m = (int) (log (s->n + 1.0) / LOG2);
  if ((1 << m) <= s->n)
    m++;
  rdpe_set (eps, s->mp_epsilon);
  rdpe_mul_eq_d (eps, (double) 4 * s->n);
  mpc_get_cdpe (ctmp, root->mvalue);
  cdpe_mod (ax, ctmp);

  mpc_set_ui (p, 1, 0);
  mpc_set_ui (pp, 0, 0);
  rdpe_set (ap, rdpe_one);
  for (i = 1; i <= m; i++)
    {
      mpc_sqr (tmp, p);
      mpc_mul (pt, root->mvalue, tmp);
      mpc_add_eq_ui (pt, 1, 0);
      mpc_mul_eq (pp, root->mvalue);
      mpc_mul_eq (pp, p);
      mpc_mul_eq_ui (pp, 2);
      mpc_add_eq (pp, tmp);
      mpc_set (p, pt);
      rdpe_mul_eq (ap, ax);
      mpc_get_cdpe (ctmp, p);
      cdpe_mod (atmp, ctmp);
      rdpe_add_eq (ap, atmp);
    }
  rdpe_mul_eq (ap, ax);
  mpc_div (corr, p, pp);

  mpc_get_cdpe (ctmp, p);
  cdpe_mod (temp, ctmp);
  rdpe_mul (apeps, ap, eps);
  rdpe_mul_eq_d (apeps, 3.0);
  root->again = rdpe_gt (temp, apeps);

  rdpe_add (root->drad, temp, apeps);
  rdpe_mul_eq_d (root->drad, (double) s->n);
  mpc_get_cdpe (ctmp, pp);
  cdpe_mod (temp, ctmp);
  rdpe_div_eq (root->drad, temp);
  if (rdpe_eq (root->drad, rdpe_zero))
    rdpe_mul (root->drad, ax, eps);

  mpc_clear (tmp);
  mpc_clear (pt);
  mpc_clear (pp);
  mpc_clear (p);
}

mps_boolean
mps_feval_usr (mps_context * ctx, mps_polynomial * p, cplx_t x, cplx_t value, double * error)
{
  int i;
  int m = (int) (log (p->degree + 1.0) / LOG2);
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
  int m = (int) (log (p->degree + 1.0) / LOG2);
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
mps_meval_usr (mps_context * ctx, mps_polynomial * p, mpc_t x, mpc_t value, rdpe_t error)
{
  int i;
  int m = (int) (log (p->degree + 1.0) / LOG2);
  rdpe_t ax, rtmp;
  mpc_t tmp;
  long int wp = mpc_get_prec (x);

  if ((1 << m) <= p->degree)
    m++;

  mpc_init2 (tmp, wp);

  mpc_rmod (ax, x);
  mpc_set_ui (value, 1U, 0U);
  mpc_rmod (error, value);

  for (i = 1; i <= m; i++)
    {
      mpc_sqr (tmp, value);
      mpc_mul (value, x, tmp);
      mpc_add_eq_ui (value, 1U, 0U);

      rdpe_mul_eq (error, ax);
      mpc_rmod (rtmp, value);
      rdpe_add_eq (error, rtmp);
    }

  rdpe_set_2dl (rtmp, 1.0, -wp);
  rdpe_mul_eq (error, rtmp);

  mpc_clear (tmp);

  return true;
}
