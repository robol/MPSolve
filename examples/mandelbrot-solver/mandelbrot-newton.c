/*
 * This file is part of MPSolve 3.0
 *
 * Copyright (C) 2001-2013, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors: 
 *   Leonardo Robol <robol@mail.dm.unipi.it>
 */

#include "mandelbrot-poly.h"

#include <float.h>

#define LOG2 0.301029995666

/**
 * @param root The approximation whose Newton correction shall be computed. 
 * @param corr The output value where the newton correction will be stored. 
 *
 * This sample computes the 'Mandelbrot polynomial by  
 * means of the relation: p=1+x*p**2, starting with p=1
 */
void
mps_mandelbrot_poly_fnewton (mps_context * ctx, mps_polynomial * poly, mps_approximation * root, cplx_t corr)
{
  cplx_t p, pp, pt, tmp;
  double ap, ax, eps;
  int i, m, n = poly->degree;

  cplx_t x;
  mps_approximation_get_fvalue (ctx, root, x);
  
  double rad = mps_approximation_get_frad (ctx, root);
  mps_boolean again = mps_approximation_get_again (ctx, root);

  m = (int) (log (n + 1.0) / LOG2);
  if ((1 << m) <= n)
    m++;
  eps = (DBL_EPSILON * 4.0) * n;
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
  
  mps_approximation_set_again (ctx, root, cplx_mod (p) > eps * ap * 3);
  mps_approximation_set_frad (ctx, root, n * (cplx_mod (p) + 3 * ap * eps) / cplx_mod (pp));
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
void
mps_mandelbrot_poly_dnewton (mps_context * ctx, mps_polynomial * poly, mps_approximation * root, cdpe_t corr)
{
  cdpe_t p, pp, pt, tmp;
  rdpe_t ap, ax, eps, temp, apeps, atmp;
  int i, m, n = poly->degree;
  cdpe_t x;
  rdpe_t drad;
  mps_boolean again;

  mps_approximation_get_dvalue (ctx, root, x);
  mps_approximation_get_drad (ctx, root, drad);
  again = mps_approximation_get_again (ctx, root);

  m = (int) (log (n + 1.0) / LOG2);
  if ((1 << m) <= n)
    m++;
  rdpe_set_d (eps, DBL_EPSILON);
  rdpe_mul_eq_d (eps, (double) 4 * n);
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

  mps_approximation_set_again (ctx, root, rdpe_gt (temp, apeps));

  rdpe_add (drad, temp, apeps);
  rdpe_mul_eq_d (drad, (double) n);
  cdpe_mod (temp, pp);
  rdpe_div_eq (drad, temp);
  if (rdpe_eq (drad, rdpe_zero))
    rdpe_mul (drad, ax, eps);

  mps_approximation_set_drad (ctx, root, drad);
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
void 
mps_mandelbrot_poly_mnewton (mps_context * ctx, mps_polynomial * poly, 
			     mps_approximation * root, mpc_t corr, long int wp)
{
  int i, m, n = poly->degree;
  rdpe_t ap, ax, eps, temp, apeps, atmp, epsilon, drad;
  cdpe_t ctmp;
  mpc_t p, pp, pt, tmp, x;
  mps_boolean again;

  mpc_init2 (p, wp);
  mpc_init2 (pp, wp);
  mpc_init2 (pt, wp);
  mpc_init2 (tmp, wp);

  mpc_init2 (x, wp);
  mps_approximation_get_mvalue (ctx, root, x);
  mps_approximation_get_drad (ctx, root, drad);
  again = mps_approximation_get_again (ctx, root);
  
  rdpe_set_2dl (epsilon, 1.0, 2 - wp);

  m = (int) (log (n + 1.0) / LOG2);
  if ((1 << m) <= n)
    m++;
  rdpe_set (eps, epsilon);
  rdpe_mul_eq_d (eps, (double) 4 * n);
  mpc_get_cdpe (ctmp, x);
  cdpe_mod (ax, ctmp);

  mpc_set_ui (p, 1, 0);
  mpc_set_ui (pp, 0, 0);
  rdpe_set (ap, rdpe_one);
  for (i = 1; i <= m; i++)
    {
      mpc_sqr (tmp, p);
      mpc_mul (pt, x, tmp);
      mpc_add_eq_ui (pt, 1, 0);
      mpc_mul_eq (pp, x);
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
  mps_approximation_set_again (ctx, root, rdpe_gt (temp, apeps));

  rdpe_add (drad, temp, apeps);
  rdpe_mul_eq_d (drad, (double) n);
  mpc_get_cdpe (ctmp, pp);
  cdpe_mod (temp, ctmp);
  rdpe_div_eq (drad, temp);
  if (rdpe_eq (drad, rdpe_zero))
    rdpe_mul (drad, ax, eps);

  mps_approximation_set_drad (ctx, root, drad);
  mps_approximation_set_again (ctx, root, again);

  mpc_clear (tmp);
  mpc_clear (pt);
  mpc_clear (pp);
  mpc_clear (p);
  mpc_clear (x);
}
