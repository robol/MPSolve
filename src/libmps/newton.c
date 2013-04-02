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


#include <float.h>
#include <mps/mps.h>
#include <time.h>
#include <math.h>

/**
 * @brief Compute the Newton correction, i.e. and the value \f$s\f$
 * given by:
 * \f[
 *      nwt = \frac{p(z)}{p'(z)} \qquad s =
 *      \left\lvert \frac{p(\lvert z \rvert)}{p'(\lvert z \rvert)} \right\rvert
 * \f]
 * and set the parameter <code>cont</code> to <code>true</code> if
 * the newton correction is greater than \f$4n\epsilon s\f$ and to
 * false otherwise.
 *
 * @param s The mps_context struct pointer.
 * @param poly The polynomial to evaluate, casted to a mps_polynomial.
 * @param root The approximation where the newton fraction should be evaluated.
 * @param corr The complex value of the newton correction. 
 */
void
mps_fnewton (mps_context * s, mps_polynomial * poly, mps_approximation * root, 
             cplx_t corr)
{
  int i;
  double ap, az, absp, azi, eps;
  cplx_t p, p1, zi, den, ppsp, tmp;

  mps_monomial_poly *mp = MPS_MONOMIAL_POLY (poly);
  cplx_t *fpc = mp->fpc;
  double *fap = mp->fap;
  int n = poly->degree;

  cplx_t z;
  cplx_set (z, root->fvalue);
  double * radius = &root->frad;
  mps_boolean * cont = &root->again;

  eps = 4 * n * DBL_EPSILON;
  az = cplx_mod (z);

  /* distinguish the cases |z|<=1, |z|>1 */
  if (az <= 1)
    {
      /*  case |z|<=1 */
      cplx_set (p, fpc[n]);
      cplx_set (p1, p);
      for (i = n - 1; i > 0; i--)
        {
          cplx_mul (tmp, p, z);
          cplx_add (p, tmp, fpc[i]);
          cplx_mul (tmp, p1, z);
          cplx_add (p1, tmp, p);
        }
      cplx_mul (tmp, p, z);
      cplx_add (p, tmp, fpc[0]);
      ap = fap[n];
      for (i = n - 1; i >= 0; i--)
        ap = ap * az + fap[i];
      absp = cplx_mod (p);
      *cont = (absp > ap * eps);
      *radius = n * (absp + eps * ap) / cplx_mod (p1) + DBL_MIN;
      cplx_div (corr, p, p1);
    }
  else
    {                           /* case |z|>1 */
      cplx_set (zi, z);
      cplx_inv_eq (zi);
      azi = 1.0 / az;
      cplx_set (p, fpc[0]);
      cplx_set (p1, p);
      for (i = n - 1; i > 0; i--)
        {
          cplx_mul (tmp, p, zi);
          cplx_add (p, tmp, fpc[n - i]);
          cplx_mul (tmp, p1, zi);
          cplx_add (p1, tmp, p);
        }
      cplx_mul (tmp, p, zi);
      cplx_add (p, tmp, fpc[n]);

      ap = fap[0];
      for (i = 1; i <= n; i++)
        ap = ap * azi + fap[i];
      absp = cplx_mod (p);
      *cont = (absp > ap * eps);

      cplx_mul_d (den, p, (double) n);
      cplx_mul (ppsp, p1, zi);
      cplx_sub_eq (den, ppsp);
      cplx_mul_eq (den, zi);
      if (cplx_mod (den) != 0)
        {
          cplx_div (corr, p, den);
          ap = (ap * eps + absp) * n;
          ap = ap / cplx_mod (den);
          *radius = ap;
        }
      else
        {
          cplx_mul (ppsp, p, z);
          cplx_div_eq (ppsp, p1);
          cplx_mul_d (den, ppsp, (double) n);
          cplx_sub_eq (den, cplx_one);
          cplx_div (corr, ppsp, den);
          cplx_mul_eq (corr, z);
          absp = cplx_mod (p);
          *cont = (absp > ap * eps);

          *radius = cplx_mod (ppsp) + (eps * ap * az) / cplx_mod (p1);
          *radius *= n / cplx_mod (den);
          *radius *= az;
        }

    }
}


/**
 * @brief Compute the Newton correction, i.e. and the value \f$s\f$
 * given by:
 * \f[
 *      nwt = \frac{p(z)}{p'(z)} \qquad s =
 *      \left\lvert \frac{p(\lvert z \rvert)}{p'(\lvert z \rvert)} \right\rvert
 * \f]
 * and set the parameter <code>cont</code> to <code>true</code> if
 * the newton correction is greater than \f$4n\epsilon s\f$ and to
 * false otherwise.
 *
 * This routine is the DPE version of <code>mps_fnewton()</code>.
 *
 *
 * @param s The mps_context struct pointer.
 * @param poly The polynomial to evaluate, casted to a mps_polynomial.
 * @param root The approximation where the newton fraction should be evaluated.
 * @param corr The complex value of the newton correction. 
 *
 *
 * @see mps_fnewton()
 */
void
mps_dnewton (mps_context * s, mps_polynomial * poly, mps_approximation * root, 
             cdpe_t corr)
{
  int i;
  rdpe_t ap, az, absp, rnew, apeps, rtmp;
  cdpe_t p, p1, tmp, z;
  double eps;

  mps_monomial_poly *mp = MPS_MONOMIAL_POLY (poly);
  cdpe_t * dpc = mp->dpc;
  rdpe_t * dap = mp->dap;
  int n = poly->degree;

  cdpe_set (z, root->dvalue);
  mps_boolean * cont = &root->again;

  eps = DBL_EPSILON * n * 4;
  cdpe_set (p, dpc[n]);
  cdpe_set (p1, p);
  for (i = n - 1; i > 0; i--)
    {
      cdpe_mul (tmp, p, z);
      cdpe_add (p, tmp, dpc[i]);
      cdpe_mul (tmp, p1, z);
      cdpe_add (p1, tmp, p);
    }
  cdpe_mul (tmp, p, z);
  cdpe_add (p, tmp, dpc[0]);
  if (cdpe_ne (p, cdpe_zero))
    if (cdpe_eq (p1, cdpe_zero))
      {
        if (s->DOLOG)
          fprintf (s->logstr, "%s", "NULL DERIVATIVE\n");
        cdpe_set (corr, cdpe_zero);
        *cont = false;
        return;
      }
    else
      cdpe_div (corr, p, p1);
  else
    {
      cdpe_set (corr, cdpe_zero);
      *cont = false;
    }
  cdpe_mod (az, z);
  rdpe_set (ap, dap[n]);
  for (i = n - 1; i >= 0; i--)
    {
      rdpe_mul (rtmp, ap, az);
      rdpe_add (ap, rtmp, dap[i]);
    }
  cdpe_mod (absp, p);
  rdpe_mul_d (apeps, ap, eps);
  *cont = rdpe_gt (absp, apeps);

  rdpe_add (rnew, absp, apeps);
  cdpe_mod (rtmp, p1);

  rdpe_div_eq (rnew, rtmp);
  if (*cont)
    rdpe_mul_d (root->drad, rnew, (double) n);
  else
    {
      rdpe_mul_eq_d (rnew, (double) (n + 1));
      if (rdpe_lt (rnew, root->drad))
        rdpe_set (root->drad, rnew);
    }

  rdpe_mul_d (rtmp, az, 4 * DBL_EPSILON);
  rdpe_add_eq (root->drad, rtmp);
}

/****************************************************
*          FUNCTION INTLOG2                         *
****************************************************/
int
mps_intlog2 (int n)
{
  int k;
  k = (int) (log (n) / LOG2);
  if (1 << k < n)
    k++;
  return k;
}

/**
 * @brief Compute the Newton correction, i.e. and the value \f$s\f$
 * given by:
 * \f[
 *      nwt = \frac{p(z)}{p'(z)} \qquad s =
 *      \left\lvert \frac{p(\lvert z \rvert)}{p'(\lvert z \rvert)} \right\rvert
 * \f]
 * and set the parameter <code>cont</code> to <code>true</code> if
 * the newton correction is greater than \f$4n\epsilon s\f$ and to
 * false otherwise.
 *
 * This routine is the multiprecision version of <code>mps_fnewton()</code>.
 * It differs from these routines, indeed, because it uses a more
 * sophisticated tecnique to perform the computation:
 * - If the polynomial are dense the computation is performed with
 *   the Horner's rule
 * - If the polynomial is sparse then the computation is performed
 *   with the parallel Horner algorithm.
 *
 * @param s The mps_context struct pointer.
 * @param poly The polynomial to evaluate, casted to a mps_polynomial.
 * @param root The approximation where the newton fraction should be evaluated.
 * @param corr The complex value of the newton correction. 
 *
 * @see mps_fnewton()
 * @see mps_dnewton()
 */
void
mps_mnewton (mps_context * s, mps_polynomial * poly, 
             mps_approximation * root, mpc_t corr)
{
  int i;
  rdpe_t ap, az, absp, temp, rnew, ep, apeps;
  cdpe_t temp1;
  mpc_t p, p1;

  mps_monomial_poly * mp = MPS_MONOMIAL_POLY (poly);
  mpc_t * mfpc = mp->mfpc;
  rdpe_t * dap = mp->dap;
  int n = poly->degree;

  long int wp = mpc_get_prec (root->mvalue);

  mpc_init2 (p, wp);
  mpc_init2 (p1, wp);

  rdpe_set_2dl (ep, 1.0, 2 - wp);
  rdpe_mul_eq_d (ep, n);
  
  if (MPS_DENSITY_IS_SPARSE (poly->density))
    {
      /* That's a dirty trick to setup a hackish derivative that 
       * points to the same internal fields of the polynomial, shifted
       * by one. */
      mps_monomial_poly derivative;
      mps_polynomial_init (s, MPS_POLYNOMIAL (&derivative));

      MPS_POLYNOMIAL (&derivative)->degree = MPS_POLYNOMIAL (mp)->degree - 1;
      derivative.spar = mp->spar + 1;
      derivative.prec = mp->prec;
      derivative.mfpc_mutex = mp->mfpc_mutex + 1;

      derivative.mfpc = mpc_valloc (n);
      mpc_vinit2 (derivative.mfpc, n, wp);
      for (i = 0; i < n; i++)
        mpc_mul_ui (derivative.mfpc[i], mp->mfpc[i+1], i+1);

      MPS_POLYNOMIAL (&derivative)->meval = mps_monomial_poly_meval;
      MPS_POLYNOMIAL (&derivative)->raise_data = mps_monomial_poly_raise_precision;

      mps_polynomial_meval (s, MPS_POLYNOMIAL (mp), root->mvalue, p, ap);
      mps_mhorner (s, &derivative, root->mvalue, p1);

      mpc_vclear (derivative.mfpc, n);
      mpc_vfree (derivative.mfpc);
    }
  else
    {                           /*  dense polynomial */
      /* commpute p(z) and p'(z) */
      mpc_set (p, mfpc[n]);
      mpc_set (p1, p);
      for (i = n - 1; i > 0; i--)
        {
          mpc_mul (p, p, root->mvalue);
          mpc_add (p, p, mfpc[i]);
          mpc_mul (p1, p1, root->mvalue);
          mpc_add (p1, p1, p);
        }
      mpc_mul (p, p, root->mvalue);
      mpc_add (p, p, mfpc[0]);

      /* compute bound to the error */
      rdpe_set (ap, dap[n]);
      mpc_get_cdpe (temp1, root->mvalue);
      cdpe_mod (az, temp1);
      for (i = n - 1; i >= 0; i--)
        {
          rdpe_mul (temp, ap, az);
          rdpe_add (ap, temp, dap[i]);
        }
    }

  /* common part */
  if (!mpc_eq_zero (p))
    if (mpc_eq_zero (p1))
      {
        if (s->DOLOG)
          fprintf (s->logstr, "%s", "NULL DERIVATIVE\n");
        root->again = false;
        mpc_set_ui (corr, 0U, 0U);
        goto exit_sub;
      }
    else
      mpc_div (corr, p, p1);
  else
    {
      mpc_set_ui (corr, 0U, 0U);
      root->again = false;
      rdpe_mul (apeps, ap, ep);
      mpc_get_cdpe (temp1, p1);
      cdpe_mod (temp, temp1);
      if (rdpe_eq_zero (temp))
        {
          if (s->DOLOG)
            fprintf (s->logstr, "%s", "NULL DERIVATIVE\n");
          goto exit_sub;
        }
      rdpe_div (root->drad, apeps, temp);
      rdpe_mul_eq_d (root->drad, (double) n + 1);
      goto exit_sub;
    }
  mpc_get_cdpe (temp1, p);
  cdpe_mod (absp, temp1);
  mpc_get_cdpe (temp1, p1);
  cdpe_mod (temp, temp1);
  rdpe_mul (apeps, ap, ep);
  root->again = rdpe_gt (absp, apeps);

  rdpe_add (rnew, absp, apeps);
  rdpe_div_eq (rnew, temp);
  if (root->again)
    rdpe_mul_d (root->drad, rnew, (double) n);
  else
    {
      rdpe_mul_d (root->drad, rnew, (double) (n + 1));
    }

exit_sub:
  mpc_clear (p1);
  mpc_clear (p);
}
