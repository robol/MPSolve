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

#include <float.h>
#include <mps/gmptools.h>
#include <mps/core.h>
#include <mps/threading.h>
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
 * @param s The mps_status struct pointer.
 * @param n The degree of the polynomial.
 * @param z The value of \f$z\f$ for which the computation should be
 *  performed.
 * @param radius A pointer to the value that will be set to the
 *  computed inclusion radius.
 * @param corr Value that will be set to the newton correction,
 *  once computed.
 * @param fpc Array with the floating point coefficients of the
 *  polynomial.
 * @param fap Array with the floating points moduli of the coefficient
 *  of the polynomial.
 * @param cont mps_boolean value that will be set to true if another
 *  iteration is needed, to false otherwise.
 */
void
mps_fnewton (mps_status * s, int n, cplx_t z, double *radius, cplx_t corr,
	     cplx_t fpc[], double fap[], mps_boolean * cont, 
	     mps_boolean skip_radius_computation)
{
  int i;
  double ap, az, absp, azi, eps;
  cplx_t p, p1, zi, den, ppsp, tmp;

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
      *radius = n * (absp + eps * ap) / cplx_mod (p1);
      cplx_div (corr, p, p1);
    }
  else
    {				/* case |z|>1 */
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

	   /* If not asked to compute the radius jump to the end */
	   if (skip_radius_computation)
	     return;

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
 * @param s The mps_status struct pointer.
 * @param n The degree of the polynomial.
 * @param z The value of \f$z\f$ for which the computation should be
 *  performed.
 * @param radius A pointer to the value that will be set to the
 *  computed inclusion radius.
 * @param corr Value that will be set to the newton correction,
 *  once computed.
 * @param fpc Array with the DPE coefficients of the
 *  polynomial.
 * @param fap Array with the DPE moduli of the coefficient
 *  of the polynomial.
 * @param cont mps_boolean value that will be set to true if another
 *  iteration is needed, to false otherwise.
 *
 * @see mps_fnewton()
 */
void
mps_dnewton (mps_status * s, int n, cdpe_t z, rdpe_t radius, cdpe_t corr,
	     cdpe_t dpc[], rdpe_t dap[], mps_boolean * cont, 
	     mps_boolean skip_radius_computation)
{
  int i;
  rdpe_t ap, az, absp, rnew, apeps, rtmp;
  cdpe_t p, p1, tmp;
  double eps;

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

  /* If not asked to compute the radius jump to the end */
  if (skip_radius_computation)
    return;

  rdpe_add (rnew, absp, apeps);
  cdpe_mod (rtmp, p1);

  rdpe_div_eq (rnew, rtmp);
  if (*cont)
    rdpe_mul_d (radius, rnew, (double) n);
  else
    {
      rdpe_mul_eq_d (rnew, (double) (n + 1));
      if (rdpe_lt (rnew, radius))
	rdpe_set (radius, rnew);
    }

  rdpe_mul_d (rtmp, az, 4 * DBL_EPSILON);
  rdpe_add_eq (radius, rtmp);
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
mps_parhorner (mps_status * st, int n, mpc_t x, mpc_t p[],
	       mps_boolean b[], mpc_t s, int n_thread)
{
  int m, j, i, i1, i2, q;
  mpc_t tmp, y;
  mps_boolean bi;

  /* Set the pointer for paraller horner to be thread specific
   * so there is not conflict with other threads.           */
  mps_boolean *spar2 = mps_thread_get_spar2 (st, n_thread);
  mpc_t *mfpc2 = mps_thread_get_mfpc2 (st, n_thread);

  mpc_init2 (tmp, st->mpwp);
  mpc_init2 (y, st->mpwp);

  for (i = 0; i < n + 1; i++)
    spar2[i] = b[i];
  for (i = 0; i < n; i++)
    if (b[i])
      mpc_set (mfpc2[i], p[i]);

  q = mps_intlog2 (n + 1);
  m = n;
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
  mpc_set (s, mfpc2[0]);

  mpc_clear (y);
  mpc_clear (tmp);
}

/**
 * @brief Compute \f$p(z)\f$ by means of the parallel Horner rule.
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
mps_aparhorner (mps_status * st,
		int n, rdpe_t x, rdpe_t p[], mps_boolean b[], rdpe_t s,
		int n_thread)
{
  int m, i, j, i1, i2, q;
  rdpe_t y, tmp;
  mps_boolean bi;

  /* Set the pointer for paraller horner to be thread specific
   * so there is not conflict with other threads.           */
  mps_boolean *spar2 = mps_thread_get_spar2 (st, n_thread);
  rdpe_t *dap2 = rdpe_valloc (st->deg + 1);

  for (i = 0; i < n + 1; i++)
    spar2[i] = b[i];
  for (i = 0; i < n; i++)
    if (b[i])
      rdpe_set (dap2[i], p[i]);	/* D99 */
  q = mps_intlog2 (n + 1);
  m = n;
  rdpe_set (y, x);
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
		    rdpe_mul (tmp, y, dap2[i2]);
		    rdpe_add (dap2[i], dap2[i1], tmp);
		  }
		else
		  rdpe_set (dap2[i], dap2[i1]);
	      else
		rdpe_mul (dap2[i], y, dap2[i2]);
	    }
	  spar2[i] = bi;
	}
      spar2[m] = false;
      rdpe_sqr_eq (y);
    }
  rdpe_set (s, dap2[0]);

  rdpe_vfree (dap2);
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
 * @param s The mps_status struct pointer.
 * @param n The degree of the polynomial.
 * @param z The value of \f$z\f$ for which the computation should be
 *  performed.
 * @param radius A pointer to the value that will be set to the
 *  computed inclusion radius.
 * @param corr Value that will be set to the newton correction,
 *  once computed.
 * @param fpc Array with the DPE coefficients of the
 *  polynomial.
 * @param fap Array with the DPE moduli of the coefficient
 *  of the polynomial.
 * @param cont mps_boolean value that will be set to true if another
 *  iteration is needed, to false otherwise.
 *
 * @see mps_fnewton()
 * @see mps_dnewton()
 */
void
mps_mnewton (mps_status * s, int n, mpc_t z, rdpe_t radius, mpc_t corr,
	     mpc_t mfpc[], mpc_t mfppc[], rdpe_t dap[],
	     mps_boolean spar[], mps_boolean * cont, int n_thread, 
	     mps_boolean skip_radius_computation)
{
  int i, n1, n2;
  rdpe_t ap, az, absp, temp, rnew, ep, apeps;
  cdpe_t temp1;
  mpc_t p, p1;

  /* Set the pointer for mnewton to be thread specific
   * so there is not conflict with other threads.      */
  mps_boolean *spar2 = mps_thread_get_spar2 (s, n_thread);

  mpc_init2 (p, s->mpwp);
  mpc_init2 (p1, s->mpwp);

  rdpe_mul_d (ep, s->mp_epsilon, (double) (n * 4));
  if (s->data_type[0] == 's')
    {				/* case of sparse polynomial */
      n1 = n + 1;
      n2 = n;

      /* compute p(z) */
      mps_parhorner (s, n1, z, mfpc, spar, p, n_thread);
      mpc_get_cdpe (temp1, z);
      cdpe_mod (az, temp1);

      /* compute bound to the error */
      mps_aparhorner (s, n1, az, dap, spar, ap, n_thread);
      for (i = 0; i < n2; i++)
	spar2[i] = spar[i + 1];
      spar2[n2] = false;

      /* compute p'(z) */
      mps_parhorner (s, n2, z, mfppc, spar2, p1, n_thread);

    }
  else
    {				/*  dense polynomial */

      /* commpute p(z) and p'(z) */
      mpc_set (p, mfpc[n]);
      mpc_set (p1, p);
      for (i = n - 1; i > 0; i--)
	{
	  mpc_mul (p, p, z);
	  mpc_add (p, p, mfpc[i]);
	  mpc_mul (p1, p1, z);
	  mpc_add (p1, p1, p);
	}
      mpc_mul (p, p, z);
      mpc_add (p, p, mfpc[0]);

      /* compute bound to the error */
      rdpe_set (ap, dap[n]);
      mpc_get_cdpe (temp1, z);
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
	*cont = false;
	mpc_set_ui (corr, 0U, 0U);
	goto exit_sub;
      }
    else
      mpc_div (corr, p, p1);
  else
    {
      mpc_set_ui (corr, 0U, 0U);
      *cont = false;
      rdpe_mul (apeps, ap, ep);
      mpc_get_cdpe (temp1, p1);
      cdpe_mod (temp, temp1);
      if (rdpe_eq_zero (temp))
	{
	  if (s->DOLOG)
	    fprintf (s->logstr, "%s", "NULL DERIVATIVE\n");
	  goto exit_sub;
	}
      rdpe_div (radius, apeps, temp);
      rdpe_mul_eq_d (radius, (double) n + 1);
      goto exit_sub;
    }
  mpc_get_cdpe (temp1, p);
  cdpe_mod (absp, temp1);
  mpc_get_cdpe (temp1, p1);
  cdpe_mod (temp, temp1);
  rdpe_mul (apeps, ap, ep);
  *cont = rdpe_gt (absp, apeps);

  /* If not asked to set the radius jump to the end */
  if (skip_radius_computation)
    goto exit_sub;

  rdpe_add (rnew, absp, apeps);
  rdpe_div_eq (rnew, temp);
  if (*cont)
    rdpe_mul_d (radius, rnew, (double) n);
  else
    {
      rdpe_mul_eq_d (rnew, (double) (n + 1));
      if (rdpe_lt (rnew, radius))
	rdpe_set (radius, rnew);
    }

exit_sub:
  mpc_clear (p1);
  mpc_clear (p);
}
