/***********************************************************
**       Multiprecision Polynomial Solver (MPSolve)       **
**                 Version 2.2, May 2001                  **
**                                                        **
**                      Written by                        **
**       Dario Andrea Bini and Giuseppe Fiorentino        **
**       (bini@dm.unipi.it)  (fiorent@dm.unipi.it)        **
**                                                        **
** (C) 2001, Dipartimento di Matematica, FRISCO LTR 21024 **
***********************************************************/

#include "mps.h"

/**********************************************************
*                     SUBROUTINE FNEWTON                  *
***********************************************************
* Compute the Newton correction c=p(z)/p'(z) together with*
* the value s=|p|(|z|)/|p'(z)|,                           *
* and set cont=true if |c|>4*n*EPS*s and false, otherwise.*
* Use different formulae according to the cases           *
* |z|>1 and |z|<=1 in order to avoid overflow.            *
* Real coefficients: to be implemented                    *
**********************************************************/
void
fnewton(int n, cplx_t z, double *radius, cplx_t corr,
	cplx_t fpc[], double fap[], boolean * cont)
{
  int i;
  double ap, az, absp, azi, eps;
  cplx_t p, p1, zi, den, ppsp, tmp;

  eps = 4 * n * DBL_EPSILON;
  az = cplx_mod(z);

  /* distinguish the cases |z|<=1, |z|>1 */

  if (az <= 1) {
    /*  case |z|<=1 */
    cplx_set(p, fpc[n]);
    cplx_set(p1, p);
    for (i = n - 1; i > 0; i--) {
      cplx_mul(tmp, p, z);
      cplx_add(p, tmp, fpc[i]);
      cplx_mul(tmp, p1, z);
      cplx_add(p1, tmp, p);
    }
    cplx_mul(tmp, p, z);
    cplx_add(p, tmp, fpc[0]);
    ap = fap[n];
    for (i = n - 1; i >= 0; i--)
      ap = ap * az + fap[i];
    absp = cplx_mod(p);
    *cont = (absp > ap * eps);
    *radius = n * (absp + eps * ap) / cplx_mod(p1);
    cplx_div(corr, p, p1);
  } else {			/* case |z|>1 */
    cplx_set(zi, z);
    cplx_inv_eq(zi);
    azi = 1.0 / az;
    cplx_set(p, fpc[0]);
    cplx_set(p1, p);
    for (i = n - 1; i > 0; i--) {
      cplx_mul(tmp, p, zi);
      cplx_add(p, tmp, fpc[n - i]);
      cplx_mul(tmp, p1, zi);
      cplx_add(p1, tmp, p);
    }
    cplx_mul(tmp, p, zi);
    cplx_add(p, tmp, fpc[n]);

    ap = fap[0];
    for (i = 1; i <= n; i++)
      ap = ap * azi + fap[i];
    absp = cplx_mod(p);
    *cont = (absp > ap * eps);

    cplx_mul_d(den,p,(double) n);
    cplx_mul(ppsp,p1,zi);
    cplx_sub_eq(den,ppsp);
    cplx_mul_eq(den,zi);
    if(cplx_mod(den) != 0) {
      cplx_div(corr,p,den);
      ap=(ap*eps+absp)*n;
      ap=ap/cplx_mod(den);
      *radius = ap;
    }
    else
      { 
	cplx_mul(ppsp, p, z);
	cplx_div_eq(ppsp, p1);
	cplx_mul_d(den, ppsp, (double) n);
	cplx_sub_eq(den, cplx_one);
	cplx_div(corr, ppsp, den);
	cplx_mul_eq(corr, z);
	absp = cplx_mod(p);
	*cont = (absp > ap * eps);
	*radius = cplx_mod(ppsp) + (eps * ap * az) / cplx_mod(p1);
	*radius *= n / cplx_mod(den);
	*radius *= az;
      }

  }
}

/**********************************************************
*                     SUBROUTINE DNEWTON                  *
***********************************************************
* Compute the Newton correction c=p(z)/p'(z) together with*
* the value s=|p|(|z|)/|p'(z)|,                           *
* and set cont=.true. if |c|>4*n*EPS*s .false., otherwise.*
* Real coefficients: to be implemented                    *
**********************************************************/
void
dnewton(int n, cdpe_t z, rdpe_t radius, cdpe_t corr,
	cdpe_t dpc[], rdpe_t dap[], boolean * cont)
{
  int i;
  rdpe_t ap, az, absp, rnew, apeps, rtmp;
  cdpe_t p, p1, tmp;
  double eps;

  eps = DBL_EPSILON * n * 4;
  cdpe_set(p, dpc[n]);
  cdpe_set(p1, p);
  for (i = n - 1; i > 0; i--) {
    cdpe_mul(tmp, p, z);
    cdpe_add(p, tmp, dpc[i]);
    cdpe_mul(tmp, p1, z);
    cdpe_add(p1, tmp, p);
  }
  cdpe_mul(tmp, p, z);
  cdpe_add(p, tmp, dpc[0]);
  if (cdpe_ne(p, cdpe_zero))
    if (cdpe_eq(p1, cdpe_zero)) {
      if (DOLOG)
	fprintf(logstr, "%s", "NULL DERIVATIVE\n");
      cdpe_set(corr, cdpe_zero);
      *cont = false;
      return;
    } else
      cdpe_div(corr, p, p1);
  else {
    cdpe_set(corr, cdpe_zero);
    *cont = false;
  }
  cdpe_mod(az, z);
  rdpe_set(ap, dap[n]);
  for (i = n - 1; i >= 0; i--) {
    rdpe_mul(rtmp, ap, az);
    rdpe_add(ap, rtmp, dap[i]);
  }
  cdpe_mod(absp, p);
  rdpe_mul_d(apeps, ap, eps);
  *cont = rdpe_gt(absp, apeps);
  rdpe_add(rnew, absp, apeps);
  cdpe_mod(rtmp, p1);

  rdpe_div_eq(rnew, rtmp);
  if (*cont)
    rdpe_mul_d(radius, rnew, (double) n);
  else {
    rdpe_mul_eq_d(rnew, (double) (n + 1));
    if (rdpe_lt(rnew, radius))
      rdpe_set(radius, rnew);
  }
}

/****************************************************
*          FUNCTION INTLOG2                         *
****************************************************/
int
intlog2(int n)
{
  int k;
  k = (int) (log(n) / LOG2);
  if (1 << k < n)
    k++;
  return k;
}

/****************************************************
*              SUBROUTINE PARHORNER                 *   
***************************************************** 
* Computes p(z) by means of the parallel Horner     *
* algorithm  b(i)= .false. means p(i)=0             *
****************************************************/
void
parhorner(int n, mpc_t x, mpc_t p[], boolean b[],
	  mpc_t s)
{
  int m, j, i, i1, i2, q;
  tmpc_t tmp, y;
  boolean bi;

  tmpc_init2(tmp, mpwp);
  tmpc_init2(y, mpwp);

  for (i = 0; i < n + 1; i++)
    spar2[i] = b[i];
  for (i = 0; i < n; i++)
    if(b[i]) 
      mpc_set(mfpc2[i], p[i]);

  q = intlog2(n + 1);
  m = n;
  mpc_set(y, x);
  for (j = 0; j < q; j++) {
    spar2[m] = false;
    m = (m + 1) >> 1;
    for (i = 0; i < m; i++) {
      i2 = (i << 1) + 1;
      i1 = i2 - 1;
      bi = spar2[i1] || spar2[i2];
      if (bi) {
	if (spar2[i1])
	  if (spar2[i2]) {
	    mpc_mul(tmp, y, mfpc2[i2]);
	    mpc_add(mfpc2[i], mfpc2[i1], tmp);
	  } else
	    mpc_set(mfpc2[i], mfpc2[i1]);
	else
	  mpc_mul(mfpc2[i], y, mfpc2[i2]);
      }
      spar2[i] = bi;
    }
    spar2[m] = false;
    mpc_sqr_eq(y);
  }
  mpc_set(s, mfpc2[0]);

  tmpc_clear(y);
  tmpc_clear(tmp);
}

/****************************************************
*              SUBROUTINE APARHORNER                *   
***************************************************** 
* Computes |p|(|z|) by means of the parallel Horner *
* algorithm  b(i)= .false. means p(i)=0             *
****************************************************/
void
aparhorner(int n, rdpe_t x, rdpe_t p[], boolean b[], rdpe_t s)
{
  int m, i, j, i1, i2, q;
  rdpe_t y, tmp;
  boolean bi;

  for (i = 0; i < n + 1; i++)
    spar2[i] = b[i];
  for (i = 0; i < n; i++)
   if(b[i]) rdpe_set(dap2[i], p[i]); /* D99 */
  q = intlog2(n + 1);
  m = n;
  rdpe_set(y, x);
  for (j = 0; j < q; j++) {
    spar2[m] = false;
    m = (m + 1) >> 1;
    for (i = 0; i < m; i++) {
      i2 = (i << 1) + 1;
      i1 = i2 - 1;
      bi = spar2[i1] || spar2[i2];
      if (bi) {
	if (spar2[i1])
	  if (spar2[i2]) {
	    rdpe_mul(tmp, y, dap2[i2]);
	    rdpe_add(dap2[i], dap2[i1], tmp);
	  } else
	    rdpe_set(dap2[i], dap2[i1]);
	else
	  rdpe_mul(dap2[i], y, dap2[i2]);
      }
      spar2[i] = bi;
    }
    spar2[m] = false;
    rdpe_sqr_eq(y);
  }
  rdpe_set(s, dap2[0]);
}

/******************************************************
*              SUBROUTINE MNEWTON                     *   
******************************************************* 
* Compute the Newton correction c=p(z)/p'(z) together *
* with the value s=|p|(|z|)/|p'(z)|, and set cont     *
* .true. if |c|>4*n*EPS*s, .false. otherwise,         *
* where EPS is the current epsilon-precision.         *
* For dense polynomials the computation of the values *
* p(z) and p'(z) is  performed with Horner's method,  *
* if the polynomial is sparse then the 'parallel'     *
* Horner's  algorithm is used.                        *
* Real coefficients: to be implemented                *
******************************************************/
void
mnewton(int n, mpc_t z, rdpe_t radius, mpc_t corr,
	mpc_t mfpc[], mpc_t mfppc[], rdpe_t dap[],
	boolean spar[], boolean * cont)
{
  int i, n1, n2;
  rdpe_t ap, az, absp, temp, rnew, ep, apeps;
  cdpe_t temp1;
  tmpc_t p, p1;

  tmpc_init2(p, mpwp);
  tmpc_init2(p1, mpwp);

  rdpe_mul_d(ep, mp_epsilon, (double) (n * 4));
  if (data_type[0] == 's') {	/* case of sparse polynomial */
    n1 = n + 1;
    n2 = n;

    /* compute p(z) */
    parhorner(n1, z, mfpc, spar, p);
    mpc_get_cdpe(temp1, z);
    cdpe_mod(az, temp1);

    /* compute bound to the error */
    aparhorner(n1, az, dap, spar, ap);
    for (i = 0; i < n2; i++)
      spar2[i] = spar[i + 1];
    spar2[n2] = false;

    /* compute p'(z) */
    parhorner(n2, z, mfppc, spar2, p1);

  } else {			/*  dense polynomial */

    /* commpute p(z) and p'(z) */
    mpc_set(p, mfpc[n]);
    mpc_set(p1, p);
    for (i = n - 1; i > 0; i--) {
      mpc_mul(p, p, z);
      mpc_add(p, p, mfpc[i]);
      mpc_mul(p1, p1, z);
      mpc_add(p1, p1, p);
    }
    mpc_mul(p, p, z);
    mpc_add(p, p, mfpc[0]);

    /* compute bound to the error */
    rdpe_set(ap, dap[n]);
    mpc_get_cdpe(temp1, z);
    cdpe_mod(az, temp1);
    for (i = n - 1; i >= 0; i--) {
      rdpe_mul(temp, ap, az);
      rdpe_add(ap, temp, dap[i]);
    }
  }

  /* common part */
  if (!mpc_eq_zero(p))
    if (mpc_eq_zero(p1)) {
      if (DOLOG)
	fprintf(logstr, "%s", "NULL DERIVATIVE\n");
      *cont = false;
      mpc_set_ui(corr, 0U, 0U);
      goto exit_sub;
    } else
      mpc_div(corr, p, p1);
  else {
    mpc_set_ui(corr, 0U, 0U);
    *cont = false;
    rdpe_mul(apeps, ap, ep);
    mpc_get_cdpe(temp1, p1);
    cdpe_mod(temp, temp1);
    if (rdpe_eq_zero(temp)) {
      fprintf(logstr, "%s", "NULL DERIVATIVE\n");
      goto exit_sub;
    }
    rdpe_div(radius, apeps, temp);
    rdpe_mul_eq_d(radius, (double) n + 1);
    goto exit_sub;
  }
  mpc_get_cdpe(temp1, p);
  cdpe_mod(absp, temp1);
  mpc_get_cdpe(temp1, p1);
  cdpe_mod(temp, temp1);
  rdpe_mul(apeps, ap, ep);
  *cont = rdpe_gt(absp, apeps);
  rdpe_add(rnew, absp, apeps);
  rdpe_div_eq(rnew, temp);
  if (*cont)
    rdpe_mul_d(radius, rnew, (double) n);
  else {
    rdpe_mul_eq_d(rnew, (double) (n + 1));
    if (rdpe_lt(rnew, radius))
      rdpe_set(radius, rnew);
  }

exit_sub:
  tmpc_clear(p1);
  tmpc_clear(p);
}
