/***********************************************************
**       Multiprecision Polynomial Solver (MPSolve)       **
**              Version 2.1, september 1999               **
**                                                        **
**                      Written by                        **
**       Dario Andrea Bini and Giuseppe Fiorentino        **
**       (bini@dm.unipi.it)  (fiorent@dm.unipi.it)        **
**                                                        **
** (C) 1999, Dipartimento di Matematica, FRISCO LTR 21024 **
***********************************************************/

/********************************************************
This file contains a sample of a user-defined polynomial.
In this sample we consider a class of polynomials suggested 
by Alan Sokal defined by means of the relation 
  p_{i+1}(x)=p_i(x)^r+x*p_{i-1}(x)^{r^2},
  p_0(x)=1,                                   (1)
  p_1(x)=1 
and for the derivatives by 
  p'_{i+1}(x)=r*p_i(x)^{r-1}*p_i'+p_{i-1}(x)^{r^2}+
              x*r^2*p_{i-1}(x)^{r^2-1}*p'_{i-1}(x).
The allowed degrees of the polynomials in this class are
0,1,2,5,10,21,42,85,170,...
The three programs below compute the values of p_i(x) and
p'_i(x) for r=2, by means of the above formulae in the float
version, dpe version and mp version. 
The programs compute also an error bound needed
for testing the stop condition,  and for computing
the radius of the inclusion disk.
The formulae for the latter computation are
obtained by means of a rounding error analysis of (1).
**********************************************************/

#include "mps.h"

/******************************************************
*              SUBROUTINE FNEWTON_USR                 *   
******************************************************* 
* Floating point computation                          *
******************************************************/
void
fnewton_usr(cplx_t x, double *rad, cplx_t corr, mps_boolean * again)
{
 /* floating point user polynomial not supported yet */
}

/******************************************************
*              SUBROUTINE DNEWTON_USR                 *   
******************************************************* 
* DPE computation                                     *
******************************************************/
void
dnewton_usr(cdpe_t x, rdpe_t rad, cdpe_t corr, mps_boolean * again)
{
  cdpe_t p0, p1, p2, pp0, pp1, pp2, tmp1, tmp2;
  rdpe_t d0, d1, d2, ap0, ap1, ap2, ax, eps, rtmp, rtmp1;
  int i, m;

  rdpe_set_d(eps, DBL_EPSILON);
  m = (int) (log(n)/LOG2) + 1;
  cdpe_mod(ax, x);

  /* initial conditions */
  cdpe_set(p0, cdpe_one);    /* p0=1   */
  cdpe_set(p1, cdpe_one);    /* p1=1   */
  cdpe_set(pp0, cdpe_zero);  /* p0'=0  */
  cdpe_set(pp1, cdpe_zero);  /* p1'=0  */

  rdpe_set(ap0, rdpe_one);   /* |p0|=1 */
  rdpe_set(ap1, rdpe_one);   /* |p1|=1 */
  rdpe_set(d0, rdpe_zero);
  rdpe_set(d1, d0);

  /* computation */
  for (i = 1; i <= m; i++) {
    /* polynomial */
    cdpe_mul(tmp2, p0, p0);
    cdpe_mul(pp2, tmp2, tmp2);
    cdpe_mul(p2, pp2, x);
    cdpe_mul(tmp1, p1, p1);
    cdpe_add(p2, p2, tmp1);

    /* derivative  */
    cdpe_mul(tmp1, p0, tmp2);
    cdpe_mod(rtmp, tmp1);
    cdpe_mul(tmp1, tmp1, pp0);
    cdpe_mul(tmp1, x, tmp1);

    cdpe_mul_2exp(tmp1, tmp1,2);

    cdpe_add(pp2, tmp1, pp2);
    cdpe_mul(tmp1, p1, pp1);
    cdpe_mul_2exp(tmp1, tmp1,1);
 
    cdpe_add(pp2, tmp1, pp2);
    cdpe_set(p0, p1);
    cdpe_set(p1, p2);
    cdpe_set(pp0, pp1);
    cdpe_set(pp1, pp2);

   /* error bound */
   /* d2=2 d1 |p1| + 4 |x| d0 |p0|^3 + 
         u( |p2| + 2sqrt 2 |p1|^2 + |x| |p0|^4(4sqrt 2+1) */
    cdpe_mod(ap2, p2);
    rdpe_mul(rtmp, rtmp, d0);
    rdpe_mul(rtmp, rtmp, ax);

    rdpe_mul_2exp(rtmp, rtmp,2);  /* 4 |x| |p0|^3 d0 */

    rdpe_mul(d2, d1, ap1);

    rdpe_mul_2exp(d2, d2,1);  /*  2 d1 |p1|  */
    rdpe_add(d2, d2, rtmp);

    rdpe_mul(rtmp1, ap0, ap0);
    rdpe_mul(rtmp1, rtmp1, rtmp1); /*  |p0|^4 */
    rdpe_mul(rtmp1, rtmp1, ax);
    rdpe_mul_d(rtmp1, rtmp1, (double) 6.66); /* |x|(4\sqrt 2 +1)|p0|^4 */
    rdpe_add(rtmp1, rtmp1, ap2);
    rdpe_mul(rtmp, ap1, ap1);
    rdpe_mul_d(rtmp, rtmp, (double) 2.83);
    rdpe_mul_2exp(rtmp, rtmp,1);

    rdpe_mul(rtmp, rtmp, eps);
    rdpe_add(d2, d2, rtmp);

    rdpe_set(d0, d1);
    rdpe_set(d1, d2);
    rdpe_set(ap0, ap1);
    rdpe_set(ap1, ap2);

  }

  cdpe_div(corr, p1, pp1);
  *again = rdpe_lt(d2, ap2);
  cdpe_mod(rtmp, pp1);
  rdpe_add(rad, ap2, d2);
  rdpe_div_eq(rad,  rtmp);
  rdpe_mul_eq_d(rad, (double) n);
}

/******************************************************
*              SUBROUTINE MNEWTON_USR                 *   
******************************************************* 
 multiprecision computation
******************************************************/
void
mnewton_usr(mpc_t x, rdpe_t rad, mpc_t corr, mps_boolean * again)
{
  rdpe_t d0, d1, d2, ap0, ap1, ap2, ax, eps, rtmp, rtmp1;
  int i, m;
  cdpe_t ctmp;
  tmpc_t p0, p1, p2, pp0, pp1, pp2, tmp1, tmp2;

  tmpc_init2(p0, mpwp);
  tmpc_init2(p1, mpwp);
  tmpc_init2(p2, mpwp);
  tmpc_init2(pp0, mpwp);
  tmpc_init2(pp1, mpwp);
  tmpc_init2(pp2, mpwp);
  tmpc_init2(tmp1, mpwp);
  tmpc_init2(tmp2, mpwp);

  rdpe_set(eps, mp_epsilon);
  m = (int) (log(n ) / LOG2);
  m=m+1;
  mpc_get_cdpe(ctmp, x);
  cdpe_mod(ax, ctmp);

  /* initial conditions */
  mpc_set_d(p0, 1.0, 0.0);    /* p0=1   */
  mpc_set_d(p1, 1.0, 0.0);    /* p1=1   */
  mpc_set_d(pp0, 0.0, 0.0);   /* p0'=0  */
  mpc_set_d(pp1, 0.0, 0.0);   /* p1'=0  */

  rdpe_set(ap0, rdpe_one);    /* |p0|=1 */
  rdpe_set(ap1, rdpe_one);    /* |p1|=1 */
  rdpe_set(d0, rdpe_zero);
  rdpe_set(d1, d0);

  /* computation */
  for (i = 1; i <= m; i++) {
    /* polynomial */
    mpc_mul(tmp2, p0, p0);
    mpc_mul(pp2, tmp2, tmp2);
    mpc_mul(p2, pp2, x);
    mpc_mul(tmp1, p1, p1);
    mpc_add(p2, p2, tmp1);

   /* derivative  */
    mpc_mul(tmp1, p0, tmp2);
    mpc_get_cdpe(ctmp, tmp1);
    cdpe_mod(rtmp, ctmp);
    mpc_mul(tmp1, tmp1, pp0);
    mpc_mul(tmp1, x, tmp1);

    mpc_add(tmp1, tmp1, tmp1);
    mpc_add(tmp1, tmp1, tmp1);

    mpc_add(pp2, tmp1, pp2);
    mpc_mul(tmp1, p1, pp1);
    mpc_add(tmp1, tmp1, tmp1);
 
    mpc_add(pp2, tmp1, pp2);
    mpc_set(p0, p1);
    mpc_set(p1, p2);
    mpc_set(pp0, pp1);
    mpc_set(pp1, pp2);

   /* error bound */
   /* d2=2 d1 |p1| + 4 |x| d0 |p0|^3 + 
         u( |p2| + 2sqrt 2 |p1|^2 + |x| |p0|^4(4sqrt 2+1) */
    mpc_get_cdpe(ctmp, p2);
    cdpe_mod(ap2, ctmp);
    rdpe_mul(rtmp, rtmp, d0);
    rdpe_mul(rtmp, rtmp, ax);

    rdpe_mul_2exp(rtmp, rtmp,2);  /* 4 |x| |p0|^3 d0 */

    rdpe_mul(d2, d1, ap1);

    rdpe_mul_2exp(d2, d2,1);  /*  2 d1 |p1|  */
    rdpe_add(d2, d2, rtmp);

    rdpe_mul(rtmp1, ap0, ap0);
    rdpe_mul(rtmp1, rtmp1, rtmp1); /*  |p0|^4 */
    rdpe_mul(rtmp1, rtmp1, ax);
    rdpe_mul_d(rtmp1, rtmp1, (double) 6.66); /* |x|(4\sqrt 2 +1)|p0|^4 */
    rdpe_add(rtmp1, rtmp1, ap2);
    rdpe_mul(rtmp, ap1, ap1);
    rdpe_mul_d(rtmp, rtmp, (double) 2.83);
    rdpe_mul_2exp(rtmp, rtmp,1);

    rdpe_mul(rtmp, rtmp, eps);
    rdpe_add(d2, d2, rtmp);

    rdpe_set(d0, d1);
    rdpe_set(d1, d2);
    rdpe_set(ap0, ap1);
    rdpe_set(ap1, ap2);

  }

  mpc_div(corr, p1, pp1);
  *again = rdpe_lt(d2, ap2);
  mpc_get_cdpe(ctmp, pp1);
  cdpe_mod(rtmp, ctmp);
  rdpe_add(rad, ap2, d2);
  rdpe_div_eq(rad,  rtmp);
  rdpe_mul_eq_d(rad, (double) n);

  tmpc_clear(tmp2);
  tmpc_clear(tmp1);
  tmpc_clear(pp2);
  tmpc_clear(pp1);
  tmpc_clear(pp0);
  tmpc_clear(p2);
  tmpc_clear(p1);
  tmpc_clear(p0);
}






