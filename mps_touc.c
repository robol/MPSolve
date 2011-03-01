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

/************************************************************
*              FUNCTION FTOUCHNWT                           *
*************************************************************
 false if the disks I and J are Newton-isolated            
************************************************************/
boolean
ftouchnwt(int n, int i, int j)
{
  cplx_t ctmp;
  double t;

  t = DBL_MAX/(2*n);  /*#G added 27/4/98 */
  if (frad[i] >= t || frad[j] >= t) return true;

  cplx_sub(ctmp, froot[i], froot[j]);
  return n * (frad[i] + frad[j]) >= cplx_mod(ctmp);
}

/************************************************************
*              FUNCTION DTOUCHNWT                           *
************************************************************/
boolean
dtouchnwt(int n, int i, int j)
{
  cdpe_t ctmp;
  rdpe_t dtmp1, dtmp2;

  rdpe_add(dtmp1, drad[i], drad[j]);
  rdpe_mul_eq_d(dtmp1, (double) n);
  cdpe_sub(ctmp, droot[i], droot[j]);
  cdpe_mod(dtmp2, ctmp);
  return rdpe_ge(dtmp1, dtmp2);
}

/************************************************************
*              FUNCTION MTOUCHNWT                           *
************************************************************/
boolean
mtouchnwt(int n, int i, int j)
{
  tmpc_t mtmp;
  cdpe_t ctmp;
  rdpe_t dtmp1, dtmp2;

  tmpc_init2(mtmp, mpwp);

  rdpe_add(dtmp1, drad[i], drad[j]);
  rdpe_mul_eq_d(dtmp1, (double) n);
  mpc_sub(mtmp, mroot[i], mroot[j]);
  mpc_get_cdpe(ctmp, mtmp);
  cdpe_mod(dtmp2, ctmp);

  tmpc_clear(mtmp);

  return rdpe_ge(dtmp1, dtmp2);
}

/************************************************************
*              FUNCTION FTOUCHREAL                          *
*************************************************************
 true if the disk intersects the real axis, false otherwise
************************************************************/
boolean
ftouchreal(int n, int i)
{
  if (frad[i] >= DBL_MAX/n) return true;

  return n * frad[i] >= fabs(cplx_Im(froot[i]));
}

/************************************************************
*              FUNCTION DTOUCHREAL                          *
************************************************************/
boolean
dtouchreal(int n, int i)
{
  rdpe_t tmp1, tmp2;

  rdpe_mul_d(tmp1, drad[i], (double) n);
  rdpe_abs(tmp2, cdpe_Im(droot[i]));
  return rdpe_ge(tmp1, tmp2);
}

/************************************************************
*              FUNCTION MTOUCHREAL                          *
************************************************************/
boolean
mtouchreal(int n, int i)
{
  rdpe_t tmp1, tmp2;

  rdpe_mul_d(tmp1, drad[i], (double) n);
  mpf_get_rdpe(tmp2, mpc_Im(mroot[i]));
  rdpe_abs_eq(tmp2);

  return rdpe_ge(tmp1, tmp2);
}

/************************************************************
*              FUNCTION  FTOUCHIMAG                         *
*************************************************************
 true iff the disk intersects the imaginary axis 
************************************************************/
boolean
ftouchimag(int n, int i)
{
  if (frad[i] >= DBL_MAX/n) return true;

  return n * frad[i] >= fabs(cplx_Re(froot[i]));
}

/************************************************************
*              FUNCTION  DTOUCHIMAG                         *
************************************************************/
boolean
dtouchimag(int n, int i)
{
  rdpe_t tmp1, tmp2;

  rdpe_mul_d(tmp1, drad[i], (double) n);
  rdpe_abs(tmp2, cdpe_Re(droot[i]));
  return rdpe_ge(tmp1, tmp2);
}

/************************************************************
*              FUNCTION  MTOUCHIMAG                         *
************************************************************/
boolean
mtouchimag(int n, int i)
{
  rdpe_t tmp1, tmp2;

  rdpe_mul_d(tmp1, drad[i], (double) n);
  mpf_get_rdpe(tmp2, mpc_Re(mroot[i]));
  rdpe_abs_eq(tmp2);

  return rdpe_ge(tmp1, tmp2);
}

/************************************************************
*              FUNCTION  FTOUCHUNIT                         *
*************************************************************
 true if the disk intersects the unit circle, false otherwise
  (n*drad[i]+1 >= |froot[i]|) && (n*drad[i]+|froot[i]| >= 1)
*************************************************************/
boolean
ftouchunit(int n, int i)
{
  double ab, rad;

  if (frad[i] >= DBL_MAX/n) return true;

  rad = n * frad[i];
  ab = cplx_mod(froot[i]);
  return (rad + 1 >= ab) && (rad + ab >= 1);
}

/************************************************************
*              FUNCTION  DTOUCHUNIT                         *
************************************************************/
boolean
dtouchunit(int n, int i)
{
  rdpe_t ab, rad, tmp;

  cdpe_mod(ab, droot[i]);
  rdpe_mul_d(rad, drad[i], (double) n);
  rdpe_add_d(tmp, rad, 1.0);
  if (rdpe_lt(tmp, ab))
    return false;
  rdpe_add(tmp, rad, ab);
  return rdpe_ge(tmp, rdpe_one);
}

/************************************************************
*              FUNCTION  MTOUCHUNIT                         *
************************************************************/
boolean
mtouchunit(int n, int i)
{
  tmpf_t mab;
  rdpe_t ab, rad;

  tmpf_init2(mab, mpwp);

  mpc_mod(mab, mroot[i]);
  mpf_sub_eq_ui(mab, 1);
  mpf_get_rdpe(ab, mab);

  tmpf_clear(mab);

  rdpe_mul_d(rad, drad[i], (double) n);

  if (rdpe_lt(rad, ab))
    return false;
  rdpe_neg_eq(ab);
  return rdpe_gt(rad, ab);
}
