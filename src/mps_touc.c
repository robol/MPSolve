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

#include <mps/mps.h>

/**
 * @brief Check if the i-th and the j-th discs are newton-isolated.
 *
 * More precisely, given a parameter <code>n</code>, check if
 * the roots <code>i</code> and <code>j</code> are separated
 * with circles whose radius is less than their distance divided
 * for <code>n</code>.
 *
 * If \f$n = 1\f$ this condition correspond to isolation,
 * if \f$n = 2*m\f$ where \f$m\f$ is the degree of the polynomial
 * then it correspond to newton isolation.
 *
 * @param s mps_status struct.
 * @param n See above.
 * @param i the first root.
 * @param j the second root.
 * @return false if the disc <code>i</code> and <code>j</code>
 *   are newton-isolated.
 */
boolean
mps_ftouchnwt(mps_status* s, int n, int i, int j)
{
  cplx_t ctmp;
  double t;

  t = DBL_MAX/(2*n);  /*#G added 27/4/98 */
  if (s->frad[i] >= t || s->frad[j] >= t) return true;

  cplx_sub(ctmp, s->froot[i], s->froot[j]);
  return n * (s->frad[i] + s->frad[j]) >= cplx_mod(ctmp);
}

/**
 * @brief Check if the i-th and the j-th discs are newton-isolated.
 *
 * More precisely, given a parameter <code>n</code>, check if
 * the roots <code>i</code> and <code>j</code> are separated
 * with circles whose radius is less than their distance divided
 * for <code>n</code>.
 *
 * If \f$n = 1\f$ this condition correspond to isolation,
 * if \f$n = 2*m\f$ where \f$m\f$ is the degree of the polynomial
 * then it correspond to newton isolation.
 *
 * @param s mps_status struct.
 * @param n See above.
 * @param i the first root.
 * @param j the second root.
 * @return false if the disc <code>i</code> and <code>j</code>
 *   are newton-isolated.
 */
boolean
mps_dtouchnwt(mps_status* s, int n, int i, int j)
{
  cdpe_t ctmp;
  rdpe_t dtmp1, dtmp2;

  rdpe_add(dtmp1, s->drad[i], s->drad[j]);
  rdpe_mul_eq_d(dtmp1, (double) n);
  cdpe_sub(ctmp, s->droot[i], s->droot[j]);
  cdpe_mod(dtmp2, ctmp);
  return rdpe_ge(dtmp1, dtmp2);
}

/**
 * @brief Check if the i-th and the j-th discs are newton-isolated.
 *
 * More precisely, given a parameter <code>n</code>, check if
 * the roots <code>i</code> and <code>j</code> are separated
 * with circles whose radius is less than their distance divided
 * for <code>n</code>.
 *
 * If \f$n = 1\f$ this condition correspond to isolation,
 * if \f$n = 2*m\f$ where \f$m\f$ is the degree of the polynomial
 * then it correspond to newton isolation.
 *
 * @param s mps_status struct.
 * @param n See above.
 * @param i the first root.
 * @param j the second root.
 * @return false if the disc <code>i</code> and <code>j</code>
 *   are newton-isolated.
 */
boolean
mps_mtouchnwt(mps_status* s, int n, int i, int j)
{
  tmpc_t mtmp;
  cdpe_t ctmp;
  rdpe_t dtmp1, dtmp2;

  tmpc_init2(mtmp, s->mpwp);

  rdpe_add(dtmp1, s->drad[i], s->drad[j]);
  rdpe_mul_eq_d(dtmp1, (double) n);
  mpc_sub(mtmp, s->mroot[i], s->mroot[j]);
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
mps_ftouchreal(mps_status* s, int n, int i)
{
  if (s->frad[i] >= DBL_MAX/n) return true;

  return n * s->frad[i] >= fabs(cplx_Im(s->froot[i]));
}

/************************************************************
*              FUNCTION DTOUCHREAL                          *
************************************************************/
boolean
mps_dtouchreal(mps_status* s, int n, int i)
{
  rdpe_t tmp1, tmp2;

  rdpe_mul_d(tmp1, s->drad[i], (double) n);
  rdpe_abs(tmp2, cdpe_Im(s->droot[i]));
  return rdpe_ge(tmp1, tmp2);
}

/************************************************************
*              FUNCTION MTOUCHREAL                          *
************************************************************/
boolean
mps_mtouchreal(mps_status* s, int n, int i)
{
  rdpe_t tmp1, tmp2;

  rdpe_mul_d(tmp1, s->drad[i], (double) n);
  mpf_get_rdpe(tmp2, mpc_Im(s->mroot[i]));
  rdpe_abs_eq(tmp2);

  return rdpe_ge(tmp1, tmp2);
}

/************************************************************
*              FUNCTION  FTOUCHIMAG                         *
*************************************************************
 true iff the disk intersects the imaginary axis 
************************************************************/
boolean
mps_ftouchimag(mps_status* s, int n, int i)
{
  if (s->frad[i] >= DBL_MAX/n) return true;

  return n * s->frad[i] >= fabs(cplx_Re(s->froot[i]));
}

/************************************************************
*              FUNCTION  DTOUCHIMAG                         *
************************************************************/
boolean
mps_dtouchimag(mps_status* s, int n, int i)
{
  rdpe_t tmp1, tmp2;

  rdpe_mul_d(tmp1, s->drad[i], (double) n);
  rdpe_abs(tmp2, cdpe_Re(s->droot[i]));
  return rdpe_ge(tmp1, tmp2);
}

/************************************************************
*              FUNCTION  MTOUCHIMAG                         *
************************************************************/
boolean
mps_mtouchimag(mps_status* s, int n, int i)
{
  rdpe_t tmp1, tmp2;

  rdpe_mul_d(tmp1, s->drad[i], (double) n);
  mpf_get_rdpe(tmp2, mpc_Re(s->mroot[i]));
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
mps_ftouchunit(mps_status* s, int n, int i)
{
  double ab, rad;

  if (s->frad[i] >= DBL_MAX/n) return true;

  rad = n * s->frad[i];
  ab = cplx_mod(s->froot[i]);
  return (rad + 1 >= ab) && (rad + ab >= 1);
}

/************************************************************
*              FUNCTION  DTOUCHUNIT                         *
************************************************************/
boolean
mps_dtouchunit(mps_status* s, int n, int i)
{
  rdpe_t ab, rad, tmp;

  cdpe_mod(ab, s->droot[i]);
  rdpe_mul_d(rad, s->drad[i], (double) n);
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
mps_mtouchunit(mps_status* s, int n, int i)
{
  tmpf_t mab;
  rdpe_t ab, rad;

  tmpf_init2(mab, s->mpwp);

  mpc_mod(mab, s->mroot[i]);
  mpf_sub_eq_ui(mab, 1);
  mpf_get_rdpe(ab, mab);

  tmpf_clear(mab);

  rdpe_mul_d(rad, s->drad[i], (double) n);

  if (rdpe_lt(rad, ab))
    return false;
  rdpe_neg_eq(ab);
  return rdpe_gt(rad, ab);
}
