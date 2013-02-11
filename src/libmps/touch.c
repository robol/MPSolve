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
#include <math.h>
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
 * @param s mps_context struct.
 * @param n See above.
 * @param i the first root.
 * @param j the second root.
 * @param frad The inclusion radii precomputed by some other routines. 
 * @return false if the disc <code>i</code> and <code>j</code>
 *   are newton-isolated.
 */
mps_boolean
mps_ftouchnwt (mps_context * s, double * frad, int n, int i, int j)
{
  cplx_t ctmp;
  double t;

  t = DBL_MAX / (2 * n);        /*#G added 27/4/98 */
  if (frad[i] >= t || frad[j] >= t)
    return true;

  cplx_sub (ctmp, s->root[i]->fvalue, s->root[j]->fvalue);

  return n * (frad[i] + frad[j]) >= cplx_mod (ctmp);
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
 * @param s mps_context struct.
 * @param drad The inclusion radii that should be used to perform
 * cluster analysis. 
 * @param n See above.
 * @param i the first root.
 * @param j the second root.
 * @return false if the disc <code>i</code> and <code>j</code>
 *   are newton-isolated.
 */
mps_boolean
mps_dtouchnwt (mps_context * s, rdpe_t * drad, int n, int i, int j)
{
  cdpe_t ctmp;
  rdpe_t dtmp1, dtmp2;

  rdpe_add (dtmp1, drad[i], drad[j]);
      
  rdpe_mul_eq_d (dtmp1, (double) n);
  cdpe_sub (ctmp, s->root[i]->dvalue, s->root[j]->dvalue);
  cdpe_mod (dtmp2, ctmp);
  return rdpe_ge (dtmp1, dtmp2);
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
 * @param s mps_context struct.
 * @param drad The inclusion radii that should be used to perform
 * cluster analysis.
 * @param n See above.
 * @param i the first root.
 * @param j the second root.
 * @return false if the disc <code>i</code> and <code>j</code>
 *   are newton-isolated.
 */
mps_boolean
mps_mtouchnwt (mps_context * s, rdpe_t * drad, int n, int i, int j)
{
  mpc_t mtmp;
  cdpe_t ctmp;
  rdpe_t dtmp1, dtmp2;

  mpc_init2 (mtmp, s->mpwp);

  rdpe_add (dtmp1, drad[i], drad[j]);

  /* if (rdpe_Esp (dtmp1) < rdpe_Esp (s->root[i]->drad) || */
  /*     rdpe_Esp (dtmp1) < rdpe_Esp (s->root[j]->drad)) */
  /*     return true; */

  rdpe_mul_eq_d (dtmp1, (double) n);
  mpc_sub (mtmp, s->root[i]->mvalue, s->root[j]->mvalue);
  mpc_get_cdpe (ctmp, mtmp);
  cdpe_mod (dtmp2, ctmp);

  mpc_clear (mtmp);

  return rdpe_ge (dtmp1, dtmp2);
}

/************************************************************
*              FUNCTION FTOUCHREAL                          *
*************************************************************
 true if the disk intersects the real axis, false otherwise
************************************************************/
mps_boolean
mps_ftouchreal (mps_context * s, int n, int i)
{
  if (s->root[i]->frad >= DBL_MAX / n)
    return true;

  return n * s->root[i]->frad >= fabs (cplx_Im (s->root[i]->fvalue));
}

/************************************************************
*              FUNCTION DTOUCHREAL                          *
************************************************************/
mps_boolean
mps_dtouchreal (mps_context * s, int n, int i)
{
  rdpe_t tmp1, tmp2;

  rdpe_mul_d (tmp1, s->root[i]->drad, (double) n);
  rdpe_abs (tmp2, cdpe_Im (s->root[i]->dvalue));
  return rdpe_ge (tmp1, tmp2);
}

/************************************************************
*              FUNCTION MTOUCHREAL                          *
************************************************************/
mps_boolean
mps_mtouchreal (mps_context * s, int n, int i)
{
  rdpe_t tmp1, tmp2;

  rdpe_mul_d (tmp1, s->root[i]->drad, (double) n);
  mpf_get_rdpe (tmp2, mpc_Im (s->root[i]->mvalue));
  rdpe_abs_eq (tmp2);

  return rdpe_ge (tmp1, tmp2);
}

/************************************************************
*              FUNCTION  FTOUCHIMAG                         *
*************************************************************
 true iff the disk intersects the imaginary axis 
************************************************************/
mps_boolean
mps_ftouchimag (mps_context * s, int n, int i)
{
  if (s->root[i]->frad >= DBL_MAX / n)
    return true;

  return n * s->root[i]->frad >= fabs (cplx_Re (s->root[i]->fvalue));
}

/************************************************************
*              FUNCTION  DTOUCHIMAG                         *
************************************************************/
mps_boolean
mps_dtouchimag (mps_context * s, int n, int i)
{
  rdpe_t tmp1, tmp2;

  rdpe_mul_d (tmp1, s->root[i]->drad, (double) n);
  rdpe_abs (tmp2, cdpe_Re (s->root[i]->dvalue));
  return rdpe_ge (tmp1, tmp2);
}

/************************************************************
*              FUNCTION  MTOUCHIMAG                         *
************************************************************/
mps_boolean
mps_mtouchimag (mps_context * s, int n, int i)
{
  rdpe_t tmp1, tmp2;

  rdpe_mul_d (tmp1, s->root[i]->drad, (double) n);
  mpf_get_rdpe (tmp2, mpc_Re (s->root[i]->mvalue));
  rdpe_abs_eq (tmp2);

  return rdpe_ge (tmp1, tmp2);
}

/************************************************************
*              FUNCTION  FTOUCHUNIT                         *
*************************************************************
 true if the disk intersects the unit circle, false otherwise
  (n*drad[i]+1 >= |froot[i]|) && (n*drad[i]+|froot[i]| >= 1)
*************************************************************/
mps_boolean
mps_ftouchunit (mps_context * s, int n, int i)
{
  double ab, rad;

  if (s->root[i]->frad >= DBL_MAX / n)
    return true;

  rad = n * s->root[i]->frad;
  ab = cplx_mod (s->root[i]->fvalue);
  return (rad + 1 >= ab) && (rad + ab >= 1);
}

/************************************************************
*              FUNCTION  DTOUCHUNIT                         *
************************************************************/
mps_boolean
mps_dtouchunit (mps_context * s, int n, int i)
{
  rdpe_t ab, rad, tmp;

  cdpe_mod (ab, s->root[i]->dvalue);
  rdpe_mul_d (rad, s->root[i]->drad, (double) n);
  rdpe_add_d (tmp, rad, 1.0);
  if (rdpe_lt (tmp, ab))
    return false;
  rdpe_add (tmp, rad, ab);
  return rdpe_ge (tmp, rdpe_one);
}

/************************************************************
*              FUNCTION  MTOUCHUNIT                         *
************************************************************/
mps_boolean
mps_mtouchunit (mps_context * s, int n, int i)
{
  mpf_t mab;
  rdpe_t ab, rad;

  mpf_init2 (mab, s->mpwp);

  mpc_mod (mab, s->root[i]->mvalue);
  mpf_sub_eq_ui (mab, 1);
  mpf_get_rdpe (ab, mab);

  mpf_clear (mab);

  rdpe_mul_d (rad, s->root[i]->drad, (double) n);

  if (rdpe_lt (rad, ab))
    return false;
  rdpe_neg_eq (ab);
  return rdpe_gt (rad, ab);
}
