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

/**
 * @file
 * @brief Implementation of sorting algorithms
 *
 * In this file are implemented the routines to sort roots
 * comparing their real parts.
 */

/********************************************************
*      SUBROUTINE FCMP                                  *
*********************************************************/
int
fcmp(const void *a, const void *b)
{
  double dif;

  dif = cplx_Re(cplx_Addr(a)) - cplx_Re(cplx_Addr(b));
  if (dif > 0)
    return 1;
  else if (dif == 0.0)
    return 0;
  return -1;
}

/*********************************************************
*      SUBROUTINE FSORT                                  *
*********************************************************/
void
fsort(void)
{
  int i;

  for (i = 0; i < n; i++) {
    cplx_Re(fppc[i]) = cplx_Re(froot[i]);
    cplx_Im(fppc[i]) = i;
  }
  
  qsort(fppc, n, sizeof(cplx_t), fcmp);
  
  for (i = 0; i < n; i++)
    order[i] = (int) cplx_Im(fppc[i]);
}

/*********************************************************
*      SUBROUTINE DCMP                                  *
*********************************************************/
int
dcmp(const void *a, const void *b)
{
  return rdpe_cmp(cdpe_Re(cdpe_Addr(a)), cdpe_Re(cdpe_Addr(b)));
}

/*********************************************************
*      SUBROUTINE DSORT                                  *
*********************************************************/
void
dsort(void)
{
  int i;

  for (i = 0; i < n; i++) {
    rdpe_set(cdpe_Re(dpc1[i]), cdpe_Re(droot[i]));
    rdpe_set_d(cdpe_Im(dpc1[i]), i);
  }
  
  qsort(dpc1, n, sizeof(cdpe_t), dcmp);
  
  for (i = 0; i < n; i++)
    order[i] = (int) rdpe_get_d(cdpe_Im(dpc1[i]));
}

/*********************************************************
*      SUBROUTINE MCMP                                  *
*********************************************************/
int
mcmp(const void *a, const void *b)
{
  return mpf_cmp(mpc_Re(mpc_Addr(a)), mpc_Re(mpc_Addr(b)));
}

/*********************************************************
*      SUBROUTINE MSORT                                  *
*********************************************************/
void
msort(void)
{
  int i;

  for (i = 0; i < n; i++) {
    mpf_set(mpc_Re(mfpc1[i]), mpc_Re(mroot[i]));
    mpf_set_ui(mpc_Im(mfpc1[i]), i);
  }
  
  qsort(mfpc1, n, sizeof(mpc_t), mcmp);
  
  for (i = 0; i < n; i++)
    order[i] = (int) mpf_get_d(mpc_Im(mfpc1[i]));
}
