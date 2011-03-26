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
fsort(mps_status* s)
{
  int i;

  for (i = 0; i < s->n; i++) {
    cplx_Re(s->fppc[i]) = cplx_Re(s->froot[i]);
    cplx_Im(s->fppc[i]) = i;
  }
  
  qsort(s->fppc, s->n, sizeof(cplx_t), fcmp);
  
  for (i = 0; i < s->n; i++)
    s->order[i] = (int) cplx_Im(s->fppc[i]);
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
dsort(mps_status* s)
{
  int i;

  for (i = 0; i < s->n; i++) {
    rdpe_set(cdpe_Re(s->dpc1[i]), cdpe_Re(s->droot[i]));
    rdpe_set_d(cdpe_Im(s->dpc1[i]), i);
  }
  
  qsort(s->dpc1, s->n, sizeof(cdpe_t), dcmp);
  
  for (i = 0; i < s->n; i++)
    s->order[i] = (int) rdpe_get_d(cdpe_Im(s->dpc1[i]));
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
msort(mps_status* s)
{
  int i;

  for (i = 0; i < s->n; i++) {
    mpf_set(mpc_Re(s->mfpc1[i]), mpc_Re(s->mroot[i]));
    mpf_set_ui(mpc_Im(s->mfpc1[i]), i);
  }
  
  qsort(s->mfpc1, s->n, sizeof(mpc_t), mcmp);
  
  for (i = 0; i < s->n; i++)
    s->order[i] = (int) mpf_get_d(mpc_Im(s->mfpc1[i]));
}
