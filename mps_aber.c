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

/***********************************************************
*               SUBROUTINE FABERTH                         *
*             No Selective Correction                      *
***********************************************************/
void
faberth(int j, cplx_t abcorr)
{
  int i;
  cplx_t z;

  cplx_set(abcorr, cplx_zero);
  for (i = 0; i < n; i++) {
    if (i == j)
      continue;
    cplx_sub(z, froot[j], froot[i]);
    cplx_inv_eq(z);
    cplx_add_eq(abcorr, z);
  }
}

/***********************************************************
*               SUBROUTINE DABERTH                         *
*             No Selective Correction                      *
***********************************************************/
void
daberth(int j, cdpe_t abcorr)
{
  int i;
  cdpe_t z;

  cdpe_set(abcorr, cdpe_zero);
  for (i = 0; i < n; i++) {
    if (i == j)
      continue;
    cdpe_sub(z, droot[j], droot[i]);
    cdpe_inv_eq(z);
    cdpe_add_eq(abcorr, z);
  }
}

/***********************************************************
*               SUBROUTINE MABERTH                         *
*             No Selective Correction                      *
***********************************************************/
void
maberth(int j, mpc_t abcorr)
{
  int i;
  cdpe_t z, temp;
  tmpc_t diff;

  tmpc_init2(diff, mpwp);

  cdpe_set(temp, cdpe_zero);
  for (i = 0; i < n; i++) {
    if (i == j)
      continue;
    mpc_sub(diff, mroot[j], mroot[i]);
    mpc_get_cdpe(z, diff);
    cdpe_inv_eq(z);
    cdpe_add_eq(temp, z);
  }
  mpc_set_cdpe(abcorr, temp);

  tmpc_clear(diff);
}

/***********************************************************
*               SUBROUTINE FABERTH_S                       *
*               Selective Correction                       *
***********************************************************/
void
faberth_s(int j, int jc, cplx_t abcorr)
{
  int i, k;
  cplx_t z;

  cplx_set(abcorr, cplx_zero);
  for (i = punt[jc]; i < punt[jc + 1]; i++) {
    k = clust[i];
    if (k == j)
      continue;
    cplx_sub(z, froot[j], froot[k]);
    cplx_inv_eq(z);
    cplx_add_eq(abcorr, z);
  }
}

/***********************************************************
*               SUBROUTINE DABERTH_S                       *
*               Selective Correction                       *
***********************************************************/
void
daberth_s(int j, int jc, cdpe_t abcorr)
{
  int i, k;
  cdpe_t z;

  cdpe_set(abcorr, cdpe_zero);
  for (i = punt[jc]; i < punt[jc + 1]; i++) {
    k = clust[i];
    if (k == j)
      continue;
    cdpe_sub(z, droot[j], droot[k]);
    cdpe_inv_eq(z);
    cdpe_add_eq(abcorr, z);
  }
}

/***********************************************************
*               SUBROUTINE MABERTH_S                       *
*               Selective Correction                       *
***********************************************************/
void
maberth_s(int j, int jc, mpc_t abcorr)
{
  int i, k;
  cdpe_t z, temp;
  tmpc_t diff;

  tmpc_init2(diff, mpwp);

  cdpe_set(temp, cdpe_zero);
  for (i = punt[jc]; i < punt[jc + 1]; i++) {
    k = clust[i];
    if (k == j)
      continue;
    mpc_sub(diff, mroot[j], mroot[k]);
    mpc_get_cdpe(z, diff);
    cdpe_inv_eq(z);
    cdpe_add_eq(temp, z);
  }
  mpc_set_cdpe(abcorr, temp);

  tmpc_clear(diff);
}
