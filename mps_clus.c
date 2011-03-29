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
 * @brief Implementation of the routines to perform cluster analysis
 */

/**
 * This subroutine makes cluster analysis, i.e., detects
 * overlapping disks, where two disks overlap if the distances
 * of their centers is less than the sum of their radii
 * multiplied by <code>nf</code>.
 *
 * Observe that \f$nf=1\f$ then this concept corresponds to overlapping,
 * if \f$nf =2 \cdot n\f$, this concept corresponds to Newton isolation.
 *
 * This routine set the vector <code>clust</code> so that it
 * contains the indices of the
 * disks in each overlapping group, while  <code>punt[i]</code>
 * points to the
 * index of <code>clust</code> where the i-th group starts. Moreover 
 * <code>m_clust[i]</code>
 * contains the  multiplicity of the i-th cluster. 
 * <code>nclust</code> is the
 * number of clusters.
 *
 * @param nf see above for a detailed description.
 */
void
mps_fcluster(mps_status* s, int nf)
{
  int incr, i, j, itemp;

  incr = 0;
  s->nclust = 0;
  s->punt[0] = 0;
  for (i = 0; i < s->n; i++)
    s->clust[i] = i;

  for (i = 0; i < s->n - 1; i++) {
    for (j = incr + 1; j < s->n; j++)
      if (mps_ftouchnwt(s, nf, s->clust[i], s->clust[j])) {
	incr++;
	itemp = s->clust[j];
	s->clust[j] = s->clust[incr];
	s->clust[incr] = itemp;
      }
    if (i == incr) {
      s->nclust++;
      s->punt[s->nclust] = i + 1;
      incr++;
    }
  }

  s->nclust++;
  s->punt[s->nclust] = s->n;
}

/**********************************************************
 *                 SUBROUTINE DCLUSTER                    *
 *********************************************************/
void
mps_dcluster(mps_status* s, int nf)
{
  int incr, i, j, itemp;

  incr = 0;
  s->nclust = 0;
  s->punt[0] = 0;
  for (i = 0; i < s->n; i++)
    s->clust[i] = i;

  for (i = 0; i < s->n - 1; i++) {
    for (j = incr + 1; j < s->n; j++)
      if (mps_dtouchnwt(s, nf, s->clust[i], s->clust[j])) {
	incr++;
	itemp = s->clust[j];
	s->clust[j] = s->clust[incr];
	s->clust[incr] = itemp;
      }
    if (i == incr) {
      s->nclust++;
      s->punt[s->nclust] = i + 1;
      incr++;
    }
  }

  s->nclust++;
  s->punt[s->nclust] = s->n;
}

/**********************************************************
 *                 SUBROUTINE XCLUSTER                    *
 *********************************************************/
void
mps_xcluster(mps_status* s, int n, int nf, int *nclust)
{
  int incr, i, j, itemp;

  incr = 0;
  *(nclust) = 0;
  s->punt_aux[0] = 0;

  for (i = 0; i < s->n - 1; i++) {
    for (j = incr + 1; j < s->n; j++)
      if (mps_mtouchnwt(s, nf, s->clust_aux[i], s->clust_aux[j])) {
	incr++;
	itemp = s->clust_aux[j];
	s->clust_aux[j] = s->clust_aux[incr];
	s->clust_aux[incr] = itemp;
      }
    if (i == incr) {
      (*nclust)++;
      s->punt_aux[*nclust] = i + 1;
      incr++;
    }
  }

  (*nclust)++;
  s->punt_aux[*nclust] = s->n;
}

/**********************************************************
 *                 SUBROUTINE MCLUSTER                    *
 **********************************************************
 Perform cluster analysis to each existing cluster by
 applying xcluster to each existing cluster, rebuild the
 vectors  clust, punt, and the integer nclust
 *********************************************************/
void
mps_mcluster(mps_status* s, int nf)
{
  int i, j, ind, n_aux, nclust_out, nclust_aux;

  s->punt[s->nclust] = s->n;
  nclust_out = 0;
  ind = 0;

  for (i = 0; i < s->nclust; i++) {
    n_aux = s->punt[i + 1] - s->punt[i];
    for (j = 0; j < n_aux; j++) {
      s->clust_aux[j] = s->clust[s->punt[i] + j];
      s->punt_aux[j] = j;
    }
    if (n_aux > 1) {
      mps_xcluster(s, n_aux, nf, &nclust_aux);
      for (j = 0; j < n_aux; j++)
	s->clust_out[s->punt[i] + j] = s->clust_aux[j];
      for (j = 0; j < nclust_aux; j++)
	s->punt_out[ind + j] = s->punt[i] + s->punt_aux[j];
      nclust_out = nclust_out + nclust_aux;
    } else {
      nclust_aux = 1;
      s->clust_out[s->punt[i]] = s->clust[s->punt[i]];
      s->punt_out[ind] = s->punt[i];
      nclust_out = nclust_out + nclust_aux;
    }
    ind = ind + nclust_aux;
  }

  s->nclust = nclust_out;
  for (i = 0; i < s->n; i++)
    s->punt[i] = s->punt_out[i];
  for (i = 0; i < s->n; i++)
    s->clust[i] = s->clust_out[i];
  s->punt[s->nclust] = s->n;
}










