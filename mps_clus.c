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
 *               SUBROUTINE FCLUSTER                       *
 ***********************************************************
 This subroutine makes cluster analysis, i.e., detects
 overlapping disks, where two disks overlap if the distances
 of their centers is less than the sum of their radii
 multiplied by 'nf'.
 Observe that nf=1 then this concept corresponds to overlapping,
 if nf =2*n, this concept corresponds to Newton isolation.
 On input: radius, root.
 On output: the vector clust contains the indices of the
 disks in each overlapping group, while  punt(i) points to the
 index of clust where the i-th group starts. Moreover m_clust(i)
 contains the  multiplicity of the i-th cluster. nclust is the
 number of clusters.
 *********************************************************/
void
fcluster(int nf)
{
  int incr, i, j, itemp;

  incr = 0;
  nclust = 0;
  punt[0] = 0;
  for (i = 0; i < n; i++)
    clust[i] = i;

  for (i = 0; i < n - 1; i++) {
    for (j = incr + 1; j < n; j++)
      if (ftouchnwt(nf, clust[i], clust[j])) {
	incr++;
	itemp = clust[j];
	clust[j] = clust[incr];
	clust[incr] = itemp;
      }
    if (i == incr) {
      nclust++;
      punt[nclust] = i + 1;
      incr++;
    }
  }

  nclust++;
  punt[nclust] = n;
}

/**********************************************************
 *                 SUBROUTINE DCLUSTER                    *
 *********************************************************/
void
dcluster(int nf)
{
  int incr, i, j, itemp;

  incr = 0;
  nclust = 0;
  punt[0] = 0;
  for (i = 0; i < n; i++)
    clust[i] = i;

  for (i = 0; i < n - 1; i++) {
    for (j = incr + 1; j < n; j++)
      if (dtouchnwt(nf, clust[i], clust[j])) {
	incr++;
	itemp = clust[j];
	clust[j] = clust[incr];
	clust[incr] = itemp;
      }
    if (i == incr) {
      nclust++;
      punt[nclust] = i + 1;
      incr++;
    }
  }

  nclust++;
  punt[nclust] = n;
}

/**********************************************************
 *                 SUBROUTINE XCLUSTER                    *
 *********************************************************/
void
xcluster(int n, int nf, int *nclust)
{
  int incr, i, j, itemp;

  incr = 0;
  *nclust = 0;
  punt_aux[0] = 0;

  for (i = 0; i < n - 1; i++) {
    for (j = incr + 1; j < n; j++)
      if (mtouchnwt(nf, clust_aux[i], clust_aux[j])) {
	incr++;
	itemp = clust_aux[j];
	clust_aux[j] = clust_aux[incr];
	clust_aux[incr] = itemp;
      }
    if (i == incr) {
      (*nclust)++;
      punt_aux[*nclust] = i + 1;
      incr++;
    }
  }

  (*nclust)++;
  punt_aux[*nclust] = n;
}

/**********************************************************
 *                 SUBROUTINE MCLUSTER                    *
 **********************************************************
 Perform cluster analysis to each existing cluster by
 applying xcluster to each existing cluster, rebuild the
 vectors  clust, punt, and the integer nclust
 *********************************************************/
void
mcluster(int nf)
{
  int i, j, ind, n_aux, nclust_out, nclust_aux;

  punt[nclust] = n;
  nclust_out = 0;
  ind = 0;

  for (i = 0; i < nclust; i++) {
    n_aux = punt[i + 1] - punt[i];
    for (j = 0; j < n_aux; j++) {
      clust_aux[j] = clust[punt[i] + j];
      punt_aux[j] = j;
    }
    if (n_aux > 1) {
      xcluster(n_aux, nf, &nclust_aux);
      for (j = 0; j < n_aux; j++)
	clust_out[punt[i] + j] = clust_aux[j];
      for (j = 0; j < nclust_aux; j++)
	punt_out[ind + j] = punt[i] + punt_aux[j];
      nclust_out = nclust_out + nclust_aux;
    } else {
      nclust_aux = 1;
      clust_out[punt[i]] = clust[punt[i]];
      punt_out[ind] = punt[i];
      nclust_out = nclust_out + nclust_aux;
    }
    ind = ind + nclust_aux;
  }

  nclust = nclust_out;
  for (i = 0; i < n; i++)
    punt[i] = punt_out[i];
  for (i = 0; i < n; i++)
    clust[i] = clust_out[i];
  punt[nclust] = n;
}










