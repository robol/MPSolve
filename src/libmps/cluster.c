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

#include <string.h>
#include <assert.h>
#include <mps/gmptools.h>
#include <mps/core.h>

/**
 * @brief Reset cluster structure information contained in <code>s</code>. After
 * the call to this routine the roots will be considered as a unique big cluster,
 * discarding every information present before.
 *
 * @param s the mps_status pointer.
 */
void
mps_cluster_reset (mps_status * s)
{
  int i;
  for (i = 0; i < s->n; i++)
      s->clust[i] = i;

  for(i = 0; i <= s->n; i++)
    {
      s->punt[i] = 0;
      s->oldpunt[i] = 0;
    }

  /* Reset delimiters */
  s->punt[0] = 0;
  s->punt[1] = s->n;
  s->oldpunt[0] = 0;
  s->oldpunt[1] = s->n;
  s->nclust = 1;

  /* Reset cluster status of the roots */
  for (i = 0; i < s->n; i++)
    {
      s->status[i][0] = 'c';
      s->status[i][1] = 'w';
      s->status[i][2] = 'u';
    }

}

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
 * @param s  The <code>mps_status</code> associated with the current
 *           computaion.
 * @param nf see above for a detailed description.
 */
void
mps_fcluster (mps_status * s, int nf)
{
  int incr, i, j, itemp;

  incr = 0;
  s->nclust = 0;
  s->punt[0] = 0;

  for (i = 0; i < s->n; i++)
    s->clust[i] = i;

  for (i = 0; i < s->n - 1; i++)
    {
      for (j = incr + 1; j < s->n; j++)
        if (mps_ftouchnwt (s, nf, s->clust[i], s->clust[j]))
          {
            incr++;
            itemp = s->clust[j];
            s->clust[j] = s->clust[incr];
            s->clust[incr] = itemp;
          }
      if (i == incr)
        {
          s->nclust++;
          s->punt[s->nclust] = i + 1;
          incr++;
        }


    }

  s->nclust++;
  s->punt[s->nclust] = s->n;

  if (s->DOLOG)
    mps_debug_cluster_structure (s);
}

/**
 * @brief Perform cluster analysis to each existing cluster by
 * applying <code>mps_xcluster</code> to each existing cluster.
 *
 * Rebuild the vectors <code>s->clust</code>,
 * <code>s->punt</code>, and the integer <code>s->nclust</code>.
 *
 * @see mps_xcluster
 */
void
mps_dcluster (mps_status * s, int nf)
{
  int incr, i, j, itemp;

  incr = 0;
  s->nclust = 0;
  s->punt[0] = 0;
  for (i = 0; i < s->n; i++)
    s->clust[i] = i;

  for (i = 0; i < s->n - 1; i++)
    {
      for (j = incr + 1; j < s->n; j++)
        if (mps_dtouchnwt (s, nf, s->clust[i], s->clust[j]))
          {
            incr++;
            itemp = s->clust[j];
            s->clust[j] = s->clust[incr];
            s->clust[incr] = itemp;
          }
      if (i == incr)
        {
          s->nclust++;
          s->punt[s->nclust] = i + 1;
          incr++;
        }
    }

  s->nclust++;
  s->punt[s->nclust] = s->n;

  if (s->DOLOG)
    mps_debug_cluster_structure (s);
}

void
mps_debug_cluster_structure (mps_status * s)
{
  int i, j;

  if (!(s->debug_level & MPS_DEBUG_CLUSTER))
    return;

  /* Debug the clusters found */
  if (s->DOLOG)
    {
      /* Debug isolated roots first */
      __MPS_DEBUG (s, "Isolated roots: ");
      for (i = 0; i < s->nclust; i++)
        {
          if (s->punt[i+1] == s->punt[i] + 1)
            {
              fprintf (s->logstr, "%d ", s->clust[s->punt[i]]);
            }
        }
      fprintf (s->logstr, "\n");

      for (i = 0; i < s->nclust; i++)
        {
          if (s->punt[i + 1] == s->punt[i] + 1)
            continue;
          __MPS_DEBUG (s, "Found cluster of %d roots: ",
                       s->punt[i + 1] - s->punt[i]);
          for (j = s->punt[i]; j < s->punt[i + 1]; j++)
            {
              fprintf (s->logstr, "%d ", s->clust[j]);
            }
          fprintf (s->logstr, "\n");
        }
    }
}

/**
 * @brief Cluster analysis on a cluster whose roots are store in clust_aux.
 *
 * @param s mps_status struct pointer;
 * @param n number of roots in the cluster;
 * @param nf see mps_fcluster();
 * @param nclust pointer to nclust;
 */
void
mps_xcluster (mps_status * s, int n, int nf, int *nclust)
{
  int incr, i, j, itemp;

  incr = 0;
  *(nclust) = 0;
  s->punt_aux[0] = 0;

  for (i = 0; i < n - 1; i++)
    {
      for (j = incr + 1; j < n; j++)
        if (mps_mtouchnwt (s, nf, s->clust_aux[i], s->clust_aux[j]))
          {
            incr++;
            itemp = s->clust_aux[j];
            s->clust_aux[j] = s->clust_aux[incr];
            s->clust_aux[incr] = itemp;
          }
      if (i == incr)
        {
          (*nclust)++;
          s->punt_aux[*nclust] = i + 1;
          incr++;
        }
    }

  (*nclust)++;
  s->punt_aux[*nclust] = n;

}

/**
 * @brief Perform cluster analysis to each existing cluster by
 * applying <code>mps_xcluster</code> to each existing cluster.
 *
 *
 * Rebuild the vectors <code>s->clust</code>,
 * <code>s->punt</code>, and the integer <code>s->nclust</code>.
 *
 * @see mps_xcluster
 */
void
mps_mcluster (mps_status * s, int nf)
{
  int i, j, ind, n_aux, nclust_out, nclust_aux;

  s->punt[s->nclust] = s->n;
  nclust_out = 0;
  ind = 0;

  for (i = 0; i < s->nclust; i++)
    {
      n_aux = s->punt[i + 1] - s->punt[i];
      for (j = 0; j < n_aux; j++)
        {
          s->clust_aux[j] = s->clust[s->punt[i] + j];
          s->punt_aux[j] = j;
        }
      if (n_aux > 1)
        {
          mps_xcluster (s, n_aux, nf, &nclust_aux);
          for (j = 0; j < n_aux; j++)
            s->clust_out[s->punt[i] + j] = s->clust_aux[j];
          for (j = 0; j < nclust_aux; j++)
            s->punt_out[ind + j] = s->punt[i] + s->punt_aux[j];
          nclust_out = nclust_out + nclust_aux;
        }
      else
        {
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

  if (s->DOLOG)
    mps_debug_cluster_structure (s);
}

void
mps_cluster_detachment_reset (mps_status * s)
{
  memset (s->clust_detached, -1, sizeof (int) * s->n);
}


/**
 * @brief Check if in the cluster <code>i_clust</code> there are quasi
 * approximated roots and detach them from the cluster into a new one.
 *
 * @param s the pointer to the mps_status struct that is holding
 * the current status of the computation.
 * @param i_clust The index of the cluster to analyze. The special
 * value <code>MPS_ALL_CLUSTERS</code> can be used to analyze all
 * clusters.
 */
void
mps_cluster_detach (mps_status * s, int i_clust)
{
  MPS_DEBUG_THIS_CALL;

  int i, ind, n_aux, j;
  rdpe_t precision, rtmp;
  mpf_t ftmp;

  if (s->debug_level & MPS_DEBUG_CLUSTER)
    {
      MPS_DEBUG (s, "Debugging cluster structure before root detaching");
      mps_debug_cluster_structure (s);
    }

  mpf_init2 (ftmp, s->mpwp);

  /* Reset the s->clust_detached vector */
  if (i_clust == MPS_ALL_CLUSTERS)
    memset (s->clust_detached, -1, sizeof (int) * s->n);
  else
    memset (s->clust_detached + s->punt[i_clust], -1,
            sizeof (int) * (s->punt[i_clust + 1] - s->punt[i_clust]));

  /* Try to remove approximated roots from the clusters, because they
   * are likely to be "fake" cluster elements. */
  for (i = 0; i < s->nclust; i++)
    {

      /* If this is not the cluster that we have to analyze we should
       * try with next one */
      if (i_clust != i && i_clust != MPS_ALL_CLUSTERS)
        {
          continue;
        }

      if (s->punt[i + 1] - s->punt[i] == 1)
        {
          /* If this is a single root cluster is not a cluster
           * so skip to the next one. */
          continue;
        }

      /* Else keep away approximated roots */
      for (j = s->punt[i]; j < s->punt[i + 1]; j++)
        {
          mpc_mod (ftmp, s->mroot[s->clust[j]]);
          mpf_get_rdpe (rtmp, ftmp);
          rdpe_set_dl (precision, 1, (long int) ((1 - 0.5 * s->mpwp) * LOG10_2
                                                 + rdpe_log10 (rtmp)));
          if (rdpe_lt (s->drad[s->clust[j]], precision))
            {

              MPS_DEBUG (s, "Separating root %d from the "
                         "rest of the cluster n°%d", s->clust[j], i);

              /* Save a log of the detachement in s->clust_detach
               * to make checking if the root is really outside the
               * cluster after the computing of Newton polygonal.
               * We are creating a new cluster in i+1, moving other
               * cluster ahead, so first move ahead the vector after
               * i+1, and then set the i+1 position to i.
               * In theory we should check if s->clust_detached[j] > i
               * and in that case shift it to s->clust_detached[j] + 1,
               * but that's not possible becase cluster after this are not
               * yet analized.
               */
              for (ind = i + 1; ind < s->nclust; ind++)
                {
                  s->clust_detached[ind + 1] = s->clust_detached[ind];
                }
              s->clust_detached[i + 1] = i;


              /* Move other roots back in the cluster */
              n_aux = s->clust[j];
              for (ind = j + 1; ind < s->punt[i + 1]; ind++)
                {
                  s->clust[ind - 1] = s->clust[ind];
                }

              s->clust[s->punt[i + 1] - 1] = n_aux;
              s->punt[i + 1]--;

              /* Move ahead s->punt */
              for (ind = s->nclust; ind > i + 1; ind--)
                {
                  s->punt[ind + 1] = s->punt[ind];
                }

              /* Set s->punt */
              s->punt[i + 2] = s->punt[i + 1] + 1;
	      s->nclust++;

              /* Start from the next root, that is shifted one position back */
              j--;

              /* If the cluster is now a single element cluster, let's
                 skip to the next one */
              if (s->punt[i + 1] - s->punt[i] == 1)
                {
                  break;
                }
            }
        }
    }

  mpf_clear (ftmp);

  if (s->debug_level & MPS_DEBUG_CLUSTER)
    {
      MPS_DEBUG (s, "Debugging cluster structure after root detaching");
      mps_debug_cluster_structure (s);
    }


}


void
mps_cluster_reassemble (mps_status * s, int i_clust)
{
  MPS_DEBUG_THIS_CALL;

  int i, l, j;

  if (s->debug_level & MPS_DEBUG_CLUSTER)
    {
      MPS_DEBUG (s, "Debugging cluster structure before reassembling the original one");
      mps_debug_cluster_structure (s);
    }

  if (i_clust == MPS_ALL_CLUSTERS)
    for (j = 0; j < s->nclust; j++)
      {
        mps_cluster_reassemble (s, j);
        return;
      }

  MPS_DEBUG (s, "Reassembling cluster %d", i_clust);

  for (i = 0; i < s->nclust; i++)
    {

      if (s->clust_detached[i] == i_clust)
        {
          MPS_DEBUG (s, "Recompacting cluster %d and %d", i_clust, i);

          /* We need this to be true to make the reassembling of
           * the cluster work as expected */
          assert (i_clust < i);

          l = s->clust[s->punt[i]];
          for (j = s->punt[i_clust + 1]; j < s->punt[i]; j++)
            {
              s->clust[j + 1] = s->clust[j];
            }
          s->clust[s->punt[i_clust + 1]] = l;
          for (j = i_clust + 1; j <= i; j++)
            {
              s->punt[j]++;
            }

          assert (s->punt[i] == s->punt[i + 1]);

          for (j = i + 1; j < s->nclust; j++)
            {
              s->punt[j - 1] = s->punt[j];
              s->clust_detached[j - 1] = s->clust_detached[j];
            }

          s->punt[s->nclust - 1] = s->punt[s->nclust];
          s->nclust--;
          for (j = 0; j < s->nclust; j++)
            {
              if (s->clust_detached[j] > i)
                {
                  s->clust_detached[j]--;
                }
            }

	  i--;
        }
    }

  if (s->debug_level & MPS_DEBUG_CLUSTER)
    {
      MPS_DEBUG (s, "Debugging cluster structure after restoring the original one");
      mps_debug_cluster_structure (s);
    }
}
