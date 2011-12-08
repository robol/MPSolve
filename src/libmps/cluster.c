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
#include <mps/cluster.h>

/**
 * @brief Get an empty mps_cluster, with no roots.
 */
mps_cluster *
mps_cluster_empty (mps_status * s)
{
  mps_cluster * cluster = mps_new (mps_cluster);
  cluster->first = NULL;
  cluster->n = 0;
  return cluster;
}

/**
 * @brief Create a cluster containing only the selected root.
 *
 * @param root_index The root that must be in the cluster.
 */
mps_cluster * 
mps_cluster_with_root (mps_status * s, long int root_index)
{
  mps_cluster * cluster = mps_new (mps_cluster);
  cluster->first = mps_new (mps_root);
  cluster->n = 1;

  cluster->first->k = root_index;
  cluster->first->next = NULL;
  cluster->first->prev = NULL;

  return cluster;
}

/**
 * @brief Free a previously allocated cluster with all the roots in
 * it. 
 *
 * @param cluster The cluster to free.
 */
void 
mps_cluster_free (mps_status * s, mps_cluster * cluster)
{
  mps_root * root = cluster->first;
  mps_root * old_root;

  /* Free all the roots in the cluster */
  while (root)
    {
      old_root = root;
      root = root->next;
      free (old_root);
    }

  free (cluster);
}

/**
 * @brief Insert a root in a cluster.
 *
 * @param cluster The cluster in which the root must be inserted.
 * @param root_index The index of the root to insert.
 */
mps_root *
mps_cluster_insert_root (mps_status * s, 
			 mps_cluster  * cluster, 
			 long int root_index)
{
  mps_root * root = mps_new (mps_root);

  /* Inserting root in the starting of the cluster */
  root->k = root_index;
  root->next = cluster->first;
  root->prev = NULL;

  cluster->n++;

  if (cluster->first)
    cluster->first->prev = root;

  cluster->first = root;

  return root;
}

/**
 * @brief Remove a root from a cluster.
 *
 * @param cluster The cluster from which the root must be removed.
 * @param root The root to remove.
 *
 * Please note the the root specified must be in the cluster, otherwise
 * an assertion error or segmentation fault will be triggered.
 */
void 
mps_cluster_remove_root (mps_status * s, mps_cluster * cluster, mps_root * root)
{
  /* Iterate over the cluster to find the root searched */
  mps_root * prev_root = root->prev;
  mps_root * next_root = root->next;

  if (prev_root)
    prev_root->next = next_root;
  if (next_root)
    next_root->prev = prev_root;

  if (cluster->first == root)
    cluster->first = root->next;

  /* Decrease root count */
  cluster->n--;

  /* Free the root */
  free (root);
}

/**
 * @brief Join two cluster in one big cluster containing the roots of
 * both. Please note that the cluster must not overlap. 
 * 
 * @param cluster_a The first cluster
 * @param cluster_b The second cluster
 * @return A new cluster containing the roots of both.
 */
mps_cluster *
mps_cluster_join (mps_status * s, mps_cluster * cluster_a, mps_cluster * cluster_b)
{
  mps_root * root;
  mps_cluster * small_cluster;
  mps_cluster * big_cluster;

  mps_cluster * new_cluster = mps_cluster_empty (s);

  if (cluster_a->n < cluster_b->n)
    {
      small_cluster = cluster_a;
      big_cluster = cluster_b;
    }
  else
    {
      small_cluster = cluster_b;
      big_cluster = cluster_a;
    }

  root = small_cluster->first;
  while (root->next != NULL)
    root = root->next;
  root->next = big_cluster->first;
  
  new_cluster->first = small_cluster->first;
  new_cluster->n = small_cluster->n + big_cluster->n;

  return new_cluster;
}

/**
 * @brief Create a new empty clusterization.
 */
mps_clusterization *
mps_clusterization_empty (mps_status * s)
{
  mps_clusterization * c = mps_new (mps_clusterization);
  c->n = 0;
  c->first = NULL;
  return c;
}

/**
 * @brief Insert a new cluster into a root clusterization.
 */
mps_cluster_item *
mps_clusterization_insert_cluster (mps_status * s, mps_clusterization * c, mps_cluster * cluster)
{
  mps_cluster_item * item = mps_new (mps_cluster_item);

  /* Set previous item as NULL and next item as the first now */
  item->prev = NULL;
  item->next = c->first;
  item->detached = NULL;

  /* Store the pointer to the cluster */
  item->cluster = cluster;

  /* The first cluster is the one we're inserting */
  c->first = item;
  
  /* The previous item of before was NULL, now should be the cluster that
   * we are inserting. */
  if (item->next)
    item->next->prev = item;

  c->n++;

  return item;
}

/**
 * @brief Pop out a cluster from a clusterization.
 */
void
mps_clusterization_pop_cluster (mps_status * s, mps_clusterization * c, mps_cluster_item * cluster_item)
{
  mps_cluster_item * next = cluster_item->next;
  mps_cluster_item * prev = cluster_item->prev;

  if (prev)
    prev->next = next;
  if (next)
    next->prev = prev;

  if (c->first == cluster_item)
    c->first = next;

  c->n--;
}

/**
 * @brief Remove a cluster item from a clusterization, freeing it.
 */
void 
mps_clusterization_remove_cluster (mps_status * s, mps_clusterization * c, mps_cluster_item * cluster_item)
{
  mps_clusterization_pop_cluster (s, c, cluster_item);
  mps_cluster_free (s, cluster_item->cluster);
  free (cluster_item);
}

/**
 * @brief Free a clusterization and all the cluster in it.
 */
void
mps_clusterization_free (mps_status * s, mps_clusterization * c)
{
  mps_cluster_item * item = c->first;
  mps_cluster_item * next_item;
  
  while (item != NULL)
    {
      mps_cluster_free (s, item->cluster);
      next_item = item->next;
      free (item);
      item = next_item;
    }

  free (c);
}

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
  /* Reset cluster status of the roots */ 
  int i;
  mps_cluster * cluster;

  for (i = 0; i < s->n; i++) 
    { 
      s->status[i][0] = 'c'; 
      s->status[i][1] = 'w'; 
      s->status[i][2] = 'u'; 
    }

  if (s->clusterization != NULL)
    mps_clusterization_free (s, s->clusterization);
  s->clusterization = mps_clusterization_empty (s);

  /* Fill in the roots */
  cluster = mps_cluster_empty (s);
  for (i = 0; i < s->n; i++)
    mps_cluster_insert_root (s, cluster, i);

  mps_clusterization_insert_cluster (s, s->clusterization, cluster);
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
 * @param frad The vector of radii to use for cluster analysis.
 * @param nf see above for a detailed description.
 */
void
mps_fcluster (mps_status * s, double * frad, int nf)
{
  /* We need to scan every cluster and make it in pieces, if possible */
  mps_clusterization * new_clusterization = mps_clusterization_empty (s);
  mps_cluster_item * item;
  int analyzed_roots = 0;

  /* Debug clusterization status if debugging was required */
  if (s->debug_level & MPS_DEBUG_CLUSTER)
    {
      int i;
      MPS_DEBUG (s, "Debugging the radius and approximations obtained for the roots before cluster analysis");
      for (i = 0; i < s->n; i++)
	{
	  MPS_DEBUG_CPLX (s, s->froot[i], "Root %d", i);
	  MPS_DEBUG (s, "radius for root %4d: %e", i, frad[i]);
	}

      MPS_DEBUG (s, "Debugging cluster structure before cluster analysis");
      mps_debug_cluster_structure (s);
    }

  item = s->clusterization->first;
  while (item)
    {
      mps_cluster * cluster = item->cluster;
      mps_cluster_item * next_item = item->next;

      /* Keep isolated cluster isolated, moving them in the new clusterization. */
      if (cluster->n == 1)
	{
	  mps_clusterization_insert_cluster (s, new_clusterization,
					     mps_cluster_with_root (s, cluster->first->k));
	  mps_clusterization_remove_cluster (s, s->clusterization, item);
	  analyzed_roots++;
	}

      item = next_item;
    }

  /* Now do cluster analysis with the rest of the clusters. */
  item = s->clusterization->first;
  while (analyzed_roots < s->n)
    {
      /* Create a new cluster to be inserted in the cluster analysis */
      mps_root * base_root;
      mps_cluster * cluster = item->cluster;
      mps_cluster * new_cluster = mps_cluster_empty (s);

      while (cluster->n == 0)
	{
	  item = item->next;
	  cluster = item->cluster;
	}
	  
      
      base_root = mps_cluster_insert_root (s, new_cluster, cluster->first->k);
      analyzed_roots++;
      mps_cluster_remove_root (s, cluster, cluster->first);

      /* Check if this root touches others root, and if new roots were added to the 
       * cluster repeat the checks. */
      while (base_root)
	{
	  mps_cluster_item * c_item;
	  mps_root * iter_root;

	  for (c_item = s->clusterization->first; c_item != NULL; c_item = c_item->next)
	    {
	      mps_cluster * iter_cluster = c_item->cluster;

	      iter_root = iter_cluster->first;
	      while (iter_root)
		{
		  if (mps_ftouchnwt (s, frad, nf, base_root->k, iter_root->k))
		    {
		      mps_root * next_root = iter_root->next;
		      mps_cluster_insert_root (s, new_cluster, iter_root->k);
		      mps_cluster_remove_root (s, iter_cluster, iter_root);
		      analyzed_roots++;
		      iter_root = next_root;
		    }
		  else
		    iter_root = iter_root->next;
		}
	    }

	  base_root = base_root->prev;
	}
      
      /* Now insert the cluster in the new clusterization */
      mps_clusterization_insert_cluster (s, new_clusterization, new_cluster);

      if (new_cluster->n == 1)
	{
	  int k = new_cluster->first->k;
	  double new_rad;

	  new_rad = cplx_mod (s->froot[k]) * 4.0f * DBL_EPSILON + frad[k]; 

	  /* Check if the computed radius is more convenient than the old one. 
	     If that's the case, apply it as inclusion radius */ 
	  if (new_rad < s->frad[k]) 
	    s->frad[k] = new_rad; 
	}
    }

  /* Set the new clusterizaition in the mps_status */
  mps_clusterization_free (s, s->clusterization);
  s->clusterization = new_clusterization;

  if (s->debug_level & MPS_DEBUG_CLUSTER)
    {
      MPS_DEBUG (s, "Debugging cluster structure after cluster analysis");
      mps_debug_cluster_structure (s);
    }
}

/**
 * @brief Perform cluster analysis to each existing cluster by
 * applying <code>mps_xcluster</code> to each existing cluster.
 *
 * Rebuild the vectors <code>s->clust</code>,
 * <code>s->punt</code>, and the integer <code>s->nclust</code>.
 *
 */
void
mps_dcluster (mps_status * s, rdpe_t * drad, int nf)
{
  /* We need to scan every cluster and make it in pieces, if possible */
  mps_clusterization * new_clusterization = mps_clusterization_empty (s);
  mps_cluster_item * item;

  /* Debug clusterization status if debugging was required */
  if (s->debug_level & MPS_DEBUG_CLUSTER)
    {
      int i;
      MPS_DEBUG (s, "Debugging the radius and approximations obtained for the roots before cluster analysis");
      for (i = 0; i < s->n; i++)
	{
	  MPS_DEBUG_CDPE (s, s->droot[i], "Root %d", i);
	  MPS_DEBUG_RDPE (s, drad[i], "radius for root %4d", i);
	}

      MPS_DEBUG (s, "Debugging cluster structure before cluster analysis");
      mps_debug_cluster_structure (s);
    }

  int analyzed_roots = 0;

  item = s->clusterization->first;
  while (item)
    {
      mps_cluster * cluster = item->cluster;
      mps_cluster_item * next_item = item->next;

      /* Keep isolated cluster isolated, moving them in the new clusterization. */
      if (cluster->n == 1)
	{
	  mps_clusterization_insert_cluster (s, new_clusterization,
					     mps_cluster_with_root (s, cluster->first->k));
	  mps_clusterization_remove_cluster (s, s->clusterization, item);
	  analyzed_roots++;
	}

      item = next_item;
    }

  /* Now do cluster analysis with the rest of the clusters. */
  item = s->clusterization->first;
  while (analyzed_roots < s->n)
    {
      /* Create a new cluster to be inserted in the cluster analysis */
      mps_root * base_root;
      mps_cluster * cluster = item->cluster;
      mps_cluster * new_cluster = mps_cluster_empty (s);

      while (cluster->n == 0)
	{
	  item = item->next;
	  cluster = item->cluster;
	}
	  
      base_root = mps_cluster_insert_root (s, new_cluster, cluster->first->k);
      analyzed_roots++;
      mps_cluster_remove_root (s, cluster, cluster->first);

      /* Check if this root touches others root, and if new roots were added to the 
       * cluster repeat the checks. */
      while (base_root)
	{
	  mps_cluster_item * c_item;
	  mps_root * iter_root;

	  for (c_item = s->clusterization->first; c_item != NULL; c_item = c_item->next)
	    {
	      mps_cluster * iter_cluster = c_item->cluster;

	      iter_root = iter_cluster->first;
	      while (iter_root)
		{
		  if (mps_dtouchnwt (s, drad, nf, base_root->k, iter_root->k))
		    {
		      mps_root * next_root = iter_root->next;
		      mps_cluster_insert_root (s, new_cluster, iter_root->k);
		      mps_cluster_remove_root (s, iter_cluster, iter_root);
		      analyzed_roots++;
		      iter_root = next_root;
		    }
		  else
		    iter_root = iter_root->next;
		}
	    }

	  base_root = base_root->prev;
	}
      
      /* Now insert the cluster in the new clusterization */
      mps_clusterization_insert_cluster (s, new_clusterization, new_cluster);

      if (new_cluster->n == 1)
	{
	  int k = new_cluster->first->k;
	  rdpe_t new_rad;
	  
	  cdpe_mod (new_rad, s->droot[k]);
	  rdpe_mul_eq_d (new_rad, 4 * DBL_EPSILON);
	  rdpe_add_eq (new_rad, drad[k]);

	  /* Check if the computed radius is more convenient than the old one.
	     If that's the case, apply it as inclusion radius */
	  if (rdpe_lt (new_rad, s->drad[k]))
	    rdpe_set (s->drad[k], new_rad);

	}
    }

  /* Set the new clusterizaition in the mps_status */
  mps_clusterization_free (s, s->clusterization);
  s->clusterization = new_clusterization;

  if (s->debug_level & MPS_DEBUG_CLUSTER)
    {
      MPS_DEBUG (s, "Debugging cluster structure after cluster analysis");
      mps_debug_cluster_structure (s);
    }
}


void
mps_debug_cluster_structure (mps_status * s)
{
  mps_cluster_item * cluster_item;
  mps_cluster * cluster;
  mps_root * root;
  mps_boolean isolated_roots = false;

  if (!(s->debug_level & MPS_DEBUG_CLUSTER)) 
    return; 

  for (cluster_item = s->clusterization->first; cluster_item != NULL; cluster_item = cluster_item->next)
    {
      cluster = cluster_item->cluster;

      /* Detect that there are isolated roots, so they will be listed
       * after. */
      if (cluster->n == 1)
	{
	  isolated_roots = true;
	  continue;
	}

      __MPS_DEBUG (s, "Found cluster of %ld roots: ", 
		   cluster->n); 

      for (root = cluster->first; root != NULL; root = root->next)
	{
	  fprintf (s->logstr, "%ld ", root->k);
	}
      fprintf (s->logstr, "\n");
    }

  if (isolated_roots)
    {
      __MPS_DEBUG (s, "Isolated roots: ");
      for (cluster_item = s->clusterization->first; cluster_item != NULL; cluster_item = cluster_item->next)
	{
	  cluster = cluster_item->cluster;
	  if (cluster->n == 1)
	    fprintf (s->logstr, "%ld ", cluster->first->k);
	}
      fprintf (s->logstr, "\n");
    }
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
mps_mcluster (mps_status * s, rdpe_t * drad, int nf)
{
  /* We need to scan every cluster and make it in pieces, if possible */
  mps_clusterization * new_clusterization = mps_clusterization_empty (s);
  mps_cluster_item * item;

  /* Debug clusterization status if debugging was required */
  if (s->debug_level & MPS_DEBUG_CLUSTER)
    {
      int i;
      MPS_DEBUG (s, "Debugging the radius and approximations obtained for the roots before cluster analysis");
      for (i = 0; i < s->n; i++)
	{
	  MPS_DEBUG_MPC (s, 15, s->mroot[i], "Root %d", i);
	  MPS_DEBUG_RDPE (s, drad[i], "radius for root %4d", i);
	}

      MPS_DEBUG (s, "Debugging cluster structure before cluster analysis");
      mps_debug_cluster_structure (s);
    }

  int analyzed_roots = 0;

  item = s->clusterization->first;
  while (item != NULL)
    {
      mps_cluster * cluster = item->cluster;
      mps_cluster_item * next_item = item->next;

      /* Keep isolated cluster isolated, moving them in the new clusterization. */
      if (cluster->n == 1)
	{
	  mps_clusterization_insert_cluster (s, new_clusterization,
					     mps_cluster_with_root (s, cluster->first->k));
	  mps_clusterization_remove_cluster (s, s->clusterization, item);
	  analyzed_roots++;
	}

      item = next_item;
    }

  /* Now do cluster analysis with the rest of the clusters. */
  item = s->clusterization->first;
  while (analyzed_roots < s->n)
    {
      /* Create a new cluster to be inserted in the cluster analysis */
      mps_root * base_root;
      mps_cluster * cluster = item->cluster;
      mps_cluster * new_cluster = mps_cluster_empty (s);

      while (cluster->n == 0)
	{
	  item = item->next;
	  cluster = item->cluster;
	}
	  
      base_root = mps_cluster_insert_root (s, new_cluster, cluster->first->k);
      analyzed_roots++;
      mps_cluster_remove_root (s, cluster, cluster->first);

      /* Check if this root touches others root, and if new roots were added to the 
       * cluster repeat the checks. */
      while (base_root)
	{
	  mps_cluster_item * c_item;
	  mps_root * iter_root;

	  for (c_item = s->clusterization->first; c_item != NULL; c_item = c_item->next)
	    {
	      mps_cluster * iter_cluster = c_item->cluster;

	      iter_root = iter_cluster->first;
	      while (iter_root)
		{
		  if (mps_mtouchnwt (s, drad, nf, base_root->k, iter_root->k))
		    {
		      mps_root * next_root = iter_root->next;
		      mps_cluster_insert_root (s, new_cluster, iter_root->k);
		      mps_cluster_remove_root (s, iter_cluster, iter_root);
		      analyzed_roots++;
		      iter_root = next_root;
		    }
		  else
		    iter_root = iter_root->next;
		}
	    }

	  base_root = base_root->prev;
	}
      
      /* Now insert the cluster in the new clusterization */
      mps_clusterization_insert_cluster (s, new_clusterization, new_cluster);

      if (new_cluster->n == 1)
	{
	  int k = new_cluster->first->k;
	  cdpe_t c;
	  rdpe_t new_rad;

	  /* Check if the computed radius is more convenient than the old one.
	     If that's the case, apply it as inclusion radius */
	  mpc_get_cdpe (c, s->mroot[k]);
	  cdpe_mod (new_rad, c);
	  rdpe_mul_eq (new_rad, s->mp_epsilon);
	  rdpe_mul_eq_d (new_rad, 4.0f);
	  rdpe_add_eq (new_rad, drad[k]);

	  if (rdpe_lt (new_rad, s->drad[k]))
	    rdpe_set (s->drad[k], new_rad);
	}
    }

  /* Set the new clusterizaition in the mps_status */
  mps_clusterization_free (s, s->clusterization);
  s->clusterization = new_clusterization;

  if (s->debug_level & MPS_DEBUG_CLUSTER)
    {
      MPS_DEBUG (s, "Debugging cluster structure after cluster analysis");
      mps_debug_cluster_structure (s);
    }
}


void 
mps_clusterization_detach_clusters (mps_status * s, mps_clusterization * c)
{
  mps_cluster_item * item;
  cdpe_t droot;
  rdpe_t rtmp, precision;
  int k;
  
  for (item = c->first; item != NULL; item = item->next)
    {
      mps_root * root;
      
      /* Skip isolated clusters */
      if (item->cluster->n == 1)
	continue;

      /* Scan the cluster for quasi approximated roots */
      root = item->cluster->first;
      while (root != NULL)
	{
	  k = root->k;
	  mpc_get_cdpe (droot, s->mroot[k]);
	  cdpe_mod (rtmp, droot);

          rdpe_set_dl (precision, 1, (long int) ((1 - 0.5 * s->mpwp) * LOG10_2
                                                 + rdpe_log10 (rtmp)));
          if (rdpe_lt (s->drad[k], precision))
            {
	      mps_cluster * detached_cluster = mps_cluster_with_root (s, k);
	      mps_root * next_root = root->next;
              MPS_DEBUG (s, "Separating root %d from the "
                         "rest of the cluster", k);
	      mps_cluster_remove_root (s, item->cluster, root);
	      root = next_root;

	      /* Insert the cluster in the clusterization */
	      mps_cluster_item * new_item = mps_clusterization_insert_cluster (s, 
									       s->clusterization, 
									       detached_cluster);
	      /* Set the new item as detached from the old one */
	      new_item->detached = item;
	    }
	  else
	    root = root->next;

	  /* If we have left only an isolated roots stop checking this cluster */
	  if (item->cluster->n == 1)
	    break;
	}
    }
}


void 
mps_clusterization_reassemble_clusters (mps_status * s, mps_clusterization * c)
{
  MPS_DEBUG_THIS_CALL;

  mps_cluster_item * cluster;
  
  cluster = s->clusterization->first;
  while (cluster != NULL)
    {
      if (cluster->detached)
	{
	  mps_clusterization_remove_cluster (s, s->clusterization, cluster);
	  mps_cluster_insert_root (s, cluster->detached->cluster, cluster->cluster->first->k);
	}

      cluster = cluster->next;
    }
}

void
mps_cluster_detachment_reset (mps_status * s)
{
  mps_clusterization_reassemble_clusters (s, s->clusterization);
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
/* void */
/* mps_cluster_detach (mps_status * s, int i_clust) */
/* { */
/*   MPS_DEBUG_THIS_CALL; */

/*   int i, ind, n_aux, j; */
/*   rdpe_t precision, rtmp; */
/*   mpf_t ftmp; */

/*   if (s->debug_level & MPS_DEBUG_CLUSTER) */
/*     { */
/*       MPS_DEBUG (s, "Debugging cluster structure before root detaching"); */
/*       mps_debug_cluster_structure (s); */
/*     } */

/*   mpf_init2 (ftmp, s->mpwp); */

/*   /\* Reset the s->clust_detached vector *\/ */
/*   if (i_clust == MPS_ALL_CLUSTERS) */
/*     memset (s->clust_detached, -1, sizeof (int) * s->n); */
/*   else */
/*     memset (s->clust_detached + s->punt[i_clust], -1, */
/*             sizeof (int) * (s->punt[i_clust + 1] - s->punt[i_clust])); */

/*   /\* Try to remove approximated roots from the clusters, because they */
/*    * are likely to be "fake" cluster elements. *\/ */
/*   for (i = 0; i < s->nclust; i++) */
/*     { */

/*       /\* If this is not the cluster that we have to analyze we should */
/*        * try with next one *\/ */
/*       if (i_clust != i && i_clust != MPS_ALL_CLUSTERS) */
/*         { */
/*           continue; */
/*         } */

/*       if (s->punt[i + 1] - s->punt[i] == 1) */
/*         { */
/*           /\* If this is a single root cluster is not a cluster */
/*            * so skip to the next one. *\/ */
/*           continue; */
/*         } */

/*       /\* Else keep away approximated roots *\/ */
/*       for (j = s->punt[i]; j < s->punt[i + 1]; j++) */
/*         { */
/*           mpc_mod (ftmp, s->mroot[s->clust[j]]); */
/*           mpf_get_rdpe (rtmp, ftmp); */
/*           rdpe_set_dl (precision, 1, (long int) ((1 - 0.5 * s->mpwp) * LOG10_2 */
/*                                                  + rdpe_log10 (rtmp))); */
/*           if (rdpe_lt (s->drad[s->clust[j]], precision)) */
/*             { */

/*               MPS_DEBUG (s, "Separating root %d from the " */
/*                          "rest of the cluster n°%d", s->clust[j], i); */

/*               /\* Save a log of the detachement in s->clust_detach */
/*                * to make checking if the root is really outside the */
/*                * cluster after the computing of Newton polygonal. */
/*                * We are creating a new cluster in i+1, moving other */
/*                * cluster ahead, so first move ahead the vector after */
/*                * i+1, and then set the i+1 position to i. */
/*                * In theory we should check if s->clust_detached[j] > i */
/*                * and in that case shift it to s->clust_detached[j] + 1, */
/*                * but that's not possible becase cluster after this are not */
/*                * yet analized. */
/*                *\/ */
/*               for (ind = i + 1; ind < s->nclust; ind++) */
/*                 { */
/*                   s->clust_detached[ind + 1] = s->clust_detached[ind]; */
/*                 } */
/*               s->clust_detached[i + 1] = i; */


/*               /\* Move other roots back in the cluster *\/ */
/*               n_aux = s->clust[j]; */
/*               for (ind = j + 1; ind < s->punt[i + 1]; ind++) */
/*                 { */
/*                   s->clust[ind - 1] = s->clust[ind]; */
/*                 } */

/*               s->clust[s->punt[i + 1] - 1] = n_aux; */
/*               s->punt[i + 1]--; */

/*               /\* Move ahead s->punt *\/ */
/*               for (ind = s->nclust; ind > i + 1; ind--) */
/*                 { */
/*                   s->punt[ind + 1] = s->punt[ind]; */
/*                 } */

/*               /\* Set s->punt *\/ */
/*               s->punt[i + 2] = s->punt[i + 1] + 1; */
/* 	      s->nclust++; */

/*               /\* Start from the next root, that is shifted one position back *\/ */
/*               j--; */

/*               /\* If the cluster is now a single element cluster, let's */
/*                  skip to the next one *\/ */
/*               if (s->punt[i + 1] - s->punt[i] == 1) */
/*                 { */
/*                   break; */
/*                 } */
/*             } */
/*         } */
/*     } */

/*   mpf_clear (ftmp); */

/*   if (s->debug_level & MPS_DEBUG_CLUSTER) */
/*     { */
/*       MPS_DEBUG (s, "Debugging cluster structure after root detaching"); */
/*       mps_debug_cluster_structure (s); */
/*     } */


/* } */


/* void */
/* mps_cluster_reassemble (mps_status * s, int i_clust) */
/* { */
/*   MPS_DEBUG_THIS_CALL; */

/*   int i, l, j; */

/*   if (s->debug_level & MPS_DEBUG_CLUSTER) */
/*     { */
/*       MPS_DEBUG (s, "Debugging cluster structure before reassembling the original one"); */
/*       mps_debug_cluster_structure (s); */
/*     } */

/*   if (i_clust == MPS_ALL_CLUSTERS) */
/*     for (j = 0; j < s->nclust; j++) */
/*       { */
/*         mps_cluster_reassemble (s, j); */
/*         return; */
/*       } */

/*   MPS_DEBUG (s, "Reassembling cluster %d", i_clust); */

/*   for (i = 0; i < s->nclust; i++) */
/*     { */

/*       if (s->clust_detached[i] == i_clust) */
/*         { */
/*           MPS_DEBUG (s, "Recompacting cluster %d and %d", i_clust, i); */

/*           /\* We need this to be true to make the reassembling of */
/*            * the cluster work as expected *\/ */
/*           assert (i_clust < i); */

/*           l = s->clust[s->punt[i]]; */
/*           for (j = s->punt[i_clust + 1]; j < s->punt[i]; j++) */
/*             { */
/*               s->clust[j + 1] = s->clust[j]; */
/*             } */
/*           s->clust[s->punt[i_clust + 1]] = l; */
/*           for (j = i_clust + 1; j <= i; j++) */
/*             { */
/*               s->punt[j]++; */
/*             } */

/*           assert (s->punt[i] == s->punt[i + 1]); */

/*           for (j = i + 1; j < s->nclust; j++) */
/*             { */
/*               s->punt[j - 1] = s->punt[j]; */
/*               s->clust_detached[j - 1] = s->clust_detached[j]; */
/*             } */

/*           s->punt[s->nclust - 1] = s->punt[s->nclust]; */
/*           s->nclust--; */
/*           for (j = 0; j < s->nclust; j++) */
/*             { */
/*               if (s->clust_detached[j] > i) */
/*                 { */
/*                   s->clust_detached[j]--; */
/*                 } */
/*             } */

/* 	  i--; */
/*         } */
/*     } */

/*   if (s->debug_level & MPS_DEBUG_CLUSTER) */
/*     { */
/*       MPS_DEBUG (s, "Debugging cluster structure after restoring the original one"); */
/*       mps_debug_cluster_structure (s); */
/*     } */
/* } */
