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


#include <string.h>
#include <float.h>
#include <assert.h>
#include <mps/gmptools.h>
#include <mps/mps.h>
#include <mps/cluster.h>
#include <math.h>

/**
 * @brief Get an empty mps_cluster, with no roots.
 * @param s The <code>mps_context</code> of the current computation.
 */
mps_cluster *
mps_cluster_empty (mps_context * s)
{
  mps_cluster * cluster = mps_new (mps_cluster);
  cluster->first = NULL;
  cluster->n = 0;
  return cluster;
}

/**
 * @brief Create a cluster containing only the selected root.
 *
 * @param s The <code>mps_context</code> of the current computation.
 * @param root_index The root that must be in the cluster.
 */
mps_cluster * 
mps_cluster_with_root (mps_context * s, long int root_index)
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
 * @param s The <code>mps_context</code> of the current computation.
 * @param cluster The cluster to free.
 */
void 
mps_cluster_free (mps_context * s, mps_cluster * cluster)
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
 * @param s The <code>mps_context</code> of the current computation.
 * @param cluster The cluster in which the root must be inserted.
 * @param root_index The index of the root to insert.
 */
mps_root *
mps_cluster_insert_root (mps_context * s, 
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
 * @param s The <code>mps_context</code> of the current computation.
 * @param cluster The cluster from which the root must be removed.
 * @param root The root to remove.
 *
 * Please note the the root specified must be in the cluster, otherwise
 * an assertion error or segmentation fault will be triggered.
 */
void 
mps_cluster_remove_root (mps_context * s, mps_cluster * cluster, mps_root * root)
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
 * @param s The <code>mps_context</code> of the current computation.
 * @param cluster_a The first cluster
 * @param cluster_b The second cluster
 * @return A new cluster containing the roots of both.
 */
mps_cluster *
mps_cluster_join (mps_context * s, mps_cluster * cluster_a, mps_cluster * cluster_b)
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
 *
 * @param s The <code>mps_context</code> of the current computation.
 */
mps_clusterization *
mps_clusterization_empty (mps_context * s)
{
  mps_clusterization * c = mps_new (mps_clusterization);
  c->n = 0;
  c->first = NULL;
  return c;
}

/**
 * @brief Insert a new cluster into a root clusterization.
 * @param s The <code>mps_context</code> of the current computation.
 * @param c The clusterization in which the cluster should be inserted.
 * @param cluster The cluster that should be inserted.
 */
mps_cluster_item *
mps_clusterization_insert_cluster (mps_context * s, mps_clusterization * c, mps_cluster * cluster)
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
 * @param s The <code>mps_context</code> of the current computation.
 * @param c The clusterization from which the cluster_item should be popped.
 * @param cluster_item The cluster item to remove.
 */
void
mps_clusterization_pop_cluster (mps_context * s, mps_clusterization * c, mps_cluster_item * cluster_item)
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
 * @param s The <code>mps_context</code> of the current computation.
 * @param c The clusterization from where the cluster_item should be removed.
 * @param cluster_item The cluster item to remove.
 */
void 
mps_clusterization_remove_cluster (mps_context * s, mps_clusterization * c, mps_cluster_item * cluster_item)
{
  mps_clusterization_pop_cluster (s, c, cluster_item);
  mps_cluster_free (s, cluster_item->cluster);
  free (cluster_item);
}

/**
 * @brief Free a clusterization and all the cluster in it.
 * @param s The <code>mps_context</code> of the current computation.
 * @param c The clusterization to free.
 */
void
mps_clusterization_free (mps_context * s, mps_clusterization * c)
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
 * @param s the mps_context pointer.
 */
void
mps_cluster_reset (mps_context * s)
{
  /* Reset cluster status of the roots */ 
  int i;
  mps_cluster * cluster;

  for (i = 0; i < s->n; i++) 
    { 
      s->root[i]->status = MPS_ROOT_STATUS_CLUSTERED;
      s->root[i]->attrs = MPS_ROOT_ATTRS_NONE;
      s->root[i]->inclusion = MPS_ROOT_INCLUSION_UNKNOWN;
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
 * @param s  The <code>mps_context</code> associated with the current
 *           computaion.
 * @param frad The vector of radii to use for cluster analysis.
 * @param nf see above for a detailed description.
 */
void
mps_fcluster (mps_context * s, double * frad, int nf)
{
  s->operation = MPS_OPERATION_CLUSTER_ANALYSIS;

  /* We need to scan every cluster and make it in pieces, if possible */
  mps_clusterization * new_clusterization = mps_clusterization_empty (s);
  mps_cluster_item * item;
  int analyzed_roots = 0;
  int i, j;

  /* This value is set to false if the radius are not newton isolated
   * by means of the newton radii. */
  mps_boolean newton_isolation = true;

  /* Debug clusterization status if debugging was required */
  if (s->debug_level & MPS_DEBUG_CLUSTER)
    {
      int i;
      MPS_DEBUG (s, "Debugging the radius and approximations obtained for the roots before cluster analysis");
      for (i = 0; i < s->n; i++)
        {
          MPS_DEBUG_CPLX (s, s->root[i]->fvalue, "Root %d", i);
          MPS_DEBUG (s, "radius for root %4d: %e", i, frad[i]);
        }

      MPS_DEBUG (s, "Debugging cluster structure before cluster analysis");
      mps_debug_cluster_structure (s);
    }

  /*
   * Mark newton isolated roots as newton isolated.
   */
  double * newton_radii = mps_newv (double, s->n);
  for (i = 0; i < s->n; i++)
    newton_radii[i] = s->root[i]->frad;
  
  for (i = 0; i < s->n; i++)
    {
      for (j = 0; j < s->n; j++)
        {
          if ((i != j) && mps_ftouchnwt (s, newton_radii, nf, i, j))
            {
              newton_isolation = false;
              break;
            }
          /* s->root_status[i] = MPS_ROOT_STATUS_NEWTON_ISOLATED; */
        }
    }

  free (newton_radii);

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

  /* Now do cluster analysis with the rest of the clusters, but only if
   * newton isolation is not guaranteed by means of the newton radii. */
  if (!newton_isolation)
    {
      /* if (MPS_INPUT_CONFIG_IS_USER (s->input_config))  */
      /*        {  */
      /*          mps_clusterization_free (s, new_clusterization);  */
      /*          return; */
      /*        } */

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
              /* mps_cluster_item * c_item; */
              mps_root * iter_root;

              /* for (c_item = s->clusterization->first; c_item != NULL; c_item = c_item->next) */
                {
                  /* mps_cluster * iter_cluster = c_item->cluster; */
                  mps_cluster * iter_cluster = cluster;

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

          /* Check if the new cluster is isolated and, in that case, set the gerschgorin
           * radius as inclusion radius if it's more conveniente than the old one.
           * In general the Gerschgorin radius cannot be used as inclusion radius, because
           * it may touch another radius and so it may be empty. */
            if (new_cluster->n == 1)  
              {  
                 int k = new_cluster->first->k;   
                 double new_rad;   

                 new_rad = cplx_mod (s->root[k]->fvalue) * 4.0f * DBL_EPSILON + frad[k];    

                 /* Check if the computed radius is more convenient than the old one.  */  
                 /*      If that's the case, apply it as inclusion radius   */  
                 if (new_rad < s->root[k]->frad)      
                   s->root[k]->frad = new_rad;      
              }  
        }

    }

  if (newton_isolation)
    {
      mps_clusterization_free (s, new_clusterization);
      new_clusterization = mps_clusterization_empty (s);

      if (s->debug_level & MPS_DEBUG_CLUSTER)
        MPS_DEBUG (s, "Reached isolation using Newton radii, so skipping every other check with Gerschgorin");

      for (i = 0; i < s->n; i++)
        {
          mps_clusterization_insert_cluster (s, new_clusterization, 
                                             mps_cluster_with_root (s, i));       
        }
    }

  /* Set the new clusterizaition in the mps_context */
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
mps_dcluster (mps_context * s, rdpe_t * drad, int nf)
{
  s->operation = MPS_OPERATION_CLUSTER_ANALYSIS;

  /* We need to scan every cluster and make it in pieces, if possible */
  mps_clusterization * new_clusterization = mps_clusterization_empty (s);
  mps_cluster_item * item;
  int i, j;

  /* This value is set to false if the radius are not newton isolated
   * by means of the newton radii. */
  mps_boolean newton_isolation = true;

  /* Debug clusterization status if debugging was required */
  if (s->debug_level & MPS_DEBUG_CLUSTER)
    {
      int i;
      MPS_DEBUG (s, "Debugging the radius and approximations obtained for the roots before cluster analysis");
      for (i = 0; i < s->n; i++)
        {
          MPS_DEBUG_CDPE (s, s->root[i]->dvalue, "Root %d", i);
          MPS_DEBUG_RDPE (s, drad[i], "radius for root %4d", i);
        }

      MPS_DEBUG (s, "Debugging cluster structure before cluster analysis");
      mps_debug_cluster_structure (s);
    }

  int analyzed_roots = 0;

  /* Do a first check of clusterization using the newton
   * radii. These are not valid to perform cluster analysis in
   * general, but can be used if they provide *COMPLETE* Newton
   * isolation. */  
  rdpe_t * newton_radii = rdpe_valloc (s->n);
  for (i = 0; i < s->n; i++)
    rdpe_set (newton_radii[i], s->root[i]->drad);
  
  for (i = 0; i < s->n; i++)
    {
      for (j = 0; j < s->n; j++)
        {
          if ((i != j) && mps_dtouchnwt (s, newton_radii, nf, i, j))
            {
              newton_isolation = false;
              break;
            }
        }
    }

  rdpe_vfree (newton_radii);

  /* If newton isolation has not been reached check with Gerschgorin */
    {
      /* if (MPS_INPUT_CONFIG_IS_USER (s->input_config))  */
      /*        {  */
      /*          mps_clusterization_free (s, new_clusterization);  */
      /*          return;  */
      /*        } */

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
              /* mps_cluster_item * c_item; */
              mps_root * iter_root;

              /* for (c_item = s->clusterization->first; c_item != NULL; c_item = c_item->next) */
                {
                  /* mps_cluster * iter_cluster = c_item->cluster; */
                  mps_cluster * iter_cluster = cluster;

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

          /* Check if the new cluster is isolated and, in that case, set the gerschgorin
           * radius as inclusion radius if it's more conveniente than the old one.
           * In general the Gerschgorin radius cannot be used as inclusion radius, because
           * it may touch another radius and so it may be empty. */
            if (new_cluster->n == 1)  
              {  
                int k = new_cluster->first->k;  
                rdpe_t new_rad;  
          
                cdpe_mod (new_rad, s->root[k]->dvalue);  
                rdpe_mul_eq_d (new_rad, 4 * DBL_EPSILON);  
                rdpe_add_eq (new_rad, drad[k]);  


                /* Check if the computed radius is more convenient than the old one.  
                 If that's the case, apply it as inclusion radius */  
                if (rdpe_lt (new_rad, s->root[k]->drad))    
                rdpe_set (s->root[k]->drad, new_rad);    
              }  
        }
    }


  if (newton_isolation)
    {
      mps_clusterization_free (s, new_clusterization);
      new_clusterization = mps_clusterization_empty (s);

      if (s->debug_level & MPS_DEBUG_CLUSTER)
        MPS_DEBUG (s, "Reached isolation using Newton radii, so skipping every other check with Gerschgorin");

      for (i = 0; i < s->n; i++)
        {
          mps_clusterization_insert_cluster (s, new_clusterization, 
                                             mps_cluster_with_root (s, i));
        }
    }

  /* Set the new clusterizaition in the mps_context */
  mps_clusterization_free (s, s->clusterization);
  s->clusterization = new_clusterization;

  if (s->debug_level & MPS_DEBUG_CLUSTER)
    {
      MPS_DEBUG (s, "Debugging cluster structure after cluster analysis");
      mps_debug_cluster_structure (s);
    }
}


void
mps_debug_cluster_structure (mps_context * s)
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
mps_mcluster (mps_context * s, rdpe_t * drad, int nf)
{
  MPS_DEBUG_THIS_CALL;

  s->operation = MPS_OPERATION_CLUSTER_ANALYSIS;

  /* We need to scan every cluster and make it in pieces, if possible */
  mps_clusterization * new_clusterization = mps_clusterization_empty (s);
  mps_cluster_item * item;
  int i, j;

  /* This value is set to false if the radius are not newton isolated
   * by means of the newton radii. */
  mps_boolean newton_isolation = true;

  /* Debug clusterization status if debugging was required */
  if (s->debug_level & MPS_DEBUG_CLUSTER)
    {
      int i;
      MPS_DEBUG (s, "Debugging the radius and approximations obtained for the roots before cluster analysis");
      for (i = 0; i < s->n; i++)
        {
          MPS_DEBUG_MPC (s, mpc_get_prec (s->root[i]->mvalue), s->root[i]->mvalue, "Root %d", i);
          MPS_DEBUG_RDPE (s, drad[i], "radius for root %4d", i);
        }

      MPS_DEBUG (s, "Debugging cluster structure before cluster analysis");
      mps_debug_cluster_structure (s);
    }

  int analyzed_roots = 0;

  /* Do a first check of clusterization using the newton
   * radii. These are not valid to perform cluster analysis in
   * general, but can be used if they provide *COMPLETE* Newton
   * isolation. */  
  rdpe_t * newton_radii = rdpe_valloc (s->n);
  for (i = 0; i < s->n; i++)
    rdpe_set (newton_radii[i], s->root[i]->drad);

  for (i = 0; i < s->n; i++)
    {
      for (j = 0; j < s->n; j++)
        {
          if ((i != j) && mps_mtouchnwt (s, newton_radii, nf, i, j))
            {
              if (s->debug_level & MPS_DEBUG_CLUSTER)
                MPS_DEBUG (s, "Failing newton isolation on root %d and %d", i, j);

              newton_isolation = false;
              break;
            }
          /* s->root_status[i] = MPS_ROOT_STATUS_NEWTON_ISOLATED; */
        }
    }

  rdpe_vfree (newton_radii);

  /* If newton isolation is not reached with Newton use Gerschgorin */
    {
      /* if (MPS_INPUT_CONFIG_IS_USER (s->input_config))  */
      /*        {  */
      /*          mps_clusterization_free (s, new_clusterization);  */
      /*          return;  */
      /*        } */
      
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
              /* mps_cluster_item * c_item; */
              mps_root * iter_root;

              /* for (c_item = s->clusterization->first; c_item != NULL; c_item = c_item->next) */
              { 
                /*        mps_cluster * iter_cluster = c_item->cluster; */
                mps_cluster * iter_cluster = iter_cluster = cluster;

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

          /* Check if the new cluster is isolated and, in that case, set the gerschgorin
           * radius as inclusion radius if it's more conveniente than the old one.
           * In general the Gerschgorin radius cannot be used as inclusion radius, because
           * it may touch another radius and so it may be empty. */
          if (new_cluster->n == 1)  
            {  
              int k = new_cluster->first->k;  
              cdpe_t c;  
              rdpe_t new_rad;  
              
              /* Check if the computed radius is more convenient than the old one.  
                 If that's the case, apply it as inclusion radius */  
              mpc_get_cdpe (c, s->root[k]->mvalue);  
              cdpe_mod (new_rad, c);  
              rdpe_mul_eq (new_rad, s->mp_epsilon);  
              rdpe_mul_eq_d (new_rad, 4.0f);  
              rdpe_add_eq (new_rad, drad[k]);  
              
              if (rdpe_lt (new_rad, s->root[k]->drad))      
                rdpe_set (s->root[k]->drad, new_rad);      
            } 
        }
    }

  if (newton_isolation)
    {
      mps_clusterization_free (s, new_clusterization);
      new_clusterization = mps_clusterization_empty (s);

      if (s->debug_level & MPS_DEBUG_CLUSTER)
        MPS_DEBUG (s, "Reached isolation using Newton radii, so skipping every other check with Gerschgorin");

      for (i = 0; i < s->n; i++)
        {
          mps_clusterization_insert_cluster (s, new_clusterization, 
                                             mps_cluster_with_root (s, i));
        }
    }

  /* Set the new clusterizaition in the mps_context */
  mps_clusterization_free (s, s->clusterization);
  s->clusterization = new_clusterization;


  if (s->debug_level & MPS_DEBUG_CLUSTER)
    {
      MPS_DEBUG (s, "Debugging cluster structure after cluster analysis");
      mps_debug_cluster_structure (s);
    }

  if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)
    {
      mps_dump (s);
    }
}


void 
mps_clusterization_detach_clusters (mps_context * s, mps_clusterization * c)
{
  MPS_DEBUG_THIS_CALL;

  /* Disable this function since it is not working as it should. */
  /* The problem, now, is that more than a root could be removed from
   * a cluster and will be then marked as isolated even if it is only isolated
   * from the base cluster and from the other detached roots. */
  return;

  mps_cluster_item * item;
  rdpe_t rtmp;
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
          mpc_rmod (rtmp, s->root[k]->mvalue);

          /* We need a complex condition here since the heuristic used to determine if a root
           * is a simple root in a cluster is based on Newton radii. 
           * These have different behavious based on the algorithm that has been selected, so
           * we introduce here two different guesses that work in each one. */
          if (((s->algorithm == MPS_ALGORITHM_STANDARD_MPSOLVE) && 
                ((rdpe_Esp (rtmp) - rdpe_Esp (s->root[k]->drad) > s->mpwp / sqrt(item->cluster->n) + 1) ||
                (s->root[k]->status == MPS_ROOT_STATUS_APPROXIMATED_IN_CLUSTER))) ||
              ((s->algorithm == MPS_ALGORITHM_SECULAR_GA) &&
                (rdpe_Esp (rtmp) - rdpe_Esp (s->root[k]->drad) > s->mpwp - 4)))
            {
              if (s->debug_level & MPS_DEBUG_CLUSTER)
                {
                  MPS_DEBUG (s, "Temporary removing root %d from its cluster since it is quasi approximated", k);
                }

              mps_cluster * detached_cluster = mps_cluster_with_root (s, k);
              mps_root * next_root = root->next;
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
mps_clusterization_reassemble_clusters (mps_context * s, mps_clusterization * c)
{
  MPS_DEBUG_THIS_CALL;

  mps_cluster_item * cluster;
  
  cluster = s->clusterization->first;
  while (cluster != NULL)
    {
      if (cluster->detached)
        {
          mps_cluster_insert_root (s, cluster->detached->cluster, cluster->cluster->first->k);
          mps_clusterization_remove_cluster (s, s->clusterization, cluster);
        }

      cluster = cluster->next;
    }
}

void
mps_cluster_detachment_reset (mps_context * s)
{
  mps_clusterization_reassemble_clusters (s, s->clusterization);
}
