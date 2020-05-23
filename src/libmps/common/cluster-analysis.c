/*
 * This file is part of MPSolve 3.1.9
 *
 * Copyright (C) 2001-2020, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Leonardo Robol <leonardo.robol@unipi.it>
 */

#include <mps/mps.h>
#include <string.h>

MPS_PRIVATE void
mps_cluster_analysis (mps_context * ctx, mps_polynomial * p)
{
  switch (ctx->lastphase)
    {
    case float_phase:
    {
      double * radii = double_valloc (ctx->n);

      mps_fradii (ctx, p, radii);
      mps_fcluster (ctx, radii, 2 * ctx->n);
      mps_fmodify (ctx, false);

      free (radii);
      break;
    }

    case dpe_phase:
    {
      rdpe_t * radii = rdpe_valloc (ctx->n);

      mps_dradii (ctx, p, radii);
      mps_dcluster (ctx, radii, 2 * ctx->n);
      mps_dmodify (ctx, false);

      free (radii);
      break;
    }

    case mp_phase:
    {
      rdpe_t * radii = rdpe_valloc (ctx->n);

      mps_mradii (ctx, p, radii);
      mps_mcluster (ctx, radii, 2 * ctx->n);
      mps_mmodify (ctx, false);

      free (radii);
      break;
    }

    case no_phase:
      break;
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
 * @param s  The <code>mps_context</code> associated with the current
 *           computaion.
 * @param frad The vector of radii to use for cluster analysis.
 * @param nf see above for a detailed description.
 */
MPS_PRIVATE void
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
MPS_PRIVATE void
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

struct _mps_cluster_worker_data {
  mps_context * ctx;
  mps_cluster * cluster;
  int * analyzed_roots;
  int base_root;
  int start_root;
  int end_root;
  rdpe_t * drad;
  int nf;
  pthread_mutex_t * block_mutex;
  mps_cluster ** original_clusters;
};

void *
_mps_mcluster_worker (void * data_ptr)
{
  struct _mps_cluster_worker_data *data = (struct _mps_cluster_worker_data*) data_ptr;
  int i, n = 0;
  mps_root * first = NULL;
  mps_root * last = NULL;
  mps_cluster * c = data->original_clusters[data->base_root];

  for (i = data->start_root; i < data->end_root; i++)
    {
      if (! data->analyzed_roots[i] && (data->original_clusters[i] == c))
	{
	  if (mps_mtouchnwt (data->ctx, data->drad, data->nf, data->base_root, i))
	    {
              if (! data->analyzed_roots[i])
                {
                  data->analyzed_roots[i] = true;

                  if (first == NULL)
                    {
                      last = first = mps_new (mps_root);
                      last->next = first->next = last->prev = first->prev = NULL;
                      last->k = i;
                    }
                  else
                    {
                      mps_root * new_root = mps_new (mps_root);
                      new_root->next = first;
                      first->prev = new_root;
                      first = new_root;
                      new_root->prev = NULL;
                      new_root->k = i;
                    }

                  n++;
                }
	    }
	}
    }

  if (n > 0)
    {
      pthread_mutex_lock (&data->cluster->lock);

      last->next = data->cluster->first;
      data->cluster->first->prev = last;
      data->cluster->first = first;
      data->cluster->n += n;

      pthread_mutex_unlock (&data->cluster->lock);
    }

  pthread_mutex_unlock (data->block_mutex);

  free (data);

  return NULL;
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
MPS_PRIVATE void
mps_mcluster (mps_context * s, rdpe_t * drad, int nf)
{
  MPS_DEBUG_THIS_CALL (s);

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

      if (! newton_isolation)
	break;
    }

  rdpe_vfree (newton_radii);

  /* Perform parallel analysis of the Gerschgorin disks. */
  int analyzed_roots = 0;
  int * already_analyzed_roots = mps_newv (int, s->n);
  mps_cluster ** original_clusters = mps_newv (mps_cluster*, s->n);
  mps_root * root = NULL;

  memset (already_analyzed_roots, 0, sizeof (int) * s->n);

  int block_size = 128;
  int block_number = (s->n - 1) / block_size + 1;
  pthread_mutex_t * block_mutexes = mps_newv (pthread_mutex_t, block_number);

  for (j = 0; j < block_number; j++)
    pthread_mutex_init (&block_mutexes[j], NULL);

  item = s->clusterization->first;
  while (item)
    {
      mps_cluster * c = item->cluster;
      mps_root * root = c->first;
      while (root)
        {
          original_clusters[root->k] = c;
          root = root->next;
        }
      item = item->next;
    }

  while (analyzed_roots < s->n)
    {
      int j;

      /* Find the first not analyzed root */
      j = 0;
      while (already_analyzed_roots[j])
	j++;
      // fprintf (stderr, "New base root = %d\n", j);

      if (j > s->n)
	break;

      item = mps_clusterization_insert_cluster (s, new_clusterization, mps_cluster_with_root (s, j));
      root = item->cluster->first;

      already_analyzed_roots[j] = true;

      do
	{
	  /* We need to check which other approximation touches our current base root
	   * and add it to our cluster. */
	  for (j = 0; j < block_number; j++)
	    {
	      struct _mps_cluster_worker_data * data = mps_new 
		(struct _mps_cluster_worker_data);

	      data->ctx = s;
	      data->cluster = item->cluster;
	      data->base_root = root->k;
	      data->start_root = j * block_size;
	      data->end_root = MIN ((j+1) * block_size, s->n);
	      data->analyzed_roots = already_analyzed_roots;
	      data->drad = drad;
	      data->nf = nf;
              data->block_mutex = &block_mutexes[j];
              data->original_clusters = original_clusters;

	      pthread_mutex_lock (data->block_mutex);
	      mps_thread_pool_assign (s, s->pool, _mps_mcluster_worker, data);
	    }

          sched_yield();

	  /* This is a cheap way of getting the next root in the cluster. 
	   * If it's already present (root->prev != NULL) just use it, otherwise
	   * lock the mutexes that the clusters hold while analyzing a single
	   * block. Getting lock implies that the cluster will be already
	   * extended, thus root->prev != NULL unless there were no other
	   * approximations inside the cluster in that block. In that case, try the
	   * next until they are all analyzed. */
	  for (j = 0; root->prev == NULL && j < block_number; j++)
	    {
	      pthread_mutex_lock (&block_mutexes[j]);
	      pthread_mutex_unlock (&block_mutexes[j]);
	    }

	  analyzed_roots++;

	} while ((root = root->prev) != NULL);
    }

  free (block_mutexes);
  free (already_analyzed_roots);
  free (original_clusters);
  mps_clusterization_free (s, s->clusterization);
  s->clusterization = new_clusterization;

  for (item = s->clusterization->first; item != NULL; item = item->next)
    {
      mps_cluster * new_cluster = item->cluster;      

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
            {
	      rdpe_set (s->root[k]->drad, new_rad);
              s->root[k]->frad = rdpe_get_d (s->root[k]->drad);
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

      s->clusterization = new_clusterization;
    }

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
