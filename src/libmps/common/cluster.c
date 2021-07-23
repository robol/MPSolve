/*
 * This file is part of MPSolve 3.2.1
 *
 * Copyright (C) 2001-2020, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Dario Andrea Bini <bini@dm.unipi.it>
 *   Giuseppe Fiorentino <fiorent@dm.unipi.it>
 *   Leonardo Robol <leonardo.robol@unipi.it>
 */


#include <string.h>
#include <float.h>
#include <assert.h>
#include <mps/mps.h>
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
  pthread_mutex_init (&cluster->lock, NULL);

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

  pthread_mutex_init (&cluster->lock, NULL);

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

  pthread_mutex_destroy(&cluster->lock);

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
  root->prev = NULL;

  root->next = cluster->first;
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

void
mps_clusterization_detach_clusters (mps_context * s, mps_clusterization * c)
{
  MPS_DEBUG_THIS_CALL (s);

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
          mpcf_rmod (rtmp, s->root[k]->mvalue);

          /* We need a complex condition here since the heuristic used to determine if a root
           * is a simple root in a cluster is based on Newton radii.
           * These have different behavious based on the algorithm that has been selected, so
           * we introduce here two different guesses that work in each one. */
          if (((s->algorithm == MPS_ALGORITHM_STANDARD_MPSOLVE) &&
               ((rdpe_Esp (rtmp) - rdpe_Esp (s->root[k]->drad) > s->mpwp / sqrt (item->cluster->n) + 1) ||
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
  MPS_DEBUG_THIS_CALL (s);

  mps_cluster_item * cluster;

  cluster = s->clusterization->first;
  while (cluster != NULL)
    {
      mps_cluster_item * next = cluster->next;

      if (cluster->detached)
        {
          mps_cluster_insert_root (s, cluster->detached->cluster, cluster->cluster->first->k);
          mps_clusterization_remove_cluster (s, s->clusterization, cluster);
        }

      cluster = next;
    }
}

void
mps_cluster_detachment_reset (mps_context * s)
{
  mps_clusterization_reassemble_clusters (s, s->clusterization);
}
