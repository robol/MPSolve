/*
 * This file is part of MPSolve 3.1.5
 *
 * Copyright (C) 2001-2015, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Leonardo Robol <leonardo.robol@sns.it>
 */

/**
 * @file
 *
 * @brief Data structures for cluster analysis and some accessors and
 * internal functions.
 */

#ifndef MPS_CLUSTER_H_
#define MPS_CLUSTER_H_

#define MPS_ALL_CLUSTERS -1

#include <mps/mps.h>

MPS_BEGIN_DECLS

/**
 * @brief This struct represent a root inside of a <code>mps_cluster</code>.
 */
struct mps_root {
  /**
   * @brief Index of the root that is considered.
   */
  long int k;

  /**
   * @brief Next root, or <code>NULL</code> if this is the last root of the cluster.
   */
  mps_root * next;

  /**
   * @brief Pointer to the previous root, or <code>NULL</code> if there is no previous
   * root.
   */
  mps_root * prev;
};

/**
 * @brief A cluster of <code>mps_roots</code>.
 */
struct mps_cluster {
  /**
   * @brief Number of roots in the cluster.
   */
  long int n;

  /**
   * @brief Pointer to the first root in the cluster.
   */
  mps_root * first;
  
  /**
   * @brief Internal mutex used to perform operations in a thread-safe
   * way. 
   */
  pthread_mutex_t lock;
};

/**
 * @brief Cluster held in a mps_clusterization.
 */
struct mps_cluster_item {
  /**
   * @brief Pointer to the actual cluster.
   */
  mps_cluster * cluster;

  /**
   * @brief Next cluster in the clusterization, or NULL if there is no such
   * cluster.
   */
  mps_cluster_item * next;

  /**
   * @brief Previous cluster in the clusterizaion or NULL if there is no
   * such cluster.
   */
  mps_cluster_item * prev;

  /**
   * @brief This pointer to the cluster from which the cluster were
   * detached, if any. Otherwise it is set to NULL.
   */
  mps_cluster_item * detached;
};

/**
 * @brief A list of <code>mps_cluster</code>.
 */
struct mps_clusterization {
  /**
   * @brief Number of cluster in the clusterization.
   */
  long int n;

  /**
   * @brief Pointer to the first cluster in the clusterization.
   */
  mps_cluster_item * first;
};

/*********************************************************************************
*                                   FUNCTIONS                                   *
*********************************************************************************/

void mps_cluster_reset (mps_context * s);
void mps_fcluster (mps_context * s, double * frad, int nf);
void mps_dcluster (mps_context * s, rdpe_t * drad, int nf);
void mps_mcluster (mps_context * s, rdpe_t * drad, int nf);
void mps_debug_cluster_structure (mps_context * s);
void mps_cluster_analysis (mps_context * ctx, mps_polynomial * p);

/* Functions for mps_cluster */
mps_cluster * mps_cluster_empty (mps_context * s);
mps_cluster * mps_cluster_with_root (mps_context * s, long int root_index);
void mps_cluster_free (mps_context * s, mps_cluster * cluster);
mps_root * mps_cluster_insert_root (mps_context * s, mps_cluster * cluster, long int root_index);
void mps_cluster_remove_root (mps_context * s, mps_cluster * cluster, mps_root * root);
mps_cluster * mps_cluster_join (mps_context * s, mps_cluster * cluster_a, mps_cluster * cluster_b);

/* Functions for mps_clusterization */
mps_clusterization * mps_clusterization_empty (mps_context * s);
mps_cluster_item * mps_clusterization_insert_cluster (mps_context * s, mps_clusterization * c, mps_cluster * cluster);
void mps_clusterization_pop_cluster (mps_context * s, mps_clusterization * c, mps_cluster_item * cluster_item);
void mps_clusterization_remove_cluster (mps_context * s, mps_clusterization * c, mps_cluster_item * cluster_item);
void mps_clusterization_free (mps_context * s, mps_clusterization * c);
void mps_clusterization_detach_clusters (mps_context * s, mps_clusterization * c);
void mps_clusterization_reassemble_clusters (mps_context * s, mps_clusterization * c);

void mps_cluster_detach (mps_context * s, mps_cluster * cluster);

MPS_END_DECLS

#endif /* endif MPS_CLUSTER_H_ */
