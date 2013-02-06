/**
 * @file
 * @brief Data structures for cluster analysis.
 */

#ifndef __MPS_CLUSTER
#define __MPS_CLUSTER

#include <mps/mps.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifdef _MPS_PRIVATE
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
};

/**
 * @brief Cluster hold in a mps_clusterization.
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
#endif /* #ifdef _MPS_PRIVATE */

/*********************************************************************************
 *                                   FUNCTIONS                                   *
 *********************************************************************************/
 
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

#ifdef __cplusplus
}
#endif

#endif
