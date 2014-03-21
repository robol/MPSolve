/*
 * This file is part of MPSolve 3.1.5
 *
 * Copyright (C) 2001-2014, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Leonardo Robol <leonardo.robol@sns.it>
 */

/**
 * @file
 * @brief
 */

#ifndef MPS_STARTING_CONFIGURATION_H_
#define MPS_STARTING_CONFIGURATION_H_

MPS_BEGIN_DECLS

#define MPS_STARTING_CONFIGURATION_INIT { \
    .n_radii = 0,                           \
    .fradii = NULL,                         \
    .dradii = NULL,                         \
    .partitioning = NULL                    \
}

/**
 * @brief This struct holds the information about the starting disposal of the
 * approximations that has been obtained by the computation of the Newton polygon
 * of the polynomial.
 */
struct mps_starting_configuration {
  /**
   * @brief The number of circles on which the approximations should be
   * put.
   */
  int n_radii;

  /**
   * @brief An array of integers that holds the indexes of the approximations that
   * should be put on the different circles.
   *
   * Precisely, this array represent a partition of \f$[0, n]\f$.
   */
  int * partitioning;

  /**
   * @brief An array containing the radius of the circles where the roots should be
   * placed.
   */
  double * fradii;

  /**
   * @brief The DPE version of fradii.
   */
  rdpe_t * dradii;
};

typedef struct mps_starting_configuration mps_starting_configuration;

/**
 * @brief Clear all the storage that has been allocated inside a starting
 * configuration.
 *
 * @param ctx The current mps_context.
 * @param c The mps_starting_configuration that should be cleared.
 */
void mps_starting_configuration_clear (mps_context * ctx, mps_starting_configuration * c);

mps_starting_configuration mps_fcompute_starting_radii (mps_context * s, int n,
                                                        mps_cluster_item * cluster_item,
                                                        double clust_rad, double g, rdpe_t eps,
                                                        double fap[]);


MPS_END_DECLS

#endif /* MPS_STARTING_CONFIGURATION_H_ */
