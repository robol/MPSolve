/*
 * This file is part of MPSolve 3.1.8
 *
 * Copyright (C) 2001-2019, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Leonardo Robol <leonardo.robol@sns.it>
 */

/**
 * @file
 *
 * @brief Implementation of the convex hull computation.
 */

#ifndef MPS_CONVEX_H_
#define MPS_CONVEX_H_

#include <mps/mps.h>

/**
 * @brief Generic vertex of a linear hypograph.
 */
struct mps_vertex {
  /**
   * @brief The x coordinate of the vertex.
   */
  int x;

  /**
   * @brief The y coordinate of the vertex.
   */
  double y;

  /**
   * @brief A pointer to the next vertex in the hypograph, or
   * NULL if this is the last vertex or a detached one.
   */
  struct mps_vertex * next;

  /**
   * @brief A pointer to the previous vertex in the hypograph,
   * or NULL if this is the first vertex or a detached one.
   */
  struct mps_vertex * previous;
};

typedef struct mps_vertex mps_vertex;

/**
 * @brief A set described as hypograph of a piecewise linear function.
 *
 * The explicit description of the set is given by a set of vertexes of the
 * type \f$(i, y_i)\f$ where \f$i\f$ is a positive integer.
 */
struct mps_linear_hypograph {
  int n;

  mps_vertex * last;

  mps_vertex * first;
};

typedef struct mps_linear_hypograph mps_linear_hypograph;

mps_linear_hypograph * mps_convex_hull (mps_context * s, mps_linear_hypograph * l);

int * mps_fconvex (mps_context * s, int n, double a[]);

mps_linear_hypograph * mps_linear_hypograph_new (mps_context * ctx);

void mps_linear_hypograph_free (mps_context * ctx, mps_linear_hypograph * l);

#endif /* endif MPS_CONVEX_H_ */


