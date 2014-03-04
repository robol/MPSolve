/*
 * This file is part of MPSolve 3.1.5
 *
 * Copyright (C) 2001-2014, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Leonardo Robol <robol@mail.dm.unipi.it>
 */

/**
 * @file 
 * 
 * @brief Implementation of the convex hull computation.
 */

#ifndef MPS_CONVEX_H_
#define MPS_CONVEX_H_

#include <mps/mps.h>

void mps_fconvex (mps_context * s, int n, double a[]);

void mps_fcompute_starting_radii (mps_context * s, int n, mps_cluster_item * cluster_item,
				  double clust_rad, double g, rdpe_t eps,
				  double fap[]);


#endif /* endif MPS_CONVEX_H_ */


