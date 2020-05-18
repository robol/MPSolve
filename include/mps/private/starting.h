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
 * @brief Selection of starting points and shifting of the polynomials
 * to zoom in the clusters.
 */

#ifndef MPS_STARTING_H_
#define MPS_STARTING_H_

MPS_BEGIN_DECLS

/* functions in starting.c */
void mps_fstart (mps_context * s, int n, mps_cluster_item * cluster, double clust_rad,
                 double g, rdpe_t eps_out, double fap[]);
void mps_dstart (mps_context * s, int n, mps_cluster_item * cluster, rdpe_t clust_rad,
                 rdpe_t g, rdpe_t eps_out, rdpe_t dap[]);
void mps_mstart (mps_context * s, int n, mps_cluster_item * cluster, rdpe_t clust_rad,
                 rdpe_t g, rdpe_t dap[], mpcf_t gg);
void mps_frestart (mps_context * s);
void mps_drestart (mps_context * s);
void mps_mrestart (mps_context * s);
void mps_fshift (mps_context * s, int m, mps_cluster_item * cluster, double clust_rad,
                 cplx_t g, rdpe_t eps);
void mps_dshift (mps_context * s, int m, mps_cluster_item * cluster, rdpe_t clust_rad,
                 cdpe_t g, rdpe_t eps);
void mps_mshift (mps_context * s, int m, mps_cluster_item * cluster, rdpe_t clust_rad,
                 mpcf_t g);

/* functions in recursive-starting.c */
void mps_recursive_fstart (mps_context * ctx, mps_polynomial * poly, mps_approximation ** approximations);
void mps_recursive_dstart (mps_context * ctx, mps_polynomial * poly, mps_approximation ** approximations);
void mps_recursive_mstart (mps_context * ctx, mps_polynomial * poly, mps_approximation ** approximations);

/* functions in file-starting.c */
void mps_file_fstart (mps_context * ctx, mps_polynomial * poly, mps_approximation ** approximations);
void mps_file_dstart (mps_context * ctx, mps_polynomial * poly, mps_approximation ** approximations);
void mps_file_mstart (mps_context * ctx, mps_polynomial * poly, mps_approximation ** approximations);

MPS_END_DECLS

#endif /* MPS_STARTING_H_ */

