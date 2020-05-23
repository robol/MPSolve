/*
 * This file is part of MPSolve 3.1.9
 *
 * Copyright (C) 2001-2020, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Leonardo Robol <leonardo.robol@unipi.it>
 */

/**
 * @file
 *
 * Implementation of the routines that compute the Aberth correction starting from the
 * Newton's one.
 *
 * @brief Implementation of Aberth correction computation.
 */

#ifndef MPS_ABERTH_H_
#define MPS_ABERTH_H_

#include <mps/mps.h>

void mps_faberth (mps_context * s, mps_approximation * root, cplx_t abcorr);
void mps_daberth (mps_context * s, mps_approximation * root, cdpe_t abcorr);
void mps_maberth (mps_context * s, mps_approximation * root, mpc_t abcorr);
void mps_faberth_s (mps_context * s, mps_approximation * root, mps_cluster * cluster, cplx_t abcorr);
void mps_faberth_wl (mps_context * s, int j, cplx_t abcorr, pthread_mutex_t * aberth_mutexes);
void mps_daberth_s (mps_context * s, mps_approximation * root, mps_cluster * cluster, cdpe_t abcorr);
void mps_daberth_wl (mps_context * s, int j, cdpe_t abcorr, pthread_mutex_t * aberth_mutexes);
void mps_maberth_s (mps_context * s, mps_approximation * root, mps_cluster * cluster, mpc_t abcorr);
void mps_maberth_s_wl (mps_context * s, int j, mps_cluster * cluster, mpc_t abcorr,
                       pthread_mutex_t * aberth_mutex);
void mps_mnewtis (mps_context * s);

#endif /* endif MPS_ABERTH_H_ */
