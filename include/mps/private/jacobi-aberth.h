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
 * @brief Implementation of the iterations using Jacobi-style updates.
 */

#ifndef MPS_JACOBI_ABERTH_H_
#define MPS_JACOBI_ABERTH_H_

#include <mps/mps.h>

MPS_BEGIN_DECLS

int mps_faberth_packet (mps_context * ctx, mps_polynomial * p, mps_boolean just_regenerated);
int mps_daberth_packet (mps_context * ctx, mps_polynomial * p, mps_boolean just_regenerated);
int mps_maberth_packet (mps_context * ctx, mps_polynomial * p, mps_boolean just_regenerated);

MPS_END_DECLS

#endif /* endif _MPS_JACOBI_ABERTH_H */
