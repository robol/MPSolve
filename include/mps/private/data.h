/*
 * This file is part of MPSolve 3.2.1
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
 * @brief Implementation of the routines that handle the management of data inside
 * mps_context objects.
 */

#ifndef MPS_DATA_H_
#define MPS_DATA_H_

#include <mps/mps.h>

MPS_BEGIN_DECLS

/* functions in data.c */
void mps_mp_set_prec (mps_context * s, long int prec);
void mps_allocate_data (mps_context * s);
void mps_prepare_data (mps_context * s, long int prec);
void mps_restore_data (mps_context * s);
void mps_free_data (mps_context * s);
long int mps_raise_data (mps_context * s, long int prec);
void mps_raise_data_raw (mps_context * s, long int prec);

/* functions in main.c */
void mps_setup (mps_context * s);
void mps_check_data (mps_context * s, char *which_case);
void mps_compute_sep (mps_context * s);

MPS_END_DECLS

#endif /* endif _MPS_DATA_H */
