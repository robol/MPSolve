/*
 * This file is part of MPSolve 3.1.5
 *
 * Copyright (C) 2001-2014, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Leonardo Robol <leonardo.robol@sns.it>
 */

#include <mps/mps.h>

/**
 * @brief Select appropriate starting point for the approximation of the roots
 * of the given polynomial by applying a divide-and-conquer strategy described
 * in {TODO: Reference missing}. 
 *
 * @param ctx The current mps_context. 
 * @param poly The polynomial whose roots should be approximated.
 */
void
mps_recursive_fstart (mps_context * ctx, mps_polynomial * poly)
{
  // Note: this function may take a long time to finish, since it will 
  // recursively solve polynomials of lower degree. 
}

/**
 * @brief Select appropriate starting point for the approximation of the roots
 * of the given polynomial by applying a divide-and-conquer strategy described
 * in {TODO: Reference missing}. 
 *
 * @param ctx The current mps_context. 
 * @param poly The polynomial whose roots should be approximated.
 */
void
mps_recursive_dstart (mps_context * ctx, mps_polynomial * poly)
{
  
}

/**
 * @brief Select appropriate starting point for the approximation of the roots
 * of the given polynomial by applying a divide-and-conquer strategy described
 * in {TODO: Reference missing}. 
 *
 * @param ctx The current mps_context. 
 * @param poly The polynomial whose roots should be approximated.
 */
void
mps_recursive_mstart (mps_context * ctx, mps_polynomial * poly)
{
  
}
