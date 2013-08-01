/*
 * This file is part of MPSolve 3.0
 *
 * Copyright (C) 2001-2013, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors: 
 *   Leonardo Robol <robol@mail.dm.unipi.it>
 */

#include <mps/mps.h>

mps_monomial_matrix_poly* 
mps_monomial_matrix_poly_new (mps_context * ctx, int degree, int m, mps_boolean monic)
{
  mps_monomial_matrix_poly * poly = mps_new (mps_monomial_matrix_poly);

  MPS_POLYNOMIAL (poly)->degree = degree;
  MPS_POLYNOMIAL (poly)->type_name = "mps_monomial_matrix_poly";

  /* Matrix polynomial specific fields */
  poly->monic = monic;
  poly->m = m;

  /* Allocation of the necessary memory to hold all the matrices. */
  poly->P = malloc (sizeof (double) * poly->m * poly->m * 
		    (degree + 1));

  /* Setup the overloaded methods for our matrix polynomial */
  MPS_POLYNOMIAL (poly)->free = mps_monomial_matrix_poly_free;
  
  mps_polynomial_init (ctx, MPS_POLYNOMIAL (poly));
  
  return poly;
}

void 
mps_monomial_matrix_poly_free (mps_context * ctx, mps_polynomial * poly)
{
  mps_monomial_matrix_poly * mpoly = MPS_MONOMIAL_MATRIX_POLY (poly);

  free (mpoly->P);
  free (poly);
}
