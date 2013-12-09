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
#include <string.h>

mps_monomial_matrix_poly* 
mps_monomial_matrix_poly_new (mps_context * ctx, int degree, int m, mps_boolean monic)
{
  mps_monomial_matrix_poly * poly = mps_new (mps_monomial_matrix_poly);

  MPS_POLYNOMIAL (poly)->degree = degree;

  /* Matrix polynomial specific fields */
  poly->monic = monic;
  poly->m = m;

  /* Allocation of the necessary memory to hold all the matrices. */
  poly->P = mps_malloc (sizeof (cplx_t) * poly->m * poly->m * 
			(degree + 1));

  /* Setup the overloaded methods for our matrix polynomial */
  MPS_POLYNOMIAL (poly)->free = mps_monomial_matrix_poly_free;
  
  mps_polynomial_init (ctx, MPS_POLYNOMIAL (poly));
  MPS_POLYNOMIAL (poly)->type_name = "mps_monomial_matrix_poly";
  
  return poly;
}

void 
mps_monomial_matrix_poly_free (mps_context * ctx, mps_polynomial * poly)
{
  mps_monomial_matrix_poly * mpoly = MPS_MONOMIAL_MATRIX_POLY (poly);

  free (mpoly->P);
  free (poly);
}

void mps_monomial_matrix_poly_add_flags (mps_context * ctx, 
					 mps_monomial_matrix_poly * mpoly, 
					 int flags)
{
  mpoly->flags |= flags; 
}


void mps_monomial_matrix_poly_clear_flags (mps_context * ctx, 
					   mps_monomial_matrix_poly * mpoly, 
					   int flags)
{
  mpoly->flags &= ~flags;
}


void mps_monomial_matrix_poly_set_coefficient_d (mps_context * ctx, 
						 mps_monomial_matrix_poly *mpoly, 
						 int i, 
						 cplx_t * matrix)
{
  mps_polynomial *poly = MPS_POLYNOMIAL (mpoly); 
  int degree = poly->degree; 

  if (i < 0 || i > degree) {
    mps_error (ctx, "Degree of the coefficient is out of bounds"); 
    return; 
  }

  /* Copy the data in the coefficients. Please note that we are assuming
   * row-major order in here. */
  cplx_t *ptr = mpoly->P + (mpoly->m * mpoly->m) * i;
  memcpy (ptr, matrix, sizeof (cplx_t) * (mpoly->m * mpoly->m)); 
}

mps_boolean
mps_monomial_matrix_poly_meval (mps_context * ctx, mps_polynomial * poly,
				mpc_t x, mpc_t value, rdpe_t error)
{
  mps_monomial_matrix_poly *mpoly = MPS_MONOMIAL_MATRIX_POLY (poly); 

  if (! mpoly->flags & MPS_MONOMIAL_MATRIX_POLY_HESSENBERG) {
    /* TODO: Insert the necessary steps to make a Hessenberg matrix 
     * available in place of this non-Hessenberg matrix polynomial. */
  }

  /* Otherwise just proceed with Horner and our recursive implementation. */  


  return true;
}
