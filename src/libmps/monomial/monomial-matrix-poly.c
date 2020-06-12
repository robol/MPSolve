/*
 * This file is part of MPSolve 3.2.1
 *
 * Copyright (C) 2001-2020, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Leonardo Robol <leonardo.robol@unipi.it>
 */

#include <mps/mps.h>
#include <string.h>

mps_monomial_matrix_poly*
mps_monomial_matrix_poly_new (mps_context * ctx, int degree, int m, mps_boolean monic)
{
  mps_monomial_matrix_poly * poly = mps_new (mps_monomial_matrix_poly);

  MPS_POLYNOMIAL (poly)->degree = degree * m;
  mps_polynomial_init (ctx, MPS_POLYNOMIAL (poly));

  /* Matrix polynomial specific fields */
  poly->monic = monic;
  poly->m = m;
  poly->degree = degree;

  /* Allocation of the necessary memory to hold all the matrices. */
  poly->P = mps_newv (cplx_t, poly->m * poly->m * (degree + 1));
  poly->mP = mps_newv (mpc_t, poly->m * poly->m * (degree + 1));
  poly->mpqPr = mps_newv (mpq_t, poly->m * poly->m * (degree + 1));
  poly->mpqPi = mps_newv (mpq_t, poly->m * poly->m * (degree + 1));

  mpc_vinit2 (poly->mP, poly->m * poly->m * (degree + 1), ctx->mpwp);
  mpq_vinit (poly->mpqPr, poly->m * poly->m * (degree + 1));
  mpq_vinit (poly->mpqPi, poly->m * poly->m * (degree + 1));

  MPS_POLYNOMIAL (poly)->type_name = "mps_monomial_matrix_poly";
  MPS_POLYNOMIAL (poly)->thread_safe = false;

  MPS_POLYNOMIAL (poly)->structure = MPS_STRUCTURE_UNKNOWN;

  /* Setup the overloaded methods for our matrix polynomial */
  MPS_POLYNOMIAL (poly)->free = mps_monomial_matrix_poly_free;
  MPS_POLYNOMIAL (poly)->raise_data = mps_monomial_matrix_poly_raise_data;
  MPS_POLYNOMIAL (poly)->meval = mps_monomial_matrix_poly_meval;

  return poly;
}

void
mps_monomial_matrix_poly_free (mps_context * ctx, mps_polynomial * poly)
{
  mps_monomial_matrix_poly * mpoly = MPS_MONOMIAL_MATRIX_POLY (poly);

  free (mpoly->P);

  mpc_vclear (mpoly->mP, mpoly->m * (poly->degree + mpoly->m));
  free (mpoly->mP);

  mpq_vclear (mpoly->mpqPr, mpoly->m * (poly->degree + mpoly->m));
  free (mpoly->mpqPr);

  mpq_vclear (mpoly->mpqPi, mpoly->m * (poly->degree + mpoly->m));
  free (mpoly->mpqPi);

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


void
mps_monomial_matrix_poly_set_coefficient_d (mps_context * ctx,
                                            mps_monomial_matrix_poly *mpoly,
                                            int i,
                                            cplx_t * matrix)
{
  mps_polynomial *poly = MPS_POLYNOMIAL (mpoly);
  int degree = poly->degree;
  int j;

  if (i < 0 || i > degree)
    {
      mps_error (ctx, "Degree of the coefficient is out of bounds");
      return;
    }

  if (poly->structure == MPS_STRUCTURE_UNKNOWN)
    {
      poly->structure = MPS_STRUCTURE_REAL_FP;
    }
  else if (!MPS_STRUCTURE_IS_FP (poly->structure))
    {
      mps_error (ctx, "Cannot assign floating point coefficients to a non-floating-point polynomial.");
      return;
    }

  /* Copy the data in the coefficients. Please note that we are assuming
   * row-major order in here. */
  cplx_t *ptr = mpoly->P + (mpoly->m * mpoly->m) * i;
  memmove (ptr, matrix, sizeof(cplx_t) * (mpoly->m * mpoly->m));

  for (j = 0; j < mpoly->m * mpoly->m; j++)
    {
      /* Adjust structure if needed */
      if (cplx_Im (matrix[j]) != 0.0)
        poly->structure = MPS_STRUCTURE_COMPLEX_FP;

      mpc_set_cplx (mpoly->mP[j], mpoly->P[j]);
    }
}

void
mps_monomial_matrix_poly_set_coefficient_q (mps_context * ctx,
                                            mps_monomial_matrix_poly *mpoly,
                                            int i,
                                            mpq_t * matrix_r,
                                            mpq_t * matrix_i)
{
  mps_polynomial *poly = MPS_POLYNOMIAL (mpoly);
  int degree = poly->degree;
  int j;

  if (i < 0 || i > degree)
    {
      mps_error (ctx, "Degree of the coefficient is out of bounds");
      return;
    }

  if (poly->structure == MPS_STRUCTURE_UNKNOWN)
    poly->structure = MPS_STRUCTURE_REAL_RATIONAL;

  if (MPS_STRUCTURE_IS_FP (poly->structure))
    {
      mps_error (ctx, "Cannot assign exact coefficients to a floating point polynomial.");
      return;
    }

  for (j = 0; j < mpoly->m * mpoly->m; j++)
    {
      mpq_set (mpoly->mpqPr[i], matrix_r[i]);
      mpq_set (mpoly->mpqPi[i], matrix_i[i]);

      if (mpq_cmp_ui (matrix_i[i], 0U, 1U) != 0)
        poly->structure = MPS_STRUCTURE_COMPLEX_RATIONAL;
    }
}

mps_boolean
mps_monomial_matrix_poly_meval (mps_context * ctx, mps_polynomial * poly,
                                mpc_t x, mpc_t value, rdpe_t error)
{
  mps_monomial_matrix_poly *mpoly = MPS_MONOMIAL_MATRIX_POLY (poly);

  if (!mpoly->flags & MPS_MONOMIAL_MATRIX_POLY_HESSENBERG)
    {
      /* TODO: Insert the necessary steps to make a Hessenberg matrix
       * available in place of this non-Hessenberg matrix polynomial. */
    }

  /* TODO: This is not thread safe - at all! Insert a double buffer strategy
   * in here, as for the case of mps_monomial_poly. */
  if (mpc_get_prec (mpoly->mP[0]) < mpc_get_prec (value))
    mps_monomial_matrix_poly_raise_data (ctx, poly, mpc_get_prec (value));

  /* Otherwise just proceed with Horner and our recursive implementation. */
  /* TODO: We are performing a floating point evaluation here, just for simplicity. */
  mps_mhessenberg_shifted_determinant (ctx, mpoly->mP, x, mpoly->m, value, error);

  return true;
}

long int
mps_monomial_matrix_poly_raise_data (mps_context * ctx,
                                     mps_polynomial * p,
                                     long int wp)
{
  int i;
  mps_monomial_matrix_poly *mpoly = MPS_MONOMIAL_MATRIX_POLY (p);
  size_t n_elems = (mpoly->degree + 1) * (mpoly->m * mpoly->m);

  for (i = 0; i < n_elems; i++)
    {
      mpc_set_prec (mpoly->mP[i], wp);
    }

  if (MPS_STRUCTURE_IS_INTEGER (p->structure) ||
      MPS_STRUCTURE_IS_RATIONAL (p->structure))
    {
      mpc_set_q (mpoly->mP[i], mpoly->mpqPr[i], mpoly->mpqPi[i]);
    }

  return mpc_get_prec (mpoly->mP[0]);
}
