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
#include <string.h>

/**
 * @brief This is the full implementation of the recursive determinant computation.
 *
 * @param ctx The current mps_context
 * @param hessenberg_matrix The hessenberg matrix whose determinant should be computed.
 * @param n The size of the matrix.
 * @param output The storage for the result.
 */
void
mps_fhessenberg_determinant (mps_context * ctx, cplx_t * hessenberg_matrix, size_t n, cplx_t output, 
			     long int * exponent)
{
  mps_fhessenberg_shifted_determinant (ctx, hessenberg_matrix, cplx_zero, n, output, exponent);
}

/**
 * @brief This is the full implementation of the recursive determinant computation of the
 * Hessnberg - \lambda I matrix.
 *
 * @param ctx The current mps_context
 * @param hessenberg_matrix The hessenberg matrix whose determinant should be computed.
 * @param shift The value of \f$\lambda\f$.
 * @param n The size of the matrix.
 * @param output The storage for the result.
 */
void
mps_fhessenberg_shifted_determinant (mps_context * ctx, cplx_t * hessenberg_matrix, const cplx_t shift, 
				     size_t n, cplx_t output, long int *acc_exponent)
{
  cplx_t *vec = mps_newv (cplx_t, n);
  size_t local_n = n;
  int i, j;
  *acc_exponent = 0;

  /* Copy the last column in vec */
  for (i = 0; i < n; i++)
    {
      cplx_set (vec[i], MPS_MATRIX_ELEM (hessenberg_matrix, i, n - 1, n));
    }
  cplx_sub_eq (vec[n-1], shift);

  while (local_n-- > 1)
    {
      /* Compress the last two cols of the matrix */
      cplx_t t, s;

      for (i = 0; i < local_n - 1; i++)
        {
          cplx_mul (s, MPS_MATRIX_ELEM (hessenberg_matrix, i, local_n - 1, n), 
		    vec[local_n]);
          cplx_mul (t, vec[i], MPS_MATRIX_ELEM (hessenberg_matrix, local_n, local_n - 1, n));
          cplx_sub (vec[i], s, t);
        }

      /* The last step require the extra accounting for the shifted case. */
      cplx_sub (s, MPS_MATRIX_ELEM (hessenberg_matrix, local_n - 1, local_n - 1, n), shift);
      cplx_mul_eq (s, vec[local_n]);
      cplx_mul (t, vec[local_n-1], MPS_MATRIX_ELEM (hessenberg_matrix, local_n, local_n - 1, n));
      cplx_sub (vec[local_n - 1], s, t);

      if (i % 50 == 0)
	{
	  int exponent;
	  double max = 0.0;

	  for (j = 0; j < local_n; j++)
	    {
	      double m = cplx_mod (vec[j]);
	      if (m > max)
		max = m;
	    }

	  frexp (max, &exponent);
	  max = pow(2.0, exponent);
	  
	  for (j = 0; j < local_n; j++)
	    {
	      cplx_div_eq_d (vec[j], max);
	    }

	  *acc_exponent += exponent;
	}
    }

  cplx_set (output, *vec);
  free (vec);
}

/**
 * @brief This is the full implementation of the recursive determinant computation.
 *
 * @param ctx The current mps_context
 * @param hessenberg_matrix The hessenberg matrix whose determinant should be computed.
 * @param n The size of the matrix.
 * @param output The storage for the result.
 */
void
mps_dhessenberg_determinant (mps_context * ctx, cdpe_t * hessenberg_matrix, size_t n, cdpe_t output)
{
  mps_dhessenberg_shifted_determinant (ctx, hessenberg_matrix, cdpe_zero, n, output);
}

/**
 * @brief This is the full implementation of the recursive determinant computation of the
 * Hessenberg - \lambda I matrix.
 *
 * @param ctx The current mps_context
 * @param hessenberg_matrix The hessenberg matrix whose determinant should be computed.
 * @param shift The value of \f$\lambda\f$.
 * @param n The size of the matrix.
 * @param output The storage for the result.
 */
MPS_PRIVATE void
mps_dhessenberg_shifted_determinant (mps_context * ctx, cdpe_t * hessenberg_matrix, const cdpe_t shift, 
				     size_t n, cdpe_t output)
{
  cdpe_t *vec = mps_newv (cdpe_t, n);
  size_t local_n = n;
  int i;

  /* Copy the last column in vec */
  for (i = 0; i < n; i++)
    {
      cdpe_set (vec[i], MPS_MATRIX_ELEM (hessenberg_matrix, i, n - 1, n));
    }
  cdpe_sub_eq (vec[n], shift);

  while (local_n-- > 1)
    {
      /* Compress the last two cols of the matrix */
      cdpe_t t, s;

      for (i = 0; i < local_n - 1; i++)
        {
          cdpe_mul (s, MPS_MATRIX_ELEM (hessenberg_matrix, i, local_n - 1, n), 
		    vec[local_n]);
          cdpe_mul (t, vec[i], MPS_MATRIX_ELEM (hessenberg_matrix, local_n, local_n - 1, n));
          cdpe_sub (vec[i], s, t);
        }

      /* The last step require the extra accounting for the shifted case. */
      cdpe_sub (s, MPS_MATRIX_ELEM (hessenberg_matrix, local_n - 1, local_n - 1, n), shift);
      cdpe_mul_eq (s, vec[local_n]);
      cdpe_mul (t, vec[local_n-1], MPS_MATRIX_ELEM (hessenberg_matrix, local_n, local_n - 1, n));
      cdpe_sub (vec[local_n - 1], s, t);
    }

  cdpe_set (output, *vec);
  free (vec);
}


/**
 * @brief This is the full implementation of the recursive determinant computation.
 *
 * @param ctx The current mps_context
 * @param hessenberg_matrix The hessenberg matrix whose determinant should be computed.
 * @param n The size of the matrix.
 * @param output The storage for the result.
 * @param error A bound on the error.
 */
MPS_PRIVATE void
mps_mhessenberg_determinant (mps_context * ctx, mpc_t * hessenberg_matrix, size_t n, mpc_t output, rdpe_t error)
{
  mpc_t zero;

  mpc_init2 (zero, mpc_get_prec (output));
  mpc_set_ui (zero, 0U, 0U);

  mps_mhessenberg_shifted_determinant (ctx, hessenberg_matrix, zero, n, output, error);

  mpc_clear (zero);
}


/**
 * @brief This is the full implementation of the recursive determinant computation of the
 * Hessnberg - \lambda I matrix.
 *
 * @param ctx The current mps_context
 * @param hessenberg_matrix The hessenberg matrix whose determinant should be computed.
 * @param shift The value of \f$\lambda\f$.
 * @param n The size of the matrix.
 * @param output The storage for the result.
 * @param error A bound on the error.
 */
MPS_PRIVATE void
mps_mhessenberg_shifted_determinant (mps_context * ctx, mpc_t * hessenberg_matrix, mpc_t shift, size_t n, mpc_t output, rdpe_t error)
{
  mpc_t *matrix = mps_newv (mpc_t, n * n);
  size_t local_n = n;
  long int wp = mpc_get_prec (output);
  mps_boolean shifted = !mpc_eq_zero (shift);
  int i, j;
  mpc_t t, s;
  rdpe_t mod, epsilon;

  rdpe_t *verrors = mps_newv (rdpe_t, n);

  memset (verrors, 0, sizeof(rdpe_t) * n);

  /* Init the elements in the matrix */
  mpc_init2 (t, wp);
  mpc_init2 (s, wp);
  mpc_vinit2 (matrix, n * n, wp);

  rdpe_set_2dl (epsilon, 1.0, 1 - wp);
  rdpe_set (error, rdpe_one);

  /* Make a copy of the original matrix and subtract the shift, if any*/
  for (i = 0; i < n; i++)
    {
      for (j = 0; j < n; j++)
        {
          if ((i == j) && shifted)
            mpc_sub (MPS_MATRIX_ELEM (matrix, i, j, n),
                     MPS_MATRIX_ELEM (hessenberg_matrix, i, j, n),
                     shift);
          else
            mpc_set (MPS_MATRIX_ELEM (matrix, i, j, n),
                     MPS_MATRIX_ELEM (hessenberg_matrix, i, j, n));
        }
    }

  while (local_n-- > 1)
    {
      /* Recursively compress the last two cols of the matrix */
      rdpe_t err_a_bottom, err_b_bottom, tmp;

      mpc_rmod (err_a_bottom, MPS_MATRIX_ELEM (matrix, local_n, local_n - 1, n));
      mpc_rmod (err_b_bottom, MPS_MATRIX_ELEM (matrix, local_n, local_n, n));

      rdpe_mul (tmp, err_b_bottom, epsilon);
      rdpe_add_eq (verrors[local_n], tmp);

      for (i = 0; i < local_n; i++)
        {
          rdpe_t err_a, err_b;

          mpc_rmod (err_a, MPS_MATRIX_ELEM (matrix, i, local_n - 1, n));
          mpc_rmod (err_b, MPS_MATRIX_ELEM (matrix, i, local_n, n));

          rdpe_mul_eq (err_a, verrors[local_n]);

          rdpe_mul_eq (err_b, epsilon);
          rdpe_add_eq (err_b, verrors[i]);
          rdpe_mul_eq (err_b, err_a_bottom);

          mpc_mul (s, MPS_MATRIX_ELEM (matrix, i, local_n - 1, n),
                   MPS_MATRIX_ELEM (matrix, local_n, local_n, n));
          mpc_mul (t, MPS_MATRIX_ELEM (matrix, i, local_n, n),
                   MPS_MATRIX_ELEM (matrix, local_n, local_n - 1, n));
          mpc_sub (MPS_MATRIX_ELEM (matrix, i, local_n - 1, n), s, t);

          mpc_rmod (mod, MPS_MATRIX_ELEM (matrix, i, local_n - 1, n));
          rdpe_mul_eq (mod, epsilon);

          rdpe_add_eq (verrors[i], mod);
          rdpe_add_eq (verrors[i], err_a);
          rdpe_add_eq (verrors[i], err_b);
        }
    }

  rdpe_set (error, verrors[0]);

  mpc_set (output, *matrix);

  mpc_vclear (matrix, n * n);
  free (matrix);

  mpc_clear (t);
  mpc_clear (s);
}
