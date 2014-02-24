/*
 * This file is part of MPSolve 3.1.5
 *
 * Copyright (C) 2001-2014, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Leonardo Robol <robol@mail.dm.unipi.it>
 */

#ifndef MPS_MONOMIAL_MATRIX_POLY_H_
#define MPS_MONOMIAL_MATRIX_POLY_H_

/**
 * @file
 * @brief Implementation of the monomial version of the matrix
 * polynomial.
 */

#include <mps/polynomial.h>

#define MPS_MONOMIAL_MATRIX_POLY(t) (MPS_POLYNOMIAL_CAST (mps_monomial_matrix_poly, t))
#define MPS_IS_MONOMIAL_MATRIX_POLY(t) \
  (mps_polynomial_check_type (t, "mps_monomial_matrix_poly"))

#ifdef  __cplusplus
extern "C"
{
#endif

#define MPS_MONOMIAL_MATRIX_POLY_HESSENBERG 0x0001

#ifdef _MPS_PRIVATE

/**
 * @brief This is the struct that holds all the data of the matrix
 * polynomial.
 */
struct mps_monomial_matrix_poly {
  /**
   * @brief Implementation of the overloaded methods for the
   * matrix polynomial.
   */
  mps_polynomial methods;

  /**
   * @brief If this flag is set to true then the
   * higher degree term of the polynomial is the identity
   * matrix and so doesn't need to be allocated and/or
   * accessed in any way.
   *
   * In particular, it's not guaranteed to be available for
   * the computations, so you should always check if this flag
   * is set before trying to operate on P[n].
   */
  mps_boolean monic;

  /**
   * @brief The size of the matrices that compose the matrix
   * polynomial.
   */
  int m;

  /**
   * @brief Degree of the matrix polynomial. Please note that this values if different
   * from the degree value in the polynomial structure, since that is the degree of the
   * scalar polynomial, and this is the degree of the matrix one.
   */
  int degree;

  /**
   * @brief The double version of the polynomial coefficients.
   *
   * NOTE: At this stage, this is the only type of data that is kept
   * for the matrix polynomial.
   */
  cplx_t * P;

  /**
   * @brief The multiprecision version of the polynomial coefficients.
   */
  mpc_t * mP;

  /**
   * @brief The rational version of the polynomial coefficients.
   * This is used only if the structure of the monomial matrix poly
   * is MPS_STRUCTURE_*_INTEGER OR MPS_STRUCTURE_*_RATIONAL.
   *
   * This array holds the real part of the coefficients.
   */
  mpq_t * mpqPr;

  /**
   * @brief The rational version of the polynomial coefficients.
   * This is used only if the structure of the monomial matrix poly
   * is MPS_STRUCTURE_*_INTEGER OR MPS_STRUCTURE_*_RATIONAL.
   *
   * This array holds the real part of the coefficients.
   */
  mpq_t * mpqPi;

  /**
   * @brief Additional properties of this polynomial. Examples are:
   * MPS_MONOMIAL_MATRIX_POLY_HESSENBERG, ...
   *
   * @seealso mps_monomial_matrix_poly_add_flag(),
   * @seealso mps_monomial_matrix_poly_clear_flags().
   */
  int flags;
};


#endif /* _MPS_PRIVATE */

/* From here on you will find declaration of the public method
   available for the matrix polynomials */

/**
 * @brief Create a new matrix polynomial of the given degree.
 *
 * @param ctx The current mps_context.
 * @param degree The degree of the matrix polynomial.
 * @param m The size of the matrices that compose the matrix polynomial.
 * @param monic A boolean value that, if set to true, specify that the leading
 *              coefficient of the polynomial is the identity matrix, and so
 *              should not specified explicitely.
 * @return A pointer to a newly allocated mps_monomial_matrix_poly that should
 * be subsequently free with a call to mps_monomial_matrix_poly_free().
 */
mps_monomial_matrix_poly* mps_monomial_matrix_poly_new (mps_context * ctx,
                                                        int degree,
                                                        int m,
                                                        mps_boolean monic);
/**
 * @brief Free a matrix polynomial.
 *
 * @param ctx The current mps_context.
 * @param poly The mps_monomial_matrix_poly that should be freed, casted to
 *             a mps_polynomial* pointer.
 */
void mps_monomial_matrix_poly_free (mps_context * ctx,
                                    mps_polynomial * poly);

/**
 * @brief Add some flags (some properties) to this matrix polynomial.
 *
 * @param ctx The current mps_context
 * @param mpoly The matrix polynomial.
 * @param flag The flags to add.
 */
void mps_monomial_matrix_poly_add_flags (mps_context * ctx,
                                         mps_monomial_matrix_poly * mpoly,
                                         int flag);

/**
 * @brief Clear some flags (properties) of this matrix polynomial.
 *
 * @param ctx The current mps_context
 * @param mpoly The matrix polynomial.
 * @param flag The flags to clear.
 */
void mps_monomial_matrix_poly_clear_flags (mps_context * ctx,
                                           mps_monomial_matrix_poly * mpoly,
                                           int flag);

/**
 * @brief Set the coefficient of degree i of the matrix polynomial.
 *
 * @param ctx The current mps_context
 * @param mpoly The mps_monomial_matrix_poly where the coefficients should be
 * set.
 * @param i The degree of the coeffient to set.
 * @param matrix A pointer to the first element of the matrix stored in row-major order.
 */
void mps_monomial_matrix_poly_set_coefficient_d (mps_context * ctx,
                                                 mps_monomial_matrix_poly *mpoly,
                                                 int i,
                                                 cplx_t * matrix);

/**
 * @brief Set the coefficient of degree i of the matrix polynomial.
 *
 * @param ctx The current mps_context
 * @param mpoly The mps_monomial_matrix_poly where the coefficients should be
 * set.
 * @param i The degree of the coeffient to set.
 * @param matrix_r A pointer to the first element of the matrix of the real parts
 * of the coefficients, stored in row-major order
 * @param matrix_i A pointer to the first element of the matrix of the imaginary parts
 * of the coefficients, stored in row-major order
 */
void mps_monomial_matrix_poly_set_coefficient_q (mps_context * ctx,
                                                 mps_monomial_matrix_poly *mpoly,
                                                 int i,
                                                 mpq_t * matrix_r,
                                                 mpq_t * matrix_i);

/**
 * @brief Evaluate a matrix polynomial at a point, in the sense of
 * evaluating \f$det(P(x))\f$.
 *
 * @param ctx The current mps_context
 * @param poly The matrix polynomial to evaluate
 * @param x The point in which the evaluation is requested
 * @param value The value of \f$det(P(x))\f$
 * @param error An upper bound to the absolute error that affects the result.
 * @return true if the evaluation was successful.
 */
mps_boolean mps_monomial_matrix_poly_meval (mps_context * ctx,
                                            mps_polynomial * poly,
                                            mpc_t x,
                                            mpc_t value,
                                            rdpe_t error);

/**
 * @brief Raise the working precision of this monomial matrix polynomal
 * to the required numnber of bits.
 *
 * @param ctx The current mps_context.
 * @param p The polynomial whose working precision should be set.
 * @param wp The bits of desired working precision.
 * @return The precision set in the polynomial. Note that this value
 * may be higher than wp due to the fact that, in general, not all
 * the precisions are available on the system and the first higher precision
 * of the value required will be used.
 */
long int
mps_monomial_matrix_poly_raise_data (mps_context * ctx,
                                     mps_polynomial * p,
                                     long int wp);

#ifdef  __cplusplus
}
#endif


#endif /* __MPS_MONOMIAL_MATRIX_POLY_H */
