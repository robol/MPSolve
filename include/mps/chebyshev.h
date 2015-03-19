/*
 * This file is part of MPSolve 3.1.5
 *
 * Copyright (C) 2001-2015, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Leonardo Robol <leonardo.robol@sns.it>
 */

#ifndef MPS_CHEBYSHEV_H_
#define MPS_CHEBYSHEV_H_

MPS_BEGIN_DECLS

 #define MPS_CHEBYSHEV_POLY_TYPE_NAME "mps_chebyshev_poly"
 #define MPS_CHEBYSHEV_POLY(t) ((mps_chebyshev_poly*)t)
 #define MPS_IS_CHEBYSHEV_POLY(t) mps_polynomial_check_type (t, "mps_chebyshev_poly")

typedef struct {
  /**
   * @brief Base implementation of a polynomial.
   */
  mps_polynomial super;

  /**
   * @brief Floating point coefficients of the polynomial in the Chebyshev
   * base
   */
  cplx_t * fpc;

  /**
   * @brief DPE floating point coefficients of the polynomial in the Chebyshev
   * base.
   */
  cdpe_t * dpc;

  /**
   * @brief Multiprecision complex coefficients of the polynomial in the Chebyshev
   * base.
   */
  mpc_t * mfpc;

  /**
   * @brief Rational coefficients of the polynomial. These are the real parts of the
   * coefficients.
   */
  mpq_t * rational_real_coeffs;

  /**
   * @brief Ratinonal coefficients of the polynomial. These are the imaginary parts
   * of the coefficients.
   */
  mpq_t * rational_imag_coeffs;

  /**
   * @brief Leading coefficient of the polynomial.
   */
  mpc_t lc;

  /**
   * @brief Internal mutex used to manage the change of precision.
   */
  pthread_mutex_t precision_mutex;
} mps_chebyshev_poly;


/**
 * @brief Create a new polynomial represented in the Chebyshev base with
 * degree set to n.
 */
mps_chebyshev_poly * mps_chebyshev_poly_new (mps_context * ctx, int n, mps_structure structure);

/**
 * @brief Set the coefficient relative to the i-th element of the Chebyshev
 * base.
 *
 * This function takes a rational number as input, and is usable only if the
 * Chebyshev polynomial is represented using rational coefficients.
 */
void mps_chebyshev_poly_set_coefficient_q (mps_context * ctx, mps_chebyshev_poly * poly, int i,
                                           mpq_t real_part, mpq_t imag_part);

/**
 * @brief Set the coefficient relative to the i-th element of the Chebyshev
 * base.
 *
 * This function takes a multiprecision floating point number as input.
 */
void mps_chebyshev_poly_set_coefficient_f (mps_context * ctx, mps_chebyshev_poly * poly,
                                           int i, mpc_t coeff);

/**
 * @brief Set the coefficient of the i-th element of the Chebyshev base.
 *
 * This function takes an integer value as input.
 *
 * @param ctx The current mps_context
 * @param poly The Chebyshev polynomial whose coefficient should be set.
 * @param i The degree of the coefficient to set.
 * @param real_coeff The real part of the new coefficient.
 * @param imag_coeff The imaginary part of the new coefficient.
 */
void mps_chebyshev_poly_set_coefficient_i (mps_context * ctx, mps_chebyshev_poly * poly,
                                           int i, long int real_coeff, long int imag_coeff);


mps_chebyshev_poly * mps_chebyshev_poly_read_from_stream (mps_context * ctx, mps_input_buffer * buffer,
                                                          mps_structure structure, mps_density density,
                                                          long int precision);


MPS_END_DECLS

#endif
