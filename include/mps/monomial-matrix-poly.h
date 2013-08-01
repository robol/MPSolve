#ifndef __MPS_MONOMIAL_MATRIX_POLY_H
#define __MPS_MONOMIAL_MATRIX_POLY_H

/**
 * @file
 * @brief Implementation of the monomial version of the matrix
 * polynomial. 
 */

#include <mps/polynomial.h>

#define MPS_MONOMIAL_MATRIX_POLY(t) ((mps_monomial_matrix_poly*) t)
#define MPS_IS_MONOMIAL_MATRIX_POLY(t) \
  (mps_polynomial_check_type (t, "mps_monomial_matrix_poly"))

#ifdef  __cplusplus
extern "C"
{
#endif

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
     * @brief The double version of the polynomial coefficients.
     *
     * NOTE: At this stage, this is the only type of data that is kept
     * for the matrix polynomial. 
     */
    double * P;

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

#ifdef  __cplusplus
}
#endif


#endif /* __MPS_MONOMIAL_MATRIX_POLY_H */
