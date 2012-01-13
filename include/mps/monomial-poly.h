#ifndef __MPS_MONOMIAL_POLY_H
#define __MPS_MONOMIAL_POLY_H

/**
 * @file
 * @brief Implementation of the allocation and edit functions for the 
 * handling of monomial polynomials. 
 */

#ifdef	__cplusplus
extern "C"
{
#endif

  #include <mps/mps.h>
  #include <gmp.h>
  #include <pthread.h>

  /**
   * @brief Data regarding a polynomial represented in the monomial
   * base.
   */
  struct mps_monomial_poly {
    
    /**
     * @brief The degree of the polynomial.
     */
    int n;

    /**
     * @brief Structure of the polynomial. This structure could be implicitly
     * set the first time that a coefficient is set with the appropriate routine,
     * and cannot be modified after that.
     */
    mps_structure structure;

    /**
     * @brief This array contains the structure of the sparse
     * polynomial.
     *
     * <code>spar[i]</code> is <code>true</code> if and only if
     * the i-th coefficients of the polynomial is a non-zero
     * coefficients
     */
    mps_boolean *spar;

    /**
     * @brief Standard real coefficients.
     */
    double *fpr;

    /**
     * @brief Standard complex coefficients.
     */
    cplx_t *fpc;

    /**
     * @brief Array containing standard complex coefficients
     */
    cplx_t *fppc;

    /**
     * @brief Dpe real coefficients.
     */
    rdpe_t *dpr;

    /**
     * @brief Dpe complex coefficients.
     */
    cdpe_t *dpc;

    /**
     * @brief Multiprecision real coefficients.
     */
    mpf_t *mfpr;

    /**
     * @brief Multiprecision complex coefficients.
     */
    mpc_t *mfpc;
    
    /**
     * @brief Array of mutexes that need to be locked when reading at the
     * i-th compoenent of the poly. 
     */
    pthread_mutex_t * mfpc_mutex;

    /**
     * @brief Multiprecision complex coefficients of \f$p'(x)\f$.
     */
    mpc_t *mfppc;

    /**
     * @brief Array containing moduli of the coefficients as double numbers.
     */
    double *fap;

    /**
     * @brief Array containing moduli of the coefficients as dpe numbers.
     */
    rdpe_t *dap;

    /**
     * @brief Real part of rational input coefficients.
     */
    mpq_t *initial_mqp_r;

    /**
     * @brief Imaginary part of rational input coefficients.
     */
    mpq_t *initial_mqp_i;

    /**
     * @brief This mutex must be locked while regenerating the coefficients
     * of the polynomial.
     */
    pthread_mutex_t regenerating;
        
  };

  /* These routines are thought for polynomial handling, i.e. allocating and 
   * setting coefficients of the polynomials, and setting the precision of the
   * floating point coefficients that are in there */

  mps_monomial_poly * mps_monomial_poly_new (mps_status * s, long int degree);

  void mps_monomial_poly_free (mps_status * s, mps_monomial_poly * mp);

  void mps_monomial_poly_raise_precision (mps_status * s, mps_monomial_poly * mp, long int prec);

  void mps_monomial_poly_set_coefficient_q (mps_status * s, mps_monomial_poly * mp, long int i, 
					    mpq_t real_part, mpq_t imag_part);
  void mps_monomial_poly_set_coefficient_d (mps_status * s, mps_monomial_poly * mp, long int i,
					    double real_part, double imag_part);

  void mps_monomial_poly_set_coefficient_int (mps_status * s, mps_monomial_poly * mp, long int i,
					      long int real_part, long int imag_part);



#ifdef	__cplusplus
}
#endif

#endif
