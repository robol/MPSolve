#ifndef __MPS_MONOMIAL_POLY_H
#define __MPS_MONOMIAL_POLY_H

/**
 * @file
 * @brief Implementation of the allocation and edit functions for the 
 * handling of monomial polynomials. 
 */


  #include <mps/polynomial.h>
  #include <mps/mps.h>
  #include <gmp.h>
  #include <pthread.h>

#define MPS_MONOMIAL_POLY(t) ((mps_monomial_poly*) t)
#define MPS_IS_MONOMIAL_POLY(t) (mps_polynomial_check_type (t, "mps_monomial_poly"))

#ifdef  __cplusplus
extern "C"
{
#endif

#ifdef _MPS_PRIVATE

  struct mps_monomial_poly_double_buffer {
    char active;
    mpc_t *mfpc1;
    mpc_t *mfpc2;
  };


  /**
   * @brief Data regarding a polynomial represented in the monomial
   * base.
   */
  struct mps_monomial_poly {

    /**
     * @brief Implementation of the methods. 
     */
    struct mps_polynomial methods;

    struct mps_monomial_poly_double_buffer db;
    
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

    /**
     * @brief Precision of the polynomial coefficients.
     */
    long int prec;
        
  };
#endif /* #ifdef _MPS_PRIVATE */

  /* These routines are thought for polynomial handling, i.e. allocating and 
   * setting coefficients of the polynomials, and setting the precision of the
   * floating point coefficients that are in there */

  mps_monomial_poly * mps_monomial_poly_new (mps_context * s, long int degree);

  void mps_monomial_poly_free (mps_context * s, mps_polynomial * mp);

  long int mps_monomial_poly_get_precision (mps_context * s, mps_monomial_poly * mp);
  
  long int mps_monomial_poly_raise_precision (mps_context * s, mps_polynomial * mp, long int prec);

  void mps_monomial_poly_set_coefficient_q (mps_context * s, mps_monomial_poly * mp, long int i, 
                                            mpq_t real_part, mpq_t imag_part);
  void mps_monomial_poly_set_coefficient_d (mps_context * s, mps_monomial_poly * mp, long int i,
                                            double real_part, double imag_part);
  void mps_mononomial_poly_set_coefficient_f (mps_context * s, mps_monomial_poly * p, long int i,
                                              mpc_t coeff);
  void mps_monomial_poly_set_coefficient_int (mps_context * s, mps_monomial_poly * mp, long int i,
                                              long long real_part, long long imag_part);
  mps_monomial_poly * mps_monomial_poly_derive (mps_context * s, mps_monomial_poly * p, int k, long int wp);

  mps_boolean mps_monomial_poly_feval (mps_context * ctx, mps_polynomial *p, cplx_t x, cplx_t value, double * error);

  mps_boolean mps_monomial_poly_deval (mps_context * ctx, mps_polynomial *p, cdpe_t x, cdpe_t value, rdpe_t error);

  mps_boolean mps_monomial_poly_meval (mps_context * ctx, mps_polynomial *p, mpc_t x, mpc_t value, rdpe_t error);

  void mps_monomial_poly_fstart (mps_context * ctx, mps_polynomial * p);

  void mps_monomial_poly_dstart (mps_context * ctx, mps_polynomial * p);

  void mps_monomial_poly_mstart (mps_context * ctx, mps_polynomial * p);

  void mps_monomial_poly_fnewton (mps_context * ctx, mps_polynomial * p, 
                                  mps_approximation * root, cplx_t corr);

  void mps_monomial_poly_dnewton (mps_context * ctx, mps_polynomial * p, 
                                  mps_approximation * root, cdpe_t corr);

  void mps_monomial_poly_mnewton (mps_context * ctx, mps_polynomial * p, 
                                  mps_approximation * root, mpc_t corr);

  void mps_monomial_poly_get_leading_coefficient (mps_context * ctx, mps_polynomial * p,
                                                  mpc_t leading_coefficient);

  void mps_monomial_poly_deflate (mps_context * ctx, mps_polynomial * p);


#ifdef  __cplusplus
}
#endif

#endif
