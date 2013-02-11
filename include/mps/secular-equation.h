/* 
 * File:   secular.h
 * Author: leonardo
 *
 * Created on 29 aprile 2011, 16.26
 */

/**
 * @file 
 * @brief Header file for secular-related routines.
 */

#ifndef SECULAR_H
#define  SECULAR_H

#include <mps/mps.h>
#include <float.h>

#ifdef  __cplusplus
extern "C"
{
#endif

#define MPS_SECULAR_EQUATION(t) ((mps_secular_equation *) t)
#define MPS_IS_SECULAR_EQUATION(t) (mps_polynomial_check_type (t, "mps_secular_equation"))

#ifdef _MPS_PRIVATE

  /* CONSTANTS */

  /**
   * @brief This is the number of bits used when first passed
   * in multiprecision. 
   */
#define MPS_SECULAR_STARTING_MP_PRECISION 128

  /**
   * @brief This is the higher precision supported by GMP that is 
   * lower than the precision supported by the standard floating
   * point machinery. It is used to set an "equivalent" precision
   * in s->mpwp for the step of multiprecision coefficient regeneration.
   */
#define MPS_SECULAR_EQUIVALENT_FP_PRECISION (MPS_SECULAR_STARTING_MP_PRECISION / 2)

  struct mps_secular_equation_double_buffer
  {
    char active;
    mpc_t *ampc1;
    mpc_t *ampc2;
    mpc_t *bmpc1;
    mpc_t *bmpc2;
  };

  /**
   * @brief Secular equation data.
   *
   * A secular equation is an equation in the form
   * \f[
   *   \sum_{i = 1}^{n} \frac{a_i}{z - b_i} = 1
   * \f]
   * and this struct holds the values of the parameters \f$a_i\f$
   * and \f$b_i\f$.
   */
  struct mps_secular_equation
  {
    struct mps_polynomial __base_class__;

    struct mps_secular_equation_double_buffer db;

    /**
     * @brief Vector of \f$a_i\f$ as complex floating
     * point numbers.
     */
    cplx_t *afpc;

    /**
     * @brief Same as <code>afpc</code>, but the <code>dpe</code>
     * version.
     */
    cdpe_t *adpc;

    /**
     * @brief Vector with the values of \f$b_i\f$ as complex
     * floating point numbers.
     */
    cplx_t *bfpc;

    /**
     * @brief Same as <code>bfpc</code>, but the <code>dpe</code>
     * version.
     */
    cdpe_t *bdpc;

    /**
     * @brief Same as <code>afpc</code>, but the multiprecision
     * version.
     */
    mpc_t *ampc;

    /**
     * @brief Mutexes thatn need to be locked to ensure consistent
     * access to ampc[j] variable.
     */
    pthread_mutex_t * ampc_mutex;

    /**
     * @brief Same as <code>bfpc</code>, but the multiprecision
     * version.
     */
    mpc_t *bmpc;

    /**
     * @brief Mutexes that need to be locked to ensure consistent
     * access to bmpc[j] variable.
     */
    pthread_mutex_t * bmpc_mutex;

    /**
     * @brief Moduli of the floating point a_i
     * coefficients of the secular equation.
     */
    double *aafpc;

    /**
     * @brief Moduli of the floating point b_i 
     * coefficients of the secular equation.
     */
    double *abfpc;

    /**
     * @brief DPE Moduli of the CDPE of Multiprecision a_i 
     * coefficients of the secular equation.
     */
    rdpe_t *aadpc;
    
    /**
     * @brief DPE Moduli of the CDPE of Multiprecision b_i 
     * coefficients of the secular equation.
     */
    rdpe_t *abdpc;

    /**
     * @brief Initial multiprecision coefficients saved for latter
     * regeneration in <code>mps_secular_ga_regenerate_coefficients()</code>.
     */
    mpc_t *initial_ampc;

    /**
     * @brief Initial multiprecision coefficients saved for latter
     * regeneration in <code>mps_secular_ga_regenerate_coefficients()</code>.
     */
    mpc_t *initial_bmpc;

    /**
     * @brief Initial rational coefficients, if rational input is selected.
     * This value is the real part of the \f$a_i\f$ coefficients.
     */
    mpq_t *initial_ampqrc;

    /**
     * @brief Initial rational coefficients, if rational input is selected.
     * This value is the real part of the \f$b_i\f$ coefficients.
     */    
    mpq_t *initial_bmpqrc;

    /**
     * @brief Initial rational coefficients, if rational input is selected.
     * This value is the imaginary part of the \f$a_i\f$ coefficients.
     */
    mpq_t *initial_ampqic;

    /**
     * @brief Initial rational coefficients, if rational input is selected.
     * This value is the imaginary part of the \f$b_i\f$ coefficients.
     */
    mpq_t *initial_bmpqic;

    /**
     * @brief This mutex is locked while changing precision. 
     */
    pthread_mutex_t precision_mutex;

  };       /* End of struct mps_secular_equation {... */

  /**
   * @brief This is a struct that represent an iteration on a root. It contains
   * information that could be useful for mps_secular_*iterate() routine to determine
   * some error bound and provide a method for the routine to communicate if
   * it was able to set the radius or not (by setting the <code>radius_set</code> 
   * in the right way).
   */
  struct mps_secular_iteration_data {
    /**
     * @brief The index of the roots on which the iterations
     * is being carried out.
     */
    long int k;

    /**
     * @brief The state of the iteration. This is a pointer
     * to a boolean that tells if the iterator has been
     * able to set a radius or not, because the radius
     * that was there before was better.
     */
    mps_boolean radius_set;

    /**
     * @brief Global mutex used to synchronization, but mainly
     * while testing new MP implementations.
     */
    pthread_mutex_t * gs_mutex;
    
    /**
     * @brief Thread local copy of the \f$a_i\f$ coefficients of the secular
     * equation.
     */
    mpc_t * local_ampc;

    /**
     * @brief Thread local copy of the \f$b_i\f$ coefficients of the secular
     * equation.
     */
    mpc_t * local_bmpc;

    /**
     * @brief Thread local copy of the floating point coefficients of the secular
     * equation.
     */
    cplx_t * local_afpc;

    /**
     * @brief Thread local copy of the floating point coefficients of the secular
     * equation.
     */
    cplx_t * local_bfpc;

    /**
     * @brief Thread local copy of the CDPE coefficients of the secular
     * equation.
     */
    cdpe_t * local_adpc;

    /**
     * @brief Thread local copy of the CDPE coefficients of the secular
     * equation.
     */
    cdpe_t * local_bdpc;
  };

#endif /* #ifdef _MPS_PRIVATE */    


  /* MACROS */
#define mps_secular_equation_from_status(s) (mps_secular_equation*) (s)->secular_equation

  /* Routines in secular-newton.c */
  void mps_secular_fnewton (mps_context * st, mps_polynomial * p, mps_approximation * root, cplx_t corr);
  void mps_secular_dnewton (mps_context * st, mps_polynomial * p, mps_approximation * root, cdpe_t corr);
  void mps_secular_mnewton (mps_context * st, mps_polynomial * p, mps_approximation * root, mpc_t corr);

  /* Routines in secular-regeneartion.c */
  mps_boolean * mps_secular_ga_find_changed_roots (mps_context * s, cdpe_t * old_b, mpc_t * old_mb);

  mps_boolean mps_secular_ga_regenerate_coefficients_mp (mps_context * s, cdpe_t * old_b, mpc_t * old_mb);

  mps_boolean mps_secular_ga_regenerate_coefficients_monomial (mps_context * s, cdpe_t * old_b, mpc_t * old_mb, mps_boolean * root_changed);

  mps_boolean mps_secular_ga_regenerate_coefficients_secular (mps_context * s, cdpe_t * old_b, mpc_t * old_mb, mps_boolean * root_changed);

  mps_boolean mps_secular_ga_regenerate_coefficients (mps_context * s);

  /* Routines in secular.c */
  void mps_secular_deflate (mps_context * s, mps_secular_equation * sec);

  void mps_secular_check_data (mps_context * s, char *which_case);

  void mps_secular_restart (mps_context * s);

  void mps_secular_switch_phase (mps_context * s, mps_phase phase);

  long int mps_secular_raise_coefficient_precision (mps_context * s, mps_polynomial * p, long int wp);

  void mps_secular_raise_precision (mps_context * s, int wp);

  void mps_secular_raise_root_precision (mps_context * s, int wp);

  /* Routines in secular-starting.c */
  void mps_secular_fstart (mps_context * s, mps_secular_equation * sec);
  void mps_secular_dstart (mps_context * s, mps_secular_equation * sec);
  void mps_secular_mstart (mps_context * s, mps_secular_equation * sec);

  /* Routines in secular-iteration.c */
  int mps_secular_ga_fiterate (mps_context * s, int maxit, mps_boolean just_regenerated);

  int mps_secular_ga_diterate (mps_context * s, int maxit, mps_boolean just_regenerated);

  int mps_secular_ga_miterate (mps_context * s, int maxit, mps_boolean just_regenerated);
  
  /* Routines in secular-ga.c */
  mps_boolean mps_secular_ga_check_stop (mps_context * s);

  void mps_secular_ga_mpsolve (mps_context * s);

  void mps_secular_ga_update_coefficients (mps_context * s);

  /* Interface functions in secular.c */
  mps_secular_equation *mps_secular_equation_new (mps_context * s,
                                                  cplx_t * afpc,
                                                  cplx_t * bfpc,
                                                  unsigned long int n);

  mps_secular_equation *mps_secular_equation_new_raw (mps_context * s,
                                                      unsigned long int n);

  void mps_secular_equation_free (mps_context * ctx, mps_polynomial * p);

  void mps_secular_set_radii (mps_context * s);

  mps_boolean mps_secular_poly_feval_with_error (mps_context * ctx, mps_polynomial * p, cplx_t x, cplx_t value, double * error);

  mps_boolean mps_secular_poly_deval_with_error (mps_context * ctx, mps_polynomial * p, cdpe_t x, cdpe_t value, rdpe_t error);

  mps_boolean mps_secular_poly_meval_with_error (mps_context * ctx, mps_polynomial * p, mpc_t x, mpc_t value, rdpe_t error);

  void mps_secular_poly_fstart (mps_context * ctx, mps_polynomial * p);

  void mps_secular_poly_dstart (mps_context * ctx, mps_polynomial * p);

  void mps_secular_poly_mstart (mps_context * ctx, mps_polynomial * p);

#ifdef  __cplusplus
}
#endif

#endif                          /* SECULAR_H */
