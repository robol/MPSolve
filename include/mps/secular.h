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
#define	 SECULAR_H

#ifdef	__cplusplus
extern "C"
{
#endif

#include <mps/mps.h>
#include <float.h>

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
     * @brief Size of the vectors of the coefficients of the
     * secular equation.
     */
    unsigned long int n;

    /**
     * @brief Selected starting case, can be 'd' for DPE
     * or 'f' for floating point
     */
    mps_phase starting_case;

    /**
     * @brief Set to true if the approximation are the best that
     * can be obtained with the current precision
     */
    mps_boolean best_approx;

    /**
     * @brief This vector contains the errors present in the coefficients
     * of the computed regeneration of the secular equation.
     *
     * This is the version use in DPE and MPC computation.
     */
    rdpe_t * dregeneration_epsilon;

    /**
     * @brief This vector contains the errors present in the coefficients
     * of the computed regeneration of the secular equation.
     *
     * This is the version use in floating point computation.
     */
    double * fregeneration_epsilon;    

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

    pthread_mutex_t * gs_mutex;
  };

#endif /* #ifdef _MPS_PRIVATE */    


  /* MACROS */
#define mps_secular_equation_from_status(s) (mps_secular_equation*) (s)->secular_equation

  /* Routines in secular-newton.c */
  void mps_secular_fnewton (mps_status * st, cplx_t x, double *rad, cplx_t corr,
			    mps_boolean * again, void * user_data, 
			    mps_boolean skip_radius_compuation);
  void mps_secular_dnewton (mps_status * st, cdpe_t x, rdpe_t rad, cdpe_t corr,
			    mps_boolean * again, void * user_data,
			    mps_boolean skip_radius_computation);
  void mps_secular_mnewton (mps_status * st, mpc_t x, rdpe_t rad, mpc_t corr,
			    mps_boolean * again, void * user_data,
			    mps_boolean skip_radius_computation);

  /* Routines in secular-regeneartion.c */
  mps_boolean * mps_secular_ga_find_changed_roots (mps_status * s, cdpe_t * old_b, mpc_t * old_mb);

  mps_boolean mps_secular_ga_regenerate_coefficients_mp (mps_status * s, cdpe_t * old_b, mpc_t * old_mb);

  mps_boolean mps_secular_ga_regenerate_coefficients_monomial (mps_status * s, cdpe_t * old_b, mpc_t * old_mb, mps_boolean * root_changed);

  mps_boolean mps_secular_ga_regenerate_coefficients_secular (mps_status * s, cdpe_t * old_b, mps_boolean * root_changed);

  mps_boolean mps_secular_ga_regenerate_coefficients (mps_status * s);

  /* Routines in secular.c */
  void mps_secular_deflate (mps_status * s, mps_secular_equation * sec);

  void mps_secular_check_data (mps_status * s, char *which_case);

  void mps_secular_switch_phase (mps_status * s, mps_phase phase);

  void mps_secular_raise_coefficient_precision (mps_status * s, int wp);

  void mps_secular_raise_precision (mps_status * s, int wp);

  void mps_secular_raise_root_precision (mps_status * s, int wp);

  /* Routines in secular-starting.c */
  void mps_secular_fstart (mps_status * s, int n, mps_cluster_item * cluster, double clust_rad,
			   double g, rdpe_t eps);
  void mps_secular_dstart (mps_status * s, int n, mps_cluster_item * cluster, rdpe_t clust_rad,
			   rdpe_t g, rdpe_t eps);
  void mps_secular_mstart (mps_status * s, int n, mps_cluster_item * cluster, rdpe_t clust_rad,
			   rdpe_t g, rdpe_t eps);

  /* Routines in secular-iteration.c */
  int mps_secular_ga_fiterate (mps_status * s, int maxit, mps_boolean just_regenerated);

  int mps_secular_ga_diterate (mps_status * s, int maxit, mps_boolean just_regenerated);

  int mps_secular_ga_miterate (mps_status * s, int maxit, mps_boolean just_regenerated);
  
  /* Routines in secular-ga.c */
  mps_boolean mps_secular_ga_check_stop (mps_status * s);

  void mps_secular_ga_improve (mps_status * s);

  void mps_secular_ga_mpsolve (mps_status * s);

  void mps_secular_ga_update_coefficients (mps_status * s);

  /* Interface functions in secular.c */
  mps_secular_equation *mps_secular_equation_new (mps_status * s,
                                                  cplx_t * afpc,
                                                  cplx_t * bfpc,
                                                  unsigned long int n);

  mps_secular_equation *mps_secular_equation_new_raw (mps_status * s,
                                                      unsigned long int n);

  void mps_secular_equation_free (mps_secular_equation * s);

  void mps_secular_set_radii (mps_status * s);

#ifdef	__cplusplus
}
#endif

#endif                          /* SECULAR_H */
