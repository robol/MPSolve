/* 
 * File:   secular.h
 * Author: leonardo
 *
 * Created on 29 aprile 2011, 16.26
 */

#ifndef SECULAR_H
#define	SECULAR_H

#ifdef	__cplusplus
extern "C" {
#endif

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
typedef struct {
    /**
     * @brief Vector of \f$a_i\f$ as complex floating
     * point numbers.
     */
    cplx_t* afpc;

    /**
     * @brief Same as <code>afpc</code>, but the <code>dpe</code>
     * version.
     */
    cdpe_t* adpc;

    /**
     * @brief Vector with the values of \f$b_i\f$ as complex
     * floating point numbers.
     */
    cplx_t* bfpc;

    /**
     * @brief Same as <code>bfpc</code>, but the <code>dpe</code>
     * version.
     */
    cdpe_t* bdpc;

    /**
     * @brief Same as <code>afpc</code>, but the multiprecision
     * version.
     */
    mpc_t * ampc;

    /**
     * @brief Same as <code>bfpc</code>, but the multiprecision
     * version.
     */
    mpc_t * bmpc;

    /**
     * @brief Size of the vectors of the coefficients of the
     * secular equation.
     */
    unsigned long int n;

} mps_secular_equation; /* End of typedef struct {... */


/* Routines in mps_secular.c */
void mps_secular_fnewton(mps_status* st, cplx_t x, double * rad, cplx_t corr, mps_boolean * again);
void mps_secular_dnewton(mps_status* st, cdpe_t x, rdpe_t rad, cdpe_t corr, mps_boolean * again);
void mps_secular_mnewton(mps_status* st, mpc_t x, rdpe_t rad, mpc_t corr, mps_boolean * again);
void mps_secular_check_data(mps_status* s, char* which_case);
void mps_secular_fstart(mps_status* s, int n, int i_clust, double clust_rad,
		double g, rdpe_t eps);
void mps_secular_dstart(mps_status* s, int n, int i_clust, rdpe_t clust_rad,
		rdpe_t g, rdpe_t eps);
void mps_secular_ga_mpsolve(mps_status* s, mps_phase phase);

/* Interface functions in mps_secular.c */
mps_secular_equation* mps_secular_equation_new(cplx_t* afpc, cplx_t* bfpc, unsigned long int n);
void mps_secular_equation_free(mps_secular_equation* s);

#ifdef	__cplusplus
}
#endif

#endif	/* SECULAR_H */

