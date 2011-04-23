/* 
 * File:   interface.h
 * Author: leonardo
 *
 * Created on 23 aprile 2011, 10.35
 */

#ifndef MPS_INTERFACE_H
#define	MPS_INTERFACE_H

#ifdef	__cplusplus
extern "C" {
#endif

/* Octave module workardound */
#ifdef __UNDEF_CPLUSPLUS
#undef __cplusplus
#endif

#include <gmp.h>
#include <mps/mt.h>
#include <mps/mpc.h>
#include <mps/tools.h>
#include <stdlib.h>
#include <stdio.h>


    /**
 * @brief Type representing the computation phase
 * of the algorithm we are in
 * now. It can assume the values:
 * - <code>no_phase</code>;
 * - <code>float_phase</code>;
 * - <code>dpe_phase</code>;
 * - <code>mp_phase</code>;
 */
typedef enum {
    no_phase, float_phase, dpe_phase, mp_phase
} mps_phase;

/**
 * @brief Function that computes \f$\frac{p}{p'}\f$ (floating point version)
 */
typedef void (*mps_fnewton_ptr)(void* status,cplx_t,double*,cplx_t,mps_boolean*);


/**
 * @brief Function that computes \f$\frac{p}{p'}\f$ (dpe version)
 */
typedef void (*mps_dnewton_ptr)(void* status, cdpe_t x, rdpe_t   rad, cdpe_t corr, mps_boolean* again);

/**
 * @brief Function that computes \f$\frac{p}{p'}\f$ (multiprecision version)
 */
typedef void (*mps_mnewton_ptr)(void* status, mpc_t  x, rdpe_t   rad, mpc_t  corr, mps_boolean* again);

/*
 * Macros for casting
 */
#define MPS_FNEWTON_PTR(x) (mps_fnewton_ptr) &x
#define MPS_DNEWTON_PTR(x) (mps_dnewton_ptr) &x
#define MPS_MNEWTON_PTR(x) (mps_mnewton_ptr) &x

/**
 * @brief this struct holds the state of the mps computation
 */
typedef struct {

   mps_boolean resume;		/* to complete                         */
   mps_boolean chkrad;		/* check radii after completion        */

   /**
    * @brief Newton isolation of the cluster.
    */
   int newtis;

   /**
    * @brief Old value for the newton isolation of the cluster.
    */
   int newtis_old;

  /*
   * INPUT / OUTPUT STREAMS
   */

  /**
   * @brief <code>true</code> if log are needed. They will
   * be written to <code>logstr</code>
   *
   * @see logstr
   */
   mps_boolean DOLOG;

  /**
   * @brief <code>true</code> if warning are needed.
   */
   mps_boolean DOWARN;

  /**
   * @brief <code>true</code> if root sorting is desired. It will
   * be performed with routines in <code>mps_sort.c</code>.
   */
   mps_boolean DOSORT;

  /**
   * @brief Default input stream.
   */
   FILE * instr;

  /**
   * @brief Default output stream.
   */
   FILE * outstr;

  /**
   * @brief Default log stream
   */
   FILE * logstr;

   FILE * rtstr;

  /*
   * CONSTANT, PARAMETERS
   */

  /**
   * @brief number of max packets of iterations
   */
   int max_pack;

  /**
   * @brief number of max iterations per packet
   */
   int max_it;

  /**
   * @brief Number of max newton iterations for gravity center
   * computations.
   */
   int max_newt_it;

  /**
   * @brief Maximum allowed number of bits for mp numbers: used in
   * high precision shift.
   */
   long int mpwp_max;

  /**
   * @brief stores the goal of the computation
   *
   * <code>goal</code> is an array of 5 chars with this meaning:
   * - <code>goal[0]</code> can assume the following values:
   *   - <code>a</code> that means \b approximate;
   *   - <code>c</code> that means \b count the roots;
   *   - <code>i</code> that means \b isolate the roots;
   * - <code>goal[1]</code>represents the search set for
   *   the roots. It can be:
   *   - <code>a</code> that means \f$ \mathbb{C} \f$;
   *   - <code>r</code> that means \f$ \{ z \in \mathbb{C} \ | \ \mathrm{Re}{z} > 0 \} \f$;
   *   - <code>l</code> that means \f$ \{ z \in \mathbb{C} \ | \ \mathrm{Re}{z} < 0 \} \f$;
   *   - <code>u</code> that means \f$ \{ z \in \mathbb{C} \ | \ \mathrm{Im}{z} > 0 \} \f$;
   *   - <code>d</code> that means \f$ \{ z \in \mathbb{C} \ | \ \mathrm{Im}{z} < 0 \} \f$;
   *   - <code>i</code> that means \f$ \{ z \in \mathbb{C} \ | \ \lVert z \rVert \leq 1 \} \f$;
   *   - <code>i</code> that means \f$ \{ z \in \mathbb{C} \ | \ \lVert z \rVert \geq 1 \} \f$;
   * - <code>goal[2]</code> controls if multiplicity check is enabled or not. Its value can be:
   *   - <code>m</code> if it is enabled;
   *   - <code>n</code> if that's not the case;
   * - <code>goal[3]</code> sets what properties of the roots have to be detected
   *   during the computation:
   *   - <code>n</code> means no one, and this is the default;
   *   - <code>r</code> means detect if a root is real;
   *   - <code>i</code> means detect if a root is imaginary;
   *   - <code>b</code> means both the above;
   * - <code>goal[4]</code> determines output format for the results:
   *   - <code>b</code> means \b bare and it's the default value;
   *   - <code>g</code> means \b gnuplot format;
   *   - <code>c</code> means \b compact;
   *   - <code>v</code> means \b verbose;
   *   - <code>f</code> means \b full;
   */
   char goal[5];

  /**
   * @brief digits of required output precision
   *
   */
   long int prec_out;

  /**
   * @brief mps_boolean value that determine if we should
   * use a random seed for startingd points
   */
   mps_boolean random_seed;

  /*
   * POLYNOMIAL DATA: SHARED VARIABLES
   */

  /**
   * @brief degree of zero-deflated polynomial.
   */
   int n;

  /**
   * @brief input degree and allocation size.
   */
   int deg;

  /**
   * @brief stores the input data type
   *
   * The value of <code>data_type[0]</code> can be:
   * - <code>'s'</code> if the input is a sparse polynomial;
   * - <code>'u'</code> if the input is a user polynomial;
   * - <code>'d'</code> if the input is a dense polybomial;
   *
   * while the value of <code>data_type[1]</code> can be:
   * - <code>'r'</code> if the coefficients are real;
   * - <code>'c'</code> if the coefficents are complex;
   *
   * and finally <code>data_type[2]</code> can assume the following values:
   * - <code>'i'</code> means integer coefficients;
   * - <code>'q'</code> means rationa cofficients;
   * - <code>'b'</code> means bigfloat coefficents;
   * - <code>'f'</code> means floating point coefficients;
   */
   char * data_type;

  /**
   * @brief Number of digits of input precision in its binary
   * representation.
   */
   long int prec_in;

  /**
   * @brief This array contains the structure of the sparse
   * polynomial.
   *
   * <code>spar[i]</code> is <code>true</code> if and only if
   * the i-th coefficients of the polynomial is a non-zero
   * coefficients
   */
   mps_boolean * spar;

   /**
    * @brief Standard real coefficients.
    */
   double * fpr;

   /**
    * @brief Standard complex coefficients.
    */
   cplx_t * fpc;

   /**
    * @brief Dpe real coefficients.
    */
   rdpe_t * dpr;

   /**
    * @brief Dpe complex coefficients.
    */
   cdpe_t * dpc;

   /**
    * @brief Real part of the integer input coefficients.
    */
   mpz_t * mip_r;

   /**
    * @brief Imaginary part of the integer input coefficients.
    */
   mpz_t * mip_i;

   /**
    * @brief Real part of rational input coefficients.
    */
   mpq_t * mqp_r;

   /**
    * @brief Imaginary part of rational input coefficients.
    */
   mpq_t * mqp_i;

   /**
    * @brief Multiprecision real coefficients.
    */
   mpf_t * mfpr;

   /**
    * @brief Multiprecision complex coefficients.
    */
   mpc_t * mfpc;

   /* Solution related variables */
   /**
    * @brief Last computing phase.
    */
   mps_phase lastphase;

   /**
    * @brief shift in the angle in the positioning of the
    * starting approximation for the last cluster. It will
    * be used to determine the new sigma to maximize distance
    * between starting points.
    */
   double last_sigma;

  /**
   * @brief Vector containing count of in, out and uncertaing roots.
   */
   int count[3];

  /**
   * @brief Number of zero roots.
   */
   int zero_roots;		/* number of roots = 0                 */

  /**
   * @brief Status of each approximation
   *
   * <code>status</code> is an array of char arrays
   * composed of three elements. For every <code>i</code>:
   * - <code>status[i][0]</code> can assume the values:
   *   - <code>m</code>: multiple root;
   *   - <code>i</code>: isolated root;
   *   - <code>a</code>: approximated single root (relative error
   *     is less then \f$ 10^{-d_{out}} \f$);
   *   - <code>o</code>: approximated cluster of roots (relative
   *     error is less than \f$ 10^{-d_{out}} \f$);
   *   - <code>c</code>: cluster of roots not yet approximated (relative
   *     error is greater than \f$ 10^{-d_{out}} \f$);
   *   - <code>f</code>: radius has reached the extreme values representable as <code>rdpe_t</code>;
   *   - <code>x</code>: this root cannot be represented as a <code>double</code>, i.e. it
   *     is <code>< DBL_MIN</code> or <code>> DBL_MAX</code>;
   * - <code>status[i][1]</code> can assume the values:
   *   - <code>R</code>: real root;
   *   - <code>r</code>: non real root;
   *   - <code>I</code>: imaginary root;
   *   - <code>i</code>: non imaginary root;
   *   - <code>w</code>: uncertain real/imaginary root;
   *   - <code>z</code>: non real and non imaginary root;
   * - <code>status[i][2]</code> can assume the values:
   *   - <code>i</code>: root in \f$ \mathcal{S} \f$;
   *   - <code>o</code>: root out of \f$ \mathcal{S} \f$;
   *   - <code>u</code>: root uncertain;
   */
   char (* status)[3];	/* status of each approximation        */

   /**
    * @brief Output index order
    */
   int * order;

   /**
    * @brief Root approximations as floating points complex
    * numbers.
    */
   cplx_t * froot;

   /**
    * @brief Root approximations as complex dpe numbers.
    */
   cdpe_t * droot;

   /**
    * @brief Root approsimations as complex multiprecision
    * numbers.
    */
   mpc_t * mroot;

   /**
    * @brief Radii of inclusion disks as real numbers.
    */
   double * frad;

   /**
    * @brief Radii of inclusion disks as dpe numbers.
    */
   rdpe_t * drad;

  /* lifetime global variables */

   /**
    * @brief <code>true</code> if the float phase should be skipped,
    * passing directly do dpe phase.
    */
   mps_boolean skip_float;

   /**
    * @brief Input precision of the coefficients.
    */
   rdpe_t eps_in;

   /**
    * @brief Output precision of the roots.
    */
   rdpe_t eps_out;

   /**
    * @brief Logarithm of the max modulus of the coefficients.
    */
   double lmax_coeff;

   /**
    * @brief Bits of working precision that mpsolve is using.
    */
   long int mpwp;

   /**
    * @brief Current multiprecision epsilon.
    */
   rdpe_t mp_epsilon;

   /**
    * @brief Log of the lower bound to the minumum distance of the roots
    */
   double sep;

  /**
   * @brief Number of active clusters.
   */
   int nclust;

   /**
    * @brief This <code>int</code> array keep information about
    * semi-converged roots that were removed by a cluster to
    * improve convergence speed.
    *
    * If <code>s->clust_detached[i]</code>
    * is not zero, than the <b>unique</b> root in the <code>i</code>-th
    * cluster has been detached by the cluster whose index is the value
    * of <code>s->clust_detached</code>.
    */
   int* clust_detached;

  /**
   * @brief indices of cluster components
   *
   * <code>clust</code> is an integer array containing the indexes
   * of roots in every cluster.
   *
   * The indexes of the <code>j+1</code>th cluster are
   * <code>clust[punt[j] : punt[j+1]]
   *
   * @see punt
   */
   int * clust;

  /**
   * @brief begginning of each cluster
   *
   * <code>punt</code> is a vector of <code>nclust + 1</code> integers;
   *
   * For every <code>i</code> <code>punt[i]</code> is the index in
   * the integer vector <code>clust</code> corresponding to the first
   * index of a cluster, i.e. the (j+1)-th cluster of roots is composed by
   * roots indexed on <code>clust[p[j] : p[j+1]]</code>.
   *
   * @see nclust
   * @see clust
   */
   int * punt;

  /**
   * @brief Array containing working precisions used for each root.
   */
   long int * rootwp;

  /**
   * @brief Multiprecision complex coefficients of \f$p'(x)\f$.
   */
   mpc_t * mfppc;

  /**
   * @brief Array containing moduli of the coefficients as double numbers.
   */
   double * fap;

  /**
   * @brief Array containing moduli of the coefficients as dpe numbers.
   */
   rdpe_t * dap;

  /**
   * @brief Array that whose i-th component is set to <code>true</code> if
   * the i-th root needs more iterations.
   */
   mps_boolean * again;

  /**
   * @brief Array containing standard complex coefficients
   */
   cplx_t * fppc;

   /**
    * @brief Standard complex coefficients of the polynomial.
    *
    * This is used as a temporary vector while shifting the polynomial
    * with a new gravity center in <code>mps_fshift()</code>.
    */
   cplx_t * fppc1;

   /**
    * @brief <code>dpe</code> complex coefficients of the polynomial.
    *
    * This is used as a temporary vector while shifting the polynomial
    * with a new gravity center in <code>mps_dshift()</code>.
    */
   cdpe_t * dpc1;

   /**
    * @brief <code>dpe</code> complex coefficients of the polynomial.
    *
    * This is used as a temporary vector while shifting the polynomial
    * with a new gravity center in <code>mps_dshift()</code>.
    */
   cdpe_t * dpc2;

   /**
    * @brief Multiprecision complex coefficients of the polynomial.
	*
    * This is used as a temporary vector while shifting the polynomial
    * with a new gravity center in <code>mps_mshift()</code>.
    */
   mpc_t * mfpc1;

   /**
    * @brief Multiprecision complex coefficients of the polynomial.
	*
    * This is used as a temporary vector while shifting the polynomial
    * with a new gravity center in <code>mps_mshift()</code>.
    */
   mpc_t * mfpc2;

   /**
    * @brief Multiprecision complex coefficients of the
    * first derivative of the polynomial.
	*
    * This is used as a temporary vector while shifting the polynomial
    * with a new gravity center in <code>mps_mshift()</code>.
    */
   mpc_t * mfppc1;

   /**
    * @brief Vector representing sparsity of the polynomial in the
    * same way that <code>spar</code> does.
    *
    * It is used as a temporary vector.
    *
    * @see spar
    */
   mps_boolean * spar1;

   /**
    * @brief Vector representing sparsity of the polynomial in the
    * same way that <code>spar</code> does.
    *
    * It is used as a temporary vector.
    *
    * @see spar
    */
   mps_boolean * spar2;

   /**
    * @brief Old value of <code>punt</code> (temporary vector).
    *
    * @see punt
    */
   int * oldpunt;

   /**
    * @brief Vector containing the moduli of the coefficients
    * of the polynomial as floating point numbers.
    *
    * It is used in the computation of the newton polygonal in
    * <code>mps_fcompute_starting_radii()</code>.
    *
    * @see mps_fcompute_starting_radii()
    */
   double * fap1;

   /**
    * @brief Vector containing the logarithm of the moduli of
    * the coefficients of the polynomial as floating
    * point numbers.
    *
    * It is used in the computation of the newton polygonal in
    * <code>mps_fcompute_starting_radii()</code>.
    *
    * @see mps_fcompute_starting_radii()
    * @see fap1
    */
   double * fap2;

   /**
    * @brief Vector containing the moduli of the coefficients
    * of the polynomial as <code>dpe</code> numbers.
    *
    * It is used in the computation of the newton polygonal in
    * <code>mps_dcompute_starting_radii()</code>.
    *
    * @see mps_dcompute_starting_radii()
    */
   rdpe_t * dap1;

   /**
    * @brief Vector containing the logarithms of the moduli of
    * the coefficients
    * of the polynomial as <code>dpe</code> numbers.
    *
    * It is used in the computation of the newton polygonal in
    * <code>mps_fcompute_starting_radii()</code>.
    *
    * @see mps_fcompute_starting_radii()
    */
   rdpe_t * dap2;

  /**
   * @brief Vector needed for convex hull computation.
   *
   * It is <code>true</code> in position \f$j\f$ if
   * the point \f$(j, log(x_j))\f$ is a vertex of the convex
   * hull computed by <code>fconvex()</code> and the other functions in
   * <code>mps_cnvx.c</code>
   */
   mps_boolean * h;

   /**
    * @brief Temporary vector containing the old value of
    * <code>again</code>.
    *
    * @see again
    */
   mps_boolean * again_old;

   int * clust_aux;		/* auxiliary vector                    */
   int * punt_aux;		/* auxiliary vector                    */
   int * punt_out;		/* auxiliary vector                    */
   int * clust_out;		/* auxiliary vector                    */

   /**
    * @brief The number of circles with initial approximations.
    */
   int n_radii;

   /**
    * @brief This variable is used to store the radii of the
    * circles with initial approximations.
    */
   double* fradii;

   /**
    * @brief This variable is used to store the radii of the
    * circles with initial approximations.
    */
   rdpe_t* dradii;

   /**
    * @brief This variable is used to store the partitioning
    * done when disposing initial approximations.
    */
   int* partitioning;


   /**
    * @brief Pointer to the function to perform newton in floating
    * point implemented by the user.
    */
   void (*fnewton_usr)(void* status,cplx_t,double*,cplx_t,mps_boolean*);

   /**
    * @brief Pointer to the function to perform newton in dpe
    * implemented by the user.
    */
   void (*dnewton_usr)(void* status, cdpe_t x, rdpe_t   rad, cdpe_t corr, mps_boolean* again);

   /**
    * @brief Pointer to the function to perform newton in multiprecision
    * implemented by the user.
    */
   void (*mnewton_usr)(void* status, mpc_t  x, rdpe_t   rad, mpc_t  corr, mps_boolean* again);

   /**
    * @brief Check data routine that has the task to determine if a float phase
    * can be performed or dpe are needed now.
    */
   void (*check_data_usr)(void* status, char* which_case);

	/**
	 * @brief A pointer that can be set to anything the user
	 * would like to access during computations. It is meant to be
	 * used when implementing fnewton, dnewton and mnewton
	 * functions to provide additional data for the
	 * computing of the polynomial.
	 */
	void * user_data;

#ifndef DISABLE_DEBUG
    /**
     * @brief Counter for the cycles used in debug prints,
     * that is allocated only if debug is allowed.
     */
    int debug_i;

    /**
     * @brief Counter for the cycles used in debug prints,
     * that is allocated only if debug is allowed.
     */
    int debug_j;
#endif

} mps_status; /* End of typedef struct { ... */

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


/*
 *    ====== ROUTINES EXPOSED TO THE INTERFACE ======
 */

/* functions in mps_defaults.c */
void mps_set_default_values(mps_status* s);

/* Functions in mps_main.c */
void mps_mpsolve(mps_status* s);

/* functions in mps_interface.c */
mps_status* mps_status_new();
void mps_status_free(mps_status* s);
int mps_status_set_poly_d(mps_status* s, cplx_t* coeff, long unsigned int n);
int mps_status_set_poly_i(mps_status* s, int* coeff, long unsigned int n);
int mps_status_get_roots_d(mps_status* s, cplx_t* roots, double* radius);


#ifdef	__cplusplus
}
#endif

#ifdef __UNDEF_CPLUSPLUS
}
#endif

#endif	/* MPS_INTERFACE_H */

