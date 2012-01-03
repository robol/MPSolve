#ifndef MPS_STATUS_H
#define MPS_STATUS_H

#ifdef __cplusplus
extern "C"
{
#endif

  /**
   * @file
   * @brief This file contains the definition of <code>mps_status</code> and 
   * most of its fields.
   */
#include <mps/mps.h>

  /**
   * @brief Function that computes \f$\frac{p}{p'}\f$ (floating point version)
   */
  typedef void (*mps_fnewton_ptr) (mps_status * status, cplx_t, double *, cplx_t,
                                   mps_boolean *, void * user_data,
				   mps_boolean skip_radius_check);

  /**
   * @brief Function that computes \f$\frac{p}{p'}\f$ (dpe version)
   */
  typedef void (*mps_dnewton_ptr) (mps_status * status, cdpe_t x, rdpe_t rad,
                                   cdpe_t corr, mps_boolean * again, 
				   void * user_data,
				   mps_boolean skip_radius_check);

  /**
   * @brief Function that computes \f$\frac{p}{p'}\f$ (multiprecision version)
   */
  typedef void (*mps_mnewton_ptr) (mps_status * status, mpc_t x, rdpe_t rad,
                                   mpc_t corr, mps_boolean * again,
				   void * user_data,
				   mps_boolean skip_radius_check);

  /**
   * @brief Functions that check if float phase is needed or not and set
   * which_case accordingly to <code>'f'</code> or <code>'d'</code>.
   */
  typedef void (*mps_check_data_ptr) (mps_status *status, char *which_case);

  /**
   * @brief Function to dispose starting approximations in the case of
   * floating point iterations.
   */
  typedef void (*mps_fstart_ptr) (mps_status *status, int n, mps_cluster_item * cluster,
                                  double clust_rad, double g, rdpe_t eps);

  /**
   * @brief Function to dispose starting approximations in the case of
   * DPE iterations.
   */
  typedef void (*mps_dstart_ptr) (mps_status *status, int n, mps_cluster_item * cluster,
                                  rdpe_t clust_rad, rdpe_t g, rdpe_t eps);

  /**
   * @brief Function that computes radii to perform cluster analysis on the
   * roots in the floating point iterations.
   */
  typedef void (*mps_fradii_ptr) (mps_status * status);

  /**
   * @brief Routine that performs the computation loop to solve the polynomial
   * or the secular equation
   */
  typedef void (*mps_mpsolve_ptr) (mps_status *status);


  /*
   * Macros for casting user functions
   */
#define MPS_FNEWTON_PTR(x) (mps_fnewton_ptr) &(x)
#define MPS_DNEWTON_PTR(x) (mps_dnewton_ptr) &(x)
#define MPS_MNEWTON_PTR(x) (mps_mnewton_ptr) &(x)
#define MPS_CHECK_DATA_PTR(x) (mps_check_data_ptr) &(x)
#define MPS_FSTART_PTR(x) (mps_fstart_ptr) &(x)
#define MPS_DSTART_PTR(x) (mps_dstart_ptr) &(x)
#define MPS_MPSOLVE_PTR(x) (mps_mpsolve_ptr) &(x)

  /**
   * @brief this struct holds the state of the mps computation
   */
  struct mps_status
  {

    mps_boolean resume;         /* to complete                         */
    mps_boolean chkrad;         /* check radii after completion        */

    /**
     * @brief Byte containing the flags of debug enabled.
     */
    unsigned int debug_level;

    /**
     * @brief True if the computation has reached the maximum allowed precision.
     */
    mps_boolean over_max;

    /**
     * @brief Configuration of the input of MPSolve
     */
    mps_input_configuration * input_config;

    /**
     * @brief Output configuration for MPSolve.
     */
    mps_output_configuration * output_config;

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
    FILE *instr;

    /**
     * @brief Default output stream.
     */
    FILE *outstr;

    /**
     * @brief Default log stream
     */
    FILE *logstr;

    FILE *rtstr;

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
     * @brief Maximum precision reached during the computation.
     */
    long int data_prec_max;

    /**
     * @brief True if this is the first iteration after the precision has been 
     * raised.
     */
    mps_boolean just_raised_precision;

    /**
     * @brief mps_boolean value that determine if we should
     * use a random seed for starting points
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
    int zero_roots;             /* number of roots = 0                 */

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
    char (*status)[3];          /* status of each approximation        */

    /**
     * @brief Output index order
     */
    int *order;

    /**
     * @brief Root approximations as floating points complex
     * numbers.
     */
    cplx_t *froot;

    /**
     * @brief Root approximations as complex dpe numbers.
     */
    cdpe_t *droot;

    /**
     * @brief Root approsimations as complex multiprecision
     * numbers.
     */
    mpc_t *mroot;

    /**
     * @brief Radii of inclusion disks as real numbers.
     */
    double *frad;

    /**
     * @brief Radii of inclusion disks as dpe numbers.
     */
    rdpe_t *drad;

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

    mps_clusterization * clusterization;

    /* /\** */
    /*  * @brief Number of active clusters. */
    /*  *\/ */
    /* int nclust; */

    /* /\** */
    /*  * @brief This <code>int</code> array keep information about */
    /*  * semi-converged roots that were removed by a cluster to */
    /*  * improve convergence speed. */
    /*  * */
    /*  * If <code>s->clust_detached[i]</code> */
    /*  * is not zero, than the <b>unique</b> root in the <code>i</code>-th */
    /*  * cluster has been detached by the cluster whose index is the value */
    /*  * of <code>s->clust_detached</code>. */
    /*  *\/ */
    /* int *clust_detached; */

    /* /\** */
    /*  * @brief indices of cluster components */
    /*  * */
    /*  * <code>clust</code> is an integer array containing the indexes */
    /*  * of roots in every cluster. */
    /*  * */
    /*  * The indexes of the <code>j+1</code>th cluster are */
    /*  * <code>clust[punt[j] : punt[j+1]]</code>. */
    /*  * */
    /*  * @see punt */
    /*  *\/ */
    /* int *clust; */

    /* /\** */
    /*  * @brief begginning of each cluster */
    /*  * */
    /*  * <code>punt</code> is a vector of <code>nclust + 1</code> integers; */
    /*  * */
    /*  * For every <code>i</code> <code>punt[i]</code> is the index in */
    /*  * the integer vector <code>clust</code> corresponding to the first */
    /*  * index of a cluster, i.e. the (j+1)-th cluster of roots is composed by */
    /*  * roots indexed on <code>clust[p[j] : p[j+1]]</code>. */
    /*  * */
    /*  * @see nclust */
    /*  * @see clust */
    /*  *\/ */
    /* int *punt; */

    /**
     * @brief Array containing working precisions used for each root.
     */
    long int *rootwp;

    /**
     * @brief Array that whose i-th component is set to <code>true</code> if
     * the i-th root needs more iterations.
     */
    mps_boolean *again;

    /* /\** */
    /*  * @brief Array containing standard complex coefficients */
    /*  *\/ */
    /* cplx_t *fppc; */

    /**
     * @brief Standard complex coefficients of the polynomial.
     *
     * This is used as a temporary vector while shifting the polynomial
     * with a new gravity center in <code>mps_fshift()</code>.
     */
    cplx_t *fppc1;

    /**
     * @brief <code>dpe</code> complex coefficients of the polynomial.
     *
     * This is used as a temporary vector while shifting the polynomial
     * with a new gravity center in <code>mps_dshift()</code>.
     */
    cdpe_t *dpc1;

    /**
     * @brief <code>dpe</code> complex coefficients of the polynomial.
     *
     * This is used as a temporary vector while shifting the polynomial
     * with a new gravity center in <code>mps_dshift()</code>.
     */
    cdpe_t *dpc2;

    /**
     * @brief Multiprecision complex coefficients of the polynomial.
     *
     * This is used as a temporary vector while shifting the polynomial
     * with a new gravity center in <code>mps_mshift()</code>.
     */
    mpc_t *mfpc1;

    /**
     * @brief Multiprecision complex coefficients of the polynomial.
     *
     * This is used as a temporary vector while shifting the polynomial
     * with a new gravity center in <code>mps_mshift()</code>.
     */
    mpc_t *mfpc2;

    /**
     * @brief Multiprecision complex coefficients of the
     * first derivative of the polynomial.
     *
     * This is used as a temporary vector while shifting the polynomial
     * with a new gravity center in <code>mps_mshift()</code>.
     */
    mpc_t *mfppc1;

    /**
     * @brief Vector representing sparsity of the polynomial in the
     * same way that <code>spar</code> does.
     *
     * It is used as a temporary vector.
     *
     * @see spar
     */
    mps_boolean *spar1;

    /**
     * @brief Vector representing sparsity of the polynomial in the
     * same way that <code>spar</code> does.
     *
     * It is used as a temporary vector.
     *
     * @see spar
     */
    mps_boolean *spar2;

    /**
     * @brief Old value of <code>punt</code> (temporary vector).
     *
     * @see punt
     */
    int *oldpunt;

    /**
     * @brief Vector containing the moduli of the coefficients
     * of the polynomial as floating point numbers.
     *
     * It is used in the computation of the newton polygonal in
     * <code>mps_fcompute_starting_radii()</code>.
     *
     * @see mps_fcompute_starting_radii()
     */
    double *fap1;

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
    double *fap2;

    /**
     * @brief Vector containing the moduli of the coefficients
     * of the polynomial as <code>dpe</code> numbers.
     *
     * It is used in the computation of the newton polygonal in
     * <code>mps_dcompute_starting_radii()</code>.
     *
     * @see mps_dcompute_starting_radii()
     */
    rdpe_t *dap1;

    /**
     * @brief Vector needed for convex hull computation.
     *
     * It is <code>true</code> in position \f$j\f$ if
     * the point \f$(j, log(x_j))\f$ is a vertex of the convex
     * hull computed by <code>fconvex()</code> and the other functions in
     * <code>mps_cnvx.c</code>
     */
    mps_boolean *h;

    /**
     * @brief Temporary vector containing the old value of
     * <code>again</code>.
     *
     * @see again
     */
    mps_boolean *again_old;

    int *clust_aux;             /* auxiliary vector                    */
    int *punt_aux;              /* auxiliary vector                    */
    int *punt_out;              /* auxiliary vector                    */
    int *clust_out;             /* auxiliary vector                    */

    /**
     * @brief The number of circles with initial approximations.
     */
    int n_radii;

    /**
     * @brief This variable is used to store the radii of the
     * circles with initial approximations.
     */
    double *fradii;

    /**
     * @brief This variable is used to store the radii of the
     * circles with initial approximations.
     */
    rdpe_t *dradii;

    /**
     * @brief This variable is used to store the partitioning
     * done when disposing initial approximations.
     */
    int *partitioning;

    /* SECTION -- Algorihtmm selection */

    /**
     * @brief This is used in the program to switch behavious based
     * on the algorithm that is been used now.
     */
    mps_algorithm algorithm;

    /**
     * @brief Pointer to the function to perform newton in floating
     * point implemented by the user.
     */
    /* void (*fnewton_usr) (mps_status *status, cplx_t, double *, cplx_t, */
    /*                      mps_boolean *, void * user_data,  */
    /* 			 mps_boolean * skip_radius_computation); */
    mps_fnewton_ptr fnewton_usr;

    /**
     * @brief Pointer to the function to perform newton in dpe
     * implemented by the user.
     */
    /* void (*dnewton_usr) (mps_status *status, cdpe_t x, rdpe_t rad, cdpe_t corr, */
    /*                      mps_boolean * again, void * user_data,  */
    /* 			 mps_boolean * skip_radius_computation); */
    mps_dnewton_ptr dnewton_usr;

    /**
     * @brief Pointer to the function to perform newton in multiprecision
     * implemented by the user.
     */
    /* void (*mnewton_usr) (mps_status *status, mpc_t x, rdpe_t rad, mpc_t corr, */
    /*                      mps_boolean * again, void * user_data,  */
    /* 			 mps_boolean * skip_radius_computation); */
    mps_mnewton_ptr mnewton_usr;

    /**
     * @brief Check data routine that has the task to determine if a float phase
     * can be performed or dpe are needed now.
     */
    void (*check_data_usr) (mps_status *status, char *which_case);

    /**
     * @brief Routine to dispose starting approximations provided by the user
     */
    void (*fstart_usr) (mps_status *status, int n, mps_cluster_item * cluster, double clust_rad,
                        double g, rdpe_t eps);

    /**
     * @brief Routine to dispose starting approximations provided
     * by user in the case of DPE computation.
     */
    void (*dstart_usr) (mps_status *status, int n, mps_cluster_item * cluster, rdpe_t clust_rad,
                        rdpe_t g, rdpe_t eps);

    /**
     * @brief Routine that performs the loop needed to coordinate
     * root finding. It has to be called to do the hard work.
     */
    void (*mpsolve_ptr) (mps_status *status);

    /**
     * @brief A pointer to the polynomial that is being solve or
     * NULL if there is no such monomial representation.
     */
    mps_monomial_poly * monomial_poly;

    /**
     * @brief A pointer that can be set to anything the user
     * would like to access during computations. It is meant to be
     * used when implementing fnewton, dnewton and mnewton
     * functions to provide additional data for the
     * computing of the polynomial.
     */
    mps_secular_equation * secular_equation;

    /**
     * @brief Number of threads to be spawned.
     */
    int n_threads;

    /**
     * @brief The thread pool used for the concurrent part of MPSolve.
     */
    mps_thread_pool * pool;

    /* DEBUG SECTIONS */

    unsigned long int regeneration_time;
    unsigned long int mp_iteration_time;
    unsigned long int dpe_iteration_time;
    unsigned long int fp_iteration_time;

  };                 /* End of typedef struct { ... */





#ifdef __cplusplus
}
#endif

#endif
