#ifndef MPS_STATUS_H
#define MPS_STATUS_H

#include <mps/mps.h>
#include <pthread.h>

#ifdef __cplusplus
extern "C"
{
#endif

  /**
   * @file
   * @brief This file contains the definition of <code>mps_context</code> and 
   * most of its fields.
   */

  /**
   * @brief Routine that performs the computation loop to solve the polynomial
   * or the secular equation
   */
  typedef void (*mps_mpsolve_ptr) (mps_context *status);

  /**
   * @brief Pointer to the callback for the async version of mpsolve
   */
  typedef void* (*mps_callback) (mps_context * status, void * user_data);


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

#ifdef _MPS_PRIVATE
  /**
   * @brief this struct holds the state of the mps computation
   */
  struct mps_context
  {
    /**
     * @brief true if an error has occurred during the computation.
     */
    mps_boolean error_state;

    /**
     * @brief The text describing the last error occurred.
     */
    char * last_error;

    /**
     * @brief This value is non NULL if mpsolve is launched via mps_mpsolve_async()
     * and in that case holds a pointer to the thread pool used to manage
     * asynchronous callbacks. 
     *
     * It will be automatically freed by mps_free_data().
     */
    mps_thread_pool * self_thread_pool;

    /**
     * @brief true if we are trying to resume previously interrupted.
     *
     * Not yet implemented.
     */
    mps_boolean resume;

    /**
     * @brief True if check of radius should be performed at the end
     * of the algorithm.
     *
     * Only works for algorithm MPS_ALGORITHM_SECULAR_GA.
     */
    mps_boolean chkrad;

    /**
     * @brief Callback called when the async version of mps_mpsolve(), i.e.
     * terminate the computation.
     */
    mps_callback callback;

    /**
     * @brief Pointer to user_data passed to the callback.
     */
    void * user_data;

    /**
     * @brief The operation running now. Can be used to debug what's happening
     * event if mpsolve was launched without debug enabled.
     */
    mps_operation operation;

    /**
     * @brief This value is true if the data for the computation has been allocated
     * by calling mps_allocate_data(). It is used by mps_free_data() to know what has
     * to be freed.
     */
    mps_boolean initialized;

    /**
     * @brief Bytes containing the flags of debug enabled.
     */
    mps_debug_level debug_level;

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
    mps_long_int_mt data_prec_max;

    /**
     * @brief Precision operation give best results when done one
     * thread at a time :)
     */
    pthread_mutex_t precision_mutex;

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

    /* Solution related variables */

    /**
     * @brief Last computing phase.
     */
    mps_phase lastphase;

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
    int zero_roots;

    /**
     * @brief Output index order
     */
    int *order;

    /**
     * @brief Vector of points to the 
     * current root approximations. 
     */
    mps_approximation ** root;

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
     * @brief Clusterization object that represent the clusterization
     * detected on the roots.
     *
     * This value is updated with the <code>mps_*cluster</code>
     * routines.
     *
     * @see mps_fcluster(), mps_dcluster(), mps_mcluster()
     */
    mps_clusterization * clusterization;

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
     * @brief Routine that performs the loop needed to coordinate
     * root finding. It has to be called to do the hard work.
     */
    void (*mpsolve_ptr) (mps_context *status);

    /**
     * @brief This is the polynomial that is currently being solved in MPSolve.
     */
    mps_polynomial * active_poly;

    /**
     * @brief Pointer to the secular equation used in computations.
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

    /**
     * @brief Auxiliary memory used in regeneation to avoid thread-safeness
     * issues. 
     */
    mpc_t * bmpc;

    /**
     * @brief True if Jacobi-style iterations must be used in the secular
     * algorithm.
     */
    mps_boolean jacobi_iterations;

    /**
     * @brief Char to be intersted after the with statement in the output piped to gnuplot.
     */
    const char * gnuplot_format;

    /* DEBUG SECTIONS */

    unsigned long int regeneration_time;
    unsigned long int mp_iteration_time;
    unsigned long int dpe_iteration_time;
    unsigned long int fp_iteration_time;

  };                 /* End of typedef struct { ... */

#endif /* #ifdef _MPS_PRIVATE */

  /* Allocator, deallocator, constructors.. */
  mps_context * mps_context_new (void);
  void mps_context_init (mps_context * s);
  void mps_context_free (mps_context * s);

  /* Accessor functions (setters) */
  int mps_context_set_poly_d (mps_context * s, cplx_t * coeff,
                             long unsigned int n);
  void mps_context_set_input_poly (mps_context * s, mps_polynomial * p);
  int mps_context_set_poly_i (mps_context * s, int *coeff, long unsigned int n);
  void mps_context_select_algorithm (mps_context * s, mps_algorithm algorithm);
  void mps_context_set_degree (mps_context * s, int n);

#ifdef _MPS_PRIVATE
  void mps_context_allocate_poly_inplace (mps_context * s, int n);
#endif

  /* Accessor functions */
  long int mps_context_get_data_prec_max (mps_context * s);
  int mps_context_get_degree (mps_context * s);
  int mps_context_get_roots_d (mps_context * s, cplx_t ** roots, double **radius);
  int mps_context_get_roots_m (mps_context * s, mpc_t ** roots, rdpe_t ** radius);
  int mps_context_get_zero_roots (mps_context * s);
  mps_root_status mps_context_get_root_status (mps_context * ctx, int i);
  mps_boolean mps_context_get_over_max (mps_context * s);
  mps_polynomial * mps_context_get_active_poly (mps_context * ctx);

  /* I/O options and flags */
  void mps_context_set_input_prec (mps_context * s, long int prec);
  void mps_context_set_output_prec (mps_context * s, long int prec);
  void mps_context_set_output_format (mps_context * s, mps_output_format format);
  void mps_context_set_output_goal (mps_context * s, mps_output_goal goal);
  void mps_context_set_starting_phase (mps_context * s, mps_phase phase);
  void mps_context_set_log_stream (mps_context * s, FILE * logstr);
  void mps_context_set_jacobi_iterations (mps_context * s, mps_boolean jacobi_iterations);

  /* Debugging */
  void mps_context_set_debug_level (mps_context * s, mps_debug_level level);
  void mps_context_add_debug_domain (mps_context * s, mps_debug_level level);
  
  /* Get input and output config pointers */
  mps_input_configuration * mps_context_get_input_config (mps_context * s);
  mps_output_configuration * mps_context_get_output_config (mps_context * s);

  /* Error handling */
  mps_boolean mps_context_has_errors (mps_context * s);
  char * mps_context_error_msg (mps_context * s);


#ifdef __cplusplus
}
#endif

#endif
