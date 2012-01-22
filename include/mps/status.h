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
#include <pthread.h>

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

#ifdef _MPS_PRIVATE
  /**
   * @brief this struct holds the state of the mps computation
   */
  struct mps_status
  {

    mps_boolean resume;         /* to complete                         */
    mps_boolean chkrad;         /* check radii after completion        */

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
     * @brief Status of approximation of a root. 
     */
    mps_root_status    * root_status;

    /**
     * @brief Array containing attributes that have been set on
     * the roots.
     */
    mps_root_attrs     * root_attrs;

    /**
     * @brief Array containing the inclusion status of the root
     * in the target set specified in the field
     * <code>input_config->search_set</code>.
     */
    mps_root_inclusion * root_inclusion;

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
     * @brief Array containing working precisions used for each root.
     */
    long int *rootwp;

    /**
     * @brief Array that whose i-th component is set to <code>true</code> if
     * the i-th root needs more iterations.
     */
    mps_boolean *again;

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

#endif /* #ifdef _MPS_PRIVATE */

  /* Allocator, deallocator, constructors.. */
  mps_status * mps_status_new (void);
  void mps_status_init (mps_status * s);
  void mps_status_free (mps_status * s);

  /* Accessor functions (setters) */
  int mps_status_set_poly_d (mps_status * s, cplx_t * coeff,
			     long unsigned int n);
  void mps_status_set_input_poly (mps_status * s, mps_monomial_poly * p);
  int mps_status_set_poly_i (mps_status * s, int *coeff, long unsigned int n);
  int mps_status_set_poly_u (mps_status * s, int n, mps_fnewton_ptr fnewton,
			     mps_dnewton_ptr dnewton,
			     mps_mnewton_ptr mnewton);
  void mps_status_select_algorithm (mps_status * s, mps_algorithm algorithm);
  void mps_status_set_degree (mps_status * s, int n);

#ifdef _MPS_PRIVATE
  void mps_status_allocate_poly_inplace (mps_status * s, int n);
#endif

  /* Accessor functions */
  long int mps_status_get_data_prec_max (mps_status * s);
  int mps_status_get_degree (mps_status * s);
  int mps_status_get_roots_d (mps_status * s, cplx_t * roots, double *radius);
  int mps_status_get_roots_m (mps_status * s, mpc_t * roots, rdpe_t * radius);
  int mps_status_get_zero_roots (mps_status * s);
  mps_boolean mps_status_get_over_max (mps_status * s);

  /* I/O options and flags */
  void mps_status_set_input_prec (mps_status * s, long int prec);
  void mps_status_set_output_prec (mps_status * s, long int prec);
  void mps_status_set_output_format (mps_status * s, mps_output_format format);
  void mps_status_set_output_goal (mps_status * s, mps_output_goal goal);
  void mps_status_set_starting_phase (mps_status * s, mps_phase phase);

  /* Debugging */
  void mps_status_set_debug_level (mps_status * s, mps_debug_level level);
  void mps_status_add_debug_domain (mps_status * s, mps_debug_level level);
  
  /* Get input and output config pointers */
  mps_input_configuration * mps_status_get_input_config (mps_status * s);
  mps_output_configuration * mps_status_get_output_config (mps_status * s);



#ifdef __cplusplus
}
#endif

#endif
