/***********************************************************
 **       Multiprecision Polynomial Solver (MPSolve)       **
 **                 Version 2.2, May 2001                  **
 **                                                        **
 **                      Written by                        **
 **       Dario Andrea Bini and Giuseppe Fiorentino        **
 **       (bini@dm.unipi.it)  (fiorent@dm.unipi.it)        **
 **                                                        **
 ** (C) 2001, Dipartimento di Matematica, FRISCO LTR 21024 **
 ***********************************************************/

/**
 * @file
 *
 * This file is the header for the libmps library. Including
 * this file is needed to access all the MPSolve routines by
 * MPSolve internals.
 *
 * @brief Header file for libmps
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef MPS_CORE_H_
#define MPS_CORE_H_

#ifdef __cplusplus
#define __MPS_NOT_DEFINE_BOOL
#endif

#ifdef __MPS_MATLAB_MODE
#define __MPS_NOT_DEFINE_BOOL
#endif

/* String type used for some hacks */
typedef const char * mps_string;

/* Boolean type used in MPSolve */
#ifndef __MPS_NOT_DEFINE_BOOL
  typedef enum
  { false = 0, true = 1 } mps_boolean;
#else
  /* Small workaround to make matlab module work; there is,
   * int matlab headers, already a false keyword defined, so
   * reusing it here make compilation fail. */
  typedef bool mps_boolean;
#endif                          /* mps_boolean */

#define mps_boolean_to_string(x) ((x) == true) ? "true" : "false"

/* Debug level */
typedef int mps_debug_level;

/* Handle systems where isnan and isinf are not available */
#include <math.h>
#ifndef isnan
          # define isnan(x) \
              (sizeof (x) == sizeof (long double) ? isnan_ld (x) \
               : sizeof (x) == sizeof (double) ? isnan_d (x) \
               : isnan_f (x))
          static inline int isnan_f  (float       x) { return x != x; }
          static inline int isnan_d  (double      x) { return x != x; }
          static inline int isnan_ld (long double x) { return x != x; }
          #endif

#ifndef isinf
          # define isinf(x) \
              (sizeof (x) == sizeof (long double) ? isinf_ld (x) \
               : sizeof (x) == sizeof (double) ? isinf_d (x) \
               : isinf_f (x))
          static inline int isinf_f  (float       x)
          { return !isnan (x) && isnan (x - x); }
          static inline int isinf_d  (double      x)
          { return !isnan (x) && isnan (x - x); }
          static inline int isinf_ld (long double x)
          { return !isnan (x) && isnan (x - x); }
#endif

#include <mps/mt-types.h>

#ifdef __cplusplus

  /* Forward declarations of the type used in the headers, so they can be
   * resolved indepently by the header inclusion order. */

  /* context.h */
  struct mps_context;

  /* cluster.h */
  struct mps_root;
  struct mps_cluster;
  struct mps_cluster_item;
  struct mps_clusterization;

  /* secular-equation.h */
  struct mps_secular_equation;
  struct mps_secular_iteration_data;

  /* monomial-poly.h */
  struct mps_monomial_poly;

  /* monomial-matrix-poly.h */
  struct mps_monomial_matrix_poly;

  /* polynomial.h */
  struct mps_polynomial;

  /* input-buffer.h */
  struct mps_input_buffer;

  /* approximation.h */
  struct mps_approximation;

  /* options.h */
  struct mps_opt;
  struct mps_input_option;

  struct mps_input_configuration;
  struct mps_output_configuration;

  /* threading.h */
  struct mps_thread_job;
  struct mps_thread_job_queue;
  struct mps_thread_worker_data;
  struct mps_thread;
  struct mps_thread_pool;
  struct mps_thread_pool_queue;
  struct mps_thread_pool_queue_item;

#else

  /* Forward declarations of the type used in the headers, so they can be
   * resolved indepently by the header inclusion order. */

  /* context.h */
  typedef struct mps_context mps_context;

  /* cluster.h */
  typedef struct mps_root mps_root;
  typedef struct mps_cluster mps_cluster;
  typedef struct mps_cluster_item mps_cluster_item;
  typedef struct mps_clusterization mps_clusterization;

  /* secular-equation.h */
  typedef struct mps_secular_equation mps_secular_equation;
  typedef struct mps_secular_iteration_data mps_secular_iteration_data;

  /* monomial-poly.h */
  typedef struct mps_monomial_poly mps_monomial_poly;

  /* monomial-matrix-poly.h */
  typedef struct mps_monomial_matrix_poly mps_monomial_matrix_poly;

  /* polynomial.h */
  typedef struct mps_polynomial mps_polynomial;

  /* input-buffer.h */
  typedef struct mps_input_buffer mps_input_buffer;

  /* approximation.h */
  typedef struct mps_approximation mps_approximation;

  /* options.h */
  typedef struct mps_opt mps_opt;
  typedef struct mps_input_option mps_input_option;

  typedef enum mps_root_status mps_root_status;
  typedef enum mps_root_inclusion mps_root_inclusion;
  typedef enum mps_root_attrs mps_root_attrs;

  typedef enum mps_algorithm mps_algorithm;
  typedef enum mps_operation mps_operation;
  typedef enum mps_option_key mps_option_key;
  typedef enum mps_structure mps_structure;
  typedef enum mps_representation mps_representation;
  typedef enum mps_density mps_density;
  typedef enum mps_output_format mps_output_format;
  typedef enum mps_output_goal mps_output_goal;
  typedef enum mps_search_set mps_search_set;
  typedef enum mps_phase mps_phase;

  typedef struct mps_input_configuration mps_input_configuration;
  typedef struct mps_output_configuration mps_output_configuration;

  /* threading.h */
  typedef struct mps_thread_job mps_thread_job;
  typedef struct mps_thread_job_queue mps_thread_job_queue;
  typedef struct mps_thread_worker_data mps_thread_worker_data;
  typedef struct mps_thread mps_thread;
  typedef struct mps_thread_pool mps_thread_pool;
  typedef struct mps_thread_pool_queue mps_thread_pool_queue;
  typedef struct mps_thread_pool_queue_item mps_thread_pool_queue_item;

#endif

  /**
   * @brief Type representing the computation phase
   * of the algorithm we are in
   * now. It can assume the values:
   * - <code>no_phase</code>;
   * - <code>float_phase</code>;
   * - <code>dpe_phase</code>;
   * - <code>mp_phase</code>;
   */
  enum mps_phase
    {
      no_phase, float_phase, dpe_phase, mp_phase
    };

  static const mps_string mps_phase_string [] = {
    "No phase", "Float phase", "DPE phase", "MP phase"
  };
#define MPS_PHASE_TO_STRING(phase) (mps_phase_string[phase])

  /**
   * @brief Used to label different operation inside the various
   * algorithms.
   */
  enum mps_operation {
    MPS_OPERATION_CLUSTER_ANALYSIS,
    MPS_OPERATION_ABERTH_FP_ITERATIONS,
    MPS_OPERATION_ABERTH_DPE_ITERATIONS,
    MPS_OPERATION_ABERTH_MP_ITERATIONS,
    MPS_OPERATION_REGENERATION,
    MPS_OPERATION_STARTING_POINTS_FP,
    MPS_OPERATION_STARTING_POINTS_DPE,
    MPS_OPERATION_STARTING_POINTS_MP,
    MPS_OPERATION_SHIFT,
    MPS_OPERATION_REFINEMENT
  };
  static const mps_string mps_operation_string [] = {
    "Cluster Analysis", "Aberth floating point iterations", "Aberth DPE iterations",
    "Aberth multiprecision iterations", "Regeneration", "Starting point computation in floating point",
    "Starting point computatino in DPE", "Starting point computation in multiprecision",
    "Shift of the polynomial", "Refinement of the approximation"
  };
#define MPS_OPERATION_TO_STRING(operation) (mps_operation_string[operation])

  /**
   * @brief Status of approximation of the root.
   */
  enum mps_root_status {
    MPS_ROOT_STATUS_NEW_CLUSTERED,
    MPS_ROOT_STATUS_CLUSTERED,
    MPS_ROOT_STATUS_ISOLATED,
    MPS_ROOT_STATUS_APPROXIMATED,
    MPS_ROOT_STATUS_APPROXIMATED_IN_CLUSTER,
    MPS_ROOT_STATUS_NOT_FLOAT,
    MPS_ROOT_STATUS_NOT_DPE,
    MPS_ROOT_STATUS_MULTIPLE
  };

  /* Macros to check root status */
  static const mps_boolean mps_table_of_approximated_roots [] = { false, false, false, true, true, false, false, false };
  static const mps_boolean mps_table_of_computed_roots [] = { false, false, true, true, true, false, false, false };
  static const mps_boolean mps_table_of_improvable_roots [] = { false, false, true, true, false, false, false, false };
#define MPS_ROOT_STATUS_IS_APPROXIMATED(status) (mps_table_of_approximated_roots[status]) 
#define MPS_ROOT_STATUS_IS_COMPUTED(status)     (mps_table_of_computed_roots[status]) 
#define MPS_ROOT_STATUS_IS_IMPROVABLE(status)   (mps_table_of_improvable_roots[status]) 

  /* Cast of root_status to string */
  static const mps_string mps_root_status_string[] = {
    "Clustered (pinned)",
    "Clustered",
    "Isolated",
    "Approximated",
    "Approximated in a cluster",
    "Not representable as floating point",
    "Not representable as DPE",
    "Multiple root"
  };
#define MPS_ROOT_STATUS_TO_STRING(status) (mps_root_status_string[status])

  /**
   * @brief Attributes that can be attached to a root and
   * are mostly aimed to detect reality or not. 
   */
  enum mps_root_attrs {
    MPS_ROOT_ATTRS_NONE,
    MPS_ROOT_ATTRS_REAL,
    MPS_ROOT_ATTRS_NOT_REAL,
    MPS_ROOT_ATTRS_IMAG,
    MPS_ROOT_ATTRS_NOT_IMAG,
    MPS_ROOT_ATTRS_NOT_REAL_AND_IMAG
  };

  /* Cast of root_attrs to string */
  static const mps_string mps_root_attrs_string [] = {
    "None",
    "Real",
    "Not real",
    "Imaginary",
    "Not imaginary",
    "Not Real nor imaginary"
  };
#define MPS_ROOT_ATTRS_TO_STRING(attrs) (mps_root_attrs_string[attrs])

  /**
   * @brief Status of inclusion of the root in the target
   * set. 
   */
  enum mps_root_inclusion {
    MPS_ROOT_INCLUSION_UNKNOWN,
    MPS_ROOT_INCLUSION_IN,
    MPS_ROOT_INCLUSION_OUT
  };

  /* Cast of mps_root_inclusion to string */
  static const mps_string mps_root_inclusion_string [] = {
    "Unknown",
    "In",
    "Out",
  };
#define MPS_ROOT_INCLUSION_TO_STRING(inclusion) (mps_root_inclusion_string[inclusion])

  /**
   * @brief Algorithm used to find the solution of the polynomial,
   * or of the secular equation.
   */
  enum  mps_algorithm
    {
      /**
       * @brief Standard MPsolve approach
       */
      MPS_ALGORITHM_STANDARD_MPSOLVE,

      /**
       * @brief Gemignani's approach applied to secular equations.
       */
      MPS_ALGORITHM_SECULAR_GA
    };

  /**
   * @brief Key for options parsed from the input source file.
   * Key that don't need values could exists (these are boolean flags,
   * actually).
   */
  enum mps_option_key
    {
      /* Flag for UNDEFINED Options */
      MPS_FLAG_UNDEFINED,

      /* Key without values associated */
      MPS_FLAG_INTEGER,
      MPS_FLAG_REAL,
      MPS_FLAG_COMPLEX,
      MPS_FLAG_RATIONAL,
      MPS_FLAG_FP,

      MPS_FLAG_SECULAR,
      MPS_FLAG_MONOMIAL,

      MPS_FLAG_DENSE,
      MPS_FLAG_SPARSE,

      /* Key with a value */
      MPS_KEY_DEGREE,
      MPS_KEY_PRECISION,

      /* Key introduced in MPSolve 3.1 */
      MPS_FLAG_CHEBYSHEV
    };

  /**
   * @brief Definition of various algebraic structure that
   * MPSolve can use as input.
   *
   * Precisely, integer, rational and floating point, either real or
   * complex, can be treated in input.
   */
  enum mps_structure
    {
      MPS_STRUCTURE_REAL_INTEGER,
      MPS_STRUCTURE_REAL_RATIONAL,
      MPS_STRUCTURE_REAL_FP,
      MPS_STRUCTURE_REAL_BIGFLOAT,
      MPS_STRUCTURE_COMPLEX_INTEGER,
      MPS_STRUCTURE_COMPLEX_RATIONAL,
      MPS_STRUCTURE_COMPLEX_FP,
      MPS_STRUCTURE_COMPLEX_BIGFLOAT,
      MPS_STRUCTURE_UNKNOWN
    };

  /**
   * @brief Density of the polynomial, or 
   * MPS_DENSITY_USER if density doesn't make sense
   * since user routines are provided to compute
   * the newton fraction.
   */
  enum mps_density {
    MPS_DENSITY_DENSE,
    MPS_DENSITY_SPARSE,
    MPS_DENSITY_USER,
  };

  /**
   * @brief Desired output format for the roots.
   */
  enum mps_output_format {
    MPS_OUTPUT_FORMAT_COMPACT,
    MPS_OUTPUT_FORMAT_GNUPLOT,
    MPS_OUTPUT_FORMAT_GNUPLOT_FULL,
    MPS_OUTPUT_FORMAT_BARE,
    MPS_OUTPUT_FORMAT_FULL,
    MPS_OUTPUT_FORMAT_VERBOSE
  };

  /**
   * @brief Goal to reach before returning the result.
   */
  enum mps_output_goal {
    MPS_OUTPUT_GOAL_ISOLATE,
    MPS_OUTPUT_GOAL_APPROXIMATE,
    MPS_OUTPUT_GOAL_COUNT
  };

  /** 
   * @brief Set in which the roots are searched.
   */
  enum mps_search_set {
    /**
     * @brief The whole complex plane.
     */
    MPS_SEARCH_SET_COMPLEX_PLANE,

    /**
     * @brief Complex numbers with a positive real part.
     */
    MPS_SEARCH_SET_POSITIVE_REAL_PART,

    /**
     * @brief Complex numbers with a negative real part.
     */
    MPS_SEARCH_SET_NEGATIVE_REAL_PART,

    /**
     * @brief Complex numbers with a positive imaginary part.
     */
    MPS_SEARCH_SET_POSITIVE_IMAG_PART,

    /**
     * @brief Complex numbers with a negative real part.
     */
    MPS_SEARCH_SET_NEGATIVE_IMAG_PART,

    /**
     * @brief Complex numbers in the unitary disc
     * \f$S = \{ z \: | \: \lvert z \rvert \leq 1 \}$
     */
    MPS_SEARCH_SET_UNITARY_DISC,

    /**
     * @brief Complex number out of the unitary disc
     * \d$S = \{ z \: | \: \lvert z \rvert \leq 1 \}$
     */
    MPS_SEARCH_SET_UNITARY_DISC_COMPL,

    /**
     * @brief Only real roots.
     */
    MPS_SEARCH_SET_REAL,

    /**
     * @brief Only pure imaginary roots.
     */
    MPS_SEARCH_SET_IMAG,

    /**
     * @brief Custom set specified by the user.
     */
    MPS_SEARCH_SET_CUSTOM
  };

  /**
   * @brief Representation chosen for the polynomial
   */
  enum mps_representation
    {
      MPS_REPRESENTATION_SECULAR,
      MPS_REPRESENTATION_MONOMIAL,
      MPS_REPRESENTATION_CHEBYSHEV
    };

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* Local include files that should not be included directly */
#include <mps/options.h>
#include <mps/cluster.h>
#include <mps/tools.h>
#include <mps/mt.h>
#include <mps/gmptools.h>
#include <mps/mpc.h>
#include <mps/polynomial.h>
#include <mps/link.h>
#include <mps/debug.h>
#include <mps/input-buffer.h>
#include <mps/context.h>
#include <mps/monomial-poly.h>
#include <mps/monomial-matrix-poly.h>
#include <mps/secular-equation.h>
#include <mps/chebyshev.h>
#include <mps/approximation.h>

/* Interface should be a subset of core, so what is defined
 * there should be included here. */
#include <mps/interface.h>
#include <mps/threading.h>

/* constants */

#define MPS_ALL_CLUSTERS -1

#ifdef __cplusplus
  extern "C" {
#endif


/* FUNCTIONS */

  /* functions in aberth.c */
#ifdef _MPS_PRIVATE
  void mps_faberth (mps_context * s, mps_approximation * root, cplx_t abcorr);
  void mps_daberth (mps_context * s, mps_approximation * root, cdpe_t abcorr);
  void mps_maberth (mps_context * s, mps_approximation * root, mpc_t abcorr);
  void mps_faberth_s (mps_context * s, mps_approximation * root, mps_cluster * cluster, cplx_t abcorr);
  void mps_faberth_wl (mps_context * s, int j, cplx_t abcorr, pthread_mutex_t * aberth_mutexes);
  void mps_daberth_s (mps_context * s, mps_approximation * root, mps_cluster * cluster, cdpe_t abcorr);
  void mps_daberth_wl (mps_context * s, int j, cdpe_t abcorr, pthread_mutex_t * aberth_mutexes);
  void mps_maberth_s (mps_context * s, mps_approximation * root, mps_cluster * cluster, mpc_t abcorr);
  void mps_maberth_s_wl (mps_context * s, int j, mps_cluster * cluster, mpc_t abcorr,
                         pthread_mutex_t * aberth_mutex);
  void mps_mnewtis (mps_context * s);
#endif

  /* functions in approximation.c */
  mps_approximation * mps_approximation_new (mps_context * s);
  void mps_approximation_free (mps_context * s, mps_approximation * appr);
  mps_approximation * mps_approximation_copy (mps_context * ctx, mps_approximation * original);

  /* functions in cluster.c */
  void mps_cluster_reset (mps_context * s);
  void mps_fcluster (mps_context * s, double * frad, int nf);
  void mps_dcluster (mps_context * s, rdpe_t * drad, int nf);
  void mps_mcluster (mps_context * s, rdpe_t * drad, int nf);
  void mps_debug_cluster_structure (mps_context * s);

  /* functions in convex.c */
  void mps_fconvex (mps_context * s, int n, double a[]);

  /* functions in cluster-analsys.c */
  void mps_cluster_analysis (mps_context * ctx, mps_polynomial * p);

  /* functions in data.c */
  void mps_mp_set_prec (mps_context * s, long int prec);
  void mps_allocate_data (mps_context * s);
  void mps_prepare_data (mps_context * s, long int prec);
  void mps_restore_data (mps_context * s);
  void mps_free_data (mps_context * s);
  long int mps_raise_data (mps_context * s, long int prec);
  void mps_raise_data_raw (mps_context * s, long int prec);

  /* functions in improve.c */
  void mps_improve (mps_context * s);

  /* functions in jacobi-aberth.c */
  int mps_faberth_packet (mps_context * ctx, mps_polynomial * p, mps_boolean just_regenerated);
  int mps_daberth_packet (mps_context * ctx, mps_polynomial * p, mps_boolean just_regenerated);
  int mps_maberth_packet (mps_context * ctx, mps_polynomial * p, mps_boolean just_regenerated);
  
  /* functions in main.c */
  void mps_setup (mps_context * s);
  void mps_check_data (mps_context * s, char *which_case);
  void mps_compute_sep (mps_context * s);
  void mps_standard_mpsolve (mps_context * s);

  /* functions in newton.c */
  void mps_fnewton (mps_context * st, mps_polynomial * p, 
                    mps_approximation * root, cplx_t corr);
  void mps_dnewton (mps_context * st, mps_polynomial * p,
                    mps_approximation * root, cdpe_t corr);
  void mps_mnewton (mps_context * st, mps_polynomial * p, 
                    mps_approximation * root, mpc_t corr);
  void mps_parhorner (mps_context * st, int n, mpc_t x, mpc_t p[],
                      mps_boolean b[], mpc_t s, int n_thread);
  void mps_aparhorner (mps_context * st, int n, rdpe_t x, rdpe_t p[],
                       mps_boolean b[], rdpe_t s, int n_thread);
  int mps_intlog2 (int n);

  /* Functions in general-radius.c */
  void mps_fradii (mps_context * s, mps_polynomial * p, double * fradii);
  void mps_dradii (mps_context * s, mps_polynomial * p, rdpe_t * dradii);
  void mps_mradii (mps_context * s, mps_polynomial * p, rdpe_t * dradii);

  /* Functions in monomial-radius.c */
  void mps_monomial_fradii (mps_context * s, double * fradii);
  void mps_monomial_dradii (mps_context * s, rdpe_t * dradii);
  void mps_monomial_mradii (mps_context * s, rdpe_t * dradii);

  /* Functions in secular-evaluation.c */
  mps_boolean mps_secular_feval (mps_context * s, mps_polynomial * p, cplx_t x, cplx_t value);
  mps_boolean mps_secular_feval_with_error (mps_context * s, mps_polynomial * p, cplx_t x, cplx_t value, double * error);
  mps_boolean mps_secular_deval (mps_context * s, mps_polynomial * p, cdpe_t x, cdpe_t value);
  mps_boolean mps_secular_deval_derivative (mps_context * s, mps_polynomial * p, cdpe_t x, cdpe_t value);
  mps_boolean mps_secular_deval_with_error (mps_context * s, mps_polynomial * p, cdpe_t x, cdpe_t value, rdpe_t error);
  mps_boolean mps_secular_meval (mps_context * s, mps_polynomial * p, mpc_t x, mpc_t value);
  mps_boolean mps_secular_meval_with_error (mps_context * s, mps_polynomial * p, mpc_t x, mpc_t value, rdpe_t error);
  mps_boolean mps_secular_feval_derivative (mps_context * s, mps_polynomial * p, cplx_t x, cplx_t value);
  
  /* Function in getopts.c */
  void mps_parse_opts (mps_context * s, int argc, char *argv[]);
  mps_boolean mps_getopts (mps_opt ** opt, int *argc_ptr, char ***argv_ptr,
                           const char *opt_format);



  /* functions in sort.c */
  void mps_fsort (mps_context * s);
  void mps_dsort (mps_context * s);
  void mps_msort (mps_context * s);

  /* functions in solve.c */
  void mps_update (mps_context * s);
  void mps_fsrad (mps_context * s, mps_cluster * cluster, cplx_t sc, double *sr);
  void mps_dsrad (mps_context * s, mps_cluster * cluster, cdpe_t sc, rdpe_t sr);
  void mps_msrad (mps_context * s, mps_cluster * cluster, mpc_t sc, rdpe_t sr);

  mps_boolean mps_check_stop (mps_context * s);
  void mps_fsolve (mps_context * s, mps_boolean * d_after_f);
  void mps_dsolve (mps_context * s, mps_boolean d_after_f);
  void mps_msolve (mps_context * s);
  void mps_fpolzer (mps_context * s, int *it, mps_boolean * excep);
  void mps_dpolzer (mps_context * s, int *it, mps_boolean * excep);
  void mps_mpolzer (mps_context * s, int *it, mps_boolean * excep);

  /* Functions in modify.c */
  void mps_fmodify (mps_context * s, mps_boolean track_new_cluster);
  void mps_dmodify (mps_context * s, mps_boolean track_new_cluster);
  void mps_mmodify (mps_context * s, mps_boolean track_new_cluster);

  /* Functions in inclusion.c */
  void mps_fupdate_inclusions (mps_context * s);
  void mps_dupdate_inclusions (mps_context * s);
  void mps_mupdate_inclusions (mps_context * s);
  
  /* functions in starting.c */
  double mps_maximize_distance (mps_context * s, double last_sigma,
                                mps_cluster_item * cluster, int n);
  void mps_fstart (mps_context * s, int n, mps_cluster_item * cluster, double clust_rad,
                   double g, rdpe_t eps_out, double fap[]);
  void mps_dstart (mps_context * s, int n, mps_cluster_item * cluster, rdpe_t clust_rad,
                   rdpe_t g, rdpe_t eps_out, rdpe_t dap[]);
  void mps_mstart (mps_context * s, int n, mps_cluster_item * cluster, rdpe_t clust_rad,
                   rdpe_t g, rdpe_t dap[], mpc_t gg);
  void mps_frestart (mps_context * s);
  void mps_drestart (mps_context * s);
  void mps_mrestart (mps_context * s);
  void mps_fshift (mps_context * s, int m, mps_cluster_item * cluster, double clust_rad,
                   cplx_t g, rdpe_t eps);
  void mps_dshift (mps_context * s, int m, mps_cluster_item * cluster, rdpe_t clust_rad,
                   cdpe_t g, rdpe_t eps);
  void mps_mshift (mps_context * s, int m, mps_cluster_item * cluster, rdpe_t clust_rad,
                   mpc_t g);

  /* functions in stio.c */
  void mps_readroots (mps_context * s);
  void mps_countroots (mps_context * s);
  void mps_outroot (mps_context * s, int i, int num);
  void mps_output (mps_context * s);
  void mps_copy_roots (mps_context * s);
  void mps_dump_status (mps_context * s, FILE * outstr);
  void mps_dump (mps_context * s);
  void mps_dump_cluster_structure (mps_context * s, FILE * outstr);
  mps_boolean mps_is_a_tty (FILE * stream);
  void mps_warn (mps_context * st, char *s);
  void mps_error (mps_context * st, const char * format, ...);
  void mps_print_errors (mps_context * s);



  /* functions in test.c */
  mps_boolean mps_inclusion (mps_context * s);

  /* functions in cluster.c */
  void mps_cluster_detach (mps_context * s, mps_cluster * cluster);

  /* functions in touch.c */
  mps_boolean mps_ftouchnwt (mps_context * s, double * frad, int n, int i, int j);
  mps_boolean mps_dtouchnwt (mps_context * s, rdpe_t * drad, int n, int i, int j);
  mps_boolean mps_mtouchnwt (mps_context * s, rdpe_t * drad, int n, int i, int j);
  mps_boolean mps_ftouchreal (mps_context * s, int n, int i);
  mps_boolean mps_dtouchreal (mps_context * s, int n, int i);
  mps_boolean mps_mtouchreal (mps_context * s, int n, int i);
  mps_boolean mps_ftouchimag (mps_context * s, int n, int i);
  mps_boolean mps_dtouchimag (mps_context * s, int n, int i);
  mps_boolean mps_mtouchimag (mps_context * s, int n, int i);
  mps_boolean mps_ftouchunit (mps_context * s, int n, int i);
  mps_boolean mps_dtouchunit (mps_context * s, int n, int i);
  mps_boolean mps_mtouchunit (mps_context * s, int n, int i);

  /* functions in user.c */
  void mps_fnewton_usr (mps_context * st, mps_polynomial * poly, mps_approximation * root, cplx_t corr);
  void mps_dnewton_usr (mps_context * st, mps_polynomial * poly, mps_approximation * root, cdpe_t corr);
  void mps_mnewton_usr (mps_context * st, mps_polynomial * poly, mps_approximation * root, mpc_t corr);
  mps_boolean mps_feval_usr (mps_context * ctx, mps_polynomial * p, cplx_t x, cplx_t value, double * error);
  mps_boolean mps_deval_usr (mps_context * ctx, mps_polynomial * p, cdpe_t x, cdpe_t value, rdpe_t error);
  mps_boolean mps_meval_usr (mps_context * ctx, mps_polynomial * p, mpc_t x, mpc_t value, rdpe_t error);

  /* functions in general-starting.c */
  void mps_general_fstart (mps_context * ctx, mps_polynomial * p);
  void mps_general_dstart (mps_context * ctx, mps_polynomial * p);
  void mps_general_mstart (mps_context * ctx, mps_polynomial * p);

  /* Routines of Input/Output in input-output.c */
  void mps_skip_comments (FILE * input_stream);
  void mps_raise_parsing_error (mps_context * s, mps_input_buffer * buffer, 
                         const char * token, 
                         const char * message, ...);
  mps_input_option mps_parse_option_line (mps_context * s, char *line, size_t length);

  mps_polynomial * mps_parse_stream (mps_context * s, FILE * input_stream);
  mps_polynomial * mps_parse_file   (mps_context * s, const char * path);

  mps_polynomial * mps_monomial_poly_read_from_stream_v2 (mps_context * s, mps_input_buffer * buffer);
  mps_monomial_poly * mps_monomial_poly_read_from_stream (mps_context * s, mps_input_buffer * buffer, 
    mps_structure structure, mps_density density);
  mps_chebyshev_poly * mps_chebyshev_poly_read_from_stream (mps_context * ctx, mps_input_buffer * buffer,
    mps_structure structure, mps_density density);
  mps_secular_equation * mps_secular_equation_read_from_stream (mps_context * ctx, mps_input_buffer * buffer,
    mps_structure structure, mps_density density);

  /* Functions in horner.c */
  void mps_fhorner (mps_context * s, mps_monomial_poly * p, cplx_t x, cplx_t value);
  void mps_fhorner_with_error (mps_context * s, mps_monomial_poly * p, cplx_t x, cplx_t value, double * relative_error);
  void mps_dhorner (mps_context * s, mps_monomial_poly * p, cdpe_t x, cdpe_t value);
  void mps_dhorner_with_error (mps_context * s, mps_monomial_poly * p, cdpe_t x, cdpe_t value, rdpe_t relative_error);
  void mps_mhorner (mps_context * s, mps_monomial_poly * p, mpc_t x, mpc_t value);
  void mps_mhorner_with_error (mps_context * s, mps_monomial_poly * p, mpc_t x, mpc_t value, rdpe_t relative_error, long int wp);
  void mps_mhorner_with_error2 (mps_context * s, mps_monomial_poly * p, mpc_t x, mpc_t value, rdpe_t relative_error, long int wp);

#ifdef __cplusplus
  }
#endif

#endif                          /* ndef MPSCORE_H */
