#ifndef MPS_TYPES_H_
#define MPS_TYPES_H_

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
  (sizeof(x) == sizeof(long double) ? isnan_ld (x) \
   : sizeof(x) == sizeof(double) ? isnan_d (x) \
   : isnan_f (x))
static inline int isnan_f (float x)
{
  return x != x;
}
static inline int isnan_d (double x)
{
  return x != x;
}
static inline int isnan_ld (long double x)
{
  return x != x;
}
          #endif

#ifndef isinf
          # define isinf(x) \
  (sizeof(x) == sizeof(long double) ? isinf_ld (x) \
   : sizeof(x) == sizeof(double) ? isinf_d (x) \
   : isinf_f (x))
static inline int isinf_f (float x)
{
  return !isnan (x) && isnan (x - x);
}
static inline int isinf_d (double x)
{
  return !isnan (x) && isnan (x - x);
}
static inline int isinf_ld (long double x)
{
  return !isnan (x) && isnan (x - x);
}
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
struct mps_command_line_option;
struct mps_command_line_option_configuration;

/* list.h */
struct mps_list_element;
struct mps_list;

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

/* regeneration-driver.h */
struct mps_regeneration_driver;

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
typedef struct mps_command_line_option mps_command_line_option;
typedef struct mps_command_line_option_configuration mps_command_line_option_configuration;

/* list.h */
typedef struct mps_list_element mps_list_element;
typedef struct mps_list mps_list;

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
typedef enum mps_starting_strategy mps_starting_strategy;

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

/* regeneration-driver.h */
typedef struct mps_regeneration_driver mps_regeneration_driver;

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
enum mps_phase {
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
enum  mps_algorithm {
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
enum mps_option_key {
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
enum mps_structure {
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
  MPS_OUTPUT_FORMAT_VERBOSE,
  MPS_OUTPUT_FORMAT_OPENSCAD
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
   * \f$S = \{ z \: | \: \lvert z \rvert \leq 1 \}\f$
   */
  MPS_SEARCH_SET_UNITARY_DISC,

  /**
   * @brief Complex number out of the unitary disc
   * \f$S = \{ z \: | \: \lvert z \rvert \leq 1 \}\f$
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
enum mps_representation {
  MPS_REPRESENTATION_SECULAR,
  MPS_REPRESENTATION_MONOMIAL,
  MPS_REPRESENTATION_CHEBYSHEV
};

/**
 * @brief Strategy used to select the starting approximations.
 */
enum mps_starting_strategy {
  MPS_STARTING_STRATEGY_DEFAULT,
  MPS_STARTING_STRATEGY_RECURSIVE,
  MPS_STARTING_STRATEGY_FILE
};

#endif /* endif MPS_TYPES_H_ */
