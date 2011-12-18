#ifndef MPS_OPTIONS_H_
#define MPS_OPTIONS_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <mps/mps.h>

  /**
   * @brief Struct holding the options passed on the command
   * line.
   *
   * Typical usage is something like this:
   * @code
   * mps_opt* opt;
   * while (opt = mps_getopt(&argc, &argv, format))
   *   {
   *     switch(opt->optchar)
   *       {
   *         case 'a':
   *           [...]
   *           break;
   *         case 'n':
   *           n = atoi(opt->opvalue);
   *           break;
   *       }
   *     free(opt);
   *   }
   * @endcode
   *
   * @see mps_getopts()
   */
  struct mps_opt
  {
    char optchar;
    char *optvalue;
  };

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
       * @brief Standard MPSolve approach applied to
       * secular equations.
       */
      MPS_ALGORITHM_SECULAR_MPSOLVE,

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
      MPS_KEY_DEGREE
    };

  /**
   * @brief This struct holds a key and the value associated
   * with it. It's used for options that require a value associated.
   *
   * For example the option Degree needs a numeric value associated to
   * it. The values are always saved as the string that are read.
   */
  struct mps_input_option
  {
    /**
     * @brief Key associated with the option.
     */
    mps_option_key flag;

    /**
     * @brief Value of the flag, or NULL if no value
     * is provided.
     */
    char *value;
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
      MPS_STRUCTURE_COMPLEX_INTEGER,
      MPS_STRUCTURE_COMPLEX_RATIONAL,
      MPS_STRUCTURE_COMPLEX_FP
    };

  /* STRUCTURES truth tables */
  const static short int mps_rational_structures[] = { 0, 1, 0, 0, 1, 0 };
  const static short int mps_integer_structures[] = { 1, 0, 0, 1, 0, 0 };
  const static short int mps_fp_structures[] = { 0, 0, 1, 0, 0, 1 };
  const static short int mps_real_structures[] = { 1, 1, 1, 0, 0, 0 };
  const static short int mps_complex_structures[] = { 0, 0, 0, 1, 1, 1 };

  /* STRUCTURE related macros */
#define MPS_INPUT_CONFIG_IS_RATIONAL(x) (mps_rational_structures[(x)->structure])
#define MPS_INPUT_CONFIG_IS_INTEGER(x)  (mps_integer_structures[(x)->structure])
#define MPS_INPUT_CONFIG_IS_FP(x)       (mps_fp_structures[(x)->structure])
#define MPS_INPUT_CONFIG_IS_REAL(x)     (mps_real_structures[(x)->structure])
#define MPS_INPUT_CONFIG_IS_COMPLEX(x)  (mps_complex_structures[(x)->structure])

  /**
   * @brief Representation chosen for the polynomial
   */
  enum mps_representation
    {
      MPS_REPRESENTATION_SECULAR,
      MPS_REPRESENTATION_MONOMIAL,
    };

  const static short int mps_secular_representations[]  = { 1, 0 };
  const static short int mps_monomial_representations[] = { 0, 1 };

#define MPS_INPUT_CONFIG_IS_SECULAR(x)  (mps_secular_representations[(x->representation)])
#define MPS_INPUT_CONFIG_IS_MONOMIAL(x) (mps_monomial_representations[(x->representation)])

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

  const static short int mps_user_representations[]   = { 0, 0, 1 };
  const static short int mps_sparse_representations[] = { 0, 1, 0 };
  const static short int mps_dense_representations[]  = { 1, 0, 0 };

#define MPS_INPUT_CONFIG_IS_USER(x)     (mps_user_representations[(x->density)])
#define MPS_INPUT_CONFIG_IS_SPARSE(x)   (mps_sparse_representations[(x->density)])
#define MPS_INPUT_CONFIG_IS_DENSE(x)    (mps_dense_representations[(x->density)])

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
   * @brief Configuration for an input stream; this struct
   * contains the information on how the input stream should
   * be parsed.
   */
  struct mps_input_configuration
  {
    /**
     * @brief Structure of the input data. For every
     * structure a particular format is defined.
     */
    mps_structure structure;

    /**
     * @brief Algorithm to be used, that will determine
     * the coefficients that will be passed to MPSolve.
     */
    mps_representation representation;

    /**
     * @brief Density of the coefficients, or MPS_DENSITY_USER
     * if the coefficients (or the newton fraction) is provided
     * via a user routine
     */
    mps_density density;

    /**
     * @brief Digits of guaranteed input precision.
     */
    long int prec;

    /**
     * @brief Selet the starting phase for the computation.
     *
     * Should be <code>float_phase</code> in the majority of
     * cases, and <code>dpe_phase</code> if the computation
     * is not manageable with the usual IEEE1354 limits.
     */
    mps_phase starting_phase;
  };

  /**
   * @brief Configuration for the output.
   *
   * This struct holds the information on what has to be
   * computed by MPSolve, such as the desired output precision
   * and the search set for the roots, etc.
   */
  struct mps_output_configuration
  {
    /**
     * @brief Digits of required output precision
     */
    long int prec;

    /**
     * @brief Desired output format.
     *
     * Can be one of:
     * @code
     *  MPS_OUTPUT_FORMAT_BARE
     *  MPS_OUTPUT_FORMAT_GNUPLOT
     *  MPS_OUTPUT_FORMAT_GNUPLOT_FULL
     *  MPS_OUTPUT_FORMAT_COMPACT
     *  MPS_OUTPUT_FORMAT_VERBOSE
     *  MPS_OUTPUT_FORMAT_FULL
     * @endcode
     */
    mps_output_format format;

  };


#ifdef __cplusplus
}
#endif

#endif /* Header end */
