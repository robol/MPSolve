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

  const static short int mps_secular_representations[]  = { 1, 0 };
  const static short int mps_monomial_representations[] = { 0, 1 };

#define MPS_INPUT_CONFIG_IS_SECULAR(x)  (mps_secular_representations[(x->representation)])
#define MPS_INPUT_CONFIG_IS_MONOMIAL(x) (mps_monomial_representations[(x->representation)])


  const static short int mps_user_representations[]   = { 0, 0, 1 };
  const static short int mps_sparse_representations[] = { 0, 1, 0 };
  const static short int mps_dense_representations[]  = { 1, 0, 0 };

#define MPS_INPUT_CONFIG_IS_USER(x)     (mps_user_representations[(x->density)])
#define MPS_INPUT_CONFIG_IS_SPARSE(x)   (mps_sparse_representations[(x->density)])
#define MPS_INPUT_CONFIG_IS_DENSE(x)    (mps_dense_representations[(x->density)])


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
     * @brief Condition to be reached to return the computed
     * approximations.
     */
    mps_output_goal goal;

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
