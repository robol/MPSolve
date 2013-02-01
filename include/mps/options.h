#ifndef MPS_OPTIONS_H_
#define MPS_OPTIONS_H_

/**
 * @file
 * @brief Implementation of option parsing for MPSolve.
 */

#include <mps/mps.h>

#ifdef __cplusplus
extern "C" {
#endif

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
  const static short int mps_rational_structures[] = { 0, 1, 0, 0, 0, 1, 0, 0, 0 };
  const static short int mps_integer_structures[] = { 1, 0, 0, 0, 1, 0, 0, 0, 0 };
  const static short int mps_fp_structures[] = { 0, 0, 1, 0, 0, 0, 1, 0, 0 };
  const static short int mps_real_structures[] = { 1, 1, 1, 1, 0, 0, 0, 0, 0 };
  const static short int mps_complex_structures[] = { 0, 0, 0, 0, 1, 1, 1, 1, 0 };
  const static short int mps_bigfloat_structures[] = { 0, 0, 0, 1, 0, 0, 0, 1, 0 };

  /* STRUCTURE related macros */
#define MPS_STRUCTURE_IS_RATIONAL(x) (mps_rational_structures[(x)])
#define MPS_STRUCTURE_IS_INTEGER(x)  (mps_integer_structures[(x)])
#define MPS_STRUCTURE_IS_FP(x)       (mps_fp_structures[(x)])
#define MPS_STRUCTURE_IS_REAL(x)     (mps_real_structures[(x)])
#define MPS_STRUCTURE_IS_COMPLEX(x)  (mps_complex_structures[(x)])


  const static short int mps_user_representations[]   = { 0, 0, 1 };
  const static short int mps_sparse_representations[] = { 0, 1, 0 };
  const static short int mps_dense_representations[]  = { 1, 0, 0 };

#define MPS_DENSITY_IS_SPARSE(x)   (mps_sparse_representations[(x)])
#define MPS_DENSITY_IS_DENSE(x)    (mps_dense_representations[(x)])

#ifdef _MPS_PRIVATE
  /**
   * @brief Configuration for an input stream; this struct
   * contains the information on how the input stream should
   * be parsed.
   */
  struct mps_input_configuration
  {
    /**
     * @brief Selet the starting phase for the computation.
     *
     * Should be <code>float_phase</code> in the majority of
     * cases, and <code>dpe_phase</code> if the computation
     * is not manageable with the usual IEEE1354 limits.
     */
    mps_phase starting_phase;
  };
#endif /* #ifdef _MPS_PRIVATE */


  /* Properties of the root */
#define MPS_OUTPUT_PROPERTY_NONE      (0x00     )
#define MPS_OUTPUT_PROPERTY_REAL      (0x01     )
#define MPS_OUTPUT_PROPERTY_IMAGINARY (0x01 << 1)

#ifdef _MPS_PRIVATE
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
     * @brief True if the mulitplicity check is enabled in
     * MPSolve.
     */
    mps_boolean multiplicity;

    /**
     * @brief The set in which the roots must be searched. 
     */
    mps_search_set search_set;

    /**
     * @brief These flags are used to determined which properties
     * of the roots must be determined by MPSolve. 
     *
     * Possible values are:
     * -# MPS_OUTPUT_PROPERTY_REAL
     * -# MPS_OUTPUT_PROPERTY_IMAGINARY
     */
    char root_properties;

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
#endif /* ifdef _MPS_PRIVATE */


#ifdef __cplusplus
}
#endif

#endif /* Header end */
