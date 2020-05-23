/*
 * This file is part of MPSolve 3.1.9
 *
 * Copyright (C) 2001-2020, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Leonardo Robol <leonardo.robol@unipi.it>
 */

#ifndef MPS_OPTIONS_H_
#define MPS_OPTIONS_H_

/**
 * @file
 * @brief Implementation of option parsing for MPSolve.
 */

#include <mps/mps.h>

MPS_BEGIN_DECLS

/**
 * @brief This struct holds a configuration for a command line option.
 * This is a step towards a more flexible implementation of the option parser,
 * compared to the traditional getopts() call.
 */
struct mps_command_line_option {
  /**
   * @brief This is the character that is recognized as starting the
   * option specification on the command line.
   *
   * This value may be '\0' if only the long format is provided for this
   * option.
   */
  char format;

  /**
   * @brief This value is true if an argument may be specified for the option.
   *
   * The argument may or may not be mandatory, according to the mandatory field
   * of this struct.
   */
  mps_boolean argument;

  /**
   * @brief If this value is true then the argument for the option is mandatory.
   * Note that this value should be true only if argument is true.
   */
  mps_boolean mandatory;

  /**
   * @brief An optional long format for the option, or NULL if no long format is
   * specified.
   */
  char * long_format;
};

/**
 * @brief Configuration for a command line parser.
 *
 * This struct essentialy holds a list of mps_command_line_option structs
 * that describe the options that should be parsed at command line.
 */
struct mps_command_line_option_configuration {
  /**
   * @brief A list of mps_command_option instances that have been provided
   * for this parser configuration.
   */
  mps_list * command_options;
};

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
struct mps_opt {
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
struct mps_input_option {
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


const static short int mps_user_representations[] = { 0, 0, 1 };
const static short int mps_sparse_representations[] = { 0, 1, 0 };
const static short int mps_dense_representations[] = { 1, 0, 0 };

#define MPS_DENSITY_IS_SPARSE(x)   (mps_sparse_representations[(x)])
#define MPS_DENSITY_IS_DENSE(x)    (mps_dense_representations[(x)])

/**
 * @brief Configuration for an input stream; this struct
 * contains the information on how the input stream should
 * be parsed.
 */
struct mps_input_configuration {
  /**
   * @brief Selet the starting phase for the computation.
   *
   * Should be <code>float_phase</code> in the majority of
   * cases, and <code>dpe_phase</code> if the computation
   * is not manageable with the usual IEEE1354 limits.
   */
  mps_phase starting_phase;
};

/* Properties of the root */
#define MPS_OUTPUT_PROPERTY_NONE      (0x00)
#define MPS_OUTPUT_PROPERTY_REAL      (0x01)
#define MPS_OUTPUT_PROPERTY_IMAGINARY (0x01 << 1)

/**
 * @brief Configuration for the output.
 *
 * This struct holds the information on what has to be
 * computed by MPSolve, such as the desired output precision
 * and the search set for the roots, etc.
 */
struct mps_output_configuration {
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

/* Function in getopts.c */
void mps_parse_opts (mps_context * s, int argc, char *argv[]);
mps_boolean mps_getopts (mps_opt ** opt, int *argc_ptr, char ***argv_ptr,
                         const char *opt_format);


mps_command_line_option_configuration * mps_command_line_option_configuration_new (void);

MPS_END_DECLS

#endif /* Header end */
