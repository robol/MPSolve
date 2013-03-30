#ifndef CHECK_IMPLEMENTATION_H
#define CHECK_IMPLEMENTATION_H

/**
 * @file
 * @brief This file contains the implementation of the checks and some 
 * commodify function to create them.
 *
 * Its main usefulness is to autolocate polynomial files in the tests
 * directory based on the <code>srcdir</code> environment variable
 * that is set by autotools when checking. 
 */

#ifdef __cplusplus
extern "C"
{
#endif

#include <check.h>
#include <mps/interface.h>
#include <mps/secular-equation.h>
#include <gmp.h>
#include <mps/mpc.h>
#include <mps/gmptools.h>
#include <stdlib.h>

  /**
   * @brief Test polynomials to be passed to the function <code>test_*_on_pol()</code>
   */
  typedef struct
  {
    char *pol_file;
    char *res_file;
    int out_digits;
    mps_phase phase;
    mps_boolean ga;
    mps_boolean DOLOG;
  } test_pol;

  void starting_setup (void);

  void append_slash (char *dest);

  const char * get_pol_name_from_path (const char * pol_path);

  char *get_pol_file (const char *pol_name, const char *type_name);

  char *get_res_file (const char *pol_name, const char *type_name);

  test_pol *test_pol_new (const char *name, const char *type_name,
                          int out_digits, mps_phase phase, mps_boolean ga);

  void starting_test_message (const char * pol_file);

  void failed_test_message (const char * pol_file);

  void success_test_message (const char * pol_file);

  void error_test_message (const char *  pol_file, const char * message);

  void test_pol_free (test_pol * pol);



#ifdef __cplusplus
}
#endif

#endif
