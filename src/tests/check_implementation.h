#ifndef CHECK_IMPLEMENTATION_H
#define CHECK_IMPLEMENTATION_H

#ifdef __cplusplus
extern "C"
{
#endif

#include <check.h>
#include <mps/interface.h>
#include <mps/secular.h>
#include <gmp.h>
#include <mps/mpc.h>
#include <mps/gmptools.h>
#include <stdlib.h>

  typedef struct
  {
    char *pol_file;
    char *res_file;
    int out_digits;
    mps_phase phase;
    mps_boolean ga;
  } test_pol;

  void starting_setup ();

  void append_slash (char *dest);

  char *get_pol_file (const char *pol_name, const char *type_name);

  char *get_res_file (const char *pol_name, const char *type_name);

  test_pol *test_pol_new (const char *name, const char *type_name,
                          int out_digits, mps_phase phase, mps_boolean ga);

  void test_pol_free (test_pol * pol);



#ifdef __cplusplus
}
#endif

#endif
