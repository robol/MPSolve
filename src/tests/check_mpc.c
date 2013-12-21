#include <mps/mps.h>
#include <check.h>
#include "check_implementation.h"

START_TEST (basics_addition)
{
  mpc_t a, b, c; 
  cplx_t x; 

  mpc_init2 (a, 150 * LOG2_10 + DBL_MANT_DIG);
  mpc_init2 (b, 150 * LOG2_10 + DBL_MANT_DIG);
  mpc_init2 (c, 150 * LOG2_10 + DBL_MANT_DIG);

  mpc_set_d (a, 1.0, 0.0); 
  mpc_set_d (b, 1e-150, 0.0); 

  mpc_add    (c, a, b); 
  mpc_sub_eq (c, a); 

  mpc_get_cplx (x, c); 

  fail_unless ( fabs(cplx_Re (x) - 1e-150) < DBL_EPSILON * 1e-150 * 8.0, 
	       "Loss of precision in multiprecision arithmetic: addition"); 
  
  mpc_clear (a); 
  mpc_clear (b); 
  mpc_clear (c); 
}
END_TEST

int
main (void)
{
  int number_failed;

  starting_setup ();

  Suite *s = suite_create ("Multiprecision arithmetic");
  TCase *tc_basics = tcase_create ("Simple arithmetic");

  // Basic operation
  tcase_add_test (tc_basics, basics_addition); 

  suite_add_tcase (s, tc_basics);

  SRunner *sr = srunner_create (s);
  srunner_run_all (sr, CK_NORMAL);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);

  return (number_failed != 0);
}
