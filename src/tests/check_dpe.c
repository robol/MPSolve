#include <check.h>
#include <check_implementation.h>
#include <mps/core.h>

/**
 * @brief Testing RDPE for safe comparison
 */
START_TEST (test_rdpe_comparison)
{
  rdpe_t a, b;
  
  rdpe_set (a, rdpe_one);
  rdpe_set (b, RDPE_MAX);

  fail_if (!(rdpe_lt (a, b)),
	   "1.0 < RDPE_MAX return false");
}
END_TEST

START_TEST (test_rdpe_overflow)
{
  rdpe_t a, b, c;
  
  rdpe_set (a, RDPE_MAX);
  rdpe_set (b, RDPE_MAX);

  rdpe_add (c, a, b);

  fail_if (rdpe_lt (c, a),
	   "RDPE_MAX + RDPE_MAX < RDPE_MAX");
  
}
END_TEST

/**
 * @brief Create the test suite used to test DPE values.
 */
Suite *
dpe_suite (void)
{
  Suite * s = suite_create ("DPE");
  TCase * tc_rdpe = tcase_create ("RDPE");

  tcase_add_test (tc_rdpe, test_rdpe_comparison);
  tcase_add_test (tc_rdpe, test_rdpe_overflow);
  
  suite_add_tcase (s, tc_rdpe);
  return s;
}

int
main (void)
{
  int number_failed;

  starting_setup ();
  
  Suite * s = dpe_suite ();
  SRunner * sr = srunner_create (s);

  srunner_run_all (sr, CK_VERBOSE);

  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);

  return (number_failed != 0);
}
