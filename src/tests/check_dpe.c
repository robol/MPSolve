#include <check.h>
#include <check_implementation.h>
#include <mps/mps.h>
#include <limits.h>

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

START_TEST (test_rdpe_sum_overflow)
{
  rdpe_t a, b, c;
  
  rdpe_set (a, RDPE_MAX);
  rdpe_set (b, RDPE_MAX);

  rdpe_add (c, a, b);

  fail_if (rdpe_lt (c, a),
           "Overflow: RDPE_MAX + RDPE_MAX < RDPE_MAX");
  
}
END_TEST

START_TEST (test_rdpe_mul_overflow)
{
  rdpe_t a, b, c;
  
  rdpe_set_2dl (a, 0.5, LONG_MAX - 5);
  rdpe_mul_eq_d (a, 10000000.0f);

  fail_if (!rdpe_eq (a, RDPE_MAX),
           "Overflow in rdpe_mul_eq_d: log2 (2^%ld * 10000000) = %ld", LONG_MAX - 6, rdpe_Esp (a));

  rdpe_set_2dl (a, 0.5, LONG_MAX - 5);
  rdpe_set_dl (b, 0.5, 15);
  rdpe_mul_eq (a, b);

  fail_if (!rdpe_eq (a, RDPE_MAX),
           "Overflow in rdpe_mul_eq: log2 (2^%ld * 2^15) = %ld", LONG_MAX - 6, rdpe_Esp (a));

  rdpe_set_2dl (a, 0.5, LONG_MAX - 5);
  rdpe_set_dl (b, 0.5, 15);
  rdpe_mul (c, b, a);

  fail_if (!rdpe_eq (c, RDPE_MAX),
           "Overflow in rdpe_mul: log2 (2^%ld * 2^15) = %ld", LONG_MAX - 6, rdpe_Esp (a));

  rdpe_set_2dl (a, 0.5, LONG_MAX - 5);
  rdpe_mul_d (c, a, 1000000000.0f);
  
  fail_if (!rdpe_eq (c, RDPE_MAX),
           "Overflow in rdpe_mul_d: log2 (2^%ld * 1000000000.0) = %ld", LONG_MAX - 6, rdpe_Esp (a));

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
  tcase_add_test (tc_rdpe, test_rdpe_sum_overflow);
  tcase_add_test (tc_rdpe, test_rdpe_mul_overflow);
  
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
