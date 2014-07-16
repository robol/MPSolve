#include <mps/mps.h>
#include <check.h>
#include "check_implementation.h"

START_TEST (test_strip)
{
  mps_context * ctx = mps_context_new ();

  const char * test1 = " !as230 ";
  char * test1_stripped = mps_utils_strip_string (ctx, test1);

  printf ("TEST_STRIP: mps_utils_strip_string(\"%s\") = \"%s\"\n", test1, test1_stripped);
  fail_unless (strcmp (test1_stripped, "!as230") == 0,
               "mps_utils_strip_string (\" !as230 \") != \"!as230\"");
  free (test1_stripped);

  const char * test2 = "as88 ";
  char * test2_stripped = mps_utils_strip_string (ctx, test2);
  printf ("TEST_STRIP: mps_utils_strip_string(\"%s\") = \"%s\"\n", test2, test2_stripped);
  fail_unless (strcmp (test2_stripped, "as88") == 0,
               "mps_utils_strip_string (\"as88 \") != \"as88\"");
  free (test2_stripped);

  const char * test3 = "    as88 ";
  char * test3_stripped = mps_utils_strip_string (ctx, test3);
  printf ("TEST_STRIP: mps_utils_strip_string(\"%s\") = \"%s\"\n", test3, test3_stripped);
  fail_unless (strcmp (test3_stripped, "as88") == 0,
               "mps_utils_strip_string (\"   as88 \") != \"as88\"");
  free (test3_stripped);

  const char * test4 = "Test";
  char * test4_stripped = mps_utils_strip_string (ctx, test4);
  printf ("TEST_STRIP: mps_utils_strip_string(\"%s\") = \"%s\"\n", test4, test4_stripped);
  fail_unless (strcmp (test4_stripped, "Test") == 0,
               "mps_utils_strip_string (\"Test\") != \"Test\"");
  free (test4_stripped);

  mps_context_free (ctx);
}
END_TEST

START_TEST (test_equivalent_rational_conversion)
{
  mps_context * ctx = mps_context_new ();

  mps_context_add_debug_domain (ctx, MPS_DEBUG_IO);

  mpq_t t1, t2;
  mpq_init (t1);
  mpq_init (t2);

  const char * test1 = "1.23";
  char * test1_eq = mps_utils_build_equivalent_rational_string (ctx, test1);
  printf ("TEST_EQUIVALENT_RATIONAL_CONVERSION: Converted 1.23 to %s\n", test1_eq);

  fail_unless ((test1_eq != NULL) && (strcmp (test1_eq, "123/100") == 0),
               "Failed to convert \"1.23\" to \"123/100\"");

  free (test1_eq);
  mpq_clear (t1);
  mpq_clear (t2);

  mps_context_free (ctx);
}
END_TEST

START_TEST (test_equivalent_rational_conversion2)
{
  mpq_t t1, t2;

  mpq_init (t1);
  mpq_init (t2);

  mps_context * ctx = mps_context_new ();

  mpq_set_si (t1, 123, 1);
  const char * test2 = "1.23e2";
  char * test2_eq = mps_utils_build_equivalent_rational_string (ctx, test2);
  printf ("TEST_EQUIVALENT_RATIONAL_CONVERSION: Converted 1.23e2 to %s\n", test2_eq);

  mpq_set_str (t2, test2_eq, 10);
  mpq_canonicalize (t2);

  fail_unless ((test2_eq != NULL) && mpq_equal (t1, t2),
               "Failed to convert \"1.23e12\" to \"123/1\"");

  free (test2_eq);

  mpq_clear (t1);
  mpq_clear (t2);
  mps_context_free (ctx);
}
END_TEST

START_TEST (test_equivalent_rational_conversion3)
{
  mps_context * ctx = mps_context_new ();
  mpq_t t1, t2;

  mpq_init (t1);
  mpq_init (t2);

  mpq_set_si (t1, 123, 1000);
  const char * test3 = "1.23e-1";
  char * test3_eq = mps_utils_build_equivalent_rational_string (ctx, test3);
  printf ("TEST_EQUIVALENT_RATIONAL_CONVERSION: Converted 1.23e-1 to %s\n", test3_eq);

  mpq_set_str (t2, test3_eq, 10);

  fail_unless ((test3_eq != NULL) && mpq_equal (t1, t2),
               "Failed to convert \"1.23e-1\" to \"123/1000\"");

  free (test3_eq);

  mpq_clear (t1);
  mpq_clear (t2);

  mps_context_free (ctx);
}
END_TEST

START_TEST (test_equivalent_rational_conversion4)
{
  mps_context * ctx = mps_context_new ();
  mpq_t t1, t2;

  mpq_init (t1);
  mpq_init (t2);

  mpq_set_si (t1, -8, 1);
  const char * test4 = "-8";
  char * test4_eq = mps_utils_build_equivalent_rational_string (ctx, test4);
  printf ("TEST_EQUIVALENT_RATIONAL_CONVERSION: Converted -8 to %s\n", test4_eq);

  mpq_set_str (t2, test4_eq, 10);

  fail_unless ((test4_eq != NULL) && mpq_equal (t1, t2),
               "Failed to convert \"-8\" to \"-8/1\"");

  free (test4_eq);

  mpq_clear (t1);
  mpq_clear (t2);

  mps_context_free (ctx);
}
END_TEST

int
main (void)
{
  int number_failed;

  starting_setup ();

  Suite *s = suite_create ("Utils");
  TCase *tc_strings = tcase_create ("String manipulation");

  tcase_add_test (tc_strings, test_strip);
  tcase_add_test (tc_strings, test_equivalent_rational_conversion);
  tcase_add_test (tc_strings, test_equivalent_rational_conversion2);
  tcase_add_test (tc_strings, test_equivalent_rational_conversion3);
  tcase_add_test (tc_strings, test_equivalent_rational_conversion4);

  suite_add_tcase (s, tc_strings);

  SRunner *sr = srunner_create (s);

  srunner_run_all (sr, CK_NORMAL);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);

  return(number_failed != 0);
}
