#include <mps/mps.h>
#include <check.h>
#include "check_implementation.h"


START_TEST (inline_simple1)
{
  mps_context * ctx = mps_context_new ();
  mps_monomial_poly * poly = MPS_MONOMIAL_POLY (
    mps_parse_inline_poly_from_string (ctx, "x^4 - 1"));

  /* Verify the parsing */
  fail_unless (poly != NULL, "Cannot parse x^4 - 1 correctly");
  
  fail_unless (mpq_cmp_si (poly->initial_mqp_r[0], -1, 1) == 0, 
	       "Coefficient of degree 0 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_i[0], 0, 1) == 0, 
	       "Coefficient of degree 0 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_r[1], 0, 1) == 0, 
	       "Coefficient of degree 1 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_i[1], 0, 1) == 0, 
	       "Coefficient of degree 1 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_r[2], 0, 1) == 0, 
	       "Coefficient of degree 2 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_i[2], 0, 1) == 0, 
	       "Coefficient of degree 2 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_r[3], 0, 1) == 0, 
	       "Coefficient of degree 3 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_i[3], 0, 1) == 0, 
	       "Coefficient of degree 3 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_r[4], 1, 1) == 0, 
	       "Coefficient of degree 4 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_i[4], 0, 1) == 0, 
	       "Coefficient of degree 4 has been parsed incorrectly");

  mps_context_free (ctx); 
}
END_TEST

START_TEST (inline_simple2)
{
  mps_context * ctx = mps_context_new ();
  mps_monomial_poly * poly = MPS_MONOMIAL_POLY (
    mps_parse_inline_poly_from_string (ctx, "100/100x^4 - 14/14"));

  /* Verify the parsing */
  fail_unless (poly != NULL, "Cannot parse 100/100x^4 - 14/14 correctly");
  
  fail_unless (mpq_cmp_si (poly->initial_mqp_r[0], -1, 1) == 0, 
	       "Coefficient of degree 0 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_i[0], 0, 1) == 0, 
	       "Coefficient of degree 0 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_r[1], 0, 1) == 0, 
	       "Coefficient of degree 1 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_i[1], 0, 1) == 0, 
	       "Coefficient of degree 1 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_r[2], 0, 1) == 0, 
	       "Coefficient of degree 2 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_i[2], 0, 1) == 0, 
	       "Coefficient of degree 2 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_r[3], 0, 1) == 0, 
	       "Coefficient of degree 3 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_i[3], 0, 1) == 0, 
	       "Coefficient of degree 3 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_r[4], 1, 1) == 0, 
	       "Coefficient of degree 4 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_i[4], 0, 1) == 0, 
	       "Coefficient of degree 4 has been parsed incorrectly");

  mps_context_free (ctx); 
}
END_TEST

START_TEST (inline_simple3)
{
  mps_context * ctx = mps_context_new ();
  mps_monomial_poly * poly = MPS_MONOMIAL_POLY (
    mps_parse_inline_poly_from_string (ctx, "1.00x^4 - 1.23"));

  /* Verify the parsing */
  fail_unless (poly != NULL, "Cannot parse 100/100x^4 - 1.23 correctly");
  
  fail_unless (mpq_cmp_si (poly->initial_mqp_r[0], -123, 100) == 0, 
	       "Coefficient of degree 0 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_i[0], 0, 1) == 0, 
	       "Coefficient of degree 0 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_r[1], 0, 1) == 0, 
	       "Coefficient of degree 1 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_i[1], 0, 1) == 0, 
	       "Coefficient of degree 1 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_r[2], 0, 1) == 0, 
	       "Coefficient of degree 2 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_i[2], 0, 1) == 0, 
	       "Coefficient of degree 2 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_r[3], 0, 1) == 0, 
	       "Coefficient of degree 3 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_i[3], 0, 1) == 0, 
	       "Coefficient of degree 3 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_r[4], 1, 1) == 0, 
	       "Coefficient of degree 4 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_i[4], 0, 1) == 0, 
	       "Coefficient of degree 4 has been parsed incorrectly");

  mps_context_free (ctx); 
}
END_TEST

START_TEST (inline_simple4)
{
  mps_context * ctx = mps_context_new ();
  mps_monomial_poly * poly = MPS_MONOMIAL_POLY (
    mps_parse_inline_poly_from_string (ctx, "1.00x^4 - 1.23e2"));

  /* Verify the parsing */
  fail_unless (poly != NULL, "Cannot parse 100/100x^4 - 1.23e2 correctly");
  
  fail_unless (mpq_cmp_si (poly->initial_mqp_r[0], -12300, 100) == 0, 
	       "Coefficient of degree 0 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_i[0], 0, 1) == 0, 
	       "Coefficient of degree 0 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_r[1], 0, 1) == 0, 
	       "Coefficient of degree 1 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_i[1], 0, 1) == 0, 
	       "Coefficient of degree 1 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_r[2], 0, 1) == 0, 
	       "Coefficient of degree 2 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_i[2], 0, 1) == 0, 
	       "Coefficient of degree 2 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_r[3], 0, 1) == 0, 
	       "Coefficient of degree 3 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_i[3], 0, 1) == 0, 
	       "Coefficient of degree 3 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_r[4], 1, 1) == 0, 
	       "Coefficient of degree 4 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_i[4], 0, 1) == 0, 
	       "Coefficient of degree 4 has been parsed incorrectly");

  mps_context_free (ctx); 
}
END_TEST

START_TEST (inline_simple5)
{
  mps_context * ctx = mps_context_new ();
  mps_monomial_poly * poly = MPS_MONOMIAL_POLY (
    mps_parse_inline_poly_from_string (ctx, "1.04x^4 - 1.23e-2"));

  /* Verify the parsing */
  fail_unless (poly != NULL, "Cannot parse 100/100x^4 - 1.23e-24 correctly");
  
  fail_unless (mpq_cmp_si (poly->initial_mqp_r[0], -123, 10000) == 0,
	       "Coefficient of degree 0 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_i[0], 0, 1) == 0, 
	       "Coefficient of degree 0 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_r[1], 0, 1) == 0, 
	       "Coefficient of degree 1 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_i[1], 0, 1) == 0, 
	       "Coefficient of degree 1 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_r[2], 0, 1) == 0, 
	       "Coefficient of degree 2 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_i[2], 0, 1) == 0, 
	       "Coefficient of degree 2 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_r[3], 0, 1) == 0, 
	       "Coefficient of degree 3 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_i[3], 0, 1) == 0, 
	       "Coefficient of degree 3 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_r[4], 104, 100) == 0, 
	       "Coefficient of degree 4 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_i[4], 0, 1) == 0, 
	       "Coefficient of degree 4 has been parsed incorrectly");

  mps_context_free (ctx); 
}
END_TEST

int
main (void)
{
  int number_failed;

  starting_setup ();

  Suite *s = suite_create ("Parsers");
  TCase *tc_inline = tcase_create ("Inline parser");

  /* Check inline parsing of polynomials */
  tcase_add_test (tc_inline, inline_simple1); 
  tcase_add_test (tc_inline, inline_simple2); 
  tcase_add_test (tc_inline, inline_simple3);
  tcase_add_test (tc_inline, inline_simple4);
  tcase_add_test (tc_inline, inline_simple5);

  suite_add_tcase (s, tc_inline);

  SRunner *sr = srunner_create (s);
  srunner_run_all (sr, CK_NORMAL);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);

  return (number_failed != 0);
}
