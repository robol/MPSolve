#include <mps/mps.h>
#include <check.h>
#include "check_implementation.h"

#define ALLOCATE_CONTEXT \
  mps_context * ctx = mps_context_new();	\
  mps_context_add_debug_domain (ctx, MPS_DEBUG_IO); \
  mps_context_set_log_stream (ctx, stdout);
  

START_TEST (inline_simple1)
{
  ALLOCATE_CONTEXT

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
  ALLOCATE_CONTEXT
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
  ALLOCATE_CONTEXT
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
  ALLOCATE_CONTEXT
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
  ALLOCATE_CONTEXT
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

START_TEST (inline_linear1)
{
  ALLOCATE_CONTEXT
  mps_monomial_poly * poly = MPS_MONOMIAL_POLY (
    mps_parse_inline_poly_from_string (ctx, "x-6"));

  /* Verify the parsing */
  fail_unless (poly != NULL, "Cannot parse x-6 correctly");
  
  fail_unless (mpq_cmp_si (poly->initial_mqp_r[0], -6, 1) == 0,
	       "Coefficient of degree 0 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_i[0], 0, 1) == 0, 
	       "Coefficient of degree 0 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_r[1], 1, 1) == 0, 
	       "Coefficient of degree 1 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_i[1], 0, 1) == 0, 
	       "Coefficient of degree 1 has been parsed incorrectly");

  mps_context_free (ctx); 
}
END_TEST

START_TEST (malformed_input1)
{
  ALLOCATE_CONTEXT
  mps_monomial_poly * poly = MPS_MONOMIAL_POLY (
    mps_parse_inline_poly_from_string (ctx, " 1.0/6x^4 - 2"));

  fail_unless ( (poly == NULL) && 
		 mps_context_has_errors (ctx), 
		"Error not raised on invalid polynomial" );

  MPS_DEBUG (ctx, "Error obtained parsing 1.0/6x^4 - 2: %s", 
	     mps_context_error_msg (ctx));

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

  /* Check for some simple extreme cases */
  tcase_add_test (tc_inline, inline_linear1);

  /* Check for correct error raising */
  tcase_add_test (tc_inline, malformed_input1);

  suite_add_tcase (s, tc_inline);

  SRunner *sr = srunner_create (s);
  srunner_run_all (sr, CK_NORMAL);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);

  return (number_failed != 0);
}
