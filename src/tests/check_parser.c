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

START_TEST (inline_simple6)
{
  ALLOCATE_CONTEXT

  mps_monomial_poly * poly = MPS_MONOMIAL_POLY (
    mps_parse_inline_poly_from_string (ctx, "x^2 - x"));

  /* Verify the parsing */
  fail_unless (poly != NULL, "Cannot parse x^2 - x correctly");
  
  fail_unless (mpq_cmp_si (poly->initial_mqp_r[0], 0, 1) == 0,
	       "Coefficient of degree 0 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_i[0], 0, 1) == 0, 
	       "Coefficient of degree 0 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_r[1], -1, 1) == 0, 
	       "Coefficient of degree 1 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_i[1], 0, 1) == 0, 
	       "Coefficient of degree 1 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_r[2], 1, 1) == 0, 
	       "Coefficient of degree 2 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_i[2], 0, 1) == 0, 
	       "Coefficient of degree 2 has been parsed incorrectly");

  mps_context_free (ctx); 
}
END_TEST

START_TEST (inline_simple7)
{
  ALLOCATE_CONTEXT

  mps_monomial_poly * poly = MPS_MONOMIAL_POLY (
    mps_parse_inline_poly_from_string (ctx, "x^2 - 6/7"));

  /* Verify the parsing */
  fail_unless (poly != NULL, "Cannot parse x^2 - x correctly");
  
  fail_unless (mpq_cmp_si (poly->initial_mqp_r[0], -6, 7) == 0,
	       "Coefficient of degree 0 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_i[0], 0, 1) == 0, 
	       "Coefficient of degree 0 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_r[1], 0, 1) == 0, 
	       "Coefficient of degree 1 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_i[1], 0, 1) == 0, 
	       "Coefficient of degree 1 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_r[2], 1, 1) == 0, 
	       "Coefficient of degree 2 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_i[2], 0, 1) == 0, 
	       "Coefficient of degree 2 has been parsed incorrectly");

  mps_context_free (ctx); 
}
END_TEST

START_TEST (inline_simple8)
{
  ALLOCATE_CONTEXT

  mps_monomial_poly * poly = MPS_MONOMIAL_POLY (
    mps_parse_inline_poly_from_string (ctx, "2342/12x^2 - 6/7 +x"));

  /* Verify the parsing */
  fail_unless (poly != NULL, "Cannot parse x^2 - x correctly");
  
  fail_unless (mpq_cmp_si (poly->initial_mqp_r[0], -6, 7) == 0,
	       "Coefficient of degree 0 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_i[0], 0, 1) == 0, 
	       "Coefficient of degree 0 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_r[1], 1, 1) == 0, 
	       "Coefficient of degree 1 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_i[1], 0, 1) == 0, 
	       "Coefficient of degree 1 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_r[2], 2342, 12) == 0, 
	       "Coefficient of degree 2 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_i[2], 0, 1) == 0, 
	       "Coefficient of degree 2 has been parsed incorrectly");

  mps_context_free (ctx); 
}
END_TEST

START_TEST (inline_simple9)
{
  ALLOCATE_CONTEXT

  mps_monomial_poly * poly = MPS_MONOMIAL_POLY (
    mps_parse_inline_poly_from_string (ctx, "-6/7 + x^2"));

  /* Verify the parsing */
  fail_unless (poly != NULL, "Cannot parse x^2 - x correctly");
  
  fail_unless (mpq_cmp_si (poly->initial_mqp_r[0], -6, 7) == 0,
	       "Coefficient of degree 0 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_i[0], 0, 1) == 0, 
	       "Coefficient of degree 0 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_r[1], 0, 1) == 0, 
	       "Coefficient of degree 1 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_i[1], 0, 1) == 0, 
	       "Coefficient of degree 1 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_r[2], 1, 1) == 0, 
	       "Coefficient of degree 2 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_i[2], 0, 1) == 0, 
	       "Coefficient of degree 2 has been parsed incorrectly");

  mps_context_free (ctx); 
}
END_TEST

START_TEST (inline_cancellation1)
{
  ALLOCATE_CONTEXT

  mps_monomial_poly * poly = MPS_MONOMIAL_POLY (
    mps_parse_inline_poly_from_string (ctx, "-6/7 + x^2 + 12/14"));

  /* Verify the parsing */
  fail_unless (poly != NULL, "Cannot parse x^2 - x correctly");
  
  fail_unless (mpq_cmp_si (poly->initial_mqp_r[0], 0, 1) == 0,
	       "Coefficient of degree 0 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_i[0], 0, 1) == 0, 
	       "Coefficient of degree 0 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_r[1], 0, 1) == 0, 
	       "Coefficient of degree 1 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_i[1], 0, 1) == 0, 
	       "Coefficient of degree 1 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_r[2], 1, 1) == 0, 
	       "Coefficient of degree 2 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_i[2], 0, 1) == 0, 
	       "Coefficient of degree 2 has been parsed incorrectly");

  mps_context_free (ctx); 
}
END_TEST

START_TEST (inline_total_cancellation)
{
  ALLOCATE_CONTEXT

  mps_monomial_poly *poly = MPS_MONOMIAL_POLY (
    mps_parse_inline_poly_from_string (ctx, "56/8 - x^4 + x^2 + x^4 - 2x^2 - 28/4 + x^2"));			 

  fail_unless (poly == NULL, "Not raising error on empty polynomial");

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

  printf ("TEST:malformed_input1 :: Error obtained parsing 1.0/6x^4 - 2: %s\n", 
	  mps_context_error_msg (ctx));

  mps_context_free (ctx);	       
}
END_TEST

START_TEST (malformed_input2)
{
  ALLOCATE_CONTEXT
  mps_monomial_poly * poly = MPS_MONOMIAL_POLY (
    mps_parse_inline_poly_from_string (ctx, " 1e4/12.3x^4 - 2"));

  fail_unless ( (poly == NULL) && 
		 mps_context_has_errors (ctx), 
		"Error not raised on invalid polynomial" );

  printf ("TEST:malformed_input1 :: Error obtained parsing 1e4/12.3x^4 - 2: %s\n", 
	  mps_context_error_msg (ctx));

  mps_context_free (ctx);	       
}
END_TEST

START_TEST (wellformed_input1)
{
  int i;

  ALLOCATE_CONTEXT
  mps_monomial_poly * poly = MPS_MONOMIAL_POLY (
    mps_parse_inline_poly_from_string (ctx, " 6/4x^9 - 2e4"));

  if (mps_context_has_errors (ctx))
    printf ("TEST:malformed_input1 :: Error obtained parsing 6/4x^9 - 2e4: %s\n", 
	    mps_context_error_msg (ctx));

  fail_unless ( (poly != NULL) && 
		 ! mps_context_has_errors (ctx), 
		"Error raised on valid polynomial" );

  fail_unless (mpq_cmp_si (poly->initial_mqp_r[0], -20000, 1) == 0,
	       "Coefficient of degree 0 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_i[0], 0, 1) == 0, 
	       "Coefficient of degree 0 has been parsed incorrectly");

  /* Check the middle ones */
  for (i = 1; i <= 8; i++)
    {
      fail_unless (mpq_cmp_si (poly->initial_mqp_r[i], 0, 2) == 0, 
		   "Coefficient of degree %d has been parsed incorrectly", i);
      fail_unless (mpq_cmp_si (poly->initial_mqp_i[i], 0, 1) == 0, 
		   "Coefficient of degree %d has been parsed incorrectly", i);

    }

  fail_unless (mpq_cmp_si (poly->initial_mqp_r[9], 3, 2) == 0, 
	       "Coefficient of degree 9 has been parsed incorrectly");
  fail_unless (mpq_cmp_si (poly->initial_mqp_i[9], 0, 1) == 0, 
	       "Coefficient of degree 9 has been parsed incorrectly");

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
  tcase_add_test (tc_inline, inline_simple6);
  tcase_add_test (tc_inline, inline_simple7);
  tcase_add_test (tc_inline, inline_simple8);
  tcase_add_test (tc_inline, inline_simple9);

  /* Check for correct add and sum of exponents */
  tcase_add_test (tc_inline, inline_cancellation1);
  tcase_add_test (tc_inline, inline_total_cancellation);

  /* Check for some simple extreme cases */
  tcase_add_test (tc_inline, inline_linear1);

  /* Check for correct error raising */
  tcase_add_test (tc_inline, malformed_input1);
  tcase_add_test (tc_inline, malformed_input2);
  tcase_add_test (tc_inline, wellformed_input1);

  suite_add_tcase (s, tc_inline);

  SRunner *sr = srunner_create (s);
  srunner_run_all (sr, CK_NORMAL);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);

  return (number_failed != 0);
}
