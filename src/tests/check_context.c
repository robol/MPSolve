#include <mps/mps.h>
#include <check.h>
#include "check_implementation.h"

START_TEST (basics_allocate_context)
{
	mps_context * ctx = mps_context_new ();
	mps_context_free (ctx);
}
END_TEST

START_TEST (basics_context_reuse_without_free)
{
	mps_context * ctx = mps_context_new ();

	/* Start with a simple polynomial of degree 2 */
	mps_monomial_poly *poly = mps_monomial_poly_new (ctx, 2);

	mps_monomial_poly_set_coefficient_d (ctx, poly, 0, -1, 0.0);
	mps_monomial_poly_set_coefficient_d (ctx, poly, 2, 1, 0.0);

	mps_context_set_input_poly (ctx, MPS_POLYNOMIAL (poly));
	mps_mpsolve (ctx);

	/* Then try to load a polynomial with a higher degree */
	poly = mps_monomial_poly_new (ctx, 10);
	mps_monomial_poly_set_coefficient_d (ctx, poly, 0, -1, 0.0);
	mps_monomial_poly_set_coefficient_d (ctx, poly, 10, 1, 0.0);

	mps_context_set_input_poly (ctx, MPS_POLYNOMIAL (poly));
	mps_mpsolve (ctx);
	mps_monomial_poly_free (ctx, MPS_POLYNOMIAL (poly));
}
END_TEST

START_TEST (basics_context_reuse_expand)
{
	mps_context * ctx = mps_context_new ();

	/* Start with a simple polynomial of degree 2 */
	mps_monomial_poly *poly = mps_monomial_poly_new (ctx, 2);

	mps_monomial_poly_set_coefficient_d (ctx, poly, 0, -1, 0.0);
	mps_monomial_poly_set_coefficient_d (ctx, poly, 2, 1, 0.0);

	mps_context_set_input_poly (ctx, MPS_POLYNOMIAL (poly));
	mps_mpsolve (ctx);
	mps_monomial_poly_free (ctx, MPS_POLYNOMIAL (poly));

	/* Then try to load a polynomial with a higher degree */
	poly = mps_monomial_poly_new (ctx, 10);
	mps_monomial_poly_set_coefficient_d (ctx, poly, 0, -1, 0.0);
	mps_monomial_poly_set_coefficient_d (ctx, poly, 10, 1, 0.0);

	mps_context_set_input_poly (ctx, MPS_POLYNOMIAL (poly));
	mps_mpsolve (ctx);
	mps_monomial_poly_free (ctx, MPS_POLYNOMIAL (poly));
}
END_TEST

START_TEST (basics_context_reuse_shrink)
{
	mps_context * ctx = mps_context_new ();

	/* Start with a simple polynomial of degree 10 */
	mps_monomial_poly *poly = mps_monomial_poly_new (ctx, 10);

	mps_monomial_poly_set_coefficient_d (ctx, poly, 0, -1, 0.0);
	mps_monomial_poly_set_coefficient_d (ctx, poly, 10, 1, 0.0);

	mps_context_set_input_poly (ctx, MPS_POLYNOMIAL (poly));
	mps_mpsolve (ctx);
	mps_monomial_poly_free (ctx, MPS_POLYNOMIAL (poly));

	/* Then try to load a polynomial with a smaller degree */
	poly = mps_monomial_poly_new (ctx, 2);
	mps_monomial_poly_set_coefficient_d (ctx, poly, 0, -1, 0.0);
	mps_monomial_poly_set_coefficient_d (ctx, poly, 2, 1, 0.0);

	mps_context_set_input_poly (ctx, MPS_POLYNOMIAL (poly));
	mps_mpsolve (ctx);
	mps_monomial_poly_free (ctx, MPS_POLYNOMIAL (poly));
}
END_TEST

int
main (void)
{
  int number_failed;

  starting_setup ();

  Suite *s = suite_create ("Matrices");
  TCase *tc_basics = tcase_create ("Basic operations");

  // Basic operation
  tcase_add_test (tc_basics, basics_allocate_context);
  tcase_add_test (tc_basics, basics_context_reuse_without_free);
  tcase_add_test (tc_basics, basics_context_reuse_expand);
  tcase_add_test (tc_basics, basics_context_reuse_shrink);

  suite_add_tcase (s, tc_basics);

  SRunner *sr = srunner_create (s);
  srunner_run_all (sr, CK_NORMAL);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);

  return (number_failed != 0);
}
