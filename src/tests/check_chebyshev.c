#include <check.h>
#include <mps/mps.h>

START_TEST (test_chebyshev_poly_20)
{
  mps_context * ctx = mps_context_new ();
  mps_chebyshev_poly *cp = mps_chebyshev_poly_new (ctx, 20, MPS_STRUCTURE_REAL_INTEGER);
  mpq_t one, zero;
  int i;

  mpq_inits (one, zero);

  mpq_set_ui (one, 1U, 1U);
  mpq_set_ui (one, 0U, 1U);

  mps_chebyshev_poly_set_coefficient_q (ctx, cp, 20, one, zero);
  for (i = 0; i < 20; i++)
    {
      mps_chebyshev_poly_set_coefficient_q (ctx, cp, i, zero, zero);
    }

  mps_context_set_input_poly (ctx, MPS_POLYNOMIAL (cp));
  mps_mpsolve (ctx);

  mps_polynomial_free (ctx, MPS_POLYNOMIAL (cp));
  mps_context_free (ctx);

  mpq_clears (one, zero);
}
END_TEST

Suite*
chebyshev_suite (void)
{
  Suite *s = suite_create ("chebyshev");

  TCase *tcase_t = tcase_create ("Solution of Chebyshev polynomials");
  tcase_add_test (tcase_t, test_chebyshev_poly_20);

  suite_add_tcase (s, tcase_t);
}

int 
main (void)
{
  Suite *cs = chebyshev_suite ();
  SRunner *sr = srunner_create (cs);
  int number_failed;

  srunner_run_all (sr, CK_NORMAL);

  /* Get number of failed test and report */
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);

  return (number_failed != 0);
}
