#include <mps/mps.h>
#include <check.h>
#include "check_implementation.h"

START_TEST (set_coefficient_d1)
{
  int n = 4;
  mps_context * ctx = mps_context_new ();
  mps_monomial_poly *poly = mps_monomial_poly_new (ctx, n);
  cplx_t output;

  mps_monomial_poly_set_coefficient_d (ctx, poly, 0,
                                       2.0, -.5);

  mps_monomial_poly_get_coefficient_d (ctx, poly, 0, output);

  fail_unless (cplx_Re (output) == 2.0 &&
               cplx_Im (output) == -0.5,
               "Failed to set the coefficients of the polynomial from double input");

  mps_monomial_poly_free (ctx, MPS_POLYNOMIAL (poly));
  mps_context_free (ctx);
}
END_TEST

START_TEST (set_coefficient_s1)
{
  int n = 4;
  mps_context * ctx = mps_context_new ();
  mps_monomial_poly * poly = mps_monomial_poly_new (ctx, n);
  mpq_t real_coeff, imag_coeff;

  mpq_init (real_coeff);
  mpq_init (imag_coeff);

  mps_monomial_poly_set_coefficient_s (ctx, poly, 0,
                                       "2/3", NULL);
  mps_monomial_poly_get_coefficient_q (ctx, poly, 0,
                                       real_coeff, imag_coeff);

  fail_unless (mpq_cmp_si (real_coeff, 2, 3) == 0 &&
               mpq_cmp_si (imag_coeff, 0, 1) == 0,
               "Failed to set coefficients from string");

  mpq_clear (real_coeff);
  mpq_clear (imag_coeff);
  mps_monomial_poly_free (ctx, MPS_POLYNOMIAL (poly));
  mps_context_free (ctx);
}
END_TEST

START_TEST (set_coefficient_s2)
{
  int n = 4;
  mps_context * ctx = mps_context_new ();
  mps_monomial_poly * poly = mps_monomial_poly_new (ctx, n);
  mpq_t real_coeff, imag_coeff;

  mpq_init (real_coeff);
  mpq_init (imag_coeff);

  mps_monomial_poly_set_coefficient_s (ctx, poly, 2,
                                       "7.8e2", "1.6");
  mps_monomial_poly_get_coefficient_q (ctx, poly, 2,
                                       real_coeff, imag_coeff);

  fail_unless (mpq_cmp_si (real_coeff, 780, 1) == 0 &&
               mpq_cmp_si (imag_coeff, 16, 10) == 0,
               "Failed to set coefficients from string");

  mpq_clear (real_coeff);
  mpq_clear (imag_coeff);
  mps_monomial_poly_free (ctx, MPS_POLYNOMIAL (poly));
  mps_context_free (ctx);
}
END_TEST


int
main (void)
{
  int number_failed;

  starting_setup ();

  Suite *s = suite_create ("Monomial Poly");
  TCase *tc_coefficients = tcase_create ("Coefficients manipulation");

  tcase_add_test (tc_coefficients, set_coefficient_d1);
  tcase_add_test (tc_coefficients, set_coefficient_s1);
  tcase_add_test (tc_coefficients, set_coefficient_s2);

  suite_add_tcase (s, tc_coefficients);

  SRunner *sr = srunner_create (s);

  srunner_run_all (sr, CK_NORMAL);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);

  return(number_failed != 0);
}
