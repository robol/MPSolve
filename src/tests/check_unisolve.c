/*
 * Check for Unisolve
 */

#include <check.h>
#include <mps/core.h>
#include <mps/poly.h>
#include <gmp.h>
#include <mps/mpc.h>
#include <mps/gmptools.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <check_implementation.h>

#define test_pol_new_simple(name, digits) test_pol_new((name), "unisolve", (digits), dpe_phase, false)

/**
 * @brief Array with the polynomials to test
 */
test_pol **test_polynomials;

int
test_unisolve_on_pol (test_pol * pol)
{
  mpspoly_t poly;
  mps_status *s = mps_status_new ();
  FILE *input_stream;
  FILE *check_stream;
  mps_boolean passed = true;
  mpc_t root, ctmp;
  mpf_t mroot;
  mpf_t ftmp;
  mpf_t eps;
  int i, j, prec = pol->out_digits * LOG2_10;
  int ch;

  /* Debug starting of this test */
  /*
     if (pol->ga)
     printf("Starting test on polynomial file %s, with %d output digits required, using GA approach and %s arithmetic.\n",
     pol->pol_file, pol->out_digits, (pol->phase == dpe_phase) ? "DPE" : "floating point");
     else
     printf("Starting test on polynomial file %s, with %d output digits required, using MPSolve approach and %s arithmetic.\n",
     pol->pol_file, pol->out_digits, (pol->phase == dpe_phase) ? "DPE" : "floating point");
   */

  mpc_init2 (root, prec);
  mpc_init2 (ctmp, prec);
  mpf_init2 (mroot, prec);
  mpf_init2 (eps, prec);
  mpf_init2 (ftmp, prec);

  mpf_set_2dl (eps, 1.0, -pol->out_digits);

  /* Open streams */
  input_stream = fopen (pol->pol_file, "r");
  check_stream = fopen (pol->res_file, "r");

  if (!(input_stream && check_stream))
    {
      fail ("Cannot open input files");
    }

  mps_set_default_values (s);
  s->prec_out = prec;
  strncpy (s->goal, "aannc", 5);
  mps_read_poly (s, input_stream, poly);

  mps_set_poly (s, poly);
  mps_allocate_data (s);
  mps_mpsolve (s);
  mps_copy_roots (s);

  /* Test if roots are equal to the roots provided in the check */
  for (i = 0; i < s->n; i++)
    {
      while (isspace (ch = getc (check_stream)));
      ungetc (ch, check_stream);
      mpc_inp_str (root, check_stream, 10);

      passed = false;
      for (j = 0; j < s->n; j++)
        {
          mpc_sub (ctmp, root, s->mroot[j]);
          mpc_mod (mroot, ctmp);
          mpc_mod (ftmp, s->mroot[j]);
          mpf_mul_eq (ftmp, eps);
          if (mpf_cmp (mroot, ftmp) <= 0)
            {
              passed = true;
              break;
            }
        }
    }

  mpf_clear (mroot);
  mpf_clear (eps);
  mpc_clear (root);
  mps_status_free (s);

  fclose (input_stream);
  fclose (check_stream);

  if (s->prec_in > pol->out_digits)
    {
      fail_unless (passed == true,
                   "Computed results are not exact to the required "
                   "precision.\n" "\n" " Dumping test configuration: \n"
                   "   => Polynomial file: %s;\n"
                   "   => Required digits: %d\n", pol->pol_file,
                   pol->out_digits);

    }
}


START_TEST (test_unisolve)
{
  test_unisolve_on_pol (test_polynomials[_i]);
}

END_TEST
/**
 * @brief Create the unisolve test suite
 */
  Suite * unisolve_suite (int n)
{
  Suite *s = suite_create ("unisolve");

  TCase *tc_mpsolve = tcase_create ("MPSolve standard approach");

  /* Add our tests */
  tcase_add_loop_test (tc_mpsolve, test_unisolve, 0, n);

  /* Add test case to the suite */
  suite_add_tcase (s, tc_mpsolve);

  return s;
}

int
main (void)
{
  int number_failed, n = 0, i;
  int digits[] = { 15, 100 };

  starting_setup ();

  test_polynomials = (test_pol **) malloc (sizeof (test_pol *) * 2 * 32);

  for (i = 0; i < 2; i++)
    {
      test_polynomials[n++] = test_pol_new_simple ("exp100", digits[i]);
      test_polynomials[n++] = test_pol_new_simple ("exp50", digits[i]);
      test_polynomials[n++] = test_pol_new_simple ("kam1_1", digits[i]);
      test_polynomials[n++] = test_pol_new_simple ("kam1_2", digits[i]);
      test_polynomials[n++] = test_pol_new_simple ("kam1_3", digits[i]);
      test_polynomials[n++] = test_pol_new_simple ("kam2_1", digits[i]);
      test_polynomials[n++] = test_pol_new_simple ("kam2_2", digits[i]);
      test_polynomials[n++] = test_pol_new_simple ("kam2_3", digits[i]);
      test_polynomials[n++] = test_pol_new_simple ("kam3_1", digits[i]);
      test_polynomials[n++] = test_pol_new_simple ("kam3_2", digits[i]);
      test_polynomials[n++] = test_pol_new_simple ("kam3_3", digits[i]);
      test_polynomials[n++] = test_pol_new_simple ("kam4", digits[i]);
      test_polynomials[n++] = test_pol_new_simple ("kir1_10", digits[i]);
      test_polynomials[n++] = test_pol_new_simple ("lar1_200", digits[i]);
      test_polynomials[n++] = test_pol_new_simple ("lar1", digits[i]);
      test_polynomials[n++] = test_pol_new_simple ("lar2", digits[i]);
      test_polynomials[n++] = test_pol_new_simple ("lar3", digits[i]);
      test_polynomials[n++] = test_pol_new_simple ("lsr_24", digits[i]);
      test_polynomials[n++] = test_pol_new_simple ("mand127", digits[i]);
      test_polynomials[n++] = test_pol_new_simple ("mand63", digits[i]);
      test_polynomials[n++] = test_pol_new_simple ("mand63", digits[i]);
      test_polynomials[n++] = test_pol_new_simple ("mand63", digits[i]);
      test_polynomials[n++] = test_pol_new_simple ("mig1_100", digits[i]);
      test_polynomials[n++] = test_pol_new_simple ("mig1_200", digits[i]);
      test_polynomials[n++] = test_pol_new_simple ("nroots50", digits[i]);
      test_polynomials[n++] = test_pol_new_simple ("spiral20", digits[i]);
      test_polynomials[n++] = test_pol_new_simple ("test", digits[i]);
      test_polynomials[n++] = test_pol_new_simple ("trv_m", digits[i]);
      test_polynomials[n++] = test_pol_new_simple ("umand31", digits[i]);
      test_polynomials[n++] = test_pol_new_simple ("wilk20", digits[i]);
      test_polynomials[n++] = test_pol_new_simple ("wilk40", digits[i]);
      test_polynomials[n++] = test_pol_new_simple ("toep1_128", digits[i]);
    }


  /* Create a new test suite for secsolve and run it */
  Suite *s = unisolve_suite (n);
  SRunner *sr = srunner_create (s);
  srunner_run_all (sr, CK_VERBOSE);

  /* Get number of failed test and report */
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);

  /* Cleanup */
  for (n--; n >= 0; n--)
    test_pol_free (test_polynomials[n]);
  free (test_polynomials);

  return (number_failed != 0);
}
