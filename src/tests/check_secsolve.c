/*
 * Check for Secsolve
 */

#include <check.h>
#include <mps/interface.h>
#include <mps/secular.h>
#include <gmp.h>
#include <mps/mpc.h>
#include <mps/gmptools.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <check_implementation.h>

/**
 * @brief Array with the polynomials to test
 */
test_pol **test_polynomials;

/**
 * @brief This function tests the resolution of a polynomial file
 * referenced by <code>pol</code>.
 */
int
test_secsolve_on_pol(test_pol *pol)
{
    mps_secular_equation* sec;
    mps_status* s = mps_status_new();
    FILE* input_stream;
    FILE* check_stream;
    mps_boolean passed = true;
    mpc_t root, ctmp;
    mpf_t mroot;
    mpf_t ftmp;
    mpf_t eps;

    /* Output digit to test, default values. Script should normally
     * alter this with options on the command line. */
    int i, j, prec = pol->out_digits * LOG2_10;
    int ch;

    /* Debug starting of this test */
    if (pol->ga)
      printf("Starting test on polynomial file %s, with %d output digits required, using GA approach and %s arithmetic.\n",
           pol->pol_file, pol->out_digits, (pol->phase == dpe_phase) ? "DPE" : "floating point");
    else
        printf("Starting test on polynomial file %s, with %d output digits required, using MPSolve approach and %s arithmetic.\n",
             pol->pol_file, pol->out_digits, (pol->phase == dpe_phase) ? "DPE" : "floating point");

    mpc_init2(root, prec);
    mpc_init2(ctmp, prec);
    mpf_init2(mroot, prec);
    mpf_init2(eps,   prec);
    mpf_init2(ftmp,  prec);

    mpf_set_2dl(eps, 1.0, -pol->out_digits);

    /* Open streams */
    input_stream = fopen(pol->pol_file, "r");
    check_stream = fopen(pol->res_file, "r");

    fail_if (!(input_stream && check_stream),
             "Cannot open one or more input files");

    /* Some default values */
    mps_set_default_values(s);
    s->prec_out = prec;
    strncpy(s->goal, "aannc", 5);

    /* Set secular equation and start in floating point */
    mps_parse_stream(s, input_stream);
    s->secular_equation->starting_case = pol->phase;

    mps_status_set_degree(s, s->n);

    if (!pol->ga)
      {
        mps_status_select_algorithm(s, MPS_ALGORITHM_SECULAR_MPSOLVE);
        mps_allocate_data(s);
      }
    else
        mps_status_select_algorithm(s, MPS_ALGORITHM_SECULAR_GA);

    /* Solve the polynomial */
    s->goal[0] = 'a';

    mps_mpsolve(s);
    mps_copy_roots(s);

    /* Test if roots are equal to the roots provided in the check */
    for(i = 0; i < s->n; i++)
      {
        while(isspace(ch = getc(check_stream)));
        ungetc(ch, check_stream);
        mpc_inp_str(root, check_stream, 10);

        passed = false;
        for(j = 0; j < s->n; j++)
          {
            mpc_sub(ctmp, root, s->mroot[j]);
            mpc_mod(mroot, ctmp);
            mpc_mod(ftmp, s->mroot[j]);
            mpf_mul_eq(ftmp, eps);
            if (mpf_cmp(mroot, ftmp) <= 0)
              {
                passed = true;
                break;
              }
          }
      }

    mpf_clear(mroot);
    mpf_clear(eps);
    mpc_clear(root);
    mps_status_free(s);

    fclose (input_stream);
    fclose (check_stream);

    fail_unless(passed == true, "Computed results are not exact to the required "
                "precision, that is of %d digits.\n"
                "Test configuration: pol: %s, ga: %d, starting_phase: %s",
                pol->out_digits, pol->pol_file, pol->ga,
                (pol->phase == float_phase) ? "float_phase" : "dpe_phase");
}

START_TEST (test_secsolve)
{
  test_secsolve_on_pol(test_polynomials[_i]);
}
END_TEST

START_TEST (test_secsolve_altern)
{
  /* Start with testing floating point without ga */
  test_pol *pol = test_pol_new("test1000", "secsolve", 20, float_phase, false);
  test_secsolve_on_pol(pol);

  /* then floating point with ga */
  pol->ga = true;
  test_secsolve_on_pol(pol);

  test_pol_free (pol);
}
END_TEST

/**
 * @brief Create the secsolve test suite
 */
Suite*
secsolve_suite (int standard)
{
  Suite *s = suite_create("secsolve");

  /* Create a test case for the standard MPSolve case and
   * one for the Gemignani's approach. */
  TCase *tc_mpsolve = tcase_create("Secsolve");

  /* Add our tests */
  tcase_add_loop_test(tc_mpsolve, test_secsolve, 0, standard);

  /* Case of a_i = (-1)^(i+1) , b_i = i */
  tcase_add_test (tc_mpsolve, test_secsolve_altern);

  /* Add test case to the suite */
  suite_add_tcase(s, tc_mpsolve);

  return s;
}

int
main (void)
{
  int number_failed, standard = 0;

  starting_setup();

  test_polynomials = (test_pol**) malloc(sizeof(test_pol*) * 9);

  /* Tests with rand15. pol */
  /* Standard MPSolvea approach */
  test_polynomials[standard++] = test_pol_new("rand15", "secsolve", 15, float_phase, false);
  test_polynomials[standard++] = test_pol_new("rand15", "secsolve", 200, float_phase, false);
  test_polynomials[standard++] = test_pol_new("rand15", "secsolve", 15, dpe_phase, false);
  test_polynomials[standard++] = test_pol_new("rand15", "secsolve", 200, dpe_phase, false);

  /* Gemignani's approach */
  test_polynomials[standard++] = test_pol_new("rand15", "secsolve", 15, float_phase, true);
  test_polynomials[standard++] = test_pol_new("rand15", "secsolve", 400, float_phase, true);

  /* Tests with rand120.pol */
  test_polynomials[standard++] = test_pol_new("rand120", "secsolve", 15, float_phase, false);
  test_polynomials[standard++] = test_pol_new("rand120", "secsolve", 15, dpe_phase, false);
  test_polynomials[standard++] = test_pol_new("rand120", "secsolve", 15, float_phase, true);

  /* Create a new test suite for secsolve and run it */
  Suite *s = secsolve_suite(standard);
  SRunner *sr = srunner_create(s);
  srunner_run_all(sr, CK_VERBOSE);

  /* Get number of failed test and report */
  number_failed = srunner_ntests_failed(sr);
  srunner_free(sr);

  for(standard--; standard >= 0; standard--)
    test_pol_free(test_polynomials[standard]);
  free(test_polynomials);

  return (number_failed != 0);
}
