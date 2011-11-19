/*
 * Check for Secsolve
 */

#include <check.h>
#include <mps/core.h>
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
test_secsolve_on_pol (test_pol * pol)
{
  mps_secular_equation *sec;
  mps_status *s = mps_status_new ();
  FILE *input_stream;
  FILE *check_stream;
  mps_boolean passed = true;
  mpc_t root, ctmp;

  /* Output digit to test, default values. Script should normally
   * alter this with options on the command line. */
  int i, j, prec = pol->out_digits * LOG2_10 + 10;
  int ch;

  fprintf (stderr, "Checking %-30s [\033[34;1mchecking\033[0m]", pol->pol_file);

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

  /* Open streams */
  input_stream = fopen (pol->pol_file, "r");
  check_stream = fopen (pol->res_file, "r");

  fail_if (!(input_stream && check_stream),
           "Cannot open one or more input files");

  /* Set secular equation and start in floating point */
  s->input_config->structure = MPS_STRUCTURE_COMPLEX_FP;
  s->input_config->representation = MPS_REPRESENTATION_SECULAR;
  s->input_config->density = MPS_DENSITY_DENSE;

  mps_parse_stream (s, input_stream);

  sec = s->secular_equation;
  s->input_config->starting_phase = pol->phase;
  s->DOLOG = pol->DOLOG;
  if (pol->DOLOG)
    s->debug_level |= MPS_DEBUG_TRACE;
  

  if (!pol->ga)
    mps_status_select_algorithm (s, MPS_ALGORITHM_SECULAR_MPSOLVE);
  else
    mps_status_select_algorithm (s, MPS_ALGORITHM_SECULAR_GA);

  strncpy (s->goal, "aannc", 5);
  s->output_config->prec = (int) ((pol->out_digits + 1) * LOG2_10) + 1;
  s->input_config->prec = 0;

  s->logstr = stderr;
  s->n_threads = 1;

  mps_mpsolve (s);

  /* printf("Computed results are not exact to the required " */
  /* 	 "precision.\n" "\n" " Dumping test configuration: \n" */
  /* 	 "   => Polynomial file: %s;\n" "   => Required digits: %d\n" */
  /* 	 "   => Gemignani's approach: %s;\n" */
  /* 	 "   => Starting phase: %s;\n", pol->pol_file, pol->out_digits, */
  /* 	 mps_boolean_to_string (pol->ga), */
  /* 	 (pol->phase == float_phase) ? "float_phase" : "dpe_phase");; */
  

  /* Test if roots are equal to the roots provided in the check */
  for (i = 0; i < s->n; i++)
    {
      rdpe_t rtmp, min_dist;
      cdpe_t cdtmp;
      
      rdpe_set (min_dist, RDPE_MAX);

      mpc_clear (root);
      mpc_init2 (root, mpc_get_prec (s->mroot[i]));

      while (isspace (ch = getc (check_stream)));
      ungetc (ch, check_stream);
      mpc_inp_str (root, check_stream, 10);

      passed = false;
      for (j = 0; j < s->n; j++)
        {
          mpc_sub (ctmp, root, s->mroot[j]);

	  mpc_get_cdpe (cdtmp, ctmp);
	  cdpe_mod (rtmp, cdtmp);

	  if (rdpe_le (rtmp, min_dist))
	      rdpe_set (min_dist, rtmp);
        }

      if (rdpe_le (min_dist, s->drad[i]))
	passed = true;
      else
	{
	  /* printf("Setting passed to true with root %d\n", i); */
	  /* printf ("s->mroot[%d] = ", i);  */
	  /* mpc_out_str (stdout, 10, 20, s->mroot[i]);  */
	  /* printf("\n");  */

	  /* printf("s->drad[%d] = ", i); */
	  /* rdpe_out_str (stdout, s->drad[i]);  */
	  /* // printf("%e", s->frad[i]); */
	  /* printf("\n");  */
	  
	  /* printf("min_dist[%d] = ", i);  */
	  /* rdpe_out_str (stdout, min_dist);  */
	  /* printf("\n");  */
	}
      
    }

  mpc_clear (root);
  mpc_clear (ctmp);

  mps_status_free (s);

  fclose (input_stream);
  fclose (check_stream);

  if (passed)
    fprintf (stderr, "\rChecking %-30s [\033[32;1m  done  \033[0m]\n", pol->pol_file);
  else
    fprintf (stderr, "\rChecking %-30s [\033[31;1m failed \033[0m]\n", pol->pol_file);


  fail_unless (passed == true,
               "Computed results are not exact to the required "
               "precision.\n" "\n" " Dumping test configuration: \n"
               "   => Polynomial file: %s;\n" "   => Required digits: %d\n"
               "   => Gemignani's approach: %s;\n"
               "   => Starting phase: %s;\n", pol->pol_file, pol->out_digits,
               mps_boolean_to_string (pol->ga),
               (pol->phase == float_phase) ? "float_phase" : "dpe_phase");
}

START_TEST (test_secsolve)
{
  test_secsolve_on_pol (test_polynomials[_i]);
}

END_TEST
START_TEST (test_secsolve_altern)
{
  /* Start with testing floating point without ga */
  test_pol *pol =
    test_pol_new ("test100", "secsolve", 10, float_phase, false);
  test_secsolve_on_pol (pol);

  /* then floating point with ga */
  pol->ga = true;
  test_secsolve_on_pol (pol);

  test_pol_free (pol);
}

END_TEST
START_TEST (test_secsolve_integer)
{
  /* Test integer parsing of secsolve, ga approach */
  test_pol *pol = test_pol_new ("integer", "secsolve", 250, dpe_phase, true);
  test_secsolve_on_pol (pol);
  test_pol_free (pol);
}
END_TEST

/**
 * @brief Simple test the checks if floating point exception
 * arise in this simple secular equation, that is likely to trigger
 * cancellation problems. 
 */
START_TEST (test_secsolve_simple)
{
  test_pol *pol = test_pol_new ("simple", "secsolve", 15, float_phase, true);
  test_secsolve_on_pol (pol);
  pol->ga = false;
  test_secsolve_on_pol (pol);
  test_pol_free (pol);
}
END_TEST

/**
 * @brief Test secsolve on some secular representation
 * of the wilkinson polynonmials. 
 */
START_TEST (test_secsolve_wilkinson)
{
  /* Testing the wilkinson polynomial of degree 20 */
  test_pol *pol = test_pol_new ("wilk20", "secsolve", 11, float_phase, true);
  test_secsolve_on_pol (pol);
  test_pol_free (pol);

  /* Testing the wilkinson polynomial of degree 40 */
  pol = test_pol_new ("wilk40", "secsolve", 11, float_phase, true);
  test_secsolve_on_pol (pol);
  test_pol_free (pol);

  /* Testing the wilkinson polynomial of degree 80 */
  pol = test_pol_new ("wilk80", "secsolve", 11, float_phase, true);
  test_secsolve_on_pol (pol);
  test_pol_free (pol);
}
END_TEST

START_TEST (test_secsolve_nroots)
{
  /* Testing secsolve on some polynomial of the type x^n - 1 */
  test_pol *pol = test_pol_new ("nroots50", "unisolve", 11, float_phase, true);
  test_secsolve_on_pol (pol);
  test_pol_free (pol);  
}
END_TEST

START_TEST (test_secsolve_kam)
{
  /* Testing the kam polynomials */
  test_pol *pol = test_pol_new ("kam1_1", "unisolve", 11, float_phase, true);
  test_secsolve_on_pol (pol);
  test_pol_free (pol);
}
END_TEST

START_TEST (test_secsolve_mand)
{
  /* Testing secsolve on the mandelbrot polynomials */

  /* Mandelbrot classic, degree 63 */
  test_pol *pol = test_pol_new ("mand63", "unisolve", 11, float_phase, true);
  test_secsolve_on_pol (pol);
  test_pol_free (pol);

  /* Mandelbrot classic, degree 127 */
  pol = test_pol_new ("mand127", "unisolve", 11, float_phase, true);
  test_secsolve_on_pol (pol);
  test_pol_free (pol);
  
}
END_TEST

START_TEST (test_secsolve_exp)
{
  /* Testing secsolve on truncated exponential series */
  test_pol * pol = test_pol_new ("exp50", "unisolve", 11, float_phase, true);
  test_secsolve_on_pol (pol);
  test_pol_free (pol);

  pol = test_pol_new ("exp100", "unisolve", 11, float_phase, true);
  test_secsolve_on_pol (pol);
  test_pol_free (pol);
}
END_TEST

START_TEST (test_secsolve_mignotte)
{
  test_pol * pol = test_pol_new ("mig1_100", "unisolve", 11, float_phase, true);
  test_secsolve_on_pol (pol);
  test_pol_free (pol);

  pol = test_pol_new ("mig1_200", "unisolve", 11, float_phase, true);
  test_secsolve_on_pol (pol);
  test_pol_free (pol);
}
END_TEST

/**
 * @brief Create the secsolve test suite
 */
  Suite * secsolve_suite (int standard)
{
  Suite *s = suite_create ("secsolve");

  /* Create a test case for the standard MPSolve case and
   * one for the Gemignani's approach. */
  TCase *tc_mpsolve = tcase_create ("Secsolve");

  /* Add our tests */
  tcase_add_loop_test (tc_mpsolve, test_secsolve, 0, standard);

  /* Case of a_i = (-1)^(i+1) , b_i = i */
  tcase_add_test (tc_mpsolve, test_secsolve_altern);

  /* Integer parsing */
  tcase_add_test (tc_mpsolve, test_secsolve_integer);

  /* Simple secular equation with cancellation problems */
  tcase_add_test (tc_mpsolve, test_secsolve_simple);

  /* Wilkinson polynomials */
  tcase_add_test (tc_mpsolve, test_secsolve_wilkinson);

  /* Roots of the unity */
  tcase_add_test (tc_mpsolve, test_secsolve_nroots);

  /* Kam polynomials */
  tcase_add_test (tc_mpsolve, test_secsolve_kam);

  /* Exponentials */
  tcase_add_test (tc_mpsolve, test_secsolve_exp);

  /* Mandelbrot polynomials */
  tcase_add_test (tc_mpsolve, test_secsolve_mand);

  /* Chebyshev */
  tcase_add_test (tc_mpsolve, test_secsolve_mignotte);

  /* Add test case to the suite */
  suite_add_tcase (s, tc_mpsolve);

  return s;
}

int
main (void)
{
  int number_failed, standard = 0;

  starting_setup ();

  test_polynomials = (test_pol **) malloc (sizeof (test_pol *) * 12);

  /* Tests with rand15. pol */
  /* Standard MPSolvea approach */
  test_polynomials[standard++] =
    test_pol_new ("rand15", "secsolve", 15, float_phase, false);
  test_polynomials[standard++] =
    test_pol_new ("rand15", "secsolve", 600, float_phase, false);
  test_polynomials[standard++] =
    test_pol_new ("rand15", "secsolve", 15, dpe_phase, false);
  test_polynomials[standard++] =
    test_pol_new ("rand15", "secsolve", 600, dpe_phase, false);

  /* Gemignani's approach */
  test_polynomials[standard++] =
    test_pol_new ("rand15", "secsolve", 15, float_phase, true);
  test_polynomials[standard++] =
    test_pol_new ("rand15", "secsolve", 600, float_phase, true);

  /* Tests with rand120.pol */
  test_polynomials[standard++] =
    test_pol_new ("rand120", "secsolve", 15, float_phase, false);
  test_polynomials[standard++] =
    test_pol_new ("rand120", "secsolve", 15, dpe_phase, false);
  test_polynomials[standard++] =
    test_pol_new ("rand120", "secsolve", 15, float_phase, true);

  /* Tests with deg500.pol */
  test_polynomials[standard++] =
    test_pol_new ("deg500", "secsolve", 15, float_phase, false);
  test_polynomials[standard++] =
    test_pol_new ("deg500", "secsolve", 15, dpe_phase, false);
  test_polynomials[standard++] =
    test_pol_new ("deg500", "secsolve", 15, float_phase, true);

  /* Create a new test suite for secsolve and run it */
  Suite *s = secsolve_suite (standard);
  SRunner *sr = srunner_create (s);
  srunner_run_all (sr, CK_VERBOSE);

  /* Get number of failed test and report */
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);

  for (standard--; standard >= 0; standard--)
    test_pol_free (test_polynomials[standard]);
  free (test_polynomials);

  return (number_failed != 0);
}
