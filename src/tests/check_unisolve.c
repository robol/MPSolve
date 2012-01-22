/*
 * Check for Unisolve
 */

#include <check.h>
#include <mps/mps.h>
#include <gmp.h>
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
  mps_status *s = mps_status_new ();
  FILE *input_stream;
  FILE *check_stream;
  mps_boolean passed = true;
  mpc_t root, ctmp;
  int i, j, prec = s->data_prec_max.value * 2;
  int zero_roots = 0;
  int ch;

  fprintf (stderr, "Checking \033[1m%-30s\033[0m [\033[34;1mchecking\033[0m]", pol->pol_file + 11);

  /* Debug starting of this test */
  /*
     if (pol->ga)
     printf("Starting test on polynomial file %s, with %d output digits required, using GA approach and %s arithmetic.\n",
     pol->pol_file, pol->out_digits, (pol->phase == dpe_phase) ? "DPE" : "floating point");
     else
     printf("Starting test on polynomial file %s, with %d output digits required, using MPSolve approach and %s arithmetic.\n",
     pol->pol_file, pol->out_digits, (pol->phase == dpe_phase) ? "DPE" : "floating point");
   */

  mpc_init2 (root, 2 * prec);
  mpc_init2 (ctmp, 2 * prec);

  /* Open streams */
  input_stream = fopen (pol->pol_file, "r");
  check_stream = fopen (pol->res_file, "r");

  if (!(input_stream && check_stream))
    {
      fail ("Cannot open input files");
    }

  s->output_config->goal = MPS_OUTPUT_GOAL_ISOLATE;

  mps_parse_stream (s, input_stream);
  s->output_config->prec = prec;

  mps_mpsolve (s);

  /* Test if roots are equal to the roots provided in the check */
  passed = true;
  for (i = 0; i < s->deg; i++)
    {
      rdpe_t rtmp;
      cdpe_t cdtmp;
      rdpe_t min_dist;
      int found_root = 0;

      while (isspace (ch = getc (check_stream)));
      ungetc (ch, check_stream);
      mpc_inp_str (root, check_stream, 10);

      if (mpc_eq_zero (root))
	{
	  zero_roots++;

	  /* We need to read it another time. This seems a bug in
	   * mpc_inp_str, but I don't get why is necessary. */
	  mpc_inp_str (root, check_stream, 10);
	  continue;
	}

      mpc_sub (ctmp, root, s->mroot[0]);
      mpc_get_cdpe (cdtmp, ctmp);
      cdpe_mod (rtmp, cdtmp);
      rdpe_set (min_dist, rtmp);

      for (j = 1; j < s->n; j++)
        {
          mpc_sub (ctmp, root, s->mroot[j]);
	  mpc_get_cdpe (cdtmp, ctmp);
	  cdpe_mod (rtmp, cdtmp);

	  if (rdpe_le (rtmp, min_dist))
	    {
	      rdpe_set (min_dist, rtmp);
	      found_root = j;
	    }
        }

      if (!rdpe_le (min_dist, s->drad[found_root]) && !s->over_max)
	{
	  passed = false;
	  if (getenv ("MPS_VERBOSE_TEST"))
	    {
	      printf("Setting passed to false with root %d\n", found_root); 
	      printf ("s->mroot[%d] = ", found_root);  
	      mpc_out_str (stdout, 10, prec, s->mroot[found_root]);  
	      printf("\n");  
	      
	      printf("s->drad[%d] = ", found_root); 
	      rdpe_out_str (stdout, s->drad[found_root]); printf ("\n");
	      
	      printf("min_dist[%d] = ", found_root);  
	      rdpe_out_str (stdout, min_dist);  
	      printf("\n");  
	    }
	}
    }

  mpc_clear (root);
  mpc_clear (ctmp);

  if (zero_roots != s->zero_roots)
    {
      passed = false;
    }

  fclose (input_stream);
  fclose (check_stream);

  mps_status_free (s);

  if (passed)
    fprintf (stderr, "\rChecking \033[1m%-30s\033[0m [\033[32;1m  done  \033[0m]\n", pol->pol_file + 11);
  else
    fprintf (stderr, "\rChecking \033[1m%-30s\033[0m [\033[31;1m failed \033[0m]\n", pol->pol_file + 11);
  
  if (getenv ("MPS_VERBOSE_TEST"))
    fail_unless (passed == true,
		 "Computed results are not exact to the required "
		 "precision.\n" "\n" " Dumping test configuration: \n"
		 "   => Polynomial file: %s;\n"
		 "   => Required digits: %d\n", pol->pol_file,
		 pol->out_digits);
  else
    fail_unless (passed == true,
		 "Computed results are not exact to the required precision");
  

  return passed;
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
      test_polynomials[n++] = test_pol_new_simple ("wilk40", digits[i] );
      test_polynomials[n++] = test_pol_new_simple ("toep1_128", digits[i]);
    }

  /* Create a new test suite for secsolve and run it */
  Suite *s = unisolve_suite (n);
  SRunner *sr = srunner_create (s);
  srunner_run_all (sr, CK_NORMAL);

  /* Get number of failed test and report */
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);

  /* Cleanup */
  for (n--; n >= 0; n--)
    test_pol_free (test_polynomials[n]);
  free (test_polynomials);

  return (number_failed != 0);
}
