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

int test_unisolve_on_pol_impl (test_pol *, mps_output_goal);

int
test_unisolve_on_pol (test_pol * pol)
{
  return test_unisolve_on_pol_impl (pol, MPS_OUTPUT_GOAL_ISOLATE) &&
    test_unisolve_on_pol_impl (pol, MPS_OUTPUT_GOAL_APPROXIMATE);
}

int
test_unisolve_on_pol_impl (test_pol * pol, mps_output_goal goal)
{
  mpc_t root, ctmp;
  mps_boolean passed = true;
  int i, j, zero_roots = 0;
  char ch;
  rdpe_t eps;
  mps_polynomial * poly = NULL;

  /* Check the roots */
  FILE* result_stream = fopen (pol->res_file, "r"); 
  FILE* input_stream  = fopen (pol->pol_file, "r");

  rdpe_set_2dl (eps, 1.0, - pol->out_digits);

  if (!result_stream) 
    {
      error_test_message ("no results file found", pol->pol_file);
      return EXIT_FAILURE;
    }
  if (!input_stream)
    {
      error_test_message ("no polynomial file found", pol->pol_file);
      return EXIT_FAILURE;
    }

  /* Create a new empty mps_context */
  mps_context * s = mps_context_new ();

  if (getenv ("MPS_VERBOSE_TEST") && strstr (pol->pol_file, getenv ("MPS_VERBOSE_TEST")))
    mps_context_set_debug_level (s, MPS_DEBUG_TRACE);

  /* Load the polynomial that has been given to us */
  poly = mps_parse_stream (s, input_stream);
  mps_context_set_input_poly (s, poly);

  starting_test_message (pol->pol_file);

  mps_context_set_output_goal (s, goal);
  mps_context_set_output_prec (s, pol->out_digits);

  /* Solve it */
  mps_context_select_algorithm (s, MPS_ALGORITHM_STANDARD_MPSOLVE);
  mps_mpsolve (s);
  
  mpc_init2 (root, mps_context_get_data_prec_max (s));
  mpc_init2 (ctmp, mps_context_get_data_prec_max (s));
    
  /* Test if roots are equal to the roots provided in the check */   
  passed = true;

  rdpe_t * drad = rdpe_valloc (mps_context_get_degree (s));
  mpc_t * mroot = mpc_valloc (mps_context_get_degree (s));
  mpc_vinit2 (mroot, mps_context_get_degree (s), 53);

  mps_context_get_roots_m (s, &mroot, &drad);

  for (i = 0; i < mps_context_get_degree (s); i++)   
    {   
      rdpe_t rtmp;   
      cdpe_t cdtmp;   
      rdpe_t min_dist;
      int found_root = 0;
      rdpe_t exp_drad;
      
      while (isspace (ch = getc (result_stream)));   
      ungetc (ch, result_stream);   
      mpc_inp_str (root, result_stream, 10);   

      if (mpc_eq_zero (root))
        {
          zero_roots++;

          /* We need to read it another time. This seems a bug in
           * mpc_inp_str, but I don't get why is necessary. */
          mpc_inp_str (root, result_stream, 10);
          continue;
        }
      
      mpc_sub (ctmp, root, mroot[0]);   
      mpc_get_cdpe (cdtmp, ctmp);   
      cdpe_mod (rtmp, cdtmp);   
      rdpe_set (min_dist, rtmp);   

      if (getenv ("MPS_VERBOSE_TEST") && (strstr (pol->pol_file, getenv ("MPS_VERBOSE_TEST"))))
        {
          printf ("Read root_%d = ", i);
          mpc_out_str_2 (stdout, 10, mps_context_get_data_prec_max (s), mps_context_get_data_prec_max (s),
                         root);
          printf ("\n");
        }
      
      for (j = 1; j < mps_context_get_degree (s); j++)   
        {   
          mpc_sub (ctmp, root, mroot[j]);
          mpc_get_cdpe (cdtmp, ctmp);   
          cdpe_mod (rtmp, cdtmp);   
          
          if (rdpe_le (rtmp, min_dist))
            {
              rdpe_set (min_dist, rtmp);
              found_root = j;
            }
        }
      
      mpc_get_cdpe (cdtmp, mroot[found_root]);
      cdpe_mod (rtmp, cdtmp);
      rdpe_mul (exp_drad, rtmp, eps);      
      rdpe_div_eq_d (min_dist, 1 + 4.0 * DBL_EPSILON);


      if ((!rdpe_le (min_dist, drad[found_root]) && !rdpe_gt (drad[found_root], exp_drad)) && !mps_context_get_over_max (s))
        {
          passed = false;
          
          if (getenv ("MPS_VERBOSE_TEST") && (strstr (pol->pol_file, getenv ("MPS_VERBOSE_TEST"))))
            {
              printf("Failing on root %d, with min_dist = ", found_root);
              rdpe_out_str (stdout, min_dist);
              printf("\ndrad_%d", found_root);
              rdpe_out_str (stdout, drad[found_root]);
              printf("\n");
              printf("Approximation_%d = ", found_root);
              mpc_out_str_2 (stdout, 10, -rdpe_Esp (drad[found_root]), -rdpe_Esp (drad[found_root]), mroot[found_root]);
              printf("\n");
            }
        }
    }

  if (zero_roots != mps_context_get_zero_roots (s))
    passed = false;
  
  fclose (result_stream);    
  
  mpc_clear (ctmp);   
  mpc_clear (root);
  mpc_vclear (mroot, mps_context_get_degree (s));
  
  free (mroot);
  free (drad);

  mps_polynomial_free (s, poly);
  mps_context_free (s);

  if (passed)
    success_test_message (pol->pol_file);
  else
    failed_test_message (pol->pol_file);

  if (getenv ("MPS_VERBOSE_TEST"))
    fail_unless (passed == true,
                 "Computed results are not exact to the required "
                 "precision.\n" "\n" " Dumping test configuration: \n"
                 "   => Polynomial file: %s;\n" "   => Required digits: %d\n"
                 "   => Gemignani's approach: %s;\n"
                 "   => Starting phase: %s;\n", pol->pol_file, pol->out_digits,
                 mps_boolean_to_string (pol->ga),
                 (pol->phase == float_phase) ? "float_phase" : "dpe_phase");
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

  test_polynomials = (test_pol **) malloc (sizeof (test_pol *) * 46);

  for (i = 0; i < 2; i++)
    {
      test_polynomials[n++] = test_pol_new_simple ("exp100", digits[i]);
      test_polynomials[n++] = test_pol_new_simple ("exp50", digits[i]);

      /* Floating point polynomials, don't require many digits */
      if (i == 0)
      {    
        test_polynomials[n++] = test_pol_new_simple ("kam1_1", digits[i]);
        test_polynomials[n++] = test_pol_new_simple ("kam1_2", digits[i]);
        test_polynomials[n++] = test_pol_new_simple ("kam1_3", digits[i]);
        test_polynomials[n++] = test_pol_new_simple ("kam2_1", digits[i]);
        test_polynomials[n++] = test_pol_new_simple ("kam2_2", digits[i]);
        test_polynomials[n++] = test_pol_new_simple ("kam2_3", digits[i]);
        test_polynomials[n++] = test_pol_new_simple ("kam3_1", digits[i]);
        test_polynomials[n++] = test_pol_new_simple ("kam3_2", digits[i]);
        test_polynomials[n++] = test_pol_new_simple ("kam3_3", digits[i]);
        test_polynomials[n++] = test_pol_new_simple ("lar1_200", digits[i]);
        test_polynomials[n++] = test_pol_new_simple ("lar1", digits[i]);
        test_polynomials[n++] = test_pol_new_simple ("lar2", digits[i]);
        test_polynomials[n++] = test_pol_new_simple ("lar3", digits[i]);
        test_polynomials[n++] = test_pol_new_simple ("lsr_24", digits[i]);
      }

      test_polynomials[n++] = test_pol_new_simple ("kir1_10", digits[i]);
      test_polynomials[n++] = test_pol_new_simple ("mand127", digits[i]);
      test_polynomials[n++] = test_pol_new_simple ("mand63", digits[i]);
      test_polynomials[n++] = test_pol_new_simple ("mand63", digits[i]);
      test_polynomials[n++] = test_pol_new_simple ("mand63", digits[i]);
      test_polynomials[n++] = test_pol_new_simple ("mig1_100", digits[i]);
      test_polynomials[n++] = test_pol_new_simple ("mig1_200", digits[i]);
      test_polynomials[n++] = test_pol_new_simple ("nroots50", digits[i]);
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
