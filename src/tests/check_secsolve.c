/*
 * Check for Secsolve
 */

#include <check.h>
#include <mps/mps.h>
#include <gmp.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <check_implementation.h>

/**
 * @brief Array with the polynomials to test
 */
test_pol **test_polynomials;

int test_secsolve_on_pol_impl (test_pol*, mps_output_goal, mps_boolean jacobi_iterations);

int 
test_secsolve_on_pol (test_pol * pol)
{
  return test_secsolve_on_pol_impl (pol, MPS_OUTPUT_GOAL_ISOLATE, false) &&
    test_secsolve_on_pol_impl (pol, MPS_OUTPUT_GOAL_APPROXIMATE, false); /* &&
    test_secsolve_on_pol_impl (pol, MPS_OUTPUT_GOAL_ISOLATE, true)       &&
    test_secsolve_on_pol_impl (pol, MPS_OUTPUT_GOAL_APPROXIMATE, true); */
}

/**
 * @brief This function tests the resolution of a polynomial file
 * referenced by <code>pol</code>.
 */
int
test_secsolve_on_pol_impl (test_pol * pol, mps_output_goal goal, mps_boolean jacobi_iterations)
{
  mpc_t root, ctmp;
  mps_boolean passed = true;
  int i, j, zero_roots = 0;
  char ch;
  rdpe_t eps;

  /* Check the roots */
  FILE* result_stream = fopen (pol->res_file, "r"); 
  FILE* input_stream  = fopen (pol->pol_file, "r");

  rdpe_set_2dl (eps, 1.0, - pol->out_digits);

  if (!result_stream) 
    {
      fprintf (stderr, "Checking \033[1m%-30s\033[0m \033[31;1mno results file found!\033[0m\n", 
               get_pol_name_from_path (pol->pol_file)); 
      return EXIT_FAILURE;
    }
  if (!input_stream)
    {
      fprintf (stderr, "Checking \033[1m%-30s\033[0m \033[31;1mno polinomial file found!\033[0m\n", 
               get_pol_name_from_path (pol->pol_file)); 
      return EXIT_FAILURE;
    }

  /* Create a new empty mps_context */
  mps_context * s = mps_context_new ();

  if ((getenv ("MPS_VERBOSE_TEST") && strstr (pol->pol_file, getenv ("MPS_VERBOSE_TEST"))) || pol->DOLOG)
    mps_context_set_debug_level (s, MPS_DEBUG_TRACE);

  /* Load the polynomial that has been given to us */
  mps_parse_stream (s, input_stream);
  
  fprintf (stderr, "Checking \033[1m%-30s\033[0m [\033[34;1mchecking\033[0m]", 
           get_pol_name_from_path (pol->pol_file));

  mps_context_set_output_prec (s, pol->out_digits);
  mps_context_set_output_goal (s, goal);
  mps_context_set_jacobi_iterations (s, jacobi_iterations);

  /* Solve it */
  mps_context_select_algorithm (s, (pol->ga) ? MPS_ALGORITHM_SECULAR_GA : MPS_ALGORITHM_STANDARD_MPSOLVE);
  mps_mpsolve (s);
  
  mpc_init2 (root, mps_context_get_data_prec_max (s));
  mpc_init2 (ctmp, mps_context_get_data_prec_max (s));
    
  /* Test if roots are equal to the roots provided in the check */   
  passed = true;

  rdpe_t * drad = rdpe_valloc (mps_context_get_degree (s));
  mpc_t * mroot = mpc_valloc (mps_context_get_degree (s));
  mpc_vinit2 (mroot, mps_context_get_degree (s), mps_context_get_data_prec_max (s));

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
      rdpe_mul_eq (rtmp, eps);
      rdpe_set (exp_drad, rtmp);
      
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

  if ((getenv ("MPS_VERBOSE_TEST") && strstr (pol->pol_file, getenv ("MPS_VERBOSE_TEST"))) || pol->DOLOG)
    {
      /* mps_context_set_output_format (s, MPS_OUTPUT_FORMAT_GNUPLOT_FULL); */
      mps_output (s);
    }
  
  fclose (result_stream);    
  
  mpc_clear (ctmp);   
  mpc_clear (root);
  mpc_vclear (mroot, mps_context_get_degree (s));
  
  free (mroot);
  free (drad);

  mps_context_free (s);

  if (passed)
    fprintf (stderr, "\rChecking \033[1m%-30s\033[0m [\033[32;1m  done  \033[0m]\n", 
             get_pol_name_from_path (pol->pol_file));
  else
    fprintf (stderr, "\rChecking \033[1m%-30s\033[0m [\033[31;1m failed \033[0m]\n", 
             get_pol_name_from_path (pol->pol_file));

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

START_TEST (test_secsolve)
{
  test_secsolve_on_pol (test_polynomials[_i]);
}
END_TEST

START_TEST (test_secsolve_altern)
{
  /* Start with testing floating point without ga */
  test_pol *pol =
    test_pol_new ("test100", "secsolve", 53, dpe_phase, false);
  test_secsolve_on_pol (pol);
  test_pol_free (pol);
}
END_TEST

START_TEST (test_secsolve_altern_ga)
{
  test_pol *pol = 
    test_pol_new ("test100", "secsolve", 53, dpe_phase, true);

  test_secsolve_on_pol (pol);
  test_pol_free (pol);
}
END_TEST

START_TEST (test_secsolve_integer)
{
  /* Test integer parsing of secsolve, ga approach */
  test_pol *pol = test_pol_new ("integer", "secsolve", 20, dpe_phase, true);
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
  test_pol *pol = test_pol_new ("wilk20", "secsolve", 53, float_phase, true);
  test_secsolve_on_pol (pol);
  test_pol_free (pol);

  /* Testing the wilkinson polynomial of degree 40 */
  pol = test_pol_new ("wilk40", "secsolve", 53, float_phase, true);
  test_secsolve_on_pol (pol);
  test_pol_free (pol);

  /* Testing the wilkinson polynomial of degree 80 */
  pol = test_pol_new ("wilk80", "secsolve", 53, float_phase, true);
  test_secsolve_on_pol (pol);
  test_pol_free (pol);
}
END_TEST

START_TEST (test_secsolve_wilkinson_monomial)
{
  /* Testinf the wilkinson polynomial of degree 20 */
  test_pol *pol = test_pol_new ("wilk20", "unisolve", 53, float_phase, true);
  test_secsolve_on_pol (pol);
  test_pol_free (pol);

  pol = test_pol_new ("wilk40", "unisolve", 53, float_phase, true);
  test_secsolve_on_pol (pol);
  test_pol_free (pol);

  pol = test_pol_new ("wilk80", "unisolve", 53, float_phase, true);
  test_secsolve_on_pol (pol);
  test_pol_free (pol);
}
END_TEST

START_TEST (test_secsolve_nroots)
{
  /* Testing secsolve on some polynomial of the type x^n - 1 */
  test_pol *pol = test_pol_new ("nroots50", "unisolve", 53, float_phase, true);
  test_secsolve_on_pol (pol);
  test_pol_free (pol);  
}
END_TEST

START_TEST (test_secsolve_kam1_1) 
{ 
  /* Testing the kam polynomials */ 
  test_pol *pol = test_pol_new ("kam1_1", "unisolve", 53, float_phase, true); 
  test_secsolve_on_pol (pol); 
  test_pol_free (pol); 
} 
END_TEST 

START_TEST (test_secsolve_kam1_2)
{
  test_pol *pol = test_pol_new ("kam1_2", "unisolve", 53, float_phase, true);
  test_secsolve_on_pol (pol);
  test_pol_free (pol);
}
END_TEST

START_TEST (test_secsolve_kam1_3)
{
  test_pol *pol = test_pol_new ("kam1_3", "unisolve", 53, float_phase, true);
  test_secsolve_on_pol (pol);
  test_pol_free (pol);
}
END_TEST

START_TEST (test_secsolve_kam2_1) 
{ 
  test_pol *pol = test_pol_new ("kam2_1", "unisolve", 53, float_phase, true); 
  test_secsolve_on_pol (pol); 
  test_pol_free (pol); 
} 
END_TEST 
 
START_TEST (test_secsolve_kam2_2)
{
  test_pol *pol = test_pol_new ("kam2_2", "unisolve", 53, float_phase, true);
  test_secsolve_on_pol (pol);
  test_pol_free (pol);
}
END_TEST

START_TEST (test_secsolve_kam2_3)
{
  test_pol *pol = test_pol_new ("kam2_3", "unisolve", 53, float_phase, true);
  test_secsolve_on_pol (pol);
  test_pol_free (pol);
}
END_TEST

START_TEST (test_secsolve_kam3_1) 
{ 
  test_pol *pol = test_pol_new ("kam3_1", "unisolve", 53, float_phase, true); 
  test_secsolve_on_pol (pol); 
  test_pol_free (pol); 
} 
END_TEST 

START_TEST (test_secsolve_kam3_2)
{
  test_pol *pol = test_pol_new ("kam3_2", "unisolve", 53, float_phase, true);
  test_secsolve_on_pol (pol);
  test_pol_free (pol);
}
END_TEST

START_TEST (test_secsolve_kam3_3)
{
  test_pol *pol = test_pol_new ("kam3_3", "unisolve", 53, float_phase, true);
  test_secsolve_on_pol (pol);
  test_pol_free (pol);
}
END_TEST

START_TEST (test_secsolve_mand)
{
  /* Testing secsolve on the mandelbrot polynomials */

  /* Mandelbrot classic, degree 63 */
  test_pol *pol = test_pol_new ("mand63", "unisolve", 53, float_phase, true);
  test_secsolve_on_pol (pol);
  test_pol_free (pol);

  /* Mandelbrot classic, degree 127 */
  pol = test_pol_new ("mand127", "unisolve", 53, float_phase, true);
  test_secsolve_on_pol (pol);
  test_pol_free (pol);
  
}
END_TEST

START_TEST (test_secsolve_exp)
{
  /* Testing secsolve on truncated exponential series */
  test_pol * pol = test_pol_new ("exp50", "unisolve", 53, float_phase, true);
  test_secsolve_on_pol (pol);
  test_pol_free (pol);

  pol = test_pol_new ("exp100", "unisolve", 53, float_phase, true);
  test_secsolve_on_pol (pol);
  test_pol_free (pol);
}
END_TEST

START_TEST (test_secsolve_mignotte)
{
  test_pol * pol = test_pol_new ("mig1_100", "unisolve", 53, float_phase, true);
  test_secsolve_on_pol (pol);
  test_pol_free (pol);

  pol = test_pol_new ("mig1_200", "unisolve", 53, float_phase, true);
  test_secsolve_on_pol (pol);
  test_pol_free (pol);
}
END_TEST

START_TEST (test_secsolve_mult)
{
  test_pol * pol = test_pol_new ("mult1", "unisolve", 15, float_phase, true);
  test_secsolve_on_pol (pol);
  test_pol_free (pol);
}
END_TEST

START_TEST (test_secsolve_mult_high_precision)
{
  test_pol * pol = test_pol_new ("mult1", "unisolve", 400, float_phase, true);
  test_secsolve_on_pol (pol);
  test_pol_free (pol);
}
END_TEST

START_TEST (test_secsolve_toep)
{
  test_pol * pol = test_pol_new ("toep1_128", "unisolve", 15, float_phase, true);
  test_secsolve_on_pol (pol);
  test_pol_free (pol);
}
END_TEST

START_TEST (test_secsolve_trv)
{
  test_pol * pol = test_pol_new ("trv_m", "unisolve", 150, float_phase, true);
  test_secsolve_on_pol (pol);
  test_pol_free (pol);
}
END_TEST

START_TEST (test_secsolve_kir1_10)
{
  test_pol * pol = test_pol_new ("kir1_10", "unisolve", 53, float_phase, true);
  test_secsolve_on_pol (pol);
  test_pol_free (pol);
}
END_TEST

START_TEST (test_secsolve_kir1_10_hp)
{  
  test_pol * pol = test_pol_new ("kir1_10", "unisolve", 50 * LOG2_10, float_phase, true);
  test_secsolve_on_pol (pol);
  test_pol_free (pol);
}
END_TEST

START_TEST (test_secsolve_kir1_20)
{
  test_pol * pol = test_pol_new ("kir1_20", "unisolve", 53, float_phase, true);
  test_secsolve_on_pol (pol);
  test_pol_free (pol);
}
END_TEST

START_TEST (test_secsolve_kir1_40) 
{ 
  test_pol * pol = test_pol_new ("kir1_40", "unisolve", 53, float_phase, true); 
  test_secsolve_on_pol (pol); 
  test_pol_free (pol); 
} 
END_TEST 

START_TEST (test_secsolve_spiral10)
{
  test_pol * pol = test_pol_new ("spiral10", "unisolve", 30, float_phase, true);
  test_secsolve_on_pol (pol);
  test_pol_free (pol);
}
END_TEST

START_TEST (test_secsolve_spiral20)
{
  test_pol * pol = test_pol_new ("spiral20", "unisolve", 30, float_phase, true);
  test_secsolve_on_pol (pol);
  test_pol_free (pol);
}
END_TEST

START_TEST (test_secsolve_spiral10_high_precision)
{
  test_pol * pol = test_pol_new ("spiral10", "unisolve", 500, float_phase, true);
  test_secsolve_on_pol (pol);
  test_pol_free (pol);
}
END_TEST

START_TEST (test_secsolve_mig1_200_high_precision)
{
  test_pol * pol = test_pol_new ("mig1_200", "unisolve", 1000, float_phase, true);
  test_secsolve_on_pol (pol);
  test_pol_free (pol);
}
END_TEST

START_TEST (test_secsolve_mig1_500_1)
{
  test_pol * pol = test_pol_new ("mig1_500_1", "unisolve", 53, float_phase, true);
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
  TCase *tc_secular = tcase_create ("Secular equation");

  /* Case of a_i = (-1)^(i+1) , b_i = i */
  tcase_add_test (tc_secular, test_secsolve_altern);
  tcase_add_test (tc_secular, test_secsolve_altern_ga);

  /* Integer parsing */
  tcase_add_test (tc_secular, test_secsolve_integer);

  /* Simple secular equation with cancellation problems */
  tcase_add_test (tc_secular, test_secsolve_simple);

  /* Wilkinson polynomials */
  tcase_add_test (tc_secular, test_secsolve_wilkinson);

  /* MONOMIAL TEST CASE */
  TCase *tc_monomial = tcase_create ("Monomial input");

  /* Roots of the unity */
  tcase_add_test (tc_monomial, test_secsolve_nroots);

  /* Kam polynomials */
  tcase_add_test (tc_monomial, test_secsolve_kam1_1); 
  tcase_add_test (tc_monomial, test_secsolve_kam1_2);
  tcase_add_test (tc_monomial, test_secsolve_kam1_3);
  tcase_add_test (tc_monomial, test_secsolve_kam2_1); 
  tcase_add_test (tc_monomial, test_secsolve_kam2_2);
  tcase_add_test (tc_monomial, test_secsolve_kam2_3);
  tcase_add_test (tc_monomial, test_secsolve_kam3_1);
  tcase_add_test (tc_monomial, test_secsolve_kam3_2);
  tcase_add_test (tc_monomial, test_secsolve_kam3_3);
                  
  /* Exponentials */
  tcase_add_test (tc_monomial, test_secsolve_exp);

  /* Mandelbrot polynomials */
  tcase_add_test (tc_monomial, test_secsolve_mand);

  /* Chebyshev */
  tcase_add_test (tc_monomial, test_secsolve_mignotte);

  /* Wilkinson polynomials */
  tcase_add_test (tc_monomial, test_secsolve_wilkinson_monomial);

  /* Mult* polynomials */
  tcase_add_test (tc_monomial, test_secsolve_mult);
  tcase_add_test (tc_monomial, test_secsolve_mult_high_precision);

  /* Topelitz */
  tcase_add_test (tc_monomial, test_secsolve_toep);

  /* Traverso, polynomial generated by resolution of a polynomial
   * system using Groebner elimination. */
  tcase_add_test (tc_monomial, test_secsolve_trv);

  /* Kirinnis polynomials */
  tcase_add_test (tc_monomial, test_secsolve_kir1_10);
  tcase_add_test (tc_monomial, test_secsolve_kir1_10_hp);
  tcase_add_test (tc_monomial, test_secsolve_kir1_20);
  tcase_add_test (tc_monomial, test_secsolve_kir1_40); 

  /* Spiral polynomials */
  tcase_add_test (tc_monomial, test_secsolve_spiral10);
  tcase_add_test (tc_monomial, test_secsolve_spiral20);
  tcase_add_test (tc_monomial, test_secsolve_spiral10_high_precision);

  /* Mig polynomial with high precision */
  tcase_add_test (tc_monomial, test_secsolve_mig1_200_high_precision);
  tcase_add_test (tc_monomial, test_secsolve_mig1_500_1);

  /* Add test case to the suite */
  suite_add_tcase (s, tc_secular);
  suite_add_tcase (s, tc_monomial);

  return s;
}

int
main (void)
{
  int number_failed, standard = 0;

  starting_setup ();

  /* Create a new test suite for secsolve and run it */
  Suite *s = secsolve_suite (standard);
  SRunner *sr = srunner_create (s);
  srunner_run_all (sr, CK_NORMAL);

  /* Get number of failed test and report */
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);

  return (number_failed != 0);
}
