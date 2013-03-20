#include <mps/mps.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <mcheck.h>
#include <string.h>
#include "check_implementation.h"

#define MPS_TEST_UNISOLVE           0
#define MPS_TEST_SECSOLVE_SECULAR   1
#define MPS_TEST_SECSOLVE_GA        2

#define TEST_UNISOLVE(pol_name) {\
  char * pol_file = get_pol_file (pol_name, "unisolve"); \
  char * res_file = get_res_file (pol_name, "unisolve"); \
  test_mpsolve (pol_file, res_file, MPS_ALGORITHM_STANDARD_MPSOLVE);    \
  free (pol_file); \
  free (res_file); \
  }

#define TEST_SECSOLVE_SECULAR(pol_name) {\
  char * pol_file = get_pol_file (pol_name, "secsolve"); \
  char * res_file = get_res_file (pol_name, "secsolve"); \
  test_mpsolve (pol_file, res_file, MPS_ALGORITHM_STANDARD_MPSOLVE);    \
  free (pol_file); \
  free (res_file); \
  }

#define TEST_SECSOLVE_MONOMIAL(pol_name) {\
  char * pol_file = get_pol_file (pol_name, "unisolve"); \
  char * res_file = get_res_file (pol_name, "unisolve"); \
  test_mpsolve (pol_file, res_file, MPS_ALGORITHM_SECULAR_GA);  \
  free (pol_file); \
  free (res_file); \
  }

static mps_boolean debug = false;

void
test_header (const char * header, const char * description)
{
  fprintf (stderr, "\n *** \033[1mTEST:\033[0m \033[34;1m%s\033[0m *** \n", header);
  if (description)
    {
      fprintf (stderr, " %s\n\n", description);
    }
  else
    fprintf (stderr, "\n");
}

int 
test_mpsolve (char * pol_file, char * res_file, mps_algorithm algorithm)
{
  mpc_t root, ctmp;
  mps_boolean passed = true;
  int i, j, zero_roots = 0;
  char ch;
  rdpe_t eps;

  /* Check the roots */
  FILE* result_stream = fopen (res_file, "r"); 
  FILE* input_stream  = fopen (pol_file, "r");

  if (!result_stream) 
    {
      fprintf (stderr, "Checking \033[1m%-30s\033[0m \033[31;1mno results file found!\033[0m\n", pol_file + 9); 
      return EXIT_FAILURE;
    }
  if (!input_stream)
    {
      fprintf (stderr, "Checking \033[1m%-30s\033[0m \033[31;1mno polinomial file found!\033[0m\n", pol_file + 9); 
      return EXIT_FAILURE;
    }

  /* Create a new empty mps_context */
  mps_context * s = mps_context_new ();

  if (debug)
    mps_context_set_debug_level (s, MPS_DEBUG_TRACE);

  /* Load the polynomial that has been given to us */
  mps_parse_stream (s, input_stream);
  
  fprintf (stderr, "Checking \033[1m%-30s\033[0m [\033[34;1mchecking\033[0m]", pol_file + 9);

  mps_context_set_output_goal (s, MPS_OUTPUT_GOAL_ISOLATE);
  mps_context_set_output_prec (s, 50 * LOG2_10);
  rdpe_set_dl (eps, 1.0, -15);

  /* Solve it */
  mps_context_select_algorithm (s, algorithm);
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

      if (getenv ("MPS_VERBOSE_TEST") && (strstr (pol_file, getenv ("MPS_VERBOSE_TEST"))))
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

      /* printf ("min_dist_%d = ", i); */
      /* rdpe_out_str (stdout, min_dist); */
      /* printf ("\nrad_%d", i); */
      /* rdpe_out_str (stdout, s->drad[i]); */
      /* printf ("\n"); */


      mpc_get_cdpe (cdtmp, mroot[found_root]);
      cdpe_mod (rtmp, cdtmp);
      rdpe_mul_eq (rtmp, eps);
      rdpe_set (exp_drad, rtmp);
      
      if ((!rdpe_le (min_dist, drad[found_root]) && !rdpe_gt (drad[found_root], exp_drad)) && !mps_context_get_over_max (s))
        {
          passed = false;
          
          if (getenv ("MPS_VERBOSE_TEST") && (strstr (pol_file, getenv ("MPS_VERBOSE_TEST"))))
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
  
  fclose (input_stream);
  fclose (result_stream);    
  
  mpc_clear (ctmp);   
  mpc_clear (root);
  mpc_vclear (mroot, mps_context_get_degree (s));
  
  free (mroot);
  free (drad);

  if (getenv ("MPS_VERBOSE_TEST") && (strstr (pol_file, getenv ("MPS_VERBOSE_TEST"))))
    {
      mps_context_set_output_format (s, MPS_OUTPUT_FORMAT_GNUPLOT_FULL);
      mps_output (s);
    }

  mps_context_free (s);

  if (passed)
    fprintf (stderr, "\rChecking \033[1m%-30s\033[0m [\033[32;1m  done  \033[0m]\n", pol_file + 9);
  else
    fprintf (stderr, "\rChecking \033[1m%-30s\033[0m [\033[31;1m failed \033[0m]\n", pol_file + 9);

  return passed;
}

void
abortfn (enum mcheck_status status)
{
  switch (status)
    {
    case MCHECK_DISABLED:
      break;
    case MCHECK_OK:
      printf ("Questo blocco funziona\n");
      break;
    case MCHECK_FREE:
      printf("Ok, ho beccato un errore di memoria: double free\n");
      break;
    case MCHECK_HEAD:
      printf("Ok, ho beccato un errore di memoria: hai letto prima di un blocco allocato\n");
      break;
    case MCHECK_TAIL:
      printf("Ok, ho beccato un errore di memoria: hai letto in coda ad un blocco allocato\n");
      break;
    }
  abort ();
}

void
standard_tests (void)
{
  test_header ("Unisolve - Classic approach", 
               "Testing polynomial solving with MPSolve classic approach");

  /* Some unisolve tests */
  TEST_UNISOLVE ("mand63");
  TEST_UNISOLVE ("kam2_1");
  TEST_UNISOLVE ("kir1_10");
  TEST_UNISOLVE ("nroots50");
  TEST_UNISOLVE ("wilk40");
  TEST_UNISOLVE ("lar3");
  TEST_UNISOLVE ("trv_m");
  TEST_UNISOLVE ("lar2"); 
  TEST_UNISOLVE ("toep1_128");
  TEST_UNISOLVE ("mult1");
  TEST_UNISOLVE ("exp50");
  TEST_UNISOLVE ("mand127");
  TEST_UNISOLVE ("kam3_2");
  TEST_UNISOLVE ("kam3_3");
  TEST_UNISOLVE ("umand31");
  TEST_UNISOLVE ("kam1_2");
  TEST_UNISOLVE ("lar1");
  TEST_UNISOLVE ("spiral20"); 
  TEST_UNISOLVE ("wilk20");
  TEST_UNISOLVE ("test"); 
  TEST_UNISOLVE ("lar1_200");
  TEST_UNISOLVE ("kam1_3");
  TEST_UNISOLVE ("kam3_1");
  TEST_UNISOLVE ("kam2_2");
  TEST_UNISOLVE ("exp100");
  TEST_UNISOLVE ("kam1_1");
  TEST_UNISOLVE ("kam2_3");
  TEST_UNISOLVE ("mig1_100");
  TEST_UNISOLVE ("wilk80");
  TEST_UNISOLVE ("mig1_200");
  TEST_UNISOLVE ("lsr_24");

  test_header ("Secsolve - Solving secular equation", 
               "Solving secular equation using MPSolve user-polynomial feature");

  /* Normal secular tests */
  TEST_SECSOLVE_SECULAR ("rand120");
  TEST_SECSOLVE_SECULAR ("rand15");

  test_header ("Secsolve - Solving polynomials",
               "Solving polynomials by representing them as secular equations");

  /* Roots of unity */
  TEST_SECSOLVE_MONOMIAL ("nroots50");

  /* Mandelbrot polynomials */
  TEST_SECSOLVE_MONOMIAL ("mand63");
  TEST_SECSOLVE_MONOMIAL ("mand127");

  /* Mignotte polynomials */
  TEST_SECSOLVE_MONOMIAL ("mig1_100");
  TEST_SECSOLVE_MONOMIAL ("mig1_200");

  /* Wilkinson polynomials */
  TEST_SECSOLVE_MONOMIAL ("wilk20");
  TEST_SECSOLVE_MONOMIAL ("wilk40");
  TEST_SECSOLVE_MONOMIAL ("wilk80");

  TEST_SECSOLVE_MONOMIAL ("spiral20");
  TEST_SECSOLVE_MONOMIAL ("mult1");
}

int 
main (int argc, char ** argv)
{
  if (argc == 1)
    standard_tests ();

  /* If that's not the case, parse options */
  mps_algorithm alg = MPS_ALGORITHM_STANDARD_MPSOLVE;
  mps_opt * opt = NULL;
  while (mps_getopts (&opt, &argc, &argv, "a::d"))
    {
      switch (opt->optchar)
        {
        case 'a':
          if (opt->optvalue)
            {
              if (*opt->optvalue == 'g')
                {
                  /* fprintf (stderr, "SECSOLVE -g mode\n"); */
                  alg = MPS_ALGORITHM_SECULAR_GA;
                }
              if (*opt->optvalue == 'u')
                {
                  /* fprintf (stderr, "UNISOLVE mode\n"); */
                  alg = MPS_ALGORITHM_STANDARD_MPSOLVE;
                }
            }
          break;
        case 'd':
          debug = true;
        default:
          break;
        }
    }

  if (argc == 3)
    {
      return (test_mpsolve (argv[1], argv[2], alg) != true);
    }

  return EXIT_SUCCESS;
}
