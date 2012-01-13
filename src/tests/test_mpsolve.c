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

#define TEST_UNISOLVE(pol_name) (test_mpsolve (get_pol_file (pol_name, "unisolve"), get_res_file (pol_name, "unisolve"), MPS_ALGORITHM_STANDARD_MPSOLVE))
#define TEST_SECSOLVE_SECULAR(pol_name) (test_mpsolve (get_pol_file (pol_name, "secsolve"), get_res_file (pol_name, "secsolve"), MPS_ALGORITHM_SECULAR_MPSOLVE))
#define TEST_SECSOLVE_MONOMIAL(pol_name) (test_mpsolve (get_pol_file (pol_name, "unisolve"), get_res_file (pol_name, "unisolve"), MPS_ALGORITHM_SECULAR_GA))

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

  /* Check the roots */
  FILE* result_stream = fopen (res_file, "r"); 

  /* Create a new empty mps_status */
  mps_status * s = mps_status_new ();

  /* Load the polynomial that has been given to us */
  s->instr = fopen (pol_file, "r");
  
  if (!result_stream) 
    {
      fprintf (stderr, "Checking \033[1m%-30s\033[0m \033[31;1mno results file found!\033[0m\n", pol_file + 9); 
      return 0;
    }
  if (!s->instr)
    {
      fprintf (stderr, "Checking \033[1m%-30s\033[0m \033[31;1mno polinomial file found!\033[0m\n", pol_file + 9); 
      return 0;
    }
  
  fprintf (stderr, "Checking \033[1m%-30s\033[0m [\033[34;1mchecking\033[0m]", pol_file + 9);

  mps_parse_stream (s, s->instr);

  s->output_config->goal = MPS_OUTPUT_GOAL_APPROXIMATE;
  /* s->DOLOG = true; */
  /* s->debug_level = MPS_DEBUG_TRACE; */
  /* s->output_config->prec = 16;  */

  /* Solve it */
  mps_status_select_algorithm (s, algorithm);
  mps_mpsolve (s);
  
  mpc_init2 (root, s->data_prec_max);
  mpc_init2 (ctmp, s->data_prec_max); 
    
  /* Test if roots are equal to the roots provided in the check */   
  passed = true;

  for (i = 0; i < s->deg; i++)   
    {   
      rdpe_t rtmp;   
      cdpe_t cdtmp;   
      rdpe_t min_dist;
      int found_root = 0;
      
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
      
      mpc_sub (ctmp, root, s->mroot[0]);   
      mpc_get_cdpe (cdtmp, ctmp);   
      cdpe_mod (rtmp, cdtmp);   
      rdpe_set (min_dist, rtmp);   

      if (getenv ("MPS_VERBOSE_TEST") && (strstr (pol_file, getenv ("MPS_VERBOSE_TEST"))))
	{
	  printf ("Read root_%d = ", i);
	  mpc_out_str_2 (stdout, 10, s->data_prec_max, s->data_prec_max,
			 root);
	  printf ("\n");
	}
      
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

      /* printf ("min_dist_%d = ", i); */
      /* rdpe_out_str (stdout, min_dist); */
      /* printf ("\nrad_%d", i); */
      /* rdpe_out_str (stdout, s->drad[i]); */
      /* printf ("\n"); */
      
      if (!rdpe_le (min_dist, s->drad[found_root]) && !s->over_max)
	{
	  passed = false;
	  
	  if (getenv ("MPS_VERBOSE_TEST") && (strstr (pol_file, getenv ("MPS_VERBOSE_TEST"))))
	    {
	      printf("Failing on root %d, with min_dist = ", found_root);
	      rdpe_out_str (stdout, min_dist);
	      printf("\ndrad_%d", found_root);
	      rdpe_out_str (stdout, s->drad[found_root]);
	      printf("\n");
	      printf("Approximation_%d = ", found_root);
	      mpc_out_str_2 (stdout, 10, -rdpe_Esp (s->drad[found_root]), -rdpe_Esp (s->drad[found_root]), s->mroot[found_root]);
	      printf("\n");
	    }
	}
    }

  if (zero_roots != s->zero_roots)
    passed = false;
  
  fclose (result_stream);    
  
  mpc_clear (ctmp);   
  mpc_clear (root);   

  mps_status_free (s);

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

int
main (int argc, char ** argv)
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
  TEST_UNISOLVE ("kam4");
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
}
