/*
 * secsolve-test.c
 *
 *  Created on: 15/giu/2011
 *      Author: leonardo
 */

#include <mps/interface.h>
#include <mps/core.h>
#include <mps/secular.h>
#include <gmp.h>
#include <mps/mpc.h>
#include <mps/gmptools.h>
#include <stdio.h>
#include <ctype.h>

int
usage (const char* program_name)
{
  fprintf(stderr,
          "%s [OPTIONS] secular_equation_file results_file\n"
          "\n"
          "Options:\n"
          "\n"
          " -g\tUse Gemignani's approach\n"
          " -oN\tCheck roots computation with N guaranteed digits\n"
          " -t[f,d]\tSelect floating point or dpe computation\n"
          "\n",
          program_name);
  exit(EXIT_FAILURE);
}

int main(int argc, char** argv)
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
  char starting_case = 'f';

  /* Select if we need to test the Gemignani's approach */
  mps_boolean ga = false;

  /* Output digit to test, default values. Script should normally
   * alter this with options on the command line. */
  int out_digits = 100;
  int i, j, prec = out_digits * LOG2_10;
  int ch;

  /* Parse the options */
  mps_opt* opt;
  while ((opt = mps_getopts(&argc, &argv, "o:t:g")))
    {
      switch(opt->optchar)
      {
        case 'g':
          ga = true;
          break;
        case 'o':
          out_digits = atoi(opt->optvalue);
          prec = out_digits * LOG2_10;
          break;
        case 't':
          starting_case = opt->optvalue[0];
          break;
        default:
          usage(argv[0]);
          break;
      }

      free(opt);
    }

  /* Check if we have the file with the coefficients and the one with the
   * results */
  if (argc != 3)
    {
      usage (argv[0]);
      exit(EXIT_FAILURE);
    }

  mpc_init2(root, prec);
  mpc_init2(ctmp, prec);
  mpf_init2(mroot, prec);
  mpf_init2(eps,   prec);
  mpf_init2(ftmp,  prec);

  mpf_set_2dl(eps, 1.0, -out_digits);

  /* Open streams */
  input_stream = fopen(argv[1], "r");
  check_stream = fopen(argv[2], "r");

  if (!(input_stream && check_stream))
    {
      fprintf(stderr, "Cannot open one or more input files\n");
      return -1;
    }

  /* Some default values */
  mps_set_default_values(s);
  s->prec_out = prec;
  strncpy(s->goal, "aannc", 5);

  /* Set secular equation and start in floating point */
  sec = mps_secular_equation_read_from_stream(s, input_stream);
  s->secular_equation = sec;
  sec->starting_case = (starting_case == 'f') ? float_phase : dpe_phase;

  mps_status_set_degree(s, sec->n);

  if (!ga)
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

  if (!passed) {
      fprintf(stderr, "Computed results are not exact to the required precision,\n"
              "that is of %d digits.\n", out_digits);
  }

  return !passed;
}
