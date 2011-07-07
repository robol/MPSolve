/*
 * unisolve-test.c
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
usage ()
{
  return -1;
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

  /* Select if we need to test the Gemignani's approach */
  mps_boolean ga = false;

  /* Output digit to test */
  int out_digits = 50;
  int i, j, prec = out_digits * LOG2_10;
  int ch;

  if (argc < 4 || argc > 5)
    {
      usage ();
      return -1;
    }

  /* If the fifth argument was passed activate
   * Gemignani's approach */
  if (argc == 5)
      ga = true;

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
  sec->starting_case = (argv[3][0] == 'f') ? 'f' : 'd';

  mps_status_set_degree(s, sec->n);
  mps_allocate_data(s);
  mps_select_algorithm(s, MPS_ALGORITHM_SECULAR_MPSOLVE);
  
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

  return !passed;
}
