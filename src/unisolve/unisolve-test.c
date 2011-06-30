/*
 * unisolve-test.c
 *
 *  Created on: 15/giu/2011
 *      Author: leonardo
 */

#include <mps/interface.h>
#include <mps/core.h>
#include <mps/poly.h>
#include <gmp.h>
#include <mps/mpc.h>
#include <mps/gmptools.h>
#include <stdio.h>

void
usage ()
{
  return -1;
}

int main(int argc, char** argv)
{
  mpspoly_t poly;
  mps_status* s = mps_status_new();
  FILE* input_stream;
  FILE* check_stream;
  mps_boolean passed = true;
  mpc_t root, ctmp;
  mpf_t mroot;
  mpf_t ftmp;
  mpf_t eps;
  int out_digits = 50;
  int i, j, prec = out_digits * LOG2_10;
  int ch;

  if (argc != 3)
    {
      usage ();
      return -1;
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

  mps_set_default_values(s);
  s->prec_out = prec;
  strncpy(s->goal, "aannc", 5);
  mps_read_poly(s, input_stream, poly);

  mps_set_poly(s, poly);
  mps_allocate_data(s);
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
