/*
 * secular.c
 *
 *  Created on: 11/apr/2011
 *      Author: leonardo
 */

#include <mps/interface.h>
#include <mps/secular.h>
#include <mps/core.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>

void
usage(mps_status *s, const char* program)
{
  /* If there is not an output stream do not print
   * the help */
  if (!s->outstr)
    return;

  fprintf(s->outstr,
      "Usage: %s [-dg] [-t type] [-n degree] [-o digits] [infile]\n"
      "\n"
      "Options:\n"
      " -d          Activate debug\n"
      " -g          Use Gemignani's approach\n"
      " -t type     Type can be 'f' for floating point\n"
      "             or 'd' for DPE\n"
      " -n degree   Degree of the polynomial associated\n"
      "             associated with the secular equation.\n"
      " -o digits   Exact digits of the roots given as output.\n",
      program);

  exit (EXIT_FAILURE);
}

mps_secular_equation*
read_secular_equation(mps_status* s, FILE* input_stream)
{
  mps_secular_equation* sec;
  int n, r, i;

  /* Read the number of the coefficients */
  r = fscanf(input_stream, "%d", &n);

  if (!r)
      mps_error(s, 1, "Error reading input coefficients of the secular equation.\n");

  /* Read directly the secular equation in DPE, so we don't need
   * to have a fallback case if the coefficients are bigger than
   * what is supported by the standard floating point arithmetic */
  sec = mps_secular_equation_new_raw(s, n);

  for(i = 0; i < n; i++)
    {
      rdpe_inp_str_flex(cdpe_Re(sec->adpc[i]), input_stream);
      rdpe_inp_str_flex(cdpe_Im(sec->adpc[i]), input_stream);
      rdpe_inp_str_flex(cdpe_Re(sec->bdpc[i]), input_stream);
      rdpe_inp_str_flex(cdpe_Im(sec->bdpc[i]), input_stream);
    }

  /* Deflate input, if identical b_i coefficients are found */
  mps_secular_deflate(s, sec);

  /* Copy coefficients back in other places */
  for(i = 0; i < sec->n; i++)
    {
      mpc_set_cdpe(sec->ampc[i], sec->adpc[i]);
      mpc_set_cdpe(sec->bmpc[i], sec->bdpc[i]);

      /* Get floating points coefficients */
      cdpe_get_x(sec->afpc[i], sec->adpc[i]);
      cdpe_get_x(sec->bfpc[i], sec->bdpc[i]);
    }

  return sec;
}

int
main(int argc, char** argv)
{
  mps_secular_equation* sec;
  mps_status* s;
  int i;

  /* Create a new secular equation with some random coefficients */
  unsigned int n = 5;
  s = mps_status_new();

  /* Gemignani's approach */
  mps_boolean ga = false;

  FILE* infile;
  double tmp1, tmp2;

  /* Parse options */
  mps_opt* opt;
  mps_phase phase = float_phase;
  while ((opt = mps_getopts(&argc, &argv, "gn:dt:o:")))
    {
      switch (opt->optchar)
        {
      case 'g':
        /* Gemignani's approach. Regenerate b_i after floating
         * point cycle */
        ga = true;
        break;
      case 'o':
        s->prec_out = atoi(opt->optvalue) * LOG2_10;
        break;
      case 'n':
        if (opt->optvalue)
          n = atoi(opt->optvalue);
        else
          usage(s, argv[0]);
        break;
      case 'd':
        s->DOLOG = true;
        s->logstr = stderr;
        break;
      case 't':
        switch (opt->optvalue[0])
        {
        case 'f':
          phase = float_phase;
          break;
        case 'd':
          phase = dpe_phase;
          break;
        default:
          usage(s, argv[0]);
        }
        break;
      default:
        usage(s, argv[0]);
        break;
        }

      free(opt);
    }

  if (argc > 2)
    usage(s, argv[0]);

  /* If no file is provided use standard input */
  if (argc == 1)
    infile = stdin;
  else
    infile = fopen(argv[1], "r");

  /* Create new secular equation */
  sec = read_secular_equation(s, infile);
  if (argc == 2)
    fclose (infile);

  /* Set secular equation in user data, so it will be
   * accessible by the secular equation routines. */
  s->user_data = sec;

  if (phase == dpe_phase)
      sec->starting_case = 'd';
  else
    sec->starting_case = 'f';

  /* If we choose gemignani's approach follow it, otherwise
   * use standard mpsolve approach applied implicitly to the
   * secular equation. */
  if (ga)
    {
      s->computation_style = 'g';

      /* Set degree and allocate polynomial-related variables
       * to allow initializitation to be performed. */
      s->deg = s->n = sec->n;
      mps_allocate_poly_inplace(s, sec->n);
      s->data_type = "uri";

      /* We set the selected phase */
      s->lastphase = phase;

      /* Allocate other data */
      mps_allocate_data(s);

      /* Manually set FILE* pointer for streams.
       * More refined options will be added later. */
      s->outstr = s->rtstr = stdout;

      /* Solve the secular equation */
      mps_secular_ga_mpsolve(s, phase);
    }
  else
    {
      s->computation_style = 'm';

      /* Set user polynomial with our custom functions */
      mps_status_set_poly_u(s, sec->n, MPS_FNEWTON_PTR(mps_secular_fnewton),
          MPS_DNEWTON_PTR(mps_secular_dnewton),
          MPS_MNEWTON_PTR(mps_secular_mnewton));

      /* Check data routine */
      s->check_data_usr = MPS_CHECK_DATA_PTR(mps_secular_check_data);

      /* Set starting point custom routine */
      s->fstart_usr = MPS_FSTART_PTR(mps_secular_fstart);
      s->dstart_usr = MPS_DSTART_PTR(mps_secular_dstart);

      /* Solve the polynomial */
      s->goal[0] = 'a';
      mps_mpsolve(s);

    }

  /* Output the roots */
  mps_copy_roots(s);
  mps_output(s);

  /* Free used data */
  mps_secular_equation_free(sec);
  mps_status_free(s);
}
