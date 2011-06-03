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

int
main(int argc, char** argv)
{
  int i;

  /* Create a new secular equation with some random coefficients */
  unsigned int n = 5;

  /* Gemignani's approach */
  mps_boolean ga = false;

  /* Allocate space for the complex coefficients */
  cplx_t* a_coefficients;
  cplx_t* b_coefficients;

  mps_secular_equation* sec;
  mps_status* s;

  FILE* infile;
  double tmp1, tmp2;

  s = mps_status_new();

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

  /* If a file is provided read the coefficients from there */
  if(argc == 2) {
    infile = fopen (argv[1], "r");
    fscanf(infile, "%d", &n);

    a_coefficients = cplx_valloc (n);
    b_coefficients = cplx_valloc (n);

    for (i = 0; i < n; i++) {
      fscanf (infile, "%lf", &tmp1);
      fscanf (infile, "%lf", &tmp2);
      cplx_set_d(a_coefficients[i], tmp1, tmp2);
      fscanf (infile, "%lf", &tmp1);
      fscanf (infile, "%lf", &tmp2);
      cplx_set_d(b_coefficients[i], tmp1, tmp2);
    }
  } else {

    /* Allocate space for the coefficients */
    a_coefficients = cplx_valloc(n);
    b_coefficients = cplx_valloc(n);

    /* Generate coefficients */
    srand(time(NULL));
    for (i = 0; i < n; i++)
      {
        cplx_set_d(a_coefficients[i], drand(), drand());
        cplx_set_d(b_coefficients[i], drand(), drand());
        //      cplx_set_d(a_coefficients[i], pow(-1, (double) i + 1), 0);
        //      cplx_set_d(b_coefficients[i], 1.0 / (i + 1) / (i + 1), 0);
      }

  }

  /* Create new secular equation */
  sec = mps_secular_equation_new(a_coefficients, b_coefficients, n);

  /* Set secular equation in user data, so it will be
   * accessible by the secular equation routines. */
  s->user_data = sec;

  /* If we choose gemignani's approach follow it, otherwise
   * use standard mpsolve approach applied implicitly to the
   * secular equation. */
  if (ga)
    {
      /* Set degree and allocate polynomial-related variables
       * to allow initializitation to be performed. */
      s->deg = s->n = n;
      mps_allocate_poly_inplace(s, n);
      s->data_type = "uri";

      /* We set the selected phase */
      s->lastphase = phase;

      /* Allocate other data */
      mps_allocate_data(s);

      /* Manually set FILE* pointer for streams.
       * More refined options will be added later. */
      s->logstr = s->outstr = s->rtstr = stdout;

      /* Solve the secular equation */
      mps_secular_ga_mpsolve(s, phase);
    }
  else
    {
      /* Set user polynomial with our custom functions */
      mps_status_set_poly_u(s, n, MPS_FNEWTON_PTR(mps_secular_fnewton),
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
