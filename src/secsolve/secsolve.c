/*
 * secular.c
 *
 *  Created on: 11/apr/2011
 *      Author: leonardo
 */

#include <mps/interface.h>
#include <mps/secular.h>
#include <mps/core.h>
#include <string.h>

void
usage (mps_status * s, const char *program)
{
  /* If there is not an output stream do not print
   * the help */
  if (!s->outstr)
    return;

  fprintf (s->outstr,
           "Usage: %s [-dg] [-t type] [-n degree] [-o digits] [infile]\n"
           "\n"
           "Options:\n"
           " -d[domains] Activate debug on selected domains, that can be one of:\n"
           "               t: Trace\n"
           "               a: Approximation\n"
           "               c: Cluster\n"
           "               i: Improvement\n"
           "               w: Timings\n"
           "               o: Input/Output\n"
           "               m: Memory management\n"
           "               f: Function calls\n"
           "               Example: -dfi for function calls and improvement\n"
           " -g          Use Gemignani's approach\n"
           " -t type     Type can be 'f' for floating point\n"
           "             or 'd' for DPE\n"
           " -i          Isolate the roots only, do not perform approximation\n"
           " -o digits   Exact digits of the roots given as output.\n",
           program);

  exit (EXIT_FAILURE);
}

int
main (int argc, char **argv)
{
  mps_secular_equation *sec;
  mps_status *s;

  /* Create a new status */
  s = mps_status_new ();
  s->prec_in = 0;

  strncpy (s->goal, "aannc", 5);

  /* Gemignani's approach */
  mps_boolean ga = false;

  FILE *infile;

  /* Parse options */
  mps_opt *opt;
  mps_phase phase = float_phase;

  opt = NULL;
  while ((mps_getopts (&opt, &argc, &argv, "gid::t:o:")))
    {
      switch (opt->optchar)
        {
        case 'g':
          /* Gemignani's approach. Regenerate b_i after floating
           * point cycle */
          ga = true;
          break;
        case 'o':
          s->prec_out = atoi (opt->optvalue) * LOG2_10;
          break;
        case 'i':
          s->goal[0] = 'i';
          break;
        case 'd':
          s->DOLOG = true;
          s->logstr = stderr;

          if (!opt->optvalue)
            {
              /* If no specific debug domain has been specified, trace. */
              s->debug_level != MPS_DEBUG_TRACE;
              break;
            }

          /* If debugging was enabled, parse debug_level */
          while (*opt->optvalue)
            {
              switch (*opt->optvalue++)
                {
                case 't':
                  s->debug_level |= MPS_DEBUG_TRACE;
                  break;
                case 'a':
                  s->debug_level |= MPS_DEBUG_APPROXIMATIONS;
                  break;
                case 'c':
                  s->debug_level |= MPS_DEBUG_CLUSTER;
                  break;
                case 'i':
                  s->debug_level |= MPS_DEBUG_IMPROVEMENT;
                  break;
                case 'w':
                  s->debug_level |= MPS_DEBUG_TIMINGS;
                  break;
                case 'o':
                  s->debug_level |= MPS_DEBUG_IO;
                  break;
                case 'm':
                  s->debug_level |= MPS_DEBUG_MEMORY;
                  break;
                case 'f':
                  s->debug_level |= MPS_DEBUG_FUNCTION_CALLS;
                  break;
                default:
                  mps_error (s, 2, "Unrecgnozied debug option: %c",
                             opt->optvalue);
                  break;
                }
            }
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
              usage (s, argv[0]);
            }
          break;
        default:
          usage (s, argv[0]);
          break;
        }
    }

  if (argc > 2)
    usage (s, argv[0]);

  /* If no file is provided use standard input */
  if (argc == 1)
    infile = stdin;
  else
    infile = fopen (argv[1], "r");

  if (!infile)
    {
      mps_error (s, 1, "Cannot open input file for read, aborting.");
      return -1;
    }

  /* Create new secular equation */
  mps_parsing_configuration default_configuration = {
    /* .structure */
    MPS_STRUCTURE_COMPLEX_FP,

    /* .representation */
    MPS_REPRESENTATION_SECULAR
  };
  mps_parse_stream (s, infile, default_configuration);
  sec = s->secular_equation;

  /* Close the file if it's not stdin */
  if (argc == 2)
    fclose (infile);

  /* Set secular equation in user data, so it will be
   * accessible by the secular equation routines. */
  sec->starting_case = phase;

  /* Use always DPE with non floating point input */
  if (!MPS_STRUCTURE_IS_FP (s->secular_equation->input_structure))
    {
      sec->starting_case = dpe_phase;
    }

  if (phase == dpe_phase)
    s->skip_float = true;

  /* If we choose gemignani's approach follow it, otherwise
   * use standard mpsolve approach applied implicitly to the
   * secular equation. */
  mps_status_set_degree (s, sec->n);

  /* Set user polynomial with our custom functions */
  mps_allocate_data (s);

  if (ga)
    {
      /* Select the right algorithm */
      mps_status_select_algorithm (s, MPS_ALGORITHM_SECULAR_GA);
    }
  else
    {
      /* Select the right algorithm */
      mps_status_select_algorithm (s, MPS_ALGORITHM_SECULAR_MPSOLVE);
    }

  /* Solve the polynomial */
  mps_mpsolve (s);

  /* Output the roots */
  mps_copy_roots (s);
  mps_output (s);

  /* Free used data */
  mps_secular_equation_free (sec);
  mps_status_free (s);
}
