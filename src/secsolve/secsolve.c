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
#include <mcheck.h>

void
abortfn (enum mcheck_status status)
{
  fprintf (stderr, "A memory error has occurred in MPSolve; aborting\n");
  abort ();
}


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
	   " -O format   Select format for output:\n"
	   "               f: Full output\n"
	   "               b: Bare output\n"
	   "               c: Compact output\n"
	   "               v: Verbose output\n"
	   "               g: Gnuplot-ready output\n"
	   "               gf: Gnuplot-full mode, can be piped to gnuplot -persist. For example:\n"
	   "                   %s -g -Ogf myfile.pol | gnuplot \n"
           " -g          Use Gemignani's approach\n"
           " -t type     Type can be 'f' for floating point\n"
           "             or 'd' for DPE\n"
           " -i          Isolate the roots only, do not perform approximation\n"
           " -o digits   Exact digits of the roots given as output.\n",
           program, program);

  exit (EXIT_FAILURE);
}

int
main (int argc, char **argv)
{
  mps_status *s;

  /* Create a new status */
  s = mps_status_new ();
  s->input_config->prec = 0;

  strncpy (s->goal, "aannc", 5);

  /* Gemignani's approach */
  mps_boolean ga = false;

  FILE *infile;

  /* Parse options */
  mps_opt *opt;
  mps_phase phase = float_phase;

  opt = NULL;
  while ((mps_getopts (&opt, &argc, &argv, "gid::t:o:O:")))
    {
      switch (opt->optchar)
        {
	case 'O':
	  /* Select the desired output format */
	  if (!opt->optvalue)
	    mps_error (s, 1, "An argument is needed for option 'O'");
	  
	  switch (*opt->optvalue)
	    {
	    case 'f':
	      s->output_config->format = MPS_OUTPUT_FORMAT_FULL;
	      break;
	    case 'b':
	      s->output_config->format = MPS_OUTPUT_FORMAT_BARE;
	      break;
	    case 'g':
	      s->output_config->format = MPS_OUTPUT_FORMAT_GNUPLOT;
	      if (*(opt->optvalue + 1) == 'f')
		s->output_config->format = MPS_OUTPUT_FORMAT_GNUPLOT_FULL;
	      break;
	    case 'v':
	      s->output_config->format = MPS_OUTPUT_FORMAT_VERBOSE;
	      break;
	    case 'c':
	      s->output_config->format = MPS_OUTPUT_FORMAT_COMPACT;
	      break;
	    default:
	      mps_error (s, 1, "The selected output format is not supported");
	      break;
	    }
	    
        case 'g':
          /* Gemignani's approach. Regenerate b_i after floating
           * point cycle */
          ga = true;
          break;
        case 'o':
          s->output_config->prec = (atoi (opt->optvalue) + 1) * LOG2_10 + 1;
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
              s->debug_level |= MPS_DEBUG_TRACE;
              break;
            }

          /* If debugging was enabled, parse debug_level */
          while (*opt->optvalue)
            {	      
	      char output[255];
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
		  sprintf (output, "Unrecognized debug option: %c", *(opt->optvalue - 1));
                  mps_error (s, 1, output);
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
  s->input_config->structure  = MPS_STRUCTURE_COMPLEX_FP;
  s->input_config->representation = MPS_REPRESENTATION_SECULAR;

  /* Parse the input stream and if a polynomial is given as output, 
   * allocate also a secular equation to be used in regeneration */
  mps_parse_stream (s, infile);

  /* Close the file if it's not stdin */
  if (argc == 2)
    fclose (infile);

  /* Select ga if we are in the case of monomial input, since is the only algorithm
   * that we can use. */
  if (MPS_INPUT_CONFIG_IS_MONOMIAL (s->input_config) && !ga)
    {
      ga = true;
      MPS_DEBUG_WITH_INFO (s, "Selecting algorithm MPS_ALGORITHM_SECULAR_GA since MPS_ALGORITHM_SECULAR_MPSOLVE is not available for monomial input");
    }

  /* If we choose gemignani's approach follow it, otherwise
   * use standard mpsolve approach applied implicitly to the
   * secular equation. */
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

  /* Select the starting phase according to user input */
  s->input_config->starting_phase = phase;

  /* Solve the polynomial */
  mps_mpsolve (s);

  /* Output the roots */
  mps_output (s);

  /* Free used data */
  mps_status_free (s);
}
