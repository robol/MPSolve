/************************************************************
 **                                                        **
 **             __  __ ___  ___      _                     **
 **            |  \/  | _ \/ __| ___| |_ _____             **
 **            | |\/| |  _/\__ \/ _ \ \ V / -_)            **
 **            |_|  |_|_|  |___/\___/_|\_/\___|            **
 **                                                        **
 **       Multiprecision Polynomial Solver (MPSolve)       **
 **                 Version 2.9, April 2011                **
 **                                                        **
 **                      Written by                        **
 **                                                        **
 **     Dario Andrea Bini       <bini@dm.unipi.it>         **
 **     Giuseppe Fiorentino     <fiorent@dm.unipi.it>      **
 **     Leonardo Robol          <robol@mail.dm.unipi.it>   **
 **                                                        **
 **           (C) 2011, Dipartimento di Matematica         **
 ***********************************************************/

#define _MPS_PRIVATE
#include <mps/mps.h>
#include <string.h>

mps_status * s = NULL;

#ifndef __WINDOWS
#include <signal.h>


void
status (int signal)
{
  int i;
  FILE * logstr = stderr;

  fprintf (logstr, "\n\nDumping the approximations:\n");

  fprintf (logstr, 
	   "  Phase = %s, In = %d, Out = %d, Uncertain = %d, Zero = %d, Cluster = %ld\n",
	   MPS_PHASE_TO_STRING (s->lastphase), s->count[0], s->count[1], s->count[2],
	   s->zero_roots, s->clusterization->n);

  fprintf (logstr, "Current approximations:");
  for (i = 0; i < s->n; i++)
    {
      switch (s->lastphase)
        {
        case no_phase:
        case float_phase:
	  fprintf (logstr, "  Approximation  %4d = ", i);
	  cplx_outln_str (logstr, s->froot[i]);
          break;

        case dpe_phase:
	  fprintf (logstr, "  Approximation  %4d = ", i);
	  cdpe_outln_str (logstr, s->droot[i]);
          break;

        case mp_phase:
	  fprintf (logstr, "  Approximation  %4d = ", i);
	  mpc_outln_str (logstr, 10, s->mpwp, s->mroot[i]);
          break;
        }
    }

  /* output radii */
  fprintf (logstr, "Current radii: \n");
  for (i = 0; i < s->n; i++)
    {
      switch (s->lastphase)
        {
        case no_phase:
        case float_phase:
	  fprintf (logstr, "  Radius of root %4d = %e\n", i, s->frad[i]);
          break;

        case dpe_phase:
        case mp_phase:
	  fprintf (logstr, "  Radius of root %4d", i);
	  rdpe_outln_str (logstr, s->drad[i]);
          break;
        }
    }

  fprintf (logstr, "\n\n");
  fprintf (logstr, "Dumping status:\n\n");
  fprintf (logstr, "                Approximation              Attributes       Inclusion\n");
  for (i = 0; i < s->n; i++)
    {
      fprintf (logstr, "  Status  %4d: %-25s  %-15s  %-15s\n", i,
	       MPS_ROOT_STATUS_TO_STRING (s->root_status[i]),
	       MPS_ROOT_ATTRS_TO_STRING  (s->root_attrs[i]), 
	       MPS_ROOT_INCLUSION_TO_STRING (s->root_inclusion[i]));
    }
}

#undef _MPS_PRIVATE
#endif

void
usage (mps_status * s, const char *program)
{
  fprintf (stderr,
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
	   " -G goal     Select the goal to reach. Possible values are:\n"
	   "              a: Approximate the roots\n"
	   "              i: Isolate the roots\n"
	   "              c: Count the roots in the search set\n"
           " -g          Use Gemignani's approach\n"
           " -t type     Type can be 'f' for floating point\n"
	   " -j n        Number of threads to spawn as workers\n"
           "             or 'd' for DPE\n"
           " -o digits   Exact digits of the roots given as output.\n",
           program, program);

  exit (EXIT_FAILURE);
}

int
main (int argc, char **argv)
{
  /* Create a new status */
  s = mps_status_new ();

  /* Associate the SIGUSR1 signal to the mps_dump () function */
#ifndef __WINDOWS
  signal (SIGUSR1, status);
#endif

  mps_status_set_input_prec (s, 0);

  /* Gemignani's approach */
  mps_boolean ga = false;

  FILE *infile;

  /* Parse options */
  mps_opt *opt;
  mps_phase phase = no_phase;

  opt = NULL;
  while ((mps_getopts (&opt, &argc, &argv, "gG:d::t:o:O:j:")))
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
	      mps_status_set_output_format (s, MPS_OUTPUT_FORMAT_FULL);
	      break;
	    case 'b':
	      mps_status_set_output_format (s, MPS_OUTPUT_FORMAT_BARE);
	      break;
	    case 'g':
	      mps_status_set_output_format (s, MPS_OUTPUT_FORMAT_GNUPLOT);
	      if (*(opt->optvalue + 1) == 'f')
		mps_status_set_output_format (s, MPS_OUTPUT_FORMAT_GNUPLOT_FULL);
	      break;
	    case 'v':
	      mps_status_set_output_format (s, MPS_OUTPUT_FORMAT_VERBOSE);
	      break;
	    case 'c':
	      mps_status_set_output_format (s, MPS_OUTPUT_FORMAT_COMPACT);
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
          mps_status_set_output_prec (s, (atoi (opt->optvalue) + 1) * LOG2_10 + 1);
          break;
        case 'G':
	  switch (*opt->optvalue)
	    {
	    case 'a':
	      mps_status_set_output_goal (s, MPS_OUTPUT_GOAL_APPROXIMATE);
	      break;
	    case 'i':
	      mps_status_set_output_goal (s, MPS_OUTPUT_GOAL_ISOLATE);
	      break;
	    case 'c':
	      mps_status_set_output_goal (s, MPS_OUTPUT_GOAL_COUNT);
	      break;
	    default:
	      mps_error (s, 1, "The selected goal does not exists");
	      break;
	    }
          break;
        case 'd':
          if (!opt->optvalue)
            {
              /* If no specific debug domain has been specified, trace. */
	      mps_status_set_debug_level (s, MPS_DEBUG_TRACE);
              break;
            }

          /* If debugging was enabled, parse debug_level */
          while (*opt->optvalue)
            {	      
	      char output[255];
              switch (*opt->optvalue++)
                {
                case 't':
                  mps_status_add_debug_domain (s, MPS_DEBUG_TRACE);
                  break;
                case 'a':
		  mps_status_add_debug_domain (s, MPS_DEBUG_APPROXIMATIONS);
                  break;
                case 'c':
		  mps_status_add_debug_domain (s, MPS_DEBUG_CLUSTER);
                  break;
                case 'i':
		  mps_status_add_debug_domain (s, MPS_DEBUG_IMPROVEMENT);
                  break;
                case 'w':
		  mps_status_add_debug_domain (s, MPS_DEBUG_TIMINGS);
                  break;
                case 'o':
		  mps_status_add_debug_domain (s, MPS_DEBUG_IO);
                  break;
                case 'm':
		  mps_status_add_debug_domain (s, MPS_DEBUG_MEMORY);
                  break;
                case 'f':
		  mps_status_add_debug_domain (s, MPS_DEBUG_FUNCTION_CALLS);
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

	case 'j':
	  mps_thread_pool_set_concurrency_limit (s, NULL, atoi (opt->optvalue));
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

  /* Parse the input stream and if a polynomial is given as output, 
   * allocate also a secular equation to be used in regeneration */
  mps_parse_stream (s, infile);

  /* Close the file if it's not stdin */
  if (argc == 2)
    fclose (infile);

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
  mps_status_set_starting_phase (s, phase);

  /* Solve the polynomial */
  mps_mpsolve (s);

  /* Output the roots */
  mps_output (s);

  /* Free used data */
  mps_status_free (s);
}
