/*
 * This file is part of MPSolve 3.2.1
 * 
 * Copyright (C) 2001-2020, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors: 
 *   Leonardo Robol <leonardo.robol@unipi.it>
 */

#define _MPS_PRIVATE
#include <mps/mps.h>
#include <string.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#if HAVE_GRAPHICAL_DEBUGGER
#define MPSOLVE_GETOPT_STRING "a:G:D:d::xt:o:O:j:S:O:i:vl:bp:rs:c"
#else
#define MPSOLVE_GETOPT_STRING "a:G:D:d::t:o:O:j:S:O:i:vl:bp:rs:c"
#endif

#if HAVE_GRAPHICAL_DEBUGGER
 #include <gtk/gtk.h>
 #include <iteration-logger.h>

MpsIterationLogger * logger = NULL;
static mps_boolean logger_closed = false;
#endif

mps_context * s = NULL;
mps_polynomial * poly = NULL;

#ifndef __WINDOWS
#include <signal.h>

void
status (int signal)
{
  int i;
  FILE * logstr = stderr;

  fprintf (stderr, "\nOperation running now: %s\n\n", MPS_OPERATION_TO_STRING (s->operation));
  fprintf (logstr, "Dumping the approximations:\n");

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
          cplx_outln_str (logstr, s->root[i]->fvalue);
          break;

        case dpe_phase:
          fprintf (logstr, "  Approximation  %4d = ", i);
          cdpe_outln_str (logstr, s->root[i]->dvalue);
          break;

        case mp_phase:
          fprintf (logstr, "  Approximation  %4d = ", i);
          mpc_outln_str (logstr, 10, s->mpwp, s->root[i]->mvalue);
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
          fprintf (logstr, "  Radius of root %4d = %e\n", i, s->root[i]->frad);
          break;

        case dpe_phase:
        case mp_phase:
          fprintf (logstr, "  Radius of root %4d", i);
          rdpe_outln_str (logstr, s->root[i]->drad);
          break;
        }
    }

  fprintf (logstr, "\n\n");
  fprintf (logstr, "Dumping status:\n\n");
  fprintf (logstr, "                Approximation              Attributes       Inclusion\n");
  for (i = 0; i < s->n; i++)
    {
      fprintf (logstr, "  Status  %4d: %-25s  %-15s  %-15s\n", i,
               MPS_ROOT_STATUS_TO_STRING (s->root[i]->status),
               MPS_ROOT_ATTRS_TO_STRING  (s->root[i]->attrs), 
               MPS_ROOT_INCLUSION_TO_STRING (s->root[i]->inclusion));
    }

  fprintf (stderr, "\n\nOperation running now: %s\n", MPS_OPERATION_TO_STRING (s->operation));
}

#undef _MPS_PRIVATE
#endif

void
usage (mps_context * s, const char *program)
{
  fprintf (stdout,
           "%s [-a alg] [-b] -c [-G goal] [-o digits] [-i digits] [-j n] [-t type] [-S set] \n"
"  [-D detect] [-O format] [-l] [-r] [filename | -p poly] "
#if HAVE_GRAPHICAL_DEBUGGER
          "[-x] "           
#endif
           "[-d] [-v] [infile]\n"
           "\n"
           "Options:\n"
           " -a alg      Select the algorithm used to solve the polynomial/secular equation:\n"
           "              u: Classic unisolve algorithm (Aberth iterations and dynamic precision)\n"
           "              s: Secular algorithm, using regeneration of increasingly better-conditioned\n"
           "                 secular equations with the same roots of the polynomial\n"
           " -b          Perform Aberth iterations in Jacobi-style instead of Gauss-Seidel\n"
	   " -c          Enable crude approximation mode. Fast but not always effective\n"
           " -G goal     Select the goal to reach. Possible values are:\n"
           "              a: Approximate the roots\n"
           "              i: Isolate the roots\n"
           "              c: Count the roots in the search set\n"
           " -o digits   Number of guaranteed digits of the roots\n"
           " -i digits   Digits of precision of the input coefficients\n"
           " -j n        Number of threads to spawn as workers\n"
           " -t type     Type can be 'f' for floating point or 'd' for DPE\n"
           " -S set      Restrict the search set for the roots \n"
           "             set can be one of:\n"
           "               u: upper half-plane { x | Im(x) > 0 } \n"
           "               d: lower half-plane { x | Im(x) < 0 } \n"
           "               l: left half-plane { x | Re(x) < 0 } \n"
           "               r: right half-plane { x | Re(x) > 0 } \n"
           "               i: inside the unit circle: { x | |x| < 1 } \n"
           "               o: outside the unit circle { x | |x| > 1 } \n"
	   "               R: real axis { x | Im(x) = 0 } \n"
	   "               I: imaginary axis { x | Re(x) = 0 }\n"
           " -D detect   Detect properties of the roots:\n"
           "               r: real roots\n"
           "               i: imaginary roots\n"
           "               b: both\n"
           " -O format   Select format for output:\n"
           "               f: full output\n"
           "               b: bare output\n"
           "               c: compact output\n"
           "               v: verbose output\n"
           "               g: gnuplot-ready output\n"
           "               gf: gnuplot-full mode, can be piped to gnuplot and display error bars. \n"
           "                   For example:\n"
           "                     %s -as -Ogf myfile.pol | gnuplot \n"
           "               gp: The same as gf but only with points (suitable for high degree polynomials)\n"
           " -l filename Set filename as the output for the log, instead of the tty. Use this option with\n"
           "             -d[domains] to activate the desired debug domains. \n"
#if HAVE_GRAPHICAL_DEBUGGER           
           " -x          Enable graphic visualization of convergence\n"
#endif
#ifndef DISABLE_DEBUG
          " -d[domains] Activate debug on selected domains, that can be one of:\n"
           "               t: trace\n"
           "               a: approximation\n"
           "               c: cluster\n"
           "               i: improvement\n"
           "               w: timings\n"
           "               o: input/Output\n"
           "               m: memory management\n"
           "               f: function calls\n"
           "               p: debug stop condition and development of iteration packets\n"
           "               r: regeneration\n"
           "               Example: -dfi for function calls and improvement\n"
#endif
	   " -p poly     Solve the polynomial specified on the command line. \n"
           "               Example: %s -p \"x^4-6*x^9+6/7*x + 5\" \n"
	   " -r          Use a recursive strategy to dispose the initial approximations.\n"
	   "             This option is available only for monomial polynomials. \n"
	   "             Note: this option is considered experimental.\n"
	   " -s file     Read the starting approximations from the given file, instead\n"
	   "             of relying on the internal algorithm of MPSolve.\n"
	   "             The format for this file is the same of the *.res files foun in\n"
	   "             src/unisolve/*.res in the source distribution of MPSolve.\n"
           " -v          Print the version and exit\n"
           "\n",
           program, program, program);

  exit (EXIT_FAILURE);
}

void*
cleanup_context (mps_context * ctx, void * user_data)
{
  /* Check for errors */
  if (mps_context_has_errors (ctx))
    {
      mps_print_errors (ctx);
      return NULL;
    }

  /* Output the roots */
  mps_output (ctx);

#ifdef HAVE_GRAPHICAL_DEBUGGER
  if (logger_closed)
    gtk_main_quit ();
  else if (logger)
  {
    /* In the other case copy the approximation in the right place
     * so the logger can display them again. */
    int i = 0;
    mps_approximation ** approximations = malloc (sizeof (mps_approximation*) * ctx->n);

    for (i = 0; i < ctx->n; i++)
      approximations[i] = mps_approximation_copy (ctx, ctx->root[i]);

    mps_iteration_logger_set_roots (logger, approximations, ctx->n);
  }
#endif  

  /* Free used data */
  if (poly) 
    mps_polynomial_free (ctx, poly);
  mps_context_free (ctx);

  s = NULL;

  return NULL;
}

#ifdef HAVE_GRAPHICAL_DEBUGGER
static void on_iteration_logger_destroy (MpsIterationLogger * logger, GdkEvent * event, gpointer user_data)
{
  gtk_widget_hide (GTK_WIDGET (logger));
  logger_closed = true;

  if (s == NULL)
    gtk_main_quit ();
}
#endif

int
main (int argc, char **argv)
{
#if HAVE_GRAPHICAL_DEBUGGER
  mps_boolean graphic_debug = false;
#endif

  /* Create a new status */
  s = mps_context_new ();

  /* Associate the SIGUSR1 signal to the mps_dump () function */
#ifndef __WINDOWS
  signal (SIGUSR1, status);
#endif

  /* Leave this value to -1 if not precision has been enforced
   * by the user. Otherwise, override the polynomial input
   * precision. */
  long int input_precision = -1;
  char * inline_poly = NULL;

  FILE *infile;

  /* Parse options */
  mps_opt *opt;
  mps_phase phase = no_phase;

  /* This will be true if the user has explicitely selected the algorithm, 
     otherwise we will be using our own heuristic. */
  mps_boolean explicit_algorithm_selection = false;

  opt = NULL;
  while ((mps_getopts (&opt, &argc, &argv, MPSOLVE_GETOPT_STRING)))
    {
      switch (opt->optchar)
        {
        case 'l':
          {
            FILE* logstr = fopen (opt->optvalue, "w");
            if (!logstr)
              mps_error (s, "Cannot open selected log file.");
            mps_context_set_log_stream (s, logstr);
          }
          break;

        case 'v':

#ifdef HAVE_CONFIG_H
          printf (PACKAGE_STRING "\n");
#else
          printf ("MPSolve 3.2.1\n");
#endif

          mps_context_free (s);
          exit (EXIT_SUCCESS);

	case 'r':
	  mps_context_select_starting_strategy (s, MPS_STARTING_STRATEGY_RECURSIVE);
	  break;

        case 'O':
          /* Select the desired output format */
          if (!opt->optvalue)
            {
              mps_error (s, "An argument is needed for option 'O'");
              break;
            }
          
          switch (*opt->optvalue)
            {
            case 'f':
              mps_context_set_output_format (s, MPS_OUTPUT_FORMAT_FULL);
              break;
            case 'b':
              mps_context_set_output_format (s, MPS_OUTPUT_FORMAT_BARE);
              break;
            case 'g':
              mps_context_set_output_format (s, MPS_OUTPUT_FORMAT_GNUPLOT);
              if (*(opt->optvalue + 1) == 'f')
                {
                  mps_context_set_output_format (s, MPS_OUTPUT_FORMAT_GNUPLOT_FULL);
                  s->gnuplot_format = "xyerrorbars";
                }
              else if (*(opt->optvalue + 1) == 'p')
                {
                  mps_context_set_output_format (s, MPS_OUTPUT_FORMAT_GNUPLOT_FULL);
                  s->gnuplot_format = "points";
                }
              break;
            case 'v':
              mps_context_set_output_format (s, MPS_OUTPUT_FORMAT_VERBOSE);
              break;
            case 'c':
              mps_context_set_output_format (s, MPS_OUTPUT_FORMAT_COMPACT);
              break;
            default:
              mps_error (s, "The selected output format is not supported");
              break;
            }

          break;

            /* select search set */
          case 'S':
            switch (*opt->optvalue)
              {
              case 'a':
                mps_context_set_search_set (s, MPS_SEARCH_SET_COMPLEX_PLANE);
                break;
              case 'r':
                mps_context_set_search_set (s, MPS_SEARCH_SET_POSITIVE_REAL_PART);
                break;
              case 'l':
                mps_context_set_search_set (s, MPS_SEARCH_SET_NEGATIVE_REAL_PART);
                break;
              case 'u':
                mps_context_set_search_set (s, MPS_SEARCH_SET_POSITIVE_IMAG_PART);
                break;
              case 'd':
                mps_context_set_search_set (s, MPS_SEARCH_SET_NEGATIVE_IMAG_PART);
                break;
              case 'i':
                mps_context_set_search_set (s, MPS_SEARCH_SET_UNITARY_DISC);
                break;
              case 'o':
                mps_context_set_search_set (s, MPS_SEARCH_SET_UNITARY_DISC_COMPL);
                break;
              case 'R':
                mps_context_set_search_set (s, MPS_SEARCH_SET_REAL);
                break;
              case 'I':
                mps_context_set_search_set (s, MPS_SEARCH_SET_IMAG);
                break;
              case 'U':
                mps_context_set_search_set (s, MPS_SEARCH_SET_CUSTOM);
                break;
              default:
                mps_error (s, "Bad search set switch: ", opt->optvalue,
                           ", use a|r|l|u|d|i|o|R|I|U");
              }
            if (strlen (opt->optvalue) != 1)
              mps_error (s, "Bad set: ", opt->optvalue);
            break;

            /* select multiplicity */
          case 'M':
            switch (*opt->optvalue)
              {
              case '+':
                s->output_config->multiplicity = true;
                break;
              case '-':
                s->output_config->multiplicity = false;
                break;
              default:
                mps_error (s, "Bad multiplicity switch: ", opt->optvalue,
                           ", use +|-");
              }
            if (strlen (opt->optvalue) != 3)
              mps_error (s, "Bad multiplicity option: ", opt->optvalue);
            break;

            /* detection */
          case 'D':
            switch (*opt->optvalue)
              {
              case 'n':
                mps_context_set_root_properties (s, MPS_OUTPUT_PROPERTY_NONE);
                break;
              case 'r':
                mps_context_set_root_properties (s, MPS_OUTPUT_PROPERTY_REAL);
                break;
              case 'i':
                mps_context_set_root_properties (s, MPS_OUTPUT_PROPERTY_IMAGINARY);
                break;
              case 'b':
                mps_context_set_root_properties (s, MPS_OUTPUT_PROPERTY_REAL | MPS_OUTPUT_PROPERTY_IMAGINARY);
                break;
              default:
                mps_error (s, "Bad detection switch: ", opt->optvalue,
                           ", use n|r|i|b");
              }
            if (strlen (opt->optvalue) != 1)
              mps_error (s, "Bad detection option: ", opt->optvalue);
            break;

            /* I/O streams */
          case 'R':
            s->rtstr = fopen (opt->optvalue, "r");
            if (s->rtstr == NULL)
              mps_error (s, "Cannot open roots file: ", opt->optvalue);
            s->resume = true;
            break;

            /* Additional checks */
          case 'C':
            switch (*opt->optvalue)
              {
              case 'R':
                s->chkrad = true;
                break;
              case 'r':
                s->chkrad = false;
                break;
              default:
                mps_error (s, "Bad check switch: ", opt->optvalue,
                           ", use R|r");
              }
            if (strlen (opt->optvalue) != 1)
              mps_error (s, "Bad check option: ", opt->optvalue);
            break;

            
        case 'a':
	  explicit_algorithm_selection = true;
          switch (*opt->optvalue)
            {
            case 'u':
              mps_context_select_algorithm (s, MPS_ALGORITHM_STANDARD_MPSOLVE);
              break;
            case 's':
              mps_context_select_algorithm (s, MPS_ALGORITHM_SECULAR_GA);
              break;
            default:
              mps_error (s, "The selected algorithm is not supported");
              break;
            }
          break;
        case 'b':
          mps_context_set_jacobi_iterations (s, true);
          break;
	case 'c':
	  mps_context_set_crude_approximation_mode (s, true);
	  break;
        case 'o':
          mps_context_set_output_prec (s, (atoi (opt->optvalue)) * LOG2_10 + 1);
          break;
        case 'i':
          input_precision = atoi (opt->optvalue) * LOG2_10;
          break;
        case 'G':
          switch (*opt->optvalue)
            {
            case 'a':
              mps_context_set_output_goal (s, MPS_OUTPUT_GOAL_APPROXIMATE);
              break;
            case 'i':
              mps_context_set_output_goal (s, MPS_OUTPUT_GOAL_ISOLATE);
              break;
            case 'c':
              mps_context_set_output_goal (s, MPS_OUTPUT_GOAL_COUNT);
              break;
            default:
              mps_error (s, "The selected goal does not exists");
              break;
            }
          break;

        case 'p':
          inline_poly = strdup (opt->optvalue);
          break;

#ifndef DISABLE_DEBUG
        case 'd':
          mps_context_add_debug_domain (s, MPS_DEBUG_INFO);
          
          if (!opt->optvalue)
            break;

          /* If debugging was enabled, parse debug_level */
          while (*opt->optvalue)
            {         
              char output[255];
              switch (*opt->optvalue++)
                {
                case 't':
                  mps_context_add_debug_domain (s, MPS_DEBUG_TRACE);
                  break;
                case 'a':
                  mps_context_add_debug_domain (s, MPS_DEBUG_APPROXIMATIONS);
                  break;
                case 'c':
                  mps_context_add_debug_domain (s, MPS_DEBUG_CLUSTER);
                  break;
                case 'i':
                  mps_context_add_debug_domain (s, MPS_DEBUG_IMPROVEMENT);
                  break;
                case 'w':
                  mps_context_add_debug_domain (s, MPS_DEBUG_TIMINGS);
                  break;
                case 'o':
                  mps_context_add_debug_domain (s, MPS_DEBUG_IO);
                  break;
                case 'm':
                  mps_context_add_debug_domain (s, MPS_DEBUG_MEMORY);
                  break;
                case 'f':
                  mps_context_add_debug_domain (s, MPS_DEBUG_FUNCTION_CALLS);
                  break;
                case 'p':
                  mps_context_add_debug_domain (s, MPS_DEBUG_PACKETS);
                  break;
                case 'r':
                  mps_context_add_debug_domain (s, MPS_DEBUG_REGENERATION);
                  break;
                default:
                  sprintf (output, "Unrecognized debug option: %c", *(opt->optvalue - 1));
                  mps_error (s, output);
                  break;
                }
            }
          break;
#endif

#if HAVE_GRAPHICAL_DEBUGGER
        case 'x':
          graphic_debug = true;
          break;
#endif          

	case 's':
	  if (opt->optvalue == NULL)
	    {
	      mps_error (s, "You need to specify a valid file with -s");
	      break;
	    }
	  else
	    {
	      s->rtstr = fopen (opt->optvalue, "r");
	      if (! s->rtstr)
		{
		  mps_error (s, "Cannot open the given file: %s", opt->optvalue);
		  break;
		}
	      else
		{
		  mps_context_select_starting_strategy (s, MPS_STARTING_STRATEGY_FILE);
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
          mps_context_set_n_threads (s, atoi (opt->optvalue));
          break;
        default:
          usage (s, argv[0]);
          break;
        }
    }

#if HAVE_GRAPHICAL_DEBUGGER
  if (graphic_debug)
  {
    /* Init GTK only if graphic debug is requested. In all other case is unnecessary
     * and will waste computational time that is likely to bias benchmarks. */
    gtk_init (&argc, &argv);

    logger = mps_iteration_logger_new ();
    mps_iteration_logger_set_mps_context (logger, s);

    g_signal_connect (GTK_WIDGET (logger), "delete_event", 
       G_CALLBACK (on_iteration_logger_destroy), NULL);

    gtk_widget_show_all (GTK_WIDGET (logger));
  }
#endif

  if (mps_context_has_errors (s))
    {
      mps_print_errors (s);
      return EXIT_FAILURE;
    }

  if (argc > 2)
    usage (s, argv[0]);

  /* If no file is provided use standard input */
  if (inline_poly)
    {
      mps_abstract_input_stream * stream = 
       (mps_abstract_input_stream *) mps_memory_file_stream_new (inline_poly);
      poly = mps_monomial_yacc_parser (s, stream);

      if (mps_context_has_errors (s))
       {
         mps_print_errors (s);
         return EXIT_FAILURE;
       }

      free (inline_poly);
    }
  else
    {
      if (argc == 1)
        infile = stdin;
      else
        infile = fopen (argv[1], "r");

      if (!infile)
        {
          mps_error (s, "Cannot open input file for read, aborting.");
          mps_print_errors (s);
          return EXIT_FAILURE;
        }

      /* Parse the input stream and if a polynomial is given as output, 
       * allocate also a secular equation to be used in regeneration */
       poly = mps_parse_stream (s, infile);
     }

  if (!poly)
    {
      mps_error (s, "Error while parsing the polynomial, aborting.");
      mps_print_errors (s);
      return EXIT_FAILURE;
    }
  else
    mps_context_set_input_poly (s, poly);

  /* Override input precision if needed */
  if (input_precision >= 0)
    mps_polynomial_set_input_prec (s, poly, input_precision);

  /* Perform some heuristic for the algorithm selection, but only if the user
   * hasn't explicitely selected one. */
  if (! explicit_algorithm_selection)
    {
      mps_context_select_algorithm (s, (MPS_IS_MONOMIAL_POLY (poly) && 
					MPS_DENSITY_IS_SPARSE (poly->density)) ? 
				    MPS_ALGORITHM_STANDARD_MPSOLVE : MPS_ALGORITHM_SECULAR_GA );
    }

  /* Close the file if it's not stdin */
  if (argc == 2 && ! inline_poly)
    fclose (infile);

  /* Select the starting phase according to user input */
  mps_context_set_starting_phase (s, phase);

  /* Solve the polynomial */
  #if HAVE_GRAPHICAL_DEBUGGER
  if (graphic_debug)
  {
    mps_mpsolve_async (s, cleanup_context, NULL);
    gtk_main ();
  }
  else
  {
  #endif
    mps_mpsolve (s);
    cleanup_context (s, NULL);
  #if HAVE_GRAPHICAL_DEBUGGER
  }
  #endif

}
