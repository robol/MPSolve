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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_GTK
 #include <gtk/gtk.h>
 #include <iteration-logger.h>

MpsIterationLogger * logger = NULL;
static mps_boolean logger_closed = false;
#endif

mps_context * s = NULL;

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
           "%s [-l filename] [-dv] [-S set] [-D detection] [-O format] [-G goal] [-t type] [-a alg] [-j threads] [-i digits] [-o digits] [infile]\n"
           "\n"
           "Options:\n"
           " -l filename Set filename as the output for the log, instead of the tty. Use this option with\n"
           "             -d[domains] to activate the desired debug domains. \n"
           " -d[domains] Activate debug on selected domains, that can be one of:\n"
           "               t: Trace\n"
           "               a: Approximation\n"
           "               c: Cluster\n"
           "               i: Improvement\n"
           "               w: Timings\n"
           "               o: Input/Output\n"
           "               m: Memory management\n"
           "               f: Function calls\n"
           "               p: Debug stop condition and development of iteration packets\n"
           "               r: Regeneration\n"
           "               Example: -dfi for function calls and improvement\n"
#if HAVE_GTK           
           " -x          Enable graphic visualization of convergence\n"
#endif           
           " -S set      If specified, restrict the search set for the roots to set. \n"
           "             set can be one of:\n"
           "               u: Semiplane { x | Im(x) > 0 } \n"
           "               d: Semiplane { x | Im(x) < 0 } \n"
           "               l: Semiplane { x | Re(x) < 0 } \n"
           "               r: Semiplane { x | Re(x) > 0 } \n"
           "               i: Inside the unit circle: { x | |x| < 1 } \n"
           "               o: Outside the unit circle { x | |x| > 1 } \n"
           " -D detect   Detect reality or imaginarity of the roots:\n"
           "               r: Detect real roots\n"
           "               i: Detect imaginary roots\n"
           "               b: Detect both\n"
           " -O format   Select format for output:\n"
           "               f: Full output\n"
           "               b: Bare output\n"
           "               c: Compact output\n"
           "               v: Verbose output\n"
           "               g: Gnuplot-ready output\n"
           "               gf: Gnuplot-full mode, can be piped to gnuplot and display error bars. \n"
           "                   For example:\n"
           "                     %s -as -Ogf myfile.pol | gnuplot \n"
           "               gp: The same as gf but only with points (suitable for high degree polynomials)\n"
           " -G goal     Select the goal to reach. Possible values are:\n"
           "              a: Approximate the roots\n"
           "              i: Isolate the roots\n"
           "              c: Count the roots in the search set\n"
           " -a alg      Select the algorithm used to solve the polynomial/secular equation:\n"
           "              u: Classic unisolve algorithm (Aberth iterations and dynamic precision)\n"
           "              s: Secular algorithm, using regeneration of increasingly better-conditioned\n"
           "                 secular equation with the same roots of the polynomial\n"
           " -t type     Type can be 'f' for floating point or 'd' for DPE\n"
           " -j n        Number of threads to spawn as workers\n"
           " -o digits   Exact digits of the roots given as output.\n"
           " -i digits   Digits of precision of the input coefficients\n"
           " -v          Print the version and exit\n"
           "\n",
           program, program);

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

  /* Free used data */
  mps_context_free (ctx);

#ifdef HAVE_GTK
  if (logger_closed)
    gtk_main_quit ();
#endif  

  s = NULL;

  return NULL;
}

#ifdef HAVE_GTK
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
  /* Create a new status */
  s = mps_context_new ();

  /* Associate the SIGUSR1 signal to the mps_dump () function */
#ifndef __WINDOWS
  signal (SIGUSR1, status);
#endif

  mps_context_set_input_prec (s, 0);

  FILE *infile;

  /* Parse options */
  mps_opt *opt;
  mps_phase phase = no_phase;

  mps_boolean graphic_debug = false;

#ifdef HAVE_GTK
    gtk_init (&argc, &argv);
#endif  

  opt = NULL;
  while ((mps_getopts (&opt, &argc, &argv, "a:G:D:d::xt:o:O:j:S:O:i:vl:")))
    {
      switch (opt->optchar)
        {
        case 'l':
          {
            FILE* logstr = fopen (opt->optvalue, "w");
            if (!logstr)
              mps_error (s, 1, "Cannot open selected log file.");
            mps_context_set_log_stream (s, logstr);
          }
          break;

        case 'v':

#ifdef HAVE_CONFIG_H
          printf ("MPSolve " VERSION "\n");
#else
          printf ("MPSolve 3.0\n");
#endif

          mps_context_free (s);
          exit (EXIT_SUCCESS);

        case 'O':
          /* Select the desired output format */
          if (!opt->optvalue)
            mps_error (s, 1, "An argument is needed for option 'O'");
          
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
              mps_error (s, 1, "The selected output format is not supported");
              break;
            }

          break;

            /* select search set */
          case 'S':
            switch (*opt->optvalue)
              {
              case 'a':
                s->output_config->search_set = MPS_SEARCH_SET_COMPLEX_PLANE;
                break;
              case 'r':
                s->output_config->search_set = MPS_SEARCH_SET_POSITIVE_REAL_PART;
                break;
              case 'l':
                s->output_config->search_set = MPS_SEARCH_SET_NEGATIVE_REAL_PART;
                break;
              case 'u':
                s->output_config->search_set = MPS_SEARCH_SET_POSITIVE_IMAG_PART;
                break;
              case 'd':
                s->output_config->search_set = MPS_SEARCH_SET_NEGATIVE_IMAG_PART;
                break;
              case 'i':
                s->output_config->search_set = MPS_SEARCH_SET_UNITARY_DISC;
                break;
              case 'o':
                s->output_config->search_set = MPS_SEARCH_SET_UNITARY_DISC_COMPL;
                break;
              case 'R':
                s->output_config->search_set = MPS_SEARCH_SET_REAL;
                break;
              case 'I':
                s->output_config->search_set = MPS_SEARCH_SET_IMAG;
                break;
              case 'U':
                s->output_config->search_set = MPS_SEARCH_SET_CUSTOM;
                break;
              default:
                mps_error (s, 3, "Bad search set switch: ", opt->optvalue,
                           ", use a|r|l|u|d|i|o|R|I|U");
              }
            if (strlen (opt->optvalue) != 1)
              mps_error (s, 2, "Bad set: ", opt->optvalue);
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
                mps_error (s, 3, "Bad multiplicity switch: ", opt->optvalue,
                           ", use +|-");
              }
            if (strlen (opt->optvalue) != 3)
              mps_error (s, 2, "Bad multiplicity option: ", opt->optvalue);
            break;

            /* detection */
          case 'D':
            switch (*opt->optvalue)
              {
              case 'n':
                s->output_config->root_properties = MPS_OUTPUT_PROPERTY_NONE;
                break;
              case 'r':
                s->output_config->root_properties = MPS_OUTPUT_PROPERTY_REAL;
                break;
              case 'i':
                s->output_config->root_properties = MPS_OUTPUT_PROPERTY_IMAGINARY;
                break;
              case 'b':
                s->output_config->root_properties = MPS_OUTPUT_PROPERTY_REAL | 
                  MPS_OUTPUT_PROPERTY_IMAGINARY;
                break;
              default:
                mps_error (s, 3, "Bad detection switch: ", opt->optvalue,
                           ", use n|r|i|b");
              }
            if (strlen (opt->optvalue) != 1)
              mps_error (s, 2, "Bad detection option: ", opt->optvalue);
            break;

            /* I/O streams */
          case 'R':
            s->rtstr = fopen (opt->optvalue, "r");
            if (s->rtstr == NULL)
              mps_error (s, 2, "Cannot open roots file: ", opt->optvalue);
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
                mps_error (s, 3, "Bad check switch: ", opt->optvalue,
                           ", use R|r");
              }
            if (strlen (opt->optvalue) != 1)
              mps_error (s, 2, "Bad check option: ", opt->optvalue);
            break;

            
        case 'a':
          switch (*opt->optvalue)
            {
            case 'u':
              mps_context_select_algorithm (s, MPS_ALGORITHM_STANDARD_MPSOLVE);
              break;
            case 's':
              mps_context_select_algorithm (s, MPS_ALGORITHM_SECULAR_GA);
              break;
            default:
              mps_error (s, 1, "The selected algorithm is not supported");
              break;
            }
          break;
        case 'o':
          mps_context_set_output_prec (s, (atoi (opt->optvalue)) * LOG2_10 + 1);
          break;
        case 'i':
          mps_context_set_input_prec (s, (atoi (opt->optvalue)));
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
              mps_error (s, 1, "The selected goal does not exists");
              break;
            }
          break;
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
                  mps_error (s, 1, output);
                  break;
                }
            }
          break;

        case 'x':
          graphic_debug = true;
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
          s->n_threads = atoi (opt->optvalue);
          break;
        default:
          usage (s, argv[0]);
          break;
        }
    }

#ifdef HAVE_GTK
  if (graphic_debug)
  {
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
  if (argc == 1)
    infile = stdin;
  else
    infile = fopen (argv[1], "r");

  if (!infile)
    {
      mps_error (s, 1, "Cannot open input file for read, aborting.");
      mps_print_errors (s);
      return EXIT_FAILURE;
    }

  /* Parse the input stream and if a polynomial is given as output, 
   * allocate also a secular equation to be used in regeneration */
  mps_parse_stream (s, infile);

  /* Close the file if it's not stdin */
  if (argc == 2)
    fclose (infile);

  /* Select the starting phase according to user input */
  mps_context_set_starting_phase (s, phase);

  /* Solve the polynomial */
  #if HAVE_GTK
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
  #if HAVE_GTK
  }
  #endif

}
