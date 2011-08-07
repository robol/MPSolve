/***********************************************************
 **       Multiprecision Polynomial Solver (MPSolve)       **
 **                 Version 2.2, May 2001                  **
 **                                                        **
 **                      Written by                        **
 **       Dario Andrea Bini and Giuseppe Fiorentino        **
 **       (bini@dm.unipi.it)  (fiorent@dm.unipi.it)        **
 **                                                        **
 ** (C) 2001, Dipartimento di Matematica, FRISCO LTR 21024 **
 ***********************************************************/

/**
 * @file
 * @brief Implementation of option parsing
 *
 */

#include <string.h>
#include <errno.h>
#include <stdlib.h>
#include <float.h>
#include <mps/core.h>

/* local definition */
void mps_print_help (mps_status * s);

/**
 * @brief Parse command line options in a similar way of getopts.
 *
 * @param op_ptr A pointer to a previous used mps_opts that can
 * be re-used, or NULL. If this is NULL it will be allocated.
 * When this function return false the opt_ptr is automatically
 * freed.
 * @param argc_ptr The address in memory of the argc variable
 * obtained by the operating system.
 * @param argv_ptr The address in memory of the argv variable
 * obtained by the operating system.
 * @param opt_format The option that <code>mps_getopts()</code>
 * has to look for. It is passed as a string like
 * <code>"ab:n" where the <code>":"</code> means that the preceding
 * character needs an argument, while single characters not followed
 * by <code>":"</code> don't expect one. In this case a correct call
 * to the program would be
 * @code
 *   ./myprogram -a -n -b 7
 * @endcode
 * or even
 * @code
 *   ./myprogam -a
 * @endcode
 *
 * The option are parsed from argc_ptr and argv_ptr and these
 * are modified taking options away (as the program has been
 * called wihtout these options) so after option parsing the program
 * can still parse other options such as filenames or similar, without
 * worrying about options switch.
 */
mps_boolean
mps_getopts (mps_opt ** opt_ptr, int *argc_ptr, char ***argv_ptr,
	     const char *opt_format)
{
  char **argv = *argv_ptr;
  int argc = *argc_ptr;
  char *tmp;
  mps_opt *opt;
  int i, l = strlen (opt_format), steps;

  if ((*opt_ptr) == NULL)
    {
      (*opt_ptr) = (mps_opt *) malloc (sizeof (mps_opt));
    }

  opt = *opt_ptr;

  /* Check if there are other arguments to parse */
  if (!argc)
    {
      free (opt);
      return false;
    }

  if (argc == 1)
    {
      free (opt);
      return false;
    }

  /* Scan for right offset */
  steps = 0;
  while (argv[1][0] != '-')
    {
      /* If we get here then argv[1] is a parameter string
       * that must be preserved, so we should put it on the end
       * of argv. */
      tmp = argv[1];
      for (i = 1; i < argc - 1; i++)
	{
	  (*argv_ptr)[i] = argv[i + 1];
	}
      (*argv_ptr)[argc - 1] = tmp;
      steps++;

      /* Check if argc permutation was performed */
      if (steps == argc - 1)
	{
	  free (opt);
	  return false;
	}
    }


  /* Search argv[0][1] in opt_format */
  for (i = 0; i < l; i++)
    {
      if (opt_format[i] == ':')
	continue;
      else if (argv[1][1] == opt_format[i])
	{
	  opt->optchar = opt_format[i];

	  /* If there is no argument we are done */
	  if (i == l || opt_format[i + 1] != ':')
	    {
	      opt->optvalue = NULL;
	      (*argc_ptr)--;
	      (*argv_ptr)[1] = argv[0];
	      (*argv_ptr)++;
	      return true;
	    }

	  /* If the string is not terminated than we should
	   * expect to find the parameter attached to it */
	  if (argv[1][2] != '\0')
	    {
	      if (argv[1][2] == '=')
		{
		  opt->optvalue = argv[1] + 3;
		  (*argc_ptr)--;
		  (*argv_ptr)[1] = argv[0];
		  (*argv_ptr)++;
		  return true;
		}
	      else
		{
		  opt->optvalue = argv[1] + 2;
		  (*argc_ptr)--;
		  (*argv_ptr)[1] = argv[0];
		  (*argv_ptr)++;
		  return true;
		}
	    }

	  /* If the parameter should be in argv[1] but is not
	   * there return an error */
	  if (argc == 2)
	    {
	      opt->optvalue = NULL;
	      (*argc_ptr)--;
	      (*argv_ptr)[1] = argv[0];
	      (*argv_ptr)++;
	      return true;
	    }

	  /* Otherwise we can set the argument in the struct */
	  opt->optvalue = argv[2];
	  (*argc_ptr) -= 2;
	  (*argv_ptr)[2] = argv[0];
	  (*argv_ptr) += 2;
	  return true;
	}
    }

  /* Fallback case, we haven't recognized the input */
  opt->optchar = '?';
  (*argc_ptr)--;
  (*argv_ptr)[1] = argv[0];
  (*argv_ptr)++;

  return true;
}

/**
 * This function parse command lines and sets global variables
 * defined in mps.h in an appropriate way.
 *
 * It is implemented in mps_opts.c
 *
 * @brief Parse options from command line.
 *
 * @param argc
 *   Argoment counter as obtained from the main()
 *   function.
 * @param argv
 *   Argoment values, i.e. array of char* as obtained
 *   in the main() function.
 */
void
mps_parse_opts (mps_status * s, int argc, char *argv[])
{
  unsigned int seed = 0;
  int i;

  /* set default flags */
  s->DOLOG = false;
  s->DOWARN = true;
  s->DOSORT = true;

  /* set default streams */
  s->instr = stdin;
  s->outstr = stdout;
  s->logstr = stderr;
  s->rtstr = stdin;

  /* set default values */
  s->prec_in = -1;		/* if != -1 ignore precision from file */
  s->prec_out = 2 * DBL_DIG;	/* default output precision            */
  s->random_seed = false;

  /* parse options */
  for (i = 1; i < argc; i++)
    {
      if (argv[i][0] != '-')
	{
	  s->instr = fopen (argv[i], "r");
	  goto finalcheck;
	  /* no options allowed after filename  */
	}
      else
	/* all options start with - */
	switch (argv[i][1])
	  {
	  case 'i':
	    s->prec_in = atol (argv[i] + 2);
	    if (s->prec_in < 0 || errno)
	      mps_error (s, 2, "Wrong input precision: ", argv[i] + 2);
	    s->prec_in = (long) (s->prec_in * LOG2_10);
	    break;

	  case 'r':
	    s->random_seed = true;
	    break;

	  case 'o':
	    s->prec_out = atol (argv[i] + 2);
	    if (s->prec_out <= 0 || errno)
	      mps_error (s, 2, "Wrong output precision: ", argv[i] + 2);
	    s->prec_out = (long) (s->prec_out * LOG2_10);
	    break;

	    /* goal */
	  case 'G':
	    switch (argv[i][2])
	      {
	      case 'a':
		s->goal[0] = 'a';
		break;
	      case 'c':
		s->goal[0] = 'c';
		break;
	      case 'i':
		s->goal[0] = 'i';
		break;
	      default:
		mps_error (s, 3, "Bad goal switch: ", argv[i] + 2,
			   ", use a|c|i");
	      }
	    if (strlen (argv[i]) != 3)
	      mps_error (s, 2, "Bad goal: ", argv[i]);
	    break;

	    /* select search set */
	  case 'S':
	    switch (argv[i][2])
	      {
	      case 'a':
		s->goal[1] = 'a';
		break;
	      case 'r':
		s->goal[1] = 'r';
		break;
	      case 'l':
		s->goal[1] = 'l';
		break;
	      case 'u':
		s->goal[1] = 'u';
		break;
	      case 'd':
		s->goal[1] = 'd';
		break;
	      case 'i':
		s->goal[1] = 'i';
		break;
	      case 'o':
		s->goal[1] = 'o';
		break;
	      case 'R':
		s->goal[1] = 'R';
		break;
	      case 'I':
		s->goal[1] = 'I';
		break;
	      case 'U':
		s->goal[1] = 'U';
		break;
	      default:
		mps_error (s, 3, "Bad search set switch: ", argv[i] + 2,
			   ", use a|r|l|u|d|i|o|R|I|U");
	      }
	    if (strlen (argv[i]) != 3)
	      mps_error (s, 2, "Bad set: ", argv[i]);
	    break;

	    /* select multiplicity */
	  case 'M':
	    switch (argv[i][2])
	      {
	      case '+':
		s->goal[2] = 'm';
		break;
	      case '-':
		s->goal[2] = 'n';
		break;
	      default:
		mps_error (s, 3, "Bad multiplicity switch: ", argv[i] + 2,
			   ", use +|-");
	      }
	    if (strlen (argv[i]) != 3)
	      mps_error (s, 2, "Bad multiplicity option: ", argv[i]);
	    break;

	    /* detection */
	  case 'D':
	    switch (argv[i][2])
	      {
	      case 'n':
		s->goal[3] = 'n';
		break;
	      case 'r':
		s->goal[3] = 'r';
		break;
	      case 'i':
		s->goal[3] = 'i';
		break;
	      case 'b':
		s->goal[3] = 'b';
		break;
	      default:
		mps_error (s, 3, "Bad detection switch: ", argv[i] + 2,
			   ", use n|r|i|b");
	      }
	    if (strlen (argv[i]) != 3)
	      mps_error (s, 2, "Bad detection option: ", argv[i]);
	    break;

	    /* I/O streams */
	  case 'R':
	    s->rtstr = fopen (argv[i] + 2, "r");
	    if (s->rtstr == NULL)
	      mps_error (s, 2, "Cannot open roots file: ", argv[i] + 2);
	    s->resume = true;
	    break;

	    /* Additional checks */
	  case 'C':
	    switch (argv[i][2])
	      {
	      case 'R':
		s->chkrad = true;
		break;
	      case 'r':
		s->chkrad = false;
		break;
	      default:
		mps_error (s, 3, "Bad check switch: ", argv[i] + 2,
			   ", use R|r");
	      }
	    if (strlen (argv[i]) != 3)
	      mps_error (s, 2, "Bad check option: ", argv[i]);
	    break;

	    /* output format */
	  case 'O':
	    switch (argv[i][2])
	      {
	      case 'b':
		s->goal[4] = 'b';
		break;
	      case 'g':
		s->goal[4] = 'g';
		break;
	      case 'c':
		s->goal[4] = 'c';
		break;
	      case 'v':
		s->goal[4] = 'v';
		break;
	      case 'f':
		s->goal[4] = 'f';
		break;
	      default:
		mps_error (s, 3, "Bad output format switch: ", argv[i] + 2,
			   ", use b|c|f|g|v");
	      }
	    if (strlen (argv[i]) != 3)
	      mps_error (s, 2, "Bad output option: ", argv[i]);
	    break;

	    /* iteration limits */
	  case 'L':
	    switch (argv[i][2])
	      {
	      case 'p':
		s->max_pack = atoi (argv[i] + 3);
		if (s->max_pack <= 0 || errno)
		  mps_error (s, 2, "Invalid number of packet iterations: ",
			     argv[i] + 3);
		break;
	      case 'i':
		s->max_it = atoi (argv[i] + 3);
		if (s->max_it <= 0 || errno)
		  mps_error (s, 2, "Invalid number of iterations: ",
			     argv[i] + 3);
		break;
	      default:
		mps_error (s, 3, "Bad limit switch: ", argv[i] + 2,
			   ", use (p|i)<num>");
	      }
	    break;

	    /* randomization */
	  case 'H':
	    seed = (unsigned int) atoi (argv[i] + 2);
	    if (seed == 0 || errno)
	      mps_error (s, 2, "Wrong random seed: ", argv[i] + 2);
	    break;

	    /* debug options */
	  case 'd':
	    s->DOLOG = true;
	    if (strlen (argv[i]) == 3 && argv[i][2] == '1')
	      s->logstr = s->outstr;
	    else if (strlen (argv[i]) != 2)
	      mps_error (s, 3, "Bad debug option: ", argv[i], ", use 1");
	    break;

	  case 'j':
	    if (strlen (argv[i]) < 2)
	      mps_error (s, 1,
			 "You should provide a number of threads to be spawned.");
	    else
	      {
		s->n_threads = atoi (argv[i] + 2);
		if (s->n_threads < 1 || errno)
		  mps_error (s, 1,
			     "The number of threads must be a positive integer");
	      }
	    break;

	    /* warning options */
	  case 'w':
	    if (strlen (argv[i]) == 3 && argv[i][2] == '+')
	      s->DOWARN = true;
	    else if (strlen (argv[i]) == 3 && argv[i][2] == '-')
	      s->DOWARN = false;
	    else
	      mps_error (s, 3, "Bad warning option: ", argv[i], ", use +|-");
	    break;

	    /* help and default case */
	  case 'h':
	    mps_print_help (s);
	    exit (EXIT_SUCCESS);
	    break;

	    /* handle illegal options */
	  default:
	    mps_error (s, 3, "Bad option: ", argv[i],
		       ", type 'unisolve -h' for help");
	  }

    }

finalcheck:

  /* If the goal is approximate or count then remove the multiplicity option */
  if (s->goal[0] != 'i')
    s->goal[2] = 'n';

  /* check I/O streams */
  if (s->instr == NULL)
    mps_error (s, 1, "Cannot open input file");
  if (s->outstr == NULL)
    mps_error (s, 1, "Cannot open output file");
  if (s->DOLOG && s->logstr == NULL)
    mps_error (s, 1, "Cannot open log file");

  /* randomize */
  if (s->random_seed)
    randomize (seed);
}

/**
 * @brief Print usage informations
 */
void
mps_print_help (mps_status * s)
{
  fprintf (s->outstr, "USAGE: unisolve {<options>} {input_file}\n");
  fprintf (s->outstr, " OPTIONS (defaults in square brackets):\n");
  fprintf (s->outstr, " -in\tn = input precision, in decimal digits [0]\n");
  fprintf (s->outstr, " -on\tn = output precision, in decimal digits [%ld]\n",
	   s->prec_out);
  fprintf (s->outstr, " -Gc\tc = Goal: (a)pproximate, [i]solate, (c)ount\n");
  fprintf (s->outstr, " -Sc\tc = Search set: [a]ll, (r)ight, (l)eft, (u)p, "
	   "(d)own complex plane\n");
  fprintf (s->outstr, "   \t                (i)n, (o)ut unitary circle\n");
  fprintf (s->outstr, " -Mc\tc = Multiplicity check: (+) on, [-] off\n");
  fprintf (s->outstr,
	   " -Dc\tc = Detect: [n]one, (r)eal/(i)maginary/(b)oth\n");
  fprintf (s->outstr,
	   " -Oc\tc = Output format: (b)are, (g)nuplot, [c]ompact, "
	   "(v)erbose, (f)ull\n");
  fprintf (s->outstr,
	   " -Hn\tn = random seed, taken from /dev/random if exists\n");
  fprintf (s->outstr, " -d\tprint debug information to standard error\n");
  fprintf (s->outstr, " -d1\tprint debug information to standard output\n");
  fprintf (s->outstr, " -wc\tc = print warning messages: [+] on, (-) off\n");
  fprintf (s->outstr, " -Lpn\tn = maximum number of packet iterations [%d]\n",
	   s->max_pack);
  fprintf (s->outstr, " -Lin\tn = maximum number of global iterations [%d]\n",
	   s->max_it);
  fprintf (s->outstr, " -r\tuse random seed for starting points\n");
  fprintf (s->outstr, " -j\tn = number of threads to spawn\n");

  /* undocumented switches
     C
     R
   */
}
