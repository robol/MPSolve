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
#include "mps.h"

/* local definition */
void mps_print_help(mps_status* s);

/***********************************************************
 *                 PARSE_OPTS                              *
 ***********************************************************
 * check command line options and streams                  *
 **********************************************************/
void
mps_parse_opts(mps_status* s, int argc, char *argv[])
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
  s->prec_in = -1;                 /* if != -1 ignore precision from file */
  s->prec_out = 2 * DBL_DIG;       /* default output precision            */
  s->random_seed = false;

  /* parse options */
  for (i = 1; i < argc; i++)
    if (argv[i][0] != '-') {
      s->instr = fopen(argv[i], "r");
      goto finalcheck;  /* no options allowed after filename  */
    } else      /* all options start with - */
      switch (argv[i][1]) {
      case 'i':
        s->prec_in = atol(argv[i] + 2);
        if (s->prec_in <= 0 || errno)
          mps_error(s, 2, "Wrong input precision: ", argv[i]+2);
        s->prec_in = (long) (s->prec_in * LOG2_10);
        break;

      case 'r':
	s->random_seed = true;
	break;

      case 'o':
        s->prec_out = atol(argv[i] + 2);
        if (s->prec_out <= 0 || errno)
          mps_error(s, 2, "Wrong output precision: ", argv[i]+2);
        s->prec_out = (long) (s->prec_out * LOG2_10);
        break;

      /* goal */
      case 'G':
        switch (argv[i][2]) {
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
          mps_error(s, 3, "Bad goal switch: ", argv[i]+2, ", use a|c|i");
        }
        if (strlen(argv[i]) != 3)
          mps_error(s, 2, "Bad goal: ", argv[i]);
        break;

      /* select search set */
      case 'S':
        switch (argv[i][2]) {
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
          mps_error(s, 3, "Bad search set switch: ", argv[i]+2, ", use a|r|l|u|d|i|o|R|I|U");
        }
        if (strlen(argv[i]) != 3)
          mps_error(s, 2, "Bad set: ", argv[i]);
        break;

      /* select multiplicity */
      case 'M':
        switch (argv[i][2]) {
        case '+':
          s->goal[2] = 'm';
          break;
        case '-':
          s->goal[2] = 'n';
          break;
        default:
          mps_error(s, 3, "Bad multiplicity switch: ", argv[i]+2, ", use +|-");
        }
        if (strlen(argv[i]) != 3)
          mps_error(s, 2, "Bad multiplicity option: ", argv[i]);
        break;

      /* detection */
      case 'D':
        switch (argv[i][2]) {
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
          mps_error(s, 3, "Bad detection switch: ", argv[i]+2, ", use n|r|i|b");
        }
        if (strlen(argv[i]) != 3)
          mps_error(s, 2, "Bad detection option: ", argv[i]);
        break;

      /* I/O streams */
      case 'R':
        s->rtstr = fopen(argv[i] + 2, "r");
        if (s->rtstr == NULL)
          mps_error(s, 2, "Cannot open roots file: ", argv[i]+2);
        s->resume = true;
        break;

      /* Additional checks */
      case 'C':
        switch (argv[i][2]) {
        case 'R':
          s->chkrad = true;
          break;
        case 'r':
          s->chkrad = false;
          break;
        default:
          mps_error(s, 3, "Bad check switch: ", argv[i]+2, ", use R|r");
        }
        if (strlen(argv[i]) != 3)
          mps_error(s, 2, "Bad check option: ", argv[i]);
        break;

      /* output format */
      case 'O':
        switch (argv[i][2]) {
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
          mps_error(s, 3, "Bad output format switch: ", argv[i]+2, ", use b|c|f|g|v");
        }
        if (strlen(argv[i]) != 3)
          mps_error(s, 2, "Bad output option: ", argv[i]);
        break;

      /* iteration limits */ 
      case 'L':
        switch (argv[i][2]) {
        case 'p':
          s->max_pack = atoi(argv[i] + 3);
          if (s->max_pack <= 0 || errno)
            mps_error(s, 2, "Invalid number of packet iterations: ", argv[i]+3);
          break;
        case 'i':
          s->max_it = atoi(argv[i] + 3);
          if (s->max_it <= 0 || errno)
            mps_error(s, 2, "Invalid number of iterations: ", argv[i]+3);
          break;
        default:
          mps_error(s, 3, "Bad limit switch: ", argv[i]+2, ", use (p|i)<num>");
        }
        break;

      /* randomization */
      case 'H':
        seed = (unsigned int) atoi(argv[i] + 2);
        if (seed == 0 || errno)
          mps_error(s, 2, "Wrong random seed: ", argv[i]+2);
        break;

      /* debug options */
      case 'd':
        s->DOLOG = true;
        if (strlen(argv[i]) == 3 && argv[i][2] == '1')
          s->logstr = s->outstr;
        else if (strlen(argv[i]) != 2)
          mps_error(s, 3, "Bad debug option: ", argv[i], ", use 1");
        break;

      /* warning options */
      case 'w':
	if (strlen(argv[i]) == 3 && argv[i][2] == '+')
          s->DOWARN = true;
        else if (strlen(argv[i]) == 3 && argv[i][2] == '-')
          s->DOWARN = false;
	else
          mps_error(s, 3, "Bad warning option: ", argv[i], ", use +|-");
        break;
	
      /* help and default case */
      case 'h':
        mps_print_help(s);
        exit(EXIT_SUCCESS);
	break;
	
      /* handle illegal options */ 
      default:
        mps_error(s, 3, "Bad option: ", argv[i], ", type 'unisolve -h' for help");
      }

finalcheck:

  /* If the goal is approximate or count then remove the multiplicity option */
  if (s->goal[0] != 'i')
    s->goal[2] = 'n';

  /* check I/O streams */
  if (s->instr == NULL)
    mps_error(s, 1, "Cannot open input file");
  if (s->outstr == NULL)
    mps_error(s, 1, "Cannot open output file");
  if (s->DOLOG && s->logstr == NULL)
    mps_error(s, 1, "Cannot open log file");

  /* randomize */
  randomize(seed);
}

/**
 * @brief Print usage informations
 */
void
mps_print_help(mps_status* s)
{
  fprintf(s->outstr, "USAGE: unisolve {<options>} {input_file}\n");
  fprintf(s->outstr, " OPTIONS (defaults in square brackets):\n");
  fprintf(s->outstr, " -in\tn = input precision, in decimal digits [0]\n");
  fprintf(s->outstr, " -on\tn = output precision, in decimal digits [%ld]\n",
                  s->prec_out);
  fprintf(s->outstr, " -Gc\tc = Goal: (a)pproximate, [i]solate, (c)ount\n");
  fprintf(s->outstr, " -Sc\tc = Search set: [a]ll, (r)ight, (l)eft, (u)p, "
                  "(d)own complex plane\n");
  fprintf(s->outstr, "   \t                (i)n, (o)ut unitary circle\n");
  fprintf(s->outstr, " -Mc\tc = Multiplicity check: (+) on, [-] off\n");
  fprintf(s->outstr, " -Dc\tc = Detect: [n]one, (r)eal/(i)maginary/(b)oth\n");
  fprintf(s->outstr, " -Oc\tc = Output format: (b)are, (g)nuplot, [c]ompact, "
                  "(v)erbose, (f)ull\n");
  fprintf(s->outstr, " -Hn\tn = random seed, taken from /dev/random if exists\n");
  fprintf(s->outstr, " -d\tprint debug information to standard error\n");
  fprintf(s->outstr, " -d1\tprint debug information to standard output\n");
  fprintf(s->outstr, " -wc\tc = print warning messages: [+] on, (-) off\n");
  fprintf(s->outstr, " -Lpn\tn = maximum number of packet iterations [%d]\n",
                  s->max_pack);
  fprintf(s->outstr, " -Lin\tn = maximum number of global iterations [%d]\n",
                  s->max_it);
  fprintf(s->outstr, " -r\tuse random seed for starting points\n");

  /* undocumented switches
     C
     R
  */
}
