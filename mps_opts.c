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

#include <string.h>
#include <errno.h>
#include "mps.h"

/* local definition */
void print_help(void);

/***********************************************************
 *                 PARSE_OPTS                              *
 ***********************************************************
 * check command line options and streams                  *
 **********************************************************/
void
parse_opts(int argc, char *argv[])
{
  unsigned int seed = 0;
  int i;

  /* set default flags */
  DOLOG = false;
  DOWARN = true;
  DOSORT = true;

  /* set default streams */
  instr = stdin;
  outstr = stdout;
  logstr = stderr;
  rtstr = stdin;

  /* set default values */
  prec_in = -1;                 /* if != -1 ignore precision from file */
  prec_out = 2 * DBL_DIG;       /* default output precision            */

  /* parse options */
  for (i = 1; i < argc; i++)
    if (argv[i][0] != '-') {
      instr = fopen(argv[i], "r");
      goto finalcheck;  /* no options allowed after filename  */
    } else      /* all options start with - */
      switch (argv[i][1]) {
      case 'i':
        prec_in = atol(argv[i] + 2);
        if (prec_in <= 0 || errno)
          error(2, "Wrong input precision: ", argv[i]+2);
        /* NOTAMIA: Perchè sovrascrivere prec_in che è stata appena acquisita
         * come parametro? */
        prec_in = (long) (prec_out * LOG2_10);
        break;

      case 'o':
        prec_out = atol(argv[i] + 2);
        if (prec_out <= 0 || errno)
          error(2, "Wrong output precision: ", argv[i]+2);
        prec_out = (long) (prec_out * LOG2_10);
        break;

      /* goal */
      case 'G':
        switch (argv[i][2]) {
        case 'a':
          goal[0] = 'a';
          break;
        case 'c':
          goal[0] = 'c';
          break;
        case 'i':
          goal[0] = 'i';
          break;
        default:
          error(3, "Bad goal switch: ", argv[i]+2, ", use a|c|i");
        }
        if (strlen(argv[i]) != 3)
          error(2, "Bad goal: ", argv[i]);
        break;

      /* select search set */
      case 'S':
        switch (argv[i][2]) {
        case 'a':
          goal[1] = 'a';
          break;
        case 'r':
          goal[1] = 'r';
          break;
        case 'l':
          goal[1] = 'l';
          break;
        case 'u':
          goal[1] = 'u';
          break;
        case 'd':
          goal[1] = 'd';
          break;
        case 'i':
          goal[1] = 'i';
          break;
        case 'o':
          goal[1] = 'o';
          break;
        case 'R':
          goal[1] = 'R';
          break;
        case 'I':
          goal[1] = 'I';
          break;
        case 'U':
          goal[1] = 'U';
          break;
        default:
          error(3, "Bad search set switch: ", argv[i]+2, ", use a|r|l|u|d|i|o|R|I|U");
        }
        if (strlen(argv[i]) != 3)
          error(2, "Bad set: ", argv[i]);
        break;

      /* select multiplicity */
      case 'M':
        switch (argv[i][2]) {
        case '+':
          goal[2] = 'm';
          break;
        case '-':
          goal[2] = 'n';
          break;
        default:
          error(3, "Bad multiplicity switch: ", argv[i]+2, ", use +|-");
        }
        if (strlen(argv[i]) != 3)
          error(2, "Bad multiplicity option: ", argv[i]);
        break;

      /* detection */
      case 'D':
        switch (argv[i][2]) {
        case 'n':
          goal[3] = 'n';
          break;
        case 'r':
          goal[3] = 'r';
          break;
        case 'i':
          goal[3] = 'i';
          break;
        case 'b':
          goal[3] = 'b';
          break;
        default:
          error(3, "Bad detection switch: ", argv[i]+2, ", use n|r|i|b");
        }
        if (strlen(argv[i]) != 3)
          error(2, "Bad detection option: ", argv[i]);
        break;

      /* I/O streams */
      case 'R':
        rtstr = fopen(argv[i] + 2, "r");
        if (rtstr == NULL)
          error(2, "Cannot open roots file: ", argv[i]+2);
        resume = true;
        break;

      /* Additional checks */
      case 'C':
        switch (argv[i][2]) {
        case 'R':
          chkrad = true;
          break;
        case 'r':
          chkrad = false;
          break;
        default:
          error(3, "Bad check switch: ", argv[i]+2, ", use R|r");
        }
        if (strlen(argv[i]) != 3)
          error(2, "Bad check option: ", argv[i]);
        break;

      /* output format */
      case 'O':
        switch (argv[i][2]) {
        case 'b':
          goal[4] = 'b';
          break;
        case 'g':
          goal[4] = 'g';
          break;
        case 'c':
          goal[4] = 'c';
          break;
        case 'v':
          goal[4] = 'v';
          break;
        case 'f':
          goal[4] = 'f';
          break;
        default:
          error(3, "Bad output format switch: ", argv[i]+2, ", use b|c|f|g|v");
        }
        if (strlen(argv[i]) != 3)
          error(2, "Bad output option: ", argv[i]);
        break;

      /* iteration limits */ 
      case 'L':
        switch (argv[i][2]) {
        case 'p':
          max_pack = atoi(argv[i] + 3);
          if (max_pack <= 0 || errno)
            error(2, "Invalid number of packet iterations: ", argv[i]+3);
          break;
        case 'i':
          max_it = atoi(argv[i] + 3);
          if (max_it <= 0 || errno)
            error(2, "Invalid number of iterations: ", argv[i]+3);
          break;
        default:
          error(3, "Bad limit switch: ", argv[i]+2, ", use (p|i)<num>");
        }
        break;

      /* randomization */
      case 'H':
        seed = (unsigned int) atoi(argv[i] + 2);
        if (seed == 0 || errno)
          error(2, "Wrong random seed: ", argv[i]+2);
        break;

      /* debug options */
      case 'd':
        DOLOG = true;
        if (strlen(argv[i]) == 3 && argv[i][2] == '1')
          logstr = outstr;
        else if (strlen(argv[i]) != 2)
          error(3, "Bad debug option: ", argv[i], ", use 1");
        break;

      /* warning options */
      case 'w':
	if (strlen(argv[i]) == 3 && argv[i][2] == '+')
          DOWARN = true;
        else if (strlen(argv[i]) == 3 && argv[i][2] == '-')
          DOWARN = false;
	else
          error(3, "Bad warning option: ", argv[i], ", use +|-");
        break;
	
      /* help and default case */
      case 'h':
        print_help();
        exit(EXIT_SUCCESS);
	break;
	
      /* handle illegal options */ 
      default:
        error(3, "Bad option: ", argv[i], ", type 'unisolve -h' for help");
      }

finalcheck:

  /* If the goal is approximate or count then remove the multiplicity option */
  if (goal[0] != 'i')
    goal[2] = 'n';

  /* check I/O streams */
  if (instr == NULL)
    error(1, "Cannot open input file");
  if (outstr == NULL)
    error(1, "Cannot open output file");
  if (DOLOG && logstr == NULL)
    error(1, "Cannot open log file");

  /* randomize */
  randomize(seed);
}

/***********************************************************
 *                 PRINT_HELP                              *
 ***********************************************************
 * printout a short command line help file                 *
 **********************************************************/
void
print_help(void)
{
  fprintf(outstr, "USAGE: unisolve {<options>} {input_file}\n");
  fprintf(outstr, " OPTIONS (defaults in square brackets):\n");
  fprintf(outstr, " -in\tn = input precision, in decimal digits [0]\n");
  fprintf(outstr, " -on\tn = output precision, in decimal digits [%ld]\n",
                  prec_out);
  fprintf(outstr, " -Gc\tc = Goal: (a)pproximate, [i]solate, (c)ount\n");
  fprintf(outstr, " -Sc\tc = Search set: [a]ll, (r)ight, (l)eft, (u)p, "
                  "(d)own complex plane\n");
  fprintf(outstr, "   \t                (i)n, (o)ut unitary circle\n");
  fprintf(outstr, " -Mc\tc = Multiplicity check: (+) on, [-] off\n");
  fprintf(outstr, " -Dc\tc = Detect: [n]one, (r)eal/(i)maginary/(b)oth\n");
  fprintf(outstr, " -Oc\tc = Output format: (b)are, (g)nuplot, [c]ompact, "
                  "(v)erbose, (f)ull\n");
  fprintf(outstr, " -Hn\tn = random seed, taken from /dev/random if exists\n");
  fprintf(outstr, " -d\tprint debug information to standard error\n");
  fprintf(outstr, " -d1\tprint debug information to standard output\n");
  fprintf(outstr, " -wc\tc = print warning messages: [+] on, (-) off\n");
  fprintf(outstr, " -Lpn\tn = maximum number of packet iterations [%d]\n",
                  max_pack);
  fprintf(outstr, " -Lin\tn = maximum number of global iterations [%d]\n",
                  max_it);

  /* undocumented switches
     C
     R
  */
}
