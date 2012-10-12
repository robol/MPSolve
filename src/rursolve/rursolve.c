/***********************************************************
**       Multiprecision Polynomial Solver (MPSolve)       **
**              Version 2.1, september 1999               **
**                                                        **
**                      Written by                        **
**       Dario Andrea Bini and Giuseppe Fiorentino        **
**       (bini@dm.unipi.it)  (fiorent@dm.unipi.it)        **
**                                                        **
** (C) 1999, Dipartimento di Matematica, FRISCO LTR 21024 **
***********************************************************/

#include <string.h>
#include <mps/rursolve.h>

/***********************************************************
 *                 MAIN                                    *
 **********************************************************/
int
main (int argc, char **argv)
{

  /* Allocate mps_context */
  mps_context *s = malloc (sizeof (mps_context));

  /* set default values */
  s->prec_in = -1;              /* default input precision */
  s->prec_out = 1000;           /* default output precision */
  strncpy (s->goal, "ianrv", 5);        /* default goal */

  /* set flags */
  s->DOLOG = false;
  s->DOWARN = false;
  s->DOSORT = true;

  /* set default streams */
  s->instr = stdin;
  s->outstr = stdout;
  s->logstr = stderr;

  /* check I/O streams */
  if (s->instr == NULL)
    mps_error (s, 1, "Cannot open input file");
  if (s->outstr == NULL)
    mps_error (s, 1, "Cannot open output file");
  if (s->DOLOG && s->logstr == NULL)
    mps_error (s, 1, "Cannot open log file");

  /* compute multivariate roots */
  mps_rursolve (s);

  return 0;
}
