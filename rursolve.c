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
#include "rursolve.h"

/***********************************************************
 *                 MAIN                                    *
 **********************************************************/
int
main(void)
{
  /* set default values */
  prec_in = -1;			/* default input precision */
  prec_out = 1000;		/* default output precision */
  strncpy(goal, "ianrv", 5);	/* default goal */

  /* set flags */
  DOLOG = false;
  DOWARN = false;
  DOSORT = true;

  /* set default streams */
  instr = stdin;
  outstr = stdout;
  logstr = stderr;

  /* check I/O streams */
  if (instr == NULL)
    mps_error(s, 1, "Cannot open input file");
  if (outstr == NULL)
    mps_error(s, 1, "Cannot open output file");
  if (DOLOG && logstr == NULL)
    mps_error(s, 1, "Cannot open log file");

  /* compute multivariate roots */
  rursolve();

  return 0;
}
