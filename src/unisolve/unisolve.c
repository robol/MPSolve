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
#include <stdio.h>
#include <mps/mps.h>

/* void */
/* abortfn (enum mcheck_status status) */
/* { */
/*   fprintf (stderr, "A memory error has occurred in MPSolve; aborting\n"); */
/*   abort (); */
/* } */

/***********************************************************
 *                 MAIN                                    *
 **********************************************************/
int
main (int argc, char *argv[])
{
  mps_status *s = mps_status_new ();

  /* Make stdout synchronous so the debugging is more
   * effective. */
  setvbuf (stdout, NULL, _IONBF, 0);

  /* parse command line options */
  mps_parse_opts (s, argc, argv);

  /* Read polynomial */
  mps_parse_stream (s, NULL);

  /* approximate roots */
  mps_mpsolve (s);

  /* output roots */
  mps_output (s);

  /* Free data */
  mps_status_free (s);

  /* return */
  return 0;
}
