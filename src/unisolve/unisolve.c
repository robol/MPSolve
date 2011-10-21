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
#include <mps/poly.h>
#include <mps/core.h>
#include <mps/interface.h>

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

  /* Set default values in s */
  mps_set_default_values (s);

  /* set default goal */
  strncpy (s->goal, "iannc", 5);

  /* parse command line options */
  mps_parse_opts (s, argc, argv);

  /* Read polynomial */
  mps_parse_stream (s, s->instr);
  // mps_read_poly (s, s->instr, p);

  /* Set polynomial */
  // mps_set_poly (s, p);

  /* allocate global variables */
  mps_allocate_data (s);

  /* approximate roots */
  mps_mpsolve (s);

  /* copy roots */
  mps_copy_roots (s);

  /* output roots */
  mps_output (s);

  /* Free data */
  mps_status_free (s);

  /* return */
  return 0;
}
