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
#include "mps_poly.h"

/***********************************************************
 *                 MAIN                                    *
 **********************************************************/
int
main(int argc, char *argv[])
{
  mpspoly_t p;
  
  /* set default goal */
  strncpy(goal, "iannc", 5);

  /* parse command line options */
  parse_opts(argc, argv);

  /* Read polynomial */
  read_poly(instr, p);
  
  /* Set polynomial */
  set_poly(p);

  /* allocate global variables */
  allocate_data();

  /* call free_data at exit */
  atexit(free_data);

  /* approximate roots */
  mpsolve();
  
  /* copy roots */
  copy_roots();

  /* output roots */
  output();

  /* return */
  return 0;
}
