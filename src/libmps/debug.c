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

#include <mps/mps.h>
#include <stdarg.h>
#include <time.h>

/*
 * This section contains the routines used to monitor time performance
 * of various sections of the code.
 */

clock_t *
mps_start_timer ()
{
  /* Create a new timer */
  clock_t *my_clock = (clock_t *) mps_malloc (sizeof (clock_t));

  /* Get the time */
  (*my_clock) = clock ();
  return my_clock;
}

unsigned long int
mps_stop_timer (clock_t * my_timer)
{
  /* Get the difference of two clocks and store the 
   * milliseconds. */
  unsigned long int delta = clock () - *my_timer;
  delta /= (CLOCKS_PER_SEC / 1000);

  /* Free the old clock */
  free (my_timer);
  return delta;
}
