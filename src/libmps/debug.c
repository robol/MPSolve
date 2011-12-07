#include <mps/core.h>
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
