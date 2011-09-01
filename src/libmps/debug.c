/**
 * @file
 * @brief Implementation of the debug functions.
 */

#include <mps/interface.h>
#include <stdarg.h>
#include <time.h>

#if __STDC_VERSION__ < 199901L
#ifndef DISABLE_DEBUG
void
__c_impl__MPS_DEBUG (mps_status * s, const char *templ, ...)
{
  va_list ap;
  if (!s->DOLOG)
    return;
  va_start (ap, templ);
  gmp_vfprintf (s->logstr, templ, ap);
  fprintf (s->logstr, "\n");
}

void
__c_impl____MPS_DEBUG (mps_status * s, const char *templ, ...)
{
  va_list ap;
  if (!s->DOLOG)
    return;
  va_start (ap, templ);
  gmp_vfprintf (s->logstr, templ, ap);
}
#endif
#endif

/*
 * This section contains the routines used to monitor time performance
 * of various sections of the code.
 */

clock_t *
mps_start_timer ()
{
  /* Create a new timer */
  clock_t * my_clock = (clock_t*) malloc (sizeof (clock_t));

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
