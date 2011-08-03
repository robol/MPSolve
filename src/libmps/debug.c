/**
 * @file
 * @brief Implementation of the debug functions.
 */

#include <mps/interface.h>
#include <stdarg.h>

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
