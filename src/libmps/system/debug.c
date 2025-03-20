/*
 * This file is part of MPSolve 3.2.2
 *
 * Copyright (C) 2001-2020, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Dario Andrea Bini <bini@dm.unipi.it>
 *   Giuseppe Fiorentino <fiorent@dm.unipi.it>
 *   Leonardo Robol <leonardo.robol@unipi.it>
 */


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
  clock_t *my_clock = (clock_t*)mps_malloc (sizeof(clock_t));

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

MPS_PRIVATE void
mps_debug_cluster_structure (mps_context * s)
{
  mps_cluster_item * cluster_item;
  mps_cluster * cluster;
  mps_root * root;
  mps_boolean isolated_roots = false;

  if (!(s->debug_level & MPS_DEBUG_CLUSTER))
    return;

  for (cluster_item = s->clusterization->first; cluster_item != NULL; cluster_item = cluster_item->next)
    {
      cluster = cluster_item->cluster;

      /* Detect that there are isolated roots, so they will be listed
       * after. */
      if (cluster->n == 1)
        {
          isolated_roots = true;
          continue;
        }

      __MPS_DEBUG (s, "Found cluster of %ld roots: ",
                   cluster->n);

      for (root = cluster->first; root != NULL; root = root->next)
        {
          fprintf (s->logstr, "%ld ", root->k);
        }
      fprintf (s->logstr, "\n");
    }

  if (isolated_roots)
    {
      __MPS_DEBUG (s, "Isolated roots: ");
      for (cluster_item = s->clusterization->first; cluster_item != NULL; cluster_item = cluster_item->next)
        {
          cluster = cluster_item->cluster;
          if (cluster->n == 1)
            fprintf (s->logstr, "%ld ", cluster->first->k);
        }
      fprintf (s->logstr, "\n");
    }
}

MPS_PRIVATE void
__c_impl__MPS_DEBUG (mps_context * ctx, const char * format, ...)
{
  va_list ap;

  va_start (ap, format);

  if (ctx->DOLOG)
    {
      fprintf (ctx->logstr, "DEBUG: ");
      vfprintf (ctx->logstr, format, ap);
      fprintf (ctx->logstr, "\n");
    }

  va_end (ap);
}

MPS_PRIVATE void
__c_impl____MPS_DEBUG (mps_context * ctx, const char * format, ...)
{
  va_list ap;

  va_start (ap, format);

  if (ctx->DOLOG)
    {
      fprintf (ctx->logstr, "DEBUG: ");
      vfprintf (ctx->logstr, format, ap);
    }

  va_end (ap);
}


