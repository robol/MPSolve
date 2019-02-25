/*
 * This file is part of MPSolve 3.1.7
 *
 * Copyright (C) 2001-2019, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Leonardo Robol <leonardo.robol@sns.it>
 */

#include <mps/mps.h>

mps_list_element *
mps_list_element_new (void * value)
{
  mps_list_element * el = mps_new (mps_list_element);

  el->value = value;
  el->previous = el->next = NULL;

  return el;
}

void
mps_list_element_free (mps_list_element * el)
{
  free (el);
}

mps_list_element *
mps_list_element_previous (mps_list_element * el)
{
  return el->previous;
}

mps_list_element * 
mps_list_element_next (mps_list_element * el)
{
  return el->next;
}

void *
mps_list_element_value (mps_list_element * el)
{
  return (el != NULL) ? el->value : NULL;
}
