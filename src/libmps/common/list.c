/*
 * This file is part of MPSolve 3.2.1
 *
 * Copyright (C) 2001-2020, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Leonardo Robol <leonardo.robol@unipi.it>
 */

#include <mps/mps.h>

/**
 * @brief Create a new empty list.
 */
mps_list *
mps_list_new (void)
{
  mps_list * l = mps_new (mps_list);

  l->first = l->last = NULL;
  l->size = 0;

  return l;
}

/**
 * @brief Free a list and all the elements inside it. 
 */
void
mps_list_free (mps_list * list)
{
  /* First delete all the elements in the list */
  mps_list_element * el = NULL;
  for (el = mps_list_first (list); el != NULL; el = mps_list_element_next (el))
    {
      mps_list_element_free (el);
    }

  free (list);
}

/**
 * @brief Return the number of elements in a list. 
 */
int 
mps_list_size (mps_list * list)
{
  return list->size;
}

void
mps_list_append (mps_list * l, mps_list_element * el)
{
  if (l->last)
    {
      l->last->next = el;
      el->previous = l->last;
      el->next = NULL;
      l->last = el;      
    }
  else
    {
      l->first = l->last = el;
      el->previous = el->next = NULL;
    }

  l->size++;
}

mps_list_element *
mps_list_first (mps_list * l)
{
  return l->first;
}

mps_list_element *
mps_list_last (mps_list * l)
{
  return l->last;
}

