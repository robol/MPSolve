/*
 * This file is part of MPSolve 3.1.8
 *
 * Copyright (C) 2001-2019, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Leonardo Robol <leonardo.robol@sns.it>
 */

/**
 * @file
 * @brief Custom implementation of list inside MPSolve.
 *
 * This implementation may be a custom one or simply a wrapper around something more
 * tested and proved to be working. Its main role is to abstract this choice to the internals
 * of MPSolve.
 */

#ifndef MPS_LIST_H_
#define MPS_LIST_H_

MPS_BEGIN_DECLS

struct mps_list_element {
  /**
   * @brief The value pointed by this list element. Note that holding this pointer
   * in the first field of the struct allows to painlessly cast a mps_list_element
   * to the pointer inside it.
   */
  void * value;

  /**
   * @brief The next element in the list.
   *
   * This value should be NULL if this is the last element of the list
   * or the element is detached.
   */
  struct mps_list_element * next;

  /**
   * @brief The previous element in the list.
   *
   * This value should be NULL if this is the first element of the list
   * or the element is detached.
   */
  struct mps_list_element * previous;
};

struct mps_list {
  /**
   * @brief A pointer to the first element of the list.
   */
  mps_list_element * first;

  /**
   * @brief A pointer to the last element of the list.
   */
  mps_list_element * last;

  /**
   * @brief The number of element of the list.
   */
  int size;
};

mps_list_element * mps_list_element_new (void * value);
void mps_list_element_free (mps_list_element * el);
mps_list_element * mps_list_element_next (mps_list_element * el);
mps_list_element * mps_list_element_previous (mps_list_element * el);
void * mps_list_element_value (mps_list_element * el);

mps_list * mps_list_new (void);
void mps_list_free (mps_list * list);
int mps_list_size (mps_list * list);
void mps_list_append (mps_list * list, mps_list_element * el);
mps_list_element * mps_list_first (mps_list * list);
mps_list_element * mps_list_last (mps_list * list);

/**
 * @brief Shortcut for iterating over lists. Note that this structure is
 * C99 only, since it would not compile under a strict C89 compiler.
 */
#define MPS_LIST_FOREACH(type, local_var, list) \
  for (type *__mps_local_iterator = (type*)mps_list_first (list),              \
       *local_var = (type*)mps_list_element_value ((mps_list_element*)__mps_local_iterator); \
       __mps_local_iterator != NULL;                                    \
       __mps_local_iterator = (type*)mps_list_element_next ((mps_list_element*)__mps_local_iterator), \
       local_var = (type*)mps_list_element_value ((mps_list_element*)__mps_local_iterator))

MPS_END_DECLS

#endif /* MPS_LIST_H_ */

