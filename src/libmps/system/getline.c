/*
 * This file is part of MPSolve 3.1.6
 *
 * Copyright (C) 2001-2015, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 */

/**
 * @file
 * @brief Porting of getline() freely taken from the Android sources. 
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#ifdef getline
#define HAVE_GETLINE
#endif
#endif

#include <mps/mps.h>

/* Workaround some getline detection failures on Windows */
#ifdef getline
  #define HAVE_GETLINE 1
#endif

#ifndef HAVE_GETLINE

#include <stdlib.h>
#include <errno.h>

#ifndef SIZE_MAX
# define SIZE_MAX ((size_t)-1)
#endif

#ifndef SSIZE_MAX
# define SSIZE_MAX ((ssize_t)(SIZE_MAX / 2))
#endif

#if !HAVE_FLOCKFILE
# undef flockfile
# define flockfile(x) ((void)0)
#endif

#if !HAVE_FUNLOCKFILE
# undef funlockfile
# define funlockfile(x) ((void)0)
#endif

/* Read up to (and including) a DELIMITER from FP into *LINEPTR (and
   NUL-terminate it).  *LINEPTR is a pointer returned from malloc (or
   NULL), pointing to *N characters of space.  It is realloc'ed as
   necessary.  Returns the number of characters read (not including
   the null terminator), or -1 on error or EOF.  */

static ssize_t
getdelim (char **lineptr, size_t *n, int delimiter, FILE *fp)
{
  ssize_t result = 0;
  size_t cur_len = 0;

  if (lineptr == NULL || n == NULL || fp == NULL)
    {
      errno = EINVAL;
      return -1;
    }

  flockfile (fp);

  if (*lineptr == NULL || *n == 0)
    {
      *n = 120;
      *lineptr = (char*)malloc (*n);
      if (*lineptr == NULL)
        {
          result = -1;
          goto unlock_return;
        }
    }

  for (;; )
    {
      int i;

      i = getc (fp);
      if (i == EOF)
        {
          result = -1;
          break;
        }

      /* Make enough space for len+1 (for final NUL) bytes.  */
      if (cur_len + 1 >= *n)
        {
          size_t needed_max =
            SSIZE_MAX < SIZE_MAX ? (size_t)SSIZE_MAX + 1 : SIZE_MAX;
          size_t needed = 2 * *n + 1;   /* Be generous. */
          char *new_lineptr;

          if (needed_max < needed)
            needed = needed_max;
          if (cur_len + 1 >= needed)
            {
              result = -1;
              goto unlock_return;
            }

          new_lineptr = (char*)realloc (*lineptr, needed);
          if (new_lineptr == NULL)
            {
              result = -1;
              goto unlock_return;
            }

          *lineptr = new_lineptr;
          *n = needed;
        }

      (*lineptr)[cur_len] = i;
      cur_len++;

      if (i == delimiter)
        break;
    }
  (*lineptr)[cur_len] = '\0';
  result = cur_len ? cur_len : result;

unlock_return:
  funlockfile (fp);
  return result;
}

ssize_t
getline (char **lineptr, size_t *n, FILE *stream)
{
  return getdelim (lineptr, n, '\n', stream);
}

#endif
