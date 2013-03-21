/*
 * This file is part of MPSolve 3.0
 *
 * Copyright (C) 2001-2013, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors: 
 *   Dario Andrea Bini <bini@dm.unipi.it>
 *   Giuseppe Fiorentino <fiorent@dm.unipi.it>
 *   Leonardo Robol <robol@mail.dm.unipi.it>
 */

#define _GNU_SOURCE

#include <stdio.h>
#include <mps/mps.h>
#include <string.h>
#include <ctype.h>

#ifndef HAVE_GETLINE

#include <stdlib.h>
#include <errno.h>

#ifndef SIZE_MAX
# define SIZE_MAX ((size_t) -1)
#endif

#ifndef SSIZE_MAX
# define SSIZE_MAX ((ssize_t) (SIZE_MAX / 2))
#endif

#if !HAVE_FLOCKFILE
# undef flockfile
# define flockfile(x) ((void) 0)
#endif

#if !HAVE_FUNLOCKFILE
# undef funlockfile
# define funlockfile(x) ((void) 0)
#endif

/* Read up to (and including) a DELIMITER from FP into *LINEPTR (and
   NUL-terminate it).  *LINEPTR is a pointer returned from malloc (or
   NULL), pointing to *N characters of space.  It is realloc'ed as
   necessary.  Returns the number of characters read (not including
   the null terminator), or -1 on error or EOF.  */

ssize_t
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
      *lineptr = (char *) malloc (*n);
      if (*lineptr == NULL)
        {
          result = -1;
          goto unlock_return;
        }
    }

  for (;;)
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
            SSIZE_MAX < SIZE_MAX ? (size_t) SSIZE_MAX + 1 : SIZE_MAX;
          size_t needed = 2 * *n + 1;   /* Be generous. */
          char *new_lineptr;

          if (needed_max < needed)
            needed = needed_max;
          if (cur_len + 1 >= needed)
            {
              result = -1;
              goto unlock_return;
            }

          new_lineptr = (char *) realloc (*lineptr, needed);
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

/**
 * @brief Create a new input buffer associated with the
 * given input stream.
 */
mps_input_buffer *
mps_input_buffer_new (FILE * stream)
{
  mps_input_buffer *buf;
  int i;
  buf = (mps_input_buffer *) mps_malloc (sizeof (mps_input_buffer));

  buf->last_token = NULL;

  /* Set initial values */
  buf->stream = stream;
  buf->line = NULL;
  buf->line_number = 0L;

  /* Set history size */
  buf->history_size = MPS_INPUT_BUFFER_HISTORY_DEFAULT_SIZE;

  /* Allocate space for the lines kept in history */
  buf->history = (char **) mps_malloc (sizeof (char *) * buf->history_size);
  for (i = 0; i < buf->history_size; ++i)
      buf->history[i] = NULL;
  buf->last = 0;

  return buf;
}

/**
 * @brief Free the buffer and internal data contained in it.
 * Unread lines may be lost.
 */
void
mps_input_buffer_free (mps_input_buffer * buffer)
{
  int i;
  if (buffer->line)
    free (buffer->line);

  for (i = 0; i < buffer->history_size; ++i)
    {
      if (buffer->history[i])
        free (buffer->history[i]);
    }

  free (buffer->history);
  free (buffer);
}

void
mps_input_buffer_set_history_size (mps_input_buffer * buffer, size_t size)
{
  /* TODO: Implement this */
}

/**
 * @brief Check if the whole stream has been read. This does
 * not mean that there is nothing more to read, since the line
 * buffer could be non-empty.
 */
mps_boolean
mps_input_buffer_eof (mps_input_buffer * buffer)
{
  return feof (buffer->stream);
}

/**
 * @brief Read a new line in the buffer, replacing the one
 * present now.
 */
int
mps_input_buffer_readline (mps_input_buffer * buf)
{
  int read_chars = 0;
  size_t length;
  int new_pos;

  /* Move the old line in the buffer, if it's not NULL */
  if (buf->line != NULL)
    {
      new_pos = (buf->last - 1 + buf->history_size) % buf->history_size;
      length = strlen (buf->line);
      
      /* Check if the line that is going to be overwritten
       * has something in it, and if that's the case, free it */
      if (buf->history[new_pos] != NULL)
          free (buf->history[new_pos]);
      
      /* Push the old line in history */
      buf->history[new_pos] = buf->line;
      buf->last = new_pos;
      buf->line = NULL;
    }

  /* Read a new line. On the first step buf->line is NULL
   * so a new space is allocated in there, that will be
   * reused on the subsequent calls. */
  read_chars = getline (&buf->line, &length, buf->stream);

  if (read_chars > 0)
    buf->last_token = buf->line;

  if (buf->line && read_chars > 0) 
    {
      buf->line_number++;
      char * comment = strstr (buf->line, "!");
      if (comment)
        {
          *comment = '\0';
          buf->line = realloc (buf->line, comment - buf->line + 1);
        }
      
    }
  
  return read_chars;
}

/**
 * @brief This function returns the next token that is in the buffer
 * but hasn't been read yet. 
 *
 * It will automagically read new lines if the one in the buffer does
 * not contains anything useful, and return NULL if the stream finish.
 *
 * The returned token shall be freed by the caller.
 */
char *
mps_input_buffer_next_token (mps_input_buffer * buf)
{
  char * token = NULL;
  char * ret;
  size_t token_size = 0;

  if (!buf->line)
    {
      if (mps_input_buffer_readline (buf) == -1)
        return NULL;
    }

  if (!buf->last_token)
    return NULL;

  do {
    /* See if we have found the starting of the token, selecting 
    * things that are not spaces nor end NULL characters. */
    if (!(isspace (*buf->last_token) || 
          (*buf->last_token == '\0')) && 
        (token == NULL))
      {
        if (buf->last_token == '\0')
          break;
        token = buf->last_token;
      }

    /* If we have already started parsing, then increase dimension */
    if (token)
      token_size++;
    buf->last_token++;

  } while ( ((token == NULL) || !isspace (*buf->last_token)) && 
            (*buf->last_token != '\0') );

  /* Check if we have parsed something or if we need to read another line */
  if (token == NULL)
    {
      if (mps_input_buffer_readline (buf) == -1)
        return NULL;
      return mps_input_buffer_next_token (buf);
    }

  /* Allocate the space for the token if we have found it */
  ret = (char *) mps_malloc (sizeof (char) * (token_size + 1));
  
  /* Copy the token in ret and set the NULL character in the end */
  strncpy (ret, token, token_size);
  ret[token_size] = '\0';

  return ret;
}
