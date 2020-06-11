/*
 * This file is part of MPSolve 3.2.0
 *
 * Copyright (C) 2001-2020, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Dario Andrea Bini <bini@dm.unipi.it>
 *   Giuseppe Fiorentino <fiorent@dm.unipi.it>
 *   Leonardo Robol <leonardo.robol@unipi.it>
 */

#define _GNU_SOURCE

#include <stdio.h>
#include <mps/mps.h>
#include <string.h>
#include <ctype.h>

/**
 * @brief Create a new input buffer associated with the
 * given input stream.
 */
mps_input_buffer *
mps_input_buffer_new (mps_abstract_input_stream * stream)
{
  mps_input_buffer *buf;
  int i;

  buf = (mps_input_buffer*)mps_malloc (sizeof(mps_input_buffer));

  buf->last_token = NULL;

  /* Set initial values */
  buf->stream = stream;
  buf->line = NULL;
  buf->line_number = 0L;

  /* Set history size */
  buf->history_size = MPS_INPUT_BUFFER_HISTORY_DEFAULT_SIZE;

  /* Allocate space for the lines kept in history */
  buf->history = (char**)mps_malloc (sizeof(char *) * buf->history_size);
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
  return mps_abstract_input_stream_eof (buffer->stream);
}

/**
 * @brief Read a new line in the buffer, replacing the one
 * present now.
 */
int
mps_input_buffer_readline (mps_input_buffer * buf)
{
  int read_chars = 0;
  size_t length = 0;
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
   * reused on the subsequent calls. If trimming the comment
   * results in an empty line, read a new one. */
  do {
    read_chars = mps_abstract_input_stream_readline (buf->stream, &buf->line, &length);

    if (read_chars > 0)
      buf->last_token = buf->line;

    if (buf->line && read_chars > 0)
      {
	buf->line_number++;
	char * comment = strstr (buf->line, "!");
	if (comment)
	  {
	    *comment = '\0';
	    read_chars = (comment - buf->line);
	  }
      }    
  } while (read_chars == 0 && buf->line != NULL); 

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

  /* Find the token */
  while (*buf->last_token != '\0' && isspace (*buf->last_token))
    buf->last_token++;

  if (*buf->last_token != '\0')
    token = buf->last_token;

  while (token && *token != '\0' && !isspace (*token))
    token++;

  /* Check if we have parsed something or if we need to read another line */
  if (token == NULL)
    {
      if (mps_input_buffer_readline (buf) == -1)
        return NULL;
      return mps_input_buffer_next_token (buf);
    }
  else
    token_size = token - buf->last_token;

  /* Allocate the space for the token if we have found it */
  ret = (char*)mps_malloc (sizeof(char) * (token_size + 1));

  /* Copy the token in ret and set the NULL character in the end */
  strncpy (ret, buf->last_token, token_size);
  ret[token_size] = '\0';

  buf->last_token = token + 1;

  /* Check that we haven't reached the end of the string. If that's the
   * case mark it clearly in the last_stoken field. */
  if (*token == '\0')
    *buf->last_token = '\0';

  return ret;
}
