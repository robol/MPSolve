/*
 * This file is part of MPSolve 3.2.2
 *
 * Copyright (C) 2001-2020, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Leonardo Robol <leonardo.robol@unipi.it>
 */

#ifndef MPS_INPUT_BUFFER_H_
#define MPS_INPUT_BUFFER_H_

#include <mps/private/system/abstract-input-stream.h>

/**
 * @file
 * @brief Implementation of a buffer for parsing input file for MPSolve.
 */

MPS_BEGIN_DECLS

/**
 * @brief Default size of the buffer of old lines in the input
 * buffer instances.
 */
#define MPS_INPUT_BUFFER_HISTORY_DEFAULT_SIZE 2

/**
 * @brief Buffer used to parse input files in MPSolve. It can
 * read a stream line by line.
 */
struct mps_input_buffer {
  /**
   * @brief Stream associated with the
   * mps_input_buffer
   */
  mps_abstract_input_stream *stream;

  /**
   * @brief Pointer the last line read in the
   * buffer. Another line can be read with
   * <code>mps_input_buffer_readline()</code>
   */
  char *line;

  /**
   * @brief Number of the last read line, the one that
   * is stored in line field.
   */
  long int line_number;

  /**
   * @brief Lines that have been read before this.
   *
   * The number of lines that are remembered is set
   * by the variable MPS_INPUT_BUFFER_HISTORY_DEFAULT_SIZE,
   * but can be overriden by calling
   * <code>mps_input_buffer_set_history_size()</code>.
   */
  char **history;

  /**
   * @brief Size of the history that is been kept in memory.
   *
   * This value should never be modified directly, but can
   * be tweaked using <code>mps_input_buffer_set_history_size()</code>.
   */
  size_t history_size;

  /**
   * @brief Index of the last line inserted in history.
   *
   * This is used internally by the mps_input_buffer to implement
   * a circular buffer.
   */
  int last;

  /**
   * @brief This is a pointer to the last parsed char in the buffer->line
   * string.
   *
   * It is used by mps_input_buffer_next_token() to determine the last
   * thing read and if there is the need to read another line.
   *
   * As for <code>last</code>, this pointer should never be manually
   * modified, even if you think that you know what you're doing.
   */
  char * last_token;
};

/* Function prototypes */
mps_input_buffer *mps_input_buffer_new (mps_abstract_input_stream * stream);
int mps_input_buffer_readline (mps_input_buffer * buf);
void mps_input_buffer_free (mps_input_buffer * buf);
void mps_input_buffer_set_history_size (mps_input_buffer * buf, size_t size);
mps_boolean mps_input_buffer_eof (mps_input_buffer * buf);
char * mps_input_buffer_next_token (mps_input_buffer * buf);

MPS_END_DECLS

#endif /* ifndef MPS_INPUT_BUFFER_H */
