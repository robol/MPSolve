/**
 * @file
 * @brief This file contains the implementation of an input buffer that
 * wraps around a FILE* and make parsing it line per line more pleasant
 * that handling the FILE* directly.
 */

#include <stdio.h>
#include <mps/core.h>
#include <string.h>

/**
 * @brief Create a new input buffer associated with the
 * given input stream.
 */
mps_input_buffer*
mps_input_buffer_new(FILE* stream)
{
   mps_input_buffer *buf;
   buf = (mps_input_buffer*) malloc(sizeof(mps_input_buffer));

   /* Set initial values */
   buf->stream = stream;
   buf->line = NULL;

   return buf;
}

void
mps_input_buffer_free(mps_input_buffer *buffer)
{
    if (buffer->line)
        free(buffer->line);

    free(buffer);
}

mps_boolean
mps_input_buffer_eof(mps_input_buffer *buffer)
{
    return feof(buffer->stream);
}

char*
mps_input_buffer_readline(mps_input_buffer *buf)
{
    ssize_t read_chars;
    ssize_t length;

    /* Free the old line if still present */
    if (buf->line != NULL)
        length = strlen(buf->line);

    /* Read a new line. On the first step buf->line is NULL
     * so a new space is allocated in there, that will be
     * reused on the subsequent calls. */
    read_chars = getline(&buf->line,
            &length,
            buf->stream);

    return buf->line;
}
