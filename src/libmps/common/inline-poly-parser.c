/*
 * This file is part of MPSolve 3.1.5
 *
 * Copyright (C) 2001-2015, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Leonardo Robol <leonardo.robol@sns.it>
 */

/* This is required for old glibc */
#define _GNU_SOURCE

#include <mps/mps.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <ctype.h>

/* Internal states of the parser */
#define PARSING_SIGN        1
#define PARSING_COEFFICIENT 2
#define PARSING_EXPONENT    3
#define PARSING_ERROR        4
#define PARSING_RESET        5

static char * parse_sign (mps_context * ctx, char * line, int * sign, mps_boolean *sign_found);

#ifndef HAVE_FMEMOPEN
/* Forward declaration. Implementation of this function, in case is not provided by the
 * system, is located in input-output.c */
FILE * fmemopen (void *buf, size_t size, const char *mode);;
#endif

static char * find_fp_separator (mps_context  * ctx, char * line)
{
  while (!isspace (*line) && *line != '\0')
    {
      if (*line == '.')
        return line;

      line++;
    }

  return NULL;
}

static long int
parse_fp_exponent (mps_context * ctx, char * exponent_start)
{
  char * ptr = exponent_start;
  char * copy = NULL;

  while (*ptr != 'x' && *ptr != '\0')
    ptr++;

  copy = mps_newv (char, ptr - exponent_start + 1);
  strncpy (copy, exponent_start, ptr - exponent_start);
  *(copy + (ptr - exponent_start)) = '\0';

  long int response = strtol (copy, &ptr, 10);

  if (*ptr != '\0')
    mps_error (ctx, "Error parsing exponent of coefficient: %s", copy);

  free (copy);

  return response;
}

MPS_PRIVATE
char * build_equivalent_rational_string (mps_context * ctx, const char * orig_line, long int * exponent, int * sign)
{
  int i;
  char * line = NULL;

  /* This an over estimate of the length of the _real_ token that we
   * have to parse but, given that we use mps_input_buffer to perform
   * the tokenization of the input, it will not hopefully be a "big"
   * over-estimate. */
  char * allocated_line = strdup (orig_line);
  char * copy = mps_newv (char, 2 * strlen (orig_line) + 5);
  char * copy_ptr = copy;
  char * line_ptr = allocated_line;
  char * line_end = allocated_line + strlen (orig_line);
  long int denominator = 0;
  mps_boolean dot_found = false;
  mps_boolean sign_found = false;

  line = allocated_line;

  char * sep = find_fp_separator (ctx, line);

  /* Note that a string could have a prepended sign */
  line = parse_sign (ctx, line, sign, &sign_found);

  /* Scan the string and truncate it if necessary */
  while (line_ptr++ < line_end)
    {
      if ((*line_ptr == '+' || *line_ptr == '-') &&
          (*(line_ptr - 1) != 'e' && *(line_ptr - 1) != 'E'))
        *line_ptr = '\0';
    }

  line_ptr = line;
  line_end = line_ptr + strlen (line);

  /* If we have floating point input and also a rational
   * separator raise an error. */
  if ((sep || strchr (line, 'e') || strchr (line, 'E'))
      && (strchr (line, '/') != NULL))
    {
      free (line);
      free (copy);
      return NULL;
    }

  *exponent = 0;

  {
    while (line_ptr < line_end)
      {
        /* Check that we don't meet 'e' or 'E'. Otherwise we should exit
         * and parse the exponent instead. */
        if (*line_ptr == 'e' || *line_ptr == 'E')
          {
            *exponent = parse_fp_exponent (ctx, line_ptr + 1);
            break;
          }

        if (*line_ptr == 'x' || *line_ptr == '+' || *line_ptr == '-')
          break;

        if (*line_ptr == '.')
          {
            dot_found = true;
          }
	else {
          denominator+= dot_found ? 1 : 0;
	  *copy_ptr = *line_ptr;
	  copy_ptr++;
	}

	line_ptr++;

      }

    /* TODO: Add the denominator part */
    if (denominator > 0)
      {
        *copy_ptr++ = '/';
        *copy_ptr++ = '1';

        for (i = 0; i < denominator; i++)
          {
            *copy_ptr++ = '0';
          }
      }

    /* Close the token */
    *copy_ptr = '\0';
  }

  /* Clean the x^D part */
  for (copy_ptr = copy; copy_ptr < copy + strlen (copy); copy_ptr++)
    {
      if (*copy_ptr == 'x')
        {
          *copy_ptr = '\0';
          break;
        }
    }

  /* Remove leading zeros, if any */
  copy_ptr = copy;
  while ((copy_ptr < copy + strlen(copy) - 1) &&
         (*copy_ptr == '0'))
    copy_ptr++;

  /* In this case, shift all the string back so that
   * the additional zeros go away. */
  if (copy_ptr != copy)
    {
      int shift = copy_ptr - copy;
      for (i = 0; i < strlen(copy) - shift + 1; i++)
        copy[i] = copy[i + shift];
    }

  if (allocated_line)
    free (allocated_line);

  return copy;
}

static char *
parse_sign (mps_context * ctx, char * line, int * sign, mps_boolean *sign_found)
{
  while (isspace (*line) || *line == '-' || *line == '+')
    {
      if (*line == '-')
        *sign *= -1;

      if (*line == '-' || *line == '+')
        {
          *sign_found = true;
        }

      line++;
    }

  return line;
}


/**
 * @brief Parse a polynomial described the "usual" way, i.e., written
 * as a_k x^k + a_{k-1} x^{k-1} + ... + a_0.
 *
 * @param ctx The current mps_context.
 * @param stream The input stream that shall be used to read the input
 * polynomial.
 */
MPS_PRIVATE mps_polynomial *
mps_parse_inline_poly_from_stream (mps_context *ctx, mps_abstract_input_stream * stream)
{
  return MPS_POLYNOMIAL (mps_monomial_yacc_parser (ctx, stream));
}
/**
 * @brief Parse a polynomial described the "usual" way, i.e., written
 * as a_k x^k + a_{k-1} x^{k-1} + ... + a_0.
 *
 * @param ctx The current mps_context.
 * @param handle A FILE* handle from which the polynomial should be read. 
 */
MPS_PRIVATE mps_polynomial *
mps_parse_inline_poly (mps_context *ctx, FILE * handle)
{
  mps_file_input_stream * stream = mps_file_input_stream_new (handle);
  mps_polynomial * p = mps_parse_inline_poly_from_stream (ctx, (mps_abstract_input_stream*) stream);
  mps_file_input_stream_free (stream);

  return p;
}

/**
 * @brief Parse a polynomial described the "usual" way, i.e., written
 * as a_k x^k + a_{k-1} x^{k-1} + ... + a_0.
 *
 * @param ctx The current mps_context.
 * @param handle A string from which the polynomial should be read. 
 */
mps_polynomial *
mps_parse_inline_poly_from_string (mps_context * ctx, const char * input)
{
  char * input_copy = strdup (input);

  mps_memory_file_stream * stream = mps_memory_file_stream_new (input_copy);
  mps_polynomial * poly = mps_parse_inline_poly_from_stream (ctx, (mps_abstract_input_stream*) stream);

  mps_memory_file_stream_free (stream);
  free (input_copy);
  
  return poly;
}
