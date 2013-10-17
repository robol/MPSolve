/*
 * This file is part of MPSolve 3.0
 *
 * Copyright (C) 2001-2013, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors: 
 *   Leonardo Robol <robol@mail.dm.unipi.it>
 */

#define _GNU_SOURCE

#include <mps/mps.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <locale.h>

void
mps_skip_comments (FILE * input_stream)
{
  int buf;
  while ((buf = fgetc (input_stream)) == '!' || isspace (buf))
    if (buf == '!')
      /* Skip until newline */
      while (fgetc (input_stream) != '\n');
  ungetc (buf, input_stream);
}

void
mps_raise_parsing_error (mps_context * s, mps_input_buffer * buffer, 
                         const char * token, 
                         const char * message, ...)
{
  if (!token)
    {
      mps_error (s, message);
      return;
    }

  char * output = (char *) mps_malloc (sizeof (char) * (strlen (token) + 256));
  sprintf (output, "Parsing error on line %ld near the token: %s", buffer->line_number, token);

  mps_error (s, output, message);
  free (output);
}

/**
 * @brief Check if the string given is equal to the option
 * identified by the given string
 */
mps_boolean
mps_is_option (mps_context * s, const char *option_string1,
               const char *option_string2)
{
  mps_boolean is_option = false;

  /* Skip initial spaces */
  while (isspace (*option_string1))
    option_string1++;
  while (isspace (*option_string2))
    option_string2++;

  /* Compare char by char */
  while ((tolower (*option_string1) == tolower (*option_string2)) &&
         (*option_string1 != '\0') && (*option_string2 != '\0'))
    {
      option_string1++;
      option_string2++;
    }

  /* We have arrived to the first different char or on the end of
   * one of the two string. If one is at the end, the other must
   * have only space characters left */
  if (*option_string1 == '\0')
    {
      while (isspace (*option_string2))
        option_string2++;

      is_option = (*option_string2 == '\0');
    }
  else if (*option_string2 == '\0')
    {
      while (isspace (*option_string1))
        option_string1++;

      is_option = (*option_string1 == '\0');
    }

  return is_option;
}

/**
 * @brief Parse a line of the input stream that contains the character
 * ';', so should be considered an option line.
 *
 * Valid options, recognized at the moment being are:
 */
mps_input_option
mps_parse_option_line (mps_context * s, char *line, size_t length)
{
  char *first_comment;
  char *option;
  char *c_ptr;
  char *equal_position;
  mps_input_option input_option = { MPS_FLAG_UNDEFINED, NULL };
  size_t real_length;

  if (length > 255)
  {
    mps_error (s,
               "Maximum line length exceeded (length > 255 while parsing)");
    return input_option;
  }

  /* Check if there are comments in this line */
  if ((first_comment = strchr (line, '!')) != NULL)
    real_length = (first_comment - line) / sizeof (char);
  else
    real_length = length;

  c_ptr = line;
  while (isspace (*c_ptr)
         && ((c_ptr < first_comment) || first_comment == NULL))
    {
      c_ptr++;
      real_length--;
    }
  option = c_ptr;
  c_ptr = strchr (option, ';');
  while (isspace (*--c_ptr) && real_length--);

  /* Now we have the option that is pointed by option and is
   * real_length characters long */
  *(c_ptr + 1) = '\0';
  if (s->debug_level & MPS_DEBUG_IO)
    {
      MPS_DEBUG (s, "Parsed option: %s", option);
    }

  input_option.flag = MPS_FLAG_UNDEFINED;
  input_option.value = NULL;

  /* Detect option about density-sparseness */
  if (mps_is_option (s, option, "dense"))
    input_option.flag = MPS_FLAG_DENSE;
  if (mps_is_option (s, option, "sparse"))
    input_option.flag = MPS_FLAG_SPARSE;

  /* Options on types */
  if (mps_is_option (s, option, "integer"))
    input_option.flag = MPS_FLAG_INTEGER;
  if (mps_is_option (s, option, "real"))
    input_option.flag = MPS_FLAG_REAL;
  if (mps_is_option (s, option, "complex"))
    input_option.flag = MPS_FLAG_COMPLEX;
  if (mps_is_option (s, option, "rational"))
    input_option.flag = MPS_FLAG_RATIONAL;
  if (mps_is_option (s, option, "floatingpoint"))
    input_option.flag = MPS_FLAG_FP;

  /* Options on the input type */
  if (mps_is_option (s, option, "secular"))
    input_option.flag = MPS_FLAG_SECULAR;
  if (mps_is_option (s, option, "monomial"))
    input_option.flag = MPS_FLAG_MONOMIAL;
  if (mps_is_option (s, option, "chebyshev"))
    input_option.flag = MPS_FLAG_CHEBYSHEV;

  /* Parsing keys with values. If = is not found in the
   * input string, than an error has occurred so we should
   * return. */
  equal_position = strchr (option, '=');
  if (equal_position == NULL)
    {
      if (input_option.flag == MPS_FLAG_UNDEFINED)
        {
          mps_error (s, "Unrecognized option: %s", option);
        }
      return input_option;
    }
  else
    {
      input_option.value = equal_position + 1;
      /* Make a copy of the option to parse it without
       * equal sign and anything after it */
      c_ptr = option;
      option = (char *) mps_malloc (sizeof (char) * (strlen (option) + 1));
      strcpy (option, c_ptr);
      *strchr (option, '=') = '\0';
    }

  if (mps_is_option (s, option, "degree"))
    input_option.flag = MPS_KEY_DEGREE;
  else if (mps_is_option (s, option, "precision"))
    input_option.flag = MPS_KEY_PRECISION;

  if (input_option.flag == MPS_FLAG_UNDEFINED)
  {
    mps_error (s, "Unrecognized option: %s", option);
  }

  /* Free the copy of the option */
  free (option);
  return input_option;
}


mps_polynomial *
mps_parse_file (mps_context * s, const char * path)
{
  FILE * handle = fopen(path, "r");
  if (!handle) 
    {
      mps_error (s, "Error while opening file: %s", path);
      return NULL;
    }
  else 
    {
      mps_polynomial *poly = mps_parse_stream (s, handle);
      fclose (handle);
      return poly;
    }
}


/**
 * @brief Parse a stream for input data.
 */
mps_polynomial *
mps_parse_stream (mps_context * s, FILE * input_stream)
{
  /* This is needed to avoid strange decimal separators */
  setlocale(LC_NUMERIC, "C");

  mps_boolean parsing_options = true;
  mps_input_buffer *buffer;
  mps_input_option input_option;
  char * line;
  mps_boolean first_pass = true;

  /* Variables used to hold the options given by the user */
  mps_structure structure = MPS_STRUCTURE_COMPLEX_FP;
  mps_density density = MPS_DENSITY_DENSE;
  mps_representation representation = MPS_REPRESENTATION_MONOMIAL;
  long int input_precision = 0;

  mps_polynomial * poly = NULL;

  if (!input_stream)
    input_stream = s->instr;

  /* Create a buffered line reader for the input stream
   * that has been assigned to us */
  buffer = mps_input_buffer_new (input_stream);

  /* Set values for required options so we can identify
   * their omission */
  s->n = -1;

  /* Skip initial comments in the stream */
  mps_skip_comments (input_stream);

  while (parsing_options)
    {
      mps_input_buffer_readline (buffer);
      line = buffer->line;
      if (strchr (line, ';') == NULL || mps_input_buffer_eof (buffer))
        {
          if (s->debug_level & MPS_DEBUG_IO)
            {
              MPS_DEBUG (s, "Finished parsing options");
            }

          if (first_pass)
            {
              /* This may be the case where an old format MPSolve file has been
               * given to MPSolve, since no option has been specified, so trying
               * to parse it that way */
              MPS_DEBUG_WITH_INFO (s, "This is not a MPSolve 3.0 pol file, so trying with 2.x format");
              poly = MPS_POLYNOMIAL (mps_monomial_poly_read_from_stream_v2 (s, buffer));

              if (poly)
		{
		  mps_input_buffer_free (buffer);
		  return poly;
		}
	      else
		return NULL;
            }
          parsing_options = false;
        }
      else
        {
          first_pass = false;
          input_option =
            mps_parse_option_line (s, line, strlen (line));

          if (mps_context_has_errors (s))
            {
              mps_input_buffer_free (buffer);
              return NULL;
            }

          /* Parsing of the degree */
          if (input_option.flag == MPS_KEY_DEGREE)
            {
              s->n = atoi (input_option.value);
              if (s->n <= 0)
              {
                mps_error (s, "Degree must be a positive integer");
                mps_input_buffer_free (buffer);
                return NULL;
              }
            }

          /* Parsing precision of input coefficients */
          if (input_option.flag == MPS_KEY_PRECISION)
            {
              input_precision = atoi (input_option.value) * LOG2_10;
              if (input_precision <= 0)
              {
                mps_error (s, "Precision must be a positive integer");
                mps_input_buffer_free (buffer);
                return NULL;
              }
            }

          /* Parsing of representations */
          else if (input_option.flag == MPS_FLAG_SECULAR)
            representation = MPS_REPRESENTATION_SECULAR;
          else if (input_option.flag == MPS_FLAG_MONOMIAL)
            representation = MPS_REPRESENTATION_MONOMIAL;
          else if (input_option.flag == MPS_FLAG_CHEBYSHEV)
            representation = MPS_REPRESENTATION_CHEBYSHEV;

          /* And of dense and or sparse input */
          else if (input_option.flag == MPS_FLAG_SPARSE)
            density = MPS_DENSITY_SPARSE;
          else if (input_option.flag == MPS_FLAG_DENSE)
            density = MPS_DENSITY_DENSE;              

          /* Parsing of algebraic structure of the input */
          else if (input_option.flag == MPS_FLAG_REAL)
            {
              /* Switch on algebraic structure that is already set */
              if (MPS_STRUCTURE_IS_INTEGER (structure))
                structure = MPS_STRUCTURE_REAL_INTEGER;
              else if (MPS_STRUCTURE_IS_RATIONAL (structure))
                structure = MPS_STRUCTURE_REAL_RATIONAL;
              else if (MPS_STRUCTURE_IS_FP (structure))
                structure = MPS_STRUCTURE_REAL_FP;
            }
          else if (input_option.flag == MPS_FLAG_COMPLEX)
            {
              /* Switch on algebraic structure that is already set */
              if (MPS_STRUCTURE_IS_INTEGER (structure))
                structure = MPS_STRUCTURE_COMPLEX_INTEGER;
              else if (MPS_STRUCTURE_IS_RATIONAL (structure))
                structure = MPS_STRUCTURE_COMPLEX_RATIONAL;
              else if (MPS_STRUCTURE_IS_FP (structure))
                structure = MPS_STRUCTURE_COMPLEX_FP;
            }

          /* Parsing of algebraic structure of the input, i.e.
           * Integer, Rational or floating point */
          else if (input_option.flag == MPS_FLAG_INTEGER)
            {
              if (MPS_STRUCTURE_IS_REAL (structure))
                structure = MPS_STRUCTURE_REAL_INTEGER;
              else if (MPS_STRUCTURE_IS_COMPLEX (structure))
                structure = MPS_STRUCTURE_COMPLEX_INTEGER;
            }
          else if (input_option.flag == MPS_FLAG_RATIONAL)
            {
              if (MPS_STRUCTURE_IS_REAL (structure))
                structure = MPS_STRUCTURE_REAL_RATIONAL;
              else if (MPS_STRUCTURE_IS_COMPLEX (structure))
                structure = MPS_STRUCTURE_COMPLEX_RATIONAL;
            }
          else if (input_option.flag == MPS_FLAG_FP)
            {
              if (MPS_STRUCTURE_IS_REAL (structure))
                structure = MPS_STRUCTURE_REAL_FP;
              else if (MPS_STRUCTURE_IS_COMPLEX (structure))
                structure = MPS_STRUCTURE_COMPLEX_FP;
            }
        }

    }

  /* Since the Degree is a required parameter, we ask that it is provided. */
  if (s->n == -1) 
    {
      mps_error (s,
                 "Degree of the polynomial must be provided via the Degree=%d configuration option.");
      return NULL;
    }
  else if (s->debug_level & MPS_DEBUG_IO)
    {
      MPS_DEBUG (s, "Degree: %d", s->n);
    }

  switch (representation) {
    case MPS_REPRESENTATION_SECULAR:
      if (s->debug_level & MPS_DEBUG_IO)
        {
          MPS_DEBUG (s, "Parsing secular equation from stream");
        }
      poly = MPS_POLYNOMIAL (mps_secular_equation_read_from_stream (s, buffer, structure, density));
      break;

    case MPS_REPRESENTATION_CHEBYSHEV:
      if (s->debug_level & MPS_DEBUG_IO)
        MPS_DEBUG (s, "Parsing mps_chebyshev_poly from stream");
      poly = MPS_POLYNOMIAL (mps_chebyshev_poly_read_from_stream (s, buffer, structure, density));
      break;

    case MPS_REPRESENTATION_MONOMIAL:
    default:
      if (s->debug_level & MPS_DEBUG_IO)
        MPS_DEBUG (s, "Parsing mps_monomial_poly from stream");
      poly = MPS_POLYNOMIAL (mps_monomial_poly_read_from_stream (s, buffer, structure, density));
      break;
  }

  if (poly)
    {
      MPS_POLYNOMIAL (poly)->structure = structure;
      MPS_POLYNOMIAL (poly)->density = density;
      mps_context_set_input_prec (s, input_precision);
    }

  mps_input_buffer_free (buffer);

  return poly;
}

