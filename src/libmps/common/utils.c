/*
 * This file is part of MPSolve 3.1.5
 *
 * Copyright (C) 2001-2015, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Leonardo Robol <leonardo.robol@sns.it>
 */

#include <mps/mps.h>
#include <string.h>
#include <ctype.h>

/**
 * @brief Strip a string by removing leading and trailing whitespaces.
 *
 * @param ctx The current mps_context.
 * @param input The string to strip.
 * @return A newly allocated stripped version of input. You are in charge
 * to free this string when not needed anymore.
 */
char *
mps_utils_strip_string (mps_context * ctx,
                        const char * input)
{
  const char * input_ptr = input;
  char * start = NULL;
  char * end = NULL;

  while (isspace (*input_ptr) && *input_ptr != '\0')
    input_ptr++;

  start = strdup (input_ptr);
  end = start + strlen (start) - 1;

  while (isspace (*end) && end > start)
    end--;

  *(end + 1) = '\0';
  start = mps_realloc (start, end - start + 2);

  return start;
}

char *
mps_utils_build_equivalent_rational_string (mps_context * ctx,
                                            const char * input)
{
  int sign = 1;
  long int exponent = 0;
  char * ptr = NULL;

  if (input == NULL)
    return strdup ("0");

  ptr =
    build_equivalent_rational_string (ctx, input,
                                      &exponent,
                                      &sign);

  if (!ptr)
    {
      return NULL;
    }

  char * equivalent_rational_string = mps_utils_strip_string (ctx, ptr);

  free (ptr);

  /* Change sign if needed */
  if (sign == -1)
    {
      if (*equivalent_rational_string == '-')
        {
          *equivalent_rational_string = ' ';
        }
      else
        {
          size_t length = strlen (equivalent_rational_string);

          equivalent_rational_string = mps_realloc
                                         (equivalent_rational_string,
                                         length + 2);
          memmove (equivalent_rational_string + 1,
                   equivalent_rational_string, length + 1);
          *equivalent_rational_string = '-';
        }
    }

  /* Insert the exponent into the string */
  if (exponent > 0)
    {
      size_t length = strlen (equivalent_rational_string);
      equivalent_rational_string = mps_realloc
                                     (equivalent_rational_string, length + exponent);

      ptr = strchr (equivalent_rational_string, '/');
      if (ptr)
        {
          /* Move the denominator se we have space to fill the numerator
           * with zeros. */
          memmove (ptr + exponent, ptr,
                   length - (ptr - equivalent_rational_string) + 1);
        }
      else
        {
          /* Set the pointer to the end of the string, that will be filled
           * with zeros. */
          ptr = equivalent_rational_string + length - 1;
          *(ptr + exponent) = '\0';
        }

      memset (ptr, '0', exponent);
    }
  else if (exponent < 0)
    {
      size_t length = strlen (equivalent_rational_string);
      char * slash = strchr (equivalent_rational_string, '/');

      length = length - exponent + ((slash == NULL) ? 1 : 0);

      equivalent_rational_string = mps_realloc
                                     (equivalent_rational_string, length);

      char * ptr = equivalent_rational_string + strlen (equivalent_rational_string);
      if (slash == NULL)
        *ptr++ = '/';

      memset (ptr, '0', (size_t)(-exponent));
      *(ptr - exponent) = '\0';
    }

  return equivalent_rational_string;
}
