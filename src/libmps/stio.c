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

#include <stdarg.h>
#include <string.h>
#include <mps/mps.h>
#include <ctype.h>

#ifndef __WINDOWS
#include <unistd.h>
#else
#include <io.h>
#endif

#define ISZERO -1

void
mps_skip_comments (FILE * input_stream)
{
  char buf;
  while ((buf = fgetc (input_stream)) == '!' || isspace (buf))
    if (buf == '!')
      /* Skip until newline */
      while (fgetc (input_stream) != '\n');
  ungetc (buf, input_stream);
}

void
mps_raise_parsing_error (mps_context * s, mps_input_buffer * buffer, 
                         const char * token, 
                         const char * message)
{
  if (!token)
    {
      mps_error (s, 1, message);
      return;
    }

  char * output = (char *) mps_malloc (sizeof (char) * (strlen (token) + 256));
  sprintf (output, "Parsing error on line %ld near the token: %s", buffer->line_number, token);

  mps_error (s, 2, output, message);
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
  mps_input_option input_option;
  size_t real_length;

  if (length > 255)
    mps_error (s, 1,
               "Maximum line length exceeded (length > 255 while parsing)");

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

  /* Free the copy of the option */
  free (option);
  return input_option;
}

void
mps_monomial_poly_read_from_stream (mps_context * s,
                                    mps_input_buffer * buffer, 
                                    mps_structure structure,
                                    mps_density density)
{
  mps_monomial_poly * poly;
  int i;
  mpf_t ftmp;
  char * token;

  mpf_init (ftmp);

  /* Allocate space for the polynomial, since we need this even
   * if we are trying to solve the associated secular_equation */
  poly = mps_monomial_poly_new (s, s->n);

  MPS_POLYNOMIAL (poly)->structure = structure;
  MPS_POLYNOMIAL (poly)->density = density;
  MPS_POLYNOMIAL (poly)->prec = 0;

  /* We still do not support sparse input */
  for (i = 0; i <= s->n; ++i)
    poly->spar[i] = true;

  /* Dense parsing */
  if (MPS_DENSITY_IS_DENSE (density))
    {
      if (MPS_STRUCTURE_IS_FP (structure))
        {
          for (i = 0; i < s->n + 1; ++i)
            {
              token = mps_input_buffer_next_token (buffer);
              if (!token || (mpf_set_str (mpc_Re (poly->mfpc[i]), token, 10) != 0))
                mps_raise_parsing_error (s, buffer, token, "Error parsing coefficients of the polynomial");
              free (token);

              if (MPS_STRUCTURE_IS_COMPLEX (structure))
                {
                  token = mps_input_buffer_next_token (buffer);
                  if (!token || (mpf_set_str (mpc_Im (poly->mfpc[i]), token, 10) != 0))
                    mps_raise_parsing_error (s, buffer, token, "Error parsing coefficients of the polynomial");
                  free (token);
                }
              else
                mpf_set_ui (mpc_Im (poly->mfpc[i]), 0U);
            }
        }
      else if (MPS_STRUCTURE_IS_RATIONAL (structure) ||
               MPS_STRUCTURE_IS_INTEGER (structure))
        {
          for (i = 0; i < s->n + 1; ++i)
            {
              token = mps_input_buffer_next_token (buffer);
              if (!token || (mpq_set_str (poly->initial_mqp_r[i], token, 10) != 0))
                mps_raise_parsing_error (s, buffer, token, "Error parsing coefficients of the polynomial");
              mpq_canonicalize (poly->initial_mqp_r[i]);
              free (token);
      
              if (MPS_STRUCTURE_IS_COMPLEX (structure))
                {
                  token = mps_input_buffer_next_token (buffer);
                  if (!token || (mpq_set_str (poly->initial_mqp_i[i], token, 10) != 0))
                    mps_raise_parsing_error (s, buffer, token, "Error parsing coefficients of the polynomial");
                  mpq_canonicalize (poly->initial_mqp_i[i]);
                  free (token);
                }
              else
                mpq_set_ui (poly->initial_mqp_i[i], 0U, 0U);

              /* Copy coefficients in the floating point ones */
              mpf_set_q (mpc_Re (poly->mfpc[i]), poly->initial_mqp_r[i]);
              mpf_set_q (mpc_Im (poly->mfpc[i]), poly->initial_mqp_i[i]);
            }
        }
    } /* closes if (MPS_INPUT_CONFIG_IS_DENSE (s->input_config)) */
  else if (MPS_DENSITY_IS_SPARSE (density))
    {
      /* Set all the spar to false, since we have still not read
       * any coefficient */
      for (i = 0; i <= s->n; ++i)
        poly->spar[i] = false;
      
      while ((token = mps_input_buffer_next_token (buffer)))
        {
          /* Read the index from the buffer */
          if (!sscanf (token, "%d", &i))
            mps_raise_parsing_error (s, buffer, token, "Error while parsing the degree of a monomial");

          if (i < 0 || i > s->n) 
            {
              mps_raise_parsing_error (s, buffer, token, "Degree of coefficient out of bounds");
              free (token);
              return;
            }

          if (poly->spar[i]) 
            mps_raise_parsing_error (s, buffer, token, "A monomial of the same degree has been inserted twice"); 
          else 
            poly->spar[i] = true;
          free (token);

          if (MPS_STRUCTURE_IS_FP (structure))
            {
              token = mps_input_buffer_next_token (buffer);
              if (!token || (mpf_set_str (mpc_Re (poly->mfpc[i]), token, 10) != 0))
                mps_raise_parsing_error (s, buffer, token, "Error parsing coefficients of the polynomial");
              free (token);
          
              if (MPS_STRUCTURE_IS_COMPLEX (structure))
                {
                  token = mps_input_buffer_next_token (buffer);
                  if (!token || (mpf_set_str (mpc_Im (poly->mfpc[i]), token, 10) != 0))
                    mps_raise_parsing_error (s, buffer, token, "Error parsing coefficients of the polynomial");
                  free (token);
                }
              else
                mpf_set_ui (mpc_Im (poly->mfpc[i]), 0U);
            }
          else if (MPS_STRUCTURE_IS_RATIONAL (structure) ||
                   MPS_STRUCTURE_IS_INTEGER (structure))
            {

              token = mps_input_buffer_next_token (buffer);
              if (!token || (mpq_set_str (poly->initial_mqp_r[i], token, 10) != 0))
                mps_raise_parsing_error (s, buffer, token, "Error parsing coefficients of the polynomial");
              mpq_canonicalize (poly->initial_mqp_r[i]);
              free (token);
      
              if (MPS_STRUCTURE_IS_COMPLEX (structure))
                {
                  token = mps_input_buffer_next_token (buffer);
                  if (!token || (mpq_set_str (poly->initial_mqp_i[i], token, 10) != 0))
                    mps_raise_parsing_error (s, buffer, token, "Error parsing coefficients of the polynomial");
                  mpq_canonicalize (poly->initial_mqp_i[i]);
                  free (token);
                }
              else
                mpq_set_ui (poly->initial_mqp_i[i], 0U, 0U);

              /* Copy coefficients in the floating point ones */
              mpf_set_q (mpc_Re (poly->mfpc[i]), poly->initial_mqp_r[i]);
              mpf_set_q (mpc_Im (poly->mfpc[i]), poly->initial_mqp_i[i]);
            }
        }
    } /* closes if (MPS_INPUT_CONFIG_IS_SPARSE (s->input_config)) */

  /* Copy coefficients back in other places */
  for (i = 0; i < s->n + 1; ++i)
    {
      if (poly->spar[i])
        {
          mpc_get_cplx (poly->fpc[i], poly->mfpc[i]);
          mpc_get_cdpe (poly->dpc[i], poly->mfpc[i]);

          /* Compute modules of coefficients */
          cdpe_mod (poly->dap[i], poly->dpc[i]);
          poly->fap[i] = rdpe_get_d (poly->dap[i]);

          if (i > 0)
            mpc_mul_ui (poly->mfppc[i-1], poly->mfppc[i], i);

          MPS_DEBUG_RDPE(s, poly->dap[i], "poly->dap[%d]", i);
        }
      else
        {
          cplx_set (poly->fpc[i], cplx_zero);
          cdpe_set (poly->dpc[i], cdpe_zero);
          
          rdpe_set (poly->dap[i], rdpe_zero);
          poly->fap[i] = 0.0f;
        }
    }

  mps_context_set_input_poly (s, MPS_POLYNOMIAL (poly));
  mpf_clear (ftmp);
}

void
mps_secular_equation_read_from_stream (mps_context * s,
                                       mps_input_buffer * buffer,
                                       mps_structure structure,
                                       mps_density density)
{
  mps_secular_equation *sec;
  int i;
  mpf_t ftmp;
  char * token;

  mpf_init (ftmp);

  /* Read directly the secular equation in DPE, so we don't need
   * to have a fallback case if the coefficients are bigger than
   * what is supported by the standard floating point arithmetic */
  sec = mps_secular_equation_new_raw (s, s->n);
  MPS_POLYNOMIAL (sec)->degree = s->n;
  MPS_POLYNOMIAL (sec)->structure = structure;
  MPS_POLYNOMIAL (sec)->density = density;
  MPS_POLYNOMIAL (sec)->prec = DBL_DIG;

  /* Parsing of integers and floating point is done with Multiprecision */
  if (MPS_STRUCTURE_IS_FP (structure))
    {
      for (i = 0; i < s->n; i++)
        {
          token = mps_input_buffer_next_token (buffer);
          if (!token || (mpf_set_str (mpc_Re (sec->initial_ampc[i]), token, 10) != 0))
            {
              MPS_DEBUG (s,
                         "Error reading coefficient a[%d] of the secular equation (real part)",
                         i);
              mps_raise_parsing_error (s, buffer, token, 
                         "Error reading some coefficients of the secular equation.\n"
                         "Please check your input file.");
            }
          free (token);

          /* Imaginary part, read only if the input is complex */
          if (MPS_STRUCTURE_IS_COMPLEX (structure))
            {
              token = mps_input_buffer_next_token (buffer);
              if (!token || (mpf_set_str (mpc_Im (sec->initial_ampc[i]), token, 10) != 0))
                {
                  MPS_DEBUG (s,
                             "Error reading coefficient a[%d] of the secular equation (imaginary part)",
                             i);
                  mps_raise_parsing_error (s, buffer, token,
                             "Error reading some coefficients of the secular equation.\n"
                             "Please check your input file.");
                }
              free (token);
            }
          else
            {
              mpf_set_ui (mpc_Im (sec->initial_ampc[i]), 0U);
            }

          token = mps_input_buffer_next_token (buffer);
          if (!token || (mpf_set_str (mpc_Re (sec->initial_bmpc[i]), token, 10) != 0))
            {
              MPS_DEBUG (s,
                         "Error reading coefficient b[%d] of the secular equation (real part)",
                         i);
              mps_raise_parsing_error (s, buffer, token,
                         "Error reading some coefficients of the secular equation.\n"
                         "Please check your input file.");
            }
          free (token);

          /* Again, read the imaginary part only if the input is complex */
          if (MPS_STRUCTURE_IS_COMPLEX (structure))
            {
              token = mps_input_buffer_next_token (buffer);
              if (!token || (mpf_set_str (mpc_Im (sec->initial_bmpc[i]), token, 10) != 0))
                {
                  MPS_DEBUG (s,
                             "Error reading coefficient b[%d] of the secular equation (imaginary part)",
                             i);
                  mps_raise_parsing_error (s, buffer, token,
                             "Error reading some coefficients of the secular equation.\n"
                             "Please check your input file.");
                }
              free (token);
            }
          else
            {
              mpf_set_ui (mpc_Im (sec->initial_bmpc[i]), 0U);
            }
        }
    }
  /*
   * Parsing of rational and integer input.
   * Parsing of the integer input is done assuming the coefficients
   * as a special case of the rational ones.
   */
  else if (MPS_STRUCTURE_IS_RATIONAL (structure) ||
           MPS_STRUCTURE_IS_INTEGER  (structure))
    {
      for (i = 0; i < s->n; i++)
        {
          /* Read real part of the a_i */
          token = mps_input_buffer_next_token (buffer);
          if (!token || (mpq_set_str (sec->initial_ampqrc[i], token, 10) != 0))
            {
              MPS_DEBUG (s, "Error reading the coefficients a[%d] of the secular equation (real part)", i);
              mps_raise_parsing_error (s, buffer, token, 
                                       "Error reading some coefficients of the secular equation.\nPlease check your input file");
            }
          mpq_canonicalize (sec->initial_ampqrc[i]);
          free (token);

          /* Read imaginary part of the a_i */
          if (MPS_STRUCTURE_IS_COMPLEX (structure))
            {
              token = mps_input_buffer_next_token (buffer);
              if (!token || (mpq_set_str (sec->initial_ampqic[i], token, 10) != 0))
                {             
                  MPS_DEBUG (s, "Error reading the coefficients a[%d] of the secular equation (imaginary part)", i);
                  mps_raise_parsing_error (s, buffer, token, 
                                           "Error reading some coefficients of the secular equation.\nPlease check your input file");
                }
              mpq_canonicalize (sec->initial_ampqic[i]);
              free (token);
            }
          else
            mpq_set_ui (sec->initial_ampqic[i], 0, 0);

          /* Read real part of the b_i */
          token = mps_input_buffer_next_token (buffer);
          if (!token || (mpq_set_str (sec->initial_bmpqrc[i], token, 10) != 0))
            {         
              MPS_DEBUG (s, "Error reading the coefficients b[%d] of the secular equation (real part)", i);
              mps_raise_parsing_error (s, buffer, token, 
                                       "Error reading some coefficients of the secular equation.\nPlease check your input file");
            }       
          mpq_canonicalize (sec->initial_bmpqrc[i]);
          free (token);

          /* Read imaginary part of the b_i */
          if (MPS_STRUCTURE_IS_COMPLEX (structure))
            {
              token = mps_input_buffer_next_token (buffer);
              if (!token || (mpq_set_str (sec->initial_bmpqic[i], token, 10) != 0))
                {             
                  MPS_DEBUG (s, "Error reading the coefficients b[%d] of the secular equation (imaginary part)", i);
                  mps_raise_parsing_error (s, buffer, token, 
                                           "Error reading some coefficients of the secular equation.\nPlease check your input file");
                }           
              mpq_canonicalize (sec->initial_bmpqic[i]);
              free (token);
            }
          else
            mpq_set_ui (sec->initial_bmpqic[i], 0, 0);
        }

      /* Set DPE coefficients */
      for (i = 0; i < s->n; i++)
        {
          mpf_set_q (ftmp, sec->initial_ampqrc[i]);
          mpf_set (mpc_Re (sec->initial_ampc[i]), ftmp);
          mpf_get_rdpe (cdpe_Re (sec->adpc[i]), ftmp);

          mpf_set_q (ftmp, sec->initial_ampqic[i]);
          mpf_set (mpc_Im (sec->initial_ampc[i]), ftmp);
          mpf_get_rdpe (cdpe_Im (sec->adpc[i]), ftmp);

          mpf_set_q (ftmp, sec->initial_bmpqrc[i]);
          mpf_set (mpc_Re (sec->initial_bmpc[i]), ftmp);
          mpf_get_rdpe (cdpe_Re (sec->bdpc[i]), ftmp);

          mpf_set_q (ftmp, sec->initial_bmpqic[i]);
          mpf_set (mpc_Im (sec->initial_bmpc[i]), ftmp);
          mpf_get_rdpe (cdpe_Im (sec->bdpc[i]), ftmp);
        }
    }

  /* Copy coefficients back in other places */
  for (i = 0; i < MPS_POLYNOMIAL (sec)->degree; i++)
    {
      /* Bulk copy of the MP coefficients */
      mpc_set (sec->ampc[i], sec->initial_ampc[i]);
      mpc_set (sec->bmpc[i], sec->initial_bmpc[i]);

      /* CDPE coefficients */
      mpc_get_cdpe (sec->adpc[i], sec->initial_ampc[i]);
      mpc_get_cdpe (sec->bdpc[i], sec->initial_bmpc[i]);

      /* Get floating points coefficients */
      cdpe_get_x (sec->afpc[i], sec->adpc[i]);
      cdpe_get_x (sec->bfpc[i], sec->bdpc[i]);

      /* Store mouduli of the coefficients */
      cdpe_mod (sec->aadpc[i], sec->adpc[i]);
      cdpe_mod (sec->abdpc[i], sec->bdpc[i]);
      sec->aafpc[i] = cplx_mod (sec->afpc[i]);
      sec->abfpc[i] = cplx_mod (sec->bfpc[i]);
    }

  mps_context_set_input_poly (s, MPS_POLYNOMIAL (sec));

  mpf_clear (ftmp);
}

void
mps_chebyshev_poly_read_from_stream (mps_context * ctx, mps_input_buffer * buffer, 
  mps_structure structure, mps_density density)
{
  int i, degree = -1;
  char * token;
  mps_chebyshev_poly * cpoly = mps_chebyshev_poly_new (ctx, ctx->n, structure);

  switch (density)
    {
      case MPS_DENSITY_DENSE:
        for (i = 0; i <= ctx->n; i++) 
          {
            if (MPS_STRUCTURE_IS_FP (structure))
              {
                token = mps_input_buffer_next_token (buffer);

                if (!token || (mpf_set_str (mpc_Re (cpoly->mfpc[i]), token, 10) != 0)) 
                    mps_error (ctx, 1, "Error while reading real part of coefficient %d", i);
                else
                  free (token);

                if (MPS_STRUCTURE_IS_COMPLEX (structure))
                  {
                    token = mps_input_buffer_next_token (buffer);

                    if (!token || (mpf_set_str (mpc_Im (cpoly->mfpc[i]), token, 10) != 0)) 
                        mps_error (ctx, 1, "Error while reading imaginary part of coefficient %d", i);
                    else
                      free (token);
                  }
                }
              else if (MPS_STRUCTURE_IS_RATIONAL (structure) || MPS_STRUCTURE_IS_INTEGER (structure))
                {
                  token = mps_input_buffer_next_token (buffer);

                  if (!token || (mpq_set_str (cpoly->rational_real_coeffs[i], token, 10)))
                    mps_error (ctx, 1, "Error while reading the real part of coefficient %d", i);
                  else
                    free (token);

                  if (MPS_STRUCTURE_IS_COMPLEX (structure))
                    {
                      token = mps_input_buffer_next_token (buffer);

                      if (!token || (mpq_set_str (cpoly->rational_imag_coeffs[i], token, 10)))
                        mps_error (ctx, 1, "Error while reading the imaginary part of coefficient %d", i);
                      else
                        free (token);
                    }

                    mpf_set_q (mpc_Re (cpoly->mfpc[i]), cpoly->rational_real_coeffs[i]);
                    mpf_set_q (mpc_Im (cpoly->mfpc[i]), cpoly->rational_imag_coeffs[i]);
                }

            /* Update other floating point coefficients */
            mpc_get_cdpe (cpoly->dpc[i], cpoly->mfpc[i]);
            mpc_get_cplx (cpoly->fpc[i], cpoly->mfpc[i]);

            if (ctx->debug_level & MPS_DEBUG_IO) {
              MPS_DEBUG_CPLX (ctx, cpoly->fpc[i], "Coefficient %d", i);
            }
          }
      break;

      case MPS_DENSITY_SPARSE:
        /* Set all the coefficients to zero first, so whatever is not given
         * will be assumed to be null. */
        for (i = 0; i < ctx->n; i++)
          {
            mpc_set_ui (cpoly->mfpc[i], 0U, 0U);
            mpq_set_ui (cpoly->rational_real_coeffs[i], 0U, 1U);
            mpq_set_ui (cpoly->rational_imag_coeffs[i], 0U, 1U);
          }

        /* Read the degree of the coefficient */
        token = mps_input_buffer_next_token (buffer);
        if (!token || !sscanf (token, "%d", &degree))
          mps_error (ctx, 1, "Cannot parse the degree of the coefficient.");
        else
          free (token);

        if (MPS_STRUCTURE_IS_FP (structure))
          {
            token = mps_input_buffer_next_token (buffer);

            if (!token || (mpf_set_str (mpc_Re (cpoly->mfpc[degree]), token, 10) != 0)) 
                mps_error (ctx, 1, "Error while reading real part of coefficient %d", degree);
            else
              free (token);

            if (MPS_STRUCTURE_IS_COMPLEX (structure))
              {
                token = mps_input_buffer_next_token (buffer);

                if (!token || (mpf_set_str (mpc_Im (cpoly->mfpc[degree]), token, 10) != 0)) 
                    mps_error (ctx, 1, "Error while reading imaginary part of coefficient %d", degree);
                else
                  free (token);                
              }
          }
        else
          {
            token = mps_input_buffer_next_token (buffer);

            if (!token || (mpq_set_str (cpoly->rational_real_coeffs[degree], token, 10)))
              mps_error (ctx, 1, "Error while reading the real part of coefficient %d", i);
            else
              free (token);

            if (MPS_STRUCTURE_IS_COMPLEX (structure))
              {
                token = mps_input_buffer_next_token (buffer);

                if (!token || (mpq_set_str (cpoly->rational_imag_coeffs[degree], token, 10)))
                  mps_error (ctx, 1, "Error while reading the imaginary part of coefficient %d", i);
                else
                  free (token);
              }
          }

        mpc_get_cdpe (cpoly->dpc[degree], cpoly->mfpc[degree]);
        mpc_get_cplx (cpoly->fpc[degree], cpoly->mfpc[degree]);
      break;

      default:
        mps_error (ctx, 1, "Only MPS_DENSITY_DENSE and MPS_DENSITY_SPARSE are supported in Chebyshev polynomials.");
        return;
        break;
    }

  mps_context_set_input_poly (ctx, MPS_POLYNOMIAL (cpoly));
}

void
mps_parse_stream_old (mps_context * s, mps_input_buffer * buffer)
{
  int i;
  mps_monomial_poly *poly;
  char data_type[3];
  char *token;
  mpf_t ftmp;
  mpq_t qtmp;

  mps_density density = MPS_DENSITY_DENSE;
  mps_structure structure = MPS_STRUCTURE_COMPLEX_FP;
  long int prec = 0;

  mpq_init (qtmp);
  mpf_init (ftmp);
  
  /* Here we have the data_type in the input_buffer, since the first line has been read, or at least
   * that should be the case. */
  token = mps_input_buffer_next_token (buffer);
  if (!token || !sscanf (token, "%3s", data_type))
    {
      mps_error (s, 1, "Error parsing the input file");
      return;
    }
  free (token);

  /* Parse data type converting it to the new format */
  switch (data_type[0])
    {
    case 's':
      density = MPS_DENSITY_SPARSE;
      break;
    case 'd':
      density = MPS_DENSITY_DENSE;
      break;
    case 'u':
      density = MPS_DENSITY_USER;
      break;
    default:
      mps_error (s, 1, "Found unsupported data_type in input file");
      return;
      break;
    }

  switch (data_type[1])
    {
    case 'r':
      structure = MPS_STRUCTURE_REAL_FP;
      break;
    case 'c':
      structure = MPS_STRUCTURE_COMPLEX_FP;
      break;
    default:
      mps_error (s, 1, "Found unsupported data_structure in input file");
      return;
      break;
    }

  switch (data_type[2])
    {
    case 'q':
      if (MPS_STRUCTURE_IS_REAL (structure))
        structure = MPS_STRUCTURE_REAL_RATIONAL;
      else if (MPS_STRUCTURE_IS_COMPLEX (structure))
        structure = MPS_STRUCTURE_COMPLEX_RATIONAL;
      break;
    case 'i':
      if (MPS_STRUCTURE_IS_REAL (structure))
        structure = MPS_STRUCTURE_REAL_INTEGER;
      else if (MPS_STRUCTURE_IS_COMPLEX (structure))
        structure = MPS_STRUCTURE_COMPLEX_INTEGER;
      break;
    case 'f':
      if (MPS_STRUCTURE_IS_REAL (structure))
        structure = MPS_STRUCTURE_REAL_FP;
      else if (MPS_STRUCTURE_IS_COMPLEX (structure))
        structure = MPS_STRUCTURE_COMPLEX_FP;
      break;      
    default:
      mps_error (s, 1, "Found unsupported data structure in input file");
      return;
      break;
    }

  /* Read precision and degree */
  prec = 0;
  token = mps_input_buffer_next_token (buffer);
  if (!token || !sscanf (token, "%ld", &prec))
    {
      mps_error (s, 1, "Error while reading the input precision of the coefficients");
      return;
    }
  else 
    prec *= LOG2_10; 
  free (token);

  token = mps_input_buffer_next_token (buffer);
  if (!token || !sscanf (token, "%d", &s->n))
    {
      mps_error (s, 1, "Error reading the degree of the polynomial");
      return;
    }
  free (token);
  s->deg = s->n;

  /* Hook up a compatiblity layer with the older MPSolve versions. Since it was possibile
   * to create custom Newton iteration routines and recomm*/
  if (density == MPS_DENSITY_USER)
    {
      mps_polynomial * user_poly = mps_polynomial_new (s);

      user_poly->density = MPS_DENSITY_USER;
      user_poly->structure = MPS_STRUCTURE_REAL_INTEGER;

      /* Newton iteration */
      user_poly->fnewton = mps_fnewton_usr;
      user_poly->dnewton = mps_dnewton_usr;
      user_poly->mnewton = mps_mnewton_usr;

      /* General disposition */
      user_poly->fstart = mps_general_fstart;
      user_poly->dstart = mps_general_dstart;
      user_poly->mstart = mps_general_mstart;

      /* Evaluation */
      user_poly->feval = mps_feval_usr;
      user_poly->deval = mps_deval_usr;
      user_poly->meval = mps_meval_usr;

      MPS_POLYNOMIAL (user_poly)->degree = s->n;

      mps_context_set_input_poly (s, user_poly);
      return;
    }

  /* Allocate the polynomial */
  poly = mps_monomial_poly_new (s, s->n);

  MPS_POLYNOMIAL (poly)->structure = structure;
  MPS_POLYNOMIAL (poly)->density = density;

  /* We still do not support sparse input */
  for (i = 0; i <= s->n; ++i)
      poly->spar[i] = true;

  /* Dense parsing */
  if (MPS_DENSITY_IS_DENSE (density))
    {
      if (MPS_STRUCTURE_IS_FP (structure))
        {
          for (i = 0; i < s->n + 1; ++i)
            {
              token = mps_input_buffer_next_token (buffer);
              if (!token || (mpf_set_str (mpc_Re (poly->mfpc[i]), token, 10) != 0))
                {
                  mps_raise_parsing_error (s, buffer, token, "Error parsing coefficients of the polynomial");
                  return;
                }
              free (token);

              if (MPS_STRUCTURE_IS_COMPLEX (structure))
                {
                  token = mps_input_buffer_next_token (buffer);
                  if (!token || (mpf_set_str (mpc_Im (poly->mfpc[i]), token, 10) != 0))
                    {
                      mps_raise_parsing_error (s, buffer, token, "Error parsing coefficients of the polynomial");
                      return;
                    }
                  free (token);
                }
              else
                mpf_set_ui (mpc_Im (poly->mfpc[i]), 0U);
            }
        }
      else if (MPS_STRUCTURE_IS_INTEGER (structure))
        {
          for (i = 0; i < s->n + 1; ++i)
            {
              token = mps_input_buffer_next_token (buffer);
              if (!token || (mpq_set_str (poly->initial_mqp_r[i], token, 10) != 0))
                {
                  mps_raise_parsing_error (s, buffer, token, "Error parsing coefficients of the polynomial");
                  return;
                }
              mpq_canonicalize (poly->initial_mqp_r[i]);
              free (token);
      
              if (MPS_STRUCTURE_IS_COMPLEX (structure))
                {
                  token = mps_input_buffer_next_token (buffer);
                  if (!token || (mpq_set_str (poly->initial_mqp_i[i], token, 10) != 0))
                    {
                      mps_raise_parsing_error (s, buffer, token, "Error parsing coefficients of the polynomial");
                      return;
                    }
                  mpq_canonicalize (poly->initial_mqp_i[i]);
                  free (token);
                }
              else
                mpq_set_ui (poly->initial_mqp_i[i], 0U, 0U);

              /* Copy coefficients in the floating point ones */
              mpf_set_q (mpc_Re (poly->mfpc[i]), poly->initial_mqp_r[i]);
              mpf_set_q (mpc_Im (poly->mfpc[i]), poly->initial_mqp_i[i]);
            }
        }
      else if (MPS_STRUCTURE_IS_RATIONAL (structure))
        {
          /* The old MPSolve format for the rational input is not understood
           * by GMP that expect rational to be represented as n / d, and here
           * we have two separate tokens n d
           */
          for (i = 0; i <= s->n; ++i)
            {
              /* Numerator of the real part of the coefficient */
              token = mps_input_buffer_next_token (buffer);
              if (!token || (mpq_set_str (qtmp, token, 10)) != 0)
                {
                  mps_raise_parsing_error (s, buffer, token, "Error parsing the numerator of a coefficient");
                  return;
                }
              mpq_set (poly->initial_mqp_r[i], qtmp);
              free (token);

              /* Denominator of the real part of the coefficient */
              token = mps_input_buffer_next_token (buffer);
              if (!token || (mpq_set_str (qtmp, token, 10)) != 0)
                {
                  mps_raise_parsing_error (s, buffer, token, "Error parsing the denominator of a coefficient");
                  return;
                }
              free (token);

              mpq_div (poly->initial_mqp_r[i], poly->initial_mqp_r[i], qtmp);
              mpq_canonicalize (poly->initial_mqp_r[i]);

              if (MPS_STRUCTURE_IS_COMPLEX (structure))
                {
                  /* Numerator of the real part of the coefficient */
                  token = mps_input_buffer_next_token (buffer);
                  if (!token || (mpq_set_str (qtmp, token, 10)) != 0)
                    {
                      mps_raise_parsing_error (s, buffer, token, "Error parsing the numerator of a coefficient");
                      return;
                    }
                  mpq_set (poly->initial_mqp_i[i], qtmp);
                  free (token);

                  /* Denominator of the real part of the coefficient */
                  token = mps_input_buffer_next_token (buffer);
                  if (!token || (mpq_set_str (qtmp, token, 10)) != 0)
                    {
                      mps_raise_parsing_error (s, buffer, token, "Error parsing the denominator of a coefficient");
                      return;
                    }
                  free (token);

                  mpq_div (poly->initial_mqp_i[i], poly->initial_mqp_i[i], qtmp);
                  mpq_canonicalize (poly->initial_mqp_i[i]);
                }
              else
                mpq_set_ui (poly->initial_mqp_i[i], 0U, 0U);
            }     
        }
    } /* closes if (MPS_INPUT_CONFIG_IS_DENSE (s->input_config)) */
  else if (MPS_DENSITY_IS_SPARSE (density))
    {
      /* There is another number in the config file that
       * represents the number of coefficients of the polynomial, so let's
       * read it and ignore it. */
      token = mps_input_buffer_next_token (buffer);
      free (token);

      /* Set all the spar to false, since we have still not read
       * any coefficient */
      for (i = 0; i <= s->n; ++i)
        poly->spar[i] = false;

      while ((token = mps_input_buffer_next_token (buffer)) != NULL)
        {
          /* Read the index from the buffer */
          if (!sscanf (token, "%d", &i))
            {
              mps_raise_parsing_error (s, buffer, token, "Error while parsing the degree of a monomial");
              free (token);
              return;
            }

          if (i < 0 || i > s->n) 
            {
              mps_raise_parsing_error (s, buffer, token, "Degree of coefficient out of bounds");
              free (token);
              return;
            }

          if (poly->spar[i])
            {
              mps_raise_parsing_error (s, buffer, token, "A monomial of the same degree has been inserted twice");
              return;
            }
          else
            poly->spar[i] = true;
          free (token);

          if (MPS_STRUCTURE_IS_FP (structure))
            {
              token = mps_input_buffer_next_token (buffer);
              if (!token || (mpf_set_str (mpc_Re (poly->mfpc[i]), token, 10) != 0))
                {
                  mps_raise_parsing_error (s, buffer, token, "Error parsing coefficients of the polynomial");
                  return;
                }
              free (token);
          
              if (MPS_STRUCTURE_IS_COMPLEX (structure))
                {
                  token = mps_input_buffer_next_token (buffer);
                  if (!token || (mpf_set_str (mpc_Im (poly->mfpc[i]), token, 10) != 0))
                    {
                      mps_raise_parsing_error (s, buffer, token, "Error parsing coefficients of the polynomial");
                      return;
                    }
                  free (token);
                }
              else
                mpf_set_ui (mpc_Im (poly->mfpc[i]), 0U);
            }
          else if (MPS_STRUCTURE_IS_INTEGER (structure))
            {
              token = mps_input_buffer_next_token (buffer);
              if (!token || (mpq_set_str (poly->initial_mqp_r[i], token, 10) != 0))
                {
                  mps_raise_parsing_error (s, buffer, token, "Error parsing coefficients of the polynomial");
                  return;
                }
              mpq_canonicalize (poly->initial_mqp_r[i]);
              free (token);
      
              if (MPS_STRUCTURE_IS_COMPLEX (structure))
                {
                  token = mps_input_buffer_next_token (buffer);
                  if (!token || (mpq_set_str (poly->initial_mqp_i[i], token, 10) != 0))
                    {
                      mps_raise_parsing_error (s, buffer, token, "Error parsing coefficients of the polynomial");
                      return;
                    }
                  mpq_canonicalize (poly->initial_mqp_i[i]);
                  free (token);
                }
              else
                mpq_set_ui (poly->initial_mqp_i[i], 0U, 0U);

              /* Copy coefficients in the floating point ones */
              mpf_set_q (mpc_Re (poly->mfpc[i]), poly->initial_mqp_r[i]);
              mpf_set_q (mpc_Im (poly->mfpc[i]), poly->initial_mqp_i[i]);
            }
          else if (MPS_STRUCTURE_IS_RATIONAL (structure))
            {
              /* Numerator of the real part of the coefficient */
              token = mps_input_buffer_next_token (buffer);
              if (!token || (mpq_set_str (qtmp, token, 10)) != 0)
                {
                  mps_raise_parsing_error (s, buffer, token, "Error parsing the numerator of a coefficient");
                  return;
                }
              mpq_set (poly->initial_mqp_r[i], qtmp);
              free (token);

              /* Denominator of the real part of the coefficient */
              token = mps_input_buffer_next_token (buffer);
              if (!token || (mpq_set_str (qtmp, token, 10)) != 0)
                {
                  mps_raise_parsing_error (s, buffer, token, "Error parsing the denominator of a coefficient");
                  return;
                }
              free (token);

              mpq_div (poly->initial_mqp_r[i], poly->initial_mqp_r[i], qtmp);
              mpq_canonicalize (poly->initial_mqp_r[i]);

              if (MPS_STRUCTURE_IS_COMPLEX (structure))
                {
                  /* Numerator of the real part of the coefficient */
                  token = mps_input_buffer_next_token (buffer);
                  if (!token || (mpq_set_str (qtmp, token, 10)) != 0)
                    {
                      mps_raise_parsing_error (s, buffer, token, "Error parsing the numerator of a coefficient");
                      return;
                    }
                  mpq_set (poly->initial_mqp_i[i], qtmp);
                  free (token);

                  /* Denominator of the real part of the coefficient */
                  token = mps_input_buffer_next_token (buffer);
                  if (!token || (mpq_set_str (qtmp, token, 10)) != 0)
                    {
                      mps_raise_parsing_error (s, buffer, token, "Error parsing the denominator of a coefficient");
                      return;
                    }
                  free (token);

                  mpq_div (poly->initial_mqp_i[i], poly->initial_mqp_i[i], qtmp);
                  mpq_canonicalize (poly->initial_mqp_i[i]);
                }
              else
                mpq_set_ui (poly->initial_mqp_i[i], 0U, 0U);
            }
        }
    } /* closes if (MPS_INPUT_CONFIG_IS_SPARSE (s->input_config)) */

  /* Copy coefficients back in other places */
  for (i = 0; i <= s->n; ++i)
    {
      if (poly->spar[i])
        {
           if (MPS_STRUCTURE_IS_INTEGER (structure) || 
               MPS_STRUCTURE_IS_RATIONAL (structure)) 
             { 
               mpf_set_q (mpc_Re (poly->mfpc[i]), poly->initial_mqp_r[i]); 
               mpf_set_q (mpc_Im (poly->mfpc[i]), poly->initial_mqp_i[i]); 
             }

          mpc_get_cplx (poly->fpc[i], poly->mfpc[i]);
          mpc_get_cdpe (poly->dpc[i], poly->mfpc[i]);

          /* Compute modules of coefficients */
          cdpe_mod (poly->dap[i], poly->dpc[i]);
          poly->fap[i] = rdpe_get_d (poly->dap[i]);

          if (MPS_STRUCTURE_IS_FP (structure))
            mpf_set (poly->mfpr[i], mpc_Re (poly->mfpc[i]));

          if (i > 0)
            mpc_mul_ui (poly->mfppc[i-1], poly->mfppc[i], i);

          if (s->debug_level & MPS_DEBUG_IO)
            {
              MPS_DEBUG_MPC (s, 15, poly->mfpc[i], "s->mfpc[%d]", i);
            }
        }
      else
        {
          cplx_set (poly->fpc[i], cplx_zero);
          cdpe_set (poly->dpc[i], cdpe_zero);
          
          rdpe_set (poly->dap[i], rdpe_zero);
          poly->fap[i] = 0.0f;

          if (MPS_STRUCTURE_IS_FP (structure))
            mpf_set (poly->mfpr[i], mpc_Re (poly->mfpc[i]));
        }

      if (s->debug_level & MPS_DEBUG_IO)
        {
          MPS_DEBUG_RDPE (s, poly->dap[i], "poly->dap[%d]", i);
          MPS_DEBUG (s, "poly->fap[%d] = %e", i, poly->fap[i]);
        }
    }

  poly->prec = prec;

  mps_context_set_input_poly (s, MPS_POLYNOMIAL (poly));

  mpf_clear (ftmp);
  mpq_clear (qtmp);
}


/**
 * @brief Parse a stream for input data.
 */
void
mps_parse_stream (mps_context * s, FILE * input_stream)
{
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
              mps_parse_stream_old (s, buffer);
              mps_input_buffer_free (buffer);
              return;
            }
          parsing_options = false;
        }
      else
        {
          first_pass = false;
          input_option =
            mps_parse_option_line (s, line, strlen (line));

          /* Parsing of the degree */
          if (input_option.flag == MPS_KEY_DEGREE)
            {
              s->n = atoi (input_option.value);
              if (s->n <= 0)
                mps_error (s, 1, "Degree must be a positive integer");
            }

          /* Parsing precision of input coefficients */
          if (input_option.flag == MPS_KEY_PRECISION)
            {
              mps_context_set_input_prec (s, atoi (input_option.value) * LOG2_10);
              if (input_precision <= 0)
                mps_error (s, 1, "Precision must be a positive integer");
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
    mps_error (s, 1,
               "Degree of the polynomial must be provided via the Degree=%d configuration option.");
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
      mps_secular_equation_read_from_stream (s, buffer, structure, density);
      break;

    case MPS_REPRESENTATION_CHEBYSHEV:
      if (s->debug_level & MPS_DEBUG_IO)
        MPS_DEBUG (s, "Parsing mps_chebyshev_poly from stream");
      mps_chebyshev_poly_read_from_stream (s, buffer, structure, density);
      break;

    case MPS_REPRESENTATION_MONOMIAL:
    default:
      if (s->debug_level & MPS_DEBUG_IO)
        MPS_DEBUG (s, "Parsing mps_polynomial from stream");
      mps_monomial_poly_read_from_stream (s, buffer, structure, density);
      break;
  }

  MPS_POLYNOMIAL (s->active_poly)->structure = structure;
  MPS_POLYNOMIAL (s->active_poly)->density = density;

  mps_context_set_input_prec (s, input_precision);

  mps_input_buffer_free (buffer);
}



/*********************************************************
*      SUBROUTINE READROOTS                              *
*********************************************************/
void
mps_readroots (mps_context * s)
{
  long digits;
  int i, read_elements;

  if (s->DOLOG)
    fprintf (s->logstr, "Reading roots...\n");

  read_elements = fscanf (s->rtstr, "%ld", &digits);
  if (!read_elements)
    {
      mps_error (s, 1, "Error while reading roots, aborting.");
    }

  /* precision setup code goes here */

  for (i = 0; i < s->n; i++)
    mpc_inp_str_u (s->root[i]->mvalue, s->rtstr, 10);
}

/*********************************************************
*      SUBROUTINE COUNTROOTS                             *
*********************************************************/
void
mps_countroots (mps_context * s)
{
  int k;

  if (s->DOLOG)
    fprintf (s->logstr, "Counting roots...\n");

  s->count[0] = s->count[1] = s->count[2] = 0;

  for (k = 0; k < s->n; k++)
    switch (s->root[k]->inclusion)
      {
      case MPS_ROOT_INCLUSION_IN:
        s->count[0]++;
        break;
      case MPS_ROOT_INCLUSION_OUT:
        s->count[1]++;
        break;
      default:
        s->count[2]++;
        break;
      }

  if (s->output_config->search_set == MPS_SEARCH_SET_UNITARY_DISC_COMPL)
    s->count[1] += s->zero_roots;
  else
    s->count[0] += s->zero_roots;
}

/*********************************************************
*      SUBROUTINE OUTCOUNT                               *
*********************************************************/
void
mps_outcount (mps_context * s)
{
  mps_countroots (s);

  fprintf (s->outstr, "%d roots are inside;\n", s->count[0]);
  fprintf (s->outstr, "%d roots are outside;\n", s->count[1]);
  fprintf (s->outstr, "%d roots are uncertain.\n", s->count[2]);
  if (s->DOLOG)
    {
      fprintf (s->logstr, "%d roots are inside;\n", s->count[0]);
      fprintf (s->logstr, "%d roots are outside;\n", s->count[1]);
      fprintf (s->logstr, "%d roots are uncertain.\n", s->count[2]);
    }
}

/*********************************************************
*      SUBROUTINE OUTFLOAT                               *
*********************************************************/
void
mps_outfloat (mps_context * s, mpf_t f, rdpe_t rad, long out_digit,
              mps_boolean sign)
{
  mpf_t t;
  rdpe_t r, ro;
  double d;
  long l, digit, true_digit;

  if (s->output_config->format == MPS_OUTPUT_FORMAT_FULL)
    {
      mpf_init2 (t, mpf_get_prec (f));
      mpf_set (t, f);
      mpf_out_str (s->outstr, 10, 0, t);
      mpf_clear (t);
      return;
    }

  mpf_init2 (t, s->output_config->prec);

  mpf_get_rdpe (ro, f);
  if (s->output_config->format == MPS_OUTPUT_FORMAT_GNUPLOT ||
      s->output_config->format == MPS_OUTPUT_FORMAT_GNUPLOT_FULL)
    rdpe_out_str_u (s->outstr, ro);
  else
    {
      rdpe_abs_eq (ro);
      if (rdpe_ne (ro, rdpe_zero))
        rdpe_div (r, rad, ro);
      else
        rdpe_set_d (r, 1.0e-10);
      digit = (long) (-rdpe_log10 (r) - 0.5);
      if (digit <= 0)
        {
          rdpe_get_dl (&d, &l, ro);
          fprintf (s->outstr, "0.e%ld", l);
        }
      else
        {
          true_digit = (long) (LOG10_2 * mpf_get_prec (f));
          true_digit = MIN (digit, true_digit);
          true_digit = MIN (true_digit, out_digit);
          if (sign)
            mpf_set (t, f);
          else
            mpf_abs (t, f);
          mpf_out_str (s->outstr, 10, true_digit, t);
        }
    }

  mpf_clear (t);
}

/*********************************************************
*      SUBROUTINE OUTROOT                                *
*********************************************************/
void
mps_outroot (mps_context * s, int i, int num)
{
  long out_digit;

  out_digit = (long) (LOG10_2 * s->output_config->prec) + 10;

  /* print format header */
  switch (s->output_config->format)
    {
    case MPS_OUTPUT_FORMAT_COMPACT:
    case MPS_OUTPUT_FORMAT_FULL:
      fprintf (s->outstr, "(");
      break;
    case MPS_OUTPUT_FORMAT_VERBOSE:
      fprintf (s->outstr, "Root(%d) = ", num);
      break;
    default:
      break;
    }

  /* print real part */
  if (i == ISZERO || s->root[i]->attrs == MPS_ROOT_ATTRS_IMAG)
    fprintf (s->outstr, "0");
  else
    mps_outfloat (s, mpc_Re (s->root[i]->mvalue), s->root[i]->drad, out_digit, true);

  /* print format middle part */
  switch (s->output_config->format)
    {
    case MPS_OUTPUT_FORMAT_BARE:
      fprintf (s->outstr, " ");
      break;
    case MPS_OUTPUT_FORMAT_GNUPLOT:
    case MPS_OUTPUT_FORMAT_GNUPLOT_FULL:
      fprintf (s->outstr, "\t");
      break;
    case MPS_OUTPUT_FORMAT_COMPACT:
    case MPS_OUTPUT_FORMAT_FULL:
      fprintf (s->outstr, ", ");
      break;
    case MPS_OUTPUT_FORMAT_VERBOSE:
      if (i == ISZERO || mpf_sgn (mpc_Im (s->root[i]->mvalue)) >= 0)
        fprintf (s->outstr, " + I * ");
      else
        fprintf (s->outstr, " - I * ");
      break;
    default:
      break;
    }

  /* print imaginary part */
  if (i == ISZERO || s->root[i]->attrs == MPS_ROOT_ATTRS_REAL)
    fprintf (s->outstr, "0");
  else
    mps_outfloat (s, mpc_Im (s->root[i]->mvalue), s->root[i]->drad, out_digit,
                  s->output_config->format != MPS_OUTPUT_FORMAT_VERBOSE);

  /* If the output format is GNUPLOT_FORMAT_FULL, print out also the radius */
  if (s->output_config->format == MPS_OUTPUT_FORMAT_GNUPLOT_FULL)
    {
      fprintf (s->outstr, "\t");
      rdpe_out_str_u (s->outstr, s->root[i]->drad);
      fprintf (s->outstr, "\t");
      rdpe_out_str_u (s->outstr, s->root[i]->drad);
    }

  /* print format ending */
  switch (s->output_config->format)
    {
    case MPS_OUTPUT_FORMAT_COMPACT:
      fprintf (s->outstr, ")");
      break;
    case MPS_OUTPUT_FORMAT_FULL:
      fprintf (s->outstr, ")\n");
      if (i != ISZERO)
        {
          rdpe_outln_str (s->outstr, s->root[i]->drad);
          fprintf (s->outstr, "Status: %s, %s, %s\n", 
                   MPS_ROOT_STATUS_TO_STRING (s->root[i]->status),
                   MPS_ROOT_ATTRS_TO_STRING (s->root[i]->attrs),
                   MPS_ROOT_INCLUSION_TO_STRING (s->root[i]->inclusion));
        }
      else
        fprintf (s->outstr, " 0\n ---\n");
      break;
    default:
      break;
    }
  fprintf (s->outstr, "\n");

  /* debug info */
  if (s->DOLOG)
    {
      if (i == ISZERO)
        fprintf (s->logstr, "zero root %-4d = 0", num);
      else
        {
          fprintf (s->logstr, "Root %-4d = ", i);
          mpc_out_str_2 (s->logstr, 10, 0, 0, s->root[i]->mvalue);
          fprintf (s->logstr, "\n");
          fprintf (s->logstr, "  Radius = ");
          rdpe_outln_str (s->logstr, s->root[i]->drad);
          fprintf (s->logstr, "  Prec = %ld\n",
                   (long) (mpc_get_prec (s->root[i]->mvalue) / LOG2_10));
          fprintf (s->logstr, "  Approximation = %s\n", 
                   MPS_ROOT_STATUS_TO_STRING (s->root[i]->status));
          fprintf (s->logstr, "  Attributes = %s\n",
                   MPS_ROOT_ATTRS_TO_STRING (s->root[i]->attrs));
          fprintf (s->logstr, "  Inclusion = %s\n",
                   MPS_ROOT_INCLUSION_TO_STRING (s->root[i]->inclusion));
          fprintf (s->logstr, "--------------------\n");
        }
    }
}

/*********************************************************
*      SUBROUTINE OUTPUT                                 *
*********************************************************/
void
mps_output (mps_context * s)
{
  int i, ind, num = 0;

  if (s->DOLOG)
    fprintf (s->logstr, "--------------------\n");

  if (s->output_config->format != MPS_OUTPUT_FORMAT_GNUPLOT && 
      s->output_config->format != MPS_OUTPUT_FORMAT_GNUPLOT_FULL)
    {
      if (s->over_max)
        {
          mps_warn (s, "Warning: Input precision has been reached during computation, "
                    "so not all the required digits may have been computed.");
        }
    }

  /* Start with plotting instructions in the case of 
   * MPS_OUTPUT_GNUPLOT_FULL, so the output can be
   * piped directly to gnuplot */
  if (s->output_config->format == MPS_OUTPUT_FORMAT_GNUPLOT_FULL)
    {
      fprintf (s->outstr, "# MPSolve output for GNUPLOT\n");
      fprintf (s->outstr, "# Make user that this output is piped into gnuplot using a command like\n");
      fprintf (s->outstr, "# mpsolve -Ogf | gnuplot \n");
      fprintf (s->outstr, "set pointsize 0.3\n");
      fprintf (s->outstr, "plot '-' title 'Computed roots' with %s\n", s->gnuplot_format);
    }

  if (s->output_config->goal == MPS_OUTPUT_GOAL_COUNT)
    mps_outcount (s);
  else
    {
      if (s->output_config->search_set != MPS_SEARCH_SET_UNITARY_DISC_COMPL)
        for (i = 0; i < s->zero_roots; i++)
            mps_outroot (s, ISZERO, num++);
      for (ind = 0; ind < s->n; ind++)
        {
          i = s->order[ind];
          if (s->root[i]->inclusion == MPS_ROOT_INCLUSION_OUT)
            continue;
          mps_outroot (s, i, num++);
        }
    }

  if (s->output_config->format == MPS_OUTPUT_FORMAT_GNUPLOT_FULL)
    {
      fprintf (s->outstr, "e\n");
      fprintf (s->outstr, "pause mouse close\n");
      fprintf (s->outstr, "# End of MPSolve GNUPLOT output. If you are seeing this maybe\n");
      fprintf (s->outstr, "# you forgot to pipe the ***solve command into gnuplot?\n");
    }
}

/*********************************************************
*      SUBROUTINE COPY_ROOTS                             *
*********************************************************/
void
mps_copy_roots (mps_context * s)
{
  int i;

  MPS_DEBUG_THIS_CALL;

  switch (s->lastphase)
    {
    case no_phase:
      mps_error (s, 1, "Nothing to copy");
      break;

    case float_phase:
      if (s->DOSORT)
        mps_fsort (s);
      for (i = 0; i < s->n; i++)
        {
          mpc_set_prec (s->root[i]->mvalue, DBL_MANT_DIG);
          mpc_set_cplx (s->root[i]->mvalue, s->root[i]->fvalue);
          rdpe_set_d (s->root[i]->drad, s->root[i]->frad);
        }
      break;

    case dpe_phase:
      if (s->DOSORT)
        mps_dsort (s);
      for (i = 0; i < s->n; i++)
        {
          mpc_set_prec (s->root[i]->mvalue, DBL_MANT_DIG);
          mpc_set_cdpe (s->root[i]->mvalue, s->root[i]->dvalue);
        }
      break;

    case mp_phase:
      if (s->DOSORT)
        mps_msort (s);
      break;

    }
}

/*************************************************************
 *                     SUBROUTINE DUMP                       *
 *************************************************************/
void
mps_dump (mps_context * s)
{
  int i;
  FILE * dmpstr = s->logstr;

  MPS_DEBUG (s, "Dumping the approximations:");

  /* output current status */
  /* fprintf (dmpstr, */
  /*          "Phase=%d, In=%d, Out=%d, Uncertain=%d, Zero=%d, Clusters=%ld\n", */
  /*          s->lastphase, s->count[0], s->count[1], s->count[2], s->zero_roots, */
  /*          s->clusterization->n); */

  MPS_DEBUG (s, "Phase = %s, In = %d, Out = %d, Uncertain = %d, Zero = %d, Cluster = %ld",
             MPS_PHASE_TO_STRING (s->lastphase), s->count[0], s->count[1], s->count[2],
             s->zero_roots, s->clusterization->n);

  /* output current approximations */
  /* fprintf (dmpstr, "\nCurrent approximations:\n"); */
  MPS_DEBUG (s, "Current approximations:");
  for (i = 0; i < s->n; i++)
    {
      switch (s->lastphase)
        {
        case no_phase:
        case float_phase:
          MPS_DEBUG_CPLX (s, s->root[i]->fvalue, "Approximation  %4d", i);
          break;

        case dpe_phase:
          MPS_DEBUG_CDPE (s, s->root[i]->dvalue, "Approximation  %4d", i);
          break;

        case mp_phase:
          MPS_DEBUG_MPC (s, s->mpwp, s->root[i]->mvalue, "Approximation  %4d", i);
          break;
        }
    }

  /* output radii */
  MPS_DEBUG (s, "Current radii:");
  for (i = 0; i < s->n; i++)
    {
      switch (s->lastphase)
        {
        case no_phase:
        case float_phase:
          MPS_DEBUG (s, "Radius of root %4d = %e", i, s->root[i]->frad);
          break;

        case dpe_phase:
        case mp_phase:
          MPS_DEBUG_RDPE (s, s->root[i]->drad, "Radius of root %4d", i);
          break;
        }
    }

  MPS_DEBUG (s, " ");
  mps_dump_status (s, dmpstr);
}

/**
 * @brief Dump cluster structure to <code>outstr</code>.
 *
 * @param s the mps_context struct pointer.
 * @param outstr The output stream where the cluster structure
 *  will be dumped.
 */
void
mps_dump_cluster_structure (mps_context * s, FILE * outstr)
{
  fprintf (outstr,
           "    MPS_DUMP_CLUSTER_STRUCTURE: Dumping cluster structure\n");

  mps_cluster_item * cluster_item;
  mps_cluster * cluster;
  mps_root * root;

  for (cluster_item = s->clusterization->first; cluster_item != NULL;
       cluster_item = cluster_item->next)
    {
      cluster = cluster_item->cluster;
      fprintf (outstr, "     Cluster contains %ld roots:\n", cluster->n);

      /* Dump cluster roots, but not more than 15 for line, to make
       * the output readable. */
      int j = 0;
      for (root = cluster->first; root != NULL; root = root->next)
        {
          /* Go to a newlint if 15 roots are printed out */
          if (j % 15 == 0)
            {
              fprintf (outstr, "\n       ");
            }

          fprintf (outstr, " %4ld", root->k);
          j++;
        }

      /* Make space untile the next cluster */
      fprintf (outstr, "\n\n");
    }
}

/**
 * @brief Dump status of all the root approximations
 */
void
mps_dump_status (mps_context * s, FILE * outstr)
{
  int i;
  MPS_DEBUG (s, "              Approximation              Attributes       Inclusion");
  for (i = 0; i < s->n; i++)
    {
      MPS_DEBUG (s, "Status  %4d: %-25s  %-15s  %-15s", i,
                 MPS_ROOT_STATUS_TO_STRING (s->root[i]->status),
                 MPS_ROOT_ATTRS_TO_STRING  (s->root[i]->attrs), 
                 MPS_ROOT_INCLUSION_TO_STRING (s->root[i]->inclusion));
    }
}

/*************************************************************
 *                     SUBROUTINE WARN                       *
 *************************************************************/
void
mps_warn (mps_context * st, char *s)
{
  if (st->DOWARN)
    {
      if (s[strlen (s)] == '\n')
        {
          fprintf (stderr, "%s", s);
        }
      else
        {
          fprintf (stderr, "%s\n", s);
        }
    }
}

/**
 * @brief Check if the file descriptor associated to stream
 * is bounded to a tty.
 *
 * @param stream the stream to check
 */
mps_boolean
mps_is_a_tty (FILE * stream)
{
#ifndef __WINDOWS
  return isatty (stream->_fileno);
#else
  return _isatty (_fileno (stream));
#endif
}

/*************************************************************
 *                     SUBROUTINE VAERROR                    *
 *************************************************************/
void
mps_error (mps_context * s, int args, ...)
{
  va_list ap;
  char *token;

  va_start (ap, args);
  while (args--)
    {
      token = va_arg (ap, char *);

      if (s->last_error == NULL)
        {
          s->last_error = strdup (token);
        }
      else
        {
          s->last_error = mps_realloc (s->last_error, strlen (s->last_error) + strlen (token) + 2);
          s->last_error = strcat (s->last_error, "\n");
          s->last_error = strcat (s->last_error, token);
        }
    }
  va_end (ap);

  s->error_state = true;

  /* Dump approximations, but only if they are present */
  if (s->root && s->lastphase)
    mps_dump (s);  /* dump status          */
  /* exit (EXIT_FAILURE);          /\* exit program         *\/ */
}

void
mps_print_errors (mps_context * s)
{
  if (mps_is_a_tty (s->logstr))
    mps_warn (s, "\033[31;1m!\033[0m MPSolve encountered an error:");  /* output error message */
  else
    mps_warn (s, "! MPSolve encountered an error:");

  mps_warn (s, s->last_error);
}
