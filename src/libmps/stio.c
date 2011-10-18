/***********************************************************
**       Multiprecision Polynomial Solver (MPSolve)       **
**                 Version 2.2, May 2001                  **
**                                                        **
**                      Written by                        **
**       Dario Andrea Bini and Giuseppe Fiorentino        **
**       (bini@dm.unipi.it)  (fiorent@dm.unipi.it)        **
**                                                        **
** (C) 2001, Dipartimento di Matematica, FRISCO LTR 21024 **
***********************************************************/

/**
 * @file
 * @brief Routines for I/O in MPSolve
 */

#include <stdarg.h>
#include <string.h>
#include <mps/gmptools.h>
#include <mps/core.h>
#include <mps/secular.h>
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
mps_raise_parsing_error (mps_status * s, const char * token, 
			 const char * message)
{
  char * output = (char *) malloc (sizeof (char) * (strlen (token) + strlen ("Parsing error near the token: ") + 1));
  sprintf (output, "Parsing error near the token: %s", token);

  mps_error (s, 2, output, message);
  free (output);
}

/**
 * @brief Check if the string given is equal to the option
 * identified by the given string
 */
mps_boolean
mps_is_option (mps_status * s, const char *option_string1,
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
mps_parse_option_line (mps_status * s, char *line, size_t length)
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

  /* Parsing keys with values. If = is not found in the
   * input string, than an error has occurred so we should
   * return. */
  equal_position = strchr (option, '=');
  if (equal_position == NULL)
    {
      if (s->debug_level & MPS_DEBUG_IO)
        {
          MPS_DEBUG (s, "Parsed input_option.flag = %d", input_option.flag);
        }
      return input_option;
    }
  else
    {
      input_option.value = equal_position + 1;
      /* Make a copy of the option to parse it without
       * equal sign and anything after it */
      c_ptr = option;
      option = (char *) malloc (sizeof (char) * (strlen (option) + 1));
      strcpy (option, c_ptr);
      *strchr (option, '=') = '\0';
    }

  if (mps_is_option (s, option, "degree"))
    input_option.flag = MPS_KEY_DEGREE;

  /* Free the copy of the option */
  free (option);
  return input_option;
}

void
mps_secular_equation_read_from_stream_poly (mps_status * s,
					    mps_input_buffer * buffer)
{
  mps_secular_equation *sec;
  int i, r;
  mpf_t ftmp;
  char * token;

  mpf_init (ftmp);

  /* Read directly the secular equation in DPE, so we don't need
   * to have a fallback case if the coefficients are bigger than
   * what is supported by the standard floating point arithmetic */
  sec = mps_secular_equation_new_raw (s, s->n);

  /* Preallocate the data that we need */
  s->spar = mps_boolean_valloc (s->n + 2);
  s->mfpc = mpc_valloc (s->n + 1);
  s->dpc  = cdpe_valloc (s->n + 1);
  s->fpc  = cplx_valloc (s->n + 1);
  s->dap  = rdpe_valloc (s->n + 1);
  s->fap  = double_valloc (s->n + 1);
  mpc_vinit (s->mfpc, s->n + 1);

  /* We still do not support sparse input */
  for (i = 0; i < s->n; ++i)
      s->spar[i] = true;

  if (MPS_STRUCTURE_IS_FP (s->config))
    {
      for (i = 0; i < s->n + 1; ++i)
	{
	  token = mps_input_buffer_next_token (buffer);
	  if (!token || (mpf_set_str (mpc_Re (s->mfpc[i]), token, 10) != 0))
	    mps_raise_parsing_error (s, token, "Error parsing coefficients of the polynomial");
	  free (token);

	  if (MPS_STRUCTURE_IS_COMPLEX (s->config))
	    {
	      token = mps_input_buffer_next_token (buffer);
	      if (!token || (mpf_set_str (mpc_Im (s->mfpc[i]), token, 10) != 0))
		mps_raise_parsing_error (s, token, "Error parsing coefficients of the polynomial");
	      free (token);
	    }
	  else
	    mpf_set_ui (mpc_Im (s->mfpc[i]), 0U);
	}
    }
  else if (MPS_STRUCTURE_IS_RATIONAL (s->config) ||
	   MPS_STRUCTURE_IS_INTEGER (s->config))
    {
      for (i = 0; i < s->n + 1; ++i)
	{
	  token = mps_input_buffer_next_token (buffer);
	  if (!token || (mpq_set_str (sec->initial_bmpqrc[i], token, 10) != 0))
	    mps_raise_parsing_error (s, token, "Error parsing coefficients of the polynomial");
	  mpq_canonicalize (sec->initial_bmpqrc[i]);
	  free (token);
      
	  if (MPS_STRUCTURE_IS_COMPLEX (s->config))
	    {
	      token = mps_input_buffer_next_token (buffer);
	      if (!token || (mpq_set_str (sec->initial_bmpqic[i], token, 10) != 0))
		mps_raise_parsing_error (s, token, "Error parsing coefficients of the polynomial");
	      mpq_canonicalize (sec->initial_bmpqic[i]);
	      free (token);
	    }
	  else
	    mpq_set_ui (sec->initial_bmpqic[i], 0U, 0U);

	  /* Copy coefficients in the floating point ones */
	  mpf_set_q (mpc_Re (s->mfpc[i]), sec->initial_bmpqrc[i]);
	  mpf_set_q (mpc_Im (s->mfpc[i]), sec->initial_bmpqic[i]);
	}
      
    }

  /* Copy coefficients back in other places */
  for (i = 0; i < s->n + 1; ++i)
    {
      mpc_get_cplx (s->fpc[i], s->mfpc[i]);
      mpc_get_cdpe (s->dpc[i], s->mfpc[i]);

      /* Compute modules of coefficients */
      cdpe_mod (s->dap[i], s->dpc[i]);
      s->fap[i] = rdpe_get_d (s->dap[i]);
    }

  s->secular_equation = sec;
  mpf_clear (ftmp);
}

void
mps_secular_equation_read_from_stream (mps_status * s,
                                       mps_input_buffer * buffer)
{
  mps_secular_equation *sec;
  int i, r;
  mpf_t ftmp;
  char * token;

  mpf_init (ftmp);

  /* Read directly the secular equation in DPE, so we don't need
   * to have a fallback case if the coefficients are bigger than
   * what is supported by the standard floating point arithmetic */
  sec = mps_secular_equation_new_raw (s, s->n);

  /* Parsing of integers and floating point is done with Multiprecision */
  if (MPS_STRUCTURE_IS_FP (s->config))
    {
      for (i = 0; i < s->n; i++)
        {
	  token = mps_input_buffer_next_token (buffer);
          if (!token || (mpf_set_str (mpc_Re (sec->initial_ampc[i]), token, 10) != 0))
            {
              MPS_DEBUG (s,
                         "Error reading coefficient a[%d] of the secular equation (real part)",
                         i);
              mps_raise_parsing_error (s, token, 
                         "Error reading some coefficients of the secular equation.\n"
                         "Please check your input file.");
            }
	  free (token);

          /* Imaginary part, read only if the input is complex */
          if (MPS_STRUCTURE_IS_COMPLEX (s->config))
            {
	      token = mps_input_buffer_next_token (buffer);
              if (!token || (mpf_set_str (mpc_Im (sec->initial_ampc[i]), token, 10) != 0))
                {
                  MPS_DEBUG (s,
                             "Error reading coefficient a[%d] of the secular equation (imaginary part)",
                             i);
                  mps_raise_parsing_error (s, token,
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
              mps_raise_parsing_error (s, token,
                         "Error reading some coefficients of the secular equation.\n"
                         "Please check your input file.");
            }
	  free (token);

          /* Again, read the imaginary part only if the input is complex */
          if (MPS_STRUCTURE_IS_COMPLEX (s->config))
            {
	      token = mps_input_buffer_next_token (buffer);
	      if (!token || (mpf_set_str (mpc_Im (sec->initial_bmpc[i]), token, 10) != 0))
                {
                  MPS_DEBUG (s,
                             "Error reading coefficient b[%d] of the secular equation (imaginary part)",
                             i);
                  mps_raise_parsing_error (s, token,
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
  else if (MPS_STRUCTURE_IS_RATIONAL (s->config) ||
           MPS_STRUCTURE_IS_INTEGER  (s->config))
    {
      for (i = 0; i < s->n; i++)
        {
          /* Read real part of the a_i */
	  token = mps_input_buffer_next_token (buffer);
          if (!token || (mpq_set_str (sec->initial_ampqrc[i], token, 10) != 0))
	    {
	      MPS_DEBUG (s, "Error reading the coefficients a[%d] of the secular equation (real part)", i);
	      mps_raise_parsing_error (s, token, 
				       "Error reading some coefficients of the secular equation.\nPlease check your input file");
	    }
          mpq_canonicalize (sec->initial_ampqrc[i]);
	  free (token);

          /* Read imaginary part of the a_i */
          if (MPS_STRUCTURE_IS_COMPLEX (s->config))
            {
	      token = mps_input_buffer_next_token (buffer);
              if (!token || (mpq_set_str (sec->initial_ampqic[i], token, 10) != 0))
		{	      
		  MPS_DEBUG (s, "Error reading the coefficients a[%d] of the secular equation (imaginary part)", i);
		  mps_raise_parsing_error (s, token, 
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
	      mps_raise_parsing_error (s, token, 
				       "Error reading some coefficients of the secular equation.\nPlease check your input file");
	    }	    
          mpq_canonicalize (sec->initial_bmpqrc[i]);
	  free (token);

          /* Read imaginary part of the b_i */
          if (MPS_STRUCTURE_IS_COMPLEX (s->config))
            {
	      token = mps_input_buffer_next_token (buffer);
              if (!token || (mpq_set_str (sec->initial_bmpqic[i], token, 10) != 0))
		{	      
		  MPS_DEBUG (s, "Error reading the coefficients b[%d] of the secular equation (imaginary part)", i);
		  mps_raise_parsing_error (s, token, 
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
  for (i = 0; i < sec->n; i++)
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
    }

  s->secular_equation = sec;

  /* Deflate input, if identical b_i coefficients are found */
  mps_secular_deflate (s, sec);

  mpf_clear (ftmp);
}


/**
 * @brief Parse a stream for input data.
 */
void
mps_parse_stream (mps_status * s, FILE * input_stream,
                  mps_input_configuration * default_configuration)
{
  mps_boolean parsing_options = true;
  mps_input_buffer *buffer;
  mps_input_option input_option;
  int i;
  ssize_t length;
  char * line;

  /* Set default values for the parsing configuration */
  s->config = default_configuration;

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
      line = strdup (buffer->line);
      if (strchr (line, ';') == NULL || mps_input_buffer_eof (buffer))
        {
	  if (s->debug_level & MPS_DEBUG_IO)
	    {
	      MPS_DEBUG (s, "Finished parsing options");
	    }
          parsing_options = false;
        }
      else
        {
          input_option =
            mps_parse_option_line (s, line, strlen (line));

          if (s->debug_level & MPS_DEBUG_IO)
            {
              MPS_DEBUG (s, "Parsed option %d", input_option.flag);
            }

          /* Parsing of the degree */
          if (input_option.flag == MPS_KEY_DEGREE)
            {
              s->n = atoi (input_option.value);
              if (s->n <= 0)
                mps_error (s, 1, "Degree must be a positive integer");
            }

          /* Parsing of representations */
          else if (input_option.flag == MPS_FLAG_SECULAR)
            s->config->representation = MPS_REPRESENTATION_SECULAR;
	  else if (input_option.flag == MPS_FLAG_MONOMIAL)
	    s->config->representation = MPS_REPRESENTATION_MONOMIAL;

          /* Parsing of algebraic structure of the input */
          else if (input_option.flag == MPS_FLAG_REAL)
            {
              /* Switch on algebraic structure that is already set */
              if (MPS_STRUCTURE_IS_INTEGER (s->config))
                s->config->structure = MPS_STRUCTURE_REAL_INTEGER;
              else if (MPS_STRUCTURE_IS_RATIONAL (s->config))
                s->config->structure = MPS_STRUCTURE_REAL_RATIONAL;
              else if (MPS_STRUCTURE_IS_FP (s->config))
                s->config->structure = MPS_STRUCTURE_REAL_FP;
            }
          else if (input_option.flag == MPS_FLAG_COMPLEX)
            {
              /* Switch on algebraic structure that is already set */
              if (MPS_STRUCTURE_IS_INTEGER (s->config))
                s->config->structure = MPS_STRUCTURE_COMPLEX_INTEGER;
              else if (MPS_STRUCTURE_IS_RATIONAL (s->config))
                s->config->structure = MPS_STRUCTURE_COMPLEX_RATIONAL;
              else if (MPS_STRUCTURE_IS_FP (s->config))
                s->config->structure = MPS_STRUCTURE_COMPLEX_FP;
            }

          /* Parsing of algebraic structure of the input, i.e.
           * Integer, Rational or floating point */
          else if (input_option.flag == MPS_FLAG_INTEGER)
            {
              if (MPS_STRUCTURE_IS_REAL (s->config))
                s->config->structure = MPS_STRUCTURE_REAL_INTEGER;
              else if (MPS_STRUCTURE_IS_COMPLEX (s->config))
                s->config->structure = MPS_STRUCTURE_COMPLEX_INTEGER;
            }
          else if (input_option.flag == MPS_FLAG_RATIONAL)
            {
              if (MPS_STRUCTURE_IS_REAL (s->config))
                s->config->structure = MPS_STRUCTURE_REAL_RATIONAL;
              else if (MPS_STRUCTURE_IS_COMPLEX (s->config))
                s->config->structure = MPS_STRUCTURE_COMPLEX_RATIONAL;
            }
          else if (input_option.flag == MPS_FLAG_FP)
            {
              if (MPS_STRUCTURE_IS_FP (s->config))
                s->config->structure = MPS_STRUCTURE_REAL_FP;
              else if (MPS_STRUCTURE_IS_COMPLEX (s->config))
                s->config->structure = MPS_STRUCTURE_COMPLEX_FP;
            }
        }

      /* Free the line obtained */
      free (line);
    }

  /* Since the Degree is a required parameter, we ask that it is provided. */
  if (s->n == -1)
    mps_error (s, 1,
               "Degree of the polynomial must be provided via the Degree=%d configuration option.");
  else if (s->debug_level & MPS_DEBUG_IO)
    {
      MPS_DEBUG (s, "Degree: %d", s->n);
    }

  /* Check that the stream in the buffer is only composed
   * by spaces */
  /* for (i = 0; i < strlen (buffer->line); i++) */
  /*   { */
  /*     if (!isspace (buffer->line[i])) */
  /*       { */
  /*         mps_error (s, 1, */
  /*                    "Options and input data are not separated by a newline."); */
  /*         break; */
  /*       } */
  /*   } */

  if (MPS_REPRESENTATION_IS_SECULAR (s->config))
    {
      if (s->debug_level & MPS_DEBUG_IO)
        {
          MPS_DEBUG (s, "Parsing secular equation from stream");
        }
      mps_secular_equation_read_from_stream (s, buffer);
    }
  else if (MPS_REPRESENTATION_IS_MONOMIAL (s->config))
    {
      if (s->debug_level & MPS_DEBUG_IO)
	MPS_DEBUG (s, "Parsing polynomial from stream");

      mps_secular_equation_read_from_stream_poly (s, buffer);
    }

  mps_input_buffer_free (buffer);
}



/*********************************************************
*      SUBROUTINE READROOTS                              *
*********************************************************/
void
mps_readroots (mps_status * s)
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
    mpc_inp_str_u (s->mroot[i], s->rtstr, 10);
}

/*********************************************************
*      SUBROUTINE COUNTROOTS                             *
*********************************************************/
void
mps_countroots (mps_status * s)
{
  int k;

  if (s->DOLOG)
    fprintf (s->logstr, "Counting roots...\n");

  s->count[0] = s->count[1] = s->count[2] = 0;

  for (k = 0; k < s->n; k++)
    switch (s->status[k][2])
      {
      case 'i':
        s->count[0]++;
        break;
      case 'o':
        s->count[1]++;
        break;
      default:
        s->count[2]++;
        break;
      }

  if (s->goal[1] == 'o')
    s->count[1] += s->zero_roots;
  else
    s->count[0] += s->zero_roots;
}

/*********************************************************
*      SUBROUTINE OUTCOUNT                               *
*********************************************************/
void
mps_outcount (mps_status * s)
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
mps_outfloat (mps_status * s, mpf_t f, rdpe_t rad, long out_digit,
              mps_boolean sign)
{
  tmpf_t t;
  rdpe_t r, ro;
  double d;
  long l, digit, true_digit;

  if (s->goal[4] == 'f')
    {
      tmpf_init2 (t, mpf_get_prec (f));
      mpf_set (t, f);
      mpf_out_str (s->outstr, 10, 0, t);
      tmpf_clear (t);
      return;
    }

  tmpf_init2 (t, s->prec_out);

  mpf_get_rdpe (ro, f);
  if (s->goal[4] == 'g')
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

  tmpf_clear (t);
}

/*********************************************************
*      SUBROUTINE OUTROOT                                *
*********************************************************/
void
mps_outroot (mps_status * s, int i)
{
  static int num = 0;           /* output roots count */
  long out_digit;

  out_digit = (long) (LOG10_2 * s->prec_out) + 10;
  num++;

  /* print format header */
  switch (s->goal[4])
    {
    case 'c':
    case 'f':
      fprintf (s->outstr, "(");
      break;
    case 'v':
      fprintf (s->outstr, "Root(%d) = ", num);
      break;
    }

  /* print real part */
  if (i == ISZERO || s->status[i][1] == 'I')
    fprintf (s->outstr, "0");
  else
    mps_outfloat (s, mpc_Re (s->mroot[i]), s->drad[i], out_digit, true);

  /* print format middle part */
  switch (s->goal[4])
    {
    case 'b':
      fprintf (s->outstr, " ");
      break;
    case 'g':
      fprintf (s->outstr, "\t");
      break;
    case 'c':
    case 'f':
      fprintf (s->outstr, ", ");
      break;
    case 'v':
      if (i == ISZERO || mpf_sgn (mpc_Im (s->mroot[i])) >= 0)
        fprintf (s->outstr, " + I * ");
      else
        fprintf (s->outstr, " - I * ");
      break;
    }

  /* print imaginary part */
  if (i == ISZERO || s->status[i][1] == 'R')
    fprintf (s->outstr, "0");
  else
    mps_outfloat (s, mpc_Im (s->mroot[i]), s->drad[i], out_digit,
                  s->goal[4] != 'v');

  /* print format ending */
  switch (s->goal[4])
    {
    case 'c':
      fprintf (s->outstr, ")");
      break;
    case 'f':
      fprintf (s->outstr, ")\n");
      if (i != ISZERO)
        {
          rdpe_outln_str (s->outstr, s->drad[i]);
          fprintf (s->outstr, "%4.3s\n", s->status[i]);
        }
      else
        fprintf (s->outstr, " 0\n ---\n");
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
          fprintf (s->logstr, "root %-4d = ", i);
          mpc_out_str_2 (s->logstr, 10, 0, 0, s->mroot[i]);
          fprintf (s->logstr, "\n");
          fprintf (s->logstr, "  radius = ");
          rdpe_outln_str (s->logstr, s->drad[i]);
          fprintf (s->logstr, "  prec = %ld\n",
                   (long) (mpc_get_prec (s->mroot[i]) / LOG2_10));
          fprintf (s->logstr, "  status = %4.3s\n", s->status[i]);
          fprintf (s->logstr, "--------------------\n");
        }
    }
}

/*********************************************************
*      SUBROUTINE OUTPUT                                 *
*********************************************************/
void
mps_output (mps_status * s)
{
  int i, ind;

  if (s->DOLOG)
    fprintf (s->logstr, "--------------------\n");

  if (s->goal[0] == 'c')
    mps_outcount (s);
  else
    {
      if (s->goal[1] != 'o')
        for (i = 0; i < s->zero_roots; i++)
          mps_outroot (s, ISZERO);
      for (ind = 0; ind < s->n; ind++)
        {
          i = s->order[ind];
          if (s->status[i][2] == 'o')
            continue;
          mps_outroot (s, i);
        }
    }
}

/*********************************************************
*      SUBROUTINE COPY_ROOTS                             *
*********************************************************/
void
mps_copy_roots (mps_status * s)
{
  int i;

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
          mpc_set_prec (s->mroot[i], DBL_MANT_DIG);
          mpc_set_cplx (s->mroot[i], s->froot[i]);
          rdpe_set_d (s->drad[i], s->frad[i]);
        }
      break;

    case dpe_phase:
      if (s->DOSORT)
        mps_dsort (s);
      for (i = 0; i < s->n; i++)
        {
          mpc_set_prec (s->mroot[i], DBL_MANT_DIG);
          mpc_set_cdpe (s->mroot[i], s->droot[i]);
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
mps_dump (mps_status * s, FILE * dmpstr)
{
  int i;

  fprintf (dmpstr, "\nDumping...\n");

  /* output current status */
  fprintf (dmpstr,
           "Phase=%d, In=%d, Out=%d, Uncertain=%d, Zero=%d, Clusters=%d\n",
           s->lastphase, s->count[0], s->count[1], s->count[2], s->zero_roots,
           s->nclust);

  /* output current approximations */
  fprintf (dmpstr, "\nCurrent approximations:\n");
  for (i = 0; i < s->n; i++)
    {
      fprintf (dmpstr, "%d:\t", i);

      switch (s->lastphase)
        {
        case no_phase:
        case float_phase:
          cplx_outln_str (dmpstr, s->froot[i]);
          break;

        case dpe_phase:
          cdpe_outln_str (dmpstr, s->droot[i]);
          break;

        case mp_phase:
          mpc_outln_str (dmpstr, 10, 0, s->mroot[i]);
          break;
        }
    }

  /* output radii */
  fprintf (dmpstr, "\nCurrent radii:\n");
  for (i = 0; i < s->n; i++)
    {
      fprintf (dmpstr, "%d:\t", i);

      switch (s->lastphase)
        {
        case no_phase:
        case float_phase:
          fprintf (dmpstr, "%e\n", s->frad[i]);
          break;

        case dpe_phase:
        case mp_phase:
          rdpe_outln_str (dmpstr, s->drad[i]);
          break;
        }
    }

  /* output position */
  fprintf (dmpstr, "\nPos:\t");
  for (i = 0; i < s->n; i++)
    fprintf (dmpstr, "%4d", i);

  /* output status information */
  fprintf (dmpstr, "\nStatus:\t");
  for (i = 0; i < s->n; i++)
    fprintf (dmpstr, "%4.3s", s->status[i]);

  /* output cluster information */
  fprintf (dmpstr, "\nClust:\t");
  for (i = 0; i < s->n; i++)
    fprintf (dmpstr, "%4d", s->clust[i]);

  fprintf (dmpstr, "\nPunt:\t");
  for (i = 0; i < s->nclust; i++)
    fprintf (dmpstr, "%4d", s->punt[i]);

  fprintf (dmpstr, "\n\n");
}

/**
 * @brief Dump cluster structure to <code>outstr</code>.
 *
 * @param s the mps_status struct pointer.
 * @param outstr The output stream where the cluster structure
 *  will be dumped.
 */
void
mps_dump_cluster_structure (mps_status * s, FILE * outstr)
{
  int i, j;
  fprintf (outstr,
           "    MPS_DUMP_CLUSTER_STRUCTURE: Dumping cluster structure\n");

  for (i = 0; i < s->nclust; i++)
    {
      fprintf (outstr, "     Cluster %d contains %d roots:\n", i,
               s->punt[i + 1] - s->punt[i]);

      /* Dump cluster roots, but not more than 15 for line, to make
       * the output readable. */
      for (j = s->punt[i]; j < s->punt[i + 1]; j++)
        {
          /* Go to a newlint if 15 roots are printed out */
          if ((j - s->punt[i]) % 15 == 0)
            {
              fprintf (outstr, "\n       ");
            }

          printf (" %d", s->clust[j]);
        }

      /* Make space untile the next cluster */
      fprintf (outstr, "\n\n");
    }
}

/**
 * @brief Dump status of all the root approximations
 */
void
mps_dump_status (mps_status * s, FILE * outstr)
{
  int i;
  for (i = 0; i < s->n; i++)
    {
      fprintf (outstr, "s->status[%d] = ", i);
      fprintf (outstr, "'%c' '%c' '%c'\n", s->status[i][0],
               s->status[i][1], s->status[i][2]);
    }
}

/*************************************************************
 *                     SUBROUTINE WARN                       *
 *************************************************************/
void
mps_warn (mps_status * st, char *s)
{
  if (st->DOWARN)
    {
      if (s[strlen (s)] == '\n')
        {
          fprintf (st->logstr, "%s", s);
        }
      else
        {
          fprintf (st->logstr, "%s\n", s);
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
mps_error (mps_status * st, int args, ...)
{
  va_list ap;
  char *s;

  if (mps_is_a_tty (st->logstr))
    mps_warn (st, "\033[31;1m!\033[0m MPSolve encountered an error:");  /* output error message */
  else
    mps_warn (st, "! MPSolve encountered an error:");
  va_start (ap, args);
  while (args--)
    {
      s = va_arg (ap, char *);
      mps_warn (st, s);         /* output error message */
    }
  va_end (ap);

  /* Dump approximations, but only if they are present */
  if (st->froot && st->lastphase)
    mps_dump (st, st->logstr);  /* dump status          */
  exit (EXIT_FAILURE);          /* exit program         */
}
