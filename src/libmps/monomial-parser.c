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

#include <mps/mps.h>

/**
 * @brief Parse the stream that has been loaded into buffer and that
 * describe a mps_monomial_poly. 
 *
 * @param s The current mps_context
 * @param buffer The buffer that needs to be parsed
 * @param The structure of the polynomial 
 * @param The density configuration of the polynomial.
 *
 * @return A newly allocated mps_polynomial, or NULL if the parsing fails.
 */
mps_monomial_poly *
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
                {
                  mps_raise_parsing_error (s, buffer, token, "Error parsing coefficients of the polynomial");
                  free (token); mps_polynomial_free (s, MPS_POLYNOMIAL (poly));
                  return NULL;
                }
              free (token);

              if (MPS_STRUCTURE_IS_COMPLEX (structure))
                {
                  token = mps_input_buffer_next_token (buffer);
                  if (!token || (mpf_set_str (mpc_Im (poly->mfpc[i]), token, 10) != 0))
                    {
                      mps_raise_parsing_error (s, buffer, token, "Error parsing coefficients of the polynomial");
                      free (token); mps_polynomial_free (s, MPS_POLYNOMIAL (poly));
                      return NULL;
                    }
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
                {
                  mps_raise_parsing_error (s, buffer, token, "Error parsing coefficients of the polynomial");
                  free (token); mps_polynomial_free (s, MPS_POLYNOMIAL (poly));
                  return NULL;
                }
              mpq_canonicalize (poly->initial_mqp_r[i]);
              free (token);
      
              if (MPS_STRUCTURE_IS_COMPLEX (structure))
                {
                  token = mps_input_buffer_next_token (buffer);
                  if (!token || (mpq_set_str (poly->initial_mqp_i[i], token, 10) != 0))
                    {
                      mps_raise_parsing_error (s, buffer, token, "Error parsing coefficients of the polynomial");
                      free (token); mps_polynomial_free (s, MPS_POLYNOMIAL (poly));
                      return NULL;
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
          {
            mps_raise_parsing_error (s, buffer, token, "Error while parsing the degree of a monomial");
            free (token); mps_polynomial_free (s, MPS_POLYNOMIAL (poly));
            return NULL;
          }

          if (i < 0 || i > s->n) 
            {
              mps_raise_parsing_error (s, buffer, token, "Degree of coefficient out of bounds");
              free (token); mps_polynomial_free (s, MPS_POLYNOMIAL (poly));
              return NULL;
            }

          if (poly->spar[i]) 
            {
              mps_raise_parsing_error (s, buffer, token, "A monomial of the same degree has been inserted twice"); 
              free (token); mps_polynomial_free (s, MPS_POLYNOMIAL (poly));
              return NULL;
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
                  free (token); mps_polynomial_free (s, MPS_POLYNOMIAL (poly));
                  return NULL;
                }
              free (token);
          
              if (MPS_STRUCTURE_IS_COMPLEX (structure))
                {
                  token = mps_input_buffer_next_token (buffer);
                  if (!token || (mpf_set_str (mpc_Im (poly->mfpc[i]), token, 10) != 0))
                    {
                      mps_raise_parsing_error (s, buffer, token, "Error parsing coefficients of the polynomial");
                      free (token); mps_polynomial_free (s, MPS_POLYNOMIAL (poly));
                      return NULL;
                    }
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
                {
                  mps_raise_parsing_error (s, buffer, token, "Error parsing coefficients of the polynomial");
                  free (token); mps_polynomial_free (s, MPS_POLYNOMIAL (poly));
                  return NULL;
                }
              mpq_canonicalize (poly->initial_mqp_r[i]);
              free (token);
      
              if (MPS_STRUCTURE_IS_COMPLEX (structure))
                {
                  token = mps_input_buffer_next_token (buffer);
                  if (!token || (mpq_set_str (poly->initial_mqp_i[i], token, 10) != 0))
                    {
                      mps_raise_parsing_error (s, buffer, token, "Error parsing coefficients of the polynomial");
                      free (token); mps_polynomial_free (s, MPS_POLYNOMIAL (poly));
                      return NULL;
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
        }
      else
        {
          cplx_set (poly->fpc[i], cplx_zero);
          cdpe_set (poly->dpc[i], cdpe_zero);
          
          rdpe_set (poly->dap[i], rdpe_zero);
          poly->fap[i] = 0.0f;
        }
    }

  mpf_clear (ftmp);
  return poly;
}

/**
 * @brief Parse the stream that has been loaded into buffer and that
 * describe a mps_monomial_poly. This function parse polynomials described
 * in the format of MPSolve 2.2. 
 *
 * @param s The current mps_context
 * @param buffer The buffer that needs to be parsed
 * @param The structure of the polynomial 
 * @param The density configuration of the polynomial.
 *
 * @return A newly allocated mps_polynomial, or NULL if the parsing fails.
 */
mps_polynomial *
mps_monomial_poly_read_from_stream_v2 (mps_context * s, mps_input_buffer * buffer)
{
  int i;
  mps_monomial_poly *poly = NULL;
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
      mps_error (s, "Error parsing the input file");

      goto cleanup;
    }

  if (token)
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
      mps_error (s, "Found unsupported data_type in input file");
      goto cleanup;
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
      mps_error (s, "Found unsupported data_structure in input file");
      goto cleanup;
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
      mps_error (s, "Found unsupported data structure in input file");
      goto cleanup;
      break;
    }

  /* Read precision and degree */
  prec = 0;
  token = mps_input_buffer_next_token (buffer);
  if (!token || !sscanf (token, "%ld", &prec))
    {
      mps_error (s, "Error while reading the input precision of the coefficients");

      if (token)
        free (token);

      goto cleanup;
    }
  else 
    prec *= LOG2_10;
  free (token);

  token = mps_input_buffer_next_token (buffer);
  if (!token || !sscanf (token, "%d", &s->n))
    {
      mps_error (s, "Error reading the degree of the polynomial");

      if (token)
        free (token);

      goto cleanup;
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

      return MPS_POLYNOMIAL (user_poly);
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

                  /* Cleanup temporary variables and exit */
                  if (token)
                    free (token);
                  mps_polynomial_free (s, MPS_POLYNOMIAL (poly));
                  poly = NULL;

                  goto cleanup;
                }
              free (token);

              if (MPS_STRUCTURE_IS_COMPLEX (structure))
                {
                  token = mps_input_buffer_next_token (buffer);
                  if (!token || (mpf_set_str (mpc_Im (poly->mfpc[i]), token, 10) != 0))
                    {
                      mps_raise_parsing_error (s, buffer, token, "Error parsing coefficients of the polynomial");

                      /* Cleanup temporary variables and exit */
                      if (token)
                        free (token);
                      mps_polynomial_free (s, MPS_POLYNOMIAL (poly));
                      poly = NULL;

                      goto cleanup;                    }
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
                  
                  /* Cleanup temporary variables and exit */
                  if (token)
                    free (token);
                  mps_polynomial_free (s, MPS_POLYNOMIAL (poly));
                  poly = NULL;

                  goto cleanup;
                }
              mpq_canonicalize (poly->initial_mqp_r[i]);
              free (token);
      
              if (MPS_STRUCTURE_IS_COMPLEX (structure))
                {
                  token = mps_input_buffer_next_token (buffer);
                  if (!token || (mpq_set_str (poly->initial_mqp_i[i], token, 10) != 0))
                    {
                      mps_raise_parsing_error (s, buffer, token, "Error parsing coefficients of the polynomial");
                      
                      /* Cleanup temporary variables and exit */
                      if (token)
                        free (token);
                      mps_polynomial_free (s, MPS_POLYNOMIAL (poly));
                      poly = NULL;

                      goto cleanup;
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

                  /* Cleanup temporary variables and exit */
                  if (token)
                    free (token);
                  mps_polynomial_free (s, MPS_POLYNOMIAL (poly));
                  poly = NULL;

                  goto cleanup;
                }
              mpq_set (poly->initial_mqp_r[i], qtmp);
              free (token);

              /* Denominator of the real part of the coefficient */
              token = mps_input_buffer_next_token (buffer);
              if (!token || (mpq_set_str (qtmp, token, 10)) != 0)
                {
                  mps_raise_parsing_error (s, buffer, token, "Error parsing the denominator of a coefficient");

                  /* Cleanup temporary variables and exit */
                  if (token)
                    free (token);
                  mps_polynomial_free (s, MPS_POLYNOMIAL (poly));
                  poly = NULL;

                  goto cleanup;
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

                      /* Cleanup temporary variables and exit */
                      if (token)
                        free (token);
                      mps_polynomial_free (s, MPS_POLYNOMIAL (poly));
                      poly = NULL;

                      goto cleanup;
                    }
                  mpq_set (poly->initial_mqp_i[i], qtmp);
                  free (token);

                  /* Denominator of the real part of the coefficient */
                  token = mps_input_buffer_next_token (buffer);
                  if (!token || (mpq_set_str (qtmp, token, 10)) != 0)
                    {
                      mps_raise_parsing_error (s, buffer, token, "Error parsing the denominator of a coefficient");

                      /* Cleanup temporary variables and exit */
                      if (token)
                        free (token);
                      mps_polynomial_free (s, MPS_POLYNOMIAL (poly));
                      poly = NULL;

                      goto cleanup;
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
          if (!token || !sscanf (token, "%d", &i))
            {
              mps_raise_parsing_error (s, buffer, token, "Error while parsing the degree of a monomial");

              /* Cleanup temporary variables and exit */
              if (token)
                free (token);
              mps_polynomial_free (s, MPS_POLYNOMIAL (poly));
              poly = NULL;

              goto cleanup;
            }

          if (i < 0 || i > s->n) 
            {
              mps_raise_parsing_error (s, buffer, token, "Degree of coefficient out of bounds");
              
              /* Cleanup temporary variables and exit */
              if (token)
                free (token);
              mps_polynomial_free (s, MPS_POLYNOMIAL (poly));
              poly = NULL;

              goto cleanup;
            }

          if (poly->spar[i])
            {
              mps_raise_parsing_error (s, buffer, token, "A monomial of the same degree has been inserted twice");
              
              /* Cleanup temporary variables and exit */
              if (token)
                free (token);
              mps_polynomial_free (s, MPS_POLYNOMIAL (poly));
              poly = NULL;

              goto cleanup;
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

                  /* Cleanup temporary variables and exit */
                  if (token)
                    free (token);
                  mps_polynomial_free (s, MPS_POLYNOMIAL (poly));
                  poly = NULL;

                  goto cleanup;
                }
              free (token);
          
              if (MPS_STRUCTURE_IS_COMPLEX (structure))
                {
                  token = mps_input_buffer_next_token (buffer);
                  if (!token || (mpf_set_str (mpc_Im (poly->mfpc[i]), token, 10) != 0))
                    {
                      mps_raise_parsing_error (s, buffer, token, "Error parsing coefficients of the polynomial");

                      /* Cleanup temporary variables and exit */
                      if (token)
                        free (token);
                      mps_polynomial_free (s, MPS_POLYNOMIAL (poly));
                      poly = NULL;

                      goto cleanup;
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
                  
                  /* Cleanup temporary variables and exit */
                  if (token)
                    free (token);
                  mps_polynomial_free (s, MPS_POLYNOMIAL (poly));
                  poly = NULL;

                  goto cleanup;
                }
              mpq_canonicalize (poly->initial_mqp_r[i]);
              free (token);
      
              if (MPS_STRUCTURE_IS_COMPLEX (structure))
                {
                  token = mps_input_buffer_next_token (buffer);
                  if (!token || (mpq_set_str (poly->initial_mqp_i[i], token, 10) != 0))
                    {
                      mps_raise_parsing_error (s, buffer, token, "Error parsing coefficients of the polynomial");

                      /* Cleanup temporary variables and exit */
                      if (token)
                        free (token);
                      mps_polynomial_free (s, MPS_POLYNOMIAL (poly));
                      poly = NULL;

                      goto cleanup;
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

                  /* Cleanup temporary variables and exit */
                  if (token)
                    free (token);
                  mps_polynomial_free (s, MPS_POLYNOMIAL (poly));
                  poly = NULL;

                  goto cleanup;
                }
              mpq_set (poly->initial_mqp_r[i], qtmp);
              free (token);

              /* Denominator of the real part of the coefficient */
              token = mps_input_buffer_next_token (buffer);
              if (!token || (mpq_set_str (qtmp, token, 10)) != 0)
                {
                  mps_raise_parsing_error (s, buffer, token, "Error parsing the denominator of a coefficient");
                  
                  /* Cleanup temporary variables and exit */
                  if (token)
                    free (token);
                  mps_polynomial_free (s, MPS_POLYNOMIAL (poly));
                  poly = NULL;

                  goto cleanup;
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
                          
                      /* Cleanup temporary variables and exit */
                      if (token)
                        free (token);
                      mps_polynomial_free (s, MPS_POLYNOMIAL (poly));
                      poly = NULL;

                      goto cleanup;
                    }
                  mpq_set (poly->initial_mqp_i[i], qtmp);
                  free (token);

                  /* Denominator of the real part of the coefficient */
                  token = mps_input_buffer_next_token (buffer);
                  if (!token || (mpq_set_str (qtmp, token, 10)) != 0)
                    {
                      mps_raise_parsing_error (s, buffer, token, "Error parsing the denominator of a coefficient");
                          
                      /* Cleanup temporary variables and exit */
                      if (token)
                        free (token);
                      mps_polynomial_free (s, MPS_POLYNOMIAL (poly));
                      poly = NULL;

                      goto cleanup;
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
              MPS_DEBUG_MPC (s, 15, poly->mfpc[i], "Coefficient of degree %d", i);
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
    }

  MPS_POLYNOMIAL (poly)->prec = prec;

cleanup:  

  mpf_clear (ftmp);
  mpq_clear (qtmp);

  return MPS_POLYNOMIAL (poly);
}

