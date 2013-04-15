/*
 * This file is part of MPSolve 3.0
 *
 * Copyright (C) 2001-2013, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors: 
 *   Leonardo Robol <robol@mail.dm.unipi.it>
 */

#include <mps/mps.h>

/**
 * @brief Parse the stream that has been loaded into buffer and that
 * describe a mps_secular_equation. 
 *
 * @param s The current mps_context
 * @param buffer The buffer that needs to be parsed
 * @param The structure of the secular equation 
 * @param The density configuration of the secular equation
 *
 * @return A newly allocated mps_secular_equation, or NULL if the parsing fails.
 */
mps_secular_equation *
mps_secular_equation_read_from_stream (mps_context * s,
                                       mps_input_buffer * buffer,
                                       mps_structure structure,
                                       mps_density density)
{
  mps_secular_equation *sec;
  int i;
  mpf_t ftmp;
  char * token;

  mpf_init2 (ftmp, 64);

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
              free (token);
              
              mps_polynomial_free (s, MPS_POLYNOMIAL (sec));
              sec = NULL;
              goto cleanup;
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
                  free (token);
              
                  mps_polynomial_free (s, MPS_POLYNOMIAL (sec));
                  sec = NULL;
                  goto cleanup;
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
              free (token);
              
              mps_polynomial_free (s, MPS_POLYNOMIAL (sec));
              sec = NULL;
              goto cleanup;
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
                  free (token);

                  mps_polynomial_free (s, MPS_POLYNOMIAL (sec));
                  sec = NULL;
                  goto cleanup;
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
                                       "Error reading some coefficients of the secular equation.\n"
                                       "Please check your input file");
              
              /* Cleanup temporary variables and exit */
              free (token);
              mps_polynomial_free (s, MPS_POLYNOMIAL (sec));
              sec = NULL;

              goto cleanup;
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
                                           "Error reading some coefficients of the secular equation."
                                           "Please check your input file");
                  free (token);
              
                  mps_polynomial_free (s, MPS_POLYNOMIAL (sec));
                  sec = NULL;
                  goto cleanup;
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
                                       "Error reading some coefficients of the secular equation."
                                       "Please check your input file");
              free (token);
              
              mps_polynomial_free (s, MPS_POLYNOMIAL (sec));
              sec = NULL;
              goto cleanup;
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
                                           "Error reading some coefficients of the secular equation."
                                           "Please check your input file");
                  free (token);
              
                  mps_polynomial_free (s, MPS_POLYNOMIAL (sec));
                  sec = NULL;
                  goto cleanup;
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

cleanup:

  mpf_clear (ftmp);

  return sec;
}
