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
 * describe a mps_chebyshev_poly.
 *
 * @param ctx The current mps_context
 * @param buffer The buffer that needs to be parsed
 * @param structure The structure of the polynomial 
 * @param density The density configuration of the polynomial.
 *
 * @return A newly allocated mps_chebyshev_poly, or NULL if the parsing fails.
 */
mps_chebyshev_poly *
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
                  {
                    mps_raise_parsing_error (ctx, buffer, token, 
                      "Error while reading real part of coefficient");
                    free (token); mps_polynomial_free (ctx, MPS_POLYNOMIAL (cpoly));
                    return NULL;
                  }
                free (token);

                if (MPS_STRUCTURE_IS_COMPLEX (structure))
                  {
                    token = mps_input_buffer_next_token (buffer);

                    if (!token || (mpf_set_str (mpc_Im (cpoly->mfpc[i]), token, 10) != 0)) 
                      {
                        mps_raise_parsing_error (ctx, buffer, token, "Error while reading imaginary part of coefficient");
                        free (token); mps_polynomial_free (ctx, MPS_POLYNOMIAL (cpoly));
                        return NULL;
                      }
                    else
                      free (token);
                  }
                }
              else if (MPS_STRUCTURE_IS_RATIONAL (structure) || MPS_STRUCTURE_IS_INTEGER (structure))
                {
                  token = mps_input_buffer_next_token (buffer);

                  if (!token || (mpq_set_str (cpoly->rational_real_coeffs[i], token, 10)))
                    {
                      mps_raise_parsing_error (ctx, buffer, token, 
                        "Error while reading the real part of coefficient");
                      free (token); mps_polynomial_free (ctx, MPS_POLYNOMIAL (cpoly));
                      return NULL;
                    }
                  free (token);

                  if (MPS_STRUCTURE_IS_COMPLEX (structure))
                    {
                      token = mps_input_buffer_next_token (buffer);

                      if (!token || (mpq_set_str (cpoly->rational_imag_coeffs[i], token, 10)))
                        {
                          mps_raise_parsing_error (ctx, buffer, token, 
                            "Error while reading the imaginary part of coefficient");
                          free (token); mps_polynomial_free (ctx, MPS_POLYNOMIAL (cpoly));
                          return NULL;
                        }
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
          {
            mps_raise_parsing_error (ctx, buffer, token, "Cannot parse the degree of the coefficient.");
            free (token); mps_polynomial_free (ctx, MPS_POLYNOMIAL (cpoly));
            return NULL;
          }

        free (token);

        if (MPS_STRUCTURE_IS_FP (structure))
          {
            token = mps_input_buffer_next_token (buffer);

            if (!token || (mpf_set_str (mpc_Re (cpoly->mfpc[degree]), token, 10) != 0)) 
              {
                mps_raise_parsing_error (ctx, buffer, token, "Error while reading real part of coefficient");
                free (token); mps_polynomial_free (ctx, MPS_POLYNOMIAL (cpoly));
                return NULL;
              }

            free (token);

            if (MPS_STRUCTURE_IS_COMPLEX (structure))
              {
                token = mps_input_buffer_next_token (buffer);

                if (!token || (mpf_set_str (mpc_Im (cpoly->mfpc[degree]), token, 10) != 0)) 
                  {
                    mps_raise_parsing_error (ctx, buffer, token, 
                      "Error while reading imaginary part of coefficient %d", degree);
                    free (token); mps_polynomial_free (ctx, MPS_POLYNOMIAL (cpoly));
                    return NULL;
                  }

                free (token);                
              }
          }
        else
          {
            token = mps_input_buffer_next_token (buffer);

            if (!token || (mpq_set_str (cpoly->rational_real_coeffs[degree], token, 10)))
              {
                mps_raise_parsing_error (ctx, buffer, token, "Error while reading the real part of coefficient %d", i);
                free (token); mps_polynomial_free (ctx, MPS_POLYNOMIAL (cpoly));
                return NULL;
              }

            free (token);

            if (MPS_STRUCTURE_IS_COMPLEX (structure))
              {
                token = mps_input_buffer_next_token (buffer);

                if (!token || (mpq_set_str (cpoly->rational_imag_coeffs[degree], token, 10)))
                  {
                    mps_raise_parsing_error (ctx, buffer, token, 
                      "Error while reading the imaginary part of coefficient %d", i);
                    free (token); mps_polynomial_free (ctx, MPS_POLYNOMIAL (cpoly));
                    return NULL;
                  }

                free (token);
              }
          }

        mpc_get_cdpe (cpoly->dpc[degree], cpoly->mfpc[degree]);
        mpc_get_cplx (cpoly->fpc[degree], cpoly->mfpc[degree]);
      break;

      default:
        mps_error (ctx, "Only MPS_DENSITY_DENSE and MPS_DENSITY_SPARSE are supported in Chebyshev polynomials.");
        mps_polynomial_free (ctx, MPS_POLYNOMIAL (cpoly));
        return NULL;
        break;
    }

  return cpoly;
}

