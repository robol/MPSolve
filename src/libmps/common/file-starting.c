/*
 * This file is part of MPSolve 3.1.9
 *
 * Copyright (C) 2001-2020, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Leonardo Robol <leonardo.robol@unipi.it>
 */

#include <mps/mps.h>

/**
 * @brief Select the starting points for the polynomial by loading the approximations loaded
 * in the file that has been opened in the given FILE* pointer. 
 *
 * @param ctx The current mps_context. 
 * @param n The number of approximations to load.
 * @param input The result of fopen() on the desired file. 
 * @param approximations The approximations that will be set according to the data found in
 * the input file. 
 */
static void
mps_load_approximations (mps_context * ctx, int n, FILE * input, 
			 mps_approximation ** approximations)
{
  MPS_DEBUG_THIS_CALL (ctx);
  int i;

  for (i = 0; i < n; i++)
    {
      mps_approximation * appr = approximations[i];      

      if (mpc_inp_str (appr->mvalue, input, 10) == 0)
	{
	  MPS_DEBUG (ctx, "Error while trying to read the %d-th approximation. Aborting", i);
	  mps_error (ctx, "Error while trying to read the %d-th approximation. Aborting", i);
	  return;
	}
      else
	{
	  char next; 

	  mpc_get_cplx (appr->fvalue, appr->mvalue);
	  mpc_get_cdpe (appr->dvalue, appr->mvalue);

	  /* This is essentialy a workaround to handle newlines. If we leave them in 
	   * mpc_inp_str will have no clue of what they are and will fail. */
	  if ((next = fgetc (input)) != '\n')
	    ungetc (next, input);
	}
    }
}

/**
 * @brief Select the starting points for the polynomial by loading the approximations loaded
 * in the file that has been opened in ctx->rtstr. 
 *
 * @param ctx The current mps_context. 
 * @param poly The polynomial for which the approxmimations should be selected. 
 * @param approximations The approximations that will be set according to the data found in
 * the input file. 
 */
void 
mps_file_fstart (mps_context * ctx, mps_polynomial * poly, mps_approximation ** approximations)
{
  mps_load_approximations (ctx, poly->degree, ctx->rtstr, approximations);
}

/**
 * @brief Select the starting points for the polynomial by loading the approximations loaded
 * in the file that has been opened in ctx->rtstr. 
 *
 * @param ctx The current mps_context. 
 * @param poly The polynomial for which the approxmimations should be selected. 
 * @param approximations The approximations that will be set according to the data found in
 * the input file. 
 */
void
mps_file_dstart (mps_context * ctx, mps_polynomial * poly, mps_approximation ** approximations)
{
  mps_load_approximations (ctx, poly->degree, ctx->rtstr, approximations);
}

/**
 * @brief Select the starting points for the polynomial by loading the approximations loaded
 * in the file that has been opened in ctx->rtstr. 
 *
 * @param ctx The current mps_context. 
 * @param poly The polynomial for which the approxmimations should be selected. 
 * @param approximations The approximations that will be set according to the data found in
 * the input file. 
 */
void
mps_file_mstart (mps_context * ctx, mps_polynomial * poly, mps_approximation ** approximations)
{
  mps_load_approximations (ctx, poly->degree, ctx->rtstr, approximations);
}

