#include "mex.h"
#include "mps_option_parser.h"
#include <mps/mps.h>
#include <stdio.h>

void
mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
  const mwSize * dims;
  const mxArray * cellArray;
  mps_context * ctx;
  mps_monomial_poly * mp;
  cplx_t * results = NULL;
  int i, n;
  double *real_res, *imag_res;
  mxArray *roots;

  _mps_matlab_options options = mps_parse_matlab_options ( (nrhs > 1) ? prhs[1] : NULL );

  dims = mxGetDimensions (prhs[0]);
  n = dims[1];

  ctx = mps_context_new ();
  mp = mps_monomial_poly_new (ctx, n - 1);
  
  for (i = 0; i < n; i++)
    {
      int size; 
      char * chars = NULL;

      cellArray = mxGetCell (prhs[0], i);
      size = mxGetN (cellArray);

      chars = (char *) malloc (sizeof (char) * (size + 1));
      mxGetString (cellArray, chars, size + 1);

      mps_monomial_poly_set_coefficient_s (ctx, mp, i, chars, NULL);   
      free (chars);
    }

  mps_context_set_input_poly (ctx, MPS_POLYNOMIAL (mp));

  /* Uncomment this to enable debug on the MPSolve operations */
  /* mps_context_add_debug_domain (ctx, MPS_DEBUG_TRACE); */

  mps_context_select_algorithm (ctx, options.algorithm);
  mps_context_set_output_prec (ctx, options.digits / log10(2));
  mps_context_set_output_goal (ctx, options.goal);

  mps_mpsolve (ctx);

  /* Only output the roots as floating point numbers if the
   * requested precision is less than 16 digits. */
  if (options.digits <= DBL_DIG)
    {
      mps_context_get_roots_d (ctx, &results, NULL);

      roots = mxCreateDoubleMatrix (n - 1, 1, mxCOMPLEX);
      real_res = mxGetPr (roots);
      imag_res = mxGetPi (roots);

      for (i = 0; i < n - 1; i++)
	{
	  real_res[i] = cplx_Re (results[i]);
	  imag_res[i] = cplx_Im (results[i]);
	}

      /* Free used data */
      cplx_vfree (results);
    }
  else
    {
      mpc_t * mresults = NULL;
      char * buffer_r, * buffer_i;
      int ndim = 2, dims[] = { n-1, 4 };

      roots = mxCreateCellArray(ndim, dims);

      mps_context_get_roots_m (ctx, &mresults, NULL);
      
      for (i = 0; i < n - 1; i++)
	{
	  mp_exp_t rexp, iexp;
          int indices[] = { i, 0 };

	  buffer_r = buffer_i = NULL;
	  buffer_r = mpf_get_str (buffer_r, &rexp, 10, 0, mpc_Re (mresults[i]));
	  buffer_i = mpf_get_str (buffer_i, &iexp, 10, 0, mpc_Im (mresults[i]));

	  mpc_clear (mresults[i]);

          indices[1] = 0;
	  mxSetCell (roots,
		     mxCalcSingleSubscript(roots, 2, indices),
		     mxCreateString(buffer_r));

	  mxArray * rexpM = mxCreateDoubleMatrix (1, 1, mxREAL);
	  mxGetPr(rexpM)[0] = rexp;
          indices[1] = 1;
	  mxSetCell (roots,
		     mxCalcSingleSubscript(roots, 2, indices),
		     rexpM);

          indices[1] = 2;
	  mxSetCell (roots,
		     mxCalcSingleSubscript(roots, 2, indices),
		     mxCreateString(buffer_i));

	  mxArray * iexpM = mxCreateDoubleMatrix (1, 1, mxREAL);
	  mxGetPr(iexpM)[0] = iexp;

          indices[1] = 3;
	  mxSetCell (roots,
		     mxCalcSingleSubscript(roots, 2, indices),
		     iexpM);

	  free (buffer_r);
	  free (buffer_i);
	}

      free (mresults);
    }

  /* Return the roots */
  plhs[0] = roots;
      
  mps_context_free (ctx);
}
