#include "mex.h"
#include "mps_option_parser.h"
#include <mps/mps.h>
#include <stdlib.h>
#include <string.h>

void 
mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{

  int i, m, n;
  mps_context *s = NULL;
  mps_polynomial *poly = NULL;
  double *real_coeff, *imag_coeff, *real_res, *imag_res;
  double *radius = NULL, *out_radius;
  mxArray *roots, *radius_vec;
  cplx_t *coefficients, *results;

  _mps_matlab_options options = mps_parse_matlab_options ( (nrhs > 1) ? prhs[1] : NULL );

  /* Get dimension of input data, but first check consistency */
  if (nrhs < 1 || nrhs > 2)
    {
      mexErrMsgTxt ("This function takes two parameters, a vector and an optional switch.");
    }

  if (nlhs > 1 && !options.radius)
    {
      mexErrMsgTxt ("This function can return only one value, unless the option \"radius\" is specified");
    }

  if (!mxIsNumeric (prhs[0]))
    {
      mexErrMsgTxt ("The first parameter must be a numeric array.");
    }

  /* Get dimensions */
  m = mxGetM (prhs[0]);
  n = mxGetN (prhs[0]);

  /* Check that this is a single dimension array */
  if (m != 1 && n != 1)
    {
      mexErrMsgTxt ("A 1D array must be provided as argument.");
    }

  /* Get the maximum between m and n */
  n = (m < n) ? n : m;

  /* Inizialize mpsolve */
  s = mps_context_new ();

  /* Allocate coefficients and space for the results */
  results = cplx_valloc (n - 1);
  if (options.chebyshev) 
    {
      poly = MPS_POLYNOMIAL(mps_chebyshev_poly_new (s, n - 1, MPS_STRUCTURE_COMPLEX_FP));
    }
  else
    {
      poly = MPS_POLYNOMIAL(mps_monomial_poly_new (s, n - 1));
    }

  /* Set coefficients */
  real_coeff = mxGetPr (prhs[0]);
  imag_coeff = mxGetPi (prhs[0]);
  mpc_t coeff;
  mpc_init2(coeff, 64);
  for (i = 0; i < n; i++)
    {
        if (options.chebyshev) 
          {
             mpc_set_d(coeff, real_coeff[i], (imag_coeff) ? imag_coeff[i] : 0.0);
             mps_chebyshev_poly_set_coefficient_f (s, MPS_CHEBYSHEV_POLY(poly), n - i - 1, coeff);
          }
        else
          {
            mps_monomial_poly_set_coefficient_d (s, MPS_MONOMIAL_POLY(poly), n - i - 1, real_coeff[i],
                                                (imag_coeff) ? imag_coeff[i] : 0.0);
          }
    }
  mpc_clear(coeff);

  /* Solve the equation */
  mps_context_set_input_poly (s, poly);
  mps_context_set_output_goal (s, MPS_OUTPUT_GOAL_APPROXIMATE);

  mps_context_select_algorithm (s, options.algorithm);
  mps_context_set_output_prec (s, options.digits / log10 (2));

  mps_mpsolve (s);

  /* Get results back */
  mps_context_get_roots_d (s, &results, &radius);

  roots = mxCreateDoubleMatrix (n - 1, 1, mxCOMPLEX);
  real_res = mxGetPr (roots);
  imag_res = mxGetPi (roots);

  options.radius = (nlhs > 1) && options.radius;

  if (options.radius)
    {
      radius_vec = mxCreateDoubleMatrix (n - 1, 1, mxREAL);
      out_radius = mxGetPr (radius_vec);
    }

  for (i = 0; i < n - 1; i++)
    {
      real_res[i] = cplx_Re (results[i]);
      imag_res[i] = cplx_Im (results[i]);

      if (options.radius)
	out_radius[i] = radius[i];
    }

  /* Free used data */
  cplx_vfree (results);
  free (radius);

  mps_polynomial_free (s, poly);
  mps_context_free (s);

  /* Return the roots */
  plhs[0] = roots;

  if (options.radius)
    plhs[1] = radius_vec;
}
