#include "mex.h"
#include <mps/mps.h>
#include <stdlib.h>
#include <string.h>


void 
mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{

  int i, m, n;
  mps_context *s = NULL;
  mps_monomial_poly *mp = NULL;
  double *real_coeff, *imag_coeff, *real_res, *imag_res;
  mxArray *roots;
  cplx_t *coefficients, *results;

  /* Get dimension of input data, but first check consistency */
  if (nrhs < 1 || nrhs > 2)
    {
      mexErrMsgTxt ("This function takes two parameters, a vector and an optional switch.");
    }

  if (nlhs > 1)
    {
      mexErrMsgTxt ("This function can return only one value.");
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
  mp = mps_monomial_poly_new (s, n - 1);

  /* Set coefficients */
  real_coeff = mxGetPr (prhs[0]);
  imag_coeff = mxGetPi (prhs[0]);
  for (i = 0; i < n; i++)
    {
        mps_monomial_poly_set_coefficient_d (s, mp, n - i - 1, real_coeff[i],
                                            (imag_coeff) ? imag_coeff[i] : 0.0);
    }

  /* Solve the equation */
  mps_context_set_input_poly (s, MPS_POLYNOMIAL (mp));
  mps_context_set_output_goal (s, MPS_OUTPUT_GOAL_APPROXIMATE);

  /* Check if the second parameter was passed */
  if (nrhs == 2)
    {
      if (!mxIsChar(prhs[1]))
        mexErrMsgTxt("The second parameter must be of type string");
      char alg[255];
      mxGetString(prhs[1], alg, 254);
      if (strlen(alg) > 1)
        mexErrMsgTxt("The second parameter must be a single character");

      switch (alg[0])
        {
        case 's':
          mps_context_select_algorithm (s, MPS_ALGORITHM_SECULAR_GA);
          break;
        case 'u':
          mps_context_select_algorithm (s, MPS_ALGORITHM_STANDARD_MPSOLVE);
          break;
        default:
          mexErrMsgTxt("The selected algorithm is not supported, use 's' or 'u'");
          break;
        }
    }

  mps_mpsolve (s);

  /* Get results back */
  mps_context_get_roots_d (s, &results, NULL);

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
  mps_context_free (s);

  /* Return the roots */
  plhs[0] = roots;
}
