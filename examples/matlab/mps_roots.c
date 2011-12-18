#include "mex.h"
#include <mps/mps.h>
#include <stdlib.h>


void
mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{

  int i, m, n;
  mps_status *s;
  double *real_coeff, *imag_coeff, *real_res, *imag_res;
  mxArray *roots;
  cplx_t *coefficients, *results;

  /* Get dimension of input data, but first check consistency */
  if (nrhs != 1)
    {
      mexErrMsgTxt ("The first parameter must be a vector.");
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
  s = mps_status_new ();

  /* Allocate coefficients and space for the results */
  coefficients = cplx_valloc (n);
  results = cplx_valloc (n - 1);

  /* Set coefficients */
  real_coeff = mxGetPr (prhs[0]);
  imag_coeff = mxGetPi (prhs[0]);
  for (i = 0; i < n; i++)
    {
      if (!imag_coeff)
        cplx_set_d (coefficients[n - i - 1], real_coeff[i], 0.0);
      else
        cplx_set_d (coefficients[n - i - 1], real_coeff[i], imag_coeff[i]);
    }

  /* Solve the equation */
  mps_status_set_poly_d (s, coefficients, n - 1);
  mps_mpsolve (s);

  /* Get results back */
  mps_status_get_roots_d (s, results, NULL);

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
  cplx_vfree (coefficients);
  mps_status_free (s);

  /* Return the roots */
  plhs[0] = roots;
}
