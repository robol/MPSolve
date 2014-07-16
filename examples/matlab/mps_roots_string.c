#include "mex.h"
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

      chars = malloc (sizeof (char) * (size + 1));
      mxGetString (cellArray, chars, size + 1);

      mps_monomial_poly_set_coefficient_s (ctx, mp, i, chars, NULL);   
      free (chars);
    }

  mps_context_set_input_poly (ctx, MPS_POLYNOMIAL (mp));

  /* Uncomment this to enable debug on the MPSolve operations */
  /* mps_context_add_debug_domain (ctx, MPS_DEBUG_TRACE); */

  mps_mpsolve (ctx);

  /* Output the roots as string to allow not to loose 
   * digits. Still to be implemented. */
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
  mps_context_free (ctx);

  /* Return the roots */
  plhs[0] = roots;
}
