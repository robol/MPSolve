/**
 * @brief MATLAB / Octave wrapper around the rational polynomial 
 * solver implemented in rational-* files. 
 *
 * This function find the roots of a rational function defined as
 * 
 *   f(x) = \sum_{i = 1}^n \frac{c_i}{(a_i - b_ix)^p}
 *
 */

#include "mex.h"
#include "rational-poly.h"

void mexFunction(int nlhs, mxArray * plhs[], 
		 int nrhs, const mxArray * prhs[])
{
  int i;

  if (nrhs != 4)
    {
      mexErrMsgTxt("Please call this function with 4 input arguments.");
    }

  int p = *mxGetPr (prhs[3]);
  int n = mxGetM(prhs[0]);

  if (n == 1)
    n = mxGetN(prhs[0]);

  cplx_t * a = cplx_valloc (n);
  cplx_t * b = cplx_valloc (n);
  cplx_t * c = cplx_valloc (n);

  double * ma = mxGetPr (prhs[0]);
  double * mb = mxGetPr (prhs[1]);
  double * mc = mxGetPr (prhs[2]);

  for (i = 0; i < n; i++)
    {
      cplx_set_d (a[i], ma[i], 0.0);
      cplx_set_d (b[i], mb[i], 0.0);
      cplx_set_d (c[i], mc[i], 0.0);
    }

  mps_context * ctx = mps_context_new ();
  mps_rational_poly * rp = mps_rational_poly_new (ctx, a, b, c, n, p);

  mps_context_set_input_poly (ctx, MPS_POLYNOMIAL (rp));
  /* Enable this flag in case you want to see the debugging 
   * output of the algorithm.
   * mps_context_add_debug_domain (ctx, MPS_DEBUG_TRACE);*/
  mps_context_select_algorithm (ctx, MPS_ALGORITHM_SECULAR_GA);
  mps_context_set_output_goal (ctx, MPS_OUTPUT_GOAL_APPROXIMATE);
  mps_mpsolve (ctx);

  cplx_t * roots = NULL;

  mps_context_get_roots_d (ctx, &roots, NULL);

  plhs[0] = mxCreateDoubleMatrix ((n - 1) * p, 1, mxCOMPLEX);

  double * rr = mxGetPr (plhs[0]);
  double * ri = mxGetPi (plhs[0]);
  
  for (i = 0; i < (n-1) * p; i++)
    {
      rr[i] = cplx_Re (roots[i]);
      ri[i] = cplx_Im (roots[i]);
    }

  free (roots);
  
  mps_rational_poly_free (ctx, rp);
  mps_context_free (ctx);
}

