#include "mex.h"
#include <mps/mps.h>
#include <gmp.h>
#include <mpfr.h>
#include "mps_option_parser.h"

void
mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
  int n, i;

  if (nrhs < 2)
    mexErrMsgTxt ("Please select the value of n and alpha\n");

  n = *mxGetPr(prhs[0]);
  double dalpha = *mxGetPr(prhs[1]);

  mps_context * ctx = mps_context_new ();
  mps_monomial_poly * mp = mps_monomial_poly_new (ctx, n);

  mxArray * roots = NULL, *radii = NULL;
  double * real_res = NULL, *imag_res = NULL;
  cplx_t * results = NULL;

  _mps_matlab_options options = mps_parse_matlab_options ( (nrhs > 2) ? prhs[2] : NULL );

  double * perturbations = (nrhs > 3) ? mxGetPr(prhs[3]) : NULL;

  /* Compute exp(- .5 * k^\alpha) */
  mpf_t f; 
  mpc_t d;
  mpfr_t f2, alpha;

  mpf_init2 (f, n * 16);
  mpc_init2 (d, n * 16);
  mpf_set_ui (mpc_Im (d), 0U);
  mpfr_init2 (f2, n * 16);
  mpfr_init2 (alpha, n * 16);

  mpfr_set_d (alpha, dalpha, MPFR_RNDN);
  mpfr_set_ui (f2, 1U, MPFR_RNDN);
  mpfr_get_f (mpc_Re (d), f2, MPFR_RNDN);

  if (perturbations)
    {
      mpf_set_d (f, perturbations[0]);
      mpf_mul_eq (mpc_Re (d), f);
    }
  mps_monomial_poly_set_coefficient_f (ctx, mp, 0, d);
  
  for (i = 1; i <= n; i++)
    {
      double random_perturbation = perturbations ? perturbations[i] : 1.0;

      mpfr_set_ui (f2, i, MPFR_RNDN);
      mpfr_pow (f2, f2, alpha, MPFR_RNDN);
      mpfr_mul_d (f2, f2, - .5, MPFR_RNDN);
      mpfr_exp (f2, f2, MPFR_RNDN);
      mpfr_get_f (mpc_Re (d), f2, MPFR_RNDN);
      
      if (perturbations)
	{
	  mpf_set_d (f, random_perturbation);
	  mpf_mul_eq (mpc_Re (d), f);
	}
      /* This can be used as a crude debug for the computed coefficients */
      /* printf("Coeff: %f\n", pow(mpf_get_d (mpc_Re (d)), 2)); */
      
      mps_monomial_poly_set_coefficient_f (ctx, mp, i, d);
    }
    
  mpc_clear (d);
  mpf_clear (f);
  mpfr_clear (f2);
  mpfr_clear (alpha);

  mps_context_set_input_poly (ctx, MPS_POLYNOMIAL (mp));

  /* Uncomment this to enable debug on the MPSolve operations */
  /* mps_context_add_debug_domain (ctx, MPS_DEBUG_TRACE); */

  mps_context_select_algorithm (ctx, options.algorithm);
  mps_context_set_output_prec (ctx, options.digits / log10(2));
  mps_context_set_output_goal (ctx, options.goal);

  mps_mpsolve (ctx);

  double * _radii = NULL;

  /* Only output the roots as floating point numbers if the
   * requested precision is less than 16 digits. */
  if (options.digits <= DBL_DIG + 1)
    {
      mps_context_get_roots_d (ctx, &results, &_radii);

      roots = mxCreateDoubleMatrix (n, 1, mxCOMPLEX);
      radii = mxCreateDoubleMatrix (n, 1, mxREAL);
      double * radii_ptr = mxGetPr(radii);

      real_res = mxGetPr (roots);
      imag_res = mxGetPi (roots);

      for (i = 0; i < n; i++)
	{
	  real_res[i] = cplx_Re (results[i]);
	  imag_res[i] = cplx_Im (results[i]);

	  radii_ptr[i] = _radii[i];
	}

      /* Free used data */
      cplx_vfree (results);
      cplx_vfree (_radii);
    }
  else
    {
      mpc_t * mresults = NULL;
      rdpe_t * _radii = NULL;
      char * buffer_r, * buffer_i;
      int ndim = 2, dims[] = { n, 4 };

      roots = mxCreateCellArray(ndim, dims);
      radii = mxCreateDoubleMatrix (n, 1, mxREAL);
      double * radii_ptr = mxGetPr(radii);

      mps_context_get_roots_m (ctx, &mresults, &_radii);
      
      for (i = 0; i < n; i++)
	{
	  mp_exp_t rexp, iexp;

	  buffer_r = buffer_i = NULL;
	  buffer_r = mpf_get_str (buffer_r, &rexp, 10, 0, mpc_Re (mresults[i]));
	  buffer_i = mpf_get_str (buffer_i, &iexp, 10, 0, mpc_Im (mresults[i]));

	  mpc_clear (mresults[i]);

	  int i0[] = { i, 0 };
          int i1[] = { i, 1 };
	  int i2[] = { i, 2 };
	  int i3[] = { i, 3 };

	  mxSetCell (roots,
		     mxCalcSingleSubscript(roots, 2, i0),
		     mxCreateString(buffer_r));

	  mxArray * rexpM = mxCreateDoubleMatrix (1, 1, mxREAL);
	  mxGetPr(rexpM)[0] = rexp;
	  mxSetCell (roots,
		     mxCalcSingleSubscript(roots, 2, i1),
		     rexpM);

	  mxSetCell (roots, 
		     mxCalcSingleSubscript(roots, 2, i2),
		     mxCreateString(buffer_i));

	  mxArray * iexpM = mxCreateDoubleMatrix (1, 1, mxREAL);
	  mxGetPr(iexpM)[0] = iexp;

	  mxSetCell (roots, 
		     mxCalcSingleSubscript(roots, 2, i3),
		     iexpM);
	  
	  free (buffer_r);
	  free (buffer_i);

	  radii_ptr[i] = rdpe_get_d (_radii[i]);
	}
      
      free (_radii);
      free (mresults);
    }

  /* Return the roots */
  plhs[0] = roots;
  plhs[1] = radii;

  mps_monomial_poly_free (ctx, MPS_POLYNOMIAL (mp));
  mps_context_free (ctx);
}

