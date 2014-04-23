#include <octave/oct.h>
#include <octave/ov-struct.h>
#include <gmp.h>
#include <octave/error.h>
#include <mps/mps.h>

#ifndef _MPS_PRIVATE
#define _MPS_PRIVATE
#endif

extern "C" {
  cplx_t * starting_points_vec = NULL;

  void _custom_start_function (mps_context * ctx, mps_polynomial * p, mps_approximation ** apprs)
  {
    int i, n = mps_context_get_degree (ctx);

    for (i = 0; i < n; i++)
      {
	mps_approximation_set_fvalue (ctx, apprs[i], starting_points_vec[i]);
      }
  }

  void _custom_dpe_start_function (mps_context * ctx, mps_polynomial * p, mps_approximation ** apprs)
  {
    int i, n = mps_context_get_degree (ctx);

    for (i = 0; i < n; i++)
      {
	cdpe_t t;
	cdpe_set_x (t, starting_points_vec[i]);
	mps_approximation_set_dvalue (ctx, apprs[i], t);
      }
  }
}


DEFUN_DLD(mps_roots, args, nargout,
"-*- texinfo -*- \n\
@deftypefn {Loadable Function} {@var{x} =} mps_roots (@var{v})\n\
@deftypefnx {Loadable Function} {@var{x} =} mps_roots (@var{v}, @var{options})\n\
@cindex root finding of a polynomial\n\
Compute the roots of the polynomial p(z) given by\n\n\
@tex\n\
$$\n p(z) = v_1z^{n-1} + v_2 z^{n-2} + \\ldots + + v_{n-1} z + v_n\n $$ \n\
@end tex\n\
@ifnottex\n\
@example\n\
        p(z) = v(1) * z^(N-1) + ... + v(N-1) * z + v(N)\n\
@end example\n\
@end ifnottex\n\n\
and return a vector with the roots.\n\n\
The optional variable @var{options} can be set to \"s\" or \"u\" to select \n\
the secular algorithm or the standard MPSolve algorithm. The default value \n\
is \"s\".  \n\
\n\
It can also be set to a struct that can contain zero or more of the following \n\
fields: \n\n\
 algorithm: can be \"s\" or \"u\" and has the same meaning of above. \n\n\
 starting_points: must be a vector of length n with the starting points for MPSolve\n\n\
@end deftypefn")
{
    int nargin = args.length();
    octave_value_list retval;
    const char* params;
    ComplexColumnVector starting_points;
    bool customStartFunction = false;

    mps_algorithm algorithm = MPS_ALGORITHM_SECULAR_GA;

    if (nargin < 1 || nargin > 3) {
        print_usage ();
        return retval;
    }

    if (nargin == 2) {
      std::string algorithm_s;
      octave_map smap;

      /* If the string conversion triggers an error try to recover
       * a struct. */
      if (args(1).is_map())
	smap = args(1).map_value();
      else
	{
	  if (! args(1).is_string())
	    {
	      print_usage();
	      return retval;
	    }

	  Cell a(1,1);
	  a(0) = args(1);

	  smap.assign ("algorithm", a);
	}

      /* Parse the arguments that are stored in the struct */
      if (smap.contains ("algorithm"))
	{
	  algorithm_s = smap.getfield ("algorithm")(0).string_value();

	  if (algorithm_s != std::string("s") && algorithm_s != std::string("u"))
	    {
	      print_usage();
	      return retval;
	    }
	  else if (algorithm_s == std::string("u"))
	    algorithm = MPS_ALGORITHM_STANDARD_MPSOLVE;
	}

      if (smap.contains ("starting_points"))
	{
	  starting_points = smap.getfield ("starting_points")(0).complex_vector_value();
	  customStartFunction = true;
	  starting_points_vec = (cplx_t*) starting_points.fortran_vec();
	}
    }
      
    /* Check that input data is a vector */
    if (!args(0).is_complex_matrix() && !args(0).is_real_matrix() && 
	!(args(0).is_int64_type() && args(0).is_matrix_type())) {
        print_usage ();
        return retval;
    }

    /* Check input data shape */
    octave_idx_type rows = args(0).rows();
    octave_idx_type cols = args(0).columns();
    if(rows != 1 && cols != 1) {
        print_usage ();
        return retval;
    }

    if (nargout > 1) {
        print_usage ();
        return retval;
    }

    ComplexMatrix parsed_args;
    if (rows > 1) 
        parsed_args = args(0).complex_matrix_value().transpose();
    else
        parsed_args = args(0).complex_matrix_value();

    /* This is the vector of che coefficients of the polynomial */
    ComplexColumnVector v = args(0).complex_column_vector_value();

    /* Get its dimension */
    octave_idx_type n = rows < cols ? cols : rows;

    /* Create vector of the roots */
    ComplexColumnVector res(n - 1);
    cplx_t *results = cplx_valloc(n-1);

    mps_context* s = mps_context_new();
    mps_context_select_algorithm (s, algorithm);

    mps_monomial_poly * p = mps_monomial_poly_new (s, n - 1);

    /* If a custom start function was specified, replace the default one. */
    if (customStartFunction)
      {
	mps_polynomial * poly = MPS_POLYNOMIAL (p);
	poly->fstart = _custom_start_function;
	poly->dstart = _custom_dpe_start_function;

	if (starting_points.length() != n -1)
	  {
	    octave_stdout << "Warning: The length of the starting points should be equal to the degree" << std::endl;
	    return retval;
	  }
      }

    if (args(0).is_int64_type()) {
      int64NDArray real_coeffs = args(0).real().int64_array_value();
      int64NDArray imag_coeffs = args(0).imag().int64_array_value();
      for (int i = 0; i < n; i++)
	mps_monomial_poly_set_coefficient_int (s, p, (n - i - 1), (int64_t) (real_coeffs(i)), 
					       (int64_t) (imag_coeffs(i)));
   }
    else {
      for(int i = 0; i < n; i++) {
	mps_monomial_poly_set_coefficient_d (s, p, n - i - 1, v(i).real(), v(i).imag());
      }
    }

    mps_context_set_input_poly (s, MPS_POLYNOMIAL (p));
    mps_context_set_output_goal (s, MPS_OUTPUT_GOAL_APPROXIMATE);

    /* Actually solve it */
    mps_mpsolve(s);

    /* Get roots and return them */
    mps_context_get_roots_d(s, &results, NULL);
    for(int i = 0; i < n - 1; i++) {
        res(i) = Complex(cplx_Re(results[i]), cplx_Im(results[i]));
    }
    cplx_vfree (results);

    /* Free mpsolve status */
    mps_context_free (s);

    return octave_value(res);
}
