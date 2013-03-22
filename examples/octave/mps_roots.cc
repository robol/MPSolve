#include <octave/oct.h>
#include <gmp.h>
#include <octave/error.h>
#include <mps/mps.h>


DEFUN_DLD(mps_roots, args, nargout,
"-*- texinfo -*- \n\
@deftypefn {Loadable Function} {@var{x} =} mps_roots(@var{v})\n\
@deftypefnx {Loadable Function} {@var{x} =} mps_roots(@var{v}, @var{alg})\n\
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
The optional variable @var{alg} can be set to \"s\" or \"u\" to select \
the secular algorithm or the standard MPSolve algorithm. The default value \
is \"s\"\
@end deftypefn")
{
    int nargin = args.length();
    octave_value_list retval;
    const char* params;

    mps_algorithm algorithm = MPS_ALGORITHM_SECULAR_GA;

    if (nargin < 1 || nargin > 3) {
        print_usage ();
        return retval;
    }

    if (nargin == 2) {
      if (!args(1).is_string() || args(1).length() != 1 || 
	  ((args(1).string_value() != std::string("s")) && 
	   (args(1).string_value() != std::string("u")))) {
	print_usage();
	return retval;
      }
      else {
	if (args(1).string_value() == std::string("u"))
	  algorithm = MPS_ALGORITHM_STANDARD_MPSOLVE;
      }
    }

    /* Check that input data is a vector */
    if (!args(0).is_complex_matrix() && !args(0).is_real_matrix() && !(args(0).is_int64_type() && args(0).is_matrix_type())) {
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
    if (rows > 1) {
        parsed_args = args(0).complex_matrix_value().transpose();
    } else
        parsed_args = args(0).complex_matrix_value();

    /* This is the vector of che coefficients of the polynomial */
    ComplexColumnVector v(parsed_args);

    /* Get its dimension */
    octave_idx_type n = cols;

    /* Create vector of the roots */
    ComplexColumnVector res(n - 1);
    cplx_t *results = cplx_valloc(n-1);

    mps_context* s = mps_context_new();
    mps_context_select_algorithm (s, algorithm);

    mps_monomial_poly * p = mps_monomial_poly_new (s, n - 1);

    if (args(0).is_int64_type()) {
      int64NDArray real_coeffs = args(0).real().int64_array_value();
      int64NDArray imag_coeffs = args(0).imag().int64_array_value();
      for (int i = 0; i < n; i++)
	mps_monomial_poly_set_coefficient_int (s, p, (n - i - 1), (int64_t) (real_coeffs(i)), (int64_t) (imag_coeffs(i)));
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
