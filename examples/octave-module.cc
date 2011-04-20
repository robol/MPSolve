#include <octave/oct.h>

extern "C" {
	#include <mps/mps.h>
}

void
usage() {
    return;
}

DEFUN_DLD(mps_roots, args, nargout,
"-*- texinfo -*- \n\
@deftypefn {Loadable Function} {@var{x} =} mps_roots(@var{v})\n\
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
and return a vector with the roots.\n\
@end deftypefn")
{
    /* Check that input data is a vector */
    if (!args(0).is_complex_matrix() && !args(0).is_real_matrix()) {
        usage();
        error("A vector must be passed as first argument.");
    }

    /* Check input data shape */
    octave_idx_type rows = args(0).rows();
    octave_idx_type cols = args(0).columns();
    if(rows != 1 && cols != 1) {
        usage();
        error("A vector must be passed as first argument.");
    }

    if (nargout > 1) {
        usage();
        error("This function returns a single value.");
    }

    ComplexMatrix parsed_args;
    if (rows > 1) {
        parsed_args = args(0).complex_matrix_value().transpose();
    } else
        parsed_args = args(0).complex_matrix_value();

    /* This is the vector of che coefficients of the polynomial */
    ComplexColumnVector v(parsed_args);

    /* Get its dimension */
    octave_idx_type n = v.cols();

    /* Create vector of the roots */
    ComplexColumnVector res(n - 1);

    /* Allocate coefficients and space for results to be passed
     * to the function mps_mpsolve */
    cplx_t *coeff = cplx_valloc(n);
    cplx_t *results = cplx_valloc(n-1);

    for(int i = 0; i < n; i++) {
        cplx_set_d(coeff[n - i - 1], v(i).real(), v(i).imag());
    }

    /* Create mps_status* struct */
    mps_status* s = mps_status_new();
    mps_status_set_poly_d(s, coeff, n - 1);

    /* Actually solve it */
    mps_mpsolve(s);

    /* Get roots and return them */
    mps_get_roots_d(s, results, NULL);
    for(int i = 0; i < n - 1; i++) {
        res(i) = Complex(cplx_Re(results[i]), cplx_Im(results[i]));
    }

    /* Free used arrays */
    cplx_vfree (coeff);
    cplx_vfree (results);

    /* Free mpsolve status */
    mps_status_free (s);

    return octave_value(res);
}
