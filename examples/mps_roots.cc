#include <octave/oct.h>

extern "C" {
	#include <mps/mps.h>
}

DEFUN_DLD(mps_roots, args, nargout,
        "Find the roots of the polynomial v(1) + v(2)z + ... + v(n)z^{n-1}")
{
    /* This is the vector of che coefficients of the polynomial */
    ColumnVector v(args(0).vector_value());

    /* Get its dimension */
    octave_idx_type n = v.cols();

    /* Create vector of the roots */
    ComplexColumnVector res(n - 1);

    /* Allocate coefficients and space for results to be passed
     * to the function mps_mpsolve */
    cplx_t *coeff = cplx_valloc(n);
    cplx_t *results = cplx_valloc(n-1);

    for(int i = 0; i < n; i++) {
        cplx_set_d(coeff[i], v(i), 0);
    }

    /* Create mps_status* struct */
    mps_status* s = mps_status_new();
    mps_status_set_poly_d(s, coeff, n - 1);

    /* Set output precision */
    s->prec_out = 15;

    /* Actually solve it */
    mps_mpsolve(s);

    /* Get roots and return them */
    mps_get_roots_d(s, results, NULL);
    for(int i = 0; i < n - 1; i++) {
        res(i) = Complex(cplx_Re(results[i]), cplx_Im(results[i]));
    }

    return octave_value(res);
}
