#ifndef _MPS_PRIVATE
#define _MPS_PRIVATE
#endif

#include <octave/oct.h>
#include <gmp.h>
#include <octave/error.h>
#include <mps/mps.h>

DEFUN_DLD(mps_det, args, nargout,
"-*- texinfo -*- \n\
@deftypefn {Loadable Function} {@var{x} =} mps_det (@var{A})\n\
Compute the determinant of @var{A} assuming it is a Hessenberg matrix. \n\
@end deftypefn")
{
    int nargin = args.length();
    octave_value_list retval;
    const char* params;
    cplx_t det;

    if (nargin != 1) {
        print_usage ();
        return retval;
    }

    ComplexMatrix A = args(0).complex_matrix_value();

    if (A.rows() != A.cols()) {
        octave_stdout << "The input matrix must be square!" << std::endl;
        return retval;
    }

    mps_context * ctx = mps_context_new ();

    mps_fhessenberg_determinant (ctx,
    		(cplx_t*) A.transpose().fortran_vec(),
    		A.rows (), det);

    mps_context_free (ctx);

    return octave_value (Complex(cplx_Re (det), cplx_Im (det)));
}
