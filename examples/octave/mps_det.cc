#ifndef _MPS_PRIVATE
#define _MPS_PRIVATE
#endif

#include <octave/oct.h>
#include <gmp.h>
#include <octave/error.h>
#include <mps/mps.h>

DEFUN_DLD(mps_det, args, nargout,
"-*- texinfo -*- \n\
@deftypefn {Loadable Function} {@var{x} =} mps_det (@var{A}, @var{mode} = 'f')\n\
Compute the determinant of @var{A} assuming it is a Hessenberg matrix. \n\
The value of @var{mode} can be one of 'f' or 'd'. 'f' is the default and implies\n\
that standard floating point will be use, while 'd' will instruct MPSolve to perform\n\
the computation in DPE to ensure that no over/underflow is encountered. \n\
@end deftypefn")
{
    int nargin = args.length();
    octave_value_list retval;
    const char* params;
    cplx_t det;
    mps_boolean use_dpe = false;

    if (nargin > 2) {
        print_usage ();
        return retval;
    }

    if (nargin == 2 && args(1).string_value() == std::string("d"))
      use_dpe = true;

    ComplexMatrix A = args(0).complex_matrix_value();

    if (A.rows() != A.cols()) {
        octave_stdout << "The input matrix must be square!" << std::endl;
        return retval;
    }

    mps_context * ctx = mps_context_new ();

    if (! use_dpe)
      mps_fhessenberg_determinant (ctx,
				   (cplx_t*) A.transpose().fortran_vec(),
				   A.rows (), det);
    else
      {
	int i, j;
	size_t n = A.rows();
	cplx_t* fp_matrix = (cplx_t*) A.fortran_vec();
	cdpe_t output;

	/* In this case we need to copy the matrix elements, since we cannot directly
	 * pass a floating point matrix to a DPE routine. */
	cdpe_t * matrix = mps_newv (cdpe_t, n * n);
	for (i = 0; i < n; i++)
	  for (j = 0; j < n; j++)
	    {
	      cdpe_set_x (MPS_MATRIX_ELEM(matrix, i, j, n), 
			  MPS_MATRIX_ELEM(fp_matrix, j, i, n));
	    }

	mps_dhessenberg_determinant (ctx, matrix, n, output);	
	cdpe_get_x (det, output);				     
	free (matrix);
      }

    mps_context_free (ctx);

    return octave_value (Complex(cplx_Re (det), cplx_Im (det)));
}
