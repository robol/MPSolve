#include <octave/oct.h>
#include <gmp.h>
#include <octave/error.h>
#include <mps/mps.h>
#include "octave_support.h"

DEFUN_DLD(mps_polyeig_impl, args, nargout,
"-*- texinfo -*- \n\
@deftypefn {Loadable Function} {@var{LAMBDA} =} mps_polyeig_impl (@var{P0}, @var{P1}, ..., @var{Pn}, @var{alg})\n\
  This is an internal function of MPSolve, and as such is undocumented. Please\n\
  do not call this function directly, but use mps_polyeig instead. \n\
@end deftypefn")
{
  int nargin = args.length ();
  int degree = nargin - 2; 
  cplx_t *roots = NULL; 
  
  /* Create a mps_monomial_matrix_poly */
  mps_context *ctx = mps_context_new (); 
  mps_monomial_matrix_poly *mp = mps_monomial_matrix_poly_new (ctx, degree, args(0).rows(), true); 

  ComplexColumnVector res(degree * args(0).rows ());

  int rows, columns; 

  for (int i = 0; i <= degree; i++) {
    ComplexMatrix coeff = args(i).complex_matrix_value ().transpose (); 

    if (i == 0) {
      rows = coeff.rows (); 
      columns = coeff.cols (); 
    }
    else {
      if (coeff.rows () != rows || coeff.cols () != columns) {
	octave_stdout << "Dimensions mismatch in the coefficients" << std::endl; 
	goto cleanup;
      }
    }

    /* This is more or less guaranteed to work by C++01 standard, not so guaranteed by C++03, 
     * apparently, but will reasonable work on _most_ platform. */
    mps_monomial_matrix_poly_set_coefficient_d (ctx, mp, i, (cplx_t*) (coeff.fortran_vec ()));
  }

  /* If we have an Hessenberg matrix poly we can set the relative flags */
  if (args(degree+1).string_value ().find ('h') != std::string::npos) {
    mps_monomial_matrix_poly_add_flags (ctx, mp, MPS_MONOMIAL_MATRIX_POLY_HESSENBERG); 
  }

  if (args(degree+1).string_value ().find ('d') != std::string::npos)
    mps_context_add_debug_domain (ctx, MPS_DEBUG_INFO);

  mps_context_set_input_poly (ctx, MPS_POLYNOMIAL (mp)); 
  mps_context_select_algorithm (ctx, MPS_ALGORITHM_SECULAR_GA); 
  mps_mpsolve (ctx); 

  mps_context_get_roots_d (ctx, &roots, NULL); 
  for (int i = 0; i < mps_context_get_degree (ctx); i++) {
    res(i) = Complex (cplx_Re(roots[i]), cplx_Im(roots[i]));
  }

  cplx_vfree (roots);

 cleanup:
  mps_monomial_matrix_poly_free (ctx, MPS_POLYNOMIAL (mp)); 
  mps_context_free (ctx); 

  /* This means nothing implemented right now */
  return octave_value (res);
}
