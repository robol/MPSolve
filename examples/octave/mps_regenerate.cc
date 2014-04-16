#ifndef _MPS_PRIVATE
#define _MPS_PRIVATE
#endif

#include <octave/oct.h>
#include <gmp.h>
#include <octave/error.h>
#include <mps/mps.h>

DEFUN_DLD(mps_regenerate, args, nargout,
"-*- texinfo -*- \n\
@deftypefn {Loadable Function} {@var{a} =} mps_det (@var{b}, @var{P0}, ..., @var{PN})\n\
Compute a secular equation equivalent to the polynomial \n\
\n\
 P(x) = det(P0 + x*P1 + ... + x^n*Pn) \n\
\n\
This function will take that b_i chosen as nodes and return the respective a_i such that \n\
\n\
 sum (a ./ (x - b)) - 1 = 0 if and only if P(x) = 0 \n\
\n\
@end deftypefn")
{
  int nargin = args.length();
  mps_context * ctx = mps_context_new ();

  // Input sanitizing goes here

  ComplexColumnVector b = args(0).complex_column_vector_value();
  int degree = nargin - 2;
  ComplexColumnVector a(b.length());

  cplx_t * fa = (cplx_t*) a.fortran_vec();
  cplx_t * fb = (cplx_t*) b.fortran_vec();

  cdpe_t ctmp;

  // Handle some special cases here. Here we are assuming that the second argument
  // is diagonal then it's the identity. 
  if (degree == 1 && args(2).is_diag_matrix())
    {
      ComplexMatrix A = args(1).complex_matrix_value().transpose();
      cplx_t * A_data = (cplx_t*) A.fortran_vec();
      
      for (int i = 0; i < a.length(); i++)
	{
	  long int exponent = 0;
	  cplx_t prod;

	  mps_fhessenberg_shifted_determinant (ctx, A_data, fb[i], A.rows(), fa[i], &exponent);
	  cplx_set (prod, cplx_one);

	  for (int j = 0; j < b.length(); j++)
	    {
	      if (i == j)
		continue;

	      cplx_t diff;
	      cplx_sub (diff, fb[i],fb[j]);
	      cplx_mul_eq (prod, diff);

	      if (j % 50 == 0)
		{
		  int pe;
		  double mod = cplx_mod (prod);
		  double new_mod = frexp (mod, &pe);

		  exponent -= pe;
		  cplx_div_eq_d (prod, mod / new_mod);
		}
	    }

	  cplx_div_eq (fa[i], prod);
	  cdpe_set_2dl (ctmp, 1.0, exponent, 0.0, 1);
	  cdpe_get_x (prod, ctmp);
	  cplx_mul_eq (fa[i], prod);
	}
    }

  mps_context_free (ctx);
  return octave_value (- a);
}
