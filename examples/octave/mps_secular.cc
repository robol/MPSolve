#include <octave/oct.h>
#include <gmp.h>
#include <octave/error.h>
#include <mps/mps.h>

#ifndef _MPS_PRIVATE
#define _MPS_PRIVATE
#endif


DEFUN_DLD(mps_secular, args, nargout,
"-*- texinfo -*- \n\
@deftypefn {Loadable Function} {@var{x} =} mps_secular (@var{a}, @var{b}, @var{alg} = 's')\n\
Compute the roots of the secular equation \n\
@tex\n\
  \\sum \\frac {a_i}{x - b_i} - 1 = 0 \n\
@end tex \n\
@ifnotttex \n\
  a_1 / (x - b_1) + ... + a_n / (x - b_n) - 1 = 0 \n\
@end ifnottex \n\
The variable @var{a} contains the values of the numerators, while @var{b} \n\
contains the values of b_i. The n roots are returned in a the output vector \n\
@var{x}.  \n\
The variable @var{alg} can be used to select the algorithm to use. Possible \n\
values are 's' for secsolve and 'u' for unisolve.\n\
@end deftytpefn")
{
  int nargin = args.length();
  mps_algorithm alg = MPS_ALGORITHM_SECULAR_GA;
  
  if (nargin < 2 || nargin > 3)
    {
      print_usage();
      return octave_value_list();
    }

  if (nargin == 3)
    {
      if (args(2).string_value() == std::string("u"))
	alg = MPS_ALGORITHM_STANDARD_MPSOLVE;
      else if (args(2).string_value() != std::string("s"))
	{
	  octave_stdout << "Warning: alg can be one one of 's' or 'u'." << std::endl;
	  print_usage();
	  return octave_value_list();
	}
    }

  if (! (args(1).is_matrix_type() && args(0).is_matrix_type()))
    {
      octave_stdout << "Warning: A and B must be complex vectors." << std::endl;
      print_usage();
      return octave_value_list();
    }

  ComplexColumnVector a = args(0).complex_column_vector_value();
  ComplexColumnVector b = args(1).complex_column_vector_value();

  int n = a.length();

  mps_context * ctx = mps_context_new();
  mps_secular_equation * sec = mps_secular_equation_new (ctx, (cplx_t*) a.fortran_vec(),
							 (cplx_t*) b.fortran_vec(), n);

  mps_context_set_input_poly (ctx, MPS_POLYNOMIAL (sec));
  mps_context_select_algorithm (ctx, alg);

  mps_mpsolve (ctx);
  mps_approximation **apprs = mps_context_get_approximations (ctx);

  ComplexColumnVector roots(n);

  for (int i = 0; i < n; i++)
    {
      cplx_t fvalue;
      mps_approximation_get_fvalue (ctx, apprs[i], fvalue);

      roots(i) = Complex (cplx_Re (fvalue), cplx_Im (fvalue));
      free (apprs[i]);
    }

  free (apprs);

  return octave_value (roots); 
}
