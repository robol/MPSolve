#include <octave/oct.h>

extern "C" {
	#include <mps/mps.h>
}

DEFUN_DLD(mps_roots,
	  args, nargout,
	  "find roots of the polynomial")
{
  ColumnVector v(args(0).vector_value());
  octave_idx_type n = v.cols();
  ComplexColumnVector res(n - 1);

  cplx_t *coeff = cplx_valloc(n);
  cplx_t *results = cplx_valloc(n);
  
  for(int i = 0; i < n; i++) {
    cplx_set_d(coeff[i], v(i), 0);
  }

  mps_status* s = mps_status_new();
  mps_status_set_poly_d(s, coeff, n - 1);
  s->prec_out = 53;

  /* Actually solve it */
  mps_mpsolve(s);
  mps_get_roots_d(s, results, NULL);

  for(int i = 0; i < n - 1; i++) {
    res(i) = Complex(cplx_Re(results[i]), cplx_Im(results[i]));
  }

  return octave_value(res);
}
