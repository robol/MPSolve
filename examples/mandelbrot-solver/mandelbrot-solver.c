#include "mandelbrot-poly.h"


mps_mandelbrot_poly *
mps_mandelbrot_poly_new (mps_context * ctx, int level)
{
  mps_mandelbrot_poly *mp = mps_new (mps_mandelbrot_poly);
  mps_polynomial *p = MPS_POLYNOMIAL (mp);
  mps_polynomial_init (ctx, p);

  mp->level = level;
  p->thread_safe = true;
  p->degree = pow(2, level) - 1;

  p->type_name = "mps_mandelbrot_poly";
  p->structure = MPS_STRUCTURE_REAL_INTEGER;
  p->density = MPS_DENSITY_DENSE;
  
  /* Load methods */
  p->fnewton = mps_mandelbrot_poly_fnewton;
  p->dnewton = mps_mandelbrot_poly_dnewton;
  p->mnewton = mps_mandelbrot_poly_mnewton;
  p->feval = mps_mandelbrot_poly_feval;
  p->deval = mps_mandelbrot_poly_deval;
  p->meval = mps_mandelbrot_poly_meval;
  p->dstart = mps_mandelbrot_poly_dstart;
  p->fstart = mps_mandelbrot_poly_fstart;

  return mp;
}

int main (int argc, char * argv[]) {

  mps_context * ctx = mps_context_new ();
  mps_mandelbrot_poly *mp = mps_mandelbrot_poly_new (ctx, atoi(argv[1]));

  mps_context_set_input_poly (ctx, MPS_POLYNOMIAL (mp));
  mps_context_select_algorithm (ctx, MPS_ALGORITHM_SECULAR_GA);
  mps_context_set_starting_phase (ctx, float_phase);
  mps_mpsolve (ctx);

  mps_context_set_output_format (ctx, MPS_OUTPUT_FORMAT_GNUPLOT_FULL);
  mps_output (ctx);

  return EXIT_SUCCESS;
}

