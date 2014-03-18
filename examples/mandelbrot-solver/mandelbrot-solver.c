#include <stdio.h>
#include "mandelbrot-poly.h"

int main (int argc, char * argv[]) {

  /* Create a new mps_context that will be used to solve a Mandelbrot polynomial
   * of the selected degree. */
  mps_context * ctx = mps_context_new ();

  if (argc != 2)
    {
      fprintf (stderr, 
	       "Usage: %s n \n"
	       "\n"
	       "Parameters: \n"
	       " - n is the level of the Mandelbrot polynomial to solve\n\n", 
	       argv[0]);
      return EXIT_FAILURE;
    }

  int n = atoi (argv[1]);

  if (n <= 0) 
    {
      fprintf (stderr, "Please specify a positive integer as Mandelbrot level.\n");
      return EXIT_FAILURE;
    }

  mps_mandelbrot_poly *mp = mps_mandelbrot_poly_new (ctx, n);

  mps_context_set_input_poly (ctx, MPS_POLYNOMIAL (mp));
  mps_context_select_algorithm (ctx, MPS_ALGORITHM_SECULAR_GA);
  mps_context_set_starting_phase (ctx, float_phase);
  mps_mpsolve (ctx);

  mps_approximation ** apprs = mps_context_get_approximations (ctx);
  for (int i = 0; i < mps_context_get_degree (ctx); i++)
    {
      cplx_t value;
      mps_approximation_get_fvalue (ctx, apprs[i], value);
      printf (" "); cplx_out_str (stdout, value); printf ("\n");

      mps_approximation_free (apprs[i]);
    }
  free (apprs);

  mps_mandelbrot_poly_free (ctx, MPS_POLYNOMIAL (mp));
  mps_context_free (ctx);

  return EXIT_SUCCESS;
}

