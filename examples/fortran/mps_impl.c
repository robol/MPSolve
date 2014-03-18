/*
 * fortran-module.c
 *
 *  Created on: 31/mag/2011
 *      Author: leonardo
 */

#include <mps/mps.h>

/**
 * @brief Caller for mps_solve routine.
 *
 * Since fortran complex are
 * identical to internal mpsolve cplx_t complex type, i.e.
 * two double in a struct, we can simply cast silently the
 * arguments of the fortran routine.
 */
void
mps_roots_ (int * n, cplx_t * coeff, cplx_t * roots)
{
  /* Create a new mps_context and a new polynomial */
  mps_context *s = mps_context_new ();
  mps_monomial_poly * p = mps_monomial_poly_new (s, *n);
  int i;

  for (i = 0; i <= *n; i++)
    mps_monomial_poly_set_coefficient_d (s, p, i, cplx_Re (coeff[i]), cplx_Im (coeff[i]));
  mps_context_set_input_poly (s, MPS_POLYNOMIAL (p));

  /* Set the output precision to DBL_EPSILON and the default goal
   * to approximate. Try to find all the possible digits representable 
   * in floating point. */
  mps_context_set_output_prec (s, 53);
  mps_context_set_output_goal (s, MPS_OUTPUT_GOAL_APPROXIMATE);

  mps_mpsolve (s);

  mps_context_get_roots_d (s, &roots, NULL);
  mps_context_free (s);
}
