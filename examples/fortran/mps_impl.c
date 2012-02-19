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
mps_roots_impl_ (cplx_t * coeff, cplx_t * roots, int *n)
{
  /* Create a new mps_status and a new polynomial */
  mps_status *s = mps_status_new ();
  mps_monomial_poly * p = mps_monomial_poly_new (s, *n);
  int i;

  for (i = 0; i <= *n; i++)
    mps_monomial_poly_set_coefficient_d (s, p, i, cplx_Re (coeff[i]), cplx_Im (coeff[i]));
  mps_status_set_input_poly (s, p);

  /* mps_status_select_algorithm (s, MPS_ALGORITHM_SECULAR_GA); */

  /* Set the output precision to DBL_EPSILON and the default goal
   * to approximate. */
  mps_status_set_output_prec (s, 53);
  mps_status_set_output_goal (s, MPS_OUTPUT_GOAL_APPROXIMATE);

  mps_mpsolve (s);

  mps_status_get_roots_d (s, roots, NULL);
  mps_status_free (s);
}
