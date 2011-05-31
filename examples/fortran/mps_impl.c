/*
 * fortran-module.c
 *
 *  Created on: 31/mag/2011
 *      Author: leonardo
 */

#include <mps/interface.h>

/**
 * @brief Caller for mps_solve routine.
 *
 * Since fortran complex are
 * identical to internal mpsolve cplx_t complex type, i.e.
 * two double in a struct, we can simply cast silently the
 * arguments of the fortran routine.
 */
void
mps_roots_impl_ (cplx_t* coeff, cplx_t* roots, int* n)
{
  /* Create a new mps_status */
  mps_status* s = mps_status_new();

  mps_status_set_poly_d(s, coeff, *n);
  mps_mpsolve(s);
  mps_status_get_roots_d(s, roots, NULL);

  mps_status_free (s);
}
