#ifndef __MPS_MONOMIAL_POLY_H
#define __MPS_MONOMIAL_POLY_H

/**
 * @file
 * @brief Implementation of the allocation and edit functions for the 
 * handling of monomial polynomials. 
 */

#ifdef	__cplusplus
extern "C"
{
#endif

  #include <mps/core.h>
  #include <gmp.h>

  /* These routines are thought for polynomial handling, i.e. allocating and 
   * setting coefficients of the polynomials, and setting the precision of the
   * floating point coefficients that are in there */

  mps_monomial_poly * mps_monomial_poly_new (mps_status * s, long int degree);

  void mps_monomial_poly_free (mps_status * s, mps_monomial_poly * mp);

  void mps_monomial_poly_raise_precision (mps_status * s, mps_monomial_poly * mp, long int prec);

  void mps_monomial_poly_set_coefficient_q (mps_status * s, mps_monomial_poly * mp, long int i, 
					    mpq_t real_part, mpq_t imag_part);


#ifdef	__cplusplus
}
#endif

#endif
