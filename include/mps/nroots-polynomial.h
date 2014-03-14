/*
 * This file is part of MPSolve 3.1.5
 *
 * Copyright (C) 2001-2014, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Leonardo Robol <robol@mail.dm.unipi.it>
 */

 /**
  * @file
  * @brief
  */

#ifndef MPS_NROOTS_POLYNOMIAL_H_
#define MPS_NROOTS_POLYNOMIAL_H_

#ifdef __cplusplus

#include <mps/mps.h>

namespace mps {

  /**
   * @brief This class is a brief example of how one could implement a custom polynomial
   * type using C++ classes. 
   *
   * It's not really meant to be of any practical use or to be efficient. The main purpose
   * of its implementation is to be straightforward so anyone can use it a a "tutorial"
   * for creating custom polynomial types. 
   */
  class NRootsPolynomial : Polynomial {

  public : 
    /**
     * @brief Create the polynomial \f$x^n - 1\f$. 
     *
     * @param n The degree of the polynomial that should be created. 
     */
    explicit NRootsPolynomial (mps_context * ctx, int n);
    
    mps_boolean eval (mps_context * ctx, cplx_t x, cplx_t value, double * error);
    mps_boolean eval (mps_context * ctx, cdpe_t x, cdpe_t value, rdpe_t   error);
    mps_boolean eval (mps_context * ctx, mpc_t x, mpc_t value, rdpe_t error);

    void newton (mps_context * ctx, mps_approximation * a, cplx_t x);
    void newton (mps_context * ctx, mps_approximation * a, cdpe_t x);
    void newton (mps_context * ctx, mps_approximation * a, mpc_t x, long int wp);
  
  };

}

#endif


#endif /* MPS_NROOTS_POLYNOMIAL_H_ */

