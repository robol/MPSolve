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

#include <mps/mpsxx.h>

namespace mps {

  class NRootsPolynomial : Polynomial {

  public : 
    mps_boolean eval (mps_context * ctx, cplx_t x, cplx_t value, double * error);
    mps_boolean eval (mps_context * ctx, cdpe_t x, cdpe_t value, rdpe_t   error);
    mps_boolean eval (mps_context * ctx, mpc_t x, mpc_t value, rdpe_t error);
  
  };

}


#endif /* MPS_NROOTS_POLYNOMIAL_H_ */

