/*
 * This file is part of MPSolve 3.1.5
 *
 * Copyright (C) 2001-2014, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Leonardo Robol <robol@mail.dm.unipi.it>
 */

#include <mps/mpsxx.h>

using namespace mps;

Polynomial::Polynomial ()
{
  this->feval = &Polynomial::feval_wrapper;
  this->deval = &Polynomial::deval_wrapper;
  this->meval = &Polynomial::meval_wrapper;
}

long int
Polynomial::raise_data (mps_context * ctx, long int wp)
{
  /* This implementation is a no-op by default */
  return mps_polynomial::raise_data(ctx, reinterpret_cast<mps_polynomial*> (this), wp);
}

mps_boolean
Polynomial::feval_wrapper (mps_context * ctx, mps_polynomial * p, 
			   cplx_t x, cplx_t value, double * error)
{
  Polynomial *thisP = reinterpret_cast<Polynomial*> (p);
  return thisP->eval (ctx, x, value, error);
}

mps_boolean
Polynomial::deval_wrapper (mps_context * ctx, mps_polynomial * p, 
			   cdpe_t x, cdpe_t value, rdpe_t error)
{
  Polynomial *thisP = reinterpret_cast<Polynomial*> (p);
  return thisP->eval (ctx, x, value, error);
}

mps_boolean
Polynomial::meval_wrapper (mps_context * ctx, mps_polynomial * p, 
			   mpc_t x, mpc_t value, rdpe_t error)
{
  Polynomial *thisP = reinterpret_cast<Polynomial*> (p);
  return thisP->eval (ctx, x, value, error);
}
