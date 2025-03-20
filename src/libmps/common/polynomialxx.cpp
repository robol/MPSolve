/*
 * This file is part of MPSolve 3.2.2
 *
 * Copyright (C) 2001-2014, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Leonardo Robol <robol@mail.dm.unipi.it>
 */

#include <mps/mps.h>
#include <mps/private/system/memory-file-stream.h>

using namespace mps;

static void
setup_vtable (mps_context * ctx, Polynomial * self, const char * type_name)
{
  self->type_name = type_name;

  /* Hooking up our wrappers to the actual function pointer of 
   * the mps_polynomial struct. Self make possible to transparently
   * use self object as a custom polynomial type in MPSolve's sense. */
  self->feval = &Polynomial::feval_wrapper;
  self->deval = &Polynomial::deval_wrapper;
  self->meval = &Polynomial::meval_wrapper;

  self->fstart = &Polynomial::fstart_wrapper;
  self->dstart = &Polynomial::dstart_wrapper;
  self->mstart = &Polynomial::mstart_wrapper;

  self->raise_data = &Polynomial::raise_data_wrapper;
  self->free = &Polynomial::free_wrapper;

  self->fnewton = &Polynomial::fnewton_wrapper;
  self->dnewton = &Polynomial::dnewton_wrapper;
  self->mnewton = &Polynomial::mnewton_wrapper;
}

static Polynomial*
polynomial_from_stream (mps_context * ctx, mps::AbstractInputStream * stream)
{
  mps_polynomial * poly = mps_monomial_yacc_parser (ctx, reinterpret_cast<mps_abstract_input_stream *> (stream));  
  return (Polynomial*) (poly);
}

Polynomial *
Polynomial::fromString (mps_context * ctx, const char * inputString)
{
  MemoryFileStream stream(inputString);
  return polynomial_from_stream(ctx, &stream);
}


Polynomial::Polynomial (mps_context * ctx, const char * type_name)
{
  mps_polynomial * self = static_cast<mps_polynomial*> (this);
  mps_polynomial_init (ctx, self);
  setup_vtable (ctx, reinterpret_cast<Polynomial*> (self), "monomial");
}

long int
Polynomial::raise_data_wp (mps_context * ctx, long int wp)
{
  /* This implementation is a no-op by default */
  return wp;
}

Polynomial::~Polynomial()
{
}

void
Polynomial::start_fp (mps_context * ctx, mps_approximation ** approximations)
{
  mps_general_fstart (ctx, static_cast<mps_polynomial*> (this), approximations);
}

void
Polynomial::start_dpe (mps_context * ctx, mps_approximation ** approximations)
{
  mps_general_dstart (ctx, static_cast<mps_polynomial*> (this), approximations);
}

void
Polynomial::start_mp (mps_context * ctx, mps_approximation ** approximations)
{
  mps_general_mstart (ctx, static_cast<mps_polynomial*> (this), approximations);
}

void
Polynomial::get_leading_coefficient (mps_context * ctx, mpc_t lc)
{
  mpc_set_ui (lc, 1U, 0U);
}

mps_boolean
Polynomial::feval_wrapper (mps_context * ctx, mps_polynomial * p, 
			   cplx_t x, cplx_t value, double * error)
{
  Polynomial *thisP = static_cast<Polynomial*> (p);
  return thisP->eval (ctx, x, value, error);
}

mps_boolean
Polynomial::deval_wrapper (mps_context * ctx, mps_polynomial * p, 
			   cdpe_t x, cdpe_t value, rdpe_t error)
{
  Polynomial *thisP = static_cast<Polynomial*> (p);
  return thisP->eval (ctx, x, value, error);
}

mps_boolean
Polynomial::meval_wrapper (mps_context * ctx, mps_polynomial * p, 
			   mpc_t x, mpc_t value, rdpe_t error)
{
  Polynomial *thisP = static_cast<Polynomial*> (p);
  return thisP->eval (ctx, x, value, error);
}

void 
Polynomial::free_wrapper (mps_context * ctx, mps_polynomial * p)
{
  Polynomial *thisP = static_cast<Polynomial*> (p);
  delete thisP;
}
    
long int
Polynomial::raise_data_wrapper (mps_context * ctx, mps_polynomial * p, 
				long int wp)
{
  Polynomial *thisP = static_cast<Polynomial*> (p);
  return thisP->raise_data_wp (ctx, wp);
}

void
Polynomial::fstart_wrapper (mps_context * ctx, mps_polynomial * p, mps_approximation ** approximations)
{
  Polynomial *thisP = static_cast<Polynomial*> (p);
  thisP->start_fp (ctx, approximations);
}

void 
Polynomial::dstart_wrapper (mps_context * ctx, mps_polynomial * p, mps_approximation ** approximations)
{
  Polynomial *thisP = static_cast<Polynomial*> (p);
  thisP->start_dpe (ctx, approximations);
}

void
Polynomial::mstart_wrapper (mps_context * ctx, mps_polynomial * p, mps_approximation ** approximations)
{
  Polynomial *thisP = static_cast<Polynomial*> (p);
  thisP->start_mp (ctx, approximations);
}

void
Polynomial::fnewton_wrapper (mps_context * ctx, mps_polynomial * p, 
			     mps_approximation * a, cplx_t x)
{
  Polynomial *thisP = static_cast<Polynomial*> (p);
  thisP->newton (ctx, a, x);
}

void
Polynomial::dnewton_wrapper (mps_context * ctx, mps_polynomial * p,
			     mps_approximation * a, cdpe_t x)
{
  Polynomial *thisP = static_cast<Polynomial*> (p);
  thisP->newton (ctx, a, x);
}

void
Polynomial::mnewton_wrapper (mps_context * ctx, mps_polynomial * p, 
			     mps_approximation * a, mpc_t x, 
			     long int wp)
{
  Polynomial *thisP = static_cast<Polynomial*> (p);
  thisP->newton (ctx, a, x, wp);
}

void
Polynomial::get_leading_coefficient_wrapper (mps_context * ctx, mps_polynomial * p,
					     mpc_t leading_coefficient)
{
  Polynomial *thisP = static_cast<Polynomial*> (p);
  thisP->get_leading_coefficient (ctx, leading_coefficient);
}
