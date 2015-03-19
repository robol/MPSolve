/*
 * This file is part of MPSolve 3.1.5
 *
 * Copyright (C) 2001-2015, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Leonardo Robol <robol@mail.dm.unipi.it>
 */

 /**
  * @file
  * @brief
  */

#include <mps/mps.h>

#ifndef MPS_REGENERATION_DRIVER_H_
#define MPS_REGENERATION_DRIVER_H_

MPS_BEGIN_DECLS

/**
 * @brief This type represent an abstract implementation of a driver for the
 * regeneration step of the main algorithm. 
 *
 * A standard implementation is given for polynomals internally in the MPSolve
 * code but a custom implementation may be given by the user by calling
 * mps_context_set_regeneration_driver(). 
 */
struct mps_regeneration_driver {

  mps_boolean (*update_fsecular_equation) (mps_context * ctx, 
					   mps_polynomial * p, 
					   mps_approximation ** approximations, 
					   mps_secular_equation * old);

  mps_boolean (*update_dsecular_equation) (mps_context * ctx,
					   mps_polynomial * p, 
					   mps_approximation ** approximations,
					   mps_secular_equation * old);

  mps_boolean (*update_msecular_equation) (mps_context * ctx, 
					   mps_polynomial * p, 
					   mps_approximation ** approximations,
					   mps_secular_equation * old);

  /**
   * @brief Optional function that is called by the mps_regeneration_driver_free()
   * method. It is intended for custom regeneration driver that needs to free
   * additional data before having the function destroyed. 
   *
   * Note that the regeneration_driver get destroyed with the context. If you don't
   * need to free data or if you want it to survive the lifecycle of the context, set
   * this pointer to NULL. 
   */
  void (*free)(mps_context * ctx, mps_regeneration_driver * rd);

};

mps_boolean
mps_standard_regeneration_driver_update_fsecular_equation (mps_context * s, 
							   mps_polynomial * p, 
							   mps_approximation ** approximations, 
							   mps_secular_equation * sec);

mps_boolean
mps_standard_regeneration_driver_update_dsecular_equation (mps_context * s, 
							   mps_polynomial * p, 
							   mps_approximation ** approximations, 
							   mps_secular_equation * sec);

mps_boolean
mps_standard_regeneration_driver_update_msecular_equation (mps_context * s, 
							   mps_polynomial * p, 
							   mps_approximation ** approximations, 
							   mps_secular_equation * sec);

mps_regeneration_driver *
mps_regeneration_driver_new_standard (mps_context * ctx);

void
mps_regeneration_driver_free (mps_context * ctx, mps_regeneration_driver * rd);

MPS_END_DECLS

#endif /* MPS_REGENERATION_DRIVER_H_ */

