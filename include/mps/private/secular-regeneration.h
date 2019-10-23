/*
 * This file is part of MPSolve 3.1.8
 *
 * Copyright (C) 2001-2019, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Leonardo Robol <robol@mail.dm.unipi.it>
 */

 /**
  * @file
  * @brief
  */

#ifndef MPS_SECULAR_REGENERATION_H_
#define MPS_SECULAR_REGENERATION_H_

MPS_BEGIN_DECLS

mps_boolean
mps_secular_ga_regenerate_coefficients_mp (mps_context * s, cdpe_t * old_b, mpc_t * old_mb);

MPS_END_DECLS

#endif /* MPS_SECULAR_REGENERATION_H_ */

