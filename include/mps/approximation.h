/*
 * This file is part of MPSolve 3.1.5
 *
 * Copyright (C) 2001-2014, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Leonardo Robol <robol@mail.dm.unipi.it>
 */

#ifndef MPS_APPROXIMATION_H_
#define MPS_APPROXIMATION_H_

#ifdef __cplusplus
extern "C" {
#endif

#ifdef _MPS_PRIVATE

struct mps_approximation {
  cplx_t fvalue;
  cdpe_t dvalue;
  mpc_t mvalue;

  double frad;
  rdpe_t drad;

  mps_boolean approximated;
  mps_boolean again;
  long int wp;

  mps_root_status status;

  /**
   * @brief Attributes that have been set on
   * the roots.
   */
  mps_root_attrs attrs;

  /**
   * @brief Inclusion status of the root
   * in the target set specified in the field
   * <code>input_config->search_set</code>.
   */
  mps_root_inclusion inclusion;
};

#endif

/* Public accessor functions */
void mps_approximation_get_fvalue (mps_context * ctx, mps_approximation * approximation, cplx_t output);
void mps_approximation_get_dvalue (mps_context * ctx, mps_approximation * approximation, cdpe_t output);
void mps_approximation_get_mvalue (mps_context * ctx, mps_approximation * approximation, mpc_t output);
double mps_approximation_get_frad (mps_context * ctx, mps_approximation * approximation);
void mps_approximation_get_drad (mps_context * ctx, mps_approximation * approximation, rdpe_t output);
mps_root_status mps_approximation_get_status (mps_context * ctx, mps_approximation * approximation);
mps_root_attrs mps_approximation_get_attrs (mps_context * ctx, mps_approximation * approximation);
mps_root_inclusion mps_approximaiton_get_inclusion (mps_context * ctx, mps_approximation * approximation);

/* Public setters functions */
void mps_approximation_set_fvalue (mps_context * ctx, mps_approximation * approximation, cplx_t value);
void mps_approximation_set_dvalue (mps_context * ctx, mps_approximation * approximation, cdpe_t value);
void mps_approximation_set_mvalue (mps_context * ctx, mps_approximation * approximation, mpc_t value);
void mps_approximation_set_frad (mps_context * ctx, mps_approximation * approximation, double frad);
void mps_approximation_set_drad (mps_context * ctx, mps_approximation * approximation, rdpe_t drad);
void mps_approximation_set_status (mps_context * ctx, mps_approximation * approximation, mps_root_status status);
void mps_approximation_set_attrs (mps_context * ctx, mps_approximation * approximation, mps_root_attrs attrs);
void mps_approximation_set_inclusion (mps_context * ctx, mps_approximation * approximation, mps_root_inclusion inclusion);



#ifdef __cplusplus
}
#endif

#endif /* #ifndef MPS_ROOT_H */
