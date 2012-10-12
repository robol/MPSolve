/*
 * This file is part of MPSolve 3.0
 *
 * Copyright (C) 2001-2012, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors: 
 *   Dario Andrea Bini <bini@dm.unipi.it>
 *   Giuseppe Fiorentino <fiorent@dm.unipi.it>
 *   Leonardo Robol <robol@mail.dm.unipi.it>
 */


#include <mps/mps.h>

/**
 * @brief Compute the floating point inclusion radius according to the 
 * polynomial representation.
 *
 * @param s The <code>mps_context</code> of the computation.
 * @param fradii The array of double where the radii will be stored.
 */
void
mps_fradii (mps_context * s, double * fradii)
{
  int i;

  if (MPS_INPUT_CONFIG_IS_SECULAR (s->input_config))
    mps_secular_fradii (s, fradii);
  else if (MPS_INPUT_CONFIG_IS_MONOMIAL (s->input_config) && 
	   !MPS_INPUT_CONFIG_IS_USER (s->input_config))
    mps_monomial_fradii (s, fradii);
  else {
	  for (i = 0; i < s->n; i++)
	    fradii[i] = s->root[i]->frad;
  }
}

/**
 * @brief Compute the DPE inclusion radius according to the 
 * polynomial representation.
 *
 * @param s The <code>mps_context</code> of the computation.
 * @param dradii The array of DPE where the radii will be stored.
 */
void
mps_dradii (mps_context * s, rdpe_t * dradii)
{
  int i;

  if (MPS_INPUT_CONFIG_IS_SECULAR (s->input_config))
    mps_secular_dradii (s, dradii);
  else if (MPS_INPUT_CONFIG_IS_MONOMIAL (s->input_config) && 
	   !MPS_INPUT_CONFIG_IS_USER (s->input_config))
    mps_monomial_dradii (s, dradii);
  else {
    for (i = 0; i < s->n; i++)
      rdpe_set (dradii[i], s->root[i]->drad);
  }
}

/**
 * @brief Compute the Multiprecision inclusion radius according to the 
 * polynomial representation.
 *
 * @param s The <code>mps_context</code> of the computation.
 * @param dradii The array of DPE where the radii will be stored.
 */
void
mps_mradii (mps_context * s, rdpe_t * dradii)
{
  int i;

  if (MPS_INPUT_CONFIG_IS_SECULAR (s->input_config))
    mps_secular_mradii (s, dradii);
  else if (MPS_INPUT_CONFIG_IS_MONOMIAL (s->input_config) && 
	   !MPS_INPUT_CONFIG_IS_USER (s->input_config))
    mps_monomial_mradii (s, dradii);
  else {
    for (i = 0; i < s->n; i++)
      rdpe_set (dradii[i], s->root[i]->drad);
  }
}
