/************************************************************
 **                                                        **
 **             __  __ ___  ___      _                     **
 **            |  \/  | _ \/ __| ___| |_ _____             **
 **            | |\/| |  _/\__ \/ _ \ \ V / -_)            **
 **            |_|  |_|_|  |___/\___/_|\_/\___|            **
 **                                                        **
 **       Multiprecision Polynomial Solver (MPSolve)       **
 **                 Version 2.9, April 2011                **
 **                                                        **
 **                      Written by                        **
 **                                                        **
 **     Dario Andrea Bini       <bini@dm.unipi.it>         **
 **     Giuseppe Fiorentino     <fiorent@dm.unipi.it>      **
 **     Leonardo Robol          <robol@mail.dm.unipi.it>   **
 **                                                        **
 **           (C) 2011, Dipartimento di Matematica         **
 ***********************************************************/

#include <mps/mps.h>

/**
 * @brief Compute the floating point inclusion radius according to the 
 * polynomial representation.
 *
 * @param s The <code>mps_status</code> of the computation.
 * @param fradii The array of double where the radii will be stored.
 */
void
mps_fradii (mps_status * s, double * fradii)
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
 * @param s The <code>mps_status</code> of the computation.
 * @param dradii The array of DPE where the radii will be stored.
 */
void
mps_dradii (mps_status * s, rdpe_t * dradii)
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
 * @param s The <code>mps_status</code> of the computation.
 * @param dradii The array of DPE where the radii will be stored.
 */
void
mps_mradii (mps_status * s, rdpe_t * dradii)
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
