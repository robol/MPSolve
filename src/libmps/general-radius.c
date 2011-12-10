#include <mps/core.h>

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
  switch (s->algorithm)
    {
    case MPS_ALGORITHM_STANDARD_MPSOLVE:
      if (!MPS_INPUT_CONFIG_IS_USER (s->input_config))
	mps_monomial_fradii (s, fradii);
      else
	{
	  for (i = 0; i < s->n; i++)
	    fradii[i] = s->frad[i];
	}
      break;
    case MPS_ALGORITHM_SECULAR_MPSOLVE:
    case MPS_ALGORITHM_SECULAR_GA:
      mps_secular_fradii (s, fradii);
      break;
    default:
      mps_error (s, 1, "Unsupported algorithm selected, aborting.");
    }
}

/**
 * @brief Compute the DPE inclusion radius according to the 
 * polynomial representation.
 *
 * @param s The <code>mps_status</code> of the computation.
 * @param s The array of DPE where the radii will be stored.
 */
void
mps_dradii (mps_status * s, rdpe_t * dradii)
{
  int i;
  switch (s->algorithm)
    {
    case MPS_ALGORITHM_STANDARD_MPSOLVE:
      MPS_DEBUG (s, "s->input_config->density = %d", s->input_config->density);
      if (MPS_INPUT_CONFIG_IS_USER (s->input_config))
	{
	  for (i = 0; i < s->n; i++)
	    rdpe_set (dradii[i], s->drad[i]);
	}
      else
	mps_monomial_dradii (s, dradii);
      break;
    case MPS_ALGORITHM_SECULAR_MPSOLVE:
    case MPS_ALGORITHM_SECULAR_GA:
      mps_secular_dradii (s, dradii);
      break;
    default:
      mps_error (s, 1, "Unsupported algorithm selected, aborting.");
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
  switch (s->algorithm)
    {
    case MPS_ALGORITHM_STANDARD_MPSOLVE:
      mps_monomial_mradii (s, dradii);
      break;
    case MPS_ALGORITHM_SECULAR_MPSOLVE:
    case MPS_ALGORITHM_SECULAR_GA:
      mps_secular_mradii (s, dradii);
      break;
    default:
      mps_error (s, 1, "Unsupported algorithm selected, aborting.");
    }
}
