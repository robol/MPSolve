#include <mps/core.h>

/**
 * @brief Compute the floating point inclusion radius according to the 
 * polynomial representation.
 *
 * @param s The <code>mps_status</code> of the computation.
 */
void
mps_fradii (mps_status * s)
{
  if (MPS_INPUT_CONFIG_IS_MONOMIAL (s->input_config))
    mps_monomial_fradii (s);
  else if (MPS_INPUT_CONFIG_IS_SECULAR (s->input_config))
    mps_secular_fradii (s);
}
