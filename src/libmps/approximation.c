#include <mps/mps.h>

mps_approximation *
mps_approximation_new (mps_status * s)
{
  mps_approximation * appr = mps_new (mps_approximation);
  mpc_init2 (appr->mvalue, s->mpwp);
  appr->again = true;
  appr->approximated = false;
  return appr;
}

void
mps_approximation_free (mps_status * s, mps_approximation * appr)
{
  mpc_clear (appr->mvalue);
  free (appr);
}
