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

mps_approximation *
mps_approximation_new (mps_context * s)
{
  mps_approximation * appr = mps_new (mps_approximation);
  mpc_init2 (appr->mvalue, s->mpwp);
  appr->again = true;
  appr->approximated = false;
  return appr;
}

void
mps_approximation_free (mps_context * s, mps_approximation * appr)
{
  mpc_clear (appr->mvalue);
  free (appr);
}
