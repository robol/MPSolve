/*
 * This file is part of MPSolve 3.0
 *
 * Copyright (C) 2001-2013, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors: 
 *   Dario Andrea Bini <bini@dm.unipi.it>
 *   Giuseppe Fiorentino <fiorent@dm.unipi.it>
 *   Leonardo Robol <robol@mail.dm.unipi.it>
 */


#include <mps/mps.h>
#include <math.h>

#define pi2 6.283184

void
mps_secular_fstart (mps_context * s, mps_secular_equation * sec)
{
  MPS_DEBUG_THIS_CALL;

  int i;
  int n = MPS_POLYNOMIAL (sec)->degree;

  for (i = 0; i < n; i++)
    {
      if (!MPS_ROOT_STATUS_IS_COMPUTED (s, i))
        {
          cplx_set (s->root[i]->fvalue, sec->bfpc[i]);
          if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)
            MPS_DEBUG_CPLX (s, s->root[i]->fvalue, "s->froot[%d]", i);
        }
    }
}

void
mps_secular_dstart (mps_context * s, mps_secular_equation * sec)
{
  MPS_DEBUG_THIS_CALL;

  int l;
  for (l = 0; l < MPS_POLYNOMIAL (sec)->degree; l++)
    {
      if (!MPS_ROOT_STATUS_IS_COMPUTED (s, l))
        {
          cdpe_set (s->root[l]->dvalue, sec->bdpc[l]);
          
          if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)
            {
              MPS_DEBUG_CDPE (s, s->root[l]->dvalue, "s->droot[%d]", l);
            }
        }
    }
}

void
mps_secular_mstart (mps_context * s, mps_secular_equation * sec)
{
  MPS_DEBUG_THIS_CALL;

  int l;
  for (l = 0; l < MPS_POLYNOMIAL (sec)->degree; l++)
    {
      if (!MPS_ROOT_STATUS_IS_COMPUTED (s, l))
        {
          mpc_set (s->root[l]->mvalue, sec->bmpc[l]);
        }         
    }
}
