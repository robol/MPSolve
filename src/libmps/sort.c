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

/********************************************************
*      SUBROUTINE FCMP                                  *
*********************************************************/
int
mps_fcmp (const void *a, const void *b)
{
  double dif;

  dif = cplx_Re (cplx_Addr (a)) - cplx_Re (cplx_Addr (b));
  if (dif > 0)
    return 1;
  else if (dif == 0.0)
    return 0;
  return -1;
}

/*********************************************************
*      SUBROUTINE FSORT                                  *
*********************************************************/
void
mps_fsort (mps_context * s)
{
  int i;
  cplx_t *real_parts = cplx_valloc (s->n);

  for (i = 0; i < s->n; i++)
    {
      cplx_set_d (real_parts[i],
                  cplx_Re (s->root[i]->fvalue),
                  i);
    }

  qsort (real_parts, s->n, sizeof (cplx_t), mps_fcmp);

  for (i = 0; i < s->n; i++)
    s->order[i] = (int) cplx_Im (real_parts[i]);

  cplx_vfree (real_parts);
}

/*********************************************************
*      SUBROUTINE DCMP                                  *
*********************************************************/
int
mps_dcmp (const void *a, const void *b)
{
  return rdpe_cmp (cdpe_Re (cdpe_Addr (a)), cdpe_Re (cdpe_Addr (b)));
}

/*********************************************************
*      SUBROUTINE DSORT                                  *
*********************************************************/
void
mps_dsort (mps_context * s)
{
  cdpe_t * real_parts = cdpe_valloc (s->n);
  int i;

  for (i = 0; i < s->n; i++)
    {
      rdpe_set (cdpe_Re (real_parts[i]), cdpe_Re (s->root[i]->dvalue));
      rdpe_set_d (cdpe_Im (real_parts[i]), i);
    }

  qsort (real_parts, s->n, sizeof (cdpe_t), mps_dcmp);

  for (i = 0; i < s->n; i++)
    s->order[i] = (int) rdpe_get_d (cdpe_Im (real_parts[i]));

  cdpe_vfree (real_parts);
}

/*********************************************************
*      SUBROUTINE MCMP                                  *
*********************************************************/
int
mps_mcmp (const void *a, const void *b)
{
  return mpf_cmp (mpc_Re (mpc_Addr (a)), mpc_Re (mpc_Addr (b)));
}

/*********************************************************
*      SUBROUTINE MSORT                                  *
*********************************************************/
void
mps_msort (mps_context * s)
{
  int i;
  mpc_t * real_parts = mpc_valloc (s->n);
  mpc_vinit2 (real_parts, s->n, s->mpwp);

  for (i = 0; i < s->n; i++)
    {
      mpf_set (mpc_Re (real_parts[i]), mpc_Re (s->root[i]->mvalue));
      mpf_set_ui (mpc_Im (real_parts[i]), i);
    }

  qsort (real_parts, s->n, sizeof (mpc_t), mps_mcmp);

  for (i = 0; i < s->n; i++)
    s->order[i] = (int) mpf_get_d (mpc_Im (real_parts[i]));

  mpc_vclear (real_parts, s->n);
  mpc_vfree (real_parts);
}
