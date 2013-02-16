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

#include <mps/gmptools.h>
#include <mps/mps.h>


/********************************************************
*      SUBROUTINE INCLUSION                             *
*********************************************************/
mps_boolean
mps_inclusion (mps_context * s)
{
  int i, j, k, oldnclust;
  rdpe_t rad, difr;
  cdpe_t difc;
  mpc_t tmp;
  rdpe_t ap, az, temp, ep, apeps;
  cdpe_t temp1;
  mpc_t p;
  mps_monomial_poly *poly = MPS_MONOMIAL_POLY (s->active_poly);


  /* add inclusion code here */
  if (!s->chkrad || s->lastphase != mp_phase)
    {
      if (s->DOLOG)
        fprintf (s->logstr, "Skipping inclusion disks check.\n");
      return true;
    }

  if (s->DOLOG)
    fprintf (s->logstr, "Checking inclusion disks...\n");

  if (s->DOLOG)
    {
      fprintf (s->logstr, "Old radii\n");
      for (i = 0; i < s->n; i++)
        {
          fprintf (s->logstr, "r(%d)=", i);
          rdpe_outln_str (s->logstr, s->root[i]->drad);
        }
    }

  /* save old radii */
  for (i = 0; i < s->n; i++)
    rdpe_set (s->dap1[i], s->root[i]->drad);

  mpc_init2 (p, s->mpwp);
  rdpe_mul_d (ep, s->mp_epsilon, (double) (s->n * 4));

  mpc_init2 (tmp, s->mpwp);

  for (i = 0; i < s->n; i++)
    {

      /* compute denominator */
      rdpe_set (rad, rdpe_one);
      for (j = 0; j < s->n; j++)
        {
          if (i == j)
            continue;
          mpc_sub (tmp, s->root[j]->mvalue, s->root[i]->mvalue);
          mpc_get_cdpe (difc, tmp);
          cdpe_smod (difr, difc);
          rdpe_mul_eq (rad, difr);
        }
      rdpe_sqrt_eq (rad);
      rdpe_mul_eq (rad, poly->dap[s->n]);

      /* compute numerator */
      if (MPS_DENSITY_IS_SPARSE (s->active_poly->density))
        {                       /* case of sparse polynomial */
          /* compute p(mroot[i]) */
	  mps_polynomial_meval (s, MPS_POLYNOMIAL (poly), s->root[i]->mvalue, p, ap);
	  rdpe_div_eq (ap, s->mp_epsilon);
        }
      else
        {                       /*  dense polynomial */

          /* commpute p(mroot[i]) and p'(mroot[i]) */
          mpc_set (p, poly->mfpc[s->n]);
          for (k = s->n - 1; k > 0; k--)
            {
              mpc_mul (p, p, s->root[i]->mvalue);
              mpc_add (p, p, poly->mfpc[k]);
            }
          mpc_mul (p, p, s->root[i]->mvalue);
          mpc_add (p, p, poly->mfpc[0]);

          /* compute bound to the error */
          rdpe_set (ap, poly->dap[s->n]);
          mpc_get_cdpe (temp1, s->root[i]->mvalue);
          cdpe_mod (az, temp1);
          for (k = s->n - 1; k >= 0; k--)
            {
              rdpe_mul (temp, ap, az);
              rdpe_add (ap, temp, poly->dap[k]);
            }
        }

      /* common part */
      mpc_get_cdpe (difc, p);
      cdpe_mod (difr, difc);
      rdpe_mul (apeps, ap, ep);
      rdpe_add_eq (apeps, difr);
      rdpe_mul_eq_d (apeps, (double) s->n);

      /* compute ratio */
      rdpe_div (s->root[i]->drad, apeps, rad);

      if (s->DOLOG)
        {
          fprintf (s->logstr, "New r(%d)=", i);
          rdpe_outln_str (s->logstr, s->root[i]->drad);
        }
    }

  oldnclust = s->clusterization->n;

  rdpe_t * newton_radii = rdpe_valloc (s->n);
  for (i = 0; i < s->n; i++)
    rdpe_set (newton_radii[i], s->root[i]->drad);

  mps_mcluster (s, newton_radii, 2 * s->n);
  free (newton_radii);

  if (s->clusterization->n >= oldnclust)
    {
      /* choose the smallest radius */
      for (i = 0; i < s->n; i++)
        if (rdpe_lt (s->dap1[i], s->root[i]->drad))
          rdpe_set (s->root[i]->drad, s->dap1[i]);
      /* update(); */
    }
  else
    mps_warn (s, "Some roots might be not approximated");

  mpc_clear (tmp);
  mpc_clear (p);

  return true;
}
