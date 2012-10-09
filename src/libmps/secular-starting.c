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
#include <math.h>

#define pi2 6.283184

void
mps_secular_fstart (mps_context * s, int n, mps_cluster_item * cluster_item, double clust_rad,
                    double g, rdpe_t eps)
{
  MPS_DEBUG_THIS_CALL;

  mps_cluster * cluster = NULL;
  mps_root * root = NULL;
  int i, l;
  double th = pi2 / n;
  double sigma;
  mps_secular_equation *sec = (mps_secular_equation *) s->secular_equation;

  /* Get best sigma possible */
  if (s->random_seed)
    sigma = drand ();
  else
    {
      /* If this is the first cluster select sigma = 0. In the other
       * case try to maximize starting points distance. */
      if (cluster_item == NULL || cluster == s->clusterization->first->cluster)
        sigma = s->last_sigma = 0.1;
      else
        sigma = mps_maximize_distance (s, s->last_sigma, cluster_item, n);
    }

  /* The roots are set as the b_i plus a small correction that is the
   * disposition on the unit cicle scaled to DBL_EPSILON */
  if (cluster_item)
    {
      cluster = cluster_item->cluster;
      root = cluster->first;
    }
  for (i = 0; i < n; i++)
    {
      if (cluster_item)
	l = root->k;
      else
	l = i;

      if (!MPS_ROOT_STATUS_IS_COMPUTED (s, i))
	{
	  if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)
	    MPS_DEBUG_CPLX (s, s->root[l]->fvalue, "s->froot[%d] (before start)", l);

	  cplx_set_d (s->root[l]->fvalue, cos (i * th + sigma),
		      sin (i * th + sigma));
	  cplx_mul_eq_d (s->root[l]->fvalue,
			 cplx_mod (s->secular_equation->bfpc[l]) *
			 DBL_EPSILON * 4.0);
	  s->root[l]->frad += cplx_mod (s->root[l]->fvalue);
	  cplx_add_eq (s->root[l]->fvalue, sec->bfpc[l]);
	  if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)
	    MPS_DEBUG_CPLX (s, s->root[l]->fvalue, "s->froot[%d]", l);
	}

      if (cluster_item)
	{
	  root = root->next;
	  if (!root)
	    return;
	}
    }
}

void
mps_secular_dstart (mps_context * s, int n, mps_cluster_item * cluster_item, rdpe_t clust_rad,
                    rdpe_t g, rdpe_t eps)
{
  MPS_DEBUG_THIS_CALL;

  mps_cluster * cluster = NULL;
  mps_root * root = NULL;
  int i, l;
  double th = pi2 / n;
  double sigma;
  mps_secular_equation *sec = (mps_secular_equation *) s->secular_equation;
  rdpe_t rtmp;

  cdpe_t ceps;
  cdpe_set (ceps, cdpe_zero);

  if (cluster_item)
    cluster = cluster_item->cluster;

  /* Get best sigma possible */
  if (s->random_seed)
    sigma = drand ();
  else
    {
      /* If this is the first cluster select sigma = 0. In the other
       * case try to maximize starting points distance. */
      if (cluster == s->clusterization->first->cluster)
        {
          sigma = s->last_sigma = 0.1;
        }
      else
        {
          sigma = mps_maximize_distance (s, s->last_sigma, cluster_item, n);
        }
    }

  if (cluster_item)
    root = cluster->first;
  for (i = 0; i < n; i++)
    {
      if (cluster_item)
	l = root->k;
      else
	l = i;

      if (!MPS_ROOT_STATUS_IS_COMPUTED (s, l))
	{
	  cdpe_mod (cdpe_Re (ceps), s->secular_equation->bdpc[l]);
	  rdpe_mul_eq_d (cdpe_Re (ceps), 4 * DBL_EPSILON);
	  cdpe_set_d (s->root[l]->dvalue, cos (l * th + sigma),
		      sin (l * th + sigma));
	  cdpe_mul_eq (s->root[l]->dvalue, ceps);

	  if (!rdpe_eq (s->root[l]->drad, RDPE_MAX))
	    {
	      cdpe_mod (rtmp, s->root[l]->dvalue);
	      rdpe_add_eq (s->root[l]->drad, rtmp);
	    }

	  cdpe_add_eq (s->root[l]->dvalue, sec->bdpc[l]);
	  
	  if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)
	    {
	      MPS_DEBUG_CDPE (s, s->root[l]->dvalue, "s->droot[%d]", l);
	    }
	}

      /* Just an experiment to see if the new method in secular-newton
       * is working */
      /* cdpe_set (s->droot[l +i], sec->bdpc[l + i]); */

      if (cluster_item)
	root = root->next;
    }
}

void
mps_secular_mstart (mps_context * s, int n, mps_cluster_item * cluster_item, rdpe_t clust_rad,
                    rdpe_t g, rdpe_t eps)
{
  MPS_DEBUG_THIS_CALL;

  mps_cluster * cluster = NULL;
  mps_root * root = NULL;
  int i, l;
  double th = pi2 / n;
  double sigma;
  mpc_t epsilon;
  mps_secular_equation *sec = (mps_secular_equation *) s->secular_equation;

  rdpe_t r_eps;
  cdpe_t ctmp;
  rdpe_t rtmp;
  rdpe_set (r_eps, rdpe_zero);

  mpc_init2 (epsilon, mpc_get_prec (s->root[0]->mvalue));
  mpc_set_ui (epsilon, 0, 0);

  /* Get best sigma possible */
  if (s->random_seed)
    sigma = drand ();
  else
    {
      /* If this is the first cluster select sigma = 0. In the other
       * case try to maximize starting points distance. */
      if (cluster == s->clusterization->first->cluster)
        {
          sigma = s->last_sigma = 0.1;
        }
      else
        {
          sigma = mps_maximize_distance (s, s->last_sigma, cluster_item, n);
        }
    }

  if (cluster_item)
    {
      cluster = cluster_item->cluster;
      root = cluster->first;
    }
  for (i = 0; i < n; i++)
    {
      if (cluster_item)
	l = root->k;
      else
	l = i;

      if (!MPS_ROOT_STATUS_IS_COMPUTED (s, l))
	{
	  /* Compute the right epsilon */
	  mpc_get_cdpe (ctmp, s->secular_equation->bmpc[l]);
	  cdpe_mod (r_eps, ctmp);
	  rdpe_mul_eq (r_eps, s->mp_epsilon);
	  rdpe_mul_eq_d (r_eps, 4.0);
	  mpf_set_rdpe (mpc_Re (epsilon), r_eps);
	  
	  mpc_set_d (s->root[l]->mvalue, cos (i * th + sigma), sin (i * th + sigma));
	  mpc_mul_eq (s->root[l]->mvalue, epsilon);
	  
	  mpc_get_cdpe (ctmp, s->root[l]->mvalue);
	  cdpe_mod (rtmp, ctmp);
	  rdpe_add_eq (s->root[l]->drad, rtmp);

	  mpc_add_eq (s->root[l]->mvalue, sec->bmpc[l]);
	}
	  
      if (cluster_item)
	{
	  root = root->next;
	  if (!root)
	    break;
	}
    }

  mpc_clear (epsilon);
}
