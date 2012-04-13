/*
 * secular-starting.c
 *
 *  Created on: 15/giu/2011
 *      Author: leonardo
 */

#include <mps/mps.h>
#include <math.h>

#define pi2 6.283184

void
mps_secular_fstart (mps_status * s, int n, mps_cluster_item * cluster_item, double clust_rad,
                    double g, rdpe_t eps)
{
  MPS_DEBUG_THIS_CALL;

  mps_cluster * cluster = NULL;
  mps_root * root = NULL;
  int i, l;
  mps_secular_equation *sec = (mps_secular_equation *) s->secular_equation;

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
	  cplx_set (s->froot[l], sec->bfpc[l]);
	  if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)
	    MPS_DEBUG_CPLX (s, s->froot[l], "s->froot[%d]", l);
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
mps_secular_dstart (mps_status * s, int n, mps_cluster_item * cluster_item, rdpe_t clust_rad,
                    rdpe_t g, rdpe_t eps)
{
  MPS_DEBUG_THIS_CALL;

  mps_cluster * cluster = NULL;
  mps_root * root = NULL;
  int i, l;
  mps_secular_equation *sec = (mps_secular_equation *) s->secular_equation;

  cdpe_t ceps;
  cdpe_set (ceps, cdpe_zero);

  if (cluster_item)
    cluster = cluster_item->cluster;

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
	  cdpe_set (s->droot[l], sec->bdpc[l]);
	  
	  if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)
	    {
	      MPS_DEBUG_CDPE (s, s->droot[l], "s->droot[%d]", l);
	    }
	}

      if (cluster_item)
	root = root->next;
    }
}

void
mps_secular_mstart (mps_status * s, int n, mps_cluster_item * cluster_item, rdpe_t clust_rad,
                    rdpe_t g, rdpe_t eps)
{
  MPS_DEBUG_THIS_CALL;

  mps_cluster * cluster = NULL;
  mps_root * root = NULL;
  int i, l;
  mpc_t epsilon;
  mps_secular_equation *sec = (mps_secular_equation *) s->secular_equation;

  mpc_init2 (epsilon, mpc_get_prec (s->mroot[0]));
  mpc_set_ui (epsilon, 0, 0);

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
	  mpc_set (s->mroot[l], sec->bmpc[l]);
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
