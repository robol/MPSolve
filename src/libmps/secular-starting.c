/*
 * secular-starting.c
 *
 *  Created on: 15/giu/2011
 *      Author: leonardo
 */

#include <mps/debug.h>
#include <mps/secular.h>
#include <mps/core.h>
#include <math.h>

#define pi2 6.283184

void
mps_secular_fstart (mps_status * s, int n, mps_cluster_item * cluster_item, double clust_rad,
                    double g, rdpe_t eps)
{
  MPS_DEBUG_THIS_CALL;

  mps_cluster * cluster = NULL;
  mps_root * root;
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

      if (s->status[l][0] != 'a' && s->status[l][0] != 'i' && s->status[l][0] != 'o')
	{
	  cplx_set_d (s->froot[l], cos (i * th + sigma),
		      sin (i * th + sigma));
	  cplx_mul_eq_d (s->froot[l],
			 cplx_mod (s->secular_equation->bfpc[l]) *
			 DBL_EPSILON * 4.0);
	  s->frad[l] += cplx_mod (s->froot[l]);
	  cplx_add_eq (s->froot[l], sec->bfpc[l]);
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


  /* Compute the inclusion radii with Gerschgorin so we can compute
   * clusterizations for the roots. */
  if (!MPS_INPUT_CONFIG_IS_USER (s->input_config))
    mps_monomial_fradii (s);
  mps_fcluster (s, 2.0 * s->n);
  mps_fmodify (s, false);
}

void
mps_secular_dstart (mps_status * s, int n, mps_cluster_item * cluster_item, rdpe_t clust_rad,
                    rdpe_t g, rdpe_t eps)
{
  MPS_DEBUG_THIS_CALL;

  mps_cluster * cluster = cluster_item->cluster;
  int i, l;
  double th = pi2 / n;
  double sigma;
  mps_secular_equation *sec = (mps_secular_equation *) s->secular_equation;
  rdpe_t rtmp;

  cdpe_t ceps;
  cdpe_set (ceps, cdpe_zero);

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

  mps_root * root = cluster->first;
  for (i = 0; i < n; i++)
    {
      l = root->k;
      if (s->status[l][0] != 'a' && s->status[l][0] != 'i' && s->status[l][0] != 'o')
	{
	  cdpe_mod (cdpe_Re (ceps), s->secular_equation->bdpc[l]);
	  rdpe_mul_eq_d (cdpe_Re (ceps), 4 * DBL_EPSILON);
	  cdpe_set_d (s->droot[l], cos (i * th + sigma),
		      sin (i * th + sigma));
	  cdpe_mul_eq (s->droot[l], ceps);

	  if (!rdpe_eq (s->drad[l], RDPE_MAX))
	    {
	      cdpe_mod (rtmp, s->droot[l]);
	      rdpe_add_eq (s->drad[l], rtmp);
	    }

	  cdpe_add_eq (s->droot[l], sec->bdpc[l]);
	  
	  if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)
	    {
	      MPS_DEBUG_CDPE (s, s->droot[l], "s->droot[%d]", l);
	    }
	}

      /* Just an experiment to see if the new method in secular-newton
       * is working */
      /* cdpe_set (s->droot[l +i], sec->bdpc[l + i]); */

      root = root->next;
    }

  mps_dcluster (s, 2.0 * s->n);
  mps_dmodify (s, false);
}

void
mps_secular_mstart (mps_status * s, int n, mps_cluster_item * cluster_item, rdpe_t clust_rad,
                    rdpe_t g, rdpe_t eps)
{
  MPS_DEBUG_THIS_CALL;

  mps_cluster * cluster = NULL;
  mps_root * root;
  int i, l;
  double th = pi2 / n;
  double sigma;
  mpc_t epsilon;
  mps_secular_equation *sec = (mps_secular_equation *) s->secular_equation;

  rdpe_t r_eps;
  cdpe_t ctmp;
  rdpe_t rtmp;
  rdpe_set (r_eps, rdpe_zero);

  mpc_init2 (epsilon, mpc_get_prec (s->mroot[0]));
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

      /* Compute the right epsilon */
      mpc_get_cdpe (ctmp, s->mroot[l]);
      cdpe_mod (r_eps, ctmp);
      rdpe_mul_eq (r_eps, s->mp_epsilon);
      rdpe_mul_eq_d (r_eps, 4.0);
      mpf_set_rdpe (mpc_Re (epsilon), r_eps);

      mpc_set_d (s->mroot[l], cos (i * th + sigma), sin (i * th + sigma));
      mpc_mul_eq (s->mroot[l], epsilon);

      mpc_get_cdpe (ctmp, s->mroot[l]);
      cdpe_mod (rtmp, ctmp);
      rdpe_add_eq (s->drad[l], rtmp);

      mpc_add_eq (s->mroot[l], sec->bmpc[l]);

      if (cluster_item)
	{
	  root = root->next;
	  if (!root)
	    break;
	}
    }

  mps_mcluster (s, 2.0 * s->n);
  mps_mmodify (s, false);
  
  mpc_clear (epsilon);
}
