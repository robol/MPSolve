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
mps_secular_fstart (mps_status * s, int n, int i_clust, double clust_rad,
                    double g, rdpe_t eps)
{
  MPS_DEBUG_THIS_CALL;

  int i, l = s->punt[i_clust];
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
      if (i_clust == 0)
        sigma = s->last_sigma = 0.1;
      else
        sigma = mps_maximize_distance (s, s->last_sigma, i_clust, n);
    }
    
  /* Determine a suitable epsilon to move the roots */
  double a_eps = 0;
  double tmp;
  for (i = 0; i < s->n; i++)
  {
      tmp = cplx_mod (s->secular_equation->afpc[i]);
      if (tmp > a_eps)
        a_eps = tmp;
  }
  a_eps *= s->n;
  a_eps *= DBL_EPSILON;

  /* The roots are set as the b_i plus a small correction that is the
   * disposition on the unit cicle scaled to DBL_EPSILON */
  for (i = 0; i < s->n; i++)
    {
      cplx_set_d (s->froot[l + i], cos (i * th + sigma),
                  sin (i * th + sigma));
      cplx_mul_eq_d (s->froot[l + i], a_eps);
      cplx_add_eq (s->froot[l + i], sec->bfpc[l + i]);
    }
}

void
mps_secular_dstart (mps_status * s, int n, int i_clust, rdpe_t clust_rad,
                    rdpe_t g, rdpe_t eps)
{
  MPS_DEBUG_THIS_CALL;

  int i, l = s->punt[i_clust];
  double th = pi2 / n;
  double sigma;
  mps_secular_equation *sec = (mps_secular_equation *) s->secular_equation;
  cdpe_t ceps;
  cdpe_set (ceps, cdpe_zero);
  
  rdpe_set (cdpe_Re (ceps), eps);
  if (cdpe_eq_zero (ceps) || !MPS_STRUCTURE_IS_FP (s->secular_equation->input_structure))
  {
      rdpe_set_d (cdpe_Re (ceps), DBL_EPSILON);
  }

  /* Get best sigma possible */
  if (s->random_seed)
    sigma = drand ();
  else
    {
      /* If this is the first cluster select sigma = 0. In the other
       * case try to maximize starting points distance. */
      if (i_clust == 0)
        {
          sigma = s->last_sigma = 0.1;
        }
      else
        {
          sigma = mps_maximize_distance (s, s->last_sigma, i_clust, n);
        }
    }

  for (i = 0; i < s->n; i++)
    {
      cdpe_set_d (s->droot[l + i], cos (i * th + sigma),
                  sin (i * th + sigma));
      cdpe_mul_eq (s->droot[l + i], ceps);
      cdpe_add_eq (s->droot[l + i], sec->bdpc[l + i]);
      rdpe_add_eq (s->drad[i], eps);
    }
}

void
mps_secular_mstart (mps_status * s, int n, int i_clust, rdpe_t clust_rad,
                    rdpe_t g, rdpe_t eps)
{
  MPS_DEBUG_THIS_CALL;

  int i, l = s->punt[i_clust];
  double th = pi2 / n;
  double sigma;
  mpc_t epsilon;
  mps_secular_equation *sec = (mps_secular_equation *) s->secular_equation;

  rdpe_t r_eps;
  rdpe_t rtmp;
  rdpe_set (r_eps, rdpe_zero);
  for (i = 0; i < s->n; i++)
  {
      cdpe_mod (rtmp, s->droot[i]);
      if (rdpe_gt (rtmp, r_eps))
        rdpe_set (r_eps, rtmp);
  }
  rdpe_mul_eq_d (r_eps, s->n);
  rdpe_mul_eq (r_eps, s->mp_epsilon);

  mpc_init2 (epsilon, s->mpwp);
  mpc_set_ui (epsilon, 0, 0);

  mpf_set_rdpe (mpc_Re (epsilon), r_eps);
  if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)
    {
      MPS_DEBUG_MPC (s, 100, epsilon, "epsilon");
    }

  /* Get best sigma possible */
  if (s->random_seed)
    sigma = drand ();
  else
    {
      /* If this is the first cluster select sigma = 0. In the other
       * case try to maximize starting points distance. */
      if (i_clust == 0)
        {
          sigma = s->last_sigma = 0.1;
        }
      else
        {
          sigma = mps_maximize_distance (s, s->last_sigma, i_clust, n);
        }
    }

  for (i = 0; i < s->n; i++)
    {
      mpc_set_d (s->mroot[l + i], cos (i * th + sigma), sin (i * th + sigma));
      mpc_mul_eq (s->mroot[l + i], epsilon);
      mpc_add_eq (s->mroot[l + i], sec->bmpc[l + i]);
      rdpe_add_eq (s->drad[i], r_eps);
    }

  mpc_clear (epsilon);
}
