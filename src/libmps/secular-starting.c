/*
 * secular-starting.c
 *
 *  Created on: 15/giu/2011
 *      Author: leonardo
 */

#include <mps/debug.h>
#include <mps/secular.h>
#include <mps/core.h>

#define pi2 6.283184

void
mps_secular_fstart(mps_status* s, int n, int i_clust, double clust_rad,
    double g, rdpe_t eps)
{
  MPS_DEBUG_THIS_CALL

  int i, l = s->punt[i_clust];
  double th = pi2 / n;
  double sigma;
  mps_secular_equation* sec = (mps_secular_equation*) s->secular_equation;

  /* Get best sigma possible */
  if (s->random_seed)
    sigma = drand();
  else
    {
      /* If this is the first cluster select sigma = 0. In the other
       * case try to maximize starting points distance. */
      if (i_clust == 0)
        sigma = s->last_sigma = 0.1;
      else
        sigma = mps_maximize_distance(s, s->last_sigma, i_clust, n);
    }

  /* The roots are set as the b_i plus a small correction that is the
   * disposition on the unit cicle scaled to DBL_EPSILON */
  for (i = 0; i < s->n; i++)
    {
      cplx_set_d(s->froot[l + i], cos(i * th + sigma), sin(i * th + sigma));
      cplx_mul_eq_d(s->froot[l + i], 5 * DBL_EPSILON);
      cplx_add_eq(s->froot[l + i], sec->bfpc[l + i]);
    }
}

void
mps_secular_dstart(mps_status* s, int n, int i_clust, rdpe_t clust_rad,
    rdpe_t g, rdpe_t eps)
{
  MPS_DEBUG_THIS_CALL

  int i, l = s->punt[i_clust];
  double th = pi2 / n;
  double sigma;
  mps_secular_equation* sec = (mps_secular_equation*) s->secular_equation;

  /* Get best sigma possible */
  if (s->random_seed)
    sigma = drand();
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
          sigma = mps_maximize_distance(s, s->last_sigma, i_clust, n);
        }
    }

  for (i = 0; i < s->n; i++)
    {
      cdpe_set_d(s->droot[l + i], cos(i * th + sigma), sin(i * th + sigma));
      cdpe_mul_eq_d(s->droot[l + i], 5 * DBL_EPSILON);
      cdpe_add_eq(s->droot[l + i], sec->bdpc[l + i]);
    }
}

void
mps_secular_mstart(mps_status* s, int n, int i_clust, rdpe_t clust_rad,
    rdpe_t g, rdpe_t eps)
{
  MPS_DEBUG_THIS_CALL

  int i, l = s->punt[i_clust];
  double th = pi2 / n;
  double sigma;
  mpc_t epsilon;
  mps_secular_equation* sec = (mps_secular_equation*) s->secular_equation;

  mpc_init2(epsilon, s->mpwp);
  mpc_set_ui(epsilon, 0, 0);
  mpf_set_rdpe(mpc_Re(epsilon), s->mp_epsilon);
  MPS_DEBUG_MPC(s, 100, epsilon, "epsilon");

  /* Get best sigma possible */
  if (s->random_seed)
    sigma = drand();
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
          sigma = mps_maximize_distance(s, s->last_sigma, i_clust, n);
        }
    }

  for (i = 0; i < s->n; i++)
    {
      mpc_set_d(s->mroot[l + i], cos(i * th + sigma), sin(i * th + sigma));
      mpc_mul_eq(s->mroot[l + i], epsilon);
      mpc_add_eq(s->mroot[l + i], sec->bmpc[l + i]);
      rdpe_add_eq(s->drad[i], s->mp_epsilon);
    }

  mpc_clear(epsilon);
}
