/************************************************************
 **                                                        **
 **             __  __ ___  ___      _                     **
 **            |  \/  | _ \/ __| ___| |_ _____             **
 **            | |\/| |  _/\__ \/ _ \ \ V / -_)            **
 **            |_|  |_|_|  |___/\___/_|\_/\___|            **
 **                                                        **
 **       Multiprecision Polynomial Solver (MPSolve)       **
 **                 Version 2.9, April 2011                **
 **                                                        **
 **                      Written by                        **
 **                                                        **
 **     Dario Andrea Bini       <bini@dm.unipi.it>         **
 **     Giuseppe Fiorentino     <fiorent@dm.unipi.it>      **
 **     Leonardo Robol          <robol@mail.dm.unipi.it>   **
 **                                                        **
 **           (C) 2011, Dipartimento di Matematica         **
 ***********************************************************/

#include <mps/core.h>
#include <mps/cluster.h>
#include <pthread.h>

/**
 * @brief Compute Aberth correction for j-th root, without
 * selective correction.
 */
void
mps_faberth (mps_status * s, int j, cplx_t abcorr)
{
  int i;
  cplx_t z;

  cplx_set (abcorr, cplx_zero);
  for (i = 0; i < s->n; i++)
    {
      if (i == j)
        continue;
      cplx_sub (z, s->froot[j], s->froot[i]);
      cplx_inv_eq (z);
      cplx_add_eq (abcorr, z);
    }
}

/**
 * @brief Compute Aberth correction for j-th root, without
 * selective correction.
 */
void
mps_daberth (mps_status * s, int j, cdpe_t abcorr)
{
  int i;
  cdpe_t z;

  cdpe_set (abcorr, cdpe_zero);
  for (i = 0; i < s->n; i++)
    {
      if (i == j)
        continue;
      cdpe_sub (z, s->droot[j], s->droot[i]);
      cdpe_inv_eq (z);
      cdpe_add_eq (abcorr, z);
    }
}

/**
 * @brief Compute Aberth correction for j-th root, without
 * selective correction.
 */
void
mps_maberth (mps_status * s, int j, mpc_t abcorr)
{
  int i;
  cdpe_t z, temp;
  mpc_t diff;

  mpc_init2 (diff, s->mpwp);

  cdpe_set (temp, cdpe_zero);
  for (i = 0; i < s->n; i++)
    {
      if (i == j)
        continue;
      mpc_sub (diff, s->mroot[j], s->mroot[i]);
      mpc_get_cdpe (z, diff);
      cdpe_inv_eq (z);
      cdpe_add_eq (temp, z);
    }
  mpc_set_cdpe (abcorr, temp);

  mpc_clear (diff);
}

/**
 * @brief Compute Aberth correction for the j-th root,
 * but only with other roots of the <code>jc</code>-th
 * cluster.
 */
void
mps_faberth_s (mps_status * s, int j, mps_cluster * cluster, cplx_t abcorr)
{
  int k;
  cplx_t z;
  mps_root * root;

  cplx_set (abcorr, cplx_zero);
  for (root = cluster->first; root != NULL; root = root->next)
    {
      k = root->k;
      if (k == j)
        continue;
      cplx_sub (z, s->froot[j], s->froot[k]);
      cplx_inv_eq (z);
      cplx_add_eq (abcorr, z);
    }
}

/**
 * @brief Compute Aberth correction for the j-th root,
 * but only with other roots of the <code>jc</code>-th
 * cluster.
 */
void
mps_daberth_s (mps_status * s, int j, mps_cluster * cluster, cdpe_t abcorr)
{
  int k;
  mps_root * root;
  cdpe_t z;

  cdpe_set (abcorr, cdpe_zero);
  for (root = cluster->first; root != NULL; root = root->next)
    {
      k = root->k;
      if (k == j)
        continue;
      cdpe_sub (z, s->droot[j], s->droot[k]);
      cdpe_inv_eq (z);
      cdpe_add_eq (abcorr, z);
    }
}

/**
 * @brief Compute Aberth correction for the j-th root,
 * but only with other roots of the <code>jc</code>-th
 * cluster.
 */
void
mps_maberth_s (mps_status * s, int j, mps_cluster * cluster, mpc_t abcorr)
{
  int k;
  mps_root * root;
  cdpe_t z, temp;
  mpc_t diff;

  mpc_init2 (diff, s->mpwp);

  cdpe_set (temp, cdpe_zero);
  for (root = cluster->first; root != NULL; root = root->next)
    {
      k = root->k;
      if (k == j)
        continue;
      mpc_sub (diff, s->mroot[j], s->mroot[k]);
      mpc_get_cdpe (z, diff);
      cdpe_inv_eq (z);
      cdpe_add_eq (temp, z);
    }
  mpc_set_cdpe (abcorr, temp);

  mpc_clear (diff);
}

void
mps_maberth_s_wl (mps_status * s, int j, mps_cluster * cluster, mpc_t abcorr,
                  pthread_mutex_t * aberth_mutexes)
{
  int k;
  mps_root * root;
  cdpe_t z, temp;
  mpc_t diff;

  mpc_init2 (diff, s->mpwp);

  cdpe_set (temp, cdpe_zero);
  for (root = cluster->first; root != NULL; root = root->next)
    {
      k = root->k;
      if (k == j)
        continue;
      pthread_mutex_lock (&aberth_mutexes[k]);
      mpc_sub (diff, s->mroot[j], s->mroot[k]);
      pthread_mutex_unlock (&aberth_mutexes[k]);
      mpc_get_cdpe (z, diff);

      if (cdpe_eq_zero(z))
	continue;
      cdpe_inv_eq (z);
      cdpe_add_eq (temp, z);
    }
  mpc_set_cdpe (abcorr, temp);

  mpc_clear (diff);
}
