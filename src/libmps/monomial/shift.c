/*
 * This file is part of MPSolve 3.1.6
 *
 * Copyright (C) 2001-2015, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Giuseppe Fiorentino <fiorent@dm.unipi.it>
 *   Leonardo Robol <leonardo.robol@sns.it>
 */

#include <mps/mps.h>

static void
mps_raisetemp (mps_context * s, unsigned long int digits)
{
  int i;

  for (i = 0; i <= s->n; i++)
    {
      mpc_set_prec (s->mfpc1[i], digits);
      mpc_set_prec (s->mfppc1[i], digits);
    }
}

/**
 * @brief This routine computes the first \f$m+1\f$ coefficients of the shifted
 * polynomial \f$p(x+g)\f$, by performing \f$m+1\f$ Horner divisions.
 * This if the floating point version of this function.
 *
 * @param s The current mps_context.
 * @param m The size of the cluster.
 * @param mps_cluster_item A pointer to the cluster that shall be shifted.
 * @param clust_rad A bound for the radius of the cluster.
 * @param g The gravity center of the cluster.
 * @parma eps The current value of epsilon that should be used as a treshold.
 *
 * Then it computes the new starting approximations for the
 * cluster selected by applying mps_fstart() and by updating the approximations.
 */
MPS_PRIVATE void
mps_fshift (mps_context * s, int m, mps_cluster_item * cluster_item, double clust_rad,
            cplx_t g, rdpe_t eps)
{
  int i, j;
  double ag;
  cplx_t t;
  mps_monomial_poly *p = MPS_MONOMIAL_POLY (s->active_poly);

  /* Perform divisions */
  ag = cplx_mod (g);
  for (i = 0; i <= s->n; i++)
    cplx_set (s->fppc1[i], p->fpc[i]);
  for (i = 0; i <= m; i++)
    {
      cplx_set (t, s->fppc1[s->n]);
      for (j = s->n - 1; j >= i; j--)
        {
          cplx_mul_eq (t, g);
          cplx_add_eq (t, s->fppc1[j]);
          cplx_set (s->fppc1[j], t);
        }
      cplx_set (p->fppc[i], t);
    }

  /* start */
  for (i = 0; i <= m; i++)
    s->fap1[i] = cplx_mod (p->fppc[i]);

  /* If there is a custom starting point function use it, otherwise
   * use the default one */
  mps_fstart (s, m, cluster_item, clust_rad, ag, eps, s->fap1);
}

/**
 * @brief This routine computes the first \f$m+1\f$ coefficients of the shifted
 * polynomial \f$p(x+g)\f$, by performing \f$m+1\f$ Horner divisions.
 * This if the DPE version of this function.
 *
 * @param s The current mps_context.
 * @param m The size of the cluster.
 * @param mps_cluster_item A pointer to the cluster that shall be shifted.
 * @param clust_rad A bound for the radius of the cluster.
 * @param g The gravity center of the cluster.
 * @parma eps The current value of epsilon that should be used as a treshold.
 *
 * Then it computes the new starting approximations for the
 * cluster selected by applying mps_fstart() and by updating the approximations.
 */
MPS_PRIVATE void
mps_dshift (mps_context * s, int m, mps_cluster_item * cluster_item, rdpe_t clust_rad,
            cdpe_t g, rdpe_t eps)
{
  int i, j;
  rdpe_t ag;
  cdpe_t t;
  mps_monomial_poly * p = MPS_MONOMIAL_POLY (s->active_poly);

  cdpe_mod (ag, g);
  for (i = 0; i <= s->n; i++)
    cdpe_set (s->dpc1[i], p->dpc[i]);
  for (i = 0; i <= m; i++)
    {
      cdpe_set (t, s->dpc1[s->n]);
      for (j = s->n - 1; j >= i; j--)
        {
          cdpe_mul_eq (t, g);
          cdpe_add_eq (t, s->dpc1[j]);
          cdpe_set (s->dpc1[j], t);
        }
      cdpe_set (s->dpc2[i], t);
    }

  /* start */
  for (i = 0; i <= m; i++)
    cdpe_mod (s->dap1[i], s->dpc2[i]);

  mps_dstart (s, m, cluster_item, clust_rad, ag, eps, s->dap1);
}

/**
 * @brief This routine computes the first \f$m+1\f$ coefficients of the shifted
 * polynomial \f$p(x+g)\f$, by performing \f$m+1\f$ Horner divisions.
 * This if the MP version of this function.
 *
 * @param s The current mps_context.
 * @param m The size of the cluster.
 * @param mps_cluster_item A pointer to the cluster that shall be shifted.
 * @param clust_rad A bound for the radius of the cluster.
 * @param g The gravity center of the cluster.
 * @parma eps The current value of epsilon that should be used as a treshold.
 *
 * Then it computes the new starting approximations for the
 * cluster selected by applying mps_fstart() and by updating the approximations.
 */
MPS_PRIVATE void
mps_mshift (mps_context * s, int m, mps_cluster_item * cluster_item, rdpe_t clust_rad, mpc_t g)
{
  int i, j, k;
  long int mpwp_temp, mpwp_max;
  rdpe_t ag, ap, abp, as, mp_ep;
  cdpe_t abd;
  mpc_t t;
  mps_monomial_poly *p = MPS_MONOMIAL_POLY (s->active_poly);

  /* mps_cluster * cluster = cluster_item->cluster;   */

  mpc_init2 (t, s->mpwp);

  /* Perform divisions
   * In the mp version of the shift stage the computation
   * is performed with increasing levels of working precision
   * until the coefficients of the shifted polynomial have at
   * least one correct bit. */
  rdpe_set (mp_ep, s->mp_epsilon);
  mpc_get_cdpe (abd, g);
  cdpe_mod (ag, abd);
  for (i = 0; i <= s->n; i++)
    mpc_set (s->mfpc1[i], p->mfpc[i]);
  rdpe_set (as, rdpe_zero);
  rdpe_set (ap, rdpe_one);
  mpc_set_ui (t, 0, 0);
  k = 0;

  /* store the current working precision mpnw into mpnw_tmp */
  mpwp_temp = s->mpwp;
  mpwp_max = m * s->mpwp;

  do
    {                           /* loop */
      mpc_set (t, s->mfpc1[MPS_POLYNOMIAL (p)->degree]);
      mpc_get_cdpe (abd, p->mfpc[s->n]);
      cdpe_mod (ap, abd);
      for (j = s->n - 1; j >= 0; j--)
        {
          mpc_get_cdpe (abd, p->mfpc[j]);
          cdpe_mod (abp, abd);
          rdpe_mul_eq (ap, ag);
          rdpe_mul_eq_d (abp, (double)j);
          rdpe_add_eq (ap, abp);
          mpc_mul_eq (t, g);
          mpc_add_eq (t, s->mfpc1[j]);
          mpc_set (s->mfpc1[j], t);
        }

      mpc_set (s->mfppc1[0], t);
      mpc_get_cdpe (abd, t);
      cdpe_mod (as, abd);
      rdpe_mul_eq (ap, mp_ep);
      rdpe_mul_eq_d (ap, 4.0 * (s->n + 1));
      k++;

      if (rdpe_lt (as, ap))
        {
          mpwp_temp += s->mpwp;

          /* if ((mpwp_temp > mpwp_max || mpwp_temp > s->output_config->prec * m * 2))    */
          /*   {    */
          /*     MPS_DEBUG (s, "Reached the maximum allowed precision in mshift");    */
          /*     break;    */
          /*   } */

          rdpe_set_2dl (mp_ep, 1.0, 1 - mpwp_temp);
          mps_raisetemp (s, mpwp_temp);
          mpc_set_prec (t, (unsigned long int)mpwp_temp);
          mpc_set_prec (g, (unsigned long int)mpwp_temp);
          if (mpwp_max < mpwp_temp)
            mpwp_max = mpwp_temp;

          for (j = 0; j <= s->n; j++)
            mpc_set (s->mfpc1[j], p->mfpc[j]);
        }
    } while (rdpe_lt (as, ap) && (k <= m)); /* loop */

  mps_raisetemp (s, 1 * mpwp_temp);

  for (i = 1; i <= m; i++)
    {
      /* mpwp_temp = MAX (mpwp_temp - s->mpwp, s->mpwp); */
      /* mps_raisetemp (s, mpwp_temp); */
      /* mpc_set_prec (t, (unsigned long int) mpwp_temp); */
      /* mpc_set_prec (g, (unsigned long int) mpwp_temp); */
      mpc_set (t, s->mfpc1[s->n]);

      for (j = s->n - 1; j >= i; j--)
        {
          mpc_mul_eq (t, g);
          mpc_add_eq (t, s->mfpc1[j]);
          mpc_set (s->mfpc1[j], t);
        }
      mpc_set (s->mfppc1[i], t);
    }
  /*
     raisetemp_raw(mpwp);
     mpc_set_prec_raw(s, (unsigned long int) mpwp);
     mpc_set_prec_raw(g, (unsigned long int) mpwp);

     segue alternativa
   */
  /* mps_raisetemp (s, 2 * mpwp_max); */
  /* mpc_set_prec (t, (unsigned long int) 2 * mpwp_max); */
  /* mpc_set_prec (g, (unsigned long int) 2 * mpwp_max); */
  mps_raisetemp (s, 2 * mpwp_temp);
  mpc_set_prec (t, (unsigned long int)s->mpwp);
  mpc_set_prec (g, (unsigned long int)s->mpwp);

  if (rdpe_lt (as, ap))
    {
      for (j = 0; j < m; j++)
        rdpe_set (s->dap1[j], ap);
      mpc_get_cdpe (abd, s->mfppc1[m]);
      cdpe_mod (s->dap1[m], abd);
    }
  else
    for (i = 0; i <= m; i++)
      {
        mpc_get_cdpe (abd, s->mfppc1[i]);
        cdpe_mod (s->dap1[i], abd);
      }

  /* Debug the coefficients of the shifted polynomial */
  if (s->debug_level & MPS_DEBUG_CLUSTER)
    for (i = 0; i <= m; i++)
      MPS_DEBUG_MPC (s, mpc_get_prec (s->mfppc1[i]), s->mfppc1[i],
                     "P(x + g), coefficient of degree %d", i);

  mps_mstart (s, m, cluster_item, clust_rad, ag, s->dap1, g);

  mpc_clear (t);
}
