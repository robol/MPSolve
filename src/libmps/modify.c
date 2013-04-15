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

void
mps_cluster_detect_properties (mps_context * ctx, mps_cluster * cluster, mps_phase phase)
{
  mps_root * root = cluster->first;
  rdpe_t log_rad;
  mps_boolean (*touch_check) (mps_context * , int, int);

  if (ctx->output_config->root_properties & MPS_OUTPUT_PROPERTY_REAL)
    {
      /* Select the correct touch routine for us */
      switch (phase)
      {
        case float_phase:
          touch_check = mps_ftouchreal;
          break;
        case dpe_phase:
          touch_check = mps_dtouchreal;
          break;
        case mp_phase:
          touch_check = mps_mtouchreal;
          break;
        default:
          return;
      }

      /* For the moment we handle only the case of isolated roots. */
      if (cluster->n == 1)
        {
          mps_boolean touch_real_axis = touch_check (ctx, ctx->n, root->k);
          if (MPS_STRUCTURE_IS_REAL (ctx->active_poly->structure))
            ctx->root[root->k]->attrs = touch_real_axis ? MPS_ROOT_ATTRS_REAL : MPS_ROOT_ATTRS_NOT_REAL;

          if (phase == float_phase)
            rdpe_set_d (log_rad, ctx->root[root->k]->frad);
          else
            rdpe_set (log_rad, ctx->root[root->k]->drad);

          /* General case */
          if (touch_real_axis && rdpe_log (log_rad) < ctx->sep - ctx->n * ctx->lmax_coeff)
            ctx->root[root->k]->attrs = MPS_ROOT_ATTRS_REAL;
        }
     }

  if (ctx->output_config->root_properties & MPS_OUTPUT_PROPERTY_IMAGINARY)
    {
      /* Select the correct touch routine for us */
      switch (phase)
        {
          case float_phase:
            touch_check = mps_ftouchimag;
            break;
          case dpe_phase:
            touch_check = mps_dtouchimag;
            break;
          case mp_phase:
            touch_check = mps_mtouchimag;
            break;
          default:
            return;
        }
      
      for (root = cluster->first; root != NULL; root = root->next)
        {        
          if (phase == float_phase)
            rdpe_set_d (log_rad, ctx->root[root->k]->frad);
          else
            rdpe_set (log_rad, ctx->root[root->k]->drad);

          /* General case */
          if (touch_check (ctx, ctx->n, root->k) && rdpe_log (log_rad) < ctx->sep - ctx->n * ctx->lmax_coeff)
            ctx->root[root->k]->attrs = MPS_ROOT_ATTRS_IMAG;
        }
    }
}

/**
 * @brief Modify the vector 'status' according to the goal, and
 * to the location of the roots.
 *
 * @param s The mps_context associated to the current computation.
 * @param track_new_cluster true if old clusters should be marked
 * with 'C' instead of 'c', so they are recognizable (for shifting).
 *
 * The subroutine is used also for marking the new cluster
 * that have been detected between two consecutive packets
 * of Aberth's iteration.
 * 
 * -# The subroutine changes into 'C' the components of
 * status[:1) corresponding to old clusters, keeping
 * status[:,1)='c' for the new formed clusters.
 * In this way applying restart selects new starting
 * approximations only for the new detected clusters.
 *
 * -# For the components for which 
 * status[:,1]!='C', 'f', 'x' performs the following
 * analysis:
 * If the cluster has mult=1 mark it with status[:1)='i'
 * if is also approximated mark it with status(:1)='a'
 * Check if c*u and i*u (i.e., uncertain set) can 
 * be made certain according to goal[1]
 *
 * -# Perform the same with options, that is,
 * If multiplicity is on then check if a cluster
 * corresponds to a  multiple root
 * If detect real then detect real roots
 * if detect imaginary then detect imaginary roots
 * If detect both then detect both imaginary and
 * real roots
 */
void
mps_fmodify (mps_context * s, mps_boolean track_new_cluster)
{
  s->operation = MPS_OPERATION_CLUSTER_ANALYSIS;

  int l, i;
  rdpe_t rtmp;
  double eps_out = rdpe_get_d (s->eps_out);

  /* Set isolation factor */
  /* nf = 2 * s->n; */

  /* If tracking of new cluster is enabled we need to mark old clusters
   * with 'C' to distinguish them from the new ones. */
  if (track_new_cluster)
    {
      for (i = 0; i < s->n; i++)
        if (s->root[i]->status == MPS_ROOT_STATUS_CLUSTERED)
          s->root[i]->status = MPS_ROOT_STATUS_NEW_CLUSTERED;
    }

  /* Iterate over the cluster to update the status of the roots */
  mps_cluster_item * c_item;
  mps_cluster * cluster;

  c_item = s->clusterization->first;
  while (c_item != NULL)
    {
      mps_root * root;
      cluster = c_item->cluster;

      mps_cluster_detect_properties (s, cluster, float_phase);
      
      /* Pick the first root in the cluster */
      root = cluster->first;
      l = root->k;

      /* Check if this is an isolated cluster */
      if (cluster->n == 1)
        {
          /* Check if the root is already approximated; if that's not the
           * case set it at least as isolated. */
          if (s->root[l]->status != MPS_ROOT_STATUS_APPROXIMATED)
            {
              s->root[l]->status = MPS_ROOT_STATUS_ISOLATED;

              /* Check if we need to mark this root as approximated */
              if (s->root[l]->frad < cplx_mod (s->root[l]->fvalue) * eps_out)
                s->root[root->k]->status = MPS_ROOT_STATUS_APPROXIMATED;
            }

          /* Grab the next cluster and continue scanning */
          c_item = c_item->next;
          continue;
        }

      /* If it's not the case scan the roots in the cluster and set them
       * to 'c'. */
      while (root != NULL)
        {
          l = root->k;

          /* If track_new_cluster is false then we may directly set here the
           * approximation status of the roots. */
          if (!track_new_cluster)
            {     
              s->root[l]->status = MPS_ROOT_STATUS_CLUSTERED;
            }

          rdpe_set_d (rtmp, s->root[l]->frad);
          rdpe_div_eq_d (rtmp, cplx_mod (s->root[l]->fvalue));
          if (rdpe_le (rtmp, s->eps_out))
            s->root[l]->status = MPS_ROOT_STATUS_APPROXIMATED_IN_CLUSTER;

          root = root->next;
        }

      /* TODO: Implement checking of the zone where the roots are. */

      c_item = c_item->next;
    }

  mps_fupdate_inclusions (s);
}


/**
 * @brief The DPE version of <code>mps_fmodify()</code>.
 *
 * @param s The mps_context associated to the current computation.
 * @param track_new_cluster true if old clusters should be marked
 * with 'C' instead of 'c', so they are recognizable (for shifting).
 *
 * @see mps_fmodify()
 */
void
mps_dmodify (mps_context * s, mps_boolean track_new_cluster)
{
  s->operation = MPS_OPERATION_CLUSTER_ANALYSIS;

  int l, i;
  rdpe_t tmpr, tmpr2;

  /* Set isolation factor */
  /* nf = 2 * s->n; */

  /* If tracking of new cluster is enabled we need to mark old clusters
   * with 'C' to distinguish them from the new ones. */
  if (track_new_cluster)
    {
      for (i = 0; i < s->n; i++)
        if (s->root[i]->status == MPS_ROOT_STATUS_CLUSTERED)
          s->root[i]->status = MPS_ROOT_STATUS_NEW_CLUSTERED;
    }

  /* Iterate over the cluster to update the status of the roots */
  mps_cluster_item * c_item;
  mps_cluster * cluster;

  c_item = s->clusterization->first;
  while (c_item != NULL)
    {
      mps_root * root;
      cluster = c_item->cluster;

      mps_cluster_detect_properties (s, cluster, dpe_phase);
      
      /* Pick the first root in the cluster */
      root = cluster->first;
      l = root->k;

      /* Check if this is an isolated cluster */
      if (cluster->n == 1)
        {
          /* Check if the root is already approximated; if that's not the
           * case set it at least as isolated. */
          if (s->root[l]->status != MPS_ROOT_STATUS_APPROXIMATED)
            s->root[l]->status = MPS_ROOT_STATUS_ISOLATED;

          /* Grab the next cluster and continue scanning */
          c_item = c_item->next;
          continue;
        }

      /* If it's not the case scan the roots in the cluster and set them
       * to 'c'. */
      while (root != NULL)
        {
          l = root->k;

          /* If track_new_cluster is false then we may directly set here the
           * approximation status of the roots. */
          if (!track_new_cluster)
            {
              s->root[l]->status = MPS_ROOT_STATUS_CLUSTERED;
            }

          rdpe_set (tmpr, s->root[l]->drad);
          cdpe_mod (tmpr2, s->root[l]->dvalue);
          rdpe_div_eq (tmpr, tmpr2);
          if (rdpe_le (tmpr, s->eps_out)) 
            {
              s->root[l]->status = MPS_ROOT_STATUS_APPROXIMATED_IN_CLUSTER;
            }

          root = root->next;
        }

      /* TODO: Implement checking of the zone where the roots are. */
      c_item = c_item->next;
    }

  mps_dupdate_inclusions (s);
}


/**
 * @brief The multiprecision version of the routine
 * <code>mps_fmodify()</code>. 
 *
 * @param s The mps_context associated to the current computation.
 * @param track_new_cluster true if old clusters should be marked
 * with 'C' instead of 'c', so they are recognizable (for shifting).
 *
 * @see mps_fmodify()
 */
void
mps_mmodify (mps_context * s, mps_boolean track_new_cluster)
{
  s->operation = MPS_OPERATION_CLUSTER_ANALYSIS;

  int l, i;
  rdpe_t tmpr, tmpr2;
  cdpe_t cdtmp;

  /* If tracking of new cluster is enabled we need to mark old clusters
   * with 'C' to distinguish them from the new ones. */
  if (track_new_cluster)
    {
      for (i = 0; i < s->n; i++)
        if (s->root[i]->status == MPS_ROOT_STATUS_CLUSTERED)
          s->root[i]->status = MPS_ROOT_STATUS_NEW_CLUSTERED;
    }

  /* Iterate over the cluster to update the status of the roots */
  mps_cluster_item * c_item;
  mps_cluster * cluster;

  c_item = s->clusterization->first;
  while (c_item != NULL)
    {
      mps_root * root;
      cluster = c_item->cluster;

      mps_cluster_detect_properties (s, cluster, mp_phase);
      
      /* Pick the first root in the cluster */
      root = cluster->first;
      l = root->k;

      /* Check if this is an isolated cluster */
      if (cluster->n == 1)
        {
          /* Check if the root is already approximated; if that's not the
           * case set it at least as isolated. */
          if (s->root[l]->status != MPS_ROOT_STATUS_APPROXIMATED)
            s->root[l]->status = MPS_ROOT_STATUS_ISOLATED;

          /* Grab the next cluster and continue scanning */
          c_item = c_item->next;
          continue;
        }

      /* If it's not the case scan the roots in the cluster and set them
       * to 'c'. */
      while (root != NULL)
        {
          l = root->k;

          /* If track_new_cluster is false then we may directly set here the
           * approximation status of the roots. */
          if (!track_new_cluster)
            {
              s->root[l]->status = MPS_ROOT_STATUS_CLUSTERED;
            }

          rdpe_set (tmpr, s->root[l]->drad);
          mpc_get_cdpe (cdtmp, s->root[l]->mvalue);
          cdpe_mod (tmpr2, cdtmp);
          rdpe_div_eq (tmpr, tmpr2);
          if (rdpe_le (tmpr, s->eps_out)) 
            s->root[l]->status = MPS_ROOT_STATUS_APPROXIMATED_IN_CLUSTER; 

          root = root->next;
        }

      /* TODO: Implement checking of the zone where the roots are. */
      c_item = c_item->next;
    }

  mps_mupdate_inclusions (s);
}
