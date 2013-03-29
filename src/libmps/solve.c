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


#include <math.h>
#include <mps/mps.h>
#include <float.h>

/**
 * @brief Set <code>again[i]</code> to <code>true</code> or to <code>false</code> 
 * according to the values of <code>status</code> and <code>inclusion</code>
 *  in the mps_approximation and the current <code>goal</code>.
 *
 * More precisely:
 *
 * - If goal is <code>MPS_OUTPUT_GOAL_COUNT</code>: <code>true</code> only 
 *   if inclusion is <code>MPS_ROOT_INCLUSION_UNKNOWN</code> and status is not 
 *   <code>MPS_ROOT_STATUS_APPROXIMATED</code>, <code>MPS_ROOT_STATUS_APPROXIMATED_IN_CLUSTER</code>
 *   or <code>MPS_ROOT_STATUS_NOT_FLOAT</code>.
 *
 * - If goal is <code>MPS_OUTPUT_GOAL_ISOLATE</code> or <code>MPS_OUTPUT_GOAL_APPROXIMATE</code>:
 *   <code>true</code> if <code>status</code> is <code>MPS_ROOT_STATUS_CLUSTERED</code> or <code>inclusion</code>
 *   is <code>MPS_ROOT_INCLUSION_UNKNOWN</code>.
 */
void
mps_update (mps_context * s)
{
  int i;

  for (i = 0; i < s->n; i++)
    s->root[i]->again = false;
  switch (s->output_config->goal)
    {

    case MPS_OUTPUT_GOAL_COUNT:                  /*  count */
      for (i = 0; i < s->n; i++)
        {
          if (s->root[i]->inclusion == MPS_ROOT_INCLUSION_UNKNOWN)
            if (!MPS_ROOT_STATUS_IS_APPROXIMATED (s->root[i]->status) && 
                s->root[i]->status != MPS_ROOT_STATUS_NOT_DPE)
              s->root[i]->again = true;
          
          if (s->output_config->multiplicity && 
              s->root[i]->status == MPS_ROOT_STATUS_CLUSTERED &&
              s->root[i]->inclusion != MPS_ROOT_INCLUSION_OUT)
            s->root[i]->again = true;

            if (s->output_config->root_properties &&
                (s->root[i]->attrs == MPS_ROOT_ATTRS_NONE &&
                 (s->root[i]->inclusion != MPS_ROOT_INCLUSION_UNKNOWN ||
                  (!MPS_ROOT_STATUS_IS_APPROXIMATED (s->root[i]->status) && 
                   s->root[i]->status != MPS_ROOT_STATUS_NOT_DPE))))
              s->root[i]->again = true;
        }
      break;
          
    case MPS_OUTPUT_GOAL_ISOLATE:                  /* isolate */
      for (i = 0; i < s->n; i++)
        {
          if (s->root[i]->inclusion == MPS_ROOT_INCLUSION_UNKNOWN ||
              (s->root[i]->status == MPS_ROOT_STATUS_CLUSTERED && 
               s->root[i]->inclusion == MPS_ROOT_INCLUSION_IN))
            if (!MPS_ROOT_STATUS_IS_APPROXIMATED (s->root[i]->status) &&
                (s->root[i]->status != MPS_ROOT_STATUS_ISOLATED ||
                 s->root[i]->inclusion != MPS_ROOT_INCLUSION_IN))
              s->root[i]->again = true;

          if (s->output_config->multiplicity && 
              s->root[i]->status == MPS_ROOT_STATUS_CLUSTERED &&
              s->root[i]->inclusion != MPS_ROOT_INCLUSION_OUT)
            s->root[i]->again = true;

          if (s->output_config->root_properties)
            {
              if (s->root[i]->attrs == MPS_ROOT_ATTRS_NONE &&
                  (!MPS_ROOT_STATUS_IS_APPROXIMATED (s->root[i]->status) &&
                   s->root[i]->status != MPS_ROOT_STATUS_NOT_DPE))
                s->root[i]->again = true;
            }
        }

      break;

    case MPS_OUTPUT_GOAL_APPROXIMATE:                  /* approximate (the same as isolate) */
      for (i = 0; i < s->n; i++)
        {
          if (s->root[i]->inclusion == MPS_ROOT_INCLUSION_UNKNOWN ||
              (s->root[i]->status == MPS_ROOT_STATUS_CLUSTERED &&
               s->root[i]->inclusion == MPS_ROOT_INCLUSION_IN))
            if (!MPS_ROOT_STATUS_IS_APPROXIMATED (s->root[i]->status))
              s->root[i]->again = true;

          if (s->output_config->multiplicity &&
              s->root[i]->status == MPS_ROOT_STATUS_CLUSTERED &&
              s->root[i]->inclusion != MPS_ROOT_INCLUSION_OUT)
            s->root[i]->again = true;

          if (s->output_config->root_properties)
            {
              if (s->root[i]->attrs == MPS_ROOT_ATTRS_NONE && 
                  MPS_ROOT_STATUS_IS_APPROXIMATED (s->root[i]->status))
                s->root[i]->again = true;
            }
        }
      break;
    }
}

/**
 * @brief Compute super center and super radius
 * 
 * This routines the super radius of the <code>i</code>-th cluster,
 * i.e. the radius of the inclusion disc for the whole cluster
 *
 * @param s The <code>mps_context</code> associated with the current computation.
 * @param cluster cluster of which the super center and radius should be computed.
 * @param sc Center of the cluster;
 * @param sr Double that will be set to the super radius of the cluster;
 */
void
mps_fsrad (mps_context * s, mps_cluster * cluster, cplx_t sc, double *sr)
{
  cplx_t ctmp;
  double sum;
  int l;

  mps_root * root;

  sum = 0.0;
  for (root = cluster->first; root != NULL; root = root->next)
    {
      l = root->k;
      sum += s->root[l]->frad;
    }
  cplx_set (sc, cplx_zero);
  for (root = cluster->first; root != NULL; root = root->next)
    {
      l = root->k;
      cplx_mul_d (ctmp, s->root[l]->fvalue, s->root[l]->frad);
      cplx_add_eq (sc, ctmp);
    }
  cplx_div_eq_d (sc, sum);
  *sr = 0.0;
  for (root = cluster->first; root != NULL; root = root->next)
    {
      l = root->k;
      cplx_sub (ctmp, sc, s->root[l]->fvalue);
      *sr = MAX (*sr, s->root[l]->frad + cplx_mod (ctmp));
    }
}

/**
 * @brief <code>dpe</code> version of <code>fsrad()</code>
 */
void
mps_dsrad (mps_context * s, mps_cluster * cluster, cdpe_t sc, rdpe_t sr)
{
  cdpe_t ctmp;
  rdpe_t sum, rtmp;
  int l;
  mps_root * root;

  rdpe_set (sum, rdpe_zero);
  for (root = cluster->first; root != NULL; root = root->next)
    {
      l = root->k;
      rdpe_add_eq (sum, s->root[l]->drad);
    }
  cdpe_set (sc, cdpe_zero);
  for (root = cluster->first; root != NULL; root = root->next)
    {
      l = root->k;
      cdpe_mul_e (ctmp, s->root[l]->dvalue, s->root[l]->drad);
      cdpe_add_eq (sc, ctmp);
    }
  cdpe_div_eq_e (sc, sum);
  rdpe_set (sr, rdpe_zero);
  for (root = cluster->first; root != NULL; root = root->next)
    {
      l = root->k;
      cdpe_sub (ctmp, sc, s->root[l]->dvalue);
      cdpe_mod (rtmp, ctmp);
      rdpe_add_eq (rtmp, s->root[l]->drad);
      if (rdpe_lt (sr, rtmp))
        rdpe_set (sr, rtmp);
    }
}

/**
 * @brief Multiprecision versione of <code>fsrad()</code>
 */
void
mps_msrad (mps_context * s, mps_cluster * cluster, mpc_t sc, rdpe_t sr)
{
  int l;
  rdpe_t rtmp;
  cdpe_t cdtmp;
  mpf_t ftmp, sum;
  mpc_t ctmp;
  mps_root * root;

  mpc_init2 (ctmp, s->mpwp);
  mpf_init2 (ftmp, s->mpwp);
  mpf_init2 (sum, s->mpwp);

  mpf_set_ui (sum, 0);
  for (root = cluster->first; root != NULL; root = root->next)
    {
      l = root->k;
      mpf_set_rdpe (ftmp, s->root[l]->drad);
      mpf_add (sum, sum, ftmp);
    }

  mpc_set_ui (sc, 0, 0);
  for (root = cluster->first; root != NULL; root = root->next)
    {
      l = root->k;
      mpf_set_rdpe (ftmp, s->root[l]->drad);
      mpc_mul_f (ctmp, s->root[l]->mvalue, ftmp);
      mpc_add_eq (sc, ctmp);
    }

  mpc_div_eq_f (sc, sum);
  rdpe_set (sr, rdpe_zero);
  for (root = cluster->first; root != NULL; root = root->next)
    {
      l = root->k;
      mpc_sub (ctmp, sc, s->root[l]->mvalue);
      mpc_get_cdpe (cdtmp, ctmp);
      cdpe_mod (rtmp, cdtmp);
      rdpe_add_eq (rtmp, s->root[l]->drad);
      if (rdpe_lt (sr, rtmp))
          rdpe_set (sr, rtmp);
      else
        {
          if (s->debug_level & MPS_DEBUG_CLUSTER) 
            {
              MPS_DEBUG_RDPE (s, sr, "Super radius of the cluster");
              // MPS_DEBUG_RDPE (s, rtmp, "rtmp");
            }
        }
    }

  mpf_clear (sum);
  mpf_clear (ftmp);
  mpc_clear (ctmp);
}



/**
 * @brief Check if the roots are computed with the required
 * precision.
 *
 * Set <code>computed</code> to <code>true</code> 
 * if the stop condition is satisfied,
 * otherwise set <code>computed</code> to <code>false</code>.
 *
 * The stop condition is obtained from the vector <code>status</code> 
 * as follows:
 *
 * If the <code>goal</code> is count stop if
 *  -  <code>**u</code> does not exist, except for 
 *     <code>a*u</code>, <code>o*u</code>, <code>f*u</code>
 *  - Mult. and does not exist <code>c**</code>
 *  - Real. and does not exist <code>*u*</code>, except for <code>au*</code>, 
 *    <code>ou*</code>;
 *  - Imag  and does not exist <code>*v*</code>, except for 
 *    <code>av*</code>, <code>ov*</code>;
 *
 * If the <code>goal</code> is isolate or approximate stop if:
 * - <code>**u</code> does not exist, except for  <code>a*u</code>, 
 *   <code>o*u</code>, <code>f*u</code>
 *   and if <code>c*i</code> does not exist;
 * - Mult. and does not exist <code>c*i</code>, <code>o*i</code>
 * - Real. and does not exist <code>*ui</code>, except for <code>aui</code>, 
 *   <code>oui</code>;
 * - Imag  and does not exist <code>*vi</code>, 
 *   except for <code>avi</code>, <code>ovi</code>;
 * 
 * @see status
 */
mps_boolean
mps_check_stop (mps_context * s)
{
  MPS_DEBUG_THIS_CALL;

  int i;
  mps_boolean computed;

  computed = false;
  /* count */
  if (s->output_config->goal == MPS_OUTPUT_GOAL_COUNT)
    {
      for (i = 0; i < s->n; i++)
        {
          if (!MPS_ROOT_STATUS_IS_APPROXIMATED (s->root[i]->status) &&
              s->root[i]->inclusion == MPS_ROOT_INCLUSION_UNKNOWN)
            return computed;

          if (s->output_config->multiplicity && 
              s->root[i]->status == MPS_ROOT_STATUS_CLUSTERED &&
              s->root[i]->inclusion != MPS_ROOT_INCLUSION_OUT)
            return computed;

          if (s->output_config->root_properties &&
              s->root[i]->attrs == MPS_ROOT_ATTRS_NONE &&
              s->root[i]->inclusion != MPS_ROOT_INCLUSION_OUT &&
              !MPS_ROOT_STATUS_IS_APPROXIMATED (s->root[i]->status) && 
              s->root[i]->status != MPS_ROOT_STATUS_MULTIPLE)
            return computed;
        }
      computed = true;
    }

  /* isolate or approximate */
  if (s->output_config->goal == MPS_OUTPUT_GOAL_ISOLATE ||
      s->output_config->goal == MPS_OUTPUT_GOAL_APPROXIMATE)
    {
      for (i = 0; i < s->n; i++)
        {
          if ((s->root[i]->inclusion == MPS_ROOT_INCLUSION_UNKNOWN || s->root[i]->inclusion == MPS_ROOT_INCLUSION_IN) &&
              !MPS_ROOT_STATUS_IS_COMPUTED (s->root[i]->status))
            {
              if (s->debug_level & MPS_DEBUG_PACKETS)
                MPS_DEBUG (s, "Cannot stop because root %d is not approximated, and its inclusion is not certain", i);
              return computed;
            }

          if (s->root[i]->status == MPS_ROOT_STATUS_CLUSTERED &&
              s->root[i]->inclusion != MPS_ROOT_INCLUSION_OUT)
            {
              if (s->debug_level & MPS_DEBUG_PACKETS)
                MPS_DEBUG (s, "Cannot stop because root %d is clustered and not certainly out of the target set", i);
              return computed;
            }

          if (s->output_config->multiplicity &&
              s->root[i]->inclusion != MPS_ROOT_INCLUSION_OUT &&
              s->root[i]->status == MPS_ROOT_STATUS_CLUSTERED)
            {
              if (s->debug_level & MPS_DEBUG_PACKETS)
                MPS_DEBUG (s, "Cannot stop because root %d is not out of the target set and is clustered", i);
              return computed;
            }

          if (s->output_config->root_properties &&
              s->root[i]->attrs == MPS_ROOT_ATTRS_NONE && 
              s->root[i]->inclusion != MPS_ROOT_INCLUSION_OUT &&
              !MPS_ROOT_STATUS_IS_APPROXIMATED (s->root[i]->status) &&
              s->root[i]->status != MPS_ROOT_STATUS_MULTIPLE)
            {
              if (s->debug_level & MPS_DEBUG_PACKETS)
                MPS_DEBUG (s, "Cannot stop because properties of root %d have not been detected, it's not out of the target set, nor approximated or multiple", i);
              return computed;
            }
        }

      MPS_DEBUG (s, "All roots are computed, stopping Aberth");
      computed = true;
    }

  return computed;
}

/**
 * @brief Actually solve the polynomial
 *
 * This routine performs the following computations:
 * -# Select starting approximations and check if some of them need dpe.
 *    Initialize the vector again which is true if the corresponding
 *    approximation is out of the root neighbourhood.
 * -# Performs max_pack packets of Aberth iterations on all the
 *    components out of the root neighbourhood belonging to the set S
 *    and on which it is possible to iterate with float.
 *    More precisely, each packet performs max_it iterations on all the
 *    components where again is true. At each iteration check if the
 *    current approximation is in the root neighbourhood; in this case
 *    set 'again' to false.  
 * -# If at the end of the general packet all the approximations
 *    are inside the root neighbourhood, i.e., 'again' is false in all
 *    the components then return.
 *    else, perform cluster analysis, select new starting approximations
 *    update the vector 'statu', update the vector 'again' that selects 
 *    the components on which to iterate, according to the goal, and 
 *    repeat until the max number of allowed packets is reached. 
 *    In the latter case output FAILURE.
 *
 * The local variable <code>again</code> controls the iteration: i.e., 
 *   <code>again[i]=true</code> means iterate on the <code>i</code>-th 
 *   component
 *
 * @param s The <code>mps_context</code> associated with the current computation.
 * @param d_after_f this variable is <code>true</code> if dpe
 * are needed after the floating point pass. 
 */
void
mps_fsolve (mps_context * s, mps_boolean * d_after_f)
{
  mps_boolean excep;
  int it_pack, iter, nit, oldnclust, i, j, required_zeros = s->n;
  mps_monomial_poly *p = MPS_MONOMIAL_POLY (s->active_poly);
  double * frad = double_valloc (s->n);

  /* == 1 ==  Initialize variables */
  it_pack = 0;
  mps_cluster_reset (s);
  for (i = 0; i < s->n; i++)
    {
      s->root[i]->again = true;
      cplx_set (s->root[i]->fvalue, cplx_zero);
      s->root[i]->frad = DBL_MAX;
    }

  /* choose starting approximations */
  if (s->DOLOG)
    fprintf (s->logstr, "FSOLVE: call fstart");

  mps_polynomial_fstart (s, MPS_POLYNOMIAL (p));

        /***************
         this part of code performs shift in the gravity center of the roots
         In order to use it, uncomment the part below and comment the
         instruction above. Dangerous for overflow.
         ************/
  /*
     {
     cplx_t ft;
     cplx_mul_d(ft, fpc[n], -n);
     cplx_div(ft, fpc[n-1], ft);
     fshift(n, 0, 100, ft, eps_out);
     }
   *//* till here */

  if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)
    mps_dump (s);

  /* Check if there are too large or too small approximations */
  *d_after_f = false;
  for (i = 0; i < s->n; i++)
    if (s->root[i]->status == MPS_ROOT_STATUS_NOT_FLOAT)
      {
        s->root[i]->again = false;
        *d_after_f = true;
      }

  /* == 2 ==  Perform max_pack packets of Aberth's iterations */
  if (s->DOLOG)
    fprintf (s->logstr, "   FSOLVE:  call fpolzer\n");
  for (iter = 0; iter < s->max_pack; iter++)
    {                           /* floop: */

       /* mps_fpolzer(s, &nit, &excep);   */
      mps_thread_fpolzer (s, &nit, &excep, required_zeros--);
      it_pack += nit;

      if (s->DOLOG)
        fprintf (s->logstr, "Packet %d  iterations= %d\n", iter, nit);

      /* perform cluster analysis, shift, restart, update 'statu', and
       * update 'again'
       */
      if (excep)
        {
          oldnclust = s->clusterization->n;

          if (s->DOLOG)
            fprintf (s->logstr, "   FSOLVE: call fcluster\n");
          /* cluster analysis */

          /* Compute the inclusion radii with Gerschgorin so we can compute
           * clusterizations for the roots. */
          mps_fradii (s, s->active_poly, frad);
          mps_fcluster (s, frad, 2 * s->n);   /* Isolation factor */

          MPS_DEBUG (s, "oldncluster = %d, ncluster = %ld", oldnclust, s->clusterization->n);

          if (oldnclust == s->clusterization->n)
            {
              if (s->DOLOG)
                fprintf (s->logstr, "   FSOLVE: cycle\n");
              continue;
            }
          else
            {
              /* modify the vector status and mark also the old
               * clusters with 'C'
               */
              if (s->DOLOG)
                fprintf (s->logstr, "   FSOLVE: call modify\n");
              mps_fmodify (s, true);

              if (iter == 0)
                for (i = 0; i < s->n; i++)
                  if (s->root[i]->status == MPS_ROOT_STATUS_NEW_CLUSTERED)
                    s->root[i]->status = MPS_ROOT_STATUS_CLUSTERED;

              /* If the polynomial is not given in terms of its coeff. then
               * skip the restart stage */
              if (MPS_IS_MONOMIAL_POLY (s->active_poly))
                {
                  /* choose new starting approximations only for new clusters */
                  if (s->DOLOG)
                    fprintf (s->logstr, "   FSOLVE: call frestart\n");
                  mps_frestart (s);
                }
              /* reset the status vector */
              for (j = 0; j < s->n; j++)
                {
                  if (s->root[j]->status == MPS_ROOT_STATUS_NEW_CLUSTERED)
                    s->root[j]->status = MPS_ROOT_STATUS_CLUSTERED;
                  s->again_old[j] = s->root[j]->again;
                }

              /* update 'again' */
              if (s->DOLOG)
                fprintf (s->logstr, "   FSOLVE: call update\n");
              mps_update (s);

              /* adjust 'again' This is needed since we are
               * between two packets */
              for (i = 0; i < s->n; i++)
                if (!s->again_old[i])
                  s->root[i]->again = false;
              if (s->DOLOG)
                fprintf (s->logstr, "   FSOLVE: call checkstop\n");
              /* Check the stop condition */
              if (mps_check_stop (s))
                goto fsolve_final_cleanup;
            }
        }
      else
        break;
    }

  /* The 'floop' has been completed:
   * If the max number of iteration has been reached then output FAILURE */
  if (iter == s->max_pack)
    {
      mps_dump (s);
      mps_error (s, 1,
                 "Float: reached the maximum number of packet iterations");
    }
  /* Otherwise exit since all the approximations are
   * in the root neighbourhood, except for the ones that cannot be
   * represented as double. */

  if (s->DOLOG)
    fprintf (s->logstr, "FLOAT: nit= %d\n", it_pack);

  /* Update */
  if (s->DOLOG)
    fprintf (s->logstr, "   FSOLVE: call fcluster\n");
  oldnclust = s->clusterization->n;

  /* Compute the inclusion radii with Gerschgorin so we can compute
   * clusterizations for the roots. */
  mps_fradii (s, s->active_poly, frad);
  mps_fcluster (s, frad, 2 * s->n);   /* Isolation factor */

  if (s->DOLOG)
    fprintf (s->logstr, "   FSOLVE: call modify\n");
  mps_fmodify (s, true);

  /* reset the status vector */
  for (j = 0; j < s->n; j++)
    {
      if (s->root[j]->status == MPS_ROOT_STATUS_NEW_CLUSTERED)
        s->root[j]->status = MPS_ROOT_STATUS_CLUSTERED;
      if (s->root[j]->status == MPS_ROOT_STATUS_NOT_FLOAT)
        *d_after_f = true;
    }


 fsolve_final_cleanup:
  double_vfree (frad);
}

/**
 * @brief This routine applies <code>nit</code> iterations of 
 * Aberth's method.
 * 
 * The method is applied to
 * the <code>i</code>-th component of the approximations for which 
 * <code>again[i]</code>
 * is <code>true</code>. Set <code>again[i]=false</code> if the <code>i</code>-th 
 * approximation is in
 * the root neighbourhood. Stop if <code>again[i]=false</code> for any <code>i</code>.
 *
 * @param s The <code>mps_context</code> associated with the current computation.
 * @param it Index of the component on which the iteration is needed. 
 * @param excep This variable is set to <code>true</code> if after <code>nit</code> 
 * iterations some approximation is still out of the root neighbourhood.
 */
void
mps_fpolzer (mps_context * s, int *it, mps_boolean * excep)
{
  int i, iter, nzeros;
  cplx_t corr, abcorr;
  double rad1, modcorr;
  mps_monomial_poly * p = MPS_MONOMIAL_POLY (s->active_poly);

  /* initialize the iteration counter */
  *it = 0;
  *excep = false;

  /* count the number of approximations in the root neighbourhood */
  nzeros = 0;
  for (i = 0; i < s->n; i++)
    if (!s->root[i]->again)
      nzeros++;
  if (nzeros == s->n)
    return;

  /* Start Aberth's iterations */
  if (s->DOLOG)
    fprintf (s->logstr, "FPOLZER: starts aberth it\n");

  for (iter = 0; iter < s->max_it; iter++)
    {                           /* do_iter : DO iter=1,nit */

      if (s->DOLOG)
        {
          fprintf (s->logstr, "FPOLZER: iteration %d\n", iter);
          mps_dump (s);
        }

      for (i = 0; i < s->n; i++)
        {                       /* do_index */
          if (s->root[i]->again)
            {
              (*it)++;
              rad1 = s->root[i]->frad;
              mps_polynomial_fnewton (s, MPS_POLYNOMIAL (p), s->root[i], corr);
              if (iter == 0 && !s->root[i]->again && s->root[i]->frad > rad1 && rad1
                  != 0)
                s->root[i]->frad = rad1;
              /***************************************
               The above condition is needed to cope with the case
               where at the first iteration the starting point
               is already in the root neighbourhood and the actually
               computed radius is too big since the value of the first
               derivative is too small.
               In this case the previous radius bound, obtained by
               means of Rouche' is more reliable and strict
               **************************************/

              if (s->root[i]->again ||
                  /* the correction is performed only if iter!=1 or rad(i)!=rad1 */
                  iter != 0 || s->root[i]->frad != rad1)
                {
                  mps_faberth (s, s->root[i], abcorr);
                  cplx_mul_eq (abcorr, corr);
                  cplx_sub (abcorr, cplx_one, abcorr);
                  cplx_div (abcorr, corr, abcorr);
                  cplx_sub_eq (s->root[i]->fvalue, abcorr);
                  modcorr = cplx_mod (abcorr);
                  s->root[i]->frad += modcorr;
                }

              /* check for new approximated roots */
              if (!s->root[i]->again)
                {
                  nzeros++;
                  if (nzeros == s->n)
                    return;
                }

            }
        }
    }
  *excep = true;
}

/**
 * @brief <code>dpe</code> version of <code>fpolzer()</code>.
 */
void
mps_dpolzer (mps_context * s, int *it, mps_boolean * excep)
{
  int iter, i, nzeros;
  rdpe_t rad1, rtmp;
  cdpe_t corr, abcorr;
  mps_monomial_poly * p = MPS_MONOMIAL_POLY (s->active_poly);

  /* initialize the iteration counter */
  *it = 0;
  *excep = false;

  /* count the number of approximations in the root neighbourhood */
  nzeros = 0;
  for (i = 0; i < s->n; i++)
    if (!s->root[i]->again)
      nzeros++;
  if (nzeros == s->n)
    return;

  /* Start Aberth's iterations */
  if (s->DOLOG)
    fprintf (s->logstr, "DPOLZER: starts aberth\n");
  for (iter = 0; iter < s->max_it; iter++)
    {                           /* do_iter: */

      for (i = 0; i < s->n; i++)
        {                       /* do_index: */

          if (s->root[i]->again)
            {
              (*it)++;
              rdpe_set (rad1, s->root[i]->drad);
              mps_polynomial_dnewton (s, MPS_POLYNOMIAL (p), s->root[i], corr);
              if (iter == 0 && !s->root[i]->again && rdpe_gt (s->root[i]->drad, rad1)
                  && rdpe_ne (rad1, rdpe_zero))
                rdpe_set (s->root[i]->drad, rad1);
                                /************************************************
                                 The above condition is needed to manage with the case where
                                 at the first iteration the starting point is already in the
                                 root neighbourhood and the actually computed radius is too
                                 big since the value of the first derivative is too small.
                                 In this case the previous radius bound, obtained by means of
                                 Rouche' is more reliable and strict
                                 **********************************************/

              if (s->root[i]->again ||
                  /* the correction is performed only if iter!=1 or rad(i)!=rad1 */
                  iter != 0
                  || rdpe_ne (s->root[i]->drad, rad1))
                {
                  mps_daberth (s, s->root[i], abcorr);
                  cdpe_mul_eq (abcorr, corr);
                  cdpe_sub (abcorr, cdpe_one, abcorr);
                  cdpe_div (abcorr, corr, abcorr);
                  cdpe_sub_eq (s->root[i]->dvalue, abcorr);
                  cdpe_mod (rtmp, abcorr);
                  rdpe_add_eq (s->root[i]->drad, rtmp);
                }

              /* check for new approximated roots */
              if (!s->root[i]->again)
                {
                  nzeros++;
                  if (nzeros == s->n)
                    return;
                }

            }
        }
    }
  *excep = true;
}

/**
 * @brief <code>dpe</code> version of <code>fsolve()</code>.
 */
void
mps_dsolve (mps_context * s, mps_boolean d_after_f)
{
  int it_pack, iter, nit, oldnclust, i, j, required_zeros = s->n;
  mps_boolean excep;
  mps_monomial_poly * p = MPS_MONOMIAL_POLY (s->active_poly);
  rdpe_t * drad = rdpe_valloc (s->n);

  if (s->DOLOG)
    {
      fprintf (s->logstr, "   DSOLVE: d_after_f= %d\n", d_after_f);
    }

  /* == 1 == Initialize variables */
  it_pack = 0;

  if (d_after_f)
    for (i = 0; i < s->n; i++)
      if (s->root[i]->status == MPS_ROOT_STATUS_NOT_FLOAT)
        {
          s->root[i]->again = true;
          rdpe_set_d (s->root[i]->drad, DBL_MAX);
        }
      else
        s->root[i]->again = false;
  else
    {
      mps_cluster_reset (s);
      for (i = 0; i < s->n; i++)
        {
          s->root[i]->again = true;
          rdpe_set_d (s->root[i]->drad, DBL_MAX);
          cdpe_set (s->root[i]->dvalue, cdpe_zero);
        }
    }

  /* Choose starting approximations */
  if (s->DOLOG)
    fprintf (s->logstr, "   DSOLVE: call dstart con again=\n");

  mps_polynomial_dstart (s, MPS_POLYNOMIAL (p));

  /* Now adjust the status vector */
  if (d_after_f)
    for (i = 0; i < s->n; i++)
      if (s->root[i]->status == MPS_ROOT_STATUS_NOT_FLOAT)
        s->root[i]->status = MPS_ROOT_STATUS_CLUSTERED;

  if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)
    mps_dump (s);

  /* == 2 == Perform s->max_pack  packets of Aberth's iterations */
  if (s->DOLOG)
    fprintf (s->logstr, "   DSOLVE: call dpolzero\n");

  for (iter = 0; iter < s->max_pack; iter++)
    {                           /* dloop : DO iter=1,s->max_pack */

       /* mps_dpolzer(s, &nit, &excep);  */
      mps_thread_dpolzer (s, &nit, &excep, required_zeros--);
      it_pack += nit;

      MPS_DEBUG (s, "DPE packet completed in %d iterations", nit);

      if (s->DOLOG)
        fprintf (s->logstr, "Packet %d iterations= %d\n", iter, nit);

      if (excep)
        {
          oldnclust = s->clusterization->n;

          /* cluster analysis */
          if (s->DOLOG)
            fprintf (s->logstr, "   DSOLVE: call dcluster\n");

          mps_dradii (s, s->active_poly, drad);
          mps_dcluster (s, drad, 2 * s->n);   /* Isolation factor */
          if (oldnclust == s->clusterization->n)
            {
              if (s->DOLOG)
                fprintf (s->logstr, "   DSOLVE:  CYCLE\n");
              continue;
            }
          else
            {
              if (s->DOLOG)
                fprintf (s->logstr, "   DSOLVE: call dmodify\n");
              mps_dmodify (s, true);

              if (iter == 0 && !d_after_f)
                for (i = 0; i < s->n; i++)
                  if (s->root[i]->status == MPS_ROOT_STATUS_NEW_CLUSTERED)
                    s->root[i]->status = MPS_ROOT_STATUS_CLUSTERED;

              /* If the polynomial is not given in terms of its
               * coeff. then skip the restart stage */
              if (MPS_IS_MONOMIAL_POLY (s->active_poly))
                {
                  /* choose new starting approximations only for new clusters */
                  if (s->DOLOG)
                    fprintf (s->logstr, "   DSOLVE: call drestart\n");
                  mps_drestart (s);
                }
              /* reset the status vector */
              for (j = 0; j < s->n; j++)
                if (s->root[j]->status == MPS_ROOT_STATUS_NEW_CLUSTERED)
                  s->root[j]->status = MPS_ROOT_STATUS_CLUSTERED;
              for (j = 0; j < s->n; j++)
                s->again_old[j] = s->root[j]->again;

              /* update 'again' */
              if (s->DOLOG)
                fprintf (s->logstr, "   DSOLVE: call update\n");
              mps_update (s);
              /* adjust 'again'
               * This is needed since we are between two packets
               */
              for (i = 0; i < s->n; i++)
                if (!s->again_old[i])
                  s->root[i]->again = false;
              if (s->DOLOG)
                fprintf (s->logstr, "   DSOLVE: call checkstop\n");
              if (mps_check_stop (s))
                goto dsolve_final_cleanup;
            }
        }
      else
        break;
    }
  if (iter == s->max_pack)
    {
      mps_dump (s);
      mps_error (s, 1,
                 "DPE: reached the maximum number of packet iterations");
    }
  /* Otherwise exit since all the approximations are
   * in the root neighbourhood, except for the ones that cannot be
   * represented as double.
   */

  if (s->DOLOG)
    fprintf (s->logstr, "DPE: nit=%d\n", nit);

  /* Update */
  if (s->DOLOG)
    fprintf (s->logstr, "   DSOLVE: now update: call dcluster\n");
  oldnclust = s->clusterization->n;

  mps_dradii (s, s->active_poly, drad);
  mps_dcluster (s, drad, 2 * s->n);   /* Isolation factor */

  if (s->DOLOG)
    fprintf (s->logstr, "   DSOLVE: now call dmodify\n");
  mps_dmodify (s, true);

  /* reset the status vector */
  for (j = 0; j < s->n; j++)
    if (s->root[j]->status == MPS_ROOT_STATUS_NEW_CLUSTERED)
      s->root[j]->status = MPS_ROOT_STATUS_CLUSTERED;

 dsolve_final_cleanup:
  rdpe_vfree (drad);
}

/**
 * @brief Multiprecision version of <code>fsolve()</code>.
 */
void 
mps_msolve (mps_context * s)
{
  int iter, nit, oldnclust, i, j, it_pack, required_zeros = s->n;
  mps_boolean excep;
  int nzc;
  rdpe_t * drad = rdpe_valloc (s->n);

  /* == 1 == Initialize variables */
  it_pack = 0;

  if (s->DOLOG)
    fprintf (s->logstr, "  MSOLVE: call restart\n");
  if (MPS_IS_MONOMIAL_POLY (s->active_poly))
    mps_mrestart (s);
  if (s->DOLOG)
    fprintf (s->logstr, "  MSOLVE: call update1\n");
  mps_update (s);
  if (s->DOLOG)
    {
      fprintf (s->logstr, "  MSOLVE: again after update = ");
      for (i = 0; i < s->n; i++)
        fprintf (s->logstr, "%d", s->root[i]->again);
      fprintf (s->logstr, "\n");
    }

  for (i = 0; i < s->n; i++)
    if (s->root[i]->again)
      s->root[i]->wp = s->mpwp;

  if (s->DOLOG)
    fprintf (s->logstr, "  MSOLVE: call checkstop\n");
  if (mps_check_stop (s))
    {
      oldnclust = s->clusterization->n;

      mps_mmodify (s, true);

      /* reset the status vector */
      for (j = 0; j < s->n; j++)
        if (s->root[j]->status == MPS_ROOT_STATUS_NEW_CLUSTERED)
          s->root[j]->status = MPS_ROOT_STATUS_CLUSTERED;

      goto msolve_final_cleanup;
    }
  nzc = 0;
  if (s->output_config->search_set == MPS_SEARCH_SET_COMPLEX_PLANE && 
      !s->output_config->multiplicity && 
      !s->output_config->root_properties)
    for (i = 0; i < s->n; i++)
      if (MPS_ROOT_STATUS_IS_COMPUTED (s->root[i]->status))
        nzc++;
  if (s->DOLOG)
    fprintf (s->logstr, "  MSOLVE: nzc=%d\n", nzc);

  if (nzc == s->n)
    {
      if (s->DOLOG)
        fprintf (s->logstr, "  MSOLVE: call mmodify and return\n");
      mps_mmodify (s, true);

      /* reset the status vector */
      for (j = 0; j < s->n; j++)
        if (s->root[j]->status == MPS_ROOT_STATUS_NEW_CLUSTERED)
          s->root[j]->status = MPS_ROOT_STATUS_CLUSTERED;

      goto msolve_final_cleanup;
    }

  /* Perform s->max_pack  packets of Aberth's iterations */
  if (s->DOLOG)
    fprintf (s->logstr, "  MSOLVE: Perform packets of Aberth\n");

  for (iter = 0; iter < s->max_pack; iter++)
    {                           /* mloop : DO iter=1,s->max_pack */
      if (s->DOLOG)
        {
          fprintf (s->logstr, "  MSOLVE: packet= %d\n", iter);
          fprintf (s->logstr, "  MSOLVE: again before mpolzer =");
          for (i = 0; i < s->n; i++)
            fprintf (s->logstr, "%d", s->root[i]->again);
          fprintf (s->logstr, "\n");
          fprintf (s->logstr, "  MSOLVE: call mpolzer\n");
        }
      /* mps_mpolzer(s, &nit, &excep); */
      mps_thread_mpolzer (s, &nit, &excep, required_zeros--);

      if (s->debug_level & MPS_DEBUG_APPROXIMATIONS)
        mps_dump (s);

        if (s->DOLOG)
          fprintf (s->logstr, "  MSOLVE: Packet %d: iterations= %d\n", iter, nit);

      it_pack += nit;
      nzc = 0;
      if (s->output_config->search_set == MPS_SEARCH_SET_COMPLEX_PLANE && 
          !s->output_config->multiplicity &&
          !s->output_config->root_properties)  /* DARIO APRILE 98 */
        for (i = 0; i < s->n; i++)
          if (MPS_ROOT_STATUS_IS_COMPUTED (s->root[i]->status))
            nzc++;
      if (s->DOLOG)
        fprintf (s->logstr, "  MSOLVE: check again nzc=%d\n", nzc);
      if (nzc == s->n)
        {
          if (s->DOLOG)
            fprintf (s->logstr, "  MSOLVE: call mmodify and return\n");
          mps_mmodify (s, true);

          /* reset the status vector */
          for (j = 0; j < s->n; j++)
            if (s->root[j]->status == MPS_ROOT_STATUS_NEW_CLUSTERED)
              s->root[j]->status = MPS_ROOT_STATUS_CLUSTERED;
          if (s->DOLOG)
            {
              MPS_DEBUG (s, "Dumping status");
              for (j = 0; j < s->n; j++)
                MPS_DEBUG (s, " %4d: %s", j, 
                           MPS_ROOT_STATUS_TO_STRING (s->root[j]->status));
            }

          goto msolve_final_cleanup;
        }

      if (s->DOLOG)
        fprintf (s->logstr, "  MSOLVE: isolated %d roots excep=%d\n", nzc,
                 excep);

      if (excep)
        {
          oldnclust = s->clusterization->n;

          /* cluster analysis */
          if (s->DOLOG)
            fprintf (s->logstr, "  MSOLVE: call mcluster\n");

          mps_mradii (s, s->active_poly, drad);
          mps_mcluster (s, drad, 2 * s->n);   /* Isolation factor */

          s->newtis_old = s->newtis;
          if (s->newtis == 0)
            mps_mnewtis (s);
          if (s->DOLOG)
            fprintf (s->logstr,
                     "  MSOLVE: newtis_old=%d, newtis=%d, oldncl=%d, s->nclust=%ld\n",
                     s->newtis_old, s->newtis, oldnclust, s->clusterization->n);

          if (oldnclust == s->clusterization->n && !(s->newtis == 1 && s->newtis_old == 0))
            /*#             if(&& iter != 0) AGO99 */
            {
              /*#D !newtis */
              if (s->DOLOG)
                fprintf (s->logstr, "  MSOLVE: CYCLE\n");
              continue;
            }
          else
            {
              if (s->DOLOG)
                fprintf (s->logstr, "  MSOLVE: call modify\n");
              mps_mmodify (s, true);
              if (iter == 0)
                /* if first packet: reset the status vector */
                for (j = 0; j < s->n; j++)
                  if (s->root[j]->status == MPS_ROOT_STATUS_NEW_CLUSTERED)
                    s->root[j]->status = MPS_ROOT_STATUS_CLUSTERED;
              if (s->DOLOG)
                {
                  MPS_DEBUG (s, "Dumping status");
                  for (j = 0; j < s->n; j++)
                    MPS_DEBUG (s, " %4d: %s", j, 
                               MPS_ROOT_STATUS_TO_STRING (s->root[j]->status));
                }
              /* If the polynomial is not given in terms of its coeff. then
               * skip the restart stage */
              if (MPS_IS_MONOMIAL_POLY (s->active_poly))
                {
                  /* choose new starting approximations only for new clusters */
                  if (s->DOLOG)
                    fprintf (s->logstr,
                             "  MSOLVE: call mrestart for new clusters\n");
                  mps_mrestart (s);
                }
              /* reset the s->status vector */
              for (j = 0; j < s->n; j++)
                {
                  if (s->root[j]->status == MPS_ROOT_STATUS_NEW_CLUSTERED)
                    s->root[j]->status = MPS_ROOT_STATUS_CLUSTERED;
                  s->again_old[j] = s->root[j]->again;
                }
              /* update 'again' */
              if (s->DOLOG)
                fprintf (s->logstr, "  MSOLVE: call update2 : ");
              mps_update (s);
              if (s->DOLOG)
                {
                  fprintf (s->logstr, "  MSOLVE: again = ");
                  for (j = 0; j < s->n; j++)
                    fprintf (s->logstr, "%d", s->root[j]->again);
                  fprintf (s->logstr, "\n");
                }
              /* adjust 'again'  This is needed since we are between two packets */
              for (i = 0; i < s->n; i++)
                if (!s->again_old[i])
                  s->root[i]->again = false;
              if (s->DOLOG)
                {
                  fprintf (s->logstr, "  MSOLVE: adjusted again = ");
                  for (j = 0; j < s->n; j++)
                    fprintf (s->logstr, "%d", s->root[j]->again);
                  fprintf (s->logstr, "\n");
                }
              if (s->DOLOG)
                fprintf (s->logstr, "  MSOLVE: call checkstop\n");
              if (mps_check_stop (s))
                {
                  mps_mmodify (s, true);

                  /* reset the s->status vector */
                  for (j = 0; j < s->n; j++)
                    if (s->root[j]->status == MPS_ROOT_STATUS_NEW_CLUSTERED)
                      s->root[j]->status = MPS_ROOT_STATUS_CLUSTERED;

                  goto msolve_final_cleanup;
                }

              nzc = 0;
              if (s->output_config->search_set == MPS_SEARCH_SET_COMPLEX_PLANE && 
                  s->output_config->multiplicity &&
                  !s->output_config->root_properties)
                for (i = 0; i < s->n; i++)
                  if (MPS_ROOT_STATUS_IS_COMPUTED (s->root[i]->status))
                    nzc++;
              if (s->DOLOG)
                fprintf (s->logstr, "  MSOLVE: check again nzc=%d\n", nzc);
              if (nzc == s->n)
                {
                  if (s->DOLOG)
                    fprintf (s->logstr, "  MSOLVE: call mmodify and return");
                  mps_mmodify (s, true);

                  /* reset the s->status vector */
                  for (j = 0; j < s->n; j++)
                    if (s->root[j]->status == MPS_ROOT_STATUS_NEW_CLUSTERED)
                      s->root[j]->status = MPS_ROOT_STATUS_CLUSTERED;

                  goto msolve_final_cleanup;
                }
            }
        }
      else
        break;
    }

  if (iter == s->max_pack)
    {
      mps_error (s, 1, "MP: reached the maximum number of packet iteration");
    }

  if (s->DOLOG)
    {
      fprintf (s->logstr, "  MSOLVE: MP: nit= %d\n", nit);
      fprintf (s->logstr, "  MSOLVE: call mcluster\n");
    }

  mps_mradii (s, s->active_poly, drad);
  mps_mcluster (s, drad, 2 * s->n);   /* Isolation factor */

  if (s->DOLOG)
    fprintf (s->logstr, "  MSOLVE:  call mmodify\n");

  oldnclust = s->clusterization->n;

  mps_mmodify (s, true);

  for (j = 0; j < s->n; j++)
    if (s->root[j]->status == MPS_ROOT_STATUS_NEW_CLUSTERED)
      s->root[j]->status = MPS_ROOT_STATUS_CLUSTERED;

  if (s->DOLOG)
    {
      MPS_DEBUG (s, "Dumping status: ");
      for (j = 0; j < s->n; j++)
        MPS_DEBUG (s, " %4d: %s", j,
                   MPS_ROOT_STATUS_TO_STRING (s->root[j]->status));
    }
  mps_update (s);

  if (s->DOLOG)
    {
      for (j = 0; j < s->n; j++)
        fprintf (s->logstr, "%d", s->root[j]->again);
      fprintf (s->logstr, "\n");
    }

 msolve_final_cleanup:
  rdpe_vfree (drad);
}

/**
 * @brief Multiprecision versione of <code>fpolzer()</code>.
 */
void
mps_mpolzer (mps_context * s, int *it, mps_boolean * excep)
{
  int nzeros, i, iter, l;
  mpc_t corr, abcorr;
  rdpe_t eps, rad1, rtmp;
  cdpe_t ctmp;
  mps_monomial_poly * p = MPS_MONOMIAL_POLY (s->active_poly);

  mpc_init2 (abcorr, s->mpwp);
  mpc_init2 (corr, s->mpwp);

  rdpe_mul_d (eps, s->mp_epsilon, (double) 4 * s->n);

  /* initialize the iteration counter */
  *it = 0;
  *excep = false;

  /* count the number of approximations in the root neighbourhood */
  nzeros = 0;
  for (i = 0; i < s->n; i++)
    if (!s->root[i]->again)
      nzeros++;
  if (nzeros == s->n)
    goto endfun;

  mps_cluster_item * c_item;
  mps_cluster * cluster;
  mps_root * root;

  /* Start Aberth's iterations */
  for (iter = 0; iter < s->max_it; iter++)
    {                           /* do_iter: */
      for (c_item = s->clusterization->first; c_item != NULL; c_item = c_item->next)
        {                       /* do_clust: */
          cluster = c_item->cluster;
          for (root = cluster->first; root != NULL; root = root->next)
            {                   /* do_indice: */
              l = root->k;

              MPS_DEBUG (s, "Iterating on root %d, iter %d", l, iter);

              if (s->root[l]->again)
                {
                  (*it)++;
                      /* sparse/dense polynomial */
                      rdpe_set (rad1, s->root[l]->drad);
                      mps_polynomial_mnewton (s, MPS_POLYNOMIAL (p), s->root[l], corr);
                      if (iter == 0 && !s->root[l]->again && rdpe_gt (s->root[l]->drad,
                                                                rad1)
                          && rdpe_ne (rad1, rdpe_zero))
                        rdpe_set (s->root[l]->drad, rad1);

                                                /************************************************
                                                 The above condition is needed to cope with the case
                                                 where at the first iteration the starting point is
                                                 already in the root neighbourhood and the actually
                                                 computed radius is too big since the value of the
                                                 first derivative is too small.
                                                 In this case the previous radius bound, obtained by
                                                 means of Rouche' is more reliable and strict
                                                 ***********************************************/

                  if (s->root[l]->again ||
                      /* the correction is performed only if iter!=1 or rad[l]!=rad1 */
                      iter != 0
                      || rdpe_ne (s->root[l]->drad, rad1))
                    {
                      mps_maberth_s (s, s->root[l], cluster, abcorr);
                      mpc_mul_eq (abcorr, corr);
                      mpc_neg_eq (abcorr);
                      mpc_add_eq_ui (abcorr, 1, 0);
                      mpc_div (abcorr, corr, abcorr);
                      mpc_sub_eq (s->root[l]->mvalue, abcorr);
                      mpc_get_cdpe (ctmp, abcorr);
                      cdpe_mod (rtmp, ctmp);
                      rdpe_add_eq (s->root[l]->drad, rtmp);
                    }

                  /* check for new approximated roots */
                  if (!s->root[l]->again)
                    {
                      nzeros++;
                      if (nzeros == s->n)
                        goto endfun;
                    }

                  MPS_DEBUG_MPC (s, 15, s->root[l]->mvalue, "s->mroot[%d]", l);
                  MPS_DEBUG_RDPE (s, s->root[l]->drad, "s->drad[%d]", l);

                }
            }
        }
      if (nzeros == s->n)
        goto endfun;
    }
  *excep = true;

endfun:                        /* free local MP variables */
  mpc_clear (corr);
  mpc_clear (abcorr);
}
