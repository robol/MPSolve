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

/**
 * @brief Check if the target set has been reached or not, and update
 * the field s->root_inclusion[i] for every root.
 */
void
mps_fupdate_inclusions (mps_context * s)
{
  mps_cluster_item * cluster_item;
  mps_cluster * cluster;
  mps_root * root;
  int i, nf = 2 * s->n;

  MPS_DEBUG_THIS_CALL;

  /* Scan the inclusion depending on the selected search set. */
  for (cluster_item = s->clusterization->first; cluster_item != NULL;
       cluster_item = cluster_item->next)
    {
      cluster = cluster_item->cluster;

      for (root = cluster->first; root != NULL; root = root->next)
        {
          i = root->k;
          
          /* First check if the root has already recongnized as part of
           * a set (or out of it) and if that's true skip to the next one. */
          if (s->root[i]->inclusion == MPS_ROOT_INCLUSION_UNKNOWN)
            switch (s->output_config->search_set)
              {
              case MPS_SEARCH_SET_COMPLEX_PLANE:
                s->root[i]->inclusion = MPS_ROOT_INCLUSION_IN;
                break;
            
              case MPS_SEARCH_SET_UNITARY_DISC:
                if (!mps_ftouchunit (s, nf, i))
                  s->root[i]->inclusion = (cplx_mod (s->root[i]->fvalue) < 1) ? MPS_ROOT_INCLUSION_IN : 
                    MPS_ROOT_INCLUSION_OUT;
                break;
            
              case MPS_SEARCH_SET_UNITARY_DISC_COMPL:
                if (!mps_ftouchunit (s, nf, i))
                  s->root[i]->inclusion = (cplx_mod (s->root[i]->fvalue) > 1) ? 
                    MPS_ROOT_INCLUSION_IN : MPS_ROOT_INCLUSION_OUT;
                break;

              case MPS_SEARCH_SET_NEGATIVE_REAL_PART:
                if (!mps_ftouchimag (s, nf, i))
                  s->root[i]->inclusion = (cplx_Re (s->root[i]->fvalue) < 0) ? 
                    MPS_ROOT_INCLUSION_IN : MPS_ROOT_INCLUSION_OUT;
                break;

              case MPS_SEARCH_SET_POSITIVE_REAL_PART:
                if (!mps_ftouchimag (s, nf, i))
                  s->root[i]->inclusion = (cplx_Re (s->root[i]->fvalue) > 0) ? 
                    MPS_ROOT_INCLUSION_IN : MPS_ROOT_INCLUSION_OUT;
                break;

              case MPS_SEARCH_SET_NEGATIVE_IMAG_PART:
                if (!mps_ftouchreal (s, nf, i))
                  s->root[i]->inclusion = (cplx_Im (s->root[i]->fvalue) < 0) ? 
                    MPS_ROOT_INCLUSION_IN : MPS_ROOT_INCLUSION_OUT;
                break;

              case MPS_SEARCH_SET_POSITIVE_IMAG_PART:
                if (!mps_ftouchreal (s, nf, i))
                  s->root[i]->inclusion = (cplx_Im (s->root[i]->fvalue) > 0) ?
                    MPS_ROOT_INCLUSION_IN : MPS_ROOT_INCLUSION_OUT;
                break;

              case MPS_SEARCH_SET_REAL:
                /* Quite a particular case since we are requiring detection
                 * of the reality of the roots. */
                if (cluster->n == 1)
                  {
                    if (mps_ftouchreal (s, 1, i))
                      {
                        if (MPS_STRUCTURE_IS_REAL (s->active_poly->structure) ||
                            (log (s->root[i]->frad) < s->sep - s->n * s->lmax_coeff))
                          {
                            s->root[i]->inclusion = MPS_ROOT_INCLUSION_IN;
                            s->root[i]->attrs = MPS_ROOT_ATTRS_REAL;
                          }
                      }
                    else 
                      {
                        s->root[i]->inclusion = MPS_ROOT_INCLUSION_OUT;
                        s->root[i]->attrs = MPS_ROOT_ATTRS_NONE;
                      }   
                  }
                
                break;

              case MPS_SEARCH_SET_IMAG:
                /* The same as the real case, a part from the fact that
                 * we don't support a pure imaginary coefficients polynomial. */
                if (cluster->n == 1)
                  {
                    if (mps_ftouchimag (s, 1, i))
                      {
                        if (log (s->root[i]->frad) < s->sep - s->n * s->lmax_coeff)
                          {
                            s->root[i]->inclusion = MPS_ROOT_INCLUSION_IN;
                            s->root[i]->attrs = MPS_ROOT_ATTRS_IMAG;
                          }
                        else
                          {
                            s->root[i]->inclusion = MPS_ROOT_INCLUSION_OUT;
                            s->root[i]->attrs = MPS_ROOT_ATTRS_NONE;
                          }
                      }
                  }
                break;

              case MPS_SEARCH_SET_CUSTOM:
                break;
              }
        }
    }

  /* Recheck all the clusters and if a cluster with an uncertain root is found reset
   * all the roots in it as uncertaing. */
  for (cluster_item = s->clusterization->first; cluster_item != NULL;
       cluster_item = cluster_item->next)
    {
      cluster = cluster_item->cluster;

      for (root = cluster->first; root != NULL; root = root->next)
        {
          i = root->k;
          if (s->root[i]->inclusion == MPS_ROOT_INCLUSION_UNKNOWN)
            {
              for (root = cluster->first; root != NULL; root = root->next)
                s->root[root->k]->inclusion = MPS_ROOT_INCLUSION_UNKNOWN;
              break;
            }
        }
    }
}

/**
 * @brief Check if the target set has been reached or not, and update
 * the field s->root[i]->inclusion for every root.
 */
void
mps_dupdate_inclusions (mps_context * s)
{
  mps_cluster_item * cluster_item;
  mps_cluster * cluster;
  mps_root * root;
  int i, nf = 2 * s->n;
  rdpe_t mod;

  MPS_DEBUG_THIS_CALL;

  /* Scan the inclusion depending on the selected search set. */
  for (cluster_item = s->clusterization->first; cluster_item != NULL;
       cluster_item = cluster_item->next)
    {
      cluster = cluster_item->cluster;

      for (root = cluster->first; root != NULL; root = root->next)
        {
          i = root->k;
          
          /* First check if the root has already recongnized as part of
           * a set (or out of it) and if that's true skip to the next one. */
          if (s->root[i]->inclusion == MPS_ROOT_INCLUSION_UNKNOWN)
            switch (s->output_config->search_set)
              {
              case MPS_SEARCH_SET_COMPLEX_PLANE:
                s->root[i]->inclusion = MPS_ROOT_INCLUSION_IN;
                break;
            
              case MPS_SEARCH_SET_UNITARY_DISC:
                if (!mps_dtouchunit (s, nf, i))
                  {
                    cdpe_mod (mod, s->root[i]->dvalue);
                    s->root[i]->inclusion = (rdpe_le (mod, rdpe_one)) ? 
                      MPS_ROOT_INCLUSION_IN : MPS_ROOT_INCLUSION_OUT;
                  }
                break;
            
              case MPS_SEARCH_SET_UNITARY_DISC_COMPL:
                if (!mps_dtouchunit (s, nf, i))
                  {
                    cdpe_mod (mod, s->root[i]->dvalue);
                    s->root[i]->inclusion = (rdpe_ge (mod, rdpe_one)) ?
                      MPS_ROOT_INCLUSION_IN : MPS_ROOT_INCLUSION_OUT;
                  }
                break;

              case MPS_SEARCH_SET_NEGATIVE_REAL_PART:
                if (!mps_dtouchimag (s, nf, i))
                  {
                    rdpe_set (mod, cdpe_Re (s->root[i]->dvalue));
                    s->root[i]->inclusion = (rdpe_le (mod, rdpe_zero)) ?
                      MPS_ROOT_INCLUSION_IN : MPS_ROOT_INCLUSION_OUT;
                  }
                break;

              case MPS_SEARCH_SET_POSITIVE_REAL_PART:
                if (!mps_dtouchimag (s, nf, i))
                  {
                    rdpe_set (mod, cdpe_Re (s->root[i]->dvalue));
                    s->root[i]->inclusion = (rdpe_ge (mod, rdpe_zero)) ?
                      MPS_ROOT_INCLUSION_IN : MPS_ROOT_INCLUSION_OUT;
                  }
                break;

              case MPS_SEARCH_SET_NEGATIVE_IMAG_PART:
                if (!mps_dtouchreal (s, nf, i))
                  {
                    rdpe_set (mod, cdpe_Im (s->root[i]->dvalue));
                    s->root[i]->inclusion = (rdpe_le (mod, rdpe_zero)) ?
                      MPS_ROOT_INCLUSION_IN : MPS_ROOT_INCLUSION_OUT;
                  }
                break;

              case MPS_SEARCH_SET_POSITIVE_IMAG_PART:
                {
                  rdpe_set (mod, cdpe_Im (s->root[i]->dvalue));
                  if (!mps_dtouchreal (s, nf, i))
                    s->root[i]->inclusion = (rdpe_ge (mod, rdpe_zero)) ?
                      MPS_ROOT_INCLUSION_IN : MPS_ROOT_INCLUSION_OUT;
                }
                break;

              case MPS_SEARCH_SET_REAL:
                /* Quite a particular case since we are requiring detection
                 * of the reality of the roots. */
                if (cluster->n == 1)
                  {
                    if (mps_dtouchreal (s, 1, i))
                      {
                        if (MPS_STRUCTURE_IS_REAL (s->active_poly->structure) ||
                            (rdpe_log (s->root[i]->drad) < s->sep - s->n * s->lmax_coeff))
                          {
                            s->root[i]->inclusion = MPS_ROOT_INCLUSION_IN;
                            s->root[i]->attrs = MPS_ROOT_ATTRS_REAL;
                          }
                      }
                    else 
                      {
                        s->root[i]->inclusion = MPS_ROOT_INCLUSION_OUT;
                        s->root[i]->attrs = MPS_ROOT_ATTRS_NONE;
                      }   
                  }
                
                break;

              case MPS_SEARCH_SET_IMAG:
                /* The same as the real case, a part from the fact that
                 * we don't support a pure imaginary coefficients polynomial. */
                if (cluster->n == 1)
                  {
                    if (mps_dtouchimag (s, 1, i))
                      {
                        if (rdpe_log (s->root[i]->drad) < s->sep - s->n * s->lmax_coeff)
                          {
                            s->root[i]->inclusion = MPS_ROOT_INCLUSION_IN;
                            s->root[i]->attrs = MPS_ROOT_ATTRS_IMAG;
                          }
                        else
                          {
                            s->root[i]->inclusion = MPS_ROOT_INCLUSION_OUT;
                            s->root[i]->attrs = MPS_ROOT_ATTRS_NONE;
                          }
                      }
                  }
                break;

              case MPS_SEARCH_SET_CUSTOM:
                break;
              }
        }
    }

  /* Recheck all the clusters and if a cluster with an uncertain root is found reset
   * all the roots in it as uncertaing. */
  for (cluster_item = s->clusterization->first; cluster_item != NULL;
       cluster_item = cluster_item->next)
    {
      cluster = cluster_item->cluster;

      for (root = cluster->first; root != NULL; root = root->next)
        {
          i = root->k;
          if (s->root[i]->inclusion == MPS_ROOT_INCLUSION_UNKNOWN)
            {
              for (root = cluster->first; root != NULL; root = root->next)
                s->root[root->k]->inclusion = MPS_ROOT_INCLUSION_UNKNOWN;
              break;
            }
        }
    }
}

/**
 * @brief Check if the target set has been reached or not, and update
 * the field s->root[i]->inclusion for every root.
 */
void
mps_mupdate_inclusions (mps_context * s)
{
  mps_cluster_item * cluster_item;
  mps_cluster * cluster;
  mps_root * root;
  int i, nf = 2 * s->n;
  cdpe_t cmod;
  rdpe_t mod;

  MPS_DEBUG_THIS_CALL;

  /* Scan the inclusion depending on the selected search set. */
  for (cluster_item = s->clusterization->first; cluster_item != NULL;
       cluster_item = cluster_item->next)
    {
      cluster = cluster_item->cluster;

      for (root = cluster->first; root != NULL; root = root->next)
        {
          i = root->k;

          /* Get a CDPE representation of s->root[i]->mvalue */
          mpc_get_cdpe (cmod, s->root[i]->mvalue);
          
          /* First check if the root has already recongnized as part of
           * a set (or out of it) and if that's true skip to the next one. */
          if (s->root[i]->inclusion == MPS_ROOT_INCLUSION_UNKNOWN)
            switch (s->output_config->search_set)
              {
              case MPS_SEARCH_SET_COMPLEX_PLANE:
                s->root[i]->inclusion = MPS_ROOT_INCLUSION_IN;
                break;
            
              case MPS_SEARCH_SET_UNITARY_DISC:
                if (!mps_mtouchunit (s, nf, i))
                  {
                    cdpe_mod (mod, cmod);
                    s->root[i]->inclusion = (rdpe_le (mod, rdpe_one)) ? 
                      MPS_ROOT_INCLUSION_IN : MPS_ROOT_INCLUSION_OUT;
                  }
                break;
            
              case MPS_SEARCH_SET_UNITARY_DISC_COMPL:
                if (!mps_mtouchunit (s, nf, i))
                  {
                    cdpe_mod (mod, cmod);
                    s->root[i]->inclusion = (rdpe_ge (mod, rdpe_one)) ?
                      MPS_ROOT_INCLUSION_IN : MPS_ROOT_INCLUSION_OUT;
                  }
                break;

              case MPS_SEARCH_SET_NEGATIVE_REAL_PART:
                if (!mps_mtouchimag (s, nf, i))
                  {
                    rdpe_set (mod, cdpe_Re (cmod));
                    s->root[i]->inclusion = (rdpe_le (mod, rdpe_zero)) ?
                      MPS_ROOT_INCLUSION_IN : MPS_ROOT_INCLUSION_OUT;
                  }
                break;

              case MPS_SEARCH_SET_POSITIVE_REAL_PART:
                if (!mps_mtouchimag (s, nf, i))
                  {
                    rdpe_set (mod, cdpe_Re (cmod));
                    s->root[i]->inclusion = (rdpe_ge (mod, rdpe_zero)) ?
                      MPS_ROOT_INCLUSION_IN : MPS_ROOT_INCLUSION_OUT;
                  }
                break;

              case MPS_SEARCH_SET_NEGATIVE_IMAG_PART:
                if (!mps_mtouchreal (s, nf, i))
                  {
                    rdpe_set (mod, cdpe_Im (cmod));
                    s->root[i]->inclusion = (rdpe_le (mod, rdpe_zero)) ?
                      MPS_ROOT_INCLUSION_IN : MPS_ROOT_INCLUSION_OUT;
                  }
                break;

              case MPS_SEARCH_SET_POSITIVE_IMAG_PART:
                {
                  rdpe_set (mod, cdpe_Im (cmod));
                  if (!mps_mtouchreal (s, nf, i))
                    s->root[i]->inclusion = (rdpe_ge (mod, rdpe_zero)) ?
                      MPS_ROOT_INCLUSION_IN : MPS_ROOT_INCLUSION_OUT;
                }
                break;

              case MPS_SEARCH_SET_REAL:
                /* Quite a particular case since we are requiring detection
                 * of the reality of the roots. */
                if (cluster->n == 1)
                  {
                    if (mps_mtouchreal (s, 1, i))
                      {
                        if (MPS_STRUCTURE_IS_REAL (s->active_poly->structure) ||
                            (rdpe_log (s->root[i]->drad) < s->sep - s->n * s->lmax_coeff))
                          {
                            s->root[i]->inclusion = MPS_ROOT_INCLUSION_IN;
                            s->root[i]->attrs = MPS_ROOT_ATTRS_REAL;
                          }
                      }
                    else 
                      {
                        s->root[i]->inclusion = MPS_ROOT_INCLUSION_OUT;
                        s->root[i]->attrs = MPS_ROOT_ATTRS_NONE;
                      }
                  }
                
                break;

              case MPS_SEARCH_SET_IMAG:
                /* The same as the real case, a part from the fact that
                 * we don't support a pure imaginary coefficients polynomial. */
                if (cluster->n == 1)
                  {
                    if (mps_mtouchimag (s, 1, i))
                      {
                        if (rdpe_log (s->root[i]->drad) < s->sep - s->n * s->lmax_coeff)
                          {
                            s->root[i]->inclusion = MPS_ROOT_INCLUSION_IN;
                            s->root[i]->attrs = MPS_ROOT_ATTRS_IMAG;
                          }
                        else
                          {
                            s->root[i]->inclusion = MPS_ROOT_INCLUSION_OUT;
                            s->root[i]->attrs = MPS_ROOT_ATTRS_NONE;
                          }
                      }
                  }
                break;

              case MPS_SEARCH_SET_CUSTOM:
                break;
              }
        }
    }

  /* Recheck all the clusters and if a cluster with an uncertain root is found reset
   * all the roots in it as uncertaing. */
  for (cluster_item = s->clusterization->first; cluster_item != NULL;
       cluster_item = cluster_item->next)
    {
      cluster = cluster_item->cluster;

      for (root = cluster->first; root != NULL; root = root->next)
        {
          i = root->k;
          if (s->root[i]->inclusion == MPS_ROOT_INCLUSION_UNKNOWN)
            {
              for (root = cluster->first; root != NULL; root = root->next)
                s->root[root->k]->inclusion = MPS_ROOT_INCLUSION_UNKNOWN;
              break;
            }
        }
    }
}
