/*
 * This file is part of MPSolve 3.1.8
 *
 * Copyright (C) 2001-2019, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Dario Andrea Bini <bini@dm.unipi.it>
 *   Giuseppe Fiorentino <fiorent@dm.unipi.it>
 *   Leonardo Robol <leonardo.robol@sns.it>
 */


#include <mps/mps.h>
#include <math.h>

typedef struct {
  mps_context * ctx;
  mps_polynomial * p;
  double * fradii;
  int i;
} _mps_fradii_worker_data;

static void *
_mps_fradii_worker (void * data_ptr)
{
  _mps_fradii_worker_data * data = (_mps_fradii_worker_data*) data_ptr;
  int i = data->i;
  double * fradii = data->fradii;
  mps_context * s = data->ctx;
  mps_polynomial * p = data->p; 

  cplx_t pol, diff;
  double new_rad, relative_error;
  int j;

  mpc_t lc;
  cplx_t ctmp;

  mpc_init2 (lc, DBL_MANT_DIG);
  mps_polynomial_get_leading_coefficient (s, p, lc);

  /* Compute the value of the polynomial in this point */
  mps_polynomial_feval (s, p, s->root[i]->fvalue, pol, &relative_error);
  
  /* If we got a floating point exception, we need to switch to DPE on this component */
  if (cplx_check_fpe (pol))
    {
      s->root[i]->status = MPS_ROOT_STATUS_NOT_FLOAT;
      fradii[i] = DBL_MAX;
      free (data);
      return NULL;
    }

  new_rad = cplx_mod (pol) + relative_error + cplx_mod (s->root[i]->fvalue) * 4.0 * DBL_EPSILON;
  new_rad *= s->n;
  
  for (j = 0; j < s->n; j++)
    {
      if (i == j)
	continue;
      
      cplx_sub (diff, s->root[i]->fvalue, s->root[j]->fvalue);
      
      /* Check for floating point exceptions in here */
      if (cplx_eq_zero (diff))
	{
	  new_rad = DBL_MAX;
	  break;
	}
      
      new_rad /= cplx_mod (diff);
    }

  /* Check for possible underflows. In case they happen, perform the
   * computation using DPE arithmetic.*/
  if (new_rad == 0.0) 
    {
      rdpe_t radius, rtmp;
      cdpe_t cdtmp;

      new_rad = cplx_mod (pol) + relative_error + cplx_mod (s->root[i]->fvalue) * 4.0 * DBL_EPSILON;
      new_rad *= s->n;
      
      rdpe_set_d (radius, new_rad);

      for (j = 0; j < s->n; j++)
	{
	  if (i == j)
	    continue;
      
	  cplx_sub (diff, s->root[i]->fvalue, s->root[j]->fvalue);
      
	  /* Check for floating point exceptions in here */
	  if (cplx_eq_zero (diff))
	    {
	      rdpe_set (radius, RDPE_MAX);
	      break;
	    }

	  rdpe_div_eq_d (radius, cplx_mod (diff));
	}

      mpc_get_cdpe (cdtmp, lc);
      cdpe_mod     (rtmp,  cdtmp);

      rdpe_div_eq (radius, rtmp);
      new_rad = rdpe_get_d (radius);
    }
  else
    {
      mpc_get_cplx (ctmp, lc);
      new_rad /= cplx_mod (ctmp);
    }
    
  fradii[i] = p->degree * new_rad * 
    (1 + p->degree * DBL_EPSILON * 2.0 * sqrt(2)) + 
    cplx_mod (s->root[i]->fvalue) * DBL_EPSILON * 2.0
    + DBL_MIN;

  mpc_clear (lc);
  free (data);

  return NULL;
}

/**
 * @brief Compute the floating point inclusion radius according to the
 * polynomial representation.
 *
 * @param s The <code>mps_context</code> of the computation.
 * @param p The mps_polynomial whose Gerschgorin radii should be evaluated.
 * @param fradii The array of double where the radii will be stored.
 */
MPS_PRIVATE void
mps_fradii (mps_context * s, mps_polynomial * p, double * fradii)
{
  MPS_DEBUG_THIS_CALL (s);
  int i;

  if (!p->feval)
    {
      for (i = 0; i < s->n; i++)
        fradii[i] = s->root[i]->frad;
      return;
    }

  for (i = 0; i < s->n; i++)
    {
      _mps_fradii_worker_data * data = mps_new (_mps_fradii_worker_data);
      
      data->ctx = s;
      data->p = p;
      data->i = i;
      data->fradii = fradii;

      mps_thread_pool_assign (s, s->pool, _mps_fradii_worker, data);
    }

  mps_thread_pool_wait (s, s->pool);
}

/**
 * @brief Compute the DPE inclusion radius according to the
 * polynomial representation.
 *
 * @param s The <code>mps_context</code> of the computation.
 * @param p The mps_polynomial whose Gerschgorin radii should be evaluated.
 * @param dradii The array of DPE where the radii will be stored.
 */
MPS_PRIVATE void
mps_dradii (mps_context * s, mps_polynomial * p, rdpe_t * dradii)
{
  MPS_DEBUG_THIS_CALL (s);

  cdpe_t pol;
  rdpe_t new_rad, relative_error, rtmp;
  int i, j;

  if (!p->deval)
    {
      for (i = 0; i < s->n; i++)
        rdpe_set (dradii[i], s->root[i]->drad);
      return;
    }

  for (i = 0; i < s->n; i++)
    {
      cdpe_t diff;
      mps_polynomial_deval (s, p, s->root[i]->dvalue, pol, relative_error);

      cdpe_mod (new_rad, pol);
      rdpe_add_eq (new_rad, relative_error);
      cdpe_mod (rtmp, s->root[i]->dvalue);
      rdpe_mul_eq_d (rtmp, 4.0 * DBL_EPSILON);
      rdpe_add_eq (new_rad, rtmp);
      rdpe_mul_eq_d (new_rad, s->n);

      for (j = 0; j < s->n; j++)
        {
          if (i == j)
            continue;

          cdpe_sub (diff, s->root[i]->dvalue, s->root[j]->dvalue);

          /* Check for floating point exceptions in here */
          if (cdpe_eq_zero (diff))
            {
              rdpe_set (new_rad, RDPE_MAX);
              break;
            }

          cdpe_mod (rtmp, diff);
          rdpe_div_eq (new_rad, rtmp);
        }

      {
        mpc_t lc;
        mpc_init2 (lc, DBL_MANT_DIG);
        mps_polynomial_get_leading_coefficient (s, p, lc);
        mpc_rmod (rtmp, lc);
        rdpe_div_eq (new_rad, rtmp);
        mpc_clear (lc);
      }

      rdpe_set (dradii[i], new_rad);
    }
}

/**
 * @brief Compute the Multiprecision inclusion radius according to the
 * polynomial representation.
 *
 * @param s The <code>mps_context</code> of the computation.
 * @param p The mps_polynomial whose Gerschgorin radii should be evaluated.
 * @param dradii The array of DPE where the radii will be stored.
 */
MPS_PRIVATE void
mps_mradii (mps_context * s, mps_polynomial * p, rdpe_t * dradii)
{
  MPS_DEBUG_THIS_CALL (s);

  mpc_t pol, mdiff;
  cdpe_t cpol, diff, cdtmp;
  rdpe_t new_rad, relative_error, rtmp;
  int i, j;

  if (!p->meval)
    {
      for (i = 0; i < s->n; i++)
        rdpe_set (dradii[i], s->root[i]->drad);
      return;
    }


  mpc_init2 (pol, s->mpwp);
  mpc_init2 (mdiff, s->mpwp);

  for (i = 0; i < s->n; i++)
    {
      mps_polynomial_meval (s, p, s->root[i]->mvalue, pol, relative_error);

      mpc_get_cdpe (cpol, pol);
      cdpe_mod (new_rad, cpol);
      rdpe_add_eq (new_rad, relative_error);
      mpc_get_cdpe (cdtmp, s->root[i]->mvalue);
      cdpe_mod (rtmp, cdtmp);
      rdpe_mul_eq (rtmp, s->mp_epsilon);
      rdpe_add_eq (new_rad, rtmp);
      rdpe_mul_eq_d (new_rad, s->n);

      rdpe_set (relative_error, rdpe_zero);

      for (j = 0; j < s->n; j++)
        {
          if (i == j)
            continue;

          mpc_sub (mdiff, s->root[i]->mvalue, s->root[j]->mvalue);
          mpc_get_cdpe (diff, mdiff);

          /* Check for floating point exceptions in here */
          if (mpc_eq_zero (mdiff))
            {
              rdpe_set (dradii[i], RDPE_MAX);
              goto mradius_cleanup;
            }

          mpc_rmod (rtmp, mdiff);
          rdpe_div_eq (new_rad, rtmp);
        }

      rdpe_mul_eq_d (new_rad, 1 + 2 * s->n * sqrt (2) * DBL_EPSILON);
      rdpe_mul_eq_d (new_rad, p->degree);

      {
        mpc_t lc;
        mpc_init2 (lc, s->mpwp);
        mps_polynomial_get_leading_coefficient (s, p, lc);
        mpc_rmod (rtmp, lc);
        rdpe_div_eq (new_rad, rtmp);
        mpc_clear (lc);
      }

      rdpe_set (dradii[i], new_rad);
    }

mradius_cleanup:

  mpc_clear (pol);
  mpc_clear (mdiff);
}
