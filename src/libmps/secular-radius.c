/*
 * This file is part of MPSolve 3.0
 *
 * Copyright (C) 2001-2012, Dipartimento di Matematica "L. Tonelli", Pisa.
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
 * @brief Compute the radius of inclusions for the roots using Gerschgorin
 * to perform cluster analysis. This routine compute Gerschgorin with the
 * implicit secular representation.
 *
 * A Gerschgorin radius shall be computed for every root and set
 * in <code>s->root[i]->frad</code>, where <code>i</code> is the index of
 * the considered component.
 *
 * @param s The <code>mps_context</code> of the computation.
 * @param fradii A vector where fradii will be stored. 
 */
void
mps_secular_fradii (mps_context * s, double * fradii)
{
  MPS_DEBUG_THIS_CALL;

  mps_secular_equation * sec = s->secular_equation;
  cplx_t sec_ev, diff;
  double prod_b, error;
  int i, j;

  for (i = 0; i < s->n; ++i)
    {
      prod_b = 1.0;

      /* If we have that the root is isolated we can simply ignore it, performing
       * a sort of cluster analysis deflation. */
      if (MPS_ROOT_STATUS_IS_COMPUTED (s, i))
	{
	  fradii[i] = DBL_MAX;
	  continue;
	}

      /* Evaluate the secular equation on root i */
      mps_secular_feval_with_error (s, sec, s->root[i]->fvalue, sec_ev, &error);
      fradii[i] = cplx_mod (sec_ev) + error;

      if (isnan (fradii[i]))
	{
	  fradii[i] = DBL_MAX;
	  continue;
	}

      /* Compute the product of (x - b_i) and p(x) = S(x) * prod(x_b_i) */
      for (j = 0; j < s->n; j++)
	{
	  double adiff;
	  cplx_sub (diff, s->root[i]->fvalue, sec->bfpc[j]);
	  adiff = cplx_mod (diff);

	  prod_b *= adiff;

	  if (i == j) 
	    continue;

	  cplx_sub (diff, s->root[i]->fvalue, s->root[j]->fvalue);
	  prod_b /= cplx_mod (diff);
	}

      /* MPS_DEBUG_CPLX (s, sec_ev, "sec_ev at %d", i); */
      /* MPS_DEBUG (s, "prod_b at %d = %e", i, prod_b); */

      fradii[i] *= prod_b * s->n;
      fradii[i] += cplx_mod (s->root[i]->fvalue) * 4.0 * DBL_EPSILON;

      if (s->root_status[i] == MPS_ROOT_STATUS_ISOLATED && 
	  (fradii[i] < s->root[i]->frad))
	s->root[i]->frad = fradii[i];
    }
}


/**
 * @brief Compute the radius of inclusions for the roots using Gerschgorin
 * to perform cluster analysis. This routine compute Gerschgorin with the
 * implicit secular representation.
 *
 * A Gerschgorin radius shall be computed for every root and set
 * in <code>s->root[i]->frad</code>, where <code>i</code> is the index of
 * the considered component.
 *
 * @param s The <code>mps_context</code> of the computation.
 * @param dradii A vector where fradii will be stored. 
 */
void
mps_secular_dradii (mps_context * s, rdpe_t * dradii)
{
  MPS_DEBUG_THIS_CALL;

  mps_secular_equation * sec = s->secular_equation;
  cdpe_t sec_ev, diff;
  rdpe_t prod_b, rtmp, error;
  int i, j;

  for (i = 0; i < s->n; ++i)
    {
      rdpe_set (prod_b, rdpe_one);

      /* If we have that the root is isolated we can simply ignore it, performing
       * a sort of cluster analysis deflation. */
      if (MPS_ROOT_STATUS_IS_COMPUTED (s, i))
	{
	  rdpe_set (dradii[i], RDPE_MAX);
	  continue;
	}

      /* Evaluate the secular equation on root i */
      mps_secular_deval_with_error (s, sec, s->root[i]->dvalue, sec_ev, error);
      cdpe_mod (dradii[i], sec_ev);
      rdpe_add_eq (dradii[i], error);

      if (isnan (dradii[i]->m))
	{
	  rdpe_set (dradii[i], RDPE_MAX);
	  continue;
	}

      /* Compute the product of (x - b_i) and p(x) = S(x) * prod(x_b_i) */
      for (j = 0; j < s->n; j++)
	{
	  cdpe_sub (diff, s->root[i]->dvalue, sec->bdpc[j]);
	  cdpe_mod (rtmp, diff);
	  rdpe_mul_eq (prod_b, rtmp);

	  if (i == j)
	    continue;

	  cdpe_sub (diff, s->root[i]->dvalue, s->root[j]->dvalue);
	  cdpe_mod (rtmp, diff);
	  rdpe_div_eq (prod_b, rtmp);

	}

      rdpe_mul_eq (dradii[i], prod_b);
      rdpe_mul_eq_d (dradii[i], s->n);

      cdpe_mod (rtmp, s->root[i]->dvalue);
      rdpe_mul_eq_d (rtmp, s->n * 4.0 * DBL_EPSILON);
      rdpe_add_eq (dradii[i], rtmp);

      if (s->root_status[i] == MPS_ROOT_STATUS_ISOLATED && 
	  (rdpe_lt (dradii[i], s->root[i]->drad)))
	rdpe_set (s->root[i]->drad, dradii[i]);
    }
}


/**
 * @brief Compute the radius of inclusions for the roots using Gerschgorin
 * to perform cluster analysis. This routine compute Gerschgorin with the
 * implicit secular representation.
 *
 * A Gerschgorin radius shall be computed for every root and set
 * in <code>s->root[i]->frad</code>, where <code>i</code> is the index of
 * the considered component.
 *
 * @param s The <code>mps_context</code> of the computation.
 * @param dradii A vector where fradii will be stored. 
 */
void
mps_secular_mradii (mps_context * s, rdpe_t * dradii)
{
  MPS_DEBUG_THIS_CALL;

  mps_secular_equation * sec = s->secular_equation;
  mpc_t mdiff, msec_ev, mprod_b;
  cdpe_t diff;
  rdpe_t rtmp, error;
  int i, j;

  mpf_t prod_b, ftmp, merror;

  if (s->lastphase != mp_phase)
    {
      for (i = 0; i < s->n; i++) 
	mpc_set_cdpe (s->root[i]->mvalue, s->root[i]->dvalue);
    }

  mpf_init2 (prod_b, s->mpwp);
  mpf_init2 (ftmp, s->mpwp);
  mpf_init2 (merror, s->mpwp);

  mpc_init2 (mdiff, s->mpwp);
  mpc_init2 (msec_ev, s->mpwp);
  mpc_init2 (mprod_b, s->mpwp);

  for (i = 0; i < s->n; ++i)
    {
      mpf_set_ui (prod_b, 1U);

      /* If we have that the root is isolated we can simply ignore it, performing
       * a sort of cluster analysis deflation. */
      if (MPS_ROOT_STATUS_IS_COMPUTED (s, i))
	{
	  rdpe_set (dradii[i], RDPE_MAX);
	  continue;
	}

      mpc_sub (mdiff, s->root[i]->mvalue, sec->bmpc[i]);
      if (mpc_eq_zero (mdiff))
	{
	  rdpe_set (dradii[i], s->root[i]->drad);
	  continue;
	}

      /* Evaluate the secular equation on root i */
      mps_secular_meval_with_error (s, sec, s->root[i]->mvalue, msec_ev, error);
      mpf_set_rdpe (merror, error);
      
      mpc_mod (ftmp, msec_ev);
      mpf_add_eq (ftmp, merror);
      mpf_get_rdpe (dradii[i], ftmp);
      rdpe_mul_eq_d (dradii[i], s->n);

      MPS_DEBUG_MPC (s, 25, msec_ev, "msec_ev"); 
      MPS_DEBUG_RDPE (s, dradii[i], "S(x)");
       
      if (isnan (dradii[i]->m))
	{
	  rdpe_set (dradii[i], RDPE_MAX);
	  continue;
	}

      /* Compute the product of (x - b_i) and p(x) = S(x) * prod(x_b_i) */
      mpc_set_ui (mprod_b, 1U, 0U);
      for (j = 0; j < s->n; j++)
	{
	  mpc_sub (mdiff, s->root[i]->mvalue, sec->bmpc[j]);
	  mpc_mul_eq (mprod_b, mdiff);

	  if (i == j)
	    continue;

	  mpc_sub (mdiff, s->root[i]->mvalue, s->root[j]->mvalue);
	  mpc_div_eq (mprod_b, mdiff);
	}

      mpc_get_cdpe (diff, mprod_b);
      mpc_mod (prod_b, mprod_b);

      mpf_get_rdpe (rtmp, prod_b);
      rdpe_mul_eq (dradii[i], rtmp);

      mpc_get_cdpe (diff, s->root[i]->mvalue);
      cdpe_mod (rtmp, diff);
      rdpe_mul_eq_d (rtmp, 4.0 * s->n);
      rdpe_mul_eq (rtmp, s->mp_epsilon);
      rdpe_add_eq (dradii[i], rtmp);
      
      if (s->root_status[i] == MPS_ROOT_STATUS_ISOLATED && 
	  (rdpe_lt (dradii[i], s->root[i]->drad)))
	rdpe_set (s->root[i]->drad, dradii[i]);

      MPS_DEBUG_RDPE (s, dradii[i], "dradii[%d]", i); 
    }

  mpf_clear (ftmp);
  mpf_clear (merror);
  mpf_clear (prod_b);

  mpc_clear (mdiff);
  mpc_clear (msec_ev);
  mpc_clear (mprod_b);
}
