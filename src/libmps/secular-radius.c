#include <mps/mps.h>
#include <math.h>

/**
 * @brief Compute the radius of inclusions for the roots using Gerschgorin
 * to perform cluster analysis. This routine compute Gerschgorin with the
 * implicit secular representation.
 *
 * A Gerschgorin radius shall be computed for every root and set
 * in <code>s->frad[i]</code>, where <code>i</code> is the index of
 * the considered component.
 *
 * @param s The <code>mps_status</code> of the computation.
 * @param fradii A vector where fradii will be stored. 
 */
void
mps_secular_fradii (mps_status * s, double * fradii)
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
      mps_secular_feval_with_error (s, sec, s->froot[i], sec_ev, &error);
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
	  cplx_sub (diff, s->froot[i], sec->bfpc[j]);
	  adiff = cplx_mod (diff);

	  prod_b *= adiff;

	  if (i == j) 
	    continue;

	  cplx_sub (diff, s->froot[i], s->froot[j]);
	  prod_b /= cplx_mod (diff);
	}

      /* MPS_DEBUG_CPLX (s, sec_ev, "sec_ev at %d", i); */
      /* MPS_DEBUG (s, "prod_b at %d = %e", i, prod_b); */

      fradii[i] *= prod_b * s->n;
      fradii[i] += cplx_mod (s->froot[i]) * 4.0 * DBL_EPSILON;

      if (s->root_status[i] == MPS_ROOT_STATUS_ISOLATED && 
	  (fradii[i] < s->frad[i]))
	s->frad[i] = fradii[i];
    }
}


/**
 * @brief Compute the radius of inclusions for the roots using Gerschgorin
 * to perform cluster analysis. This routine compute Gerschgorin with the
 * implicit secular representation.
 *
 * A Gerschgorin radius shall be computed for every root and set
 * in <code>s->frad[i]</code>, where <code>i</code> is the index of
 * the considered component.
 *
 * @param s The <code>mps_status</code> of the computation.
 * @param dradii A vector where fradii will be stored. 
 */
void
mps_secular_dradii (mps_status * s, rdpe_t * dradii)
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
      mps_secular_deval_with_error (s, sec, s->droot[i], sec_ev, error);
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
	  cdpe_sub (diff, s->droot[i], sec->bdpc[j]);
	  cdpe_mod (rtmp, diff);
	  rdpe_mul_eq (prod_b, rtmp);

	  if (i == j)
	    continue;

	  cdpe_sub (diff, s->droot[i], s->droot[j]);
	  cdpe_mod (rtmp, diff);
	  rdpe_div_eq (prod_b, rtmp);

	}

      rdpe_mul_eq (dradii[i], prod_b);
      rdpe_mul_eq_d (dradii[i], s->n);

      cdpe_mod (rtmp, s->droot[i]);
      rdpe_mul_eq_d (rtmp, s->n * 4.0 * DBL_EPSILON);
      rdpe_add_eq (dradii[i], rtmp);

      if (s->root_status[i] == MPS_ROOT_STATUS_ISOLATED && 
	  (rdpe_lt (dradii[i], s->drad[i])))
	rdpe_set (s->drad[i], dradii[i]);
    }
}


/**
 * @brief Compute the radius of inclusions for the roots using Gerschgorin
 * to perform cluster analysis. This routine compute Gerschgorin with the
 * implicit secular representation.
 *
 * A Gerschgorin radius shall be computed for every root and set
 * in <code>s->frad[i]</code>, where <code>i</code> is the index of
 * the considered component.
 *
 * @param s The <code>mps_status</code> of the computation.
 * @param dradii A vector where fradii will be stored. 
 */
void
mps_secular_mradii (mps_status * s, rdpe_t * dradii)
{
  MPS_DEBUG_THIS_CALL;

  mps_secular_equation * sec = s->secular_equation;
  mpc_t mdiff, msec_ev;
  cdpe_t sec_ev, diff;
  rdpe_t prod_b, rtmp, error;
  int i, j;

  mpc_init2 (mdiff, s->mpwp);
  mpc_init2 (msec_ev, s->mpwp);

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

      mpc_sub (mdiff, s->mroot[i], sec->bmpc[i]);
      if (mpc_eq_zero (mdiff))
	{
	  rdpe_set (dradii[i], s->drad[i]);
	  continue;
	}

      /* Evaluate the secular equation on root i */
      mpc_set_si (mdiff, -3, 1);
      mps_secular_meval_with_error (s, sec, mdiff, msec_ev, error);
      mps_secular_meval_with_error (s, sec, s->mroot[i], msec_ev, error);
      mpc_get_cdpe (sec_ev, msec_ev);
      cdpe_mod (dradii[i], sec_ev);
      rdpe_add_eq (dradii[i], error);
      rdpe_mul_eq_d (dradii[i], s->n);

      MPS_DEBUG_MPC (s, 25, msec_ev, "msec_ev");
      MPS_DEBUG_RDPE (s, dradii[i], "S(x)");

      if (isnan (dradii[i]->m))
	{
	  rdpe_set (dradii[i], RDPE_MAX);
	  continue;
	}

      /* Compute the product of (x - b_i) and p(x) = S(x) * prod(x_b_i) */
      for (j = 0; j < s->n; j++)
	{
	  mpc_sub (mdiff, s->mroot[i], sec->bmpc[j]);
	  mpc_get_cdpe (diff, mdiff);
	  cdpe_mod (rtmp, diff);
	  rdpe_mul_eq (prod_b, rtmp);

	  if (i == j)
	    continue;

	  mpc_sub (mdiff, s->mroot[i], s->mroot[j]);
	  mpc_get_cdpe (diff, mdiff);
	  cdpe_mod (rtmp, diff);
	  rdpe_div_eq (prod_b, rtmp);
	}

      rdpe_mul_eq (dradii[i], prod_b);

      mpc_get_cdpe (diff, s->mroot[i]);
      cdpe_mod (rtmp, diff);
      rdpe_mul_eq_d (rtmp, 4.0 * s->n);
      rdpe_mul_eq (rtmp, s->mp_epsilon);
      rdpe_add_eq (dradii[i], rtmp);
      
      if (s->root_status[i] == MPS_ROOT_STATUS_ISOLATED && 
	  (rdpe_lt (dradii[i], s->drad[i])))
	rdpe_set (s->drad[i], dradii[i]);

      MPS_DEBUG_RDPE (s, dradii[i], "dradii[%d]", i);
    }

  mpc_clear (mdiff);
  mpc_clear (msec_ev);
}
