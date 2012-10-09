#include <mps/mps.h>
#include <math.h>

/**
 * @brief Compute the radius of inclusions for the roots using Gerschgorin
 * to perform cluster analysis. 
 *
 * A Gerschgorin radius shall be computed for every root and set
 * in <code>s->root[i]->frad</code>, where <code>i</code> is the index of
 * the considered component.
 *
 * @param s The <code>mps_context</code> of the computation.
 * @param fradii The array where the computed radii will be stored.
 */
void
mps_monomial_fradii (mps_context * s, double * fradii)
{
  MPS_DEBUG_THIS_CALL;

  cplx_t pol;
  double new_rad, relative_error;
  mps_monomial_poly * p = s->monomial_poly;
  int i, j;

  for (i = 0; i < s->n; i++)
    {
      cplx_t diff;
      
      /* Compute the value of the polynomial in this point */
      mps_fhorner_with_error (s, s->monomial_poly, s->root[i]->fvalue, pol, &relative_error);

      /* If we got a floating point exception, we need to switch to DPE on this component */
      if (cplx_check_fpe (pol))
	{
	  s->root_status[i] = MPS_ROOT_STATUS_NOT_FLOAT;
	  fradii[i] = DBL_MAX;
	  continue;
	}
      new_rad = cplx_mod (pol) + relative_error + cplx_mod (s->root[i]->fvalue) * 4.0 * DBL_EPSILON;

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
      new_rad /= p->fap[s->n];
      
      fradii[i] = new_rad;
    }
}

/**
 * @brief Compute the radius of inclusions for the roots using Gerschgorin
 * to perform cluster analysis. 
 *
 * A Gerschgorin radius shall be computed for every root and set
 * in <code>s->root[i]->frad</code>, where <code>i</code> is the index of
 * the considered component.
 *
 * @param s The <code>mps_context</code> of the computation.
 * @param dradii The array of DPE where the radii will be stored.
 */
void
mps_monomial_dradii (mps_context * s, rdpe_t * dradii)
{
  MPS_DEBUG_THIS_CALL;

  cdpe_t pol;
  rdpe_t new_rad, relative_error, rtmp;
  mps_monomial_poly * p = s->monomial_poly;
  int i, j;

  for (i = 0; i < s->n; i++)
    {
      cdpe_t diff;
      mps_dhorner_with_error (s, s->monomial_poly, s->root[i]->dvalue, pol, relative_error);

      cdpe_mod (new_rad, pol);
      rdpe_add_eq (new_rad, relative_error);
      cdpe_mod (rtmp, s->root[i]->dvalue);
      rdpe_mul_eq_d (rtmp, 4.0 * DBL_EPSILON);
      rdpe_add_eq (new_rad, rtmp);

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

      rdpe_div_eq (new_rad, p->dap[s->n]);
      rdpe_set (dradii[i], new_rad);
    }
}

/**
 * @brief Compute the radius of inclusions for the roots using Gerschgorin
 * to perform cluster analysis. 
 *
 * A Gerschgorin radius shall be computed for every root and set
 * in <code>s->root[i]->frad</code>, where <code>i</code> is the index of
 * the considered component.
 *
 * @param s The <code>mps_context</code> of the computation.
 * @param dradii The array of DPE where the radii will be stored.
 */
void
mps_monomial_mradii (mps_context * s, rdpe_t * dradii)
{
  MPS_DEBUG_THIS_CALL;

  mpc_t pol, mdiff;
  cdpe_t cpol, diff, cdtmp;
  rdpe_t new_rad, relative_error, rtmp;
  mps_monomial_poly * p = s->monomial_poly;
  int i, j;

  mpc_init2 (pol, s->mpwp);
  mpc_init2 (mdiff, s->mpwp);

  for (i = 0; i < s->n; i++)
    {
      mps_mhorner_with_error2 (s, s->monomial_poly, s->root[i]->mvalue, pol, relative_error, s->mpwp);

      mpc_get_cdpe (cpol, pol);
      cdpe_mod (new_rad, cpol);
      rdpe_add_eq (new_rad, relative_error);
      mpc_get_cdpe (cdtmp, s->root[i]->mvalue);
      cdpe_mod (rtmp, cdtmp);
      rdpe_mul_eq (rtmp, s->mp_epsilon);
      rdpe_add_eq (new_rad, rtmp);

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

      rdpe_mul_eq_d (new_rad, 1 + 2 * s->n * sqrt(2) * DBL_EPSILON);
      rdpe_mul_eq_d (new_rad, p->n);
      rdpe_div_eq (new_rad, p->dap[s->n]);
      rdpe_set (dradii[i], new_rad);
    }

 mradius_cleanup:

  mpc_clear (pol);
  mpc_clear (mdiff);
}


