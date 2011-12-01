#include <mps/core.h>

/**
 * @brief Compute the radius of inclusions for the roots using Gerschgorin
 * to perform cluster analysis. 
 *
 * A Gerschgorin radius shall be computed for every root and set
 * in <code>s->frad[i]</code>, where <code>i</code> is the index of
 * the considered component.
 *
 * @param s The <code>mps_status</code> of the computation.
 */
void
mps_monomial_fradii (mps_status * s)
{
  MPS_DEBUG_THIS_CALL;

  cplx_t pol;
  double new_rad, relative_error;
  mps_monomial_poly * p = s->monomial_poly;
  int i, j;

  for (i = 0; i < s->n; i++)
    {
      cplx_t diff;
      mps_fhorner_with_error (s, s->monomial_poly, s->froot[i], pol, &relative_error);
      new_rad = cplx_mod (pol) + relative_error + cplx_mod (s->froot[i]) * 4.0 * DBL_EPSILON;

      for (j = 0; j < s->n; j++)
	{
	  if (i == j)
	    continue;

	  cplx_sub (diff, s->froot[i], s->froot[j]);
	      
	  /* Check for floating point exceptions in here */
	  if (cplx_eq_zero (diff))
	    {
	      new_rad = DBL_MAX;
	      break;
	    }

	  new_rad /= cplx_mod (diff);
	}
      new_rad /= p->fap[s->n];
      
      s->frad[i] = new_rad;
    }
}

/**
 * @brief Compute the radius of inclusions for the roots using Gerschgorin
 * to perform cluster analysis. 
 *
 * A Gerschgorin radius shall be computed for every root and set
 * in <code>s->frad[i]</code>, where <code>i</code> is the index of
 * the considered component.
 *
 * @param s The <code>mps_status</code> of the computation.
 */
void
mps_monomial_dradii (mps_status * s)
{
  MPS_DEBUG_THIS_CALL;

  cdpe_t pol;
  rdpe_t new_rad, relative_error, rtmp;
  mps_monomial_poly * p = s->monomial_poly;
  int i, j;

  for (i = 0; i < s->n; i++)
    {
      cdpe_t diff;
      mps_dhorner_with_error (s, s->monomial_poly, s->droot[i], pol, relative_error);

      cdpe_mod (new_rad, pol);
      rdpe_add_eq (new_rad, relative_error);
      cdpe_mod (rtmp, s->droot[i]);
      rdpe_mul_eq_d (rtmp, 4.0 * DBL_EPSILON);
      rdpe_add_eq (new_rad, rtmp);

      for (j = 0; j < s->n; j++)
	{
	  if (i == j)
	    continue;

	  cdpe_sub (diff, s->droot[i], s->droot[j]);
	      
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
      rdpe_set (s->drad[i], new_rad);
    }
}
