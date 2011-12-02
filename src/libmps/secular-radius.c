#include <mps/core.h>

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
 */
void
mps_secular_fradii (mps_status * s)
{
  mps_secular_equation * sec = s->secular_equation;
  cplx_t sec_ev, diff;
  double prod_b;
  int i, j;

  for (i = 0; i < s->n; ++i)
    {
      /* Evaluate the secular equation on root i */
      mps_secular_feval (s, sec, s->froot[i], sec_ev);
      s->frad[i] = cplx_mod (sec_ev);

      /* Compute the product of (x - b_i) and p(x) = S(x) * prod(x_b_i) */
      for (j = 0; j < s->n; j++)
	{
	  cplx_sub (diff, s->froot[i], sec->bfpc[j]);
	  prod_b *= cplx_mod (diff);
	}

      s->frad[i] /= prod_b;
    }
}
