#include <mps/core.h>

/**
 * @brief Evaluate a secular equation <code>sec</code> in the point <code>x</code>
 *
 * @param s The <code>mps_status</code> of the computation.
 * @param sec The secular equation to evaluate.
 * @param x The point in which the secular equation must be evaluated.
 * @param value The value of the secular equation in the pointer <code>x</code>.
 */
void
mps_secular_feval (mps_status * s, mps_secular_equation * sec, cplx_t x, cplx_t value)
{
  cplx_t ctmp;
  int i;
  
  cplx_set (value, cplx_zero);
  
  for (i = 0; i < s->n; i++)
    {
      cplx_sub (ctmp, x, sec->bfpc[i]);
      cplx_div (ctmp, sec->afpc[i], ctmp);
      cplx_add_eq (value, ctmp);
    }

  cplx_sub_eq (value, cplx_one);
}

/**
 * @brief Evaluate a secular equation <code>sec</code> in the point <code>x</code>
 *
 * @param s The <code>mps_status</code> of the computation.
 * @param sec The secular equation to evaluate.
 * @param x The point in which the secular equation must be evaluated.
 * @param value The value of the secular equation in the pointer <code>x</code>.
 */
void
mps_secular_deval (mps_status * s, mps_secular_equation * sec, cdpe_t x, cdpe_t value)
{
  cdpe_t ctmp;
  int i;

  cdpe_set (value, cdpe_zero);

  for (i = 0; i < s->n; i++)
    {
      cdpe_sub (ctmp, x, sec->bdpc[i]);
      cdpe_div (ctmp, sec->adpc[i], ctmp);
      cdpe_add_eq (value, ctmp);
    }

  cdpe_sub_eq (value, cdpe_one);
}
