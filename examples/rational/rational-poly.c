#include "rational-poly.h"
#include <string.h>

void 
mps_rational_poly_free (mps_context * ctx, mps_rational_poly * p)
{
  free (p->a);
  free (p->b);
  free (p->c);
  mpcf_clear (p->lc);
  free (p);
}

mps_rational_poly * 
mps_rational_poly_new (mps_context * ctx, cplx_t * a, cplx_t * b,
		       cplx_t * c, int nn, int pp)
{
  int i;
  mps_rational_poly * p = mps_new (mps_rational_poly);
  mps_polynomial_init (ctx, MPS_POLYNOMIAL (p));

  MPS_POLYNOMIAL (p)->type_name = "mps_rational_poly";
  MPS_POLYNOMIAL (p)->density = MPS_DENSITY_DENSE;
  MPS_POLYNOMIAL (p)->structure = MPS_STRUCTURE_REAL_FP;

  MPS_POLYNOMIAL (p)->feval = mps_rational_poly_feval;
  MPS_POLYNOMIAL (p)->deval = mps_rational_poly_deval;
  MPS_POLYNOMIAL (p)->meval = mps_rational_poly_meval;
  MPS_POLYNOMIAL (p)->mnewton = mps_rational_poly_mnewton;
  MPS_POLYNOMIAL (p)->thread_safe = false;
  MPS_POLYNOMIAL (p)->get_leading_coefficient = mps_rational_poly_get_leading_coefficient;
  
  /* TODO: Check if the degree is lower. */
  MPS_POLYNOMIAL (p)->degree = (nn-1) * pp;

  p->a = mps_newv (cplx_t, nn);
  p->b = mps_newv (cplx_t, nn);
  p->c = mps_newv (cplx_t, nn);
  
  memcpy (p->a, a, sizeof (cplx_t) * nn);
  memcpy (p->b, b, sizeof (cplx_t) * nn);
  memcpy (p->c, c, sizeof (cplx_t) * nn);

  p->nn = nn;
  p->pp = pp;

  mpcf_init2 (p->lc, 64);

  /* Compute the leading coefficient */
  cplx_t lc, prod;
  cplx_set (lc, cplx_zero);
  cplx_set (prod, cplx_one);
  for (i = 0; i < nn; i++)
    {
      cplx_t xx;      
      cplx_pow_si (xx, b[i], pp);
      cplx_mul_eq (prod, xx);
      cplx_div (xx, c[i], xx);
      cplx_add_eq (lc, xx);
    }

  cplx_mul_eq (lc, prod);
  mpcf_set_d (p->lc, cplx_Re (lc), cplx_Im (lc));

  return p;
}

void mps_rational_poly_get_leading_coefficient (mps_context * ctx, 
						mps_polynomial * p, 
						mpcf_t lc)
{
  mps_rational_poly * rp = MPS_RATIONAL_POLY (p);
  mpcf_set (lc, rp->lc);
}
