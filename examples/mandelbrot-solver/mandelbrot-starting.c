#include "mandelbrot-poly.h"

void 
mps_mandelbrot_poly_dstart (mps_context *ctx, mps_polynomial *p)
{
  mps_mandelbrot_poly *mp = MPS_MANDELBROT_POLY (p);

  if (mp->level <= 4)
    {
      mps_general_dstart (ctx, p);
    }
  else
    {
      int i;
      cplx_t *roots = NULL;
      cdpe_t rot;

      cdpe_set_d (rot, 1 + p->degree * DBL_EPSILON, DBL_EPSILON);

      mps_context * new_ctx = mps_context_new ();
      mps_mandelbrot_poly *new_mp = mps_mandelbrot_poly_new (new_ctx, mp->level - 1);
      
      mps_context_set_input_poly (new_ctx, MPS_POLYNOMIAL (new_mp));
      mps_context_select_algorithm (new_ctx, MPS_ALGORITHM_SECULAR_GA);
      mps_context_set_starting_phase (new_ctx, float_phase);
      mps_mpsolve (new_ctx);

      mps_context_get_roots_d (new_ctx, &roots, NULL);

      for (i = 0; i < ctx->n; i++)
	{
	  int j = i / 2;

	  cdpe_set_x (ctx->root[i]->dvalue, roots[j]);

	  if (i % 2 == 0)
	    cdpe_mul_eq (ctx->root[i]->dvalue, rot);
	}

      cdpe_set (ctx->root[ctx->n - 1]->dvalue, cdpe_one);

      cplx_vfree (roots);
      mps_polynomial_free (new_ctx, MPS_POLYNOMIAL (new_mp));
      mps_context_free (new_ctx);
    }
}

void
mps_mandelbrot_poly_fstart (mps_context *ctx, mps_polynomial *p)
{
  mps_mandelbrot_poly *mp = MPS_MANDELBROT_POLY (p);

  if (mp->level <= 5)
    {
      mps_general_fstart (ctx, p);
    }
  else
    {
      int i;
      cplx_t *roots = cplx_valloc (p->degree / 2);
      cplx_t rot;

      cplx_set_d (rot, 1 + p->degree * DBL_EPSILON, 0.0);

      mps_context * new_ctx = mps_context_new ();
      mps_mandelbrot_poly *new_mp = mps_mandelbrot_poly_new (new_ctx, mp->level - 1);
      
      mps_context_set_input_poly (new_ctx, MPS_POLYNOMIAL (new_mp));
      mps_context_select_algorithm (new_ctx, MPS_ALGORITHM_SECULAR_GA);
      mps_mpsolve (new_ctx);

      mps_context_get_roots_d (new_ctx, &roots, NULL);

      for (i = 0; i < ctx->n; i++)
	{
	  int j = i / 2;
	  cplx_set (ctx->root[i]->fvalue, roots[j]);

	  if (i % 2 == 0)
	    cplx_mul_eq (ctx->root[i]->fvalue, rot);

	  ctx->root[i]->frad = DBL_MAX;
	}

      cplx_set (ctx->root[ctx->n - 1]->fvalue, cplx_one);

      cplx_vfree (roots);
      mps_polynomial_free (new_ctx, MPS_POLYNOMIAL (new_mp));
      mps_context_free (new_ctx);
    }
}
