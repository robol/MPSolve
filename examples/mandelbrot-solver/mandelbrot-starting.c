/*
 * This file is part of MPSolve 3.0
 *
 * Copyright (C) 2001-2013, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors: 
 *   Leonardo Robol <robol@mail.dm.unipi.it>
 */

#include "mandelbrot-poly.h"

void 
mps_mandelbrot_poly_dstart (mps_context *ctx, mps_polynomial *p, mps_approximation ** approximations)
{
  mps_mandelbrot_poly *mp = MPS_MANDELBROT_POLY (p);
  int n = mps_context_get_degree (ctx);

  if (mp->level <= 4)
    {
      mps_general_dstart (ctx, p, approximations);
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

      for (i = 0; i < n - 1; i++)
	{
	  cdpe_t starting_point;
	  int j = i / 2;

	  cdpe_set_x (starting_point, roots[j]);

	  if (i % 2 == 0)
	    cdpe_mul_eq (starting_point, rot);

	  mps_approximation_set_dvalue (ctx, approximations[i], starting_point);
	}

      mps_approximation_set_dvalue (ctx, approximations[n - 1], cdpe_one);

      cplx_vfree (roots);
      mps_mandelbrot_poly_free (new_ctx, MPS_POLYNOMIAL (new_mp));
      mps_context_free (new_ctx);
    }
}

void
mps_mandelbrot_poly_fstart (mps_context *ctx, mps_polynomial *p, mps_approximation ** approximations)
{
  mps_mandelbrot_poly *mp = MPS_MANDELBROT_POLY (p);
  int n = mps_context_get_degree (ctx);

  if (mp->level <= 7)
    {
      mps_general_fstart (ctx, p, approximations);
    }
  else
    {
      int i;
      cplx_t *roots = cplx_valloc (p->degree / 2);
      cplx_t rot_up, rot_down;

      int eps = (mp->level > DBL_DIG) ? mp->level : DBL_DIG + 1;

      cplx_set_d (rot_up, 1 + pow (.5, eps), -4.0 * DBL_EPSILON);
      cplx_set_d (rot_down, 1 + pow (.5, eps), 4.0 * DBL_EPSILON);

      mps_context * new_ctx = mps_context_new ();
      mps_mandelbrot_poly *new_mp = mps_mandelbrot_poly_new (new_ctx, mp->level - 1);

      mps_context_set_input_poly (new_ctx, MPS_POLYNOMIAL (new_mp));
      mps_context_select_algorithm (new_ctx, MPS_ALGORITHM_SECULAR_GA);
      mps_context_set_output_prec (new_ctx, 16);
      mps_mpsolve (new_ctx);

      mps_context_get_roots_d (new_ctx, &roots, NULL);

      for (i = 0; i < n - 1; i++)
	{
	  cplx_t starting_point;
	  int j = i / 2;
	  cplx_set (starting_point, roots[j]);

	  if (i % 2 == 0)
	    cplx_mul_eq (starting_point, (cplx_Im (roots[j]) > 0) ? rot_up : rot_down);
	  else
	    cplx_div_eq (starting_point, (cplx_Im (roots[j]) > 0) ? rot_up : rot_down);

	  mps_approximation_set_fvalue (ctx, approximations[i], starting_point);
	  mps_approximation_set_frad (ctx, approximations[i], DBL_MAX);
	}

      if (mp->level % 2 == 1)
	{
	  cplx_set_d (rot_up, -1.0, 0.0);
	  mps_approximation_set_fvalue (ctx, approximations[i], rot_up);
	}
      else
	mps_approximation_set_fvalue (ctx, approximations[i], cplx_one);

      cplx_vfree (roots);
      mps_mandelbrot_poly_free (new_ctx, MPS_POLYNOMIAL (new_mp));
      mps_context_free (new_ctx);
    }
}
