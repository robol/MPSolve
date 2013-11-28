/*
 * This file is part of MPSolve 3.0
 *
 * Copyright (C) 2001-2013, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors: 
 *   Leonardo Robol <robol@mail.dm.unipi.it>
 */

#include <mps/mps.h>

#define MPS_MANDELBROT_POLY(t) (MPS_POLYNOMIAL_CAST (mps_mandelbrot_poly, t))
#define MPS_IS_MANDELBROT_POLY(t) (mps_polynomial_check_type (t, "mps_mandelbrot_poly"))

struct mps_mandelbrot_poly { 
  mps_polynomial p;
  int level; 
}; 

typedef struct mps_mandelbrot_poly mps_mandelbrot_poly;

mps_mandelbrot_poly *mps_mandelbrot_poly_new (mps_context * ctx, int level);
void mps_mandelbrot_poly_free (mps_context * ctx, mps_polynomial *p);

void mps_mandelbrot_poly_fstart (mps_context *ctx, mps_polynomial *p);
void mps_mandelbrot_poly_dstart (mps_context *ctx, mps_polynomial *p);

void mps_mandelbrot_poly_fnewton (mps_context * s, mps_polynomial * poly, mps_approximation * root, cplx_t corr);
void mps_mandelbrot_poly_dnewton (mps_context * s, mps_polynomial * poly, mps_approximation * root, cdpe_t corr);
void mps_mandelbrot_poly_mnewton (mps_context * s, mps_polynomial * poly, mps_approximation * root, mpc_t corr, long int wp);

mps_boolean mps_mandelbrot_poly_feval (mps_context * ctx, mps_polynomial * p, cplx_t x, cplx_t value, double * error);
mps_boolean mps_mandelbrot_poly_deval (mps_context * ctx, mps_polynomial * p, cdpe_t x, cdpe_t value, rdpe_t error);
mps_boolean mps_mandelbrot_poly_meval (mps_context * ctx, mps_polynomial * p, mpc_t x, mpc_t value, rdpe_t error);

