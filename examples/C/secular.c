/*
 * Example code for libmps
 *
 * This code computes the zeros of the secular equation 
 *
 *  S(x) = 1 / (x - 2) + 3 / (x - 4) + 5 / (x - 6) - 1
 *
 * using its representation as floating point, high precision
 * floating point, and rational coefficients. 
 *
 * Can be compiled with:
 *   gcc -o secular -lm -lgmp -lmps secular.c
 *
 * Author: Leonardo Robol <leonardo.robol@unipi.it>
 */


#include <mps/mps.h>
#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include <string.h>

int
main (int argc, char **argv)
{
  mps_context * ctx = NULL;
  int i;
  long wp = 128;

  cplx_t *af = NULL, *bf = NULL;
  mpcf_t *am = NULL, *bm = NULL;
  mpq_t *aqr = NULL, *bqr = NULL;
  mpq_t qzero;

  cplx_t *roots = NULL;
  mpcf_t  *mroots = NULL;
  mps_secular_equation * p = NULL;

  /* First version: floating point coefficients. */
  printf("== Computing the roots using floating point coefficients \n");
  af = cplx_valloc(3);
  bf = cplx_valloc(3);

  for (i = 0; i < 3; i++) {
    cplx_set_d(af[i], 2.0 * i + 1, 0.0);
    cplx_set_d(bf[i], 2.0 * i + 2.0, 0.0);
  }

  ctx = mps_context_new();
  p = mps_secular_equation_new (ctx, af, bf, 3);

  mps_context_set_input_poly (ctx, MPS_POLYNOMIAL (p));

  mps_mpsolve(ctx);

  mps_context_get_roots_d (ctx, &roots, NULL);

  for (i = 0; i < 3; i++) {
    printf("Root %d: %f + %f I\n", i, cplx_Re(roots[i]), cplx_Im(roots[i]));
  }

  free (roots);
  mps_polynomial_free (ctx, MPS_POLYNOMIAL (p));
  mps_context_free(ctx);
  
  free(af);
  free(bf);

  /* Version 2: higher precision coefficients and approximations */
  printf("\n\n== Computing the roots using multiprecision floating-point coefficients (asking 40 digits)\n");
  am = mpcf_valloc(3);
  bm = mpcf_valloc(3);

  for (i = 0; i < 3; i++) {
    mpcf_init2(am[i], wp);
    mpcf_set_ui(am[i], 2U * i + 1U, 0U);
    mpcf_init2(bm[i], wp);
    mpcf_set_ui(bm[i], 2U * i + 2U, 0U);    
  }

  ctx = mps_context_new();
  p = mps_secular_equation_new_raw (ctx, 3);

  for (i = 0; i < 3; i++) {
    mps_secular_equation_set_coefficient_f (ctx, p, i, am[i], bm[i]);
  }  

  mps_context_set_input_poly (ctx, MPS_POLYNOMIAL (p));
  mps_context_set_output_goal (ctx, MPS_OUTPUT_GOAL_APPROXIMATE);
  mps_context_set_output_prec (ctx, 40);
  mps_mpsolve(ctx);

  mps_context_get_roots_m (ctx, &mroots, NULL);

  for (i = 0; i < 3; i++) {
    printf("Root %d: ", i);
    mpcf_out_str(stdout, 10, wp, mroots[i]);
    printf("\n");
    mpcf_clear(mroots[i]);

    mpcf_clear(am[i]);
    mpcf_clear(bm[i]);    
  }

  free (mroots); mroots = NULL;
  mps_polynomial_free (ctx, MPS_POLYNOMIAL (p));
  mps_context_free(ctx);
  
  free(am);
  free(bm);

  /* Version 3: higher precision coefficients and approximations */
  printf("\n\n== Computing the roots using rational coefficients (asking 240 digits)\n");
  aqr = mpq_valloc(3);
  bqr = mpq_valloc(3);
  mpq_init(qzero);

  mpq_set_ui (qzero, 0U, 1U);

  for (i = 0; i < 3; i++) {
    mpq_init(aqr[i]);
    mpq_set_ui(aqr[i], 2U * i + 1U, 1U);
    mpq_init(bqr[i]);
    mpq_set_ui(bqr[i], 2U * i + 2U, 1U);    
  }

  ctx = mps_context_new();
  p = mps_secular_equation_new_raw (ctx, 3);

  for (i = 0; i < 3; i++) {
    mps_secular_equation_set_coefficient_q (ctx, p, i, aqr[i], qzero, bqr[i], qzero);
  }

  mpq_clear(qzero);

  mps_context_set_input_poly (ctx, MPS_POLYNOMIAL (p));
  mps_context_set_output_goal (ctx, MPS_OUTPUT_GOAL_APPROXIMATE);
  mps_context_set_output_prec (ctx, 240);
  mps_mpsolve(ctx);

  mps_context_get_roots_m (ctx, &mroots, NULL);

  for (i = 0; i < 3; i++) {
    printf("Root %d: ", i);
    mpcf_out_str(stdout, 10, wp, mroots[i]);
    printf("\n");
    
    mpcf_clear(mroots[i]);

    mpq_clear(aqr[i]);
    mpq_clear(bqr[i]);
  }

  free (mroots);
  mps_polynomial_free (ctx, MPS_POLYNOMIAL (p));
  mps_context_free(ctx);
  
  free(aqr);
  free(bqr);    
  
  return EXIT_SUCCESS;
}
