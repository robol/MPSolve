/*
 * Example code for libmps
 *
 * This code computes the roots of Wilkinson polynomial using
 * mps_mpsolve() routine and print them to stdout.
 *
 * Can be compiled with:
 *   gcc wilkinson.c  -o wilkinson -lm -lgmp -lmps -lpthread
 *
 * Author: Mikhail Kagalenko <kagalenko.m.b@rsreu.ru>
 */


#include <mps/mps.h>
#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include <string.h>

#define W_DEGREE 20

int
main ()
{
  /* n is the degree of the polynomial,
   * i is used as counter */
  long int m = 0, k = 0, M = W_DEGREE;
  long int phi[W_DEGREE + 1] = {0}, phi0[W_DEGREE + 1] = {0};
  mps_monomial_poly *p;
  mps_context *ctx;
  cplx_t cf;
  mpq_t one, m_one, zero, cf_q;

  phi[0] = 1;
  for(m = 0; m < M; m++) {
    for(k = 0; k < m + 1; k++) {
      phi0[k] = phi[k];
    }
    phi[0] = 0;
    for(k = 0; k < m + 1; k++) {
      phi[k + 1] = phi0[k];
    }
    for(k = 0; k < m + 1; k++) {
      phi[k] = phi[k] - (m + 1)*phi0[k];
    }
  }
  /* for (k = 0; k < M + 1; k++) {printf("%ld\n", phi[k]);};  return 0; */
  mpq_init (one);
  mpq_init (m_one);
  mpq_init (zero);
  mpq_init (cf_q);
  
  mpq_set_si (one, 1, 1);
  mpq_set_si (m_one, -1, 1);
  mpq_set_si (zero, 0, 1);

  
  ctx = mps_context_new ();
  p = mps_monomial_poly_new (ctx, M);


  mps_context_select_algorithm(ctx, MPS_ALGORITHM_SECULAR_GA);

  for(m = 0; m < M + 1; m++) {
    /* mps_monomial_poly_set_coefficient_d (ctx, p, m, phi[m], 0); */
    mpq_set_d(cf_q, (double)phi[m]);
    mps_monomial_poly_set_coefficient_q (ctx, p, m, cf_q, zero);
  }

  /*
  for(m = 0; m < M + 1; m++) {
    mps_monomial_poly_get_coefficient_d(ctx, p, m, cf);
    printf("%ld\n",  (long)(cplx_Re(cf)));
  }
  return 0;
  */
  /* Set the input polynomial */
  mps_context_set_input_poly (ctx, MPS_POLYNOMIAL (p));

  /* Allocate space to hold the results. We check only floating point results
   * in here */
  cplx_t *results = cplx_valloc (m);

  /* Actually solve the polynomial */
  mps_mpsolve (ctx);

  /* Save roots computed in the vector results */
  mps_context_get_roots_d (ctx, &results, NULL);

  /* Print out roots */
  for (k = 0; k < m; k++)
    {
      printf("%.17e %.17e\n", cplx_Re(results[k]), cplx_Im(results[k]));
    }

  free (results);

  return EXIT_SUCCESS;
}
