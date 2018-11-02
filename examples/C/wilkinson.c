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
#define STR_BUF_SZ 256

int
main ()
{
  long int m = 0, k = 0, M = W_DEGREE;
  mpz_t *phi, *phi0;
  mps_monomial_poly *p;
  mps_context *ctx;
  cplx_t cf;
  mpq_t one, m_one, zero, cf_q;
  char  str_buf[STR_BUF_SZ];
  int ret_val;
  
  phi = malloc((M + 1)*sizeof(mpz_t));
  phi0 = malloc((M + 1)*sizeof(mpz_t));
  
  for (m = 0; m < M + 1; m++) {
    mpz_init(phi[m]);
    mpz_init(phi0[m]);
  }
  mpz_set_si(phi[0], 1);
  for(m = 0; m < M; m++) {
    for(k = 0; k < m + 1; k++) {
      mpz_set(phi0[k],phi[k]);
    }
    mpz_set_si(phi[0], 0);
    for(k = 0; k < m + 1; k++) {
       mpz_set(phi[k + 1], phi0[k]);
    }
    for(k = 0; k < m + 1; k++) {
      mpz_mul_si (phi0[k], phi0[k], m + 1);
      mpz_sub (phi[k], phi[k], phi0[k]);
    }
  }
  /* for (k = 0; k < M + 1; k++) {gmp_printf("%Zd\n", phi[k]);};  return 0; */
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
    ret_val = gmp_snprintf (str_buf, STR_BUF_SZ, "%Zd", phi[m]);
    if(ret_val > STR_BUF_SZ) {
      fprintf(stderr, "String buffer too short: %d instead of %d",
        STR_BUF_SZ, ret_val);
      return -1;
    }
    mps_monomial_poly_set_coefficient_s(ctx, p, m, str_buf, "0");
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
  cplx_t *results = cplx_valloc (M);

  /* Actually solve the polynomial */
  mps_mpsolve (ctx);

  /* Save roots computed in the vector results */
  mps_context_get_roots_d (ctx, &results, NULL);

  /* Print out roots */
  for (k = 0; k < M; k++)
    {
      printf("%.17e %.17e\n", cplx_Re(results[k]), cplx_Im(results[k]));
    }

  free (results);

  return EXIT_SUCCESS;
}
