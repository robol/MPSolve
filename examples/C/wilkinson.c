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
#include <complex.h>
#include <time.h>
#include <math.h>

#define W_DEGREE 20


void wilkinson_coefficients(long int M, mpz_t **phi_out);
void sort_roots(long M, mpc_t *roots, rdpe_t *radii);

int
main (int argc, char **argv)
{
  long int m = 0, M = 0;
  double err_max = 0, err_d = 0, t_elapsed = 0;
  clock_t t_begin, t_end;
  mpz_t *phi;
  mpq_t cf, zero;
  mps_monomial_poly *p;
  mps_context *ctx;
  mpc_t *roots = NULL, err;
  rdpe_t *radii = NULL;
  cplx_t err_c, root_c;


  if (argc > 1) {
    M = atoi (argv[1]);
  } else {
    M = W_DEGREE;
  }
  wilkinson_coefficients(M, &phi);
  
  /* for (m = 0; m < M + 1; m++) {gmp_printf("%Zd\n", phi[m]);};  return 0; */
  /* Create the context and set the monomial coefficients */
  ctx = mps_context_new ();
  p = mps_monomial_poly_new (ctx, M);

  mps_context_select_algorithm(ctx, MPS_ALGORITHM_SECULAR_GA);

  mpq_init(cf);
  mpq_init(zero);
  for(m = 0; m < M + 1; m++) {
    mpq_set_z ( cf, phi[m]);
    mps_monomial_poly_set_coefficient_q(ctx, p, m, cf, zero);
  }
  /* Set the input polynomial */
  mps_context_set_input_poly (ctx, MPS_POLYNOMIAL (p));
  
  /* Set the output precision */
  mps_context_set_output_prec (ctx, 53);
  mps_context_set_output_goal (ctx, MPS_OUTPUT_GOAL_APPROXIMATE);
  /* Find the roots */
  t_begin = clock();
  mps_mpsolve (ctx);
  t_end = clock();
  /* Read the roots from context */
  mps_context_get_roots_m (ctx, &roots, &radii);
  /* Sort roots in the order of increasing real part */
  sort_roots(M, roots, radii);
  /* Print out roots */
  printf("Roots:\t\t\t\t\tError\t\tInclusion\n");
  printf("real part\t imaginary part\t\t\t\tradius\n");
  printf("=================================================================\n");
  mpc_init2(err, mpc_get_prec (roots[0]));
  for (m = 0; m < M; m++) {
    mpc_sub_ui(err, roots[m], m + 1, 0);
    mpc_get_cplx(err_c, err);
    err_d = sqrt(cplx_Re(err_c)*cplx_Re(err_c) + cplx_Im(err_c)*cplx_Im(err_c));
    err_max = err_d > err_max ? err_d : err_max;
    mpc_get_cplx(root_c, roots[m]);
    printf("%.3e\t%.3e\t\t%.3e\t%.3e\n", cplx_Re(root_c), cplx_Im(root_c),
      err_d, rdpe_get_d(radii[m]));
  }
  printf("=================================================================\n");
  printf("\t\t\t\tmax(Error) = %.3e\n",err_max);
  t_elapsed = (double)(t_end - t_begin) / CLOCKS_PER_SEC;
  printf("\n\tRoot finding took %ld min %.3f sec\n\n",
    (long)floor(t_elapsed/60), t_elapsed - 60*floor(t_elapsed/60));
  /* 
     Free the allocated memory 
  */
  mpq_clear(cf);
  mpq_clear(zero);
  mpc_clear(err);
  mps_monomial_poly_free(ctx, (mps_polynomial *)p);
  mps_context_free(ctx);
  for (m = 0; m < M + 1; m++) {mpz_clear(phi[m]);}
  free(phi);
  free (roots);
  /* free(roots); */
  free(radii);
  /* mpz_clear() */
  return EXIT_SUCCESS;
}

void
wilkinson_coefficients(long int M, mpz_t **phi_out)
{
  long int k = 0, m = 0;
  mpz_t  *phi0 = NULL, *phi = NULL;

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
  for (m = 0; m < M + 1; m++) {mpz_clear(phi0[m]);}
  free(phi0);
  *phi_out = phi;
}

/*
  Sort complex roots and the corresponding error radii by the
  real part of the roots
 */
typedef struct {
  mpc_t root;
  rdpe_t radius;
} val_wrapper;
int
compar(const void *a, 
  const void *b) {
  cplx_t r1, r2;
  mpc_get_cplx(r1, ((val_wrapper*)a)->root);
  mpc_get_cplx(r2, ((val_wrapper*)b)->root);
  int res = cplx_Re(r1) > cplx_Re(r2) ? 1 : -1;
  if (cplx_Re(r1) == cplx_Re(r2)) {res = 0;}
  return res;
}
void
sort_roots(long M, mpc_t *roots, rdpe_t *radii)
{
  long m;
  val_wrapper *vals = malloc(M * sizeof(val_wrapper));
  long int wp = mpc_get_prec (roots[0]);

  for (m = 0; m < M; m++) {
    mpc_init2(vals[m].root, wp);
    mpc_set(vals[m].root, roots[m]);
    rdpe_set(vals[m].radius, radii[m]);
  }
  qsort(vals, M, sizeof(val_wrapper), compar);
  for (m = 0; m < M; m++) {
    mpc_set(roots[m], vals[m].root);
    rdpe_set(radii[m], vals[m].radius);
  }
  for (m = 0; m < M; m++) {
    mpc_clear(vals[m].root);
  }
  free(vals);
}
