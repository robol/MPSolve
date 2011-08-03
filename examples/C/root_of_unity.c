/*
 * Example code for libmps
 *
 * This code computes the n-th roots of unity using the
 * mps_mpsolve() routine and print them to stdout.
 *
 * Can be compiled with:
 *   gcc -o root_of_unity -lm -lgmp -lmps root_of_unity.c
 *
 * Author: Leonardo Robol <robol@mail.dm.unipi.it>
 */


#include <mps/interface.h>
#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>

int
main (int argc, char **argv)
{

  /* n is the degree of the polynomial,
   * i is used as counter */
  int n, i;

  /* Get n from command line */
  if (argc > 1)
    {
      n = atoi (argv[1]);
    }

  /* If parsing failed set n = 5 */
  if (n == 0)
    {
      n = 5;
    }

  /* Allocate space for the coefficients and fill it */
  cplx_t *coeff = cplx_valloc (n + 1);

  /* Set first and last coefficients */
  cplx_set_d (coeff[0], (double) -1, 0);
  cplx_set_d (coeff[n], (double) 1, 0);
  for (i = 1; i < n; i++)
    {
      cplx_set (coeff[i], cplx_zero);
    }

  /* Allocate space to hold the results */
  cplx_t *results = cplx_valloc (n);

  /* Create a new mps_status and set the polynomial */
  mps_status *s = mps_status_new ();
  mps_status_set_poly_d (s, coeff, n);

  /* Actually solve the polynomial */
  mps_mpsolve (s);

  /* Save roots computed in the vector results */
  mps_status_get_roots_d (s, results, NULL);

  /* Print out roots */
  for (i = 0; i < n; i++)
    {
      cplx_out_str (stdout, results[i]);
      printf ("\n");
    }

  return EXIT_SUCCESS;
}
