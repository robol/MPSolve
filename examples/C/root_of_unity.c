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


#include <mps/mps.h>
#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include <string.h>

int
main (int argc, char **argv)
{
  /* n is the degree of the polynomial,
   * i is used as counter */
  long int n = 0, i;
  mps_monomial_poly *p;
  mps_context *s;

  mpq_t one, m_one, zero;

  mpq_init (one);
  mpq_init (m_one);
  mpq_init (zero);

  mpq_set_si (one, 1, 1);
  mpq_set_si (m_one, -1, 1);
  mpq_set_si (zero, 0, 1);

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
  
  s = mps_context_new ();
  p = mps_monomial_poly_new (s, n);

  mps_monomial_poly_set_coefficient_q (s, p, 0, m_one, zero); 
  mps_monomial_poly_set_coefficient_q (s, p, n, one, zero); 

  /* mps_monomial_poly_set_coefficient_d (s, p, 0, 1, 0); */
  /* mps_monomial_poly_set_coefficient_d (s, p, n, -1, 0); */

  /* Set the input polynomial */
  mps_context_set_input_poly (s, p);

  /* Allocate space to hold the results. We check only floating point results
   * in here */
  cplx_t *results = cplx_valloc (n);

  /* Actually solve the polynomial */
  mps_mpsolve (s);

  /* Save roots computed in the vector results */
  mps_context_get_roots_d (s, results, NULL);

  /* Print out roots */
  for (i = 0; i < n; i++)
    {
      cplx_out_str (stdout, results[i]);
      printf ("\n");
    }

  return EXIT_SUCCESS;
}
