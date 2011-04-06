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
#include <string.h>
#include <gmp.h>

int main(int argc, char** argv) {

    /* n is the degree of the polynomial,
     * i is used as counter */
    int n, i;

    /* Get n from command line */
    if (argc > 1) {
    	n = atoi(argv[1]);
    }

    /* If parsing failed set n = 5 */
    if (n == 0) { n = 5; }

    /* Allocate space for the coefficients and fill it */
    double* coeff = (double*) malloc(sizeof(double) * (n + 1));
    coeff[0] = -1;
    coeff[n] = 1;
    for(i = 1; i < n; i++) {
    	coeff[i] = 0;
    }

    /* Allocate space to hold the results */
    cplx_t* results = cplx_valloc(n);
  
    /* Create a new mps_status and set the polynomial */
    mps_status* s = mps_status_new();
    mps_status_set_poly_f(s, coeff, n);

    /* Actually solve the polynomial */
    mps_mpsolve(s);

    /* Save roots computed in the vector results */
    mps_get_roots_f(s, results, NULL);

    /* Print out roots */
    for(i = 0; i < n; i++) {
    	cplx_out_str(stdout, results[i]); printf("\n");
    }

    return EXIT_SUCCESS;
}
