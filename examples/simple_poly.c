/*
 * Example code for libmps
 *
 * This code computes the 5-th roots of unity using the
 * mps_mpsolve() routine and print them to stdout.
 *
 * Can be compiled with:
 *   gcc -o simple_poly -lm -lgmp -lmps simple_poly.c
 *
 * Author: Leonardo Robol <robol@mail.dm.unipi.it>
 */


#include <mps/mps.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gmp.h>

int main(int argc, char** argv) {

	/* Declare variables with coefficients of the polynomial */
    double coeff[] = { -1, 0, 0, 0, 0, 1 };

    /* n is the degree and we need a cplx_t vector with dimension
     * n to hold the roots.
     */
    int n = 5, i;
    cplx_t* results = cplx_valloc(5);
  
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
