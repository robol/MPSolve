/*
 * secular.c
 *
 *  Created on: 11/apr/2011
 *      Author: leonardo
 */

#include <mps/interface.h>
#include <mps/secular.h>
#include <mps/core.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>


int main(int argc, char** argv) {

	int i;

	/* Create a new empty mps_status */
	mps_status* s = mps_status_new();

	/* Create a new secular equation with some random coefficients */
	unsigned int n;
        if (argc < 2)
            n = 5;
        else
            n = atoi(argv[1]);


        
	cplx_t* a_coefficients = cplx_valloc(n);
	cplx_t* b_coefficients = cplx_valloc(n);

	srand(time(NULL));
	for(i = 0; i < n; i++) {
		cplx_set_d(a_coefficients[i], drand(), drand());
		cplx_set_d(b_coefficients[i], drand(), drand());
//		cplx_set_d(a_coefficients[i], (double) i+4, 0);
//		cplx_set_d(b_coefficients[i], ((double) i) + 1.5, 0);
                cplx_set_d(a_coefficients[i], pow(-1, (double) i + 1), 0);
                cplx_set_d(b_coefficients[i], 1.0 / (i+1) / (i+1), 0);
	}

	/* Dump coefficients */
	printf("Coefficients: \n");
	printf("a = [ ", n);
	for(i = 0; i < n; i++) {
		printf("%f + %fi ", cplx_Re(a_coefficients[i]), cplx_Im(a_coefficients[i]));
	}
	printf("];\n");

	printf("b = [ ", n);
	for(i=0; i < n; i++) {
		printf("%f + %fi ", cplx_Re(b_coefficients[i]), cplx_Im(b_coefficients[i]));
	}
	printf(" ]; \n");

	mps_secular_equation* sec = mps_secular_equation_new(a_coefficients, b_coefficients, n);

	/* Set user polynomial with our custom functions */
	mps_status_set_poly_u(s, n,
			MPS_FNEWTON_PTR(mps_secular_fnewton),
			MPS_DNEWTON_PTR(mps_secular_dnewton),
			MPS_MNEWTON_PTR(mps_secular_mnewton));

	s->check_data_usr = MPS_CHECK_DATA_PTR(mps_secular_check_data);

	/* Set secular equation in user data */
	s->user_data = sec;

	/* Set DOLOG to true to see the output */
        if (argc < 3)
            s->DOLOG = false;
        else
            s->DOLOG = true;


	/* Solve the polynomial */
        s->goal[0] = 'a';
	s->prec_out = 17;
	mps_mpsolve(s);

	/* Output the roots */
	mps_copy_roots(s);
	mps_output(s);

        mps_status_free (s);
}
