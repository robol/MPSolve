/*
 * mps_secular.c
 *
 *  Created on: 10/apr/2011
 *      Author: leonardo
 */

#include <mps/mps.h>
#include <stdlib.h>
#include <string.h>

/**
 * @brief Create a new secular equation struct
 */
mps_secular_equation*
mps_secular_equation_new(cplx_t* afpc, cplx_t* bfpc, unsigned long int n) {

	/* Allocate the space for the new struct */
	mps_secular_equation* s = (mps_secular_equation*) malloc(sizeof(mps_secular_equation));

	/* Copy data in the struct, so the user shall not worry about the scope of
	 * its input data.
	 */
	s->afpc = cplx_valloc(n);
	s->bfpc = cplx_valloc(n);

	memcpy(s->afpc, afpc, sizeof(cplx_t) * n);
	memcpy(s->bfpc, bfpc, sizeof(cplx_t) * n);

	/* Allocate temporary variables */
	s->sum_bz = cplx_valloc(n);
	s->sum_ab = cplx_valloc(n);

	return s;
}

void
mps_secular_equation_free(mps_secular_equation* s) {
	/* Free internal data */
	cplx_vfree(s->afpc);
	cplx_vfree(s->bfpc);

	/* ...and temporary variables */
	cplx_vfree(s->sum_bz);
	cplx_vfree(s->sum_ab);

	/* ...and then release it */
	free(s);
}

void
mps_secular_fnewton(mps_status* s, cplx_t x, double *rad, cplx_t corr, boolean * again) {

	int i;
	cplx_t ctmp, pol, fp;

	/* Get pointer to the mps_secular_equation */
	mps_secular_equation* sec = (mps_secular_equation*) s->user_data;

	for(i = 0; i < sec->n; i++) {
		/* Compute (z - b_i) and store it in sec->sum_bz[i] */
		cplx_sub(ctmp, sec->bfpc[i], x);
		cplx_add_eq(sec->sum_bz[i], ctmp);

		/* Compute a_i / (z - b_i) and store it in sec->sum_ab[i] */
		cplx_inv(ctmp, sec->sum_bz[i]);
		cplx_mul(sec->sum_ab[i], sec->afpc[i], ctmp);

		/* Overwrite (z - b_i) with a_i / (z-b_i)^2. Now ctmp
		 * is 1 / z - b_i */
		cplx_mul_eq(sec->sum_bz[i], ctmp);
	}

	/* Compute polynomial divided for the product of (b_j - z) */
	cplx_set(pol, cplx_zero);
	for(i = 0; i < sec->n; i++) {
		cplx_add_eq(pol, sec->sum_ab[i]);
	}
	cplx_sub_eq(pol, cplx_one);

	/* Compute the first derivative of the polynomial divided for
	 * (b_j - z) */
	cplx_set(fp, cplx_zero);
	for(i = 0; i < sec->n; i++) {
		cplx_sub_eq(fp, sec->sum_bz[i]);
	}

	/* Compute newton correction */
	cplx_div(corr, pol, fp);

	/* Compute radius of inclusion
	 * TODO: Check the right way to compute this */

	/* Start with computing p(|z|), remembering that
	 * sec->sum_ad contains the values a_i / (z-b_i).
	 * Compute their modulus and then sum them together. */
	*rad = 0;
	for(i = 0; i < sec->n; i++) {
		*rad += cplx_mod(sec->sum_ab[i]);
	}

	/* Radius is eps * n * p(|z|) */
	*rad = (*rad) * (double) sec->n * 4 * DBL_EPSILON;

	/* Check if we need to continue */
	if ((*rad < cplx_mod(pol)) &&
			(cplx_mod(fp) > DBL_EPSILON * cplx_mod(x))) {
		*again = false;
	} else { *again = true; }
}
