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
	s->dsum_ab = cdpe_valloc(n);
	s->dsum_bz = cdpe_valloc(n);

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
	cdpe_vfree(s->dsum_ab);
	cdpe_vfree(s->dsum_bz);

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
		cplx_sub(sec->sum_bz[i], x, sec->bfpc[i]);

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


void
mps_secular_dnewton(mps_status* s, cdpe_t x, rdpe_t rad, cdpe_t corr, boolean * again) {

	int i;
	cdpe_t ctmp, pol, fp;
	rdpe_t rtmp, rtmp2, rtmp3;

	/* Get pointer to the mps_secular_equation */
	mps_secular_equation* sec = (mps_secular_equation*) s->user_data;

	for(i = 0; i < sec->n; i++) {
		/* Compute (z - b_i) and store it in sec->sum_bz[i] */
		cdpe_sub(sec->dsum_bz[i], sec->bdpc[i], x);

		/* Compute a_i / (z - b_i) and store it in sec->sum_ab[i] */
		cdpe_inv(ctmp, sec->dsum_bz[i]);
		cdpe_mul(sec->dsum_ab[i], sec->adpc[i], ctmp);

		/* Overwrite (z - b_i) with a_i / (z-b_i)^2. Now ctmp
		 * is 1 / z - b_i */
		cdpe_mul_eq(sec->dsum_bz[i], ctmp);
	}

	/* Compute polynomial divided for the product of (b_j - z) */
	cdpe_set(pol, cdpe_zero);
	for(i = 0; i < sec->n; i++) {
		cdpe_add_eq(pol, sec->dsum_ab[i]);
	}
	cdpe_sub_eq(pol, cdpe_one);

	/* Compute the first derivative of the polynomial divided for
	 * (b_j - z) */
	cdpe_set(fp, cdpe_zero);
	for(i = 0; i < sec->n; i++) {
		cdpe_sub_eq(fp, sec->dsum_bz[i]);
	}

	/* Compute newton correction */
	cdpe_div(corr, pol, fp);

	/* Compute radius of inclusion
	 * TODO: Check the right way to compute this */

	/* Start with computing p(|z|), remembering that
	 * sec->sum_ad contains the values a_i / (z-b_i).
	 * Compute their modulus and then sum them together. */
	rdpe_set(rad, rdpe_zero);
	for(i = 0; i < sec->n; i++) {
		cdpe_mod(rtmp, sec->dsum_ab[i]);
		rdpe_add_eq(rad, rtmp);
	}

	/* Radius is eps * 4n * p(|z|) */
	rdpe_mul_eq_d(rad, (double) sec->n * 4 * DBL_EPSILON);

	/* Compute |p(z)| , |p'(z)| and |z| and check if we need to continue */
	cdpe_mod(rtmp, pol);
	cdpe_mod(rtmp2, fp);
	cdpe_mod(rtmp3, x);
	rdpe_mul_eq_d(rtmp3, DBL_EPSILON);
	if( (rdpe_lt(rad, rtmp) ) ||
			(rdpe_gt(rtmp2, rtmp3)) ) {
		*again = false;
	} else { *again = true; }
}

