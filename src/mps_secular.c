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

	int i;

	/* Allocate the space for the new struct */
	mps_secular_equation* s = (mps_secular_equation*) malloc(sizeof(mps_secular_equation));

	/* Copy data in the struct, so the user shall not worry about the scope of
	 * its input data.
	 */
	s->afpc = cplx_valloc(n);
	s->bfpc = cplx_valloc(n);

	/* Allocate complex dpe coefficients of the secular equation */
	s->adpc = cdpe_valloc(n);
	s->bdpc = cdpe_valloc(n);

	/* Allocate multiprecision complex coefficients of the secular equation */
	s->ampc = mpc_valloc(n);
	s->bmpc = mpc_valloc(n);


	/* Copy the complex coefficients passed as argument */
	for(i = 0; i < n; i++) {
		/* a_i coefficients */
		cplx_set(s->afpc[i], afpc[i]);

		cdpe_init(s->adpc[i]);
		cdpe_set_x(s->adpc[i], afpc[i]);

		mpc_init2(s->ampc[i], 53);
		mpc_set_cplx(s->ampc[i], afpc[i]);

		/* b_i coefficients */
		cplx_set(s->bfpc[i], bfpc[i]);

		cdpe_init(s->bdpc[i]);
		cdpe_set_x(s->bdpc[i], bfpc[i]);

		mpc_init2(s->bmpc[i], 53);
		mpc_set_cplx(s->bmpc[i], bfpc[i]);
	}

	s->n = n;

	return s;
}

void
mps_secular_equation_free(mps_secular_equation* s) {
	/* Free internal data */
	cplx_vfree(s->afpc);
	cplx_vfree(s->bfpc);

	/* ...and then release it */
	free(s);
}

void
mps_secular_fnewton(mps_status* s, cplx_t x, double *rad, cplx_t corr, boolean * again) {

	int i;
	cplx_t ctmp, ctmp2, pol, fp, sumb;

	mps_secular_equation* sec = (mps_secular_equation*) s->user_data;

	cplx_set(pol,  cplx_zero);
	cplx_set(fp,   cplx_zero);
	cplx_set(sumb, cplx_zero);
	*rad = 0;

	for(i = 0; i < sec->n; i++) {
		/* Compute z - b_i */
		cplx_sub(ctmp, x, sec->bfpc[i]);

		/* Compute (z-b_i)^{-1} */
		cplx_inv_eq(ctmp);

		/* Compute sum of (z-b_i)^{-1} */
		cplx_add_eq(sumb, ctmp);

		/* Compute a_i / (z - b_i) */
		cplx_mul(ctmp2, sec->afpc[i], ctmp);

		/* Add a_i / (z - b_i) to pol */
		cplx_add_eq(pol, ctmp2);

		/* Compute a_i / (z - b_i)^2 */
		cplx_mul_eq(ctmp2, ctmp);

		/* Add it to fp */
		cplx_sub_eq(fp, ctmp2);
	}

	/* Compute secular function */
	cplx_sub_eq(pol, cplx_one);

//	printf("|pol| = %f\n", cplx_mod(pol));

//	if (cplx_mod(pol) < DBL_EPSILON * cplx_mod(fp) * cplx_mod(x)) {
	if(cplx_mod(pol) < DBL_EPSILON * cplx_mod(x)) {
		*again = false;
		return;
	}

	/* Compute newton correction */
	// cplx_div(corr, pol, fp);
	cplx_mul(ctmp, pol, sumb);
	cplx_add_eq(ctmp, fp);

	if (cplx_ne(ctmp, cplx_zero)) {
		cplx_div(corr, pol, ctmp);
	} else {
		cplx_div(corr, pol, fp);
	}

//	printf("Correction = "); cplx_outln(corr);

	/* Compute radius of inclusion
	 * TODO: Check the right way to compute this */

	/* Radius is eps * n * p(|z|) */
	*rad = (*rad) * (double) sec->n * 4 * DBL_EPSILON;

	*again  = true;
}

void
mps_secular_dnewton(mps_status* s, cdpe_t x, rdpe_t rad, cdpe_t corr, boolean * again) {

	int i;

	mps_secular_equation* sec = (mps_secular_equation*) s->user_data;

	cdpe_t pol, fp, sumb, ctmp, ctmp2;
	rdpe_t rtmp, rtmp2, xmod;

	cdpe_set(pol,  cdpe_zero);
	cdpe_set(fp,   cdpe_zero);
	cdpe_set(sumb, cdpe_zero);
	rdpe_set(rad,  rdpe_zero);

	for(i = 0; i < sec->n; i++) {
		/* Compute z - b_i */
		cdpe_sub(ctmp, x, sec->bdpc[i]);

		/* Compute modulus of x and |x| - b_i */
		cdpe_mod(xmod, x);

		/* Invert it, i.e. compute 1 / (z - b_i) */
		cdpe_inv_eq(ctmp);

		/* Compute sum of 1 / (z - b_i) */
		cdpe_add_eq(sumb, ctmp);

		/* Compute a / (z - b_i) and its modulus */
		cdpe_mul(ctmp2, sec->adpc[i], ctmp);
		cdpe_add_eq(pol, ctmp2);
		cdpe_mod(rtmp, ctmp2);
		rdpe_add_eq(rad, rtmp);

		/* Compute a / (z - b_i)^2 and add it to the first derivative */
		cdpe_mul_eq(ctmp2, ctmp);
		cdpe_sub_eq(fp, ctmp2);
	}

	/* Compute poly */
	cdpe_sub_eq(pol, cdpe_one);

	/* Compute radius */
	rdpe_add_eq(rad, rdpe_one);
	rdpe_mul_eq_d(rad, 4 * sec->n * DBL_EPSILON);

	/* Check if |p(z)| < |p'(z)| * |x| * epsilon */
	cdpe_mod(rtmp, fp);
	cdpe_mod(rtmp2, x);
	rdpe_mul_eq(rtmp, rtmp2);
	rdpe_mul_eq_d(rtmp, DBL_EPSILON);
	cdpe_mod(rtmp2, pol);

//	printf("Polynomial = "); cdpe_outln(pol);

	if(rdpe_lt(rtmp2, rtmp)) {
		cdpe_set(corr, cdpe_zero);
		*again = false;
		return;
	}


	/* Compute correction */
	cdpe_div(corr, fp, pol);
	cdpe_add_eq(corr, sumb);
	cdpe_inv_eq(corr);

//	printf("Correction = "); cdpe_outln(corr);


	/* Compute poly modulus. If it is less than
	 * epsilon * p(|z|) than stop */
	cdpe_mod(rtmp, pol);

	/* Stop if radius get small */
	if (rdpe_lt(rtmp, rad)) {
		*again = false;
	} else {
		*again = true;
	}

	*again = true;
}

void mps_secular_mnewton(mps_status* s, mpc_t x, rdpe_t rad, mpc_t corr, boolean * again) {

	int i;

	/* Get a pointer to the secular equation */
	mps_secular_equation* sec = (mps_secular_equation*) s->user_data;

	/* Declare temporary variables */
	mpc_t sumb, pol, fp, ctmp, ctmp2;
	mpf_t ftmp;
	rdpe_t rtmp, rtmp2;

	/* Set working precision */
	mpc_init2(sumb,  s->mpwp);
	mpc_init2(pol,   s->mpwp);
	mpc_init2(fp,    s->mpwp);
	mpc_init2(ctmp,  s->mpwp);
	mpc_init2(ctmp2, s->mpwp);
	mpf_init2(ftmp , s->mpwp);

	/* Adjust precision of coefficients */
	if (s->mpwp != mpc_get_prec(sec->ampc[0])) {
		for(i = 0; i < sec->n; i++) {
			mpc_set_prec(sec->ampc[i], s->mpwp);
			mpc_set_prec(sec->bmpc[i], s->mpwp);
		}
	}

	/* Set some starting values */
	mpc_set_d(sumb, 0, 0);
	mpc_set_d(pol,  0, 0);
	mpc_set_d(fp,   0, 0);
	rdpe_set(rad, rdpe_zero);

	for(i = 0; i < sec->n; i++) {
		/* Compute z - b_i */
		mpc_sub(ctmp, x, sec->bmpc[i]);

		/* Compute (z-b_i)^{-1} */
		mpc_inv_eq(ctmp);

		/* Multiply sum of (z-b_i)^{-1} */
		mpc_add_eq(sumb, ctmp);

		/* Compute a_i / (z - b_i) and its modulus */
		mpc_mul(ctmp2, sec->ampc[i], ctmp);
		mpc_mod(ftmp, ctmp2);
		mpf_get_rdpe(rtmp, ftmp);
		rdpe_add_eq(rad, rtmp);

		/* Add a_i / (z - b_i) to pol */
		mpc_add_eq(pol, ctmp2);

		/* Compute a_i / (z - b_i)^2 */
		mpc_mul_eq(ctmp2, ctmp);

		/* Add it to fp */
		mpc_sub_eq(fp, ctmp2);
	}

	/* Subtract one from pol */
	mpc_sub_eq_ui(pol, 1, 0);

	/* Compute modulus of |p(z)| */
	mpc_mod(ftmp, fp);
	mpf_get_rdpe(rtmp, ftmp);
	mpc_mod(ftmp, x);
	mpf_get_rdpe(rtmp2, ftmp);
	rdpe_mul_eq(rtmp, rtmp2);
	rdpe_mul_eq_d(rtmp, DBL_EPSILON);
	mpc_mod(ftmp, pol);
	mpf_get_rdpe(rtmp2, ftmp);
	if (rdpe_le(rtmp2, rtmp)) {
		*again = false;
		mpc_set_d(corr, 0, 0);
		return;
	} else {
		printf("Polynomial mod. = "); rdpe_outln(rtmp2);
	}

	/* Compute correction */
	mpc_div(corr, fp, pol);
	mpc_add_eq(corr, sumb);
	mpc_inv_eq(corr);



	/* Compute radius */
	rdpe_mul_eq_d(rad, (double) sec->n * DBL_EPSILON * 4);

	/* Compute modulus of pol */
	mpc_mod(ftmp, pol);
	mpf_add_eq_ui(ftmp, 1);
	mpf_get_rdpe(rtmp, ftmp);

	if(rdpe_lt(rad, rtmp)) {
		*again = false;
	} else {
		*again = true;
	}

}

void mps_secular_check_data(mps_status* s, char* which_case) {
	*which_case = 'f';
}
