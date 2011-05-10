/*
 * mps_secular.c
 *
 *  Created on: 10/apr/2011
 *      Author: leonardo
 */

#include <mps/core.h>
#include <mps/secular.h>
#include <mps/debug.h>
#include <mps/mt.h>
#include <float.h>
#include <mps/mpc.h>

#define pi2 6.283184

/**
 * @brief Create a new secular equation struct
 */
mps_secular_equation*
mps_secular_equation_new(cplx_t* afpc, cplx_t* bfpc, unsigned long int n) {

	int i;

	/* Allocate the space for the new struct */
	mps_secular_equation* s = (mps_secular_equation*) malloc(
			sizeof(mps_secular_equation));

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
	for (i = 0; i < n; i++) {
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

void mps_secular_equation_free(mps_secular_equation* s) {
	/* Free internal data */
	cplx_vfree(s->afpc);
	cplx_vfree(s->bfpc);

	/* ...and then release it */
	free(s);
}

void mps_secular_fstart(mps_status* s, int n, int i_clust, double clust_rad,
		double g, rdpe_t eps) {

	int i, l = s->punt[i_clust];
	double th = pi2 / n;
	double sigma;
	mps_secular_equation* sec = (mps_secular_equation*) s->user_data;

	/* Get best sigma possible */
	if (s->random_seed)
		sigma = drand();
	else {
		/* If this is the first cluster select sigma = 0. In the other
		 * case try to maximize starting points distance. */
		if (i_clust == 0) {
			sigma = s->last_sigma = 0;
		} else {
			sigma = mps_maximize_distance(s, s->last_sigma, i_clust, n);
		}
	}

	for (i = 0; i < s->n; i++) {
		cplx_set_d(s->froot[l + i], cos(i * th + sigma), sin(i * th + sigma));
		cplx_mul_eq_d(s->froot[l + i], 10 * DBL_EPSILON);
		cplx_add_eq(s->froot[l + i], sec->bfpc[l + i]);
	}

}

void mps_secular_dstart(mps_status* s, int n, int i_clust, rdpe_t clust_rad,
		rdpe_t g, rdpe_t eps) {

	int i, l = s->punt[i_clust];
	double th = pi2 / n;
	double sigma;
	mps_secular_equation* sec = (mps_secular_equation*) s->user_data;

	/* Get best sigma possible */
	if (s->random_seed)
		sigma = drand();
	else {
		/* If this is the first cluster select sigma = 0. In the other
		 * case try to maximize starting points distance. */
		if (i_clust == 0) {
			sigma = s->last_sigma = 0;
		} else {
			sigma = mps_maximize_distance(s, s->last_sigma, i_clust, n);
		}
	}

	for (i = 0; i < s->n; i++) {
		cdpe_set_d(s->droot[l + i], cos(i * th + sigma), sin(i * th + sigma));
		cdpe_mul_eq_d(s->droot[l + i], DBL_EPSILON);
		cdpe_add_eq(s->droot[l + i], sec->bdpc[l + i]);
	}
}

void mps_secular_fnewton(mps_status* s, cplx_t x, double *rad, cplx_t corr,
		mps_boolean * again) {

	int i;
	cplx_t ctmp, ctmp2, pol, fp, sumb;
	double apol;
	*again = true;

	mps_secular_equation* sec = (mps_secular_equation*) s->user_data;

	cplx_set(pol, cplx_zero);
	cplx_set(fp, cplx_zero);
	cplx_set(sumb, cplx_zero);
	*rad = 0;
	apol = 0;

	for (i = 0; i < sec->n; i++) {
		/* Compute z - b_i */
		cplx_sub(ctmp, x, sec->bfpc[i]);

		/* Compute (z-b_i)^{-1} */
		cplx_inv_eq(ctmp);

		if (cplx_eq_zero(ctmp)) {
			cplx_set(corr, cplx_zero);
			*rad = 0;
			return;
		}

		/* Compute sum of (z-b_i)^{-1} */
		cplx_add_eq(sumb, ctmp);

		/* Compute a_i / (z - b_i) */
		cplx_mul(ctmp2, sec->afpc[i], ctmp);

		/* Add a_i / (z - b_i) to pol */
		cplx_add_eq(pol, ctmp2);
		apol += cplx_mod(ctmp2);

		/* Compute a_i / (z - b_i)^2 */
		cplx_mul_eq(ctmp2, ctmp);

		/* Add it to fp */
		cplx_sub_eq(fp, ctmp2);
	}

	/* Compute secular function */
	cplx_sub_eq(pol, cplx_one);

	/* Compute newton correction */
	cplx_mul(ctmp, pol, sumb);
	cplx_add_eq(fp, ctmp);

	if (!cplx_eq(fp, cplx_zero)) {
		cplx_div(corr, pol, fp);
	} else {
		cplx_set(corr, pol);
	}

	// apol += 1;
	apol *= (sec->n + 3) * DBL_EPSILON;

	if (cplx_mod(pol) < apol) {
		*again = false;
	}

	/* Radius is n * newt_corr */
	*rad = cplx_mod(corr) * sec->n;

	if (*again && (cplx_mod(corr) < cplx_mod(x) * DBL_EPSILON)) {
		*again = false;
	}
}

void mps_secular_dnewton(mps_status* s, cdpe_t x, rdpe_t rad, cdpe_t corr,
		mps_boolean * again) {

	int i;
	*again = true;

	mps_secular_equation* sec = (mps_secular_equation*) s->user_data;

	cdpe_t pol, fp, sumb, ctmp, ctmp2;
	rdpe_t rtmp, rtmp2, apol;

	cdpe_set(pol, cdpe_zero);
	cdpe_set(fp, cdpe_zero);
	cdpe_set(sumb, cdpe_zero);
	rdpe_set(rad, rdpe_zero);
	rdpe_set(apol, rdpe_zero);

	for (i = 0; i < sec->n; i++) {
		/* Compute z - b_i */
		cdpe_sub(ctmp, x, sec->bdpc[i]);

		/* Invert it, i.e. compute 1 / (z - b_i) */
		cdpe_inv_eq(ctmp);

		/* Compute sum of 1 / (z - b_i) */
		cdpe_add_eq(sumb, ctmp);

		/* Compute a / (z - b_i) and its modulus */
		cdpe_mul(ctmp2, sec->adpc[i], ctmp);
		cdpe_add_eq(pol, ctmp2);
		cdpe_mod(rtmp, ctmp2);
		rdpe_add_eq(apol, rtmp);

		/* Compute a / (z - b_i)^2 and add it to the first derivative */
		cdpe_mul_eq(ctmp2, ctmp);
		cdpe_sub_eq(fp, ctmp2);
	}

	/* Compute poly */
	cdpe_sub_eq(pol, cdpe_one);

	/* Compute correction */
	cdpe_mul(ctmp, pol, sumb);
	cdpe_add_eq(fp, ctmp);

	if (!cdpe_eq(fp, cdpe_zero)) {
		cdpe_div(corr, pol, fp);
	} else {
		cdpe_set(corr, pol);
	}

	/* Compute radius as n * | corr | */
	cdpe_mod(rad, corr);
	rdpe_mul_eq_d(rad, s->n);

	/* Compute \sum_i | a_i / (z - b_i) | + 1
	 * and check if the secular equation is smaller
	 * than this multiplied for n * epsilon */
	rdpe_add_eq(apol, rdpe_one);
	rdpe_mul_eq_d(apol, (sec->n + 3) * DBL_EPSILON);

	cdpe_mod(rtmp, pol);
	//	if (rdpe_lt(rtmp, apol))
	//		*again = false;

	/* If newton correction is less than
	 * the modules of |x| multiplied for
	 * for epsilon stop */
	if (*again) {
		/* Computation of |x| and |corr| */
		cdpe_mod(rtmp, corr);
		cdpe_mod(rtmp2, x);
		rdpe_mul_eq_d(rtmp2, sec->n * DBL_EPSILON);

		/* If |corr| < |x| * DBL_EPSILON then stop */
		if (rdpe_lt(rtmp, rtmp2))
			*again = false;
	}

}

void mps_secular_mnewton(mps_status* s, mpc_t x, rdpe_t rad, mpc_t corr,
		mps_boolean * again) {

	int i;

	/* Set again to true. If the convergence will be proved
	 * during the iteration it will be set to false */
	*again = true;

	/* Get a pointer to the secular equation */
	mps_secular_equation* sec = (mps_secular_equation*) s->user_data;

	/* Declare temporary variables */
	mpc_t sumb, pol, fp, ctmp, ctmp2;
	mpf_t ftmp, ftmp2;
	rdpe_t rtmp, rtmp2, apol;

	/* Set working precision */
	mpc_init2(sumb, s->mpwp);
	mpc_init2(pol, s->mpwp);
	mpc_init2(fp, s->mpwp);
	mpc_init2(ctmp, s->mpwp);
	mpc_init2(ctmp2, s->mpwp);
	mpf_init2(ftmp, s->mpwp);
	mpf_init2(ftmp2, s->mpwp);

	/* Adjust precision of coefficients */
	if (s->mpwp != mpc_get_prec(sec->ampc[0])) {
		for (i = 0; i < sec->n; i++) {
			mpc_set_prec(sec->ampc[i], s->mpwp);
			mpc_set_prec(sec->bmpc[i], s->mpwp);
		}
	}

	/* Set some starting values */
	mpc_set_d(sumb, 0, 0);
	mpc_set_d(pol, 0, 0);
	mpc_set_d(fp, 0, 0);
	rdpe_set(apol, rdpe_zero);

	for (i = 0; i < sec->n; i++) {
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
		rdpe_add_eq(apol, rtmp);

		/* Add a_i / (z - b_i) to pol */
		mpc_add_eq(pol, ctmp2);

		/* Compute a_i / (z - b_i)^2 */
		mpc_mul_eq(ctmp2, ctmp);

		/* Add it to fp */
		mpc_sub_eq(fp, ctmp2);
	}

	/* Subtract one from pol */
	mpc_sub_eq_ui(pol, 1, 0);

	/* Compute correction */
	mpc_mul_eq(sumb, pol);
	mpc_add_eq(fp, sumb);
	if (!mpc_eq_zero(fp)) {
		mpc_div(corr, pol, fp);
	} else {
		mpc_set(corr, pol);
	}

	/* Compute radius */
	mpc_mod(ftmp, corr);
	mpf_get_rdpe(rad, ftmp);
	rdpe_mul_eq_d(rad, (double) sec->n);

	/* Compute modulus of pol */
	mpc_mod(ftmp, pol);
	mpf_get_rdpe(rtmp, ftmp);

	/* If |pol| < (n+3) * eps * |apol|
	 * stop */
	rdpe_add_eq(apol, rdpe_one);
	rdpe_mul_eq_d(apol, (sec->n + 3) * s->prec_out);
	if (rdpe_lt(rtmp, apol))
		*again = false;

	/* Check if newton correction is less than
	 * the modules of x for s->prec_out, and if
	 * that's the case, stop. */
	if (*again) {
		mpc_mod(ftmp, x);
		mpf_get_rdpe(rtmp, ftmp);
		rdpe_mul_eq(rtmp, s->eps_out);
		mpc_mod(ftmp, corr);
		mpf_get_rdpe(rtmp2, ftmp);

		if (rdpe_lt(rtmp2, rtmp)) {
			*again = false;
		}
	}

}

void mps_secular_check_data(mps_status* s, char* which_case) {
	/* While we can't found a good criterion to check
	 * the possibility to start in pure floating point we
	 * use the DPE version. */
	*which_case = 'f';
}

/**
 * @brief MPSolve main function for the secular equation solving.
 */
void mps_secular_mpsolve_ga(mps_status* s) {
    int i;

    /* Set initial cluster structure as no cluster structure. */
    mps_cluster_reset(s);

    /* Select initial approximations using the custom secular
     * routine. */
    mps_secular_fstart(s, s->n, 0, 0.0, 0.0, s->eps_out);
}
