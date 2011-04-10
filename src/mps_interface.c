/*
 * mps_interface.c
 *
 *  Created on: 05/apr/2011
 *      Author: Leonardo Robol <robol@poisson.phc.unipi.it>
 */

/**
 * 	@file
 * 	@brief Implementation of the routines to interact with MPSolve
 * 	as a library.
 */

#include <mps/mps.h>
#include <mps/mps_poly.h>
#include <mps/mps_link.h>
#include <gmp.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/**
 * @brief Allocate polynomial related variables directly in mps_status.
 *
 * This routine allocate
 */
void mps_allocate_poly_inplace(mps_status* s, int n) {

	int i;

	/* If n is provided, then we should allocate variables for a polynomial
	 * of degree n. If it is <= 0 then we can assume that is already set
	 * to the right vaule.
	 */
	if (n > 0) {
		s->deg = s->n = n;
	}

	if (s->DOLOG)
	   fprintf(s->logstr, "Allocating polynomial in place\n");

	s->data_type = (char*) malloc(sizeof(char) * 3);

	s->spar = boolean_valloc(s->deg + 2);

	s->fpr = double_valloc(s->deg + 1);
	s->fpc = cplx_valloc(s->deg + 1);

	s->dpr = rdpe_valloc(s->deg + 1);
	s->dpc = cdpe_valloc(s->deg + 1);

	s->mip_r = mpz_valloc(s->deg + 1);
	s->mip_i = mpz_valloc(s->deg + 1);
	for (i = 0; i <= s->deg; i++) {
		mpz_init(s->mip_r[i]);
		mpz_init(s->mip_i[i]);
	}

	s->mqp_r = mpq_valloc(s->deg + 1);
	s->mqp_i = mpq_valloc(s->deg + 1);
	for (i = 0; i <= s->deg; i++) {
		mpq_init(s->mqp_r[i]);
		mpq_init(s->mqp_i[i]);
	}

	s->mfpr = mpf_valloc(s->deg + 1);
	for (i = 0; i <= s->deg; i++)
		mpf_init2(s->mfpr[i], s->prec_in);
	s->mfpc = mpc_valloc(s->deg + 1);
	for (i = 0; i <= s->deg; i++)
		mpc_init2(s->mfpc[i], s->prec_in);

}

/**
 * @brief Allocate a new mps_status struct with default
 * options.
 */
mps_status* mps_status_new() {
	/* Allocate the new mps_status and load default options */
	mps_status* s = (mps_status*) malloc(sizeof(mps_status));
	mps_set_default_values (s);

	/* Set default streams */
	s->instr = stdin;
	s->outstr = stdout;
	s->logstr = stdout;

	/* Set standard precision */
	s->prec_out = 53;
	s->prec_in = 0;

	return s;
}

/**
 * @brief Set active polynomial as a real floating point coefficient
 * polynomial of degree <code>n</code> with coefficient exactly
 * determined by components of vector coeff.
 *
 * Precisely, if \f${\mathrm coeff}\f$ is a vector of \f$n+1\f$ components,
 * \f[
 *   p(x) = \sum_{i = 0}^{n} {\mathrm coeff}_i  x^i
 * \f]
 */
int mps_status_set_poly_d(mps_status* s, cplx_t* coeff, long unsigned int n) {

	int i;

	/* Allocate space for a polynomial of degree n */
	mps_allocate_poly_inplace(s, n);

	/* Set type to a dense, real, floating point polynomial */
	s->data_type[0] = 'd';
	s->data_type[1] = 'c';
	s->data_type[2] = 'f';

	/* Fill polynomial coefficients */
	for(i = 0; i <= n; i++) {
		mpc_set_cplx(s->mfpc[i], coeff[i]);
	}

	/* Allocate space for computation related data */
	mps_allocate_data(s);

	return 0;
}

/**
 * @brief Set active polynomial as a integer coefficient
 * polynomial of degree <code>n</code> with coefficient exactly
 * determined by components of vector coeff.
 *
 * Precisely, if \f${\mathrm coeff}\f$ is a vector of \f$n+1\f$ components,
 * \f[
 *   p(x) = \sum_{i = 0}^{n} {\mathrm coeff}_i  x^i
 * \f]
 */
int mps_status_set_poly_i(mps_status* s, int* coeff, long unsigned int n) {

	int i;

	/* Allocate data in mps_status to hold the polynomial of degree n */
	mps_allocate_poly_inplace(s, n);

	/* Dense, real, integer coefficients */
	s->data_type[0] = 'd';
	s->data_type[1] = 'r';
	s->data_type[2] = 'i';

	/* Fill polynomial */
	for(i = 0; i <= n; i++) {
		mpz_set_si(s->mip_r[i], coeff[i]);
	}

	/* Allocate data for the computation */
	mps_allocate_data(s);

	return 0;
}

/**
 * @brief Set <code>roots[i]</code> to the i-th root of the polynomial
 * and (if it is not <code>NULL</code>) <code>radius[i]</code>
 * to the i-th inclusion radius.
 */
int mps_get_roots_d(mps_status* s, cplx_t* roots, double* radius) {
	int i;
	for(i = 0; i < s->n; i++) {

	if (radius != NULL) {
		  if (s->lastphase == float_phase ||
				  s->lastphase == dpe_phase) {
			  radius[i] = s->frad[i];
		  }
		  else {
			  radius[i] = rdpe_get_d(s->drad[i]);
		  }

	}

	if (s->lastphase == mp_phase) {
		  mpc_get_cplx(roots[i], s->mroot[i]);
	  } else if (s->lastphase == float_phase) {
		  cplx_set(roots[i], s->froot[i]);
	  } else if (s->lastphase == dpe_phase) {
		  cdpe_get_x(roots[i], s->droot[i]);
	  }
  }
  return 0;
}
