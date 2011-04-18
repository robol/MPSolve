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

/**
 * @mainpage Using MPSolve as a library
 *
 * @section Installation Installing MPSolve system-wide
 * First, you need to get MPSolve. You can get the latest release via <code>git</code>
 * or download it via <code>http</code> grabbing it at http://www.dm.unipi.it/...
 * If you downloaded the source tarball this operation is pretty straightforward.
 * You can simply unpack it and then
 * @code
 *   make
 *   [sudo] make install
 * @endcode
 *
 * These commands will install the library <code>libmps.so</code> in your system library
 * directory. In this way you will be able compile your source file using a command similar
 * to
 * @code
 * gcc -o myprogram -lmps -lgmp -lm myprogram.c.
 * @endcode
 *
 * @section Interface Using the libmps interface
 *
 * The library provides some useful routine to interact with the polynomial solver. Most of
 * them are designed to handle polynomial definition and are implemented in mps_interface.c
 *
 */

#include <mps/mps.h>
#include <mps/mps_poly.h>
#include <mps/mps_link.h>
#include <gmp.h>
#include <stdlib.h>
#include <stdio.h>

/**
 * @brief Allocate polynomial related variables directly in mps_status.
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
 * @brief Set active poly as a user poly, providing routines to compute
 * newton corrections.
 *
 * This is an example of call to this function:
 * @code
 * // Set a polynomial of degree n with associated mps_status* s
 * // and use the provided routines to compute newton corrections.
 * mps_status_set_poly_u(s, n,
 *   MPS_FNEWTON_PTR(mps_secular_fnewton),
 *	 MPS_DNEWTON_PTR(mps_secular_dnewton),
 *	 MPS_MNEWTON_PTR(mps_secular_mnewton));
 * @endcode
 *
 * @param s The <code>mps_status</code> struct;
 * @param n The degree of the polynomial;
 * @param fnewton The routine that performs the computation of the newton correction
 *   in floating point. It must be of the type
 *   <code>(void*)(mps_status* s, cplx_t x, double *rad, cplx_t corr, boolean * again)</code>
 *   and can be passed to the function with the right casting using the macro
 *   <code>MPS_FNEWTON_PTR</code>.
 * @param dnewton The routine that performs the computation of the newton correction in
 *   <code>dpe</code> precision. It must be of the type
 *   <code>(void*)(mps_status* s, cdpe_t x, rdpe_t rad, cdpe_t corr, boolean * again)</code>
 *   and can be passed to the function with the right casting using the macro
 *   <code>MPS_DNEWTON_PTR</code>.
 * @param fnewton The routine that performs the computation of the newton correction in
 *   multiprecision. It must be of the type
 *   <code>(void*)(mps_status* s, mpc_t x, rdpe_t rad, mpc_t corr, boolean * again)</code>
 *   and can be passed to the function with the right casting using the macro
 *   <code>MPS_MNEWTON_PTR</code>.
 */
int mps_status_set_poly_u(mps_status* s, int n, mps_fnewton_ptr fnewton,
		mps_dnewton_ptr dnewton, mps_mnewton_ptr mnewton) {

	/* Set degree and allocate data */
	s->deg = s->n = n;
	mps_allocate_data(s);

	/* TODO: Apart from u, what should be set here? */
	s->data_type = "uri";

	/* Set functions */
	s->fnewton_usr = fnewton;
	s->dnewton_usr = dnewton;
	s->mnewton_usr = mnewton;

	return 0;
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

/**
 * @brief Get the roots computed as multiprecision complex numbers.
 */
int
mps_get_roots_m(mps_status* s, mpc_t* roots, rdpe_t* radius) {
	/* TODO: Implement mps_get_roots_d() */
	return 0;
}
