/*
 * This file is part of MPSolve 3.0
 *
 * Copyright (C) 2001-2013, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors: 
 *   Leonardo Robol <robol@mail.dm.unipi.it>
 */

#include <mps/mps.h>
#include <mps/chebyshev.h>

 void mps_chebyshev_poly_free(mps_context *, mps_polynomial *);
 long int mps_chebyshev_poly_raise_data (mps_context * ctx, mps_polynomial * poly, long int wp);
 void mps_chebyshev_get_leading_coefficient (mps_context * ctx, mps_polynomial * poly, mpc_t lc);

 /* This is implemented in chebyshev-evaluation.c */
 mps_boolean mps_chebyshev_poly_meval (mps_context * ctx, mps_polynomial * poly, mpc_t x, mpc_t value, rdpe_t error);

 mps_chebyshev_poly * 
 mps_chebyshev_poly_new (mps_context * ctx, int n, mps_structure structure)
 {
        mps_chebyshev_poly * poly = mps_new (mps_chebyshev_poly);

        /* Call the base constructor */
        MPS_POLYNOMIAL (poly)->degree = n;
        mps_polynomial_init (ctx, MPS_POLYNOMIAL (poly));

        /* Our implementation is still not thread safe */
        MPS_POLYNOMIAL (poly)->thread_safe = false;

        MPS_POLYNOMIAL (poly)->structure = structure;

        /* Chances are that these are not really needed. In this case we can 
         * simply leave this values set to NULL. They will be allocated on
         * demand by mps_chebyshev_poly_set_coefficient_q(). */
        poly->rational_real_coeffs = poly->rational_imag_coeffs = NULL;

        if (MPS_STRUCTURE_IS_INTEGER (structure) || MPS_STRUCTURE_IS_RATIONAL (structure))
        {
          poly->rational_real_coeffs = mps_newv (mpq_t, n + 1);
          poly->rational_imag_coeffs = mps_newv (mpq_t, n + 1);

          mpq_vinit (poly->rational_real_coeffs, n + 1);
          mpq_vinit (poly->rational_imag_coeffs, n + 1);
        }

        poly->fpc = cplx_valloc (n + 1);
        poly->dpc = cdpe_valloc (n + 1);
        poly->mfpc = mpc_valloc (n + 1);

        mpc_vinit2 (poly->mfpc, n + 1, ctx->mpwp);

        /* Construct the polynomial vtable in a proper way */
        MPS_POLYNOMIAL (poly)->free = mps_chebyshev_poly_free;
        MPS_POLYNOMIAL (poly)->raise_data = mps_chebyshev_poly_raise_data;
        MPS_POLYNOMIAL (poly)->meval = mps_chebyshev_poly_meval;
        MPS_POLYNOMIAL (poly)->get_leading_coefficient = mps_chebyshev_get_leading_coefficient;

        /* Attach the typename to this polynomial */
        MPS_POLYNOMIAL (poly)->type_name = MPS_CHEBYSHEV_POLY_TYPE_NAME;

        /* Compute the leading coefficient of the polynomial */
        mpc_init2 (poly->lc, ctx->mpwp);
        if (n > 0) {
                mpc_set_ui (poly->lc, 2U, 0U);
                mpc_pow_si (poly->lc, poly->lc, n - 1);
        }
        else {
                mpc_set_ui (poly->lc, 1U, 0U);
        }

        pthread_mutex_init (&poly->precision_mutex, NULL);

        return poly;
 }

 void
 mps_chebyshev_poly_free (mps_context * ctx, mps_polynomial * poly)
 {
        mps_chebyshev_poly * cpoly = MPS_CHEBYSHEV_POLY (poly);

        mpc_vclear (cpoly->mfpc, poly->degree + 1);
        mpc_clear (cpoly->lc);

        cplx_vfree (cpoly->fpc);
        cdpe_vfree (cpoly->dpc);
        mpc_vfree (cpoly->mfpc);

        if (MPS_STRUCTURE_IS_INTEGER (poly->structure) || MPS_STRUCTURE_IS_RATIONAL (poly->structure))
          {
            mpq_vclear (cpoly->rational_real_coeffs, poly->degree + 1);
            mpq_vclear (cpoly->rational_imag_coeffs, poly->degree + 1);

            free (cpoly->rational_real_coeffs);
            free (cpoly->rational_imag_coeffs);
          }

        free (poly);
 }

 long int
 mps_chebyshev_poly_raise_data (mps_context * ctx, mps_polynomial * poly, long int wp)
 {
        int i;
        mps_chebyshev_poly * cpoly = MPS_CHEBYSHEV_POLY (poly);

        pthread_mutex_lock (&cpoly->precision_mutex);

        /* Check if raising precision is a worth operation, or if we are already
         * to a higher preicison. */
        if (wp < mpc_get_prec (cpoly->mfpc[0])) {
            pthread_mutex_unlock (&cpoly->precision_mutex);
            return mpc_get_prec (cpoly->mfpc[0]);
        }

        mpc_set_prec (cpoly->lc, wp);
        mpc_set_ui (cpoly->lc, 2U, 0U);
        mpc_pow_si (cpoly->lc, cpoly->lc, MAX(poly->degree - 1, 0));

        /* Otherwise really increase it */
        for (i = 0; i <= poly->degree; i++) {
                mpc_set_prec (cpoly->mfpc[i], wp);
        }

        if (MPS_STRUCTURE_IS_INTEGER (poly->structure) ||
            MPS_STRUCTURE_IS_RATIONAL (poly->structure))
          {
            /* Regenerate the coefficients starting from the exact input
             * of the user. */
            for (i = 0; i <= poly->degree; i++)
              {
                mpf_set_q (mpc_Re (cpoly->mfpc[i]), cpoly->rational_real_coeffs[i]);
                mpf_set_q (mpc_Im (cpoly->mfpc[i]), cpoly->rational_imag_coeffs[i]);
              }
          }

        pthread_mutex_unlock (&cpoly->precision_mutex);

        return mpc_get_prec (cpoly->mfpc[0]);
 }

void mps_chebyshev_get_leading_coefficient (mps_context * ctx, mps_polynomial * poly, mpc_t lc)
{
        mps_chebyshev_poly * cpoly = MPS_CHEBYSHEV_POLY (poly);
        mpc_set (lc, cpoly->lc);
        mpc_mul_eq (lc, cpoly->mfpc[poly->degree]);
}


 void
 mps_chebyshev_poly_set_coefficient_f (mps_context * ctx, mps_chebyshev_poly * poly,
        int i, mpc_t coeff)
 {
        if (mpc_get_prec (coeff) > mpc_get_prec (poly->mfpc[0])) {
                mps_chebyshev_poly_raise_data (ctx, MPS_POLYNOMIAL (poly), mpc_get_prec (coeff));
        }

        mpc_set (poly->mfpc[i], coeff);
        mpc_get_cdpe (poly->dpc[i], coeff);
        mpc_get_cplx (poly->fpc[i], coeff);
 }

void 
mps_chebyshev_poly_set_coefficient_q (mps_context * ctx, mps_chebyshev_poly * poly, 
       int i, mpq_t real_part, mpq_t imag_part)
{
  mps_chebyshev_poly * cpoly = MPS_CHEBYSHEV_POLY (poly);

  mpq_set (cpoly->rational_real_coeffs[i], real_part);
  mpq_set (cpoly->rational_imag_coeffs[i], imag_part);

  mpf_set_q (mpc_Re (cpoly->mfpc[i]), real_part);
  mpf_set_q (mpc_Im (cpoly->mfpc[i]), imag_part);
}
