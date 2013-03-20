/***********************************************************
**                       GMP Tools                        **
**                      Version 2.0                       **
**                                                        **
**             Written by Giuseppe Fiorentino             **
**                 (fiorent@dm.unipi.it)                  **
***********************************************************/

/**
 * @file
 * @brief Additional functions used to complete the GMP package with what
 * is needed in MPSolve.
 */

#ifndef __GMPTOOLS_H__
#define __GMPTOOLS_H__

#include <gmp.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C"
{
#endif

/**********************************************
*                  MPZ_T                      *
**********************************************/

#define mpz_Val(Z)              (*Z)
#define mpz_Move(Z1, Z2)        (*Z1 = *Z2)

/* missing functions */
#ifndef mpz_swap
  void mpz_swap (mpz_t z1, mpz_t z2);
#endif
#ifndef mpz_tstbit
  int mpz_tstbit (mpz_t z, unsigned long int pos);
#endif
#define mpz_get_bit(Z, N)       mpz_tstbit(Z, N)

#define mpz_mul_eq(Z1, Z2)      mpz_mul(Z1, Z1, Z2)
#define mpz_add_eq(Z1, Z2)      mpz_add(Z1, Z1, Z2)

/* vector support functions */
#define mpz_valloc(N)         (mpz_t *) malloc((N) * sizeof(mpz_t))
  void mpz_vinit (mpz_t v[], unsigned long int size);
  void mpz_vclear (mpz_t v[], unsigned long int size);
#define mpz_vfree(V)          free(V)

/**********************************************
*                  MPQ_T                      *
**********************************************/

#define mpq_Val(Q)            (*Q)
#define mpq_Move(Q1, Q2)      (*Q1 = *Q2)

/* missing functions */
#ifndef mpq_swap
  void mpq_swap (mpq_t q1, mpq_t q2);
#endif

/* I/O */
#ifndef mpq_out_str
  void mpq_out_str (FILE * stream, int base, mpq_t q);
#endif

/* vector support functions */
#define mpq_valloc(N)         (mpq_t *) malloc((N) * sizeof(mpq_t))
  void mpq_vinit (mpq_t v[], unsigned long int size);
  void mpq_vclear (mpq_t v[], unsigned long int size);
#define mpq_vfree(V)          free(V)

/**********************************************
*                  MPF_T                      *
**********************************************/

#define mpf_Val(F)            (*F)
#define mpf_Move(F1, F2)      (*F1 = *F2)

/* missing functions */
#ifndef mpf_swap
  void mpf_swap (mpf_t f1, mpf_t f2);
#endif
  void mpf_set_2dl (mpf_t f, double d, long int l);
  void mpf_get_2dl (double *d, long int *l, mpf_t f);
  long int mpf_size_2 (mpf_t f);

/* missing operators */
#define mpf_inv(R, F)         mpf_ui_div(R, 1, F)
#define mpf_sqr(R, F)         mpf_mul(R, F, F)
  void mpf_add_si (mpf_t r, mpf_t f, long int i);
  void mpf_sub_si (mpf_t r, mpf_t f, long int i);
  void mpf_si_sub (mpf_t r, long int i, mpf_t f);
  void mpf_mul_si (mpf_t r, mpf_t f, long int i);
  void mpf_div_si (mpf_t r, mpf_t f, long int i);
#ifndef mpf_pow_ui
  void mpf_pow_ui (mpf_t r, mpf_t f, unsigned long int i);
#endif
  void mpf_pow_si (mpf_t r, mpf_t f, long int i);

/* op= style operators for mpf_t */
#define mpf_neg_eq(F)         mpf_neg(F, F)
#define mpf_inv_eq(F)         mpf_ui_div(F, 1, F)
#define mpf_sqr_eq(F)         mpf_mul(F, F, F)
#define mpf_sqrt_eq(F)        mpf_sqrt(F, F)
#define mpf_add_eq(F1, F2)    mpf_add(F1, F1, F2)
#define mpf_add_eq_ui(F, I)   mpf_add_ui(F, F, I)
#define mpf_add_eq_si(F, I)   mpf_add_si(F, F, I)
#define mpf_sub_eq(F1, F2)    mpf_sub(F1, F1, F2)
#define mpf_sub_eq_ui(F, I)   mpf_sub_ui(F, F, I)
#define mpf_sub_eq_si(F, I)   mpf_sub_si(F, F, I)
#define mpf_ui_sub_eq(F, I)   mpf_ui_sub(F, I, F)
#define mpf_si_sub_eq(F, I)   mpf_si_sub(F, I, F)
#define mpf_mul_eq(F1, F2)    mpf_mul(F1, F1, F2)
#define mpf_mul_eq_ui(F, I)   mpf_mul_ui(F, F, I)
#define mpf_mul_eq_si(F, I)   mpf_mul_si(F, F, I)
#define mpf_mul_eq_2exp(F, I) mpf_mul_2exp(F, F, I)
#define mpf_div_eq(F1, F2)    mpf_div(F1, F1, F2)
#define mpf_div_eq_ui(F, I)   mpf_div_ui(F, F, I)
#define mpf_div_eq_si(F, I)   mpf_div_si(F, F, I)
#define mpf_ui_div_eq(F, I)   mpf_ui_div(F, I, F)
#define mpf_si_div_eq(F, I)   mpf_si_div(F, I, F)
#define mpf_div_eq_2exp(F, I) mpf_div_2exp(F, F, I)
#define mpf_pow_eq_si(F, I)   mpf_pow_si(F, F, I)

#define mpf_is_zero_p(F)      (mpf_sgn(F) ? 0 : 1)

/* vector support functions */
#define mpf_valloc(N)         (mpf_t *) malloc((N) * sizeof(mpf_t))
  void mpf_vinit (mpf_t v[], unsigned long int size);
  void mpf_vinit2 (mpf_t v[], unsigned long int size, unsigned long int prec);
  void mpf_vclear (mpf_t v[], unsigned long int size);
#define mpf_vfree(V)          free(V)



/*
 * End of extern "C" {
 *   ...
 * }
 */
#ifdef __cplusplus
}
#endif

#endif
